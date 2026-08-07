//! Extracted cohesive responsibility from the command façade.

use super::genotype::{equivalent_gt, parse_gt_alleles};
use super::metrics::{add_variant_stats, add_variant_stats_subtype};
use super::rows::{
    almismatch_same_locus, bk_for_row, cluster_query_filter, fn_fp_combined_row, fn_row,
    fp_like_row, split_query_primitives_with_neighbors, tp_combined_row, tp_single_side_row,
    trim_variant, unk_combined_row, unk_truth_row,
};
use super::{
    AnnotatedRow, Cluster, ComparisonConfig, Entry, Event, MAX_CLUSTER_VARIANTS, RegionState,
    SPLIT_LEFT_SHIFT_WINDOW, Side, XCMP_ENUMERATION_THRESHOLD, row_matches_variant_key,
};
use crate::adapters::vcf::{Variant, VariantKey};
use crate::domain::{Interval, TypeCounts};
use crate::engines::partial_credit;
use anyhow::{Result, bail};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub(super) fn build_clusters_with_gap(
    truth: &[Variant],
    query: &[Variant],
    cluster_gap: usize,
) -> Vec<Cluster> {
    let mut entries = Vec::new();
    entries.extend(truth.iter().cloned().map(|variant| Entry {
        side: Side::Truth,
        variant,
    }));
    entries.extend(query.iter().cloned().map(|variant| Entry {
        side: Side::Query,
        variant,
    }));
    entries.sort_by(|left, right| {
        left.variant
            .key
            .chrom
            .cmp(&right.variant.key.chrom)
            .then(left.variant.key.pos.cmp(&right.variant.key.pos))
            .then(left.variant.end_pos().cmp(&right.variant.end_pos()))
    });

    let mut clusters = Vec::new();
    let mut current: Option<Cluster> = None;
    for entry in entries {
        let start = entry.variant.key.pos;
        let end = entry.variant.end_pos();
        match &mut current {
            Some(cluster)
                if cluster.chrom == entry.variant.key.chrom
                    && start <= cluster.end.saturating_add(cluster_gap)
                    && cluster.truth.len() + cluster.query.len() < MAX_CLUSTER_VARIANTS =>
            {
                cluster.end = cluster.end.max(end);
                match entry.side {
                    Side::Truth => cluster.truth.push(entry.variant),
                    Side::Query => cluster.query.push(entry.variant),
                }
            }
            _ => {
                if let Some(cluster) = current.take() {
                    clusters.push(cluster);
                }
                let mut cluster = Cluster {
                    chrom: entry.variant.key.chrom.clone(),
                    start,
                    end,
                    truth: Vec::new(),
                    query: Vec::new(),
                };
                match entry.side {
                    Side::Truth => cluster.truth.push(entry.variant),
                    Side::Query => cluster.query.push(entry.variant),
                }
                current = Some(cluster);
            }
        }
    }

    if let Some(cluster) = current {
        clusters.push(cluster);
    }
    clusters
}

pub(super) fn process_cluster(
    cluster: &Cluster,
    reference_sequences: &BTreeMap<String, String>,
    conf_bed: Option<&[Interval]>,
    config: ComparisonConfig,
    counts: &mut BTreeMap<String, TypeCounts>,
    subtype_counts: &mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    rows: &mut Vec<AnnotatedRow>,
) -> Result<()> {
    let output_start = rows.len();
    let reference = reference_sequences
        .get(&cluster.chrom)
        .ok_or_else(|| anyhow::anyhow!("reference contig {} not found", cluster.chrom))?;
    // Class F: when a multi-allelic deletion query primitive slides left of
    // the cluster's original anchor, the BS column and BED-derived Region
    // tags must reflect the extended cluster span. Pre-compute the
    // post-shift primitive positions, extend cluster.start (and end) to
    // cover them, and proceed with the widened span. Without this both
    // chr21:44413756 / chr21:47906001 emit `BS=44413761` / `BS=47906004`
    // (the original parent anchor) while legacy emits `BS=44413756` /
    // `BS=47906001`.
    let mut min_primitive_start = cluster.start;
    let mut max_primitive_end = cluster.end;
    for query in &cluster.query {
        for primitive in split_query_primitives_with_neighbors(
            query,
            reference,
            cluster.start,
            &cluster.query,
            &cluster.truth,
        ) {
            if primitive.key.pos < min_primitive_start {
                min_primitive_start = primitive.key.pos;
            }
            let p_end = primitive.key.pos + primitive.key.ref_allele.len().max(1) - 1;
            if p_end > max_primitive_end {
                max_primitive_end = p_end;
            }
        }
    }
    let cluster: Cluster = if min_primitive_start < cluster.start || max_primitive_end > cluster.end
    {
        Cluster {
            chrom: cluster.chrom.clone(),
            start: min_primitive_start,
            end: max_primitive_end,
            truth: cluster.truth.clone(),
            query: cluster.query.clone(),
        }
    } else {
        cluster.clone()
    };
    let cluster = &cluster;
    let region_state = RegionState::from_cluster(cluster, reference, conf_bed);
    let mut truth_remaining = cluster.truth.clone();
    let mut query_remaining = cluster.query.clone();
    // Capture the row range emitted by `exact_match_pairs` so the post-
    // hap_mismatch pass can degrade `BK=lm` on UNK combined rows when the
    // cluster turns out to have hap_mismatch=false. legacy emits BK=lm
    // on outside-CONF exact-match pairs only when the cluster's block-
    // level hap-compare also disagrees; otherwise it stays at BK=`.`.
    let exact_match_pre_count = rows.len();
    exact_match_pairs(
        cluster,
        reference,
        &region_state,
        ComparisonOutputs {
            counts,
            subtype_counts,
            rows,
        },
        &mut truth_remaining,
        &mut query_remaining,
    );
    let exact_match_post_count = rows.len();
    // Prepared truth can carry the same unphased GT spelling as query. For
    // those byte-identical exact indel pairs legacy's block comparison leaves
    // an outside-CONF row at BK=`.` when the rest of the block reconciles. The
    // ordinary (phased) truth stream does not take this path.
    let identical_exact_keys = identical_gt_exact_indel_keys(cluster);

    if truth_remaining.is_empty() && query_remaining.is_empty() {
        if !config.no_hc {
            degrade_identical_exact_unk_rows(
                &mut rows[exact_match_pre_count..exact_match_post_count],
                &identical_exact_keys,
            );
        }
        set_xcmp_context(&mut rows[output_start..], "simple:match", false);
        return Ok(());
    }

    // Legacy xcmp only runs block-level haplotype comparison when ALL
    // three conditions of `hap_run` fire in `xcmp.cpp:finish_block`:
    // (a) `n_nonsnp > 0` — block contains at least one GT-selected non-
    // SNP allele on either side; (b) `calls_1 > 0` — truth side has any
    // variant in the block; (c) `calls_2 > 0` — query side has any
    // variant in the block. A SNP-only block, a truth-only block, or a
    // query-only block never reaches hap-compare: each unmatched record
    // keeps its per-variant FN/FP/UNK classification from simple compare
    // and BK stays `.` regardless of block shape. Rust's
    // `cluster_signature` would otherwise synthesize a spurious "block
    // mismatch" against an empty counterpart side and promote every
    // unmatched indel to BK=lm.
    let allow_haplotype_match = !config.no_hc
        && cluster_has_gt_selected_nonsnp(cluster)
        && !cluster.truth.is_empty()
        && !cluster.query.is_empty();
    // Truth variants that were exact-matched (removed from truth_remaining).
    // These are the only variants for which a matching deletion in truth
    // should suppress the spurious BK=lm via drain restoration in
    // enumerate_haplotype_assignments. Variants still present in
    // truth_remaining were not matched, so BK=lm remains correct.
    let truth_matched: Vec<Variant> = cluster
        .truth
        .iter()
        .filter(|tv| !truth_remaining.iter().any(|r| r.key == tv.key))
        .cloned()
        .collect();
    // Class C narrow relaxation: at positions where the query has both an
    // Insert and a Subst record AND truth has a multi-allelic record at
    // the same anchor whose alts cover BOTH the query's insert allele and
    // its subst allele, allow the Insert+Subst conflict in the query
    // enumeration. This mirrors legacy's reinterpretation of overlapping
    // query records as a single hetalt multi-allelic — chr21:30374435
    // query `G→GT 1/1` + `G→T 0/1` reconciles against truth `G→GT,T 2|1`
    // by placing the insert+sub combo on one hap and the insert alone on
    // the other (apply_events emits sub_alt + inserted at the shared
    // anchor). The truth-counterpart guard prevents over-firing on shapes
    // like chr21:16328989 (truth has `G→GA` only, no SNP allele in alts)
    // where legacy keeps the strict drain semantics → BK=`.`.
    let relax_positions = compute_class_c_relaxation_positions(&cluster.query, &cluster.truth);
    let signature_cluster = Cluster {
        chrom: cluster.chrom.clone(),
        start: cluster.start.saturating_sub(config.hb_expand).max(1),
        end: cluster
            .end
            .saturating_add(config.hb_expand)
            .min(reference.len()),
        truth: cluster.truth.clone(),
        query: cluster.query.clone(),
    };
    let (truth_sig, query_sig) = if allow_haplotype_match {
        (
            cluster_signature_with_limit(
                &signature_cluster,
                &cluster.truth,
                reference,
                None,
                &BTreeSet::new(),
                config.max_enum,
            )?,
            cluster_signature_with_limit(
                &signature_cluster,
                &cluster.query,
                reference,
                Some(&truth_matched),
                &relax_positions,
                config.max_enum,
            )?,
        )
    } else {
        (None, None)
    };

    // Legacy's `DiploidCompare::setRegion` flags `hap_match = true` as
    // soon as ANY enumerated (h1, h2) pair from truth's di_haps equals
    // ANY pair from query's di_haps — cross-product membership, not set
    // equality. Replicate that rule here: mismatch only when the two
    // signature sets are disjoint. `None` on either side (enumeration
    // budget exceeded) still falls through to the mismatch path because
    // legacy's `hap_fail = true` takes the same `ctype != "hap:mismatch"`
    // branch we want for non-evaluable blocks.
    let is_match = matches!(
        (&truth_sig, &query_sig),
        (Some(left), Some(right)) if left.intersection(right).next().is_some()
    );
    // Legacy's `ctype == "hap:mismatch"` fires only when the block-level
    // haplotype comparator actually ran (both signatures computed) AND
    // the two sides disagreed. A None signature — hapcmp skipped on the
    // n_nonsnp gate OR state-count cap exceeded — corresponds to
    // legacy's "simple" / "hapfail" ctypes, neither of which promotes to
    // BK=lm.
    let hap_mismatch = if allow_haplotype_match && truth_sig.is_some() && query_sig.is_some() {
        !is_match
    } else if allow_haplotype_match
        && truth_sig.is_some()
        && query_sig.is_none()
        && estimated_state_count_with_limit(&cluster.query, config.max_enum) <= config.max_enum
    {
        // Query states drained due to Insert+Subst conflict or ref-overlap (not
        // budget overflow). Budget overflow → hapfail → BK=.; state drain →
        // check truth counterpart.
        //   Some(true)  — Insert+Subst drain AND truth has the insert → BK=.
        //   Some(false) — Insert+Subst drain, no truth counterpart  → BK=lm
        //                 OR deletion-covers-insert drain            → BK=lm
        //   None        — other drain (overlapping deletions, etc.)  → BK=.
        matches!(
            query_insert_conflict_has_truth_counterpart(
                &cluster.query,
                &cluster.truth,
                &truth_remaining,
            ),
            Some(false)
        )
    } else {
        false
    };

    // Class E (chr21:35384302 chr21 case) — degrade BK=lm to `.` on UNK
    // exact-match-pair rows when the cluster's haplotype comparator also
    // matches (hap_mismatch=false). `unk_combined_row` hardcodes BK=lm at
    // emission time because exact_match_pairs runs before cluster_signature
    // is computed. Once we know the cluster has no hap-level mismatch, the
    // legacy verdict is BK=`.` (same as a TP combined row would carry
    // BK=gm). Restrict to clusters where hap-compare actually ran
    // (`allow_haplotype_match`) so SNP-only or single-side clusters keep
    // their pre-fix behaviour — those never had hap_mismatch evaluated and
    // their BK=lm hardcode is what legacy emits in those shapes.
    if allow_haplotype_match && !hap_mismatch {
        for row in &mut rows[exact_match_pre_count..exact_match_post_count] {
            if row.record.sample_values_contain(":UNK:lm:") {
                row.record.replace_sample_values(":UNK:lm:", ":UNK:.:");
            }
        }
    }

    let remainder = Cluster {
        chrom: cluster.chrom.clone(),
        start: cluster.start,
        end: cluster.end,
        truth: truth_remaining,
        query: query_remaining,
    };
    // Legacy's graph hapcmp can reconcile a hom-alt query indel with two
    // nearby heterozygous truth copies of the same edit inside a repetitive
    // block, even when the copies use different VCF anchors. The linear
    // signature enumerator intentionally does not generally collapse such
    // anchors (that would be unsound for heterogeneous insertions), so retain
    // this narrowly evidenced promotion path. chr21_preprocess_controls at
    // 15181523/15181526 is the fixture: the earlier truth copy is outside
    // CONF, the same-key 15181526 pair is a GT mismatch, and an exact shared
    // SNP anchors the block. Legacy emits the mismatch pair as TP/gm with
    // HapMatch rather than FN/FP am.
    let legacy_hap_promotions = if is_match {
        BTreeSet::new()
    } else {
        legacy_repetitive_indel_hap_promotions(cluster, &region_state)
    };
    let (mut xcmp_ctype, mut xcmp_hap_match) = if !allow_haplotype_match {
        ("simple:mismatch", false)
    } else {
        match (&truth_sig, &query_sig) {
            (Some(_), Some(_)) if is_match => ("hap:match", true),
            (Some(_), Some(_)) => ("hap:mismatch", false),
            _ if hap_mismatch => ("hap:mismatch", false),
            _ => ("hapfail:mismatch", false),
        }
    };
    if is_match {
        mark_cluster_match(
            &remainder,
            cluster,
            reference,
            &region_state,
            ComparisonOutputs {
                counts,
                subtype_counts,
                rows,
            },
        );
    } else {
        mark_cluster_mismatch(
            &remainder,
            cluster,
            hap_mismatch,
            reference,
            &region_state,
            ComparisonOutputs {
                counts,
                subtype_counts,
                rows,
            },
        );
    }
    if !legacy_hap_promotions.is_empty() {
        for row in &mut rows[exact_match_post_count..] {
            if row_matches_variant_key(row, &legacy_hap_promotions)
                && row.record.sample_values_contain(":FN:am:")
                && row.record.sample_values_contain(":FP:am:")
            {
                let key = legacy_hap_promotions
                    .iter()
                    .find(|key| row_matches_variant_key(row, &BTreeSet::from([(*key).clone()])))
                    .expect("promoted row has a matching variant key");
                let truth = cluster
                    .truth
                    .iter()
                    .find(|variant| variant.key == *key)
                    .expect("promoted truth variant exists");
                let query = cluster
                    .query
                    .iter()
                    .find(|variant| variant.key == *key)
                    .expect("promoted query variant exists");
                *row = tp_combined_row(
                    truth,
                    query,
                    reference,
                    cluster.start,
                    &region_state.row_tags(Some(truth), Some(query)),
                );
            }
        }
        xcmp_ctype = "hap:match";
        xcmp_hap_match = true;
    }
    // The provisional exact-pair `UNK:lm` becomes `UNK:.` only when no
    // residual row proves that the enclosing block remains unreconciled.
    // Preprocess-controls dot blocks contain only TP/gm sibling rows; true
    // local-mismatch blocks retain at least one residual lm/am row.
    if !config.no_hc
        && !rows[exact_match_post_count..]
            .iter()
            .any(annotated_row_has_unreconciled_allele)
    {
        degrade_identical_exact_unk_rows(
            &mut rows[exact_match_pre_count..exact_match_post_count],
            &identical_exact_keys,
        );
    }
    if xcmp_hap_match
        && rows[output_start..]
            .iter()
            .any(annotated_row_has_unreconciled_allele)
    {
        xcmp_ctype = "hap:mismatch";
        xcmp_hap_match = false;
    }
    // pre.py emits unphased, decomposed truth records. When a decomposed SNP
    // and indel share an anchor and the SNP has an exact query counterpart,
    // legacy's VariantLocationAggregator orders the SNP record first. The
    // phased non-preprocessed stream follows the ordinary decision/type sort,
    // so gate this on the all-unphased same-anchor truth shape.
    let snp_first_positions = legacy_preprocessed_snp_first_positions(cluster);
    for row in &mut rows[output_start..] {
        if snp_first_positions.contains(&row.sort_key.1) {
            row.sort_key.2 = usize::from(!annotated_row_is_snp(row));
        }
    }
    set_xcmp_context(&mut rows[output_start..], xcmp_ctype, xcmp_hap_match);
    Ok(())
}

pub(super) fn identical_gt_exact_indel_keys(cluster: &Cluster) -> BTreeSet<VariantKey> {
    cluster
        .truth
        .iter()
        .filter(|truth| {
            truth.primary_type() == "INDEL"
                && cluster
                    .query
                    .iter()
                    .any(|query| truth.key == query.key && truth.gt == query.gt)
        })
        .map(|truth| truth.key.clone())
        .collect()
}

pub(super) fn degrade_identical_exact_unk_rows(
    rows: &mut [AnnotatedRow],
    keys: &BTreeSet<VariantKey>,
) {
    for row in rows {
        if row_matches_variant_key(row, keys) && row.record.sample_values_contain(":UNK:lm:") {
            row.record.replace_sample_values(":UNK:lm:", ":UNK:.:");
        }
    }
}

pub(super) fn legacy_preprocessed_snp_first_positions(cluster: &Cluster) -> BTreeSet<usize> {
    let positions = cluster
        .truth
        .iter()
        .map(|variant| variant.key.pos)
        .collect::<BTreeSet<_>>();
    positions
        .into_iter()
        .filter(|pos| {
            let truth_at_pos = cluster
                .truth
                .iter()
                .filter(|variant| variant.key.pos == *pos)
                .collect::<Vec<_>>();
            truth_at_pos
                .iter()
                .all(|variant| variant.gt.contains('/') && !variant.gt.contains('|'))
                && truth_at_pos
                    .iter()
                    .any(|variant| variant.primary_type() == "SNP")
                && truth_at_pos
                    .iter()
                    .any(|variant| variant.primary_type() == "INDEL")
                && truth_at_pos.iter().any(|truth| {
                    truth.primary_type() == "SNP"
                        && cluster.query.iter().any(|query| {
                            query.key == truth.key && equivalent_gt(&query.gt, &truth.gt)
                        })
                })
        })
        .collect()
}

pub(super) fn annotated_row_is_snp(row: &AnnotatedRow) -> bool {
    row.record.ref_allele.len() == 1
        && row
            .record
            .alt_allele
            .split(',')
            .all(|allele| allele.len() == 1)
}

pub(super) fn legacy_repetitive_indel_hap_promotions(
    cluster: &Cluster,
    region_state: &RegionState,
) -> BTreeSet<VariantKey> {
    let has_exact_shared_anchor = cluster.truth.iter().any(|truth| {
        cluster.query.iter().any(|query| {
            truth.key == query.key
                && equivalent_gt(&truth.gt, &query.gt)
                && truth.primary_type() == "SNP"
        })
    });
    if !has_exact_shared_anchor {
        return BTreeSet::new();
    }

    cluster
        .truth
        .iter()
        .filter(|truth| {
            truth.primary_type() == "INDEL"
                && truth.is_het()
                && region_state.truth_is_conf(truth)
                && cluster.query.iter().any(|query| {
                    query.key == truth.key
                        && query.is_homalt()
                        && region_state.query_is_conf(query)
                        && cluster.truth.iter().any(|other_truth| {
                            other_truth.key != truth.key
                                && other_truth.primary_type() == "INDEL"
                                && other_truth.is_het()
                                && !region_state.truth_is_conf(other_truth)
                                && same_single_alt_indel_edit(other_truth, truth)
                        })
                })
        })
        .map(|truth| truth.key.clone())
        .collect()
}

pub(super) fn same_single_alt_indel_edit(left: &Variant, right: &Variant) -> bool {
    fn edit(variant: &Variant) -> Option<(String, String)> {
        if variant.key.alt_allele.contains(',') {
            return None;
        }
        let (_, reference, alternate) = trim_variant(
            variant.key.pos,
            &variant.key.ref_allele,
            &variant.key.alt_allele,
        );
        Some((reference, alternate))
    }

    edit(left) == edit(right)
}

pub(super) fn set_xcmp_context(rows: &mut [AnnotatedRow], ctype: &'static str, hap_match: bool) {
    for row in rows {
        row.xcmp_ctype = Some(ctype);
        row.xcmp_hap_match = hap_match;
    }
}

pub(super) fn annotated_row_has_unreconciled_allele(row: &AnnotatedRow) -> bool {
    let Some(format) = row.record.format.as_deref() else {
        return false;
    };
    let format = format.split(':').collect::<Vec<_>>();
    let Some(bk_index) = format.iter().position(|key| *key == "BK") else {
        return false;
    };
    for sample in &row.record.samples {
        let values = sample.split(':').collect::<Vec<_>>();
        if values
            .get(bk_index)
            .is_some_and(|bk| matches!(*bk, "lm" | "am"))
        {
            return true;
        }
    }
    false
}

#[cfg(test)]
pub(super) fn cluster_signature(
    cluster: &Cluster,
    variants: &[Variant],
    reference: &str,
    truth_variants: Option<&[Variant]>,
    relax_positions: &BTreeSet<usize>,
) -> Result<Option<BTreeSet<String>>> {
    cluster_signature_with_limit(
        cluster,
        variants,
        reference,
        truth_variants,
        relax_positions,
        XCMP_ENUMERATION_THRESHOLD,
    )
}

pub(super) fn cluster_signature_with_limit(
    cluster: &Cluster,
    variants: &[Variant],
    reference: &str,
    truth_variants: Option<&[Variant]>,
    relax_positions: &BTreeSet<usize>,
    max_enum: usize,
) -> Result<Option<BTreeSet<String>>> {
    // Skip enumeration when the predicted state space exceeds the legacy
    // xcmp cap. A `None` result forces the caller to treat the cluster as
    // mismatch — the same conservative fallback legacy applies when the
    // hap-block enumeration budget is blown.
    if estimated_state_count_with_limit(variants, max_enum) > max_enum {
        return Ok(None);
    }

    let segment = reference_segment(reference, cluster.start, cluster.end)?;
    let mut signatures = BTreeSet::new();
    for (hap1_events, hap2_events) in enumerate_haplotype_assignments(
        variants,
        reference,
        cluster.start,
        cluster.end,
        truth_variants,
        relax_positions,
        max_enum,
    )? {
        let Some(hap1) = apply_events(reference, cluster.start, cluster.end, &hap1_events)? else {
            continue;
        };
        let Some(hap2) = apply_events(reference, cluster.start, cluster.end, &hap2_events)? else {
            continue;
        };
        let mut pair = [hap1, hap2];
        pair.sort();
        signatures.insert(format!("{}|{}", pair[0], pair[1]));
    }
    if signatures.is_empty() {
        // An empty cluster has a well-defined signature (the reference
        // segment on both haps). A non-empty cluster with zero valid
        // haplotype strings means every assignment was invalidated by
        // apply_events — falling back to `{segment}|{segment}` here would
        // spuriously match an empty-side cluster's reference signature and
        // emit phantom TP rows. Treat as enumeration failure → mismatch.
        if !variants.is_empty() {
            return Ok(None);
        }
        signatures.insert(format!("{segment}|{segment}"));
    }
    Ok(Some(signatures))
}

/// Legacy xcmp's `finish_block` only triggers block-level haplotype
/// comparison when at least one variant in the block selects a non-SNP
/// alt via its GT (the `n_nonsnp > 0` gate in xcmp.cpp). Returns true
/// when any cluster variant's GT picks an alt whose trimmed ref/alt
/// lengths aren't both 1. Unselected alts (declared but GT=0 on both
/// haps) don't count — legacy counts only selected alleles.
/// The GT-selected non-reference alt strings for a single row.
///
/// Legacy's `SimpleDiploidCompare` computes its allele-set comparison on
/// the allele STRINGS each side's GT actually picks, not on the raw ALT
/// column entries. A multi-allelic truth record `T→C,TATC` with GT `1|0`
/// selects only `C`; the `TATC` alt declared in the row never enters the
/// per-row comparison set. Missing (`.` / `*`) and empty alts are ignored
/// to match legacy's skip of those indices in `alleles_seen`.
pub(super) fn gt_selected_nonref_alts(variant: &Variant) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    let alts: Vec<&str> = variant.key.alt_allele.split(',').collect();
    for allele_index in parse_gt_alleles(&variant.gt) {
        if allele_index == 0 {
            continue;
        }
        let Some(alt) = alts.get(allele_index.saturating_sub(1)).copied() else {
            continue;
        };
        if alt.is_empty() || alt == "." || alt == "*" {
            continue;
        }
        out.insert(alt.to_string());
    }
    out
}

pub(super) fn cluster_has_gt_selected_nonsnp(cluster: &Cluster) -> bool {
    for variant in cluster.truth.iter().chain(cluster.query.iter()) {
        let alts: Vec<&str> = variant.key.alt_allele.split(',').collect();
        let ref_allele = &variant.key.ref_allele;
        let gt_indices: std::collections::HashSet<usize> =
            parse_gt_alleles(&variant.gt).into_iter().collect();

        // Check GT-selected alleles for any non-SNP (insertion OR deletion).
        for &allele_index in &gt_indices {
            if allele_index == 0 {
                continue;
            }
            let Some(alt) = alts.get(allele_index.saturating_sub(1)).copied() else {
                continue;
            };
            if alt.is_empty() || alt == "." || alt == "*" {
                continue;
            }
            // Mirror the prefix/suffix trim `normalize_ref_alt` applies so
            // equal-length substitution chains (e.g. TAT→CAT) count as SNPs.
            let (trimmed_ref_len, trimmed_alt_len) = trimmed_primitive_lens(ref_allele, alt);
            if trimmed_ref_len != 1 || trimmed_alt_len != 1 {
                return true;
            }
        }

        // Also check non-GT-selected alleles for deletions: a multi-allelic
        // truth record like TA→AA,T (GT=1|0) has a non-selected deletion T
        // that legacy's xcmp still counts toward n_nonsnp, enabling hapcmp.
        // Non-selected insertions (e.g. T→C,TATC GT=1|0) do NOT trigger it.
        for (idx, &alt) in alts.iter().enumerate() {
            let allele_index = idx + 1;
            if gt_indices.contains(&allele_index) {
                continue; // handled above
            }
            if alt.is_empty() || alt == "." || alt == "*" {
                continue;
            }
            let (trimmed_ref_len, trimmed_alt_len) = trimmed_primitive_lens(ref_allele, alt);
            if trimmed_ref_len > trimmed_alt_len {
                // Non-GT-selected deletion allele.
                return true;
            }
        }
    }
    false
}

pub(super) fn trimmed_primitive_lens(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();
    let mut prefix = 0usize;
    while prefix < ref_bytes.len()
        && prefix < alt_bytes.len()
        && ref_bytes[prefix] == alt_bytes[prefix]
    {
        prefix += 1;
    }
    let mut suffix = 0usize;
    while suffix < ref_bytes.len().saturating_sub(prefix)
        && suffix < alt_bytes.len().saturating_sub(prefix)
        && ref_bytes[ref_bytes.len() - 1 - suffix] == alt_bytes[alt_bytes.len() - 1 - suffix]
    {
        suffix += 1;
    }
    (
        ref_bytes.len().saturating_sub(prefix + suffix),
        alt_bytes.len().saturating_sub(prefix + suffix),
    )
}

/// When query haplotype enumeration drains (states → empty), determine whether
/// the cause is an Insert+Subst same-anchor conflict AND whether the Insert at
/// that position has a matching truth variant (same pos/ref/alt).
///
/// Returns:
///   `None`        — no Insert+Subst conflict detected (budget exhaustion or
///                   ref-overlap drain, not the Insert+Subst case).
///   `Some(true)`  — conflict found AND truth has a counterpart for the query
///                   Insert at that position → expected hapfail (the query has
///                   an extra FP SNP beside a shared insert) → hap_mismatch=FALSE.
///   `Some(false)` — conflict found but no truth counterpart for any conflict-
///                   position insert → genuine mismatch → hap_mismatch=TRUE.
/// Class D BK classifier for `fn_fp_combined_row` pairs (truth+query at
/// same byte-equal alt column, GT multisets disagree). Mirrors legacy's
/// XCmpQuantify branching:
///   * Same selected allele set, different multiset (zygosity diff —
///     truth het vs query homalt of the same allele) → `am`.
///   * Different selected sets with at least one shared allele,
///     INDEL → `lm`. chr21:38861935 fixture: truth `T→TAA,TA 1|1` vs
///     query `T→TAA,TA 1/2`.
///   * Different selected sets with at least one shared allele,
///     SNP → `.` (legacy quirk; chr21:9922359 truth `T→A,C 1|0` vs
///     query `T→A,C 2/1`).
///   * Disjoint selected sets → `lm` (almismatch path).
pub(super) fn compute_paired_bk(truth: &Variant, query: &Variant) -> &'static str {
    let truth_set = gt_selected_nonref_alts(truth);
    let query_set = gt_selected_nonref_alts(query);
    if truth_set == query_set {
        return "am";
    }
    let intersect_count = truth_set.intersection(&query_set).count();
    if intersect_count == 0 {
        return "lm";
    }
    // Overlapping but unequal selected sets — type-dependent.
    if truth.primary_type() == "SNP" && query.primary_type() == "SNP" {
        return ".";
    }
    "lm"
}

/// Compute the set of positions where the query side's `Insert+Subst at
/// same anchor` conflict should be RELAXED during haplotype enumeration.
///
/// A position qualifies when:
///   * the query has at least one Insert (alt longer than ref) AND at
///     least one Subst (alt same length as ref) record AT THAT POSITION,
///   * truth has a record AT THE SAME (chrom, pos, ref) whose declared
///     alts include BOTH the query's insert allele AND the query's subst
///     allele — i.e. truth's multi-allelic mirrors the union of query's
///     overlapping single-allelic records.
///
/// chr21:30374435 example: query has `G→GT 1/1` (insert) + `G→T 0/1`
/// (subst); truth has `G→GT,T 2|1` (alts `{GT, T}`). Truth's alt set
/// covers both query alleles → relaxation fires at 30374435 → query's
/// hap pair reconciles into `(GT, T+T)` matching truth's `(GT, T)` over
/// the cluster span.
///
/// chr21:16328989 (negative test): query has `G→GA 1/1` (insert) +
/// `G→A 0/1` (subst); truth has `G→GA 1|1` (alts `{GA}`). Truth's alt
/// set lacks the SNP `A` → no relaxation → strict drain → BK=`.`.
pub(super) fn compute_class_c_relaxation_positions(
    query: &[Variant],
    truth: &[Variant],
) -> BTreeSet<usize> {
    use std::collections::BTreeMap;
    // Group query selected alts by (pos, ref_allele), splitting into
    // insert and subst buckets.
    let mut insert_alts: BTreeMap<(usize, String), BTreeSet<String>> = BTreeMap::new();
    let mut subst_alts: BTreeMap<(usize, String), BTreeSet<String>> = BTreeMap::new();
    for v in query {
        let key = (v.key.pos, v.key.ref_allele.clone());
        for alt in selected_alt_sequences(v) {
            if alt.len() > v.key.ref_allele.len() {
                insert_alts.entry(key.clone()).or_default().insert(alt);
            } else if alt.len() == v.key.ref_allele.len() {
                subst_alts.entry(key.clone()).or_default().insert(alt);
            }
        }
    }
    let mut out: BTreeSet<usize> = BTreeSet::new();
    let chrom_opt = query.first().map(|v| v.key.chrom.clone());
    let Some(chrom) = chrom_opt else {
        return out;
    };
    for (key, ialts) in &insert_alts {
        let Some(salts) = subst_alts.get(key) else {
            continue;
        };
        let combined: BTreeSet<String> = ialts.iter().chain(salts.iter()).cloned().collect();
        for tv in truth {
            if tv.key.chrom != chrom || tv.key.pos != key.0 || tv.key.ref_allele != key.1 {
                continue;
            }
            let t_alts: BTreeSet<String> = tv
                .key
                .alt_allele
                .split(',')
                .map(|s| s.to_string())
                .collect();
            if combined.is_subset(&t_alts) {
                out.insert(key.0);
                break;
            }
        }
    }
    out
}

pub(super) fn query_insert_conflict_has_truth_counterpart(
    query_variants: &[Variant],
    truth_variants: &[Variant],
    truth_remaining: &[Variant],
) -> Option<bool> {
    // Classify each query variant position as Insert (any alt longer than ref)
    // or Subst (all alts same length as ref, i.e. SNP/substitution).
    let mut insert_positions: BTreeSet<usize> = BTreeSet::new();
    let mut subst_positions: BTreeSet<usize> = BTreeSet::new();
    for v in query_variants {
        let max_alt_len = v
            .key
            .alt_allele
            .split(',')
            .map(|a| a.len())
            .max()
            .unwrap_or(0);
        if max_alt_len > v.key.ref_allele.len() {
            insert_positions.insert(v.key.pos);
        } else {
            subst_positions.insert(v.key.pos);
        }
    }
    // Positions where both an Insert and a Subst coexist.
    let conflict_positions: Vec<usize> = insert_positions
        .intersection(&subst_positions)
        .copied()
        .collect();
    if conflict_positions.is_empty() {
        // Check for deletion-covers-insert drain: a deletion with no ref (0) allele
        // in its GT spans a position where an insert also exists. Both haplotypes
        // carry the deletion, leaving no valid haplotype for the insert. Legacy
        // hapcmp can handle this and finds a mismatch → BK=lm, not BK=.
        //
        // Guard: only fire when the truth side has genuine mismatch context, i.e.
        //   (a) truth_remaining is non-empty — truth has variants unaccounted for, OR
        //   (b) cluster.truth (original) has a variant at the blocked insert position
        //       (even if it was exact-matched, its presence indicates truth expected
        //        something there).
        // If neither holds (e.g. truth's only variant was an exact-match TP with no
        // truth call at ipos), the drain is an internal query conflict unrelated to
        // truth → return None instead of a spurious Some(false).
        if insert_positions.iter().any(|&ipos| {
            // Find a deletion present on BOTH query haplotypes that covers the
            // insert anchor `ipos`. Coverage means `ipos ∈ (deletion.pos,
            // deletion.pos + ref_len - 1]` (the anchor itself is not consumed
            // by the deletion in VCF normalization, so we use strict `<`).
            let blocking_deletion = query_variants.iter().find(|v| {
                let ref_len = v.key.ref_allele.len();
                let max_alt = v
                    .key
                    .alt_allele
                    .split(',')
                    .map(|a| a.len())
                    .max()
                    .unwrap_or(0);
                ref_len > max_alt // deletion
                    && v.key.pos < ipos
                    && ipos < v.key.pos + ref_len
                    // deletion must be on both haplotypes (no 0/ref allele in GT)
                    && !v.gt.split(['/', '|']).any(|a| a == "0")
            });
            let Some(del) = blocking_deletion else {
                return false;
            };
            // Genuine-mismatch evidence must be POSITIONALLY RELATED to the
            // blocked insert: either the original truth side declared a variant
            // exactly at the blocked anchor (truth predicted the insert), or an
            // unmatched truth variant overlaps the deletion's claimed range
            // (truth disagreement falls inside the deleted span). chr21:44049606
            // (chr21_passonly) is a counter-example for the prior loose gate:
            // truth has TGATA→T at 44049663, well outside the AATGATAGATAG→A
            // deletion at 44049606..44049617 — legacy emits BK=`.` because the
            // drain is an internal query conflict unrelated to truth.
            let del_end = del.key.pos + del.key.ref_allele.len() - 1;
            let truth_at_ipos = truth_variants.iter().any(|tv| tv.key.pos == ipos);
            let truth_in_del_range = truth_remaining
                .iter()
                .any(|tv| tv.key.pos >= del.key.pos && tv.key.pos <= del_end);
            truth_at_ipos || truth_in_del_range
        }) {
            return Some(false); // deletion on both haps covers insert → genuine mismatch
        }
        return None; // not an Insert+Subst conflict drain
    }
    // For each conflict position, check whether any query Insert there has an
    // exact-ALT-match truth variant (same pos/ref/alt string).
    for pos in &conflict_positions {
        for qv in query_variants {
            if qv.key.pos != *pos {
                continue;
            }
            let max_alt_len = qv
                .key
                .alt_allele
                .split(',')
                .map(|a| a.len())
                .max()
                .unwrap_or(0);
            if max_alt_len <= qv.key.ref_allele.len() {
                continue; // not an insert at this position
            }
            // Use alt-set equality comparison: "CAA,CA" matches truth "CA,CAA"
            // (same set, different order), but "TAA" does NOT match "TAA,TA"
            // (strict subset → genuine mismatch, legacy gives BK=lm).
            let q_alts: std::collections::BTreeSet<&str> = qv.key.alt_allele.split(',').collect();
            if truth_variants.iter().any(|tv| {
                if tv.key.pos != qv.key.pos || tv.key.ref_allele != qv.key.ref_allele {
                    return false;
                }
                let t_alts: std::collections::BTreeSet<&str> =
                    tv.key.alt_allele.split(',').collect();
                q_alts == t_alts
            }) {
                return Some(true); // expected hapfail: truth has identical insert allele set
            }
        }
    }
    Some(false) // genuine mismatch: no truth counterpart for any conflict-pos insert
}

/// Upper-bound the number of haplotype assignments the enumeration will
/// produce for a variant list without actually allocating them. Heterozygous
/// unphased variants double the state count; everything else leaves it alone.
/// Returns `usize::MAX` as soon as the cap is exceeded, so callers can branch
/// without waiting for overflow.
#[cfg(test)]
pub(super) fn estimated_state_count(variants: &[Variant]) -> usize {
    estimated_state_count_with_limit(variants, XCMP_ENUMERATION_THRESHOLD)
}

pub(super) fn estimated_state_count_with_limit(variants: &[Variant], max_enum: usize) -> usize {
    let mut total: usize = 1;
    for variant in variants {
        let alleles = parse_gt_alleles(&variant.gt);
        let phased = variant.gt.contains('|');
        let multiplier = if !phased && alleles.len() == 2 && alleles[0] != alleles[1] {
            2
        } else {
            1
        };
        total = total.saturating_mul(multiplier);
        if total > max_enum {
            return usize::MAX;
        }
    }
    total
}

pub(super) fn reference_segment(reference: &str, start: usize, end: usize) -> Result<String> {
    // ASCII-only reference (DNA): byte slicing avoids allocating a full Vec<char>
    // over a 48 Mbp chromosome each time the segment is requested.
    let bytes = reference.as_bytes();
    if start == 0 || end > bytes.len() || start > end {
        bail!("invalid reference slice {start}-{end}");
    }
    std::str::from_utf8(&bytes[start - 1..end])
        .map(|s| s.to_string())
        .map_err(|err| anyhow::anyhow!("non-UTF8 reference slice: {err}"))
}

pub(super) fn enumerate_haplotype_assignments(
    variants: &[Variant],
    reference: &str,
    cluster_start: usize,
    cluster_end: usize,
    truth_variants: Option<&[Variant]>,
    relax_positions: &BTreeSet<usize>,
    max_enum: usize,
) -> Result<Vec<(Vec<Event>, Vec<Event>)>> {
    // State per haplotype: (claimed_end, subst_end, insert_end).
    //
    // claimed_end — max ref pos claimed by any Subst or Delete (using var_end =
    //               var_start + ref_allele.len() - 1, matching legacy's VCF-level
    //               comparison). Blocks new Subst/Delete that overlap.
    // subst_end   — max var_start of any Subst placed (for Insert conflict check).
    //               A new Insert at anchor A conflicts if A <= subst_end: in the
    //               DiploidReference graph you cannot traverse both a Subst edge
    //               and an Insert edge at the same ref node.
    // insert_end  — max var_start of any Insert placed (for Subst conflict check).
    //               A new Subst at pos P conflicts if P <= insert_end: same single-
    //               node constraint, opposite direction.
    //
    // Conflict rules:
    //   Subst/Delete: conflict if effective_start <= claimed_end  (ref overlap)
    //                         OR  effective_start <= insert_end   (prior Insert at anchor)
    //     where effective_start = min(Subst.pos, Delete.start) from events.
    //     Using events-based start rather than VCF var_start allows a deletion
    //     like CAAAA→C at pos P (Delete{P+1..P+4}) to coexist with a Subst at P
    //     (anchor is shared; the deletion does not consume the anchor base).
    //   Insert:       conflict if var_start <= subst_end    (prior Subst at anchor)
    //               — NOT gated by claimed_end: an Insert at anchor P can coexist
    //                 with a prior Delete whose claimed_end >= P (Insert inside a
    //                 deleted region is valid; legacy emits it unconditionally).
    //
    // Update rules (all non-ref assignments update their respective counter):
    //   Subst/Delete: claimed_end = max(claimed_end, var_end)
    //                 subst_end   = max(subst_end, var_start)  [Subst only]
    //   Insert:       claimed_end = max(claimed_end, var_end)   [keeps BK=lm for
    //                               Insert-before-Subst-at-same-pos]
    //                 insert_end  = max(insert_end, var_start)
    type HaplotypeState = (
        Vec<Event>,
        Vec<Event>,
        usize,
        usize,
        usize,
        usize,
        usize,
        usize,
    );
    let mut states: Vec<HaplotypeState> = vec![(Vec::new(), Vec::new(), 0, 0, 0, 0, 0, 0)];
    for variant in variants {
        let var_start = variant.key.pos;
        let var_end = var_start + variant.key.ref_allele.len() - 1;
        let is_insert_only = |events: &[Event]| {
            !events.is_empty() && events.iter().all(|e| matches!(e, Event::Insert { .. }))
        };
        let variant_assignments =
            variant_haplotype_assignments(variant, reference, cluster_start, cluster_end)?;
        // Pre-size `next` so we do not thrash the allocator when a variant
        // doubles the state count. Bail out early if the product blows past
        // the enumeration cap — this is a belt-and-braces check behind the
        // `estimated_state_count` gate in `cluster_signature`.
        let projected = states.len().saturating_mul(variant_assignments.len());
        if projected > max_enum {
            bail!(
                "xcmp enumeration exceeded threshold ({} projected states)",
                projected
            );
        }
        let mut next = Vec::with_capacity(projected);
        for (
            hap1_events,
            hap2_events,
            h1_claimed,
            h1_subst,
            h1_insert,
            h2_claimed,
            h2_subst,
            h2_insert,
        ) in &states
        {
            for (left_events, right_events) in &variant_assignments {
                let left_nonref = !left_events.is_empty();
                let right_nonref = !right_events.is_empty();
                let left_insert_only = is_insert_only(left_events);
                let right_insert_only = is_insert_only(right_events);
                // Effective modification start from events: min(Subst.pos, Delete.start).
                // For VCF-normalized deletions the anchor base is NOT consumed, so the
                // effective start is pos+prefix_len, not pos. Using this lets a deletion
                // CAAAA→C (Delete{P+1..P+4}) coexist with a Subst at P on the same hap.
                let event_min_start = |evts: &[Event]| -> usize {
                    evts.iter()
                        .filter_map(|e| match e {
                            Event::Subst { pos, .. } => Some(*pos),
                            Event::Delete { start, .. } => Some(*start),
                            Event::Insert { .. } => None,
                        })
                        .min()
                        .unwrap_or(var_start)
                };
                let left_eff_start = event_min_start(left_events);
                let right_eff_start = event_min_start(right_events);
                // When truth also has this exact deletion (same pos/ref/alt), every
                // query hap state would include deletion+SNP while truth has deletion
                // only — no intersection → spurious BK=lm. Restore legacy drain by
                // falling back to var_start (pre-fix behaviour) for this variant.
                let truth_has_variant = truth_variants.is_some_and(|tv| {
                    tv.iter().any(|t| {
                        t.key.pos == variant.key.pos
                            && t.key.ref_allele == variant.key.ref_allele
                            && t.key.alt_allele == variant.key.alt_allele
                    })
                });
                let left_conflict_start = if left_eff_start > var_start && truth_has_variant {
                    var_start
                } else {
                    left_eff_start
                };
                let right_conflict_start = if right_eff_start > var_start && truth_has_variant {
                    var_start
                } else {
                    right_eff_start
                };
                // Class C narrow relaxation: at relaxation positions, allow
                // Insert+Subst at the same anchor on the same hap. The
                // single-node "ref overlap with prior Subst/Delete" check
                // still applies, but the Insert↔Subst single-node block is
                // skipped — this lets the chr21:30374435 query records
                // (G→GT 1/1 + G→T 0/1) compose into truth-matching
                // haplotypes.
                let left_relax_subst = left_nonref
                    && !left_insert_only
                    && relax_positions.contains(&left_conflict_start);
                let right_relax_subst = right_nonref
                    && !right_insert_only
                    && relax_positions.contains(&right_conflict_start);
                let left_relax_insert =
                    left_nonref && left_insert_only && relax_positions.contains(&var_start);
                let right_relax_insert =
                    right_nonref && right_insert_only && relax_positions.contains(&var_start);
                // Subst/Delete conflict: ref overlap OR prior Insert at same
                // anchor. The Insert-at-same-anchor part is suppressed at
                // relaxation positions.
                if left_nonref
                    && !left_insert_only
                    && (left_conflict_start <= *h1_claimed
                        || (!left_relax_subst && left_conflict_start <= *h1_insert))
                {
                    continue;
                }
                if right_nonref
                    && !right_insert_only
                    && (right_conflict_start <= *h2_claimed
                        || (!right_relax_subst && right_conflict_start <= *h2_insert))
                {
                    continue;
                }
                // Insert conflict: prior Subst at same anchor (suppressed at
                // relaxation positions).
                if left_nonref && left_insert_only && !left_relax_insert && var_start <= *h1_subst {
                    continue;
                }
                if right_nonref
                    && right_insert_only
                    && !right_relax_insert
                    && var_start <= *h2_subst
                {
                    continue;
                }
                let mut next_hap1 = hap1_events.clone();
                next_hap1.extend_from_slice(left_events);
                let mut next_hap2 = hap2_events.clone();
                next_hap2.extend_from_slice(right_events);
                // Update claimed_end for ALL non-ref (including Insert) so that
                // a later Subst at the same pos sees conflict via claimed_end.
                // Exception: at relaxation positions, an Insert must NOT bump
                // h_claimed; otherwise a follow-on Subst at the same anchor
                // would still be blocked by the ref-overlap rule.
                let left_inserts_skip_claim = left_nonref && left_insert_only && left_relax_insert;
                let right_inserts_skip_claim =
                    right_nonref && right_insert_only && right_relax_insert;
                let new_h1_claimed = if left_nonref && !left_inserts_skip_claim {
                    (*h1_claimed).max(var_end)
                } else {
                    *h1_claimed
                };
                let new_h2_claimed = if right_nonref && !right_inserts_skip_claim {
                    (*h2_claimed).max(var_end)
                } else {
                    *h2_claimed
                };
                // Update subst_end for Subst/Delete only.
                let new_h1_subst = if left_nonref && !left_insert_only {
                    (*h1_subst).max(var_start)
                } else {
                    *h1_subst
                };
                let new_h2_subst = if right_nonref && !right_insert_only {
                    (*h2_subst).max(var_start)
                } else {
                    *h2_subst
                };
                // Update insert_end for Insert only.
                let new_h1_insert = if left_nonref && left_insert_only {
                    (*h1_insert).max(var_start)
                } else {
                    *h1_insert
                };
                let new_h2_insert = if right_nonref && right_insert_only {
                    (*h2_insert).max(var_start)
                } else {
                    *h2_insert
                };
                next.push((
                    next_hap1,
                    next_hap2,
                    new_h1_claimed,
                    new_h1_subst,
                    new_h1_insert,
                    new_h2_claimed,
                    new_h2_subst,
                    new_h2_insert,
                ));
            }
        }
        states = next;
    }
    Ok(states
        .into_iter()
        .map(|(h1, h2, _, _, _, _, _, _)| (h1, h2))
        .collect())
}

pub(super) fn variant_haplotype_assignments(
    variant: &Variant,
    reference: &str,
    cluster_start: usize,
    cluster_end: usize,
) -> Result<Vec<(Vec<Event>, Vec<Event>)>> {
    let alleles = parse_gt_alleles(&variant.gt);
    let phased = variant.gt.contains('|');
    if alleles.len() != 2 {
        bail!("unsupported genotype {}", variant.gt);
    }
    let left =
        normalized_events_for_allele(variant, alleles[0], reference, cluster_start, cluster_end)?;
    let right =
        normalized_events_for_allele(variant, alleles[1], reference, cluster_start, cluster_end)?;
    if phased || alleles[0] == alleles[1] {
        return Ok(vec![(left, right)]);
    }
    Ok(vec![(left.clone(), right.clone()), (right, left)])
}

/// Sorted multi-set of ALT-allele SEQUENCES that `variant`'s GT selects
/// (reference alleles are dropped). Legacy's `SimpleDiploidCompare`
/// builds a per-index bitmask and compares it across samples; because
/// the loader unifies the `variation[]` table across truth + query
/// samples at each locus, the per-index bitmask is implicitly a per-
/// SEQUENCE bitmask. Rust stores truth and query separately with
/// independent ALT-column ordering, so we reconstruct the legacy
/// comparison by mapping each GT index back to its alt string.
pub(super) fn selected_alt_sequences(variant: &Variant) -> Vec<String> {
    let alts: Vec<&str> = variant.key.alt_allele.split(',').collect();
    let mut out = Vec::new();
    for idx in parse_gt_alleles(&variant.gt) {
        if idx == 0 {
            continue;
        }
        let Some(alt) = alts.get(idx - 1).copied() else {
            continue;
        };
        if alt.is_empty() || alt == "." || alt == "*" {
            continue;
        }
        out.push(alt.to_string());
    }
    out.sort();
    out
}

/// Legacy-equivalent per-record simple-compare match: two variants at
/// the same (chrom, pos, ref) whose *declared* ALT sets are equal and
/// whose GT-selected sub-multisets agree. Mirrors the
/// `alleles_seen_1 == alleles_seen_2 && gtt1 == gtt2` branch of
/// `SimpleDiploidCompare.cpp:232-245` once the VariantReader-level
/// allele unification is factored out.
///
/// Requiring equal declared ALT sets (not merely equal selected sets)
/// keeps us from over-pairing records that legacy emits as two separate
/// VCF lines — e.g. truth `A→ATCTC,ATC` GT `1|0` vs query `A→ATCTC` GT
/// `0/1`. Both select `{ATCTC}`, but legacy's reader keeps them as two
/// records with different ALT columns and emits one truth-only TP row +
/// one query-only TP row.
pub(super) fn simple_compare_pairs_match(
    truth: &Variant,
    query: &Variant,
    reference: &str,
    cluster_start: usize,
    cluster_truth: &[Variant],
    cluster_query: &[Variant],
) -> bool {
    if truth.key.chrom != query.key.chrom
        || truth.key.pos != query.key.pos
        || truth.key.ref_allele != query.key.ref_allele
    {
        return false;
    }
    let truth_alts: BTreeSet<&str> = truth.key.alt_allele.split(',').collect();
    let query_alts: BTreeSet<&str> = query.key.alt_allele.split(',').collect();
    if truth_alts != query_alts {
        return false;
    }
    let truth_selected = selected_alt_sequences(truth);
    let query_selected = selected_alt_sequences(query);
    if truth_selected.is_empty() || query_selected.is_empty() {
        return false;
    }
    if truth_selected != query_selected {
        return false;
    }
    // If EITHER side's alleles primitive-split into distinct trimmed
    // anchors (mixed-length indels like `GAA→GA,G` or `TTG→TTGTG,T`),
    // legacy emits per-primitive rows at the primitives' canonical
    // positions. The TP-vs-FN classification then depends on whether
    // the *block-level* haplotype comparator reconciles the two sides —
    // so we leave those pairs unmatched here and let
    // `cluster_signature` / `mark_cluster_match` / `mark_cluster_mismatch`
    // drive the classification.
    if query_primitive_splits(
        query,
        reference,
        cluster_start,
        cluster_truth,
        cluster_query,
    ) || query_primitive_splits(
        truth,
        reference,
        cluster_start,
        cluster_truth,
        cluster_query,
    ) {
        return false;
    }
    // `gttype` equality is implied by the selected multi-set equality —
    // homalt `1|1` selects `[X, X]` vs het `0|1` selects `[X]`; these
    // differ as multi-sets even when the nonref-set matches.
    true
}

/// True iff `split_query_primitives(query)` would fan the record into
/// multiple per-primitive rows. Two paths qualify:
///
/// * Distinct trimmed anchors — the `all_same_anchor` negation at the top of
///   `split_query_primitives_with_neighbors`.
/// * Same-anchor multi-allelic insertion that splits via per-primitive
///   left-shift to orphan anchors (Class B fan-out — `try_split_same_anchor_via_shift`).
pub(super) fn query_primitive_splits(
    query: &Variant,
    reference: &str,
    cluster_start: usize,
    cluster_truth: &[Variant],
    cluster_neighbors: &[Variant],
) -> bool {
    let alleles = parse_gt_alleles(&query.gt);
    let used: BTreeSet<usize> = alleles.into_iter().filter(|a| *a > 0).collect();
    if used.is_empty() {
        return false;
    }
    let alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    let mut trimmed: Vec<(usize, String, String)> = Vec::new();
    for idx in &used {
        let Some(alt) = alts.get(*idx - 1).copied() else {
            continue;
        };
        trimmed.push(trim_variant(query.key.pos, &query.key.ref_allele, alt));
    }
    if trimmed.len() <= 1 {
        return false;
    }
    let first = (trimmed[0].0, trimmed[0].1.clone());
    let all_same_anchor = trimmed
        .iter()
        .all(|(pos, r, _)| *pos == first.0 && r == &first.1);
    if !all_same_anchor {
        return true;
    }
    // Same-anchor: only "splits" when Class B insertion fan-out would fire
    // (at least one primitive shifts to a distinct anchor that has no truth
    // representation in the cluster).
    let bases = reference.as_bytes();
    let pos_min = cluster_start.saturating_sub(SPLIT_LEFT_SHIFT_WINDOW).max(1);
    try_split_same_anchor_via_shift(
        &query.key.chrom,
        query.key.pos,
        &trimmed,
        bases,
        pos_min,
        cluster_truth,
        cluster_neighbors,
    )
    .is_some()
}

/// Compute the post-`partial_credit::left_shift` (pos, ref, alt) for a
/// trimmed primitive, falling back to the input on out-of-bounds reference
/// access. Mirrors the slide block inside
/// `split_query_primitives_with_neighbors` so the fan-out decision and the
/// emission share identical semantics.
pub(super) fn compute_shift_target(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    bases: &[u8],
    pos_min: usize,
) -> (usize, String, String) {
    if pos == 0 || ref_allele.is_empty() || pos + ref_allele.len() - 1 > bases.len() {
        return (pos, ref_allele.to_string(), alt_allele.to_string());
    }
    let mut rv = partial_credit::RefVar {
        start: pos,
        end: pos + ref_allele.len() - 1,
        alt: alt_allele.to_string(),
    };
    partial_credit::left_shift(bases, &mut rv, pos_min.max(1), true);
    let ref_start = rv.start.saturating_sub(1);
    let ref_end = rv.end;
    if ref_end < ref_start || ref_end > bases.len() || rv.start < 1 {
        return (pos, ref_allele.to_string(), alt_allele.to_string());
    }
    let new_ref: String = bases[ref_start..ref_end]
        .iter()
        .map(|&b| b.to_ascii_uppercase() as char)
        .collect();
    if new_ref.is_empty() {
        return (pos, ref_allele.to_string(), alt_allele.to_string());
    }
    (rv.start, new_ref, rv.alt)
}

/// Class B same-anchor multi-allelic insertion fan-out. When all trimmed
/// primitives share the same `(pos, ref)` anchor (legacy keeps these
/// multi-allelic by default), this helper checks whether per-primitive
/// `partial_credit::left_shift` would canonicalize at least one alt to a
/// DISTINCT anchor while leaving its sibling at the original anchor. If
/// that shifted alt lands at a position where the cluster's truth has no
/// matching record, return the per-primitive shifted (pos, ref, alt)
/// triples — the caller fans the record out into per-row primitives.
///
/// The truth-orphan gate distinguishes:
/// * chr21:21690513 — truth declares `C→CACAT` only; CACAC slides to
///   (21690501, T, TACAC) which has no truth match → fan out.
/// * chr21:40096658 — truth declares `T→TAGATAGAG` AND `C→CAGATAGAT,...`
///   at 40096650; TAGATAGAT slides to (40096650, C, CAGATAGAT) which IS
///   represented in truth → keep multi-allelic and let block-level
///   haplotype matching reconcile.
///
/// Returns `None` when the multi-allelic should stay intact (insertions
/// don't shift, all shifts collapse to the same anchor, or every shifted
/// alt has truth representation at the shifted anchor).
pub(super) fn try_split_same_anchor_via_shift(
    chrom: &str,
    variant_pos: usize,
    trimmed: &[(usize, String, String)],
    bases: &[u8],
    pos_min: usize,
    cluster_truth: &[Variant],
    cluster_neighbors: &[Variant],
) -> Option<Vec<(usize, String, String)>> {
    if trimmed.len() < 2 {
        return None;
    }
    // Insertion-only: ref length 1, alt length > 1. Same-anchor deletions
    // are handled separately by the existing all_same_anchor=false fan-out.
    if !trimmed
        .iter()
        .all(|(_, r, a)| r.len() == 1 && a.len() > r.len())
    {
        return None;
    }
    // Per-primitive slide floor: the highest cluster-query record position
    // strictly below the primitive's anchor (excluding the variant being
    // shifted). Sliding onto an already-occupied locus would duplicate
    // that record (chr21:21690513 chr21 case has a neighboring SNP
    // T→C at 21690501; the CACAC primitive would otherwise slide on
    // top of it — legacy stops one position above at 21690502 A→ACACA
    // because xcmp doesn't permit two records to share an anchor).
    let neighbor_floor = cluster_neighbors
        .iter()
        .filter(|n| n.key.pos != variant_pos)
        .map(|n| n.key.pos)
        .max();
    let effective_pos_min = match neighbor_floor {
        Some(n) => pos_min.max(n),
        None => pos_min,
    };
    let shifted: Vec<(usize, String, String)> = trimmed
        .iter()
        .map(|(pos, r, a)| compute_shift_target(*pos, r, a, bases, effective_pos_min))
        .collect();
    let first = (shifted[0].0, shifted[0].1.clone());
    let distinct = !shifted
        .iter()
        .all(|(pos, r, _)| *pos == first.0 && r == &first.1);
    if !distinct {
        return None;
    }
    // Require at least one stayer to have truth representation at the
    // ORIGINAL (un-shifted) anchor. Without this gate a multi-allelic
    // query whose alts shift apart but neither side matches truth (e.g.
    // chr21:32767041 `T→TCTCACA,TCTCTCT` in a truth-less cluster) would
    // fan out into two orphan rows, while legacy keeps it as a single
    // hetalt UNK row. The truth match at the original anchor is what
    // licenses the truth_subset_match-style emit (combined TP at
    // truth's repr plus residual at the shifted anchor).
    let any_truth_at_original = trimmed.iter().any(|(pos, r, a)| {
        cluster_truth.iter().any(|t| {
            t.key.chrom == chrom
                && t.key.pos == *pos
                && t.key.ref_allele == *r
                && t.key.alt_allele.split(',').any(|alt| alt == a)
        })
    });
    if !any_truth_at_original {
        return None;
    }
    // Block fan-out when ANY shifted alt has truth representation at the
    // shifted anchor. Legacy keeps the query as multi-allelic in those
    // cases and lets the block-level haplotype matcher reconcile.
    let any_truth_at_shifted = shifted.iter().enumerate().any(|(i, (p, r, a))| {
        // Only count primitives that actually moved.
        if *p == trimmed[i].0 && *r == trimmed[i].1 {
            return false;
        }
        cluster_truth.iter().any(|t| {
            t.key.chrom == chrom
                && t.key.pos == *p
                && t.key.ref_allele == *r
                && t.key.alt_allele.split(',').any(|alt| alt == a)
        })
    });
    if any_truth_at_shifted {
        return None;
    }
    Some(shifted)
}

/// Truth-subset match: truth selects ALL of its declared alts AND those
/// alts are a (proper) subset of query's GT-selected alleles. Legacy
/// emits a single TP/gm row at truth's representation, with query's GT
/// remapped against truth's allele indices (alleles missing from truth's
/// column collapse to ref `0`); the unmatched portion of the multi-
/// allelic query then emits as a separate per-primitive FP row.
///
/// Pinning case: chr21:27249918 — truth `CTAAATAAA→C` GT `1|0` selects
/// `[C]`; query `CTAAATAAA→C,CTAAA` GT `1/2` selects `[C, CTAAA]`.
/// Truth's `[C]` ⊊ query's `[C, CTAAA]` and the shared C allele is what
/// legacy haplotype-matches against. Output: combined TP/gm row at
/// `CTAAATAAA→C` (truth GT `1|0`, query GT `0/1`) plus orphan FP at
/// pos+4 `ATAAA→A` for the unmatched CTAAA primitive.
///
/// This path is intentionally narrow:
///   * truth declared alts ⊆ query declared alts (proper subset);
///   * truth selects all of its declared alts (otherwise the unselected
///     truth alt would not have a haplotype counterpart);
///   * truth's selected MULTISET ⊆ query's selected MULTISET — guards
///     against zygosity mismatches like truth `1|1` (homalt G×2) vs
///     query `2/1` ({G, A}) where set-subset would over-match a
///     hetalt query into a homalt truth;
///   * `query_primitive_splits(query)` is true — the multi-allelic
///     query trims into per-primitive rows at distinct anchors.
///     Otherwise legacy emits two separate rows (truth-only TP at
///     truth's repr, query-only TP at query's repr — chr21:40096658).
pub(super) fn truth_subset_match(
    truth: &Variant,
    query: &Variant,
    reference: &str,
    cluster_start: usize,
    cluster_truth: &[Variant],
    cluster_query: &[Variant],
) -> bool {
    if truth.key.chrom != query.key.chrom
        || truth.key.pos != query.key.pos
        || truth.key.ref_allele != query.key.ref_allele
    {
        return false;
    }
    let truth_alts: BTreeSet<&str> = truth.key.alt_allele.split(',').collect();
    let query_alts: BTreeSet<&str> = query.key.alt_allele.split(',').collect();
    if truth_alts == query_alts {
        // Equal alt sets are handled by simple_compare_pairs_match.
        return false;
    }
    if !truth_alts.is_subset(&query_alts) {
        return false;
    }
    let mut truth_selected = selected_alt_sequences(truth);
    if truth_selected.is_empty() {
        return false;
    }
    let truth_selected_set: BTreeSet<&str> = truth_selected.iter().map(String::as_str).collect();
    // Truth must select all of its declared alts.
    if truth_selected_set != truth_alts {
        return false;
    }
    let mut query_selected = selected_alt_sequences(query);
    if query_selected.is_empty() {
        return false;
    }
    // Multiset subset: every truth-selected occurrence must have a
    // corresponding occurrence in query. selected_alt_sequences returns
    // sorted vectors so we can sweep both with a two-pointer scan.
    truth_selected.sort();
    query_selected.sort();
    let mut qi = 0;
    for t in &truth_selected {
        loop {
            if qi >= query_selected.len() {
                return false;
            }
            if &query_selected[qi] == t {
                qi += 1;
                break;
            }
            if &query_selected[qi] > t {
                return false;
            }
            qi += 1;
        }
    }
    // Legacy only fans into a combined row + orphan primitive when the
    // multi-allelic query trims into distinct per-primitive anchors
    // (after left-shift) AND no shifted alt has truth representation at
    // the shifted anchor. Without this gate we'd over-match same-anchor
    // multi-allelics like chr21:40096658 (truth `T→TAGATAGAG` vs query
    // `T→TAGATAGAG,TAGATAGAT`) — TAGATAGAT slides to (40096650, C,
    // CAGATAGAT) which IS in cluster_truth at 40096650, so the shift
    // fan-out is suppressed and the multi-allelic stays intact.
    query_primitive_splits(
        query,
        reference,
        cluster_start,
        cluster_truth,
        cluster_query,
    )
}

/// Remap query's GT into truth's allele-index space, with alleles that
/// aren't present in truth's ALT column collapsing to ref (`0`). Unphased
/// hetalt-with-ref forms are normalised so the smaller index prints
/// first (e.g. unphased `1/2` whose second allele drops to `0` becomes
/// `0/1` rather than `1/0`). Phased GTs preserve their original order.
pub(super) fn remap_query_gt_subset(truth: &Variant, query: &Variant) -> String {
    let truth_alts: Vec<&str> = truth.key.alt_allele.split(',').collect();
    let query_alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    let separator = if query.gt.contains('|') {
        '|'
    } else if query.gt.contains('/') {
        '/'
    } else {
        return query.gt.clone();
    };
    let parts: Vec<String> = query
        .gt
        .split(['/', '|'])
        .map(|tok| {
            if tok == "." {
                return ".".to_string();
            }
            let Ok(idx) = tok.parse::<usize>() else {
                return tok.to_string();
            };
            if idx == 0 {
                return "0".to_string();
            }
            let Some(alt) = query_alts.get(idx - 1).copied() else {
                return ".".to_string();
            };
            match truth_alts.iter().position(|ta| *ta == alt) {
                Some(pos) => (pos + 1).to_string(),
                None => "0".to_string(),
            }
        })
        .collect();
    if separator == '/'
        && parts.len() == 2
        && let (Ok(a), Ok(b)) = (parts[0].parse::<i32>(), parts[1].parse::<i32>())
        && a > b
    {
        return format!("{}/{}", b, a);
    }
    parts.join(&separator.to_string())
}

/// Render a query hetalt GT in legacy xcmp's canonical "alpha-later /
/// alpha-earlier" form, expressed against the OUTPUT record's ALT
/// ordering (`output_alts` — usually truth's). For every observed
/// hetalt case the rule is:
///   - identify the two called alt strings from the query's own alt
///     column + raw GT;
///   - order them as (alpha-later, alpha-earlier);
///   - emit their positions within `output_alts`.
///
/// Covers both same-alt (`C→C,G 1/2` → `2/1`) and reordered
/// (`truth C→CA,CAA` vs `query C→CAA,CA 1/2` → `2/1`) paths. Hom
/// / het-with-ref / no-call genotypes fall through unchanged.
pub(super) fn canonical_hetalt_gt(output_alts: &str, query: &Variant) -> String {
    let gt = query.gt.as_str();
    let separator = if gt.contains('|') {
        '|'
    } else if gt.contains('/') {
        '/'
    } else {
        return gt.to_string();
    };
    let tokens: Vec<&str> = gt.split(separator).collect();
    if tokens.len() != 2 {
        return gt.to_string();
    }
    let (Some(a), Some(b)) = (
        tokens[0].parse::<usize>().ok(),
        tokens[1].parse::<usize>().ok(),
    ) else {
        return gt.to_string();
    };
    if a == 0 || b == 0 || a == b {
        return gt.to_string();
    }
    let q_alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    let (Some(&alt_a), Some(&alt_b)) = (q_alts.get(a - 1), q_alts.get(b - 1)) else {
        return gt.to_string();
    };
    // "Later" = longer allele; ties broken alphabetically. Legacy's
    // addAlleleToVariant canonicalises hetalt GTs by length-then-alpha,
    // not purely alphabetically.
    let (later, earlier) = if (alt_a.len(), alt_a) > (alt_b.len(), alt_b) {
        (alt_a, alt_b)
    } else {
        (alt_b, alt_a)
    };
    let out_alts: Vec<&str> = output_alts.split(',').collect();
    let Some(p_later) = out_alts.iter().position(|alt| *alt == later) else {
        return gt.to_string();
    };
    let Some(p_earlier) = out_alts.iter().position(|alt| *alt == earlier) else {
        return gt.to_string();
    };
    format!("{}{}{}", p_later + 1, separator, p_earlier + 1)
}

pub(super) fn is_multi_allelic(variant: &Variant) -> bool {
    variant.key.alt_allele.contains(',')
}

pub(super) fn normalized_events_for_allele(
    variant: &Variant,
    allele_index: usize,
    reference: &str,
    cluster_start: usize,
    cluster_end: usize,
) -> Result<Vec<Event>> {
    if allele_index == 0 {
        return Ok(Vec::new());
    }
    let alts: Vec<&str> = variant.key.alt_allele.split(',').collect();
    let alt = alts
        .get(allele_index.saturating_sub(1))
        .copied()
        .ok_or_else(|| {
            anyhow::anyhow!(
                "missing alt allele {} for {}",
                allele_index,
                variant.key.alt_allele
            )
        })?;
    Ok(normalize_ref_alt(
        variant.key.pos,
        &variant.key.ref_allele,
        alt,
        reference,
        cluster_start,
        cluster_end,
    ))
}

pub(super) fn normalize_ref_alt(
    pos: usize,
    ref_allele: &str,
    alt_allele: &str,
    reference: &str,
    cluster_start: usize,
    cluster_end: usize,
) -> Vec<Event> {
    let ref_chars: Vec<char> = ref_allele.chars().collect();
    let alt_chars: Vec<char> = alt_allele.chars().collect();

    let mut prefix = 0usize;
    while prefix < ref_chars.len()
        && prefix < alt_chars.len()
        && ref_chars[prefix] == alt_chars[prefix]
    {
        prefix += 1;
    }

    let mut suffix = 0usize;
    while suffix < ref_chars.len().saturating_sub(prefix)
        && suffix < alt_chars.len().saturating_sub(prefix)
        && ref_chars[ref_chars.len() - 1 - suffix] == alt_chars[alt_chars.len() - 1 - suffix]
    {
        suffix += 1;
    }

    let ref_mid = &ref_chars[prefix..ref_chars.len() - suffix];
    let alt_mid = &alt_chars[prefix..alt_chars.len() - suffix];
    let mut events = Vec::new();
    let overlap = ref_mid.len().min(alt_mid.len());
    for index in 0..overlap {
        if ref_mid[index] != alt_mid[index] {
            events.push(Event::Subst {
                pos: pos + prefix + index,
                alt: alt_mid[index],
            });
        }
    }

    if ref_mid.len() > overlap {
        events.push(Event::Delete {
            start: pos + prefix + overlap,
            end: pos + prefix + ref_mid.len() - 1,
        });
    } else if alt_mid.len() > overlap {
        let inserted: String = alt_mid[overlap..].iter().collect();
        let mut anchor = pos + prefix + overlap - 1;
        if ref_allele.len() == 1 {
            // Homopolymer anchor slide — only meaningful when the inserted
            // sequence is a homopolymer extension of the anchor's base
            // (e.g. inserting "TT" into a "T..." run is position-invariant).
            // Sliding heterogeneous insertions like "ATATA" through a
            // T-homopolymer is UNSOUND — the resulting sequence depends on
            // the exact insertion site, and two distinct variants in the
            // same run would collapse to a single anchor and conflict in
            // `apply_events` (see chr21:15671076 cluster regression).
            let bases = reference.as_bytes();
            let anchor_base = bases[anchor - 1];
            let is_homopolymer_extension = !inserted.is_empty()
                && inserted
                    .bytes()
                    .all(|b| b.eq_ignore_ascii_case(&anchor_base));
            if is_homopolymer_extension {
                while anchor < cluster_end && bases[anchor].eq_ignore_ascii_case(&anchor_base) {
                    anchor += 1;
                }
                if anchor < cluster_start {
                    anchor = cluster_start;
                }
            }
        }
        events.push(Event::Insert {
            anchor,
            seq: inserted,
        });
    }

    events
}

pub(super) fn apply_events(
    reference: &str,
    start: usize,
    end: usize,
    events: &[Event],
) -> Result<Option<String>> {
    // DNA is ASCII — byte-slice the reference instead of allocating a chr-scale
    // Vec<char> on every call.
    let bases = reference.as_bytes();
    let mut substitutions: BTreeMap<usize, char> = BTreeMap::new();
    let mut insertions: BTreeMap<usize, String> = BTreeMap::new();
    let mut deleted = BTreeSet::new();

    for event in events {
        match event {
            Event::Subst { pos, alt } => {
                // An earlier Delete covering this pos would silently swallow
                // the substitution under the output-walk loop below (the
                // deleted-pos base is skipped, so the substitution never
                // emits). That masks genuine event conflicts as "no-op"
                // haplotypes that collide with empty-query's reference
                // signature. Treat the combination as an invalid assignment.
                if deleted.contains(pos) {
                    return Ok(None);
                }
                if let Some(existing) = substitutions.get(pos) {
                    if *existing != *alt {
                        return Ok(None);
                    }
                } else {
                    substitutions.insert(*pos, *alt);
                }
            }
            Event::Insert { anchor, seq } => {
                // An insertion anchored inside a deleted span is semantically
                // "insert these bases at the junction where the deleted
                // region was" — the output-walk loop below emits
                // `insertions[pos]` regardless of whether `pos` itself is
                // deleted, so the bases survive the deletion. This mirrors
                // legacy `Haplotype::seq`, which applies a downstream
                // insertion relative to the *shifted* string after the
                // upstream delete.
                //
                // Rejecting here (prior behavior) masked legitimate
                // delete+insert pairs on the same haplotype — see the
                // chr21:16997949/16997951 cluster where query hap1 has
                // `GCA→G` followed by `A→ACG`.
                if let Some(existing) = insertions.get(anchor) {
                    if existing != seq {
                        return Ok(None);
                    }
                } else {
                    insertions.insert(*anchor, seq.clone());
                }
            }
            Event::Delete { start, end } => {
                for pos in *start..=*end {
                    if substitutions.contains_key(&pos) {
                        return Ok(None);
                    }
                    deleted.insert(pos);
                }
            }
        }
    }

    let mut output = String::new();
    for pos in start..=end {
        if !deleted.contains(&pos) {
            let base = substitutions
                .get(&pos)
                .copied()
                .unwrap_or_else(|| bases[pos - 1] as char);
            output.push(base);
        }
        if let Some(inserted) = insertions.get(&pos) {
            output.push_str(inserted);
        }
    }
    // FASTA may mix uppercase and lowercase bases (soft-masked repeats);
    // Insert events carry VCF ALT sequences which are uppercase. The byte-
    // wise signature compare then diverges when truth and query place the
    // inserted bases at different anchors inside the same homopolymer run,
    // even though the final sequence is biologically identical. Canonicalize
    // case so signatures compare biologically rather than lexically.
    Ok(Some(output.to_ascii_uppercase()))
}

pub(super) fn exact_match_pairs(
    cluster: &Cluster,
    reference: &str,
    region_state: &RegionState,
    outputs: ComparisonOutputs<'_>,
    truth_remaining: &mut Vec<Variant>,
    query_remaining: &mut Vec<Variant>,
) {
    let ComparisonOutputs {
        counts,
        subtype_counts,
        rows,
    } = outputs;
    let mut matched_truth = BTreeSet::new();
    let mut matched_query = BTreeSet::new();
    // Residual query records produced by truth-subset matches: the matched
    // alleles are emitted in the combined TP row, but the unmatched alts
    // need to fall through to the downstream per-primitive emission. We
    // collect them here and append after the matched-index filter.
    let mut residual_queries: Vec<Variant> = Vec::new();

    for (truth_index, truth) in truth_remaining.iter().enumerate() {
        // Match precedence: prefer the strict simple-compare match
        // (equal declared alt sets, equal GT-selected sub-multisets);
        // fall back to truth-subset match (truth's alts ⊊ query's
        // alts, truth selects all of its alts, truth's selected
        // alleles ⊆ query's selected). The subset case emits using
        // truth's representation with query's GT remapped.
        let exact = query_remaining.iter().enumerate().find(|(_, query)| {
            simple_compare_pairs_match(
                truth,
                query,
                reference,
                cluster.start,
                &cluster.truth,
                &cluster.query,
            )
        });
        let subset = if exact.is_none() {
            query_remaining.iter().enumerate().find(|(_, query)| {
                truth_subset_match(
                    truth,
                    query,
                    reference,
                    cluster.start,
                    &cluster.truth,
                    &cluster.query,
                )
            })
        } else {
            None
        };
        if let Some((query_index, query, is_subset)) = exact
            .map(|(i, q)| (i, q, false))
            .or_else(|| subset.map(|(i, q)| (i, q, true)))
        {
            // Per-record CONF gate (legacy `count_unk && !CONF → UNK`):
            // a same-key truth+query pair outside the confident region
            // collapses to UNK/lm on both sides, not TP/gm. The stats
            // counters mirror the row's BD: outside-CONF query becomes
            // query_unk, and the truth side drops out of truth_tp.
            let truth_in_conf = region_state.truth_is_conf(truth);
            let query_in_conf = region_state.query_is_conf(query);
            let combined_unk = !truth_in_conf && !query_in_conf;
            if truth_in_conf {
                add_variant_stats(
                    &mut counts
                        .entry(truth.primary_type().to_string())
                        .or_default()
                        .truth_tp,
                    truth,
                );
                add_variant_stats_subtype(subtype_counts, truth.primary_type(), truth, |stats| {
                    &mut stats.truth_tp
                });
            }
            if query_in_conf {
                add_variant_stats(
                    &mut counts
                        .entry(query.primary_type().to_string())
                        .or_default()
                        .query_tp,
                    query,
                );
                add_variant_stats_subtype(subtype_counts, query.primary_type(), query, |stats| {
                    &mut stats.query_tp
                });
            } else {
                add_variant_stats(
                    &mut counts
                        .entry(query.primary_type().to_string())
                        .or_default()
                        .query_unk,
                    query,
                );
                add_variant_stats_subtype(subtype_counts, query.primary_type(), query, |stats| {
                    &mut stats.query_unk
                });
            }
            // When truth and query share a byte-equal ALT column we can
            // emit a single combined row. Legacy's loader canonicalises
            // unphased multi-allelic hetalt GTs into `<later>/<earlier>`
            // form (the `addAlleleToVariant` rule in
            // `VariantLocationAggregator` with `MAX_GT=2`), so apply that
            // swap to the query GT before emission. Truth keeps its
            // source (phased) GT verbatim.
            //
            // When ALTs differ (same allele *set*, column-reordered) we
            // reuse truth's representation for the row and re-index the
            // query GT via the shared allele table — legacy does this
            // implicitly. In that branch the re-index alone produces the
            // final ordering, so NO extra swap is applied.
            //
            // simple_compare_pairs_match only fires when neither side
            // primitive-splits, so the combined or remap paths always
            // emit a single TP row.
            if is_subset {
                // Truth-subset path: truth's alts ⊊ query's selected.
                // Emit at truth's representation; remap query GT so
                // alleles missing from truth's column collapse to ref.
                let canonicalized = Variant {
                    key: VariantKey {
                        chrom: query.key.chrom.clone(),
                        pos: query.key.pos,
                        ref_allele: query.key.ref_allele.clone(),
                        alt_allele: truth.key.alt_allele.clone(),
                    },
                    qual: query.qual.clone(),
                    filter: query.filter.clone(),
                    gt: remap_query_gt_subset(truth, query),
                };
                let regions = region_state.row_tags(Some(truth), Some(query));
                if combined_unk {
                    rows.push(unk_combined_row(
                        truth,
                        &canonicalized,
                        reference,
                        cluster.start,
                        &regions,
                    ));
                } else {
                    rows.push(tp_combined_row(
                        truth,
                        &canonicalized,
                        reference,
                        cluster.start,
                        &regions,
                    ));
                }
                // Build residual queries for the unmatched alleles so the
                // downstream emission path produces orphan FP/UNK rows at
                // each primitive's natural anchor. Legacy keeps the
                // unmatched primitive (e.g. CTAAA at chr21:27249918 →
                // ATAAA→A at pos+4) as a separate row.
                //
                // Pre-trim each unmatched alt with `trim_variant` so the
                // residual's VariantKey matches the per-primitive key
                // RegionState already registered in `covered_query` (line
                // 158). Without the trim the residual lands at the parent
                // anchor with an un-trimmed alt and misses the CONF
                // registration → BD downgrades from FP to UNK.
                //
                // Class B same-anchor multi-allelic insertion (chr21:21690513
                // — `C→CACAC,CACAT`): when truth declares only CACAT and the
                // unmatched CACAC primitive shifts via `partial_credit::
                // left_shift` into a microsat-canonical anchor (21690501
                // T→TACAC), the residual must be created at the SHIFTED
                // anchor. Otherwise the orphan emits at the parent's 21690513
                // anchor — losing the byte-equality with legacy's
                // `T→TACAC` row. Apply `compute_shift_target` for
                // insertion residuals so they land where
                // `split_query_primitives_with_neighbors` placed the
                // matching primitive in `covered_query`.
                let truth_alts_set: BTreeSet<&str> = truth.key.alt_allele.split(',').collect();
                let bases_for_residual = reference.as_bytes();
                // Neighbor floor mirrors `try_split_same_anchor_via_shift`:
                // the residual must NOT slide onto an anchor occupied by
                // another raw cluster query. chr21:21690513 chr21 case has
                // a filtered SNP `T→C` at 21690501; the CACAC residual
                // would otherwise slide there and clobber that record's
                // anchor — legacy stops one position above (21690502
                // A→ACACA) because xcmp doesn't permit two records to
                // share an anchor.
                let neighbor_floor_for_residual = cluster
                    .query
                    .iter()
                    .filter(|n| n.key.pos != query.key.pos)
                    .map(|n| n.key.pos)
                    .max();
                let pos_min_for_residual = match neighbor_floor_for_residual {
                    Some(n) => cluster.start.saturating_sub(SPLIT_LEFT_SHIFT_WINDOW).max(n),
                    None => cluster.start.saturating_sub(SPLIT_LEFT_SHIFT_WINDOW).max(1),
                };
                for unmatched_alt in query
                    .key
                    .alt_allele
                    .split(',')
                    .filter(|a| !truth_alts_set.contains(a))
                {
                    let (tpos, tref, talt) =
                        trim_variant(query.key.pos, &query.key.ref_allele, unmatched_alt);
                    // Shift insertion residuals through homopolymer /
                    // microsat reference patterns. Deletions and SNP
                    // primitives are returned unchanged by
                    // `compute_shift_target` because `partial_credit::
                    // left_shift` slides only when the deleted segment's
                    // last base matches the alt's last base — which is
                    // why the existing `apply_slide` gate in
                    // `split_query_primitives_with_neighbors` is
                    // restricted to deletion primitives.
                    let (rpos, rref, ralt) = if talt.len() > tref.len() && tref.len() == 1 {
                        compute_shift_target(
                            tpos,
                            &tref,
                            &talt,
                            bases_for_residual,
                            pos_min_for_residual,
                        )
                    } else {
                        (tpos, tref, talt)
                    };
                    residual_queries.push(Variant {
                        key: VariantKey {
                            chrom: query.key.chrom.clone(),
                            pos: rpos,
                            ref_allele: rref,
                            alt_allele: ralt,
                        },
                        qual: query.qual.clone(),
                        filter: query.filter.clone(),
                        gt: "0/1".to_string(),
                    });
                }
            } else if truth.key.alt_allele == query.key.alt_allele
                && equivalent_gt(&truth.gt, &query.gt)
            {
                let query_for_row = if is_multi_allelic(query) {
                    // For a phased query (GT has '|'), legacy mirrors truth's
                    // haplotype assignment unphased: truth 1|2 → query 1/2.
                    // canonical_hetalt_gt (alphabetical "later/earlier") is wrong
                    // for phased cases like CACACACAT,CAT where it gives 2|1.
                    // For an unphased query (GT has '/'), legacy uses canonical
                    // alphabetical ordering (e.g. C,G unphased 1/2 → 2/1).
                    let gt = if query.gt.contains('|') {
                        truth.gt.replace('|', "/")
                    } else {
                        canonical_hetalt_gt(&truth.key.alt_allele, query)
                    };
                    Variant {
                        key: query.key.clone(),
                        qual: query.qual.clone(),
                        filter: query.filter.clone(),
                        gt,
                    }
                } else {
                    query.clone()
                };
                let regions = region_state.row_tags(Some(truth), Some(query));
                if combined_unk {
                    rows.push(unk_combined_row(
                        truth,
                        &query_for_row,
                        reference,
                        cluster.start,
                        &regions,
                    ));
                } else {
                    rows.push(tp_combined_row(
                        truth,
                        &query_for_row,
                        reference,
                        cluster.start,
                        &regions,
                    ));
                }
            } else {
                // Multi-allelic pair whose ALT columns are reordered
                // (e.g. truth `C→CA,CAA` vs query `C→CAA,CA`). Legacy's
                // xcmp output keeps TRUTH's ALT ordering in the alt
                // column but re-indexes the QUERY's GT against the
                // *alphabetically-sorted* version of query's own ALTs —
                // not against truth's order. That's why query `CAA,CA
                // GT=1/2` emits `2/1` (CA sorts first, so CAA→2, CA→1)
                // while `AT,ATT GT=1/2` emits `1/2` verbatim (AT,ATT is
                // already alpha-sorted). Truth GT is always printed
                // verbatim under truth's displayed alt ordering.
                let canonicalized = Variant {
                    key: VariantKey {
                        chrom: query.key.chrom.clone(),
                        pos: query.key.pos,
                        ref_allele: query.key.ref_allele.clone(),
                        alt_allele: truth.key.alt_allele.clone(),
                    },
                    qual: query.qual.clone(),
                    filter: query.filter.clone(),
                    gt: canonical_hetalt_gt(&truth.key.alt_allele, query),
                };
                let regions = region_state.row_tags(Some(truth), Some(query));
                if combined_unk {
                    rows.push(unk_combined_row(
                        truth,
                        &canonicalized,
                        reference,
                        cluster.start,
                        &regions,
                    ));
                } else {
                    rows.push(tp_combined_row(
                        truth,
                        &canonicalized,
                        reference,
                        cluster.start,
                        &regions,
                    ));
                }
            }
            matched_truth.insert(truth_index);
            matched_query.insert(query_index);
        }
    }

    *truth_remaining = truth_remaining
        .iter()
        .enumerate()
        .filter(|(index, _)| !matched_truth.contains(index))
        .map(|(_, variant)| variant.clone())
        .collect();
    let mut new_query_remaining: Vec<Variant> = query_remaining
        .iter()
        .enumerate()
        .filter(|(index, _)| !matched_query.contains(index))
        .map(|(_, variant)| variant.clone())
        .collect();
    new_query_remaining.append(&mut residual_queries);
    *query_remaining = new_query_remaining;
}

pub(super) fn mark_cluster_match(
    cluster: &Cluster,
    // Full (pre-exact-match) cluster — legacy's per-record QQ comes from
    // the block's ORIGINAL representative query qual, not the post-exact-
    // match remainder. Without this the truth-only TP row at a multi-
    // allelic split locus picks up the wrong QQ when an earlier SNP pair
    // in the same block already consumed its byte-equal counterpart.
    full_cluster: &Cluster,
    reference: &str,
    region_state: &RegionState,
    outputs: ComparisonOutputs<'_>,
) {
    let ComparisonOutputs {
        counts,
        subtype_counts,
        rows,
    } = outputs;
    // Per-record CONF gate: legacy's `XCmpQuantify::countVariants`
    // unconditionally rewrites any output record's BD to UNK when the
    // record's Regions tag set lacks "CONF":
    //   if (count_unk && !tag_contains("CONF")) type = "UNK";
    // This applies even inside a haplotype-matched cluster, so a query
    // SNP at a TS_boundary-but-not-CONF position is emitted as UNK with
    // BK=. while the in-CONF indels around it stay TP/gm. Without this,
    // rust over-credits outside-CONF query records as TP and inflates
    // shared_qq calculation by including their qual values.
    for truth in &cluster.truth {
        if region_state.truth_is_conf(truth) {
            add_variant_stats(
                &mut counts
                    .entry(truth.primary_type().to_string())
                    .or_default()
                    .truth_tp,
                truth,
            );
            add_variant_stats_subtype(subtype_counts, truth.primary_type(), truth, |stats| {
                &mut stats.truth_tp
            });
        }
    }
    for query in &cluster.query {
        if region_state.query_is_conf(query) {
            add_variant_stats(
                &mut counts
                    .entry(query.primary_type().to_string())
                    .or_default()
                    .query_tp,
                query,
            );
            add_variant_stats_subtype(subtype_counts, query.primary_type(), query, |stats| {
                &mut stats.query_tp
            });
        } else {
            add_variant_stats(
                &mut counts
                    .entry(query.primary_type().to_string())
                    .or_default()
                    .query_unk,
                query,
            );
            add_variant_stats_subtype(subtype_counts, query.primary_type(), query, |stats| {
                &mut stats.query_unk
            });
        }
    }

    if cluster.truth.len() == 1
        && cluster.query.len() == 1
        && query_matches_truth_key(&cluster.query[0], &cluster.truth[0])
        && !query_primitive_splits(
            &cluster.truth[0],
            reference,
            cluster.start,
            &cluster.truth,
            &cluster.query,
        )
        && !query_primitive_splits(
            &cluster.query[0],
            reference,
            cluster.start,
            &cluster.truth,
            &cluster.query,
        )
    {
        let truth = &cluster.truth[0];
        let query = &cluster.query[0];
        let regions = region_state.row_tags(Some(truth), Some(query));
        // Combined-row emit also obeys the per-record CONF gate: a single
        // truth+query exact match outside the confident region collapses
        // to UNK/lm on both sides rather than TP/gm. BK=lm is the legacy
        // verdict whenever the same-locus counterpart shares an alt
        // allele (which is by construction here — `query_matches_truth_key`
        // demands ref+alt equality).
        if !region_state.truth_is_conf(truth) && !region_state.query_is_conf(query) {
            rows.push(unk_combined_row(
                truth,
                query,
                reference,
                cluster.start,
                &regions,
            ));
        } else {
            rows.push(tp_combined_row(
                truth,
                query,
                reference,
                cluster.start,
                &regions,
            ));
        }
        return;
    }

    // Haplotype-matched clusters: truth and query agree at the haplotype
    // level even when their decomposed records differ (e.g. a 2-bp insert
    // as TA→TAA,T truth vs T→TA query). Legacy's xcmp stamps each output
    // record's IQQ with the block-level minimum query QUAL — i.e. the
    // most-conservative call quality in the superlocus — and
    // XCmpQuantify propagates that IQQ into every sample column's FORMAT/QQ.
    // Picking the minimum (rather than the first) keeps the truth-side QQ
    // on truth-only TP rows byte-equal to legacy when the cluster carries
    // multiple query variants at different quals (e.g. chr21:15246157
    // TA→TAA,T → QQ=174.59, not the earlier SNP's 817.09).
    //
    // Restrict shared_qq to IN-CONF query records: outside-CONF queries
    // are reclassified to UNK below and never contribute their qual to
    // the truth-side TP rows. Without this filter, a chr21:18827409
    // cluster (truth indels in CONF, query SNPs outside CONF) picks up
    // the SNP min-qual instead of the matched indel's qual.
    let shared_qq = full_cluster
        .query
        .iter()
        .filter(|q| region_state.query_is_conf(q))
        .filter_map(|q| {
            q.qual
                .parse::<f64>()
                .ok()
                .filter(|v| *v > 0.0)
                .map(|v| (v, q.qual.as_str()))
        })
        .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(_, s)| s);
    // Truth-side TP split rows inherit the cluster's aggregated query
    // filter: legacy's bcftools-merged input stamps these rows with the
    // union of non-PASS filter tokens from every query record in the
    // cluster, not with the truth's (always PASS/`.`) filter. Compute
    // once per cluster against `full_cluster` so exact-matched queries
    // that were already consumed still contribute.
    let truth_cluster_filter = cluster_query_filter(full_cluster);
    for truth in &cluster.truth {
        if region_state.truth_is_conf(truth) {
            rows.push(tp_single_side_row(
                truth,
                reference,
                cluster.start,
                &region_state.row_tags(Some(truth), None),
                Side::Truth,
                shared_qq,
                &truth_cluster_filter,
            ));
        } else {
            // Truth outside CONF in a gm-matched cluster: emit UNK,
            // matching legacy's unconditional `count_unk && !CONF → UNK`
            // rewrite. BK=lm if a same-locus query counterpart exists,
            // else BK=. via bk_for_row. mark_cluster_match is reached
            // only when hapcmp returned match, so hap_mismatch=false.
            let bk = bk_for_row(truth, &full_cluster.query, false);
            rows.push(unk_truth_row(
                truth,
                reference,
                cluster.start,
                &region_state.row_tags(Some(truth), None),
                bk,
            ));
        }
    }
    for query in &cluster.query {
        // pre.py decomposes a complex indel into adjacent insertion/deletion
        // records with the same QUAL. In a hap:match block legacy still
        // stamps those unmatched sibling primitives BK=lm. Keep this on the
        // hap-match path only: `--unhappy` skips hapcmp and leaves the same
        // query-only primitives at BK=`.`.
        let split_sibling = has_nonconf_split_sibling(query, full_cluster, region_state);
        // Multi-allelic query records emit one VCF row per per-allele
        // primitive, matching legacy's per-primitive output grain. CONF
        // status is evaluated per primitive: a multi-allelic record's
        // deletion allele can sit fully inside CONF while its insertion
        // allele straddles a CONF edge, in which case the deletion gets
        // TP/gm but the insertion gets UNK/. (legacy's
        // `count_unk && !CONF → UNK` rewrite applies per output row).
        let primitives = split_query_primitives_with_neighbors(
            query,
            reference,
            cluster.start,
            &full_cluster.query,
            &full_cluster.truth,
        );
        let any_primitive_in_conf = primitives
            .iter()
            .any(|primitive| region_state.query_is_conf(primitive));
        let fanned_out = primitives.len() > 1;
        for primitive in primitives {
            let primitive_in_conf = region_state.query_is_conf(&primitive);
            let regions = region_state.row_tags(None, Some(&primitive));
            if primitive_in_conf {
                rows.push(tp_single_side_row(
                    &primitive,
                    reference,
                    cluster.start,
                    &regions,
                    Side::Query,
                    None,
                    ".",
                ));
            } else {
                let fallback = bk_for_row(&primitive, &full_cluster.truth, false);
                let bk = if split_sibling {
                    "lm"
                } else {
                    matched_query_unk_bk(fanned_out, any_primitive_in_conf, fallback)
                };
                rows.push(fp_like_row(
                    &primitive,
                    reference,
                    cluster.start,
                    &regions,
                    "UNK",
                    None,
                    bk,
                ));
            }
        }
    }
}

pub(super) fn has_nonconf_split_sibling(
    query: &Variant,
    cluster: &Cluster,
    region_state: &RegionState,
) -> bool {
    if query.primary_type() != "INDEL" || region_state.query_is_conf(query) {
        return false;
    }
    let query_is_insertion = query.key.alt_allele.len() > query.key.ref_allele.len();
    let query_is_deletion = query.key.ref_allele.len() > query.key.alt_allele.len();
    if !query_is_insertion && !query_is_deletion {
        return false;
    }

    cluster.query.iter().any(|sibling| {
        if sibling.key == query.key
            || sibling.qual != query.qual
            || sibling.primary_type() != "INDEL"
            || region_state.query_is_conf(sibling)
        {
            return false;
        }
        let sibling_is_insertion = sibling.key.alt_allele.len() > sibling.key.ref_allele.len();
        let sibling_is_deletion = sibling.key.ref_allele.len() > sibling.key.alt_allele.len();
        let complementary = (query_is_insertion && sibling_is_deletion)
            || (query_is_deletion && sibling_is_insertion);
        let split_span = query
            .key
            .ref_allele
            .len()
            .max(sibling.key.ref_allele.len())
            .max(1);
        complementary && query.key.pos.abs_diff(sibling.key.pos) < split_span
    })
}

pub(super) fn matched_query_unk_bk(
    fanned_out: bool,
    any_primitive_in_conf: bool,
    fallback: &'static str,
) -> &'static str {
    if fanned_out && !any_primitive_in_conf {
        "lm"
    } else {
        fallback
    }
}

pub(super) fn mark_cluster_mismatch(
    cluster: &Cluster,
    // Full (pre-exact-match) cluster — counterpart pool scanned by
    // `almismatch_same_locus`. Must include variants already paired as
    // TP so that same-locus allele-disjoint pairings survive the exact-
    // match sweep; legacy's simple-compare sees every call in the block
    // regardless of pairing outcome.
    full_cluster: &Cluster,
    // `ctype == "hap:mismatch"` verdict from `compare_cluster`: true iff
    // the block-level haplotype comparator ran (n_nonsnp gate fired on
    // the cluster) AND the two haplotype signatures disagreed. Feeds
    // Path B of `bk_for_row`.
    hap_mismatch: bool,
    reference: &str,
    region_state: &RegionState,
    outputs: ComparisonOutputs<'_>,
) {
    let ComparisonOutputs {
        counts,
        subtype_counts,
        rows,
    } = outputs;
    // Same-locus FN+FP pairs: when a truth record and a query record
    // share chrom/pos/ref/alt set but disagree on genotype (e.g. truth
    // 0|1 het vs query 1/1 homalt), legacy emits a single combined row
    // with BD=FN on the truth sample, BD=FP on the query sample, and
    // BK=am (allele-match, genotype-mismatch). Rust previously split
    // these into two rows, inflating the VCF body line count by one
    // per pair and losing the BK=am annotation.
    let mut paired_truth: BTreeSet<usize> = BTreeSet::new();
    let mut paired_query: BTreeSet<usize> = BTreeSet::new();
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for (ti, truth) in cluster.truth.iter().enumerate() {
        if paired_truth.contains(&ti) {
            continue;
        }
        for (qi, query) in cluster.query.iter().enumerate() {
            if paired_query.contains(&qi) {
                continue;
            }
            if query_matches_truth_allele_set(query, truth)
                && selected_alt_sequences(truth) != selected_alt_sequences(query)
            {
                pairs.push((ti, qi));
                paired_truth.insert(ti);
                paired_query.insert(qi);
                break;
            }
        }
    }
    for (ti, qi) in pairs {
        let truth = &cluster.truth[ti];
        let query = &cluster.query[qi];
        let truth_in_conf = region_state.truth_is_conf(truth);
        let query_in_conf = region_state.query_is_conf(query);
        // Only count as truth_fn when truth is inside CONF — outside CONF
        // the row emits BD=UNK so it must not contribute to the FN tally.
        if truth_in_conf {
            add_variant_stats(
                &mut counts
                    .entry(truth.primary_type().to_string())
                    .or_default()
                    .truth_fn,
                truth,
            );
            add_variant_stats_subtype(subtype_counts, truth.primary_type(), truth, |stats| {
                &mut stats.truth_fn
            });
        }
        let query_bd: &'static str = if query_in_conf { "FP" } else { "UNK" };
        // Per legacy's `count_unk && !CONF → UNK` rewrite, a paired
        // truth+query in a non-CONF region downgrades both samples to
        // BD=UNK. Truth-side BD respects this even when query is FP/CONF:
        // a hetalt truth in TS_boundary keeps `UNK:am` while the matched
        // query inside CONF emits `FP:am`. Without this, chr21:38484260
        // emits FN where legacy says UNK on truth.
        let truth_bd: &'static str = if truth_in_conf { "FN" } else { "UNK" };
        if query_in_conf {
            add_variant_stats(
                &mut counts
                    .entry(query.primary_type().to_string())
                    .or_default()
                    .query_fp,
                query,
            );
            add_variant_stats_subtype(subtype_counts, query.primary_type(), query, |stats| {
                &mut stats.query_fp
            });
        } else {
            add_variant_stats(
                &mut counts
                    .entry(query.primary_type().to_string())
                    .or_default()
                    .query_unk,
                query,
            );
            add_variant_stats_subtype(subtype_counts, query.primary_type(), query, |stats| {
                &mut stats.query_unk
            });
        }
        let bk = compute_paired_bk(truth, query);
        // Canonicalise unphased hetalt query GT to legacy's
        // `<later>/<earlier>` ordering when the row is hetalt and the
        // declared alts have a non-canonical order. SNP cases like
        // chr21:9922359 (`T→A,C` query GT `1/2`) emit `2/1` per
        // VariantLocationAggregator's MAX_GT=2 rule. Phased / hom /
        // het-with-ref pass through unchanged.
        let query_gt = if is_distinct_hetalt(&query.gt) {
            canonical_hetalt_gt(&truth.key.alt_allele, query)
        } else {
            remap_query_gt_subset(truth, query)
        };
        let query_for_row = Variant {
            key: truth.key.clone(),
            qual: query.qual.clone(),
            filter: query.filter.clone(),
            gt: query_gt,
        };
        rows.push(fn_fp_combined_row(
            truth,
            &query_for_row,
            reference,
            cluster.start,
            &region_state.row_tags(Some(truth), Some(query)),
            truth_bd,
            query_bd,
            bk,
        ));
    }
    for (ti, truth) in cluster.truth.iter().enumerate() {
        if paired_truth.contains(&ti) {
            continue;
        }
        let emit_unk = !region_state.truth_is_conf(truth);
        if emit_unk {
            // Truth outside the confident region is always emitted BD=UNK
            // in germline mode, matching legacy's
            // `XCmpQuantify::countVariants` unconditional rewrite:
            //   if (count_unk && !tag_contains("CONF")) type = "UNK";
            // The prior rust gate also required `!full_cluster.query.is_empty()`
            // but that is stricter than legacy — a truth-only cluster in a
            // TS_boundary-but-not-CONF region still becomes UNK. BK=lm
            // fires when a nearby query variant anchors the local-match
            // heuristic via `bk_for_row`.
            let bk = bk_for_row(truth, &full_cluster.query, hap_mismatch);
            rows.push(unk_truth_row(
                truth,
                reference,
                cluster.start,
                &region_state.row_tags(Some(truth), None),
                bk,
            ));
        } else {
            add_variant_stats(
                &mut counts
                    .entry(truth.primary_type().to_string())
                    .or_default()
                    .truth_fn,
                truth,
            );
            add_variant_stats_subtype(subtype_counts, truth.primary_type(), truth, |stats| {
                &mut stats.truth_fn
            });
            let bk = bk_for_row(truth, &full_cluster.query, hap_mismatch);
            rows.push(fn_row(
                truth,
                reference,
                cluster.start,
                &region_state.row_tags(Some(truth), None),
                bk,
            ));
        }
    }
    for (qi, query) in cluster.query.iter().enumerate() {
        if paired_query.contains(&qi) {
            continue;
        }
        if let Some(split_rows) =
            split_query_mismatch_rows(query, reference, region_state, cluster.start)
        {
            for (pseudo, bd, regions) in split_rows {
                match bd {
                    "FP" => {
                        add_variant_stats(
                            &mut counts
                                .entry(pseudo.primary_type().to_string())
                                .or_default()
                                .query_fp,
                            &pseudo,
                        );
                        add_variant_stats_subtype(
                            subtype_counts,
                            pseudo.primary_type(),
                            &pseudo,
                            |stats| &mut stats.query_fp,
                        );
                    }
                    "UNK" => {
                        add_variant_stats(
                            &mut counts
                                .entry(pseudo.primary_type().to_string())
                                .or_default()
                                .query_unk,
                            &pseudo,
                        );
                        add_variant_stats_subtype(
                            subtype_counts,
                            pseudo.primary_type(),
                            &pseudo,
                            |stats| &mut stats.query_unk,
                        );
                    }
                    _ => {}
                }
                add_variant_stats(
                    &mut counts
                        .entry(pseudo.primary_type().to_string())
                        .or_default()
                        .query_total,
                    &pseudo,
                );
                add_variant_stats_subtype(
                    subtype_counts,
                    pseudo.primary_type(),
                    &pseudo,
                    |stats| &mut stats.query_total,
                );
                let pseudo_bk = bk_for_row(&pseudo, &full_cluster.truth, hap_mismatch);
                let pseudo_fp_class = if bd == "FP" {
                    fp_class_from_bk(pseudo_bk)
                } else {
                    None
                };
                rows.push(fp_like_row(
                    &pseudo,
                    reference,
                    cluster.start,
                    &regions,
                    bd,
                    pseudo_fp_class,
                    pseudo_bk,
                ));
            }
            continue;
        }
        // Even same-type multi-allelic queries fan out one primitive row
        // per active allele in legacy's output — matches the decomposed
        // per-primitive representation xcmp emits after classification.
        // Each primitive picks up its OWN Regions classification: a
        // multi-allelic record whose deletion allele sits inside CONF
        // but whose insertion allele straddles a CONF edge emits one
        // row with FP/CONF and a second with UNK/no-CONF. The parent's
        // any-allele coverage decision is too coarse for that.
        let primitives = split_query_primitives_with_neighbors(
            query,
            reference,
            cluster.start,
            &full_cluster.query,
            &full_cluster.truth,
        );
        let any_primitive_in_conf = primitives
            .iter()
            .any(|primitive| region_state.query_is_conf(primitive));
        let fanned_out = primitives.len() > 1;
        for primitive in primitives {
            let is_conf = region_state.query_is_conf(&primitive);
            let bd: &'static str = if is_conf { "FP" } else { "UNK" };
            let regions = region_state.row_tags(None, Some(&primitive));
            if is_conf {
                add_variant_stats(
                    &mut counts
                        .entry(primitive.primary_type().to_string())
                        .or_default()
                        .query_fp,
                    &primitive,
                );
                add_variant_stats_subtype(
                    subtype_counts,
                    primitive.primary_type(),
                    &primitive,
                    |stats| &mut stats.query_fp,
                );
            } else {
                add_variant_stats(
                    &mut counts
                        .entry(primitive.primary_type().to_string())
                        .or_default()
                        .query_unk,
                    &primitive,
                );
                add_variant_stats_subtype(
                    subtype_counts,
                    primitive.primary_type(),
                    &primitive,
                    |stats| &mut stats.query_unk,
                );
            }
            let fallback = bk_for_row(&primitive, &full_cluster.truth, hap_mismatch);
            let primitive_bk = matched_query_unk_bk(fanned_out, any_primitive_in_conf, fallback);
            let fp_class = if bd == "FP" {
                fp_class_from_bk(primitive_bk)
            } else {
                None
            };
            rows.push(fp_like_row(
                &primitive,
                reference,
                cluster.start,
                &regions,
                bd,
                fp_class,
                primitive_bk,
            ));
        }
    }
}

/// Map the per-row `BK` (block kind) tag to the FP classification used by
/// the summary / extended / ROC `FP.gt` / `FP.al` columns. Mirrors legacy
/// hap.py's quantify aggregation: FP rows tagged `BK=am` (allele match,
/// genotype mismatch) bump `FP.gt`; FP rows tagged `BK=lm` (locus match,
/// allele mismatch) bump `FP.al`; everything else (novel FPs with `BK=.`)
/// stays unclassified and is omitted from both columns.
pub(super) fn fp_class_from_bk(bk: &str) -> Option<&'static str> {
    match bk {
        "am" => Some("gt"),
        "lm" => Some("al"),
        _ => None,
    }
}

/// Two VCF keys describe the same locus when they agree on chrom/pos/ref
/// and their alt lists are equal as sets (order-independent). Legacy's
/// xcmp exact-match pairing treats `A ATT,AT` truth and `A AT,ATT` query
/// as the same record with reshuffled alt indices; rust's byte-equal
/// check would miss this and fall through to the haplotype-match path
/// (which emits two rows instead of the combined row legacy expects).
pub(super) fn query_matches_truth_key(query: &Variant, truth: &Variant) -> bool {
    // Byte-equal ALT column, not BTreeSet-of-alts. Legacy xcmp's
    // `simpleCompare` pairs records on raw (chrom, pos, ref, alt)
    // string equality, so truth `GAA→GA,G` (GT 2|1) does NOT pair with
    // query `GAA→G,GA` (GT 1/2) — the ALT column-order flip is a
    // genuine non-match at the simple-compare layer even though both
    // encode the same unordered allele set. Under set-semantic rust
    // used to over-pair them, bypassing the block-level hap engine
    // that would have emitted FN/FP once a neighbouring variant
    // broke the hap signature (see chr21:18743964-18743965).
    // GT multi-set equivalence is checked separately by
    // `equivalent_gt` at the call site.
    query.key.chrom == truth.key.chrom
        && query.key.pos == truth.key.pos
        && query.key.ref_allele == truth.key.ref_allele
        && query.key.alt_allele == truth.key.alt_allele
}

/// Same physical VCF locus with the same declared alleles, allowing ALT
/// columns to use different index orders. Legacy's shared allele table uses
/// this equivalence when it emits a combined FN/FP genotype-mismatch row;
/// the query GT is remapped into truth's ALT order before serialization.
pub(super) fn query_matches_truth_allele_set(query: &Variant, truth: &Variant) -> bool {
    query.key.chrom == truth.key.chrom
        && query.key.pos == truth.key.pos
        && query.key.ref_allele == truth.key.ref_allele
        && query.key.alt_allele.split(',').collect::<BTreeSet<_>>()
            == truth.key.alt_allele.split(',').collect::<BTreeSet<_>>()
}

pub(super) fn is_distinct_hetalt(gt: &str) -> bool {
    let alleles = parse_gt_alleles(gt);
    alleles.len() == 2 && alleles[0] > 0 && alleles[1] > 0 && alleles[0] != alleles[1]
}

pub(super) fn split_query_mismatch_rows(
    query: &Variant,
    _reference: &str,
    region_state: &RegionState,
    _block_start: usize,
) -> Option<Vec<(Variant, &'static str, String)>> {
    if !query.key.alt_allele.contains(',') {
        return None;
    }
    let alleles = parse_gt_alleles(&query.gt);
    let used: BTreeSet<usize> = alleles.into_iter().filter(|allele| *allele > 0).collect();
    if used.is_empty() {
        return None;
    }

    let alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    let mut components: Vec<(Variant, bool)> = Vec::new();
    let mut types = BTreeSet::new();
    for allele in used {
        let alt = *alts.get(allele - 1)?;
        let pseudo = Variant {
            key: VariantKey {
                chrom: query.key.chrom.clone(),
                pos: query.key.pos,
                ref_allele: query.key.ref_allele.clone(),
                alt_allele: alt.to_string(),
            },
            qual: query.qual.clone(),
            filter: query.filter.clone(),
            gt: "0/1".to_string(),
        };
        types.insert(pseudo.primary_type().to_string());
        let covered = component_is_conf(&pseudo, region_state);
        components.push((pseudo, covered));
    }

    if types.len() <= 1 {
        return None;
    }

    let any_conf = components.iter().any(|(_, covered)| *covered);
    let any_nonconf = components.iter().any(|(_, covered)| !*covered);
    let mut rows = Vec::new();
    for (pseudo, covered) in components {
        let mut tags = Vec::new();
        if covered {
            tags.push("CONF");
        }
        if any_conf {
            if any_nonconf {
                tags.push("TS_boundary");
            } else {
                tags.push("TS_contained");
            }
        } else if region_state.any_conf {
            tags.push("TS_boundary");
        }
        let regions = if tags.is_empty() {
            String::new()
        } else {
            format!(";Regions={}", tags.join(","))
        };
        let bd = if covered { "FP" } else { "UNK" };
        rows.push((pseudo, bd, regions));
    }
    Some(rows)
}

pub(super) fn component_is_conf(variant: &Variant, region_state: &RegionState) -> bool {
    if !region_state.conf_enabled {
        return true;
    }
    if variant.primary_type() == "SNP" {
        return region_state
            .conf_intervals
            .iter()
            .any(|interval| interval.matches(&variant.key.chrom, variant.key.pos));
    }
    if variant.key.alt_allele.starts_with(&variant.key.ref_allele)
        && variant.key.alt_allele.len() > variant.key.ref_allele.len()
    {
        let anchor = variant.key.pos;
        return region_state
            .conf_intervals
            .iter()
            .any(|interval| interval.matches(&variant.key.chrom, anchor))
            && region_state
                .conf_intervals
                .iter()
                .any(|interval| interval.matches(&variant.key.chrom, anchor + 1));
    }
    region_state
        .conf_intervals
        .iter()
        .any(|interval| interval.matches(&variant.key.chrom, variant.key.pos))
}

/// Legacy's `QuantifyRegions::annotate` computes an "effective" reference
/// span per record by per-allele trimming (trimRight then trimLeft with
/// refpadding=false), then:
///   - substitutions / deletions: [start, end] from the post-trim remnant
///   - pure insertions (all trimmed alts empty ref): [anchor, anchor+1]
///     — the anchor base and the one after it bracket the insertion point
///
/// The record then requires `!is_pure_insertion || fully_covered` to land
/// in a region. Port that rule here so insertion / SNP records get the
/// same CONF / TS_boundary classification legacy emits.
///
/// Returns `(refstart_1b, refend_1b, is_pure_insertion)` in 1-based
/// inclusive VCF coordinates (matching `variant.key.pos`). `None` when
/// the record has no nucleotide alts (all symbolic / missing).
pub(super) fn effective_refrange(variant: &Variant) -> Option<(usize, usize, bool)> {
    let ref_bytes = variant.key.ref_allele.as_bytes();
    let pos_0b = variant.key.pos.saturating_sub(1) as i64;
    let mut updated_start = i64::MAX;
    let mut updated_end = i64::MIN;
    let mut is_pure_insertion = false;
    let mut has_nuc = false;

    for alt in variant.key.alt_allele.split(',') {
        if alt.is_empty() || alt == "." {
            continue;
        }
        if alt.starts_with('<') {
            continue;
        }
        // trimRight (refpadding=false): strip common suffix from ref/alt.
        let alt_bytes = alt.as_bytes();
        let mut reflen = ref_bytes.len();
        let mut altlen = alt_bytes.len();
        while reflen > 0 && altlen > 0 && ref_bytes[reflen - 1] == alt_bytes[altlen - 1] {
            reflen -= 1;
            altlen -= 1;
        }
        // trimLeft (refpadding=false): strip common prefix.
        let mut rel_start = 0usize;
        while rel_start < reflen
            && rel_start < altlen
            && ref_bytes[rel_start] == alt_bytes[rel_start]
        {
            rel_start += 1;
        }
        let al_start = pos_0b + rel_start as i64;
        let al_end = pos_0b + reflen as i64 - 1;

        if !has_nuc {
            is_pure_insertion = true;
        }
        has_nuc = true;

        if al_end >= al_start {
            // subst/del — record still has ref bases after trim
            updated_start = updated_start.min(al_start);
            updated_end = updated_end.max(al_end);
            is_pure_insertion = false;
        } else {
            // insertion — effective range brackets the insertion point
            updated_start = updated_start.min(al_start - 1);
            updated_end = updated_end.max(al_start);
        }
    }

    if !has_nuc {
        return None;
    }
    // Convert 0-based inclusive back to 1-based inclusive (add 1 to both).
    Some((
        (updated_start + 1) as usize,
        (updated_end + 1) as usize,
        is_pure_insertion,
    ))
}

/// Replicate legacy `gvcf2bed`'s sequentially-merged padding output from
/// the truth VCF (src/c++/main/gvcf2bed.cpp). For each truth record we
/// compute the effective refstart/refend per legacy rules (the same as
/// `effective_refrange`, but done in 0-based coords). Records are merged
/// into contiguous bed intervals using legacy's idiosyncratic rule:
///   `current.end = max(refstart, current.end)` — **refstart, not refend**.
/// That intentional under-extension means an insertion folded into a
/// prior overlapping interval loses its following-base coverage, which
/// reproduces legacy's TS_boundary behavior for insertions straddling a
/// CONF edge. Replicate exactly — don't "fix" the legacy merge.
///
/// When `target_bed` is `Some`, mirror legacy's `bcf_sr_set_targets(..,
/// target.c_str(), 1, 0)` filter on the gvcf2bed call: the htslib synced
/// reader gates on the **start position only** (single-point overlap of
/// `pos_0b` against any half-open `[s, e)` target interval), not the
/// full ref span. Records whose start position falls outside the raw
/// CONF bed are dropped before emission. This filter is critical for
/// `Subset.IS_CONF.Size` parity (legacy sums per-file BED lengths
/// without cross-file merging).
///
/// The returned half-open `[start, end)` intervals are later merged with
/// the raw CONF bed (standard overlap/touching merge) before
/// `variant_is_conf` consumes them, but the **un-merged sum** is what
/// feeds `Subset.IS_CONF.Size`.
pub(super) fn gvcf2bed_padding(
    truth: &[Variant],
    target_bed: Option<&[Interval]>,
) -> Vec<Interval> {
    struct ActiveInterval {
        chrom: String,
        start: i64,
        end: i64,
    }

    // Build per-chrom sorted target ranges for fast point-in-set lookup.
    let target_index: Option<BTreeMap<String, Vec<(usize, usize)>>> = target_bed.map(|tb| {
        let mut by_chrom: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
        for iv in tb {
            by_chrom
                .entry(iv.chrom.clone())
                .or_default()
                .push((iv.start, iv.end));
        }
        for v in by_chrom.values_mut() {
            v.sort_unstable();
        }
        by_chrom
    });

    let pos_in_target = |chrom: &str, pos_0b: i64| -> bool {
        let Some(idx) = target_index.as_ref() else {
            return true;
        };
        let Some(ranges) = idx.get(chrom) else {
            return false;
        };
        if pos_0b < 0 {
            return false;
        }
        let p = pos_0b as usize;
        // Binary search for the first interval with start > p; check the
        // candidate predecessor for half-open containment [s, e).
        let idx = ranges.partition_point(|(s, _)| *s <= p);
        if idx == 0 {
            return false;
        }
        let (s, e) = ranges[idx - 1];
        s <= p && p < e
    };

    let mut sorted: Vec<&Variant> = truth.iter().collect();
    sorted.sort_by(|a, b| {
        a.key
            .chrom
            .cmp(&b.key.chrom)
            .then(a.key.pos.cmp(&b.key.pos))
    });

    let mut out = Vec::new();
    let mut active: Option<ActiveInterval> = None;

    for variant in sorted {
        let ref_bytes = variant.key.ref_allele.as_bytes();
        if ref_bytes.is_empty() {
            continue;
        }
        let pos_0b = variant.key.pos.saturating_sub(1) as i64;
        if !pos_in_target(&variant.key.chrom, pos_0b) {
            continue;
        }
        // Legacy default refstart/refend (from getLocation): the raw
        // 0-based-inclusive ref-allele span. These persist when the alt
        // loop produces no NUC alleles (all symbolic / `<NON_REF>` /
        // missing) — the record is still emitted with this raw range.
        let mut rec_start = pos_0b;
        let mut rec_end = pos_0b + ref_bytes.len() as i64 - 1;
        let ref_is_nuc = ref_bytes.iter().all(|b| {
            matches!(
                *b,
                b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
            )
        });

        if ref_is_nuc {
            let mut updated_start = i64::MAX;
            let mut updated_end = i64::MIN;
            let mut nuc_alleles = false;
            // Mirror legacy's behavior: a non-NUC alt **breaks** the
            // loop (it does not just skip — `break` in the C++ source).
            // Subsequent NUC alts after a symbolic one are ignored, so
            // updated_start/end track only alts processed up to the
            // first symbolic.
            for alt in variant.key.alt_allele.split(',') {
                if alt.is_empty() || alt == "." {
                    // MISSING — legacy treats as NUC with empty alt
                    // string, which after trim resolves to a full
                    // ref-deletion span [pos, pos+reflen-1].
                    nuc_alleles = true;
                    let al_start = pos_0b;
                    let al_end = pos_0b + ref_bytes.len() as i64 - 1;
                    if al_end >= al_start {
                        updated_start = updated_start.min(al_start);
                        updated_end = updated_end.max(al_end);
                    }
                    continue;
                }
                if alt.starts_with('<')
                    || alt.bytes().any(|b| {
                        !matches!(
                            b,
                            b'A' | b'C' | b'G' | b'T' | b'N' | b'a' | b'c' | b'g' | b't' | b'n'
                        )
                    })
                {
                    break;
                }
                let alt_bytes = alt.as_bytes();
                let mut reflen = ref_bytes.len();
                let mut altlen = alt_bytes.len();
                while reflen > 0 && altlen > 0 && ref_bytes[reflen - 1] == alt_bytes[altlen - 1] {
                    reflen -= 1;
                    altlen -= 1;
                }
                let mut rel_start = 0usize;
                while rel_start < reflen
                    && rel_start < altlen
                    && ref_bytes[rel_start] == alt_bytes[rel_start]
                {
                    rel_start += 1;
                }
                let al_start = pos_0b + rel_start as i64;
                let al_end = pos_0b + reflen as i64 - 1;
                nuc_alleles = true;
                if al_end >= al_start {
                    updated_start = updated_start.min(al_start);
                    updated_end = updated_end.max(al_end);
                } else {
                    updated_start = updated_start.min(al_start - 1);
                    updated_end = updated_end.max(al_start);
                }
            }
            if nuc_alleles {
                rec_start = updated_start;
                rec_end = updated_end;
            }
        }

        let flush = match active.as_ref() {
            Some(iv) => {
                !(iv.chrom == variant.key.chrom && rec_start <= iv.end && rec_end >= iv.start)
            }
            None => true,
        };

        if flush {
            if let Some(iv) = active.take() {
                out.push(Interval {
                    chrom: iv.chrom,
                    start: iv.start as usize,
                    end: (iv.end + 1) as usize,
                });
            }
            active = Some(ActiveInterval {
                chrom: variant.key.chrom.clone(),
                start: rec_start,
                end: rec_end,
            });
        } else if let Some(iv) = active.as_mut() {
            iv.start = iv.start.min(rec_start);
            // NOTE: legacy uses refstart (not refend) here — replicate.
            iv.end = iv.end.max(rec_start);
        }
    }
    if let Some(iv) = active.take() {
        out.push(Interval {
            chrom: iv.chrom,
            start: iv.start as usize,
            end: (iv.end + 1) as usize,
        });
    }
    out
}

/// Merge overlapping or touching half-open bed intervals. Matches
/// legacy `IntervalList::add` semantics: two intervals [a,b) and [c,d)
/// merge when `c <= b` (overlap or end-touch). This runs on the union
/// of the raw CONF bed and `gvcf2bed_padding` output — the result is
/// the effective CONF lane that legacy's `QuantifyRegions::annotate`
/// queries against.
pub(super) fn merge_bed_intervals(intervals: &[Interval]) -> Vec<Interval> {
    let mut by_chrom: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    for iv in intervals {
        by_chrom
            .entry(iv.chrom.clone())
            .or_default()
            .push((iv.start, iv.end));
    }
    let mut out = Vec::new();
    for (chrom, mut ranges) in by_chrom {
        ranges.sort_unstable();
        let mut current: Option<(usize, usize)> = None;
        for (start, end) in ranges {
            match current {
                Some((cs, ce)) if start <= ce => {
                    current = Some((cs, ce.max(end)));
                }
                Some((cs, ce)) => {
                    out.push(Interval {
                        chrom: chrom.clone(),
                        start: cs,
                        end: ce,
                    });
                    current = Some((start, end));
                }
                None => current = Some((start, end)),
            }
        }
        if let Some((cs, ce)) = current {
            out.push(Interval {
                chrom: chrom.clone(),
                start: cs,
                end: ce,
            });
        }
    }
    out
}

pub(super) fn variant_is_conf(
    variant: &Variant,
    reference: &str,
    cluster_start: usize,
    cluster_end: usize,
    conf_intervals: &[Interval],
) -> bool {
    if conf_intervals.is_empty() {
        return false;
    }
    // Legacy rule (QuantifyRegions::annotate):
    //   compute effective [refstart, refend] per gvcf2bed trim rules;
    //   add CONF if hasOverlap && (!is_pure_insertion || fully_covered).
    // `conf_intervals` at this point is the raw BED merged with the
    // gvcf2bed-style padding derived from truth, so that insertion-
    // bridging at CONF edges matches legacy.
    if let Some((refstart_1b, refend_1b, is_pure_insertion)) = effective_refrange(variant) {
        let has_overlap = (refstart_1b..=refend_1b).any(|pos| {
            conf_intervals
                .iter()
                .any(|interval| interval.matches(&variant.key.chrom, pos))
        });
        if !has_overlap {
            return false;
        }
        if !is_pure_insertion {
            return true;
        }
        // Pure insertion — require full coverage of the bracketed range.
        let fully = (refstart_1b..=refend_1b).all(|pos| {
            conf_intervals
                .iter()
                .any(|interval| interval.matches(&variant.key.chrom, pos))
        });
        if fully {
            return true;
        }
        return false;
    }
    // Symbolic-alt fallback: reuse the normalized-events path so <DEL>
    // etc still classify against the adjusted CONF lane.
    let alleles = parse_gt_alleles(&variant.gt);
    for allele in alleles {
        if allele == 0 {
            continue;
        }
        let Ok(events) =
            normalized_events_for_allele(variant, allele, reference, cluster_start, cluster_end)
        else {
            continue;
        };
        let mut covered = true;
        for event in events {
            match event {
                Event::Subst { pos, .. } => {
                    covered &= conf_intervals
                        .iter()
                        .any(|interval| interval.matches(&variant.key.chrom, pos));
                }
                Event::Delete { start, end } => {
                    covered &= (start..=end).all(|pos| {
                        conf_intervals
                            .iter()
                            .any(|interval| interval.matches(&variant.key.chrom, pos))
                    });
                }
                Event::Insert { anchor, .. } => {
                    covered &= conf_intervals
                        .iter()
                        .any(|interval| interval.matches(&variant.key.chrom, anchor))
                        && conf_intervals
                            .iter()
                            .any(|interval| interval.matches(&variant.key.chrom, anchor + 1));
                }
            }
            if !covered {
                break;
            }
        }
        if covered {
            return true;
        }
    }
    false
}
