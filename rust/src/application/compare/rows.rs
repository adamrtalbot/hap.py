//! Extracted cohesive responsibility from the command façade.

use super::{
    AnnotatedRow, Cluster, Side, Variant, comparison_info, fp_class_from_bk, genotype_label,
    query_type_rank,
};
use crate::adapters::vcf::VariantKey;
use crate::domain::{ComparisonRecord, RawVcfRecord};
use crate::engines::partial_credit;
use std::collections::BTreeSet;

fn comparison_record(
    variant: &Variant,
    reference: String,
    alternate: String,
    quality: String,
    filter: String,
    info: String,
    truth_sample: String,
    query_sample: String,
) -> ComparisonRecord {
    RawVcfRecord {
        chrom: variant.key.chrom.clone(),
        pos: variant.key.pos,
        id: ".".to_string(),
        ref_allele: reference,
        alt_allele: alternate,
        qual: quality,
        filter,
        info,
        format: Some("GT:BD:BK:BI:BVT:BLT:QQ".to_string()),
        samples: vec![truth_sample, query_sample],
    }
    .into()
}

pub(super) fn combined_record_qual<'a>(truth: &'a Variant, query: &'a Variant) -> &'a str {
    let numeric = |qual: &str| {
        qual.parse::<f32>()
            .ok()
            .filter(|value| !value.is_nan())
            .unwrap_or(0.0)
    };
    if numeric(&truth.qual) >= numeric(&query.qual) {
        &truth.qual
    } else {
        &query.qual
    }
}

pub(super) fn tp_combined_row(
    truth: &Variant,
    query: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
) -> AnnotatedRow {
    let info = comparison_info(truth, reference);
    AnnotatedRow {
        sort_key: (truth.key.chrom.clone(), truth.key.pos, 1, 0),
        query_pass: filter_is_pass(&query.filter),
        fp_class: None,
        xcmp_ctype: None,
        xcmp_hap_match: false,
        record: comparison_record(
            truth,
            display_ref(truth, reference),
            display_alt(truth),
            combined_record_qual(truth, query).to_string(),
            filter_for_output(&query.filter).to_string(),
            format!("BS={block_start}{regions}"),
            format!(
                "{}:TP:gm:{info}:{}:{}:{}",
                truth.gt,
                truth.primary_type(),
                genotype_label(truth),
                query.qual
            ),
            format!(
                "{}:TP:gm:{info}:{}:{}:{}",
                query.gt,
                truth.primary_type(),
                genotype_label(query),
                query.qual
            ),
        ),
    }
}

/// Combined UNK+UNK row for a same-key truth+query pair where the locus
/// falls outside the confident region. Legacy emits BD=UNK on both
/// samples with BK=lm — the locus-match heuristic fires because the two
/// sides share the variant exactly, and the unconditional non-CONF →
/// UNK rewrite trumps the gm verdict from xcmp. Truth-side QQ stays `.`
/// (truth's input qual is always "0") while query-side carries the
/// query's own qual.
pub(super) fn unk_combined_row(
    truth: &Variant,
    query: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
) -> AnnotatedRow {
    let info = comparison_info(truth, reference);
    AnnotatedRow {
        sort_key: (truth.key.chrom.clone(), truth.key.pos, 1, 0),
        query_pass: filter_is_pass(&query.filter),
        fp_class: None,
        xcmp_ctype: None,
        xcmp_hap_match: false,
        record: comparison_record(
            truth,
            display_ref(truth, reference),
            display_alt(truth),
            combined_record_qual(truth, query).to_string(),
            filter_for_output(&query.filter).to_string(),
            format!("BS={block_start}{regions}"),
            format!(
                "{}:UNK:lm:{info}:{}:{}:.",
                truth.gt,
                truth.primary_type(),
                genotype_label(truth)
            ),
            format!(
                "{}:UNK:lm:{info}:{}:{}:{}",
                query.gt,
                truth.primary_type(),
                genotype_label(query),
                query.qual
            ),
        ),
    }
}

/// Legacy treats an empty, `.`, or `PASS` filter column as the passing set.
/// Anything else (`OverlapConflict`, `SuspiciousHomAlt`, user-defined filter
/// names) is excluded from the PASS-tier rollups in summary/extended.
pub(super) fn filter_is_pass(filter: &str) -> bool {
    filter.is_empty() || filter == "." || filter == "PASS"
}

/// Render a query variant's FILTER column for the output VCF. Legacy
/// collapses empty / missing / PASS to `.` and passes everything else
/// through verbatim, preserving the original ;-separated filter tags
/// (e.g. `TruthSensitivityTranche99.00to99.90;LowQD`). Truth-only rows
/// always print `.` — the legacy emitter only forwards the query call's
/// filter, never the truth call's.
pub(super) fn filter_for_output(filter: &str) -> &str {
    if filter.is_empty() || filter == "PASS" {
        "."
    } else {
        filter
    }
}

/// Legacy emits a cluster-level filter column for truth-only TP rows (the
/// truth side of a split TP match where the matching query primitive is
/// emitted on a separate row). After bcftools merges truth+query into a
/// single VCF, each merged record's FILTER is the union of non-PASS
/// filter tokens from every query record in the cluster, with original
/// positional ordering preserved and de-duplicated. Empty → `.`.
pub(super) fn cluster_query_filter(cluster: &Cluster) -> String {
    let mut tokens: Vec<String> = Vec::new();
    for q in &cluster.query {
        for tok in q.filter.split(';') {
            let tok = tok.trim();
            if tok.is_empty() || tok == "PASS" || tok == "." {
                continue;
            }
            if !tokens.iter().any(|t| t == tok) {
                tokens.push(tok.to_string());
            }
        }
    }
    if tokens.is_empty() {
        ".".to_string()
    } else {
        // Legacy's bcftools-merged input deterministically orders the
        // FILTER column when stamping the cluster's union onto truth-side
        // TP rows: tokens land in byte-wise ascending order regardless of
        // the source order seen in any single query record. (Single-source
        // single-record rows preserve their own filter order; this path
        // only fires when the filter is the cluster aggregate.)
        tokens.sort();
        tokens.join(";")
    }
}

pub(super) fn tp_single_side_row(
    variant: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
    side: Side,
    shared_qq: Option<&str>,
    truth_filter: &str,
) -> AnnotatedRow {
    let info = comparison_info(variant, reference);
    match side {
        Side::Truth => AnnotatedRow {
            sort_key: (variant.key.chrom.clone(), variant.key.pos, 1, 0),
            // Truth-side TP split rows inherit the cluster's aggregated query
            // filter (see `cluster_query_filter` and call sites). When the
            // matching query is non-PASS, legacy demotes the truth-TP to FN
            // in the PASS-tier rollup; mirror that by deriving query_pass
            // from the same filter string that already lands in the FILTER
            // column.
            query_pass: filter_is_pass(truth_filter),
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
            record: comparison_record(
                variant,
                display_ref(variant, reference),
                display_alt(variant),
                variant.qual.clone(),
                truth_filter.to_string(),
                format!("BS={block_start}{regions}"),
                format!(
                    "{}:TP:gm:{info}:{}:{}:{}",
                    variant.gt,
                    variant.primary_type(),
                    genotype_label(variant),
                    shared_qq.unwrap_or(variant.qual.as_str())
                ),
                "./.:.:.:.:NOCALL:nocall:0".to_string(),
            ),
        },
        Side::Query => AnnotatedRow {
            sort_key: (
                variant.key.chrom.clone(),
                variant.key.pos,
                1,
                // +1 so query-only TP rows sort after truth-only TP rows at
                // the same position, matching legacy's truth-before-query
                // output order when both sides have a single-side TP row.
                query_type_rank(variant) + 1,
            ),
            query_pass: filter_is_pass(&variant.filter),
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
            record: comparison_record(
                variant,
                display_ref(variant, reference),
                display_alt(variant),
                variant.qual.clone(),
                filter_for_output(&variant.filter).to_string(),
                format!("BS={block_start}{regions}"),
                "./.:.:.:.:NOCALL:nocall:.".to_string(),
                format!(
                    "{}:TP:gm:{info}:{}:{}:{}",
                    variant.gt,
                    variant.primary_type(),
                    genotype_label(variant),
                    variant.qual
                ),
            ),
        },
    }
}

/// Combined FN+FP row for a truth+query pair that share chrom/pos/ref/alt
/// as a set but disagree on GT (e.g. truth 0|1 het vs query 1/1 homalt).
/// Legacy emits BD=FN on truth, BD=FP on query, BK=am on both. Record QUAL
/// is the maximum call QUAL, while query QQ retains the query's own score
/// for downstream ROC enumeration.
#[allow(clippy::too_many_arguments)] // Mirrors the two-sample legacy VCF row contract.
pub(super) fn fn_fp_combined_row(
    truth: &Variant,
    query: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
    truth_bd: &'static str,
    query_bd: &'static str,
    bk: &str,
) -> AnnotatedRow {
    // Per-sample BI: legacy emits each sample's `BI` based on the alleles
    // ITS own GT selects, not the merged record's overall subtype. The
    // chr21:9922359 fixture (truth `T→A,C 1|0` selects {A}=tv vs query
    // `T→A,C 1/2` selects {A,C}=ti,tv) demands different BI strings on
    // the two samples even on a combined fn_fp row.
    let truth_info = comparison_info(truth, reference);
    let query_info = comparison_info(query, reference);
    let fp_class = if query_bd == "FP" {
        fp_class_from_bk(bk)
    } else {
        None
    };
    AnnotatedRow {
        sort_key: (truth.key.chrom.clone(), truth.key.pos, 0, 0),
        // Combined rows inherit the query's PASS status for the ALL vs
        // PASS bifurcation — query filter drives whether this row
        // contributes to the PASS-tier summary, consistent with
        // legacy's derived-from-row contract.
        query_pass: filter_is_pass(&query.filter),
        // FP class follows the same rule as `fp_like_row`: BK=am→gt,
        // BK=lm→al, otherwise unclassified. Earlier this was hardcoded
        // to `Some("gt")` on the assumption that combined truth+query
        // rows are always genotype mismatches; that's true for the
        // common case but not for hetalt-vs-different-hetalt pairings
        // where compute_paired_bk returns `lm` (chr21:38861935 truth
        // `T→TAA,TA` 1|2 vs query `T→TAA` 0/1) — those should land in
        // FP.al, not FP.gt. Only emit the class when the query side is
        // actually FP (truth-side UNK pairs leave the row unclassified).
        fp_class,
        xcmp_ctype: None,
        xcmp_hap_match: false,
        record: comparison_record(
            truth,
            display_ref(truth, reference),
            display_alt(truth),
            combined_record_qual(truth, query).to_string(),
            filter_for_output(&query.filter).to_string(),
            format!("BS={block_start}{regions}"),
            format!(
                "{}:{truth_bd}:{bk}:{truth_info}:{}:{}:.",
                truth.gt,
                truth.primary_type(),
                genotype_label(truth)
            ),
            format!(
                "{}:{query_bd}:{bk}:{query_info}:{}:{}:{}",
                query.gt,
                truth.primary_type(),
                genotype_label(query),
                query.qual
            ),
        ),
    }
}

pub(super) fn fn_row(
    truth: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
    bk: &'static str,
) -> AnnotatedRow {
    let info = comparison_info(truth, reference);
    AnnotatedRow {
        sort_key: (truth.key.chrom.clone(), truth.key.pos, 0, 0),
        // FN rows count in both ALL and PASS tiers (legacy truth is already
        // pass-filtered upstream; the query-filter flag is not a gating
        // condition on truth-side FN counts).
        query_pass: true,
        fp_class: None,
        xcmp_ctype: None,
        xcmp_hap_match: false,
        record: comparison_record(
            truth,
            display_ref(truth, reference),
            display_alt(truth),
            truth.qual.clone(),
            ".".to_string(),
            format!("BS={block_start}{regions}"),
            format!(
                "{}:FN:{bk}:{info}:{}:{}:.",
                truth.gt,
                truth.primary_type(),
                genotype_label(truth)
            ),
            "./.:.:.:.:NOCALL:nocall:0".to_string(),
        ),
    }
}

/// Truth variant outside the confident region emits BD=UNK (not FN) in
/// legacy; legacy additionally sets BK=lm when a cluster-adjacent query
/// variant shares at least one alternate allele at byte level. Mirrors
/// the legacy classifier in `XCmpQuantify::countVariants` where
/// non-CONF truth paired with a local query match short-circuits to
/// UNK+lm rather than claiming an FN against a non-evaluable region.
pub(super) fn unk_truth_row(
    truth: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
    bk: &'static str,
) -> AnnotatedRow {
    let info = comparison_info(truth, reference);
    AnnotatedRow {
        sort_key: (truth.key.chrom.clone(), truth.key.pos, 0, 0),
        // UNK truth rows are not counted against the pass-only tier (same
        // semantics as FN rows, per legacy's Counts.cpp treatment).
        query_pass: true,
        fp_class: None,
        xcmp_ctype: None,
        xcmp_hap_match: false,
        record: comparison_record(
            truth,
            display_ref(truth, reference),
            display_alt(truth),
            truth.qual.clone(),
            ".".to_string(),
            format!("BS={block_start}{regions}"),
            format!(
                "{}:UNK:{bk}:{info}:{}:{}:.",
                truth.gt,
                truth.primary_type(),
                genotype_label(truth)
            ),
            "./.:.:.:.:NOCALL:nocall:0".to_string(),
        ),
    }
}

/// Split a multi-allelic query variant into per-active-allele primitives,
/// trimming shared prefix/suffix so the emitted VCF records land at the
/// natural anchor position (e.g. `GAAGA → GAAGAAAGA,G` with GT=1/2 splits
/// into `A → AAAGA` at pos+4 plus `GAAGA → G` at pos). Each primitive
/// keeps the cluster-level classification (TP/FP/UNK) that was decided
/// for the parent record. Single-allelic variants pass through unchanged.
/// Legacy emits one VCF row per primitive; without this, rust under-counts
/// QUERY.TOTAL by ~40 records per chr21 case.
/// Variant of `split_query_primitives` that knows about sibling records
/// in the cluster. The neighbor list lets the per-primitive left-shift
/// (Class F) avoid sliding onto a position already occupied by another
/// raw record — sliding chr21:28720259 alt 2 onto 28720258 would
/// duplicate the raw `28720258 AAAG→A` record because that single-alt
/// deletion sits at the same locus the slide would target.
///
/// `cluster_truth` enables the Class B same-anchor insertion fan-out
/// (`try_split_same_anchor_via_shift`) — when a multi-allelic insertion's
/// shifted primitive lands at a position represented in truth, the
/// fan-out is suppressed so block-level haplotype matching can reconcile
/// the multi-allelic record (chr21:40096658 case). Pass `&[]` when truth
/// context is unavailable.
pub(super) fn split_query_primitives_with_neighbors(
    variant: &Variant,
    reference: &str,
    cluster_start: usize,
    cluster_neighbors: &[Variant],
    cluster_truth: &[Variant],
) -> Vec<Variant> {
    if !variant.key.alt_allele.contains(',') {
        return vec![variant.clone()];
    }
    let alleles = parse_gt_alleles(&variant.gt);
    let used: Vec<usize> = alleles
        .iter()
        .copied()
        .filter(|a| *a > 0)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    if used.is_empty() {
        return vec![variant.clone()];
    }
    let alts: Vec<&str> = variant.key.alt_allele.split(',').collect();

    // Trim each active allele to its minimal (pos, ref, alt) triple. If all
    // trimmed primitives share the same (pos, ref) — i.e. all alts are
    // length-1 variations or same-length anchored events at the same
    // reference position — legacy keeps the record multi-allelic. Only
    // when trimming produces distinct anchor positions does legacy fan
    // the record out into per-primitive rows.
    //
    // Class B (left-shift through microsatellites — chr21:21690513) is
    // INTENTIONALLY OUT OF SCOPE here for INSERTION primitives: applying
    // `partial_credit::left_shift` per insertion over-shifts microsat
    // inserts already canonicalized by pre.py (regressed
    // chr21_passonly_region from 0 → 60 in a prior session). Class F
    // (chr21:44413761 + chr21:47906004) IS handled by the `apply_slide`
    // block below for DELETION primitives only, gated on the absence of
    // a neighboring raw record at the slide target — sliding onto an
    // already-occupied locus would duplicate that record.
    let mut trimmed: Vec<(usize, String, String)> = Vec::new();
    for allele_idx in &used {
        let Some(alt) = alts.get(*allele_idx - 1).copied() else {
            continue;
        };
        trimmed.push(trim_variant(variant.key.pos, &variant.key.ref_allele, alt));
    }
    if trimmed.is_empty() {
        return vec![variant.clone()];
    }

    let first_anchor = (&trimmed[0].0, &trimmed[0].1);
    let all_same_anchor = trimmed
        .iter()
        .all(|(pos, r, _)| (pos, r) == (first_anchor.0, first_anchor.1));
    if all_same_anchor {
        // Class B same-anchor insertion fan-out. When the multi-allelic
        // primitives canonicalize to distinct anchors via `partial_credit::
        // left_shift` AND no shifted alt has truth representation at the
        // shifted anchor, fan out into per-primitive rows with the slide
        // applied. chr21:21690513 (`C→CACAC,CACAT`): CACAC slides to
        // (21690501, T, TACAC), CACAT stays at 21690513 with truth's
        // exact match → fan out. chr21:40096658 (`T→TAGATAGAG,TAGATAGAT`):
        // TAGATAGAT slides to (40096650, C, CAGATAGAT), but truth at
        // 40096650 already declares that allele → keep multi-allelic.
        let bases_for_class_b = reference.as_bytes();
        let pos_min_for_class_b = cluster_start.saturating_sub(SPLIT_LEFT_SHIFT_WINDOW).max(1);
        if let Some(shifted) = try_split_same_anchor_via_shift(
            &variant.key.chrom,
            variant.key.pos,
            &trimmed,
            bases_for_class_b,
            pos_min_for_class_b,
            cluster_truth,
            cluster_neighbors,
        ) {
            let mut sorted = shifted;
            sorted.sort_by_key(|(p, r, _)| (*p, r.len()));
            return sorted
                .into_iter()
                .map(|(p, r, a)| Variant {
                    key: VariantKey {
                        chrom: variant.key.chrom.clone(),
                        pos: p,
                        ref_allele: r,
                        alt_allele: a,
                    },
                    qual: variant.qual.clone(),
                    filter: variant.filter.clone(),
                    gt: "0/1".to_string(),
                })
                .collect();
        }
        // Keep as multi-allelic; sort alphabetically so query-only
        // records emit `T → TCA,TCACACA` (shorter first) matching the
        // canonical representation legacy pre.py generates. When a
        // paired truth exists the exact-match branch switches to
        // truth's representation, so this sort only affects query-only
        // emission.
        let (pos, ref_allele, _) = &trimmed[0];
        let mut alt_list: Vec<String> = trimmed.iter().map(|(_, _, alt)| alt.clone()).collect();
        alt_list.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
        return vec![Variant {
            key: VariantKey {
                chrom: variant.key.chrom.clone(),
                pos: *pos,
                ref_allele: ref_allele.clone(),
                alt_allele: alt_list.join(","),
            },
            qual: variant.qual.clone(),
            filter: variant.filter.clone(),
            // Canonical hetalt GT for the sorted multi-allelic. Legacy
            // emits the reversed form `<later>/<earlier>` (e.g. `2/1`)
            // because `VariantLocationAggregator::addAlleleToVariant`
            // with `MAX_GT=2` writes the second allele into the first
            // zero slot, producing reversed-index GTs. Mirror that
            // ordering so downstream byte output matches.
            gt: if alt_list.len() > 1 {
                "2/1".to_string()
            } else {
                "0/1".to_string()
            },
        }];
    }

    // Class F per-primitive left-shift (deletions only). Sort by start so
    // the accumulated `local_min` boundary covers siblings left-to-right.
    // Skip the slide entirely when a neighboring cluster record already
    // sits at or before the trim anchor — those are the chr21:28720259 /
    // chr21:30548467 shapes where legacy keeps each raw record at its own
    // anchor and a per-primitive slide would create duplicates.
    let bases = reference.as_bytes();
    let mut sorted = trimmed;
    sorted.sort_by_key(|(pos, ref_allele, _)| (*pos, ref_allele.len()));
    let neighbor_anchors: BTreeSet<usize> = cluster_neighbors
        .iter()
        .filter(|n| n.key.pos != variant.key.pos)
        .map(|n| n.key.pos)
        .collect();
    let any_neighbor_below = sorted
        .iter()
        .any(|(pos, _, _)| neighbor_anchors.iter().any(|&n| n < *pos));
    let mut local_min = cluster_start.saturating_sub(SPLIT_LEFT_SHIFT_WINDOW).max(1);
    sorted
        .into_iter()
        .map(|(pos, ref_allele, alt_allele)| {
            let apply_slide = !any_neighbor_below
                && ref_allele.len() > alt_allele.len()
                && pos > 0
                && pos + ref_allele.len() - 1 <= bases.len();
            let (final_pos, final_ref, final_alt) = if apply_slide {
                let mut rv = partial_credit::RefVar {
                    start: pos,
                    end: pos + ref_allele.len() - 1,
                    alt: alt_allele.clone(),
                };
                partial_credit::left_shift(bases, &mut rv, local_min, true);
                let ref_start = rv.start.saturating_sub(1);
                let ref_end = rv.end;
                if ref_end >= ref_start && ref_end <= bases.len() && rv.start >= 1 {
                    let new_ref: String = bases[ref_start..ref_end]
                        .iter()
                        .map(|&b| b.to_ascii_uppercase() as char)
                        .collect();
                    let new_end = rv.start + new_ref.len().max(1) - 1;
                    local_min = local_min.max(new_end);
                    (rv.start, new_ref, rv.alt)
                } else {
                    let original_end = pos + ref_allele.len().max(1) - 1;
                    local_min = local_min.max(original_end);
                    (pos, ref_allele, alt_allele)
                }
            } else {
                let original_end = pos + ref_allele.len().max(1) - 1;
                local_min = local_min.max(original_end);
                (pos, ref_allele, alt_allele)
            };
            Variant {
                key: VariantKey {
                    chrom: variant.key.chrom.clone(),
                    pos: final_pos,
                    ref_allele: final_ref,
                    alt_allele: final_alt,
                },
                qual: variant.qual.clone(),
                filter: variant.filter.clone(),
                // Each split primitive becomes a simple 0/1 het so downstream
                // is_het / genotype_label / primary_type classify it the same
                // way legacy does for its decomposed representations.
                gt: "0/1".to_string(),
            }
        })
        .collect()
}

/// Trim common prefix and suffix from (ref, alt) and re-anchor the
/// position so both sides end up non-empty. Matches VCF primitive
/// normalisation: insertions anchor on the last prefix base, deletions
/// anchor on the same, substitutions just shift `pos` by the prefix
/// length.
pub(super) fn trim_variant(pos: usize, ref_allele: &str, alt: &str) -> (usize, String, String) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt.as_bytes();
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
    let ref_mid_start = prefix;
    let ref_mid_end = ref_bytes.len() - suffix;
    let alt_mid_start = prefix;
    let alt_mid_end = alt_bytes.len() - suffix;
    let mut new_ref: String = ref_bytes[ref_mid_start..ref_mid_end]
        .iter()
        .map(|b| *b as char)
        .collect();
    let mut new_alt: String = alt_bytes[alt_mid_start..alt_mid_end]
        .iter()
        .map(|b| *b as char)
        .collect();
    let mut new_pos = pos + prefix;
    if new_ref.is_empty() || new_alt.is_empty() {
        if suffix > 0 {
            // Right-anchor: use first char of shared suffix. Legacy right-anchors
            // indels when possible (anchor appended rather than prepended).
            let anchor = ref_bytes[ref_bytes.len() - suffix] as char;
            new_ref.push(anchor);
            new_alt.push(anchor);
        } else if prefix > 0 {
            // Fallback left-anchor when no suffix is available.
            new_pos = new_pos.saturating_sub(1);
            let anchor = ref_bytes[prefix - 1] as char;
            new_ref = format!("{}{}", anchor, new_ref);
            new_alt = format!("{}{}", anchor, new_alt);
        }
    }
    (new_pos, new_ref, new_alt)
}

/// Per-row BK classifier ported from legacy hap.py's
/// `XCmpQuantify::countVariants` (lines 188–202) plus its upstream `kind`
/// and `ctype` assignments in `SimpleDiploidCompare.cpp` (lines 260–265)
/// and `xcmp.cpp` (lines 440–525). Returns `"lm"` for rows that hit
/// either of legacy's two local-match paths, `"."` otherwise:
///
/// * **Path A — `almismatch`** (`SimpleDiploidCompare.cpp:265`): a
///   counterpart-side record exists at the *exact same* (chrom, pos, ref)
///   with a GT-selected non-reference allele set that is completely
///   disjoint from this row's selected alt set. Legacy requires both
///   sides to be present and both alt sets to be non-empty; an empty set
///   intersects no set, so a hom-ref counterpart never triggers
///   `almismatch`. Filtered counterparts are excluded — simple-compare
///   operates on post-filter calls.
///
/// * **Path B — `hap:mismatch`** (`xcmp.cpp:440-525`): the block-level
///   haplotype comparator was actually invoked (the `n_nonsnp > 0` gate
///   fired on the cluster) AND the reconstructed signatures on the two
///   sides failed to reconcile. Pure-SNP clusters do not qualify because
///   the hapcmp gate won't fire on them, so this path only ever hits
///   clusters containing at least one GT-selected indel / MNP.
///
/// `gm` is produced by TP emitters (`mark_cluster_match` /
/// `exact_match_pairs`) before a row ever reaches this helper, and `am`
/// is written directly by `fn_fp_combined_row` for GT-mismatch pairs —
/// both legacy's `type == "TP"` and `kind == "gtmismatch"` branches are
/// handled outside this function. The `else kind = "."` fallthrough at
/// the bottom of legacy's chain corresponds to returning `"."` here.
pub(super) fn bk_for_row(
    row: &Variant,
    counterparts: &[Variant],
    hap_mismatch: bool,
) -> &'static str {
    if almismatch_same_locus(row, counterparts) {
        return "lm";
    }
    if hap_mismatch {
        return "lm";
    }
    "."
}

/// Does the counterpart side declare a record at the same VCF locus with
/// a completely disjoint GT-selected non-reference allele set? Mirrors
/// legacy's `(nonref_als_1 & nonref_als_2) == 0` check.
///
/// Legacy's `compareVariants` only fires the `almismatch` path when it
/// sees a **single merged `Variants` record** with non-ref calls on both
/// samples. `VariantReader` groups VCF rows that share the full
/// `(chrom, pos, ref, alt)` tuple — rows with different ALT columns
/// stay as separate Variants objects. So rust's counterpart scan must
/// also require an exact ALT-column match, not just the same
/// `(chrom, pos, ref)`. Without the alt check, rust flags cases where
/// legacy's reader kept the rows separate (e.g. `G→GA` insertion on one
/// side vs `G→A` SNP on the other at the same position) and both sides
/// end up with `kind=missing` in legacy rather than `almismatch`.
pub(super) fn almismatch_same_locus(row: &Variant, counterparts: &[Variant]) -> bool {
    let row_alts = gt_selected_nonref_alts(row);
    if row_alts.is_empty() {
        return false;
    }
    // Legacy's `compareVariants` only sees per-allele mismatch when
    // `VariantReader` merged the two records into a single `Variants`
    // object — which only happens when both sides' REF and ALT lengths
    // are compatible as a substitution (same ref length, same alt
    // length, i.e. both pure SNPs). Cross-shape pairs (e.g. truth
    // insertion `G→GA` vs query SNP `G→A` at the same pos) stay as
    // separate Variants rows and each side's compareVariants returns
    // `kind=missing`, not `almismatch` → BK=`.`.
    let row_snp_like =
        row.key.ref_allele.len() == 1 && row.key.alt_allele.split(',').all(|a| a.len() == 1);
    for c in counterparts {
        if c.key.chrom != row.key.chrom
            || c.key.pos != row.key.pos
            || c.key.ref_allele != row.key.ref_allele
            || c.key.alt_allele != row.key.alt_allele
        {
            continue;
        }
        let c_snp_like =
            c.key.ref_allele.len() == 1 && c.key.alt_allele.split(',').all(|a| a.len() == 1);
        if row_snp_like != c_snp_like {
            continue;
        }
        if !filter_is_pass(&c.filter) {
            continue;
        }
        let c_alts = gt_selected_nonref_alts(c);
        if c_alts.is_empty() {
            continue;
        }
        if c_alts.is_disjoint(&row_alts) {
            return true;
        }
    }
    false
}

pub(super) fn fp_like_row(
    query: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
    bd: &str,
    fp_class: Option<&'static str>,
    bk: &'static str,
) -> AnnotatedRow {
    let info = comparison_info(query, reference);
    AnnotatedRow {
        sort_key: (
            query.key.chrom.clone(),
            query.key.pos,
            2,
            query_type_rank(query),
        ),
        query_pass: filter_is_pass(&query.filter),
        fp_class,
        xcmp_ctype: None,
        xcmp_hap_match: false,
        record: comparison_record(
            query,
            display_ref(query, reference),
            display_alt(query),
            query.qual.clone(),
            filter_for_output(&query.filter).to_string(),
            format!("BS={block_start}{regions}"),
            "./.:.:.:.:NOCALL:nocall:.".to_string(),
            format!(
                "{}:{bd}:{bk}:{info}:{}:{}:{}",
                query.gt,
                query.primary_type(),
                genotype_label(query),
                query.qual
            ),
        ),
    }
}

pub(super) fn display_alt(variant: &Variant) -> String {
    // Preserve the source VCF's ALT ordering rather than sorting lexicographically.
    // Legacy hap.py passes the original multi-allelic ALT string through xcmp
    // and writes it back verbatim; sorting here would rewrite the ALT column
    // (e.g. `ATCTC,ATC → ATC,ATCTC`) and also invalidate the GT indices that
    // were emitted for the pre-sort allele order.
    variant.key.alt_allele.clone()
}

pub(super) fn display_ref(variant: &Variant, reference: &str) -> String {
    if variant
        .key
        .alt_allele
        .split(',')
        .any(|alt| alt.starts_with('<') && alt.ends_with('>'))
        && let Some(base) = reference.as_bytes().get(variant.key.pos.saturating_sub(1))
    {
        return (*base as char).to_ascii_uppercase().to_string();
    }
    variant.key.ref_allele.clone()
}

pub(super) fn query_type_rank(variant: &Variant) -> usize {
    if variant.primary_type() == "SNP" {
        0
    } else {
        1
    }
}

pub(super) fn genotype_label(variant: &Variant) -> &'static str {
    // BLT / location label follows the active-allele view of the GT, not the
    // literal string: two distinct non-REF alleles are hetalt; two equal
    // non-REF alleles are homalt (covers 1|1, 2|2, 3|3 …); one REF + one
    // non-REF is het (covers 0|1, 1|0, 0|3, 3|0 …); anything else nocall.
    // Kept in-sync with genotype_label without going through is_het /
    // is_homalt, which are narrower on purpose so summary counts stay
    // byte-equal with legacy's literal-GT-match bucketing.
    let all_alleles = parse_gt_alleles(&variant.gt);
    let nonzero: Vec<usize> = all_alleles.iter().copied().filter(|a| *a > 0).collect();
    if nonzero.len() == 2 && nonzero[0] != nonzero[1] {
        return "hetalt";
    }
    if all_alleles.len() == 2 && all_alleles[0] > 0 && all_alleles[0] == all_alleles[1] {
        return "homalt";
    }
    if all_alleles.len() == 2 && (all_alleles[0] == 0) != (all_alleles[1] == 0) {
        return "het";
    }
    "nocall"
}

pub(super) fn comparison_info(variant: &Variant, reference: &str) -> String {
    // Multi-allelic records — SNP or INDEL — go through subtype_label so
    // mixed ti/tv or mixed indel sizes emit the legacy comma-joined BI
    // (e.g. `A → G,T` GT=2/1 must be `ti,tv`, not a single ti/tv token).
    if variant.key.alt_allele.contains(',')
        && let Some(subtype) = subtype_label(variant)
    {
        return subtype.to_lowercase();
    }
    if variant.primary_type() == "SNP" {
        return snp_bucket_label(variant).unwrap_or("tv").to_string();
    }
    if let Some(subtype) = subtype_label(variant) {
        return subtype.to_lowercase();
    }
    let _ = reference;
    "indel".to_string()
}

pub(super) fn snp_bucket_label(variant: &Variant) -> Option<&'static str> {
    if variant.primary_type() != "SNP" {
        return None;
    }
    if variant.is_single_base_snp() {
        return Some(if variant.is_transition() { "ti" } else { "tv" });
    }
    let mut saw_change = false;
    let mut saw_ti = false;
    let mut saw_tv = false;
    for (ref_base, alt_base) in variant
        .key
        .ref_allele
        .chars()
        .zip(variant.key.alt_allele.chars())
    {
        if ref_base == alt_base {
            continue;
        }
        saw_change = true;
        if is_transition_pair(ref_base, alt_base) {
            saw_ti = true;
        } else {
            saw_tv = true;
        }
    }
    if !saw_change {
        None
    } else if saw_ti && !saw_tv {
        Some("ti")
    } else if saw_tv && !saw_ti {
        Some("tv")
    } else {
        None
    }
}

pub(super) fn subtype_label(variant: &Variant) -> Option<String> {
    // Multi-allelic records may mix SNP + INDEL alts (e.g. `A → AT,T` is
    // an insertion + an SNP). Legacy emits the BI field as a comma-
    // joined list of each allele's subtype token — `i1_5,tv` in that
    // example. Don't short-circuit on primary_type=INDEL at the top;
    // handle each allele individually first.
    if variant.key.alt_allele.contains(',') {
        let used = parse_gt_alleles(&variant.gt)
            .into_iter()
            .filter(|allele| *allele > 0)
            .collect::<BTreeSet<_>>();
        let alts: Vec<&str> = variant.key.alt_allele.split(',').collect();
        let mut subtypes: BTreeSet<String> = BTreeSet::new();
        for allele in used {
            let alt = *alts.get(allele.saturating_sub(1))?;
            let pseudo = Variant {
                key: VariantKey {
                    chrom: variant.key.chrom.clone(),
                    pos: variant.key.pos,
                    ref_allele: variant.key.ref_allele.clone(),
                    alt_allele: alt.to_string(),
                },
                qual: variant.qual.clone(),
                filter: variant.filter.clone(),
                gt: "0/1".to_string(),
            };
            // Per-allele subtype: INDELs get size-bucketed (I1_5 etc),
            // SNPs get ti/tv. Both contribute to the comma-joined BI.
            let one = if pseudo.primary_type() == "SNP" {
                snp_bucket_label(&pseudo).map(|s| s.to_string())
            } else {
                subtype_label(&pseudo)
            };
            if let Some(token) = one {
                subtypes.insert(token);
            }
        }
        if subtypes.is_empty() {
            return None;
        }
        return Some(subtypes.into_iter().collect::<Vec<_>>().join(","));
    }
    if variant.primary_type() != "INDEL" {
        return None;
    }
    let ref_allele = &variant.key.ref_allele;
    let alt_allele = &variant.key.alt_allele;

    // Symbolic alts (`<DEL>`, `<INS>`, `<DUP>`) don't align character-by-
    // character against the ref. Legacy's `realignRefVar` treats them as
    // pure del/ins of the REF length (anchor-base convention — at least
    // one REF base). Map the tag directly; size comes from ref_len.max(1)
    // (matches legacy's `T → <DEL>` → `d1_5` classification).
    if alt_allele.starts_with('<') && alt_allele.ends_with('>') {
        let inner = &alt_allele[1..alt_allele.len() - 1];
        let class = if inner.starts_with("DEL") {
            "D"
        } else if inner.starts_with("INS") || inner.starts_with("DUP") {
            "I"
        } else {
            "C"
        };
        let bucket = bucket_label(ref_allele.len().max(1));
        return Some(format!("{class}{bucket}"));
    }

    // Two-sided trim: legacy uses Smith-Waterman-style realignment to
    // decide whether a ref/alt pair is a pure del/ins after normalising
    // both shared prefix AND shared suffix. Rust previously only trimmed
    // shared prefix — that misclassifies e.g. `GGAAAGAAAAAGAAA… →
    // GGAAAGAAAGAAA` as complex when it's really a plain deletion of
    // 10 bases (legacy emits `d6_15`). Trimming both sides recovers the
    // canonical I/D/C classification.
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();
    let prefix = ref_bytes
        .iter()
        .zip(alt_bytes.iter())
        .take_while(|(r, a)| r == a)
        .count();
    let max_suffix = (ref_bytes.len() - prefix).min(alt_bytes.len() - prefix);
    let suffix = (0..max_suffix)
        .take_while(|i| ref_bytes[ref_bytes.len() - 1 - i] == alt_bytes[alt_bytes.len() - 1 - i])
        .count();
    let ref_rem = ref_bytes.len() - prefix - suffix;
    let alt_rem = alt_bytes.len() - prefix - suffix;

    let (class, size) = if alt_rem == 0 && ref_rem > 0 {
        ("D", ref_rem)
    } else if ref_rem == 0 && alt_rem > 0 {
        ("I", alt_rem)
    } else {
        // Complex: size is the larger remnant (matches legacy's
        // VariantStatistics bucketing on total_ins+total_del after
        // realignment, which for a substitution-plus-indel takes the
        // dominant indel length).
        ("C", ref_rem.max(alt_rem))
    };
    let bucket = bucket_label(size);
    let mut label = format!("{class}{bucket}");
    if class == "C" && ref_rem > 0 && alt_rem > 0 {
        // Complex variants carry a comma-joined ti/tv tag derived from
        // the first post-trim ref/alt pair, matching legacy
        // `VariantStatistics`'s additional Ti/Tv bucket on
        // substitution-plus-indel records (e.g. chr21:26105569
        // `T→AAAGAAAA` is `c6_15,tv` because T→A is a transversion;
        // chr21:32759510 allele 2 `TTTTTTTTT→C` is `c6_15,ti`). The
        // pair at index `prefix` is guaranteed mismatched — `prefix`
        // counts shared leading bases and stops on the first
        // disagreement.
        let ref_first = ref_bytes[prefix] as char;
        let alt_first = alt_bytes[prefix] as char;
        let ti_tv = if is_transition_pair(ref_first, alt_first) {
            "ti"
        } else {
            "tv"
        };
        label.push(',');
        label.push_str(ti_tv);
    }
    Some(label)
}

pub(super) fn bucket_label(size: usize) -> &'static str {
    match size {
        0..=5 => "1_5",
        6..=15 => "6_15",
        _ => "16_PLUS",
    }
}

pub(super) fn is_transition_pair(ref_base: char, alt_base: char) -> bool {
    matches!(
        (ref_base, alt_base),
        ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C')
    )
}
