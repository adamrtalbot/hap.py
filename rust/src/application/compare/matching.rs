//! Cluster nearby variants and match truth against query within each superlocus.

#[cfg(test)]
use super::XCMP_ENUMERATION_THRESHOLD;
use super::genotype::{equivalent_gt, parse_gt_alleles};
use super::legacy_graph;
use super::metrics::{add_variant_stats, add_variant_stats_subtype};
use super::rows::{
    BenchmarkDecision, bk_for_row, cluster_query_filter, compute_shift_target, fn_fp_combined_row,
    fn_row, fp_like_row, gt_selected_nonref_alts, legacy_duplicate_alt_query_output_projection,
    split_query_primitives_with_neighbors, tp_combined_row, tp_single_side_row, trim_variant,
    try_split_same_anchor_via_shift, unk_combined_row, unk_truth_row,
};
use super::{
    AnnotatedRow, Cluster, ComparisonConfig, ComparisonOutputs, Entry, Event, MAX_CLUSTER_VARIANTS,
    RegionState, SPLIT_LEFT_SHIFT_WINDOW, Side, row_matches_variant_key,
};
use crate::adapters::vcf::{Variant, VariantKey};
use crate::domain::{FpClass, Interval, TypeCounts, XcmpCtype};
use anyhow::{Result, bail};
use std::borrow::Borrow;
use std::collections::{BTreeMap, BTreeSet};

#[cfg(test)]
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
                    && start <= cluster.end.saturating_add(cluster_gap) =>
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

pub(super) struct StreamingClusters<T, Q>
where
    T: Iterator<Item = Result<Variant>>,
    Q: Iterator<Item = Result<Variant>>,
{
    truth: std::iter::Peekable<T>,
    query: std::iter::Peekable<Q>,
    pending: Option<Entry>,
    cluster_gap: usize,
    contig_ranks: BTreeMap<String, usize>,
}

impl<T, Q> StreamingClusters<T, Q>
where
    T: Iterator<Item = Result<Variant>>,
    Q: Iterator<Item = Result<Variant>>,
{
    pub(super) fn new(
        truth: T,
        query: Q,
        cluster_gap: usize,
        contig_ranks: BTreeMap<String, usize>,
    ) -> Self {
        Self {
            truth: truth.peekable(),
            query: query.peekable(),
            pending: None,
            cluster_gap,
            contig_ranks,
        }
    }

    fn next_entry(&mut self) -> Result<Option<Entry>> {
        if self.truth.peek().is_some_and(Result::is_err) {
            return self.truth.next().transpose().map(|variant| {
                variant.map(|variant| Entry {
                    side: Side::Truth,
                    variant,
                })
            });
        }
        if self.query.peek().is_some_and(Result::is_err) {
            return self.query.next().transpose().map(|variant| {
                variant.map(|variant| Entry {
                    side: Side::Query,
                    variant,
                })
            });
        }
        let take_truth = match (self.truth.peek(), self.query.peek()) {
            (Some(Ok(truth)), Some(Ok(query))) => {
                variant_stream_order(truth, query, &self.contig_ranks).is_le()
            }
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => return Ok(None),
            (Some(_), Some(_)) => unreachable!("errors handled before ordering"),
        };
        if take_truth {
            self.truth.next().transpose().map(|variant| {
                variant.map(|variant| Entry {
                    side: Side::Truth,
                    variant,
                })
            })
        } else {
            self.query.next().transpose().map(|variant| {
                variant.map(|variant| Entry {
                    side: Side::Query,
                    variant,
                })
            })
        }
    }
}

impl<T, Q> Iterator for StreamingClusters<T, Q>
where
    T: Iterator<Item = Result<Variant>>,
    Q: Iterator<Item = Result<Variant>>,
{
    type Item = Result<Cluster>;

    fn next(&mut self) -> Option<Self::Item> {
        let first = match self.pending.take() {
            Some(entry) => entry,
            None => match self.next_entry() {
                Ok(Some(entry)) => entry,
                Ok(None) => return None,
                Err(error) => return Some(Err(error)),
            },
        };
        let mut cluster = Cluster {
            chrom: first.variant.key.chrom.clone(),
            start: first.variant.key.pos,
            end: first.variant.end_pos(),
            truth: Vec::new(),
            query: Vec::new(),
        };
        push_cluster_entry(&mut cluster, first);
        loop {
            let entry = match self.next_entry() {
                Ok(Some(entry)) => entry,
                Ok(None) => return Some(Ok(cluster)),
                Err(error) => return Some(Err(error)),
            };
            let start = entry.variant.key.pos;
            if cluster.chrom == entry.variant.key.chrom
                && start <= cluster.end.saturating_add(self.cluster_gap)
            {
                if cluster.truth.len() + cluster.query.len() >= MAX_CLUSTER_VARIANTS {
                    return Some(Err(anyhow::anyhow!(
                        "connected comparison cluster {}:{}-{} exceeds the {} variant active-window limit",
                        cluster.chrom,
                        cluster.start,
                        cluster.end.max(entry.variant.end_pos()),
                        MAX_CLUSTER_VARIANTS
                    )));
                }
                cluster.end = cluster.end.max(entry.variant.end_pos());
                push_cluster_entry(&mut cluster, entry);
            } else {
                self.pending = Some(entry);
                return Some(Ok(cluster));
            }
        }
    }
}

fn variant_stream_order(
    left: &Variant,
    right: &Variant,
    contig_ranks: &BTreeMap<String, usize>,
) -> std::cmp::Ordering {
    match (
        contig_ranks.get(&left.key.chrom),
        contig_ranks.get(&right.key.chrom),
    ) {
        (Some(left_rank), Some(right_rank)) => left_rank.cmp(right_rank),
        (Some(_), None) => std::cmp::Ordering::Less,
        (None, Some(_)) => std::cmp::Ordering::Greater,
        (None, None) => left.key.chrom.cmp(&right.key.chrom),
    }
    .then(left.key.pos.cmp(&right.key.pos))
    .then(left.end_pos().cmp(&right.end_pos()))
}

pub(super) fn comparison_contig_ranks(
    truth_headers: &[String],
    query_headers: &[String],
) -> BTreeMap<String, usize> {
    let mut ranks = BTreeMap::new();
    for line in truth_headers.iter().chain(query_headers) {
        let Some(body) = line.strip_prefix("##contig=<ID=") else {
            continue;
        };
        let Some(contig) = body.split([',', '>']).next() else {
            continue;
        };
        let next_rank = ranks.len();
        ranks.entry(contig.to_string()).or_insert(next_rank);
    }
    ranks
}

fn push_cluster_entry(cluster: &mut Cluster, entry: Entry) {
    match entry.side {
        Side::Truth => cluster.truth.push(entry.variant),
        Side::Query => cluster.query.push(entry.variant),
    }
}

/// Order one side's block records the way legacy xcmp sees them.
///
/// legacy's `VariantReader` pulls records through htslib's synced reader,
/// whose `bcf_sr_sort` pass regroups every line sharing a VCF POS before
/// they reach `GraphReference::makeGraph`. With `COLLAPSE_NONE` the pairing
/// mode is `BCF_SR_PAIR_EXACT`, so each distinct REF/ALT spelling becomes one
/// variant set whose count is the number of *files* carrying it, and the
/// emission loop repeatedly pops the highest-count set (ties keep creation
/// order). A REF/ALT spelling present in both inputs therefore precedes the
/// single-file spellings at the same POS, independent of the order the lines
/// appear in either file.
///
/// The pairing is on allele spelling alone — genotypes never enter it — so
/// the key here is "does the other side of this block carry a record with the
/// same (pos, ref, alt)", not "did this record survive rust's GT-aware exact
/// match".
///
/// Inside a class the stream keeps each reader's own line order, and the two
/// readers do not agree on it. hap.py leaves the truth VCF unpreprocessed
/// (`preprocessing_truth` is off), so truth records reach xcmp in raw input
/// order. The query goes through pre.py's leftshift/decompose, which re-emits
/// same-position records by trimmed edit start — the coordinate `makeGraph`
/// derives after `trimLeft`/`trimRight`, putting a substitution at POS ahead
/// of an insertion anchored on it. Measured against the legacy result VCF that
/// holds for every one of its 7,566 multi-record query-only positions, while
/// truth-only positions keep raw order in the cases where the two disagree.
pub(super) fn legacy_graph_truth_order(
    truth: &[Variant],
    paired_keys: &BTreeSet<VariantKey>,
) -> Vec<Variant> {
    let mut ordered = truth.to_vec();
    ordered.sort_by_key(|variant| {
        (
            variant.key.pos,
            usize::from(!paired_keys.contains(&variant.key)),
        )
    });
    ordered
}

/// Query-side counterpart of [`legacy_graph_truth_order`]. Paired records take
/// their position from the truth stream — htslib creates each variant set while
/// scanning reader 0 (truth) first, so the pairs are emitted in truth's order —
/// and the query-only remainder follows pre.py's trimmed-start order.
pub(super) fn legacy_graph_query_order(
    query: &[Variant],
    paired_keys: &BTreeSet<VariantKey>,
    truth_order: &[Variant],
) -> Vec<Variant> {
    let truth_rank = truth_order
        .iter()
        .enumerate()
        .map(|(rank, variant)| (variant.key.clone(), rank as isize))
        .collect::<BTreeMap<_, _>>();
    let mut ordered = query.to_vec();
    ordered.sort_by_key(|variant| {
        let paired = paired_keys.contains(&variant.key);
        (
            variant.key.pos,
            usize::from(!paired),
            if paired {
                truth_rank.get(&variant.key).copied().unwrap_or(0)
            } else {
                legacy_graph::selected_edit_start(variant)
            },
        )
    });
    ordered
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
        set_xcmp_context(&mut rows[output_start..], XcmpCtype::SimpleMatch, false);
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
    let (linear_truth_sig, linear_query_sig) = if allow_haplotype_match {
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
    let (mut truth_sig, mut query_sig) = (linear_truth_sig.clone(), linear_query_sig.clone());
    if allow_haplotype_match {
        // htslib pairs the two inputs' lines on allele spelling alone, so the
        // "seen on both sides" class is the intersection of the block's
        // (pos, ref, alt) keys — genotypes play no part.
        let truth_keys = cluster
            .truth
            .iter()
            .map(|variant| variant.key.clone())
            .collect::<BTreeSet<_>>();
        let query_keys = cluster
            .query
            .iter()
            .map(|variant| variant.key.clone())
            .collect::<BTreeSet<_>>();
        let paired_keys = truth_keys
            .intersection(&query_keys)
            .cloned()
            .collect::<BTreeSet<_>>();
        let graph_truth = legacy_graph_truth_order(&cluster.truth, &paired_keys);
        let graph_query = legacy_graph_query_order(&cluster.query, &paired_keys, &graph_truth);
        truth_sig = legacy_graph::signatures(
            &graph_truth,
            reference,
            signature_cluster.start,
            signature_cluster.end,
            config.max_enum,
        );
        query_sig = legacy_graph::signatures(
            &graph_query,
            reference,
            signature_cluster.start,
            signature_cluster.end,
            config.max_enum,
        );
    }

    // Legacy's `DiploidCompare::setRegion` flags `hap_match = true` as
    // soon as ANY enumerated (h1, h2) pair from truth's di_haps equals
    // ANY pair from query's di_haps — cross-product membership, not set
    // equality. Replicate that rule here: mismatch only when the two
    // signature sets are disjoint. `None` on either side (enumeration
    // budget exceeded) still falls through to the mismatch path because
    // legacy's `hap_fail = true` takes the same `ctype != "hap:mismatch"`
    // branch we want for non-evaluable blocks.
    let overlapping_deletion_mismatch = legacy_overlapping_deletion_mismatch(cluster);
    let signatures_overlap = matches!(
        (&truth_sig, &query_sig),
        (Some(left), Some(right)) if left.intersection(right).next().is_some()
    );
    let graph_failed = allow_haplotype_match && (truth_sig.is_none() || query_sig.is_none());
    // The overlapping-deletion GT discordance is a mismatch-promotion signal,
    // not a match veto. When the truth and query signature sets share a diploid
    // haplotype pair, legacy's DiploidCompare returns hap:match regardless of
    // the deletion-overlap shape, so `is_match` keys off the signatures alone.
    // `overlapping_deletion_mismatch` still promotes an *unmatched* block to
    // BK=lm via `hap_mismatch` below (chr20 fixture drains to hapfail without
    // it); it just no longer overrides a genuine match (HG003 DeepVariant
    // chr11:7757959 adjacent/compound INDEL block, which legacy resolves gm).
    let is_match = signatures_overlap && !graph_failed;
    // Legacy's `ctype == "hap:mismatch"` fires only when the block-level
    // haplotype comparator actually ran (both signatures computed) AND
    // the two sides disagreed. A None signature — hapcmp skipped on the
    // n_nonsnp gate OR state-count cap exceeded — corresponds to
    // legacy's "simple" / "hapfail" ctypes, neither of which promotes to
    // BK=lm.
    let insert_conflict_counterpart = if allow_haplotype_match
        && estimated_state_count_with_limit(&cluster.query, config.max_enum) <= config.max_enum
    {
        query_insert_conflict_has_truth_counterpart(
            &cluster.query,
            &cluster.truth,
            &truth_remaining,
        )
    } else {
        None
    };
    // A query-side SNP beside an otherwise shared insertion can leave both
    // graph signatures evaluable but disjoint. Legacy still classifies that
    // narrow shape as hapfail (rather than hap:mismatch), so its unmatched SNP
    // rows retain BK=`.`. The counterpart predicate is deliberately limited to
    // a same-anchor Insert+Subst conflict with a matching truth insertion.
    //
    // The suppression only holds while the Insert+Subst conflict is still
    // *live* — i.e. at least one record at a conflict position remains
    // unmatched after `exact_match_pairs`. When both members of the conflict
    // were exact-matched (drained to TP), the conflict no longer explains the
    // disjoint signatures: the genuine mismatch lies at another position, and
    // legacy reaches `hap:mismatch` → BK=lm there. Broad germline data hits
    // this via a compound-het insertion pair (truth split, query merged into a
    // hetalt) sitting one anchor away from an already-matched Insert+Subst
    // pair; without this gate those leftover FN/FP rows lost their BK=lm.
    let conflict_positions_still_unmatched = {
        let query_remaining_positions: BTreeSet<usize> =
            query_remaining.iter().map(|v| v.key.pos).collect();
        let (insert_positions, subst_positions) = insert_subst_positions(&cluster.query);
        insert_positions
            .intersection(&subst_positions)
            .any(|pos| query_remaining_positions.contains(pos))
    };
    let shared_insert_conflict = truth_sig.is_some()
        && query_sig.is_some()
        && matches!(insert_conflict_counterpart, Some(true))
        && conflict_positions_still_unmatched
        && cluster.query.iter().any(|variant| {
            // The hapfail shape needs a genuine hetalt insertion aggregate
            // (two DISTINCT alts, e.g. `TTGG,TTGGG`). A duplicate-alt
            // aggregate (`ACCCT,ACCCT`) is a persisted HOMOZYGOUS insertion,
            // not a reciprocal het pair; legacy still reaches hap:mismatch on
            // those blocks (chr1:150042104 HG003 → BK=lm), so exclude them.
            let alts: Vec<&str> = variant.key.alt_allele.split(',').collect();
            let distinct: BTreeSet<&str> = alts.iter().copied().collect();
            distinct.len() >= 2
                && distinct
                    .iter()
                    .any(|alt| alt.len() > variant.key.ref_allele.len())
        });
    let covered_multi_conflict_mismatch = graph_failed
        && truth_remaining.is_empty()
        && legacy_covered_multi_conflict_mismatch(&cluster.query);
    // A failed graph outranks every other signal. Legacy's `finish_block`
    // arms `hap_fail` before entering the `try` and clears it only when
    // `DiploidCompare` returns match or mismatch, so any throw out of
    // `makeGraph`/`enumeratePaths`/`setRegion` leaves the block at
    // `ctype=hapfail:*`, and quantify only promotes `BK=lm` on
    // `ctype=hap:mismatch`. A block rust could not enumerate therefore
    // cannot carry a local mismatch, whatever the shape heuristics say.
    let hap_mismatch = if graph_failed || shared_insert_conflict {
        false
    } else if overlapping_deletion_mismatch || covered_multi_conflict_mismatch {
        true
    } else if allow_haplotype_match && truth_sig.is_some() && query_sig.is_some() {
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
        matches!(insert_conflict_counterpart, Some(false))
    } else {
        false
    };

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
        (XcmpCtype::SimpleMismatch, false)
    } else if graph_failed || shared_insert_conflict {
        (XcmpCtype::HapfailMismatch, false)
    } else {
        match (&truth_sig, &query_sig) {
            (Some(_), Some(_)) if is_match => (XcmpCtype::HapMatch, true),
            (Some(_), Some(_)) => (XcmpCtype::HapMismatch, false),
            _ if hap_mismatch => (XcmpCtype::HapMismatch, false),
            _ => (XcmpCtype::HapfailMismatch, false),
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
    // This graph-level verdict is independent of the linear haplotype
    // signature above. A shared insertion that reaches a later het-alt
    // aggregate is one such case: the linear signatures overlap, while
    // legacy xcmp still marks the outside-CONF aggregate block as a local
    // mismatch. Apply the correction on both match and mismatch paths.
    let unknown_aggregate_local_mismatch =
        legacy_unknown_aggregate_local_mismatch(cluster, &region_state, hap_mismatch);
    if let Some(local_mismatch) = unknown_aggregate_local_mismatch {
        for row in &mut rows[output_start..] {
            if local_mismatch && row.record.samples_contain(":UNK:.:") {
                row.record
                    .replace_sample_fragment(":UNK:.:", ":UNK:lm:")
                    .expect("comparison decision edits preserve valid records");
            } else if !local_mismatch && row.record.samples_contain(":UNK:lm:") {
                row.record
                    .replace_sample_fragment(":UNK:lm:", ":UNK:.:")
                    .expect("comparison decision edits preserve valid records");
            }
        }
    }
    correct_long_aggregate_block_kind(&mut rows[output_start..]);
    if region_state.any_conf {
        correct_reaching_insertion_aggregate_block_kind(&mut rows[output_start..]);
    }
    if !legacy_hap_promotions.is_empty() {
        for row in &mut rows[exact_match_post_count..] {
            if row_matches_variant_key(row, &legacy_hap_promotions)
                && row.record.samples_contain(":FN:am:")
                && row.record.samples_contain(":FP:am:")
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
        xcmp_ctype = XcmpCtype::HapMatch;
        xcmp_hap_match = true;
    }
    // A residual lm/am allele makes the entire outside-CONF exact-pair part
    // of the same local mismatch. Conversely, exact-only preprocessed blocks
    // retain the missing-kind verdict for byte-identical unphased records.
    let residual_is_unreconciled = rows[exact_match_post_count..]
        .iter()
        .any(annotated_row_has_unreconciled_allele);
    if should_propagate_unreconciled_exact_rows(
        xcmp_ctype,
        residual_is_unreconciled,
        covered_multi_conflict_mismatch,
    ) {
        for row in &mut rows[exact_match_pre_count..exact_match_post_count] {
            if row.record.samples_contain(":UNK:.:") {
                row.record
                    .replace_sample_fragment(":UNK:.:", ":UNK:lm:")
                    .expect("comparison decision edits preserve valid records");
            }
        }
    } else if !config.no_hc {
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
        xcmp_ctype = XcmpCtype::HapMismatch;
        xcmp_hap_match = false;
    }
    // pre.py emits unphased, decomposed truth records. When a decomposed SNP
    // and indel share an anchor and the SNP has an exact query counterpart,
    // legacy's VariantLocationAggregator orders the SNP record first. The
    // phased non-preprocessed stream follows the ordinary decision/type sort,
    // so gate this on the all-unphased same-anchor truth shape.
    apply_legacy_combined_before_truth_only_order(&mut rows[output_start..]);
    // Apply the mixed primitive ordering after the broader combined-row
    // rule. Both the SNP and INDEL can be combined TP rows, in which case
    // that rule assigns both rank zero and would otherwise fall back to ALT
    // lexical order (`AT` before `T`).
    apply_legacy_preprocessed_snp_first_order(&mut rows[output_start..], cluster);
    apply_legacy_halfcall_order(&mut rows[output_start..]);
    apply_legacy_halfcall_block_kind(&mut rows[output_start..]);
    set_xcmp_context(&mut rows[output_start..], xcmp_ctype, xcmp_hap_match);
    Ok(())
}

pub(super) fn should_propagate_unreconciled_exact_rows(
    xcmp_ctype: XcmpCtype,
    residual_is_unreconciled: bool,
    genuine_drain_mismatch: bool,
) -> bool {
    residual_is_unreconciled && (xcmp_ctype == XcmpCtype::HapMismatch || genuine_drain_mismatch)
}

pub(super) fn apply_legacy_preprocessed_snp_first_order(
    rows: &mut [AnnotatedRow],
    cluster: &Cluster,
) {
    let snp_first_positions = legacy_preprocessed_snp_first_positions(cluster);
    for row in rows {
        if snp_first_positions.contains(&row.sort_key.pos) {
            row.sort_key.side_rank = usize::from(!annotated_row_is_snp(row));
        }
    }
}

pub(super) fn apply_legacy_combined_before_truth_only_order(rows: &mut [AnnotatedRow]) {
    let combined_positions = rows
        .iter()
        .filter_map(|row| {
            let samples = &row.record.raw().samples;
            (samples.len() >= 2
                && !samples[0].contains("NOCALL:nocall")
                && !samples[1].contains("NOCALL:nocall"))
            .then_some(row.sort_key.pos)
        })
        .collect::<BTreeSet<_>>();
    for row in rows {
        if !combined_positions.contains(&row.sort_key.pos) {
            continue;
        }
        let raw = row.record.raw();
        if raw.alt_allele == "." {
            continue;
        }
        if raw.samples.len() >= 2
            && !raw.samples[0].contains("NOCALL:nocall")
            && !raw.samples[1].contains("NOCALL:nocall")
        {
            row.sort_key.side_rank = 0;
        } else if raw
            .samples
            .get(1)
            .is_some_and(|sample| sample.contains("NOCALL:nocall"))
        {
            row.sort_key.side_rank = 1;
        }
    }
}
pub(super) fn apply_legacy_halfcall_order(rows: &mut [AnnotatedRow]) {
    let halfcall_positions = rows
        .iter()
        .filter_map(|row| {
            let raw = row.record.raw();
            (raw.alt_allele == "."
                && raw
                    .samples
                    .iter()
                    .any(|sample| sample.contains(":UNK:halfcall:")))
            .then_some(raw.pos)
        })
        .collect::<BTreeSet<_>>();
    for row in rows {
        if !halfcall_positions.contains(&row.sort_key.pos) {
            continue;
        }
        let raw = row.record.raw();
        if raw.alt_allele == "." {
            // Legacy emits a matched/combined allele first, then its
            // spanning-deletion halfcall companion.
            row.sort_key.side_rank = 2;
        } else if raw
            .samples
            .get(1)
            .is_some_and(|sample| sample.contains("NOCALL:nocall"))
        {
            // When both rows are truth-only, the halfcall precedes the
            // ordinary FN/UNK allele instead.
            row.sort_key.side_rank = 3;
        } else if raw
            .samples
            .first()
            .is_some_and(|sample| sample.contains("NOCALL:nocall"))
        {
            // With all three grains present at one locus, legacy orders the
            // halfcall, then the truth-only ordinary allele, then the
            // query-only residual allele.
            row.sort_key.side_rank = 4;
        }
    }
}

/// Legacy hapcmp reports a local mismatch when a genotype-discordant
/// deletion pair is overlapped by another truth deletion. The linear event
/// enumerator can drain this shape and surface it as hapfail instead, which
/// loses both `am` on the paired allele and block-wide `lm` on companions.
pub(super) fn legacy_overlapping_deletion_mismatch(cluster: &Cluster) -> bool {
    cluster.truth.iter().any(|paired_truth| {
        if paired_truth.primary_type() != "INDEL" {
            return false;
        }
        let Some(paired_query) = cluster
            .query
            .iter()
            .find(|query| query.key == paired_truth.key)
        else {
            return false;
        };
        if equivalent_gt(&paired_truth.gt, &paired_query.gt) {
            return false;
        }
        cluster.truth.iter().any(|overlapping_truth| {
            overlapping_truth.key != paired_truth.key
                && overlapping_truth.primary_type() == "INDEL"
                && overlapping_truth.key.pos < paired_truth.key.pos
                && overlapping_truth.end_pos() >= paired_truth.end_pos()
        })
    })
}

/// Legacy's spanning-deletion halfcall keeps the local-mismatch block kind
/// of an unmatched truth sibling emitted at the same position.  This is
/// observable when a matched deletion covers both a `*` halfcall and a
/// separate unmatched SNP: the SNP is `FN/lm` and the halfcall is `N/lm`,
/// even though the halfcall has no query-side record of its own.
pub(super) fn apply_legacy_halfcall_block_kind(rows: &mut [AnnotatedRow]) {
    // BK is a block verdict. A halfcall inherits local mismatch from any
    // truth allele in the comparison block, not merely a sibling at the
    // same position. Long spanning deletions in GIAB expose this across
    // hundreds of bases and multiple `*` companions.
    let block_has_local_mismatch = rows.iter().any(|row| {
        row.record
            .raw()
            .samples
            .first()
            .is_some_and(|sample| sample.contains(":lm:"))
    });
    if !block_has_local_mismatch {
        return;
    }
    for row in rows {
        let raw = row.record.raw();
        if raw.alt_allele != "." {
            continue;
        }
        if raw
            .samples
            .first()
            .is_some_and(|sample| sample.contains(":N:.:.:UNK:halfcall:"))
        {
            row.record
                .replace_sample_fragment(":N:.:.:UNK:halfcall:", ":N:lm:.:UNK:halfcall:")
                .expect("comparison decision edits preserve valid records");
        } else if raw
            .samples
            .first()
            .is_some_and(|sample| sample.contains(":UNK:.:.:UNK:halfcall:"))
        {
            row.record
                .replace_sample_fragment(":UNK:.:.:UNK:halfcall:", ":UNK:lm:.:UNK:halfcall:")
                .expect("comparison decision edits preserve valid records");
        }
    }
}

pub(super) fn legacy_unknown_aggregate_local_mismatch(
    cluster: &Cluster,
    region_state: &RegionState,
    rust_hap_mismatch: bool,
) -> Option<bool> {
    if region_state.any_conf {
        return None;
    }
    cluster.query.iter().find_map(|query| {
        if query.gt != "2/1" || !query.key.alt_allele.contains(',') {
            return None;
        }
        let alts = query.key.alt_allele.split(',').collect::<Vec<_>>();
        if alts.len() != 2 {
            return None;
        }
        let truth_at_locus = cluster
            .truth
            .iter()
            .filter(|truth| {
                truth.key.pos == query.key.pos
                    && truth.key.ref_allele == query.key.ref_allele
                    && !truth.key.alt_allele.contains(',')
                    && alts.contains(&truth.key.alt_allele.as_str())
            })
            .collect::<Vec<_>>();
        if truth_at_locus.len() != 2 {
            return None;
        }

        if rust_hap_mismatch {
            // A GT=2/1 insertion aggregate immediately followed by a shared
            // hom-alt deletion can reconstruct the same two legacy graph
            // paths even though the linear event model reports a mismatch.
            // The deletion must consume more reference than either aggregate
            // ALT contributes, making this a graph-splice correction rather
            // than a generic nearby-variant heuristic.
            let max_alt_len = alts.iter().map(|alt| alt.len()).max().unwrap_or(0);
            let followed_by_shared_covering_deletion = cluster.truth.iter().any(|truth| {
                truth.key.pos == query.key.pos + 1
                    && truth.key.ref_allele.len() > max_alt_len
                    && truth
                        .key
                        .alt_allele
                        .split(',')
                        .all(|alt| alt.len() < truth.key.ref_allele.len())
                    && cluster.query.iter().any(|candidate| {
                        candidate.key == truth.key && equivalent_gt(&candidate.gt, &truth.gt)
                    })
            });
            return followed_by_shared_covering_deletion.then_some(false);
        }
        None
    })
}

pub(super) fn correct_long_aggregate_block_kind(rows: &mut [AnnotatedRow]) {
    for row in rows {
        let raw = row.record.raw();
        let is_long_aggregate = raw.alt_allele.contains(',')
            && raw.alt_allele.split(',').all(|alt| alt.len() > 512)
            && raw
                .samples
                .iter()
                .any(|sample| sample.starts_with("2/1:UNK:lm:"));
        if is_long_aggregate {
            row.record
                .replace_sample_fragment(":UNK:lm:", ":UNK:.:")
                .expect("comparison decision edits preserve valid records");
        }
    }
}

pub(super) fn correct_reaching_insertion_aggregate_block_kind(rows: &mut [AnnotatedRow]) {
    let aggregate_positions = rows
        .iter()
        .filter_map(|row| {
            let raw = row.record.raw();
            (raw.alt_allele.contains(',')
                && raw
                    .samples
                    .iter()
                    .any(|sample| sample.starts_with("2/1:UNK:")))
            .then_some(raw.pos)
        })
        .collect::<Vec<_>>();
    let legacy_local_mismatch = aggregate_positions.iter().any(|aggregate_pos| {
        rows.iter().any(|row| {
            let raw = row.record.raw();
            let max_alt_len = raw.alt_allele.split(',').map(str::len).max().unwrap_or(0);
            raw.pos < *aggregate_pos
                && max_alt_len > raw.ref_allele.len()
                && (16..=64).contains(&max_alt_len)
                && raw.pos.saturating_add(max_alt_len).saturating_add(2) >= *aggregate_pos
        })
    });
    if legacy_local_mismatch {
        for row in rows {
            if row.record.samples_contain(":UNK:.:") {
                row.record
                    .replace_sample_fragment(":UNK:.:", ":UNK:lm:")
                    .expect("comparison decision edits preserve valid records");
            }
        }
    }
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
        if row_matches_variant_key(row, keys) && row.record.samples_contain(":UNK:lm:") {
            row.record
                .replace_sample_fragment(":UNK:lm:", ":UNK:.:")
                .expect("comparison decision edits preserve valid records");
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
                // Decomposition products of a single truth variant carry the
                // identical genotype string, and only those get SNP-first
                // ordering from legacy's aggregator. Independent colocated
                // records on opposite haplotypes (e.g. an insertion `1/0` and
                // a SNP `0/1`) are two separate variants; legacy keeps their
                // prep-stream order, which the ordinary ALT-lexical sort
                // already reproduces. Gate on shared GT so this rule only
                // fires for the genuine decomposition shape.
                && truth_at_pos
                    .first()
                    .is_some_and(|first| truth_at_pos.iter().all(|v| v.gt == first.gt))
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

pub(super) fn set_xcmp_context(rows: &mut [AnnotatedRow], ctype: XcmpCtype, hap_match: bool) {
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

/// Partition variant positions into Insert (any alt longer than ref) and
/// Subst (all alts no longer than ref, i.e. SNP/MNP/deletion) buckets. A
/// position present in both sets carries a same-anchor Insert+Subst conflict.
/// Shared by the hapfail-suppression predicate and its still-unmatched gate so
/// both classify identically.
fn insert_subst_positions(variants: &[Variant]) -> (BTreeSet<usize>, BTreeSet<usize>) {
    let mut insert_positions: BTreeSet<usize> = BTreeSet::new();
    let mut subst_positions: BTreeSet<usize> = BTreeSet::new();
    for v in variants {
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
    (insert_positions, subst_positions)
}

pub(super) fn query_insert_conflict_has_truth_counterpart(
    query_variants: &[Variant],
    truth_variants: &[Variant],
    truth_remaining: &[Variant],
) -> Option<bool> {
    // Classify each query variant position as Insert (any alt longer than ref)
    // or Subst (all alts no longer than ref). The graph drain path groups a
    // same-anchor deletion with substitutions; the truth-counterpart check
    // below distinguishes the reciprocal aggregate that remains hapfail.
    let (insert_positions, subst_positions) = insert_subst_positions(query_variants);
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
            // A truth insertion can be one selected member of the query's
            // reciprocal insertion/deletion aggregate. Accept that direction
            // as a counterpart. The reverse remains strict: query "TAA" does
            // not match truth "TAA,TA" because the query omits a truth allele.
            let q_alts: std::collections::BTreeSet<&str> = qv.key.alt_allele.split(',').collect();
            if truth_variants.iter().any(|tv| {
                if tv.key.pos != qv.key.pos || tv.key.ref_allele != qv.key.ref_allele {
                    return false;
                }
                let t_alts: std::collections::BTreeSet<&str> =
                    tv.key.alt_allele.split(',').collect();
                !t_alts.is_empty() && t_alts.is_subset(&q_alts)
            }) {
                return Some(true); // expected hapfail: truth has identical insert allele set
            }
        }
    }
    Some(false) // genuine mismatch: no truth counterpart for any conflict-pos insert
}

/// Legacy graph comparison can complete a phased query block when one
/// deletion spans multiple later Insert+Subst conflicts. The simpler Rust
/// graph can drain that block instead. A single conflict is not sufficient:
/// the public HG001 corpus contains many such hapfail blocks with BK=`.`.
pub(super) fn legacy_covered_multi_conflict_mismatch(query_variants: &[Variant]) -> bool {
    let mut insert_positions = BTreeSet::new();
    let mut substitution_positions = BTreeSet::new();
    for variant in query_variants {
        let ref_len = variant.key.ref_allele.len();
        let alt_lens = variant
            .key
            .alt_allele
            .split(',')
            .map(str::len)
            .collect::<Vec<_>>();
        if alt_lens.iter().any(|alt_len| *alt_len > ref_len) {
            insert_positions.insert(variant.key.pos);
        }
        if alt_lens.iter().all(|alt_len| *alt_len == ref_len) {
            substitution_positions.insert(variant.key.pos);
        }
    }
    let conflicts = insert_positions
        .intersection(&substitution_positions)
        .copied()
        .collect::<Vec<_>>();
    if conflicts.len() < 2 {
        return false;
    }
    query_variants.iter().any(|variant| {
        let ref_len = variant.key.ref_allele.len();
        let max_alt_len = variant
            .key
            .alt_allele
            .split(',')
            .map(str::len)
            .max()
            .unwrap_or(0);
        ref_len > max_alt_len
            && !variant.gt.split(['/', '|']).any(|allele| allele == "0")
            && conflicts.iter().all(|position| {
                variant.key.pos < *position && *position < variant.key.pos + ref_len
            })
    })
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
    // A query duplicate-alt aggregate (`X,X`, GT 2/1) is the persisted form of
    // two source het records that both produced the same edit X — one on each
    // haplotype. Pairing it to a truth homalt claims both query haplotypes are
    // exactly X. That claim is false when a SEPARATE query record shares this
    // anchor with a GT-selected nonref allele: one of those haplotypes then
    // also carries the neighbour's edit, so it is not X alone. Legacy defers
    // such a block to hap-compare, which finds the two sides' haplotypes differ
    // and leaves the truth FN with both query alleles FP. chr1:150042104
    // (HG003 DeepVariant): truth `A>ACCCT 1/1`, query `A>ACCCT,ACCCT 2/1` +
    // `A>T 1/0` — the SNP on one hap makes it `TCCCT`, so the aggregate is not
    // a clean homozygous insertion and must not drain as a TP here.
    if query_duplicate_alt_aggregate_has_conflicting_neighbor(query, cluster_query) {
        return false;
    }
    // `gttype` equality is implied by the selected multi-set equality —
    // homalt `1|1` selects `[X, X]` vs het `0|1` selects `[X]`; these
    // differ as multi-sets even when the nonref-set matches.
    true
}

/// True when `query` is a duplicate-alt aggregate (an ALT column listing the
/// same allele more than once, e.g. `ACCCT,ACCCT`) AND another query record at
/// the same anchor position carries a GT-selected nonref allele. In that shape
/// the aggregate's implied homozygosity is contradicted by the neighbour's edit
/// sharing a haplotype, so the pair is not a clean simple-compare match and
/// must fall through to the block-level haplotype comparison.
fn query_duplicate_alt_aggregate_has_conflicting_neighbor(
    query: &Variant,
    cluster_query: &[Variant],
) -> bool {
    let alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    if alts.len() < 2 {
        return false;
    }
    let distinct: BTreeSet<&str> = alts.iter().copied().collect();
    if distinct.len() == alts.len() {
        // No repeated alt — not a duplicate-alt aggregate.
        return false;
    }
    cluster_query.iter().any(|neighbour| {
        neighbour.key.pos == query.key.pos
            && neighbour.key != query.key
            && !selected_alt_sequences(neighbour).is_empty()
    })
}

/// Legacy renders a query duplicate-alt DELETION aggregate (`X,X` selecting
/// both indices, X shorter than ref) that faces a truth homalt for the same
/// deletion as TWO het rows: one deletion copy allele-matches the truth
/// (combined FN:am/FP:am row) and the other is excess (FP:lm). Its graph
/// aligns one homozygous truth-deletion haplotype to the clean query
/// haplotype and leaves the copy sharing a haplotype with a conflicting
/// neighbour unmatched. Duplicate-alt INSERTION aggregates stay collapsed as
/// one homalt FP:lm row (chr1:150042104, chr3:30297039, chr7:68324852,
/// chr13:74256688 all measured that way), so this fires for deletions only.
///
/// Rewriting the aggregate into `[X 0/1, X 1/0]` lets the existing pairing
/// loop bind the `0/1` copy to the truth homalt as `am` (its selected
/// multiset `[X]` differs from the truth's `[X, X]`, so the genotype-mismatch
/// gate passes) while the `1/0` copy falls through to the unpaired-query
/// emitter and picks up `lm` from `bk_for_row`. The `0/1` copy is listed
/// first so the pairing loop's first-match `break` binds it, matching
/// legacy's `0/1` on the am row and `1/0` on the standalone FP.
/// chr6:91567856 (HG003 DeepVariant): truth `TC>T 1/1`, query `TC>T,T 1/2`
/// beside a `T>A 1/0` SNP.
pub(super) fn split_matched_deletion_aggregate(
    query: &Variant,
    truth: &[Variant],
) -> Option<[Variant; 2]> {
    let alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    let [left, right] = alts.as_slice() else {
        return None;
    };
    if left != right {
        return None;
    }
    let alt = *left;
    // Deletion only — an insertion aggregate stays collapsed as one homalt row.
    if alt.len() >= query.key.ref_allele.len() {
        return None;
    }
    let selected: BTreeSet<usize> = parse_gt_alleles(&query.gt).into_iter().collect();
    if selected != BTreeSet::from([1, 2]) {
        return None;
    }
    // A truth homalt for exactly this deletion must exist in the block — that
    // homozygous truth haplotype is what makes one query copy an allele match.
    let has_truth_homalt = truth.iter().any(|t| {
        t.key.chrom == query.key.chrom
            && t.key.pos == query.key.pos
            && t.key.ref_allele == query.key.ref_allele
            && t.key.alt_allele == alt
            && {
                let gt = parse_gt_alleles(&t.gt);
                !gt.is_empty() && gt.iter().all(|allele| *allele == 1)
            }
    });
    if !has_truth_homalt {
        return None;
    }
    let make = |gt: &str| Variant {
        key: VariantKey {
            chrom: query.key.chrom.clone(),
            pos: query.key.pos,
            ref_allele: query.key.ref_allele.clone(),
            alt_allele: alt.to_string(),
        },
        qual: query.qual.clone(),
        filter: query.filter.clone(),
        gt: gt.to_string(),
    };
    Some([make("0/1"), make("1/0")])
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
    // A `2/1` multi-allelic query is the persisted output of legacy's
    // VariantLocationAggregator. It represents both source calls as one
    // final hetalt record and must survive comparison at that row grain.
    // Reapplying the subset matcher consumes one allele into a truth row
    // and emits the other as a residual, turning one real query record into
    // two counted primitives (for example chr1:963700 in test_full).
    // Raw pre-aggregation subset fixtures retain their source `1/2` GT and
    // continue through the legacy combined-row-plus-residual path below.
    if query.gt == "2/1" && query.key.alt_allele.contains(',') {
        return false;
    }
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
            let projected_query = legacy_duplicate_alt_query_output_projection(query);
            let query_for_output = projected_query.as_ref();
            if is_subset {
                // Truth-subset path: truth's alts ⊊ query's selected.
                // Emit at truth's representation; remap query GT so
                // alleles missing from truth's column collapse to ref.
                let canonicalized = Variant {
                    key: VariantKey {
                        chrom: query_for_output.key.chrom.clone(),
                        pos: query_for_output.key.pos,
                        ref_allele: query_for_output.key.ref_allele.clone(),
                        alt_allele: truth.key.alt_allele.clone(),
                    },
                    qual: query_for_output.qual.clone(),
                    filter: query_for_output.filter.clone(),
                    gt: remap_query_gt_subset(truth, query_for_output),
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
                let query_for_row = if is_multi_allelic(query_for_output) {
                    // For a phased query (GT has '|'), legacy mirrors truth's
                    // haplotype assignment unphased: truth 1|2 → query 1/2.
                    // canonical_hetalt_gt (alphabetical "later/earlier") is wrong
                    // for phased cases like CACACACAT,CAT where it gives 2|1.
                    // For an unphased query (GT has '/'), legacy uses canonical
                    // alphabetical ordering (e.g. C,G unphased 1/2 → 2/1).
                    let gt = if query_for_output.gt.contains('|') {
                        truth.gt.replace('|', "/")
                    } else {
                        canonical_hetalt_gt(&truth.key.alt_allele, query_for_output)
                    };
                    Variant {
                        key: query_for_output.key.clone(),
                        qual: query_for_output.qual.clone(),
                        filter: query_for_output.filter.clone(),
                        gt,
                    }
                } else {
                    query_for_output.clone()
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
                        chrom: query_for_output.key.chrom.clone(),
                        pos: query_for_output.key.pos,
                        ref_allele: query_for_output.key.ref_allele.clone(),
                        alt_allele: truth.key.alt_allele.clone(),
                    },
                    qual: query_for_output.qual.clone(),
                    filter: query_for_output.filter.clone(),
                    gt: canonical_hetalt_gt(&truth.key.alt_allele, query_for_output),
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
        if region_state.truth_is_conf(truth)
            && (truth.primary_type() != "UNK"
                || halfcall_is_covered_by_matched_deletion(truth, full_cluster))
        {
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

    // A haplotype-matched superlocus can still contain a same-allele record
    // whose per-record genotypes differ (for example truth 0|1 vs query 1/1).
    // Legacy's VariantReader unifies that record before block comparison and
    // therefore emits one combined row even when other records share the
    // superlocus. Keep those pairs together on the hap-match path; the old
    // implementation only combined the degenerate one-truth/one-query case.
    let mut paired_truth = BTreeSet::new();
    let mut paired_query = BTreeSet::new();
    let mut paired_records = Vec::new();
    for (truth_index, truth) in cluster.truth.iter().enumerate() {
        for (query_index, query) in cluster.query.iter().enumerate() {
            if paired_query.contains(&query_index) {
                continue;
            }
            if query_matches_truth_allele_set(query, truth)
                && selected_alt_sequences(truth) != selected_alt_sequences(query)
                && !query_primitive_splits(
                    truth,
                    reference,
                    cluster.start,
                    &cluster.truth,
                    &cluster.query,
                )
                && !query_primitive_splits(
                    query,
                    reference,
                    cluster.start,
                    &cluster.truth,
                    &cluster.query,
                )
            {
                paired_truth.insert(truth_index);
                paired_query.insert(query_index);
                paired_records.push((truth_index, query_index));
                break;
            }
        }
    }
    for (truth_index, query_index) in &paired_records {
        let truth = &cluster.truth[*truth_index];
        let query = &cluster.query[*query_index];
        let regions = region_state.row_tags(Some(truth), Some(query));
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
                .filter(|v| v.is_finite() && *v >= 0.0)
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
    for (truth_index, truth) in cluster.truth.iter().enumerate() {
        if paired_truth.contains(&truth_index) {
            continue;
        }
        if region_state.truth_is_conf(truth) {
            if truth.primary_type() == "UNK"
                && !halfcall_is_covered_by_matched_deletion(truth, full_cluster)
            {
                rows.push(fn_row(
                    truth,
                    reference,
                    cluster.start,
                    &region_state.row_tags(Some(truth), None),
                    ".",
                ));
            } else {
                rows.push(tp_single_side_row(
                    truth,
                    reference,
                    cluster.start,
                    &region_state.row_tags(Some(truth), None),
                    Side::Truth,
                    shared_qq,
                    &truth_cluster_filter,
                ));
            }
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
    for (query_index, query) in cluster.query.iter().enumerate() {
        if paired_query.contains(&query_index) {
            continue;
        }
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
                let bk = bk_for_row(&primitive, &full_cluster.truth, false);
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

pub(super) fn halfcall_is_covered_by_matched_deletion(
    halfcall: &Variant,
    cluster: &Cluster,
) -> bool {
    if halfcall.primary_type() != "UNK" || !halfcall.gt.contains('.') {
        return false;
    }
    let requires_haplotype_reconciliation = cluster.query.iter().any(|query| {
        (query.gt == "2/1"
            && query.key.alt_allele.contains(',')
            && !cluster.truth.iter().any(|truth| truth.key == query.key))
            || cluster
                .truth
                .iter()
                .any(|truth| truth.key == query.key && !equivalent_gt(&truth.gt, &query.gt))
    });
    cluster.truth.iter().any(|truth| {
        let covering_deletion = truth.primary_type() == "INDEL"
            && truth
                .key
                .alt_allele
                .split(',')
                .any(|alt| alt.len() < truth.key.ref_allele.len())
            && truth.key.pos <= halfcall.key.pos
            && halfcall.key.pos <= truth.end_pos();
        if !covering_deletion {
            return false;
        }
        // In a hap:match block, a covering truth deletion with no byte-equal
        // query counterpart was necessarily rescued by haplotype
        // reconciliation, so its `*` companion is TP/gm as well. For a
        // simple exact deletion pair legacy leaves the companion neutral;
        // it is promoted only when another allele forces reconciliation.
        let has_exact_query_deletion = cluster.query.iter().any(|query| query.key == truth.key);
        !has_exact_query_deletion || requires_haplotype_reconciliation
    })
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
    // Decompose query duplicate-alt DELETION aggregates that face a truth
    // homalt into two het copies (`0/1` am-copy first, `1/0` lm-copy) so the
    // pairing loop below binds one as `am` and the other emits as an `lm` FP,
    // mirroring legacy at chr6:91567856. Insertion aggregates are untouched.
    let decomposed_cluster;
    let cluster = {
        let mut changed = false;
        let mut new_query = Vec::with_capacity(cluster.query.len());
        for query in &cluster.query {
            match split_matched_deletion_aggregate(query, &cluster.truth) {
                Some([am_copy, lm_copy]) => {
                    new_query.push(am_copy);
                    new_query.push(lm_copy);
                    changed = true;
                }
                None => new_query.push(query.clone()),
            }
        }
        if changed {
            decomposed_cluster = Cluster {
                chrom: cluster.chrom.clone(),
                start: cluster.start,
                end: cluster.end,
                truth: cluster.truth.clone(),
                query: new_query,
            };
            &decomposed_cluster
        } else {
            cluster
        }
    };
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
                && !mixed_type_same_locus_keeps_indel_rows_separate(cluster, truth, query)
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
        let projected_query = legacy_duplicate_alt_query_output_projection(query);
        let query_for_output = projected_query.as_ref();
        // Canonicalise unphased hetalt query GT to legacy's
        // `<later>/<earlier>` ordering when the row is hetalt and the
        // declared alts have a non-canonical order. SNP cases like
        // chr21:9922359 (`T→A,C` query GT `1/2`) emit `2/1` per
        // VariantLocationAggregator's MAX_GT=2 rule. Phased / hom /
        // het-with-ref pass through unchanged.
        let query_gt = if is_distinct_hetalt(&query_for_output.gt) {
            canonical_hetalt_gt(&truth.key.alt_allele, query_for_output)
        } else {
            remap_query_gt_subset(truth, query_for_output)
        };
        let query_for_row = Variant {
            key: truth.key.clone(),
            qual: query_for_output.qual.clone(),
            filter: query_for_output.filter.clone(),
            gt: query_gt,
        };
        rows.push(fn_fp_combined_row(
            truth,
            &query_for_row,
            reference,
            cluster.start,
            &region_state.row_tags(Some(truth), Some(query)),
            &BenchmarkDecision {
                truth_bd,
                query_bd,
                bk,
            },
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
            let primitive_bk = bk_for_row(&primitive, &full_cluster.truth, hap_mismatch);
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

pub(super) fn mixed_type_same_locus_keeps_indel_rows_separate(
    cluster: &Cluster,
    truth: &Variant,
    query: &Variant,
) -> bool {
    if truth.primary_type() != "INDEL" || query.primary_type() != "INDEL" {
        return false;
    }
    cluster.truth.iter().any(|snp_truth| {
        snp_truth.key.pos == truth.key.pos
            && snp_truth.primary_type() == "SNP"
            && cluster.query.iter().any(|snp_query| {
                snp_query.key.pos == query.key.pos
                    && snp_query.primary_type() == "SNP"
                    && query_matches_truth_allele_set(snp_query, snp_truth)
                    && selected_alt_sequences(snp_truth) != selected_alt_sequences(snp_query)
            })
    })
}

/// Map the per-row `BK` (block kind) tag to the FP classification used by
/// the summary / extended / ROC `FP.gt` / `FP.al` columns. Mirrors legacy
/// hap.py's quantify aggregation: FP rows tagged `BK=am` (allele match,
/// genotype mismatch) bump `FP.gt`; FP rows tagged `BK=lm` (locus match,
/// allele mismatch) bump `FP.al`; everything else (novel FPs with `BK=.`)
/// stays unclassified and is omitted from both columns.
pub(super) fn fp_class_from_bk(bk: &str) -> Option<FpClass> {
    match bk {
        "am" => Some(FpClass::Gt),
        "lm" => Some(FpClass::Al),
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
    gvcf2bed_padding_iter(truth.iter().map(Ok::<_, anyhow::Error>), target_bed)
        .expect("in-memory variants are infallible")
}

pub(super) fn gvcf2bed_padding_iter<I, V>(
    truth: I,
    target_bed: Option<&[Interval]>,
) -> Result<Vec<Interval>>
where
    I: IntoIterator<Item = Result<V>>,
    V: Borrow<Variant>,
{
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

    let mut out = Vec::new();
    let mut active: Option<ActiveInterval> = None;

    for variant in truth {
        let variant = variant?;
        let variant = variant.borrow();
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
    Ok(out)
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
    // A selected spanning-deletion (`ALT=*`) is loaded as a legacy
    // half-call (`ALT=.`, GT containing one missing allele). QuantifyRegions
    // classifies that record by its anchor, even though the source REF span
    // is retained so the emitted END field can describe the deletion it
    // overlaps. Treating the entire REF span as the confidence footprint
    // incorrectly turns confident N/TP half-calls into UNK and makes their
    // clusters TS_boundary in test_full.
    if variant.primary_type() == "UNK" && variant.gt.contains('.') {
        return conf_intervals
            .iter()
            .any(|interval| interval.matches(&variant.key.chrom, variant.key.pos));
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
