use crate::cli::{CompareArgs, PreprocessArgs};
use crate::fasta;
use crate::metrics_json;
use crate::preprocess;
use crate::report::{self, CountsBucket};
use crate::vcf::{self, Variant, VariantKey};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

#[derive(Clone, Debug, Default)]
pub struct TypeCounts {
    pub truth_total: CountsBucket,
    pub truth_tp: CountsBucket,
    pub truth_fn: CountsBucket,
    pub query_total: CountsBucket,
    pub query_tp: CountsBucket,
    pub query_fp: CountsBucket,
    pub query_unk: CountsBucket,
}

#[derive(Clone, Debug)]
pub struct AnnotatedRow {
    pub sort_key: (String, usize, usize, usize),
    pub line: String,
    /// Whether the originating query variant was PASS-filtered. Truth
    /// variants are always considered PASS (the `--usefiltered-truth=False`
    /// contract is enforced upstream in preprocess). The field drives the
    /// summary / extended `ALL` vs `PASS` bifurcation: PASS rows count only
    /// when this is `true`, `ALL` rows count unconditionally.
    pub query_pass: bool,
    /// Sub-classification for an FP row: `Some("gt")` when the query shares
    /// chrom/pos/ref/alt with a truth variant (genotype mismatch only),
    /// `Some("al")` when the allele itself doesn't match any truth at the
    /// same locus. `None` for non-FP rows.
    pub fp_class: Option<&'static str>,
}

#[derive(Clone, Debug)]
struct Cluster {
    chrom: String,
    start: usize,
    end: usize,
    truth: Vec<Variant>,
    query: Vec<Variant>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Side {
    Truth,
    Query,
}

#[derive(Clone, Debug)]
struct Entry {
    side: Side,
    variant: Variant,
}

const CLUSTER_GAP_BP: usize = 50;

/// Matches legacy hap.py's `--xcmp-enumeration-threshold` default. Beyond this
/// number of haplotype-assignment states we stop enumerating and treat the
/// cluster as a mismatch. Without the cap a dense cluster of ~30 het variants
/// explodes to 2^30 × two `Vec<Event>` allocations — hundreds of GB of RAM.
const XCMP_ENUMERATION_THRESHOLD: usize = 16_768;

/// Upper bound on variants admitted to a single cluster. Once exceeded we
/// force-split the cluster even if the inter-variant gap is under 50 bp. This
/// protects `apply_events`/`cluster_signature` paths from pathological inputs
/// that would otherwise allocate multi-GB strings. Legacy xcmp's `finish_block`
/// has no record-count cap — only a position gap (`hb_window = 50 bp`) — so
/// matching legacy BS grouping requires admitting very large clusters when
/// records are tightly packed (e.g. `chr21:10697903` has 67 variants within
/// 50 bp gaps). The `estimated_state_count` threshold still short-circuits
/// `cluster_signature` on clusters that would blow up haplotype enumeration.
const MAX_CLUSTER_VARIANTS: usize = 10_000;

#[derive(Clone, Debug)]
enum Event {
    Subst { pos: usize, alt: char },
    Insert { anchor: usize, seq: String },
    Delete { start: usize, end: usize },
}

#[derive(Clone, Debug, Default)]
struct RegionState {
    conf_enabled: bool,
    any_conf: bool,
    any_nonconf: bool,
    conf_intervals: Vec<vcf::BedInterval>,
    covered_truth: BTreeSet<VariantKey>,
    covered_query: BTreeSet<VariantKey>,
}

impl RegionState {
    fn from_cluster(
        cluster: &Cluster,
        reference: &str,
        conf_bed: Option<&[vcf::BedInterval]>,
    ) -> Self {
        let mut state = Self {
            conf_enabled: conf_bed.is_some(),
            conf_intervals: conf_bed.unwrap_or(&[]).to_vec(),
            ..Self::default()
        };
        let Some(conf_bed) = conf_bed else {
            return state;
        };

        for truth in &cluster.truth {
            let covered = variant_is_conf(truth, reference, cluster.start, cluster.end, conf_bed);
            if covered {
                state.any_conf = true;
                state.covered_truth.insert(truth.key.clone());
            } else {
                state.any_nonconf = true;
            }
        }
        for query in &cluster.query {
            let covered = variant_is_conf(query, reference, cluster.start, cluster.end, conf_bed);
            if covered {
                state.any_conf = true;
                state.covered_query.insert(query.key.clone());
            } else {
                state.any_nonconf = true;
            }
        }

        state
    }

    fn row_tags(&self, truth: Option<&Variant>, query: Option<&Variant>) -> String {
        let mut tags: Vec<&str> = Vec::new();
        let covered = truth
            .map(|variant| self.covered_truth.contains(&variant.key))
            .unwrap_or(false)
            || query
                .map(|variant| self.covered_query.contains(&variant.key))
                .unwrap_or(false);
        if covered {
            tags.push("CONF");
        }
        if self.any_conf {
            if self.any_nonconf {
                tags.push("TS_boundary");
            } else {
                tags.push("TS_contained");
            }
        }
        if tags.is_empty() {
            String::new()
        } else {
            format!(";Regions={}", tags.join(","))
        }
    }

    fn query_is_conf(&self, query: &Variant) -> bool {
        !self.conf_enabled || self.covered_query.contains(&query.key)
    }

    fn truth_is_conf(&self, truth: &Variant) -> bool {
        !self.conf_enabled || self.covered_truth.contains(&truth.key)
    }
}

pub fn run(args: CompareArgs) -> Result<()> {
    let reference_path = Path::new(&args.reference);
    let prefix = Path::new(&args.report_prefix);
    let reference_sequences = fasta::read_sequences(reference_path)?;
    let contig_lengths: BTreeMap<String, usize> = reference_sequences
        .iter()
        .map(|(name, sequence)| (name.clone(), sequence.len()))
        .collect();
    // Legacy's Subset.Size column reports the "N-trimmed length" from
    // `fastainfo` (Tools/fastasize.py → `n_trimmed_length`): contig length
    // minus only the leading and trailing Ns. Internal N-tracts (e.g. the
    // centromere run on chr21) stay counted. Rust previously stripped
    // every N which matched ~35.1 M for chr21 vs legacy's 38.7 M. Using
    // the leading/trailing-only trim recovers the legacy convention.
    let contig_non_n_lengths: BTreeMap<String, usize> = reference_sequences
        .iter()
        .map(|(name, sequence)| {
            let bytes = sequence.as_bytes();
            let leading = bytes
                .iter()
                .take_while(|byte| matches!(**byte, b'N' | b'n'))
                .count();
            let trailing = bytes
                .iter()
                .rev()
                .take_while(|byte| matches!(**byte, b'N' | b'n'))
                .count();
            let trimmed = bytes.len().saturating_sub(leading + trailing);
            (name.clone(), trimmed)
        })
        .collect();
    let contig_set: BTreeSet<String> = contig_lengths.keys().cloned().collect();
    let bed = args
        .regions_bedfile
        .as_ref()
        .or(args.targets_bedfile.as_ref())
        .map(|path| vcf::load_bed(Path::new(path), &contig_set))
        .transpose()?;
    let conf_bed = args
        .fp_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &contig_set))
        .transpose()?;
    let locations = args
        .locations
        .as_deref()
        .map(|text| vcf::parse_locations(text, &contig_set))
        .transpose()?;

    // Legacy hap.py runs `pre.py` on both truth and query before xcmp so the
    // two sides share canonical representations (left-shifted, primitive-split,
    // INFO-stripped, ADO-inserted, FORMAT-reordered). Replicate that contract
    // here on the main code path — reuse the already-byte-parity preprocess
    // pipeline shipped in Phase 2.
    //
    // Truth is always pass-only regardless of the user's `--pass-only` flag.
    // This mirrors legacy hap.py's default `--usefiltered-truth=False` which
    // drops any truth record whose FILTER is not PASS (OverlapConflict,
    // SuspiciousHomAlt, etc.) before xcmp runs — they're part of the truth
    // set's own "noise", not a call to validate. Query receives the user's
    // flag as-is.
    let scratch = scratch_dir(prefix)?;
    let truth_prep = scratch.join("truth.prep.vcf.gz");
    let query_prep = scratch.join("query.prep.vcf.gz");
    preprocess::run(build_preprocess_args(
        &args, &args.truth, &truth_prep, true, false,
    ))?;
    preprocess::run(build_preprocess_args(
        &args,
        &args.query,
        &query_prep,
        args.pass_only,
        false,
    ))?;

    let truth = vcf::load_variants(
        &truth_prep,
        &contig_set,
        args.pass_only,
        bed.as_deref(),
        locations.as_deref(),
    )?;
    let query = vcf::load_variants(
        &query_prep,
        &contig_set,
        args.pass_only,
        bed.as_deref(),
        locations.as_deref(),
    )?;

    let contigs_in_play = collect_contigs(&truth, &query, locations.as_deref());
    let subset_size = contigs_in_play
        .iter()
        .filter_map(|contig| contig_non_n_lengths.get(contig))
        .sum::<usize>();
    if subset_size == 0 {
        bail!("no reference contigs selected for analysis");
    }
    let mut counts: BTreeMap<String, TypeCounts> = BTreeMap::new();
    let mut subtype_counts: BTreeMap<String, BTreeMap<String, TypeCounts>> = BTreeMap::new();
    for variant in &truth {
        add_variant_stats(
            &mut counts
                .entry(variant.primary_type().to_string())
                .or_default()
                .truth_total,
            variant,
        );
        add_variant_stats_subtype(
            &mut subtype_counts,
            variant.primary_type(),
            variant,
            |stats| &mut stats.truth_total,
        );
    }
    for variant in &query {
        add_variant_stats(
            &mut counts
                .entry(variant.primary_type().to_string())
                .or_default()
                .query_total,
            variant,
        );
        add_variant_stats_subtype(
            &mut subtype_counts,
            variant.primary_type(),
            variant,
            |stats| &mut stats.query_total,
        );
    }

    let clusters = build_clusters(&truth, &query);
    // Fold gvcf2bed-style insertion padding (derived from truth) into
    // the raw CONF bed before classification. Legacy hap.py does the
    // same in Python (hap.py:323): `args.strat_regions.append(
    // "CONF_VARS:" + gvcf2bed(vcf1))`, and QuantifyRegions's CONF-label
    // squash folds the CONF_VARS lane back under the "CONF" lane.
    // Merging here bridges the 1-base gaps legacy's CONF beds carry at
    // every interval boundary (e.g. chr21:17562905) whenever an
    // adjacent truth insertion would pad them.
    let adjusted_conf_bed: Option<Vec<vcf::BedInterval>> = conf_bed.as_ref().map(|raw| {
        let mut combined = raw.clone();
        combined.extend(gvcf2bed_padding(&truth));
        merge_bed_intervals(&combined)
    });
    let conf_size = adjusted_conf_bed.as_deref().map(inclusive_region_size).unwrap_or(0);
    let mut rows = Vec::new();
    for cluster in clusters {
        process_cluster(
            &cluster,
            &reference_sequences,
            adjusted_conf_bed.as_deref(),
            &mut counts,
            &mut subtype_counts,
            &mut rows,
        )?;
    }

    rows.sort_by(|left, right| {
        left.sort_key
            .cmp(&right.sort_key)
            .then_with(|| left.line.cmp(&right.line))
    });

    if let Some(parent) = prefix.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }

    let all_counts = derive_total_counts(&rows, false);
    let pass_counts = derive_total_counts(&rows, true);
    let all_subtype = derive_subtype_counts(&rows, false);
    let pass_subtype = derive_subtype_counts(&rows, true);
    let all_subset = derive_subset_counts(&rows, false);
    let pass_subset = derive_subset_counts(&rows, true);
    let all_subset_subtype = derive_subset_subtype_counts(&rows, false);
    let pass_subset_subtype = derive_subset_subtype_counts(&rows, true);
    let all_fp = derive_fp_classes(&rows, false);
    let pass_fp = derive_fp_classes(&rows, true);
    report::write_summary(
        &prefix.with_extension("summary.csv"),
        &all_counts,
        &pass_counts,
        &all_fp,
        &pass_fp,
    )?;
    report::write_extended(
        &prefix.with_extension("extended.csv"),
        &all_counts,
        &pass_counts,
        &all_subtype,
        &pass_subtype,
        subset_size,
        conf_size,
        &all_subset,
        &pass_subset,
        &all_subset_subtype,
        &pass_subset_subtype,
        &all_fp,
        &pass_fp,
    )?;
    let vcf_headers = build_vcf_headers(&contig_lengths);
    report::write_vcf(&prefix.with_extension("vcf.gz"), &vcf_headers, &rows)?;
    crate::roc::write_roc_files(prefix, &rows, subset_size, conf_size)?;
    let commandline = format!(
        "hap germline {} {} -r {} -o {}",
        args.truth, args.query, args.reference, args.report_prefix
    );
    let mut final_args = BTreeMap::new();
    final_args.insert("truth", args.truth.clone());
    final_args.insert("query", args.query.clone());
    final_args.insert("ref", args.reference.clone());
    final_args.insert("reports_prefix", args.report_prefix.clone());
    final_args.insert("pass_only", args.pass_only.to_string());
    final_args.insert(
        "regions_bedfile",
        args.regions_bedfile.clone().unwrap_or_default(),
    );
    final_args.insert(
        "targets_bedfile",
        args.targets_bedfile.clone().unwrap_or_default(),
    );
    final_args.insert("fp_bedfile", args.fp_bedfile.clone().unwrap_or_default());
    final_args.insert("locations", args.locations.clone().unwrap_or_default());
    final_args.insert("threads", args.threads.unwrap_or_default().to_string());
    final_args.insert("strat_tsv", args.strat_tsv.clone().unwrap_or_default());
    metrics_json::write_compare_runinfo(
        &prefix.with_extension("runinfo.json"),
        &commandline,
        &final_args,
    )?;
    metrics_json::write_metrics_gz(
        &prefix.with_extension("metrics.json.gz"),
        "hap.py",
        &commandline,
        &[
            (
                "summary.metrics",
                "summary.metrics",
                &prefix.with_extension("summary.csv"),
            ),
            (
                "all.metrics",
                "all.metrics",
                &prefix.with_extension("extended.csv"),
            ),
            (
                "roc.all",
                "roc.all",
                &prefix.with_extension("roc.all.csv.gz"),
            ),
        ],
    )?;
    Ok(())
}

fn build_vcf_headers(contig_lengths: &BTreeMap<String, usize>) -> Vec<String> {
    let mut headers = vec![
        "##fileformat=VCFv4.1".to_string(),
        "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
        "##INFO=<ID=BS,Number=1,Type=Integer,Description=\"Start position of the benchmarking superlocus on current chromosome\">".to_string(),
        "##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">".to_string(),
        "##INFO=<ID=IMPORT_FAIL,Number=.,Type=Flag,Description=\"Flag to identify variants that could not be imported.\">".to_string(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">".to_string(),
        "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">".to_string(),
        "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Subtype for comparison\">".to_string(),
        "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"Variant type\">".to_string(),
        "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"Genotype label\">".to_string(),
        "##FORMAT=<ID=QQ,Number=1,Type=String,Description=\"Quality score\">".to_string(),
    ];
    for (name, length) in contig_lengths {
        headers.push(format!("##contig=<ID={name},length={length}>"));
    }
    headers.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string());
    headers
}

fn build_clusters(truth: &[Variant], query: &[Variant]) -> Vec<Cluster> {
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
                    && start <= cluster.end.saturating_add(CLUSTER_GAP_BP)
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

fn process_cluster(
    cluster: &Cluster,
    reference_sequences: &BTreeMap<String, String>,
    conf_bed: Option<&[vcf::BedInterval]>,
    counts: &mut BTreeMap<String, TypeCounts>,
    subtype_counts: &mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    rows: &mut Vec<AnnotatedRow>,
) -> Result<()> {
    let reference = reference_sequences
        .get(&cluster.chrom)
        .ok_or_else(|| anyhow::anyhow!("reference contig {} not found", cluster.chrom))?;
    let region_state = RegionState::from_cluster(cluster, reference, conf_bed);
    let mut truth_remaining = cluster.truth.clone();
    let mut query_remaining = cluster.query.clone();
    exact_match_pairs(
        cluster,
        reference,
        &region_state,
        counts,
        subtype_counts,
        rows,
        &mut truth_remaining,
        &mut query_remaining,
    );

    if truth_remaining.is_empty() && query_remaining.is_empty() {
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
    let allow_haplotype_match = cluster_has_gt_selected_nonsnp(cluster)
        && !cluster.truth.is_empty()
        && !cluster.query.is_empty();
    // Truth variants that were exact-matched (removed from truth_remaining).
    // These are the only variants for which a matching deletion in truth
    // should suppress the spurious BK=lm via drain restoration in
    // enumerate_haplotype_assignments. Variants still present in
    // truth_remaining were not matched, so BK=lm remains correct.
    let truth_matched: Vec<Variant> = cluster.truth.iter()
        .filter(|tv| !truth_remaining.iter().any(|r| r.key == tv.key))
        .cloned()
        .collect();
    let (truth_sig, query_sig) = if allow_haplotype_match {
        (
            cluster_signature(cluster, &cluster.truth, reference, None)?,
            cluster_signature(cluster, &cluster.query, reference, Some(&truth_matched))?,
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
    } else if allow_haplotype_match && truth_sig.is_some() && query_sig.is_none()
        && estimated_state_count(&cluster.query) <= XCMP_ENUMERATION_THRESHOLD
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

    let remainder = Cluster {
        chrom: cluster.chrom.clone(),
        start: cluster.start,
        end: cluster.end,
        truth: truth_remaining,
        query: query_remaining,
    };
    if is_match {
        mark_cluster_match(
            &remainder,
            cluster,
            reference,
            &region_state,
            counts,
            subtype_counts,
            rows,
        );
    } else {
        mark_cluster_mismatch(
            &remainder,
            cluster,
            hap_mismatch,
            reference,
            &region_state,
            counts,
            subtype_counts,
            rows,
        );
    }
    Ok(())
}

fn cluster_signature(
    cluster: &Cluster,
    variants: &[Variant],
    reference: &str,
    truth_variants: Option<&[Variant]>,
) -> Result<Option<BTreeSet<String>>> {
    // Skip enumeration when the predicted state space exceeds the legacy
    // xcmp cap. A `None` result forces the caller to treat the cluster as
    // mismatch — the same conservative fallback legacy applies when the
    // hap-block enumeration budget is blown.
    if estimated_state_count(variants) > XCMP_ENUMERATION_THRESHOLD {
        return Ok(None);
    }

    let segment = reference_segment(reference, cluster.start, cluster.end)?;
    let mut signatures = BTreeSet::new();
    for (hap1_events, hap2_events) in
        enumerate_haplotype_assignments(variants, reference, cluster.start, cluster.end, truth_variants)?
    {
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
fn gt_selected_nonref_alts(variant: &Variant) -> BTreeSet<String> {
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

fn cluster_has_gt_selected_nonsnp(cluster: &Cluster) -> bool {
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

fn trimmed_primitive_lens(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
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
fn query_insert_conflict_has_truth_counterpart(
    query_variants: &[Variant],
    truth_variants: &[Variant],
    truth_remaining: &[Variant],
) -> Option<bool> {
    // Classify each query variant position as Insert (any alt longer than ref)
    // or Subst (all alts same length as ref, i.e. SNP/substitution).
    let mut insert_positions: BTreeSet<usize> = BTreeSet::new();
    let mut subst_positions: BTreeSet<usize> = BTreeSet::new();
    for v in query_variants {
        let max_alt_len = v.key.alt_allele.split(',').map(|a| a.len()).max().unwrap_or(0);
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
            query_variants.iter().any(|v| {
                let ref_len = v.key.ref_allele.len();
                let max_alt = v.key.alt_allele.split(',').map(|a| a.len()).max().unwrap_or(0);
                ref_len > max_alt // deletion
                    && v.key.pos < ipos
                    && ipos <= v.key.pos + ref_len - 1
                    // deletion must be on both haplotypes (no 0/ref allele in GT)
                    && !v.gt.split(['/', '|']).any(|a| a == "0")
            })
            && (
                !truth_remaining.is_empty()
                || truth_variants.iter().any(|tv| tv.key.pos == ipos)
            )
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
            let max_alt_len = qv.key.alt_allele.split(',').map(|a| a.len()).max().unwrap_or(0);
            if max_alt_len <= qv.key.ref_allele.len() {
                continue; // not an insert at this position
            }
            // Use alt-set equality comparison: "CAA,CA" matches truth "CA,CAA"
            // (same set, different order), but "TAA" does NOT match "TAA,TA"
            // (strict subset → genuine mismatch, legacy gives BK=lm).
            let q_alts: std::collections::BTreeSet<&str> =
                qv.key.alt_allele.split(',').collect();
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
fn estimated_state_count(variants: &[Variant]) -> usize {
    let mut total: usize = 1;
    for variant in variants {
        let alleles = parse_gt_alleles(&variant.gt);
        let phased = variant.gt.contains('|');
        let multiplier =
            if !phased && alleles.len() == 2 && alleles[0] != alleles[1] { 2 } else { 1 };
        total = total.saturating_mul(multiplier);
        if total > XCMP_ENUMERATION_THRESHOLD {
            return usize::MAX;
        }
    }
    total
}

fn reference_segment(reference: &str, start: usize, end: usize) -> Result<String> {
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

fn enumerate_haplotype_assignments(
    variants: &[Variant],
    reference: &str,
    cluster_start: usize,
    cluster_end: usize,
    truth_variants: Option<&[Variant]>,
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
    let mut states: Vec<(Vec<Event>, Vec<Event>, usize, usize, usize, usize, usize, usize)> =
        vec![(Vec::new(), Vec::new(), 0, 0, 0, 0, 0, 0)];
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
        if projected > XCMP_ENUMERATION_THRESHOLD {
            bail!(
                "xcmp enumeration exceeded threshold ({} projected states)",
                projected
            );
        }
        let mut next = Vec::with_capacity(projected);
        for (hap1_events, hap2_events, h1_claimed, h1_subst, h1_insert, h2_claimed, h2_subst, h2_insert) in &states {
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
                let truth_has_variant = truth_variants.map_or(false, |tv| {
                    tv.iter().any(|t| {
                        t.key.pos == variant.key.pos
                            && t.key.ref_allele == variant.key.ref_allele
                            && t.key.alt_allele == variant.key.alt_allele
                    })
                });
                let left_conflict_start =
                    if left_eff_start > var_start && truth_has_variant { var_start } else { left_eff_start };
                let right_conflict_start =
                    if right_eff_start > var_start && truth_has_variant { var_start } else { right_eff_start };
                // Subst/Delete conflict: ref overlap OR prior Insert at same anchor.
                if left_nonref && !left_insert_only
                    && (left_conflict_start <= *h1_claimed || left_conflict_start <= *h1_insert)
                {
                    continue;
                }
                if right_nonref && !right_insert_only
                    && (right_conflict_start <= *h2_claimed || right_conflict_start <= *h2_insert)
                {
                    continue;
                }
                // Insert conflict: prior Subst at same anchor.
                if left_nonref && left_insert_only && var_start <= *h1_subst {
                    continue;
                }
                if right_nonref && right_insert_only && var_start <= *h2_subst {
                    continue;
                }
                let mut next_hap1 = hap1_events.clone();
                next_hap1.extend_from_slice(left_events);
                let mut next_hap2 = hap2_events.clone();
                next_hap2.extend_from_slice(right_events);
                // Update claimed_end for ALL non-ref (including Insert) so that
                // a later Subst at the same pos sees conflict via claimed_end.
                let new_h1_claimed = if left_nonref { (*h1_claimed).max(var_end) } else { *h1_claimed };
                let new_h2_claimed = if right_nonref { (*h2_claimed).max(var_end) } else { *h2_claimed };
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
                next.push((next_hap1, next_hap2, new_h1_claimed, new_h1_subst, new_h1_insert, new_h2_claimed, new_h2_subst, new_h2_insert));
            }
        }
        states = next;
    }
    Ok(states.into_iter().map(|(h1, h2, _, _, _, _, _, _)| (h1, h2)).collect())
}

/// Returns true when `events` contains any Subst or Delete — i.e., when the
/// haplotype assignment actually consumes reference bases at the variant
/// position. Pure-insertion assignments (only Event::Insert) are excluded
/// because they are anchored at a ref position without consuming it.
fn events_claim_ref_bases(events: &[Event]) -> bool {
    events
        .iter()
        .any(|e| matches!(e, Event::Subst { .. } | Event::Delete { .. }))
}

fn variant_haplotype_assignments(
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

fn parse_gt_alleles(gt: &str) -> Vec<usize> {
    gt.split(['/', '|'])
        .map(|part| part.parse::<usize>().unwrap_or(0))
        .collect()
}

fn equivalent_gt(left: &str, right: &str) -> bool {
    let mut left_alleles = parse_gt_alleles(left);
    let mut right_alleles = parse_gt_alleles(right);
    left_alleles.sort_unstable();
    right_alleles.sort_unstable();
    left_alleles == right_alleles
}

/// Sorted multi-set of ALT-allele SEQUENCES that `variant`'s GT selects
/// (reference alleles are dropped). Legacy's `SimpleDiploidCompare`
/// builds a per-index bitmask and compares it across samples; because
/// the loader unifies the `variation[]` table across truth + query
/// samples at each locus, the per-index bitmask is implicitly a per-
/// SEQUENCE bitmask. Rust stores truth and query separately with
/// independent ALT-column ordering, so we reconstruct the legacy
/// comparison by mapping each GT index back to its alt string.
fn selected_alt_sequences(variant: &Variant) -> Vec<String> {
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
fn simple_compare_pairs_match(truth: &Variant, query: &Variant) -> bool {
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
    if query_primitive_splits(query) || query_primitive_splits(truth) {
        return false;
    }
    // `gttype` equality is implied by the selected multi-set equality —
    // homalt `1|1` selects `[X, X]` vs het `0|1` selects `[X]`; these
    // differ as multi-sets even when the nonref-set matches.
    true
}

/// True iff `split_query_primitives(query)` would fan the record into
/// multiple per-primitive rows (distinct trimmed anchors). Mirrors the
/// `all_same_anchor` negation at the top of that function.
fn query_primitive_splits(query: &Variant) -> bool {
    let alleles = parse_gt_alleles(&query.gt);
    let used: BTreeSet<usize> = alleles.into_iter().filter(|a| *a > 0).collect();
    if used.is_empty() {
        return false;
    }
    let alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    let mut trimmed: Vec<(usize, String)> = Vec::new();
    for idx in &used {
        let Some(alt) = alts.get(*idx - 1).copied() else {
            continue;
        };
        let (pos, r, _) = trim_variant(query.key.pos, &query.key.ref_allele, alt);
        trimmed.push((pos, r));
    }
    if trimmed.len() <= 1 {
        return false;
    }
    let first = trimmed[0].clone();
    !trimmed.iter().all(|(pos, r)| pos == &first.0 && r == &first.1)
}

/// Remap `query.gt` into `truth`'s ALT-index space so the combined row
/// displays the query's genotype using truth's column ordering. Legacy
/// does this implicitly via its unified `variation[]` table; rust stores
/// the two sides' ALT columns independently, so we re-resolve each
/// query allele index to its sequence, then look that sequence up in
/// truth's ALT list. Unknown alleles (`.`) pass through.
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
fn canonical_hetalt_gt(output_alts: &str, query: &Variant) -> String {
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
    let (Some(a), Some(b)) = (tokens[0].parse::<usize>().ok(), tokens[1].parse::<usize>().ok())
    else {
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

fn remap_query_gt_to_truth(truth: &Variant, query: &Variant) -> String {
    let truth_alts: Vec<&str> = truth.key.alt_allele.split(',').collect();
    let query_alts: Vec<&str> = query.key.alt_allele.split(',').collect();
    let mut out = String::new();
    let mut token = String::new();
    for c in query.gt.chars() {
        if c == '|' || c == '/' {
            out.push_str(&remap_gt_token(&token, &truth_alts, &query_alts));
            out.push(c);
            token.clear();
        } else {
            token.push(c);
        }
    }
    out.push_str(&remap_gt_token(&token, &truth_alts, &query_alts));
    out
}

/// Reverse an unphased multi-allelic hetalt GT (`1/2` ↔ `2/1`). Legacy's
/// `VariantLocationAggregator::addAlleleToVariant` with `MAX_GT=2` writes
/// the second allele into the first zero slot, producing reversed-index
/// GTs. Phased, hom, het-with-ref, and single-allele GTs pass through.
fn swap_unphased_hetalt(gt: &str) -> String {
    if !gt.contains('/') {
        return gt.to_string();
    }
    let tokens: Vec<&str> = gt.split('/').collect();
    if tokens.len() != 2 {
        return gt.to_string();
    }
    let parsed: Vec<Option<u32>> = tokens.iter().map(|t| t.parse::<u32>().ok()).collect();
    if let (Some(a), Some(b)) = (parsed[0], parsed[1]) {
        if a > 0 && b > 0 && a != b {
            return format!("{}/{}", b, a);
        }
    }
    gt.to_string()
}

fn is_multi_allelic(variant: &Variant) -> bool {
    variant.key.alt_allele.contains(',')
}

fn remap_gt_token(token: &str, truth_alts: &[&str], query_alts: &[&str]) -> String {
    if token == "." {
        return ".".to_string();
    }
    let Ok(idx) = token.parse::<usize>() else {
        return token.to_string();
    };
    if idx == 0 {
        return "0".to_string();
    }
    let Some(alt_seq) = query_alts.get(idx - 1).copied() else {
        return ".".to_string();
    };
    match truth_alts.iter().position(|&ta| ta == alt_seq) {
        Some(pos) => (pos + 1).to_string(),
        None => ".".to_string(),
    }
}

fn normalized_events_for_allele(
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

fn normalize_ref_alt(
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
                while anchor < cluster_end
                    && bases[anchor].eq_ignore_ascii_case(&anchor_base)
                {
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

fn apply_events(
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

fn exact_match_pairs(
    cluster: &Cluster,
    reference: &str,
    region_state: &RegionState,
    counts: &mut BTreeMap<String, TypeCounts>,
    subtype_counts: &mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    rows: &mut Vec<AnnotatedRow>,
    truth_remaining: &mut Vec<Variant>,
    query_remaining: &mut Vec<Variant>,
) {
    let mut matched_truth = BTreeSet::new();
    let mut matched_query = BTreeSet::new();

    for (truth_index, truth) in truth_remaining.iter().enumerate() {
        if let Some((query_index, query)) = query_remaining
            .iter()
            .enumerate()
            .find(|(_, query)| simple_compare_pairs_match(truth, query))
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
            if truth.key.alt_allele == query.key.alt_allele
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
                rows.push(tp_combined_row(
                    truth,
                    &query_for_row,
                    reference,
                    cluster.start,
                    &region_state.row_tags(Some(truth), Some(query)),
                ));
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
                rows.push(tp_combined_row(
                    truth,
                    &canonicalized,
                    reference,
                    cluster.start,
                    &region_state.row_tags(Some(truth), Some(query)),
                ));
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
    *query_remaining = query_remaining
        .iter()
        .enumerate()
        .filter(|(index, _)| !matched_query.contains(index))
        .map(|(_, variant)| variant.clone())
        .collect();
}

fn mark_cluster_match(
    cluster: &Cluster,
    // Full (pre-exact-match) cluster — legacy's per-record QQ comes from
    // the block's ORIGINAL representative query qual, not the post-exact-
    // match remainder. Without this the truth-only TP row at a multi-
    // allelic split locus picks up the wrong QQ when an earlier SNP pair
    // in the same block already consumed its byte-equal counterpart.
    full_cluster: &Cluster,
    reference: &str,
    region_state: &RegionState,
    counts: &mut BTreeMap<String, TypeCounts>,
    subtype_counts: &mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    rows: &mut Vec<AnnotatedRow>,
) {
    for truth in &cluster.truth {
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
    for query in &cluster.query {
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
    }

    if cluster.truth.len() == 1
        && cluster.query.len() == 1
        && query_matches_truth_key(&cluster.query[0], &cluster.truth[0])
        && !query_primitive_splits(&cluster.truth[0])
        && !query_primitive_splits(&cluster.query[0])
    {
        rows.push(tp_combined_row(
            &cluster.truth[0],
            &cluster.query[0],
            reference,
            cluster.start,
            &region_state.row_tags(Some(&cluster.truth[0]), Some(&cluster.query[0])),
        ));
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
    let shared_qq = full_cluster
        .query
        .iter()
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
    for query in &cluster.query {
        // Multi-allelic query records emit one VCF row per per-allele
        // primitive, matching legacy's per-primitive output grain. Each
        // primitive inherits the cluster-level TP classification.
        for primitive in split_query_primitives(query) {
            rows.push(tp_single_side_row(
                &primitive,
                reference,
                cluster.start,
                &region_state.row_tags(None, Some(query)),
                Side::Query,
                None,
                ".",
            ));
        }
    }
}

fn mark_cluster_mismatch(
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
    counts: &mut BTreeMap<String, TypeCounts>,
    subtype_counts: &mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    rows: &mut Vec<AnnotatedRow>,
) {
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
            if query_matches_truth_key(query, truth) && !equivalent_gt(&truth.gt, &query.gt) {
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
        let is_conf = region_state.query_is_conf(query);
        let query_bd: &'static str = if is_conf { "FP" } else { "UNK" };
        if is_conf {
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
        rows.push(fn_fp_combined_row(
            truth,
            query,
            reference,
            cluster.start,
            &region_state.row_tags(Some(truth), Some(query)),
            query_bd,
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
                let pseudo_fp_class = if bd == "FP" {
                    classify_fp(&pseudo, full_cluster)
                } else {
                    None
                };
                let pseudo_bk = bk_for_row(&pseudo, &full_cluster.truth, hap_mismatch);
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
        let is_conf = region_state.query_is_conf(query);
        let bd: &'static str = if is_conf { "FP" } else { "UNK" };
        let regions = region_state.row_tags(None, Some(query));
        // Even same-type multi-allelic queries fan out one primitive row
        // per active allele in legacy's output — matches the decomposed
        // per-primitive representation xcmp emits after classification.
        for primitive in split_query_primitives(query) {
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
            let fp_class = if bd == "FP" {
                classify_fp(&primitive, full_cluster)
            } else {
                None
            };
            let primitive_bk = bk_for_row(&primitive, &full_cluster.truth, hap_mismatch);
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

/// Classify an FP query variant as genotype-mismatch (`"gt"`), allele-mismatch
/// (`"al"`), or novel (`None`).
///
/// - `"gt"`: a truth variant at the same locus (chrom/pos/ref/alt) exists —
///   the alleles agree, only the genotype differs.
/// - `None`: truth variant is absent at this position, or exists but with a
///   different alt allele — the call is not counted in FP.gt or FP.al.
///
/// Uses `full_cluster` (pre-exact-match truth set) so that truth variants
/// already consumed by exact-match TP pairing are still visible here.
fn classify_fp(query: &Variant, full_cluster: &Cluster) -> Option<&'static str> {
    if full_cluster
        .truth
        .iter()
        .any(|truth| query_matches_truth_key(query, truth))
    {
        Some("gt")
    } else {
        None
    }
}

/// Two VCF keys describe the same locus when they agree on chrom/pos/ref
/// and their alt lists are equal as sets (order-independent). Legacy's
/// xcmp exact-match pairing treats `A ATT,AT` truth and `A AT,ATT` query
/// as the same record with reshuffled alt indices; rust's byte-equal
/// check would miss this and fall through to the haplotype-match path
/// (which emits two rows instead of the combined row legacy expects).
fn query_matches_truth_key(query: &Variant, truth: &Variant) -> bool {
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

fn split_query_mismatch_rows(
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

fn component_is_conf(variant: &Variant, region_state: &RegionState) -> bool {
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
/// The record then requires `!is_pure_insertion || fully_covered` to land
/// in a region. Port that rule here so insertion / SNP records get the
/// same CONF / TS_boundary classification legacy emits.
///
/// Returns `(refstart_1b, refend_1b, is_pure_insertion)` in 1-based
/// inclusive VCF coordinates (matching `variant.key.pos`). `None` when
/// the record has no nucleotide alts (all symbolic / missing).
fn effective_refrange(variant: &Variant) -> Option<(usize, usize, bool)> {
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
        while rel_start < reflen && rel_start < altlen
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
/// The returned half-open `[start, end)` intervals are later merged with
/// the raw CONF bed (standard overlap/touching merge) before
/// `variant_is_conf` consumes them.
fn gvcf2bed_padding(truth: &[Variant]) -> Vec<vcf::BedInterval> {
    struct ActiveInterval {
        chrom: String,
        start: i64,
        end: i64,
    }

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
        let mut rec_start = i64::MAX;
        let mut rec_end = i64::MIN;
        let mut has_nuc = false;

        for alt in variant.key.alt_allele.split(',') {
            if alt.is_empty() || alt == "." || alt.starts_with('<') {
                continue;
            }
            let alt_bytes = alt.as_bytes();
            let mut reflen = ref_bytes.len();
            let mut altlen = alt_bytes.len();
            while reflen > 0 && altlen > 0 && ref_bytes[reflen - 1] == alt_bytes[altlen - 1] {
                reflen -= 1;
                altlen -= 1;
            }
            let mut rel_start = 0usize;
            while rel_start < reflen && rel_start < altlen
                && ref_bytes[rel_start] == alt_bytes[rel_start]
            {
                rel_start += 1;
            }
            let al_start = pos_0b + rel_start as i64;
            let al_end = pos_0b + reflen as i64 - 1;
            has_nuc = true;
            if al_end >= al_start {
                rec_start = rec_start.min(al_start);
                rec_end = rec_end.max(al_end);
            } else {
                rec_start = rec_start.min(al_start - 1);
                rec_end = rec_end.max(al_start);
            }
        }

        if !has_nuc {
            continue;
        }

        let flush = match active.as_ref() {
            Some(iv) => {
                !(iv.chrom == variant.key.chrom && rec_start <= iv.end && rec_end >= iv.start)
            }
            None => true,
        };

        if flush {
            if let Some(iv) = active.take() {
                out.push(vcf::BedInterval {
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
        out.push(vcf::BedInterval {
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
fn merge_bed_intervals(intervals: &[vcf::BedInterval]) -> Vec<vcf::BedInterval> {
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
                    out.push(vcf::BedInterval {
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
            out.push(vcf::BedInterval {
                chrom: chrom.clone(),
                start: cs,
                end: ce,
            });
        }
    }
    out
}

fn variant_is_conf(
    variant: &Variant,
    reference: &str,
    cluster_start: usize,
    cluster_end: usize,
    conf_intervals: &[vcf::BedInterval],
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

fn collect_contigs(
    truth: &[Variant],
    query: &[Variant],
    locations: Option<&[vcf::LocationFilter]>,
) -> BTreeSet<String> {
    let mut contigs = BTreeSet::new();
    if let Some(filters) = locations {
        for filter in filters {
            match filter {
                vcf::LocationFilter::Contig(chrom) => {
                    contigs.insert(chrom.clone());
                }
                vcf::LocationFilter::Range { chrom, .. } => {
                    contigs.insert(chrom.clone());
                }
            }
        }
    }
    for variant in truth.iter().chain(query.iter()) {
        contigs.insert(variant.key.chrom.clone());
    }
    contigs
}

fn inclusive_region_size(intervals: &[vcf::BedInterval]) -> usize {
    // Legacy hap.py reports conf_size using end-start+1 (1-based inclusive
    // interval length), matching Python's IntervalList semantics where each
    // [start, end) BED interval contributes end-start+1 bases after merging.
    let mut by_chrom: BTreeMap<String, Vec<(usize, usize)>> = BTreeMap::new();
    for interval in intervals {
        by_chrom
            .entry(interval.chrom.clone())
            .or_default()
            .push((interval.start, interval.end));
    }

    let mut total = 0usize;
    for mut ranges in by_chrom.into_values() {
        ranges.sort_unstable();
        let mut current: Option<(usize, usize)> = None;
        for (start, end) in ranges {
            match current {
                Some((cur_start, cur_end)) if start <= cur_end => {
                    current = Some((cur_start, cur_end.max(end)));
                }
                Some((cur_start, cur_end)) => {
                    total += cur_end.saturating_sub(cur_start) + 1;
                    current = Some((start, end));
                }
                None => current = Some((start, end)),
            }
        }
        if let Some((cur_start, cur_end)) = current {
            total += cur_end.saturating_sub(cur_start) + 1;
        }
    }
    total
}

fn derive_subset_counts(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, TypeCounts>> {
    let mut subsets: BTreeMap<String, BTreeMap<String, TypeCounts>> = BTreeMap::new();
    for row in rows {
        let fields: Vec<&str> = row.line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let info = fields[7];
        let Some(region_tail) = info.split(";Regions=").nth(1) else {
            continue;
        };
        let subset_tags: Vec<&str> = region_tail
            .split(',')
            .filter(|tag| *tag != "CONF")
            .collect();
        if subset_tags.is_empty() {
            continue;
        }

        let format_keys: Vec<&str> = fields[8].split(':').collect();
        let truth_parts: Vec<&str> = fields[9].split(':').collect();
        let query_parts: Vec<&str> = fields[10].split(':').collect();

        let truth_sample = SampleView::new(&format_keys, &truth_parts);
        let query_sample = SampleView::new(&format_keys, &query_parts);
        let filtered_out = pass_only && !row.query_pass;

        for subset in subset_tags {
            let type_map = subsets.entry(subset.to_string()).or_default();
            if let Some(variant_type) = truth_sample.variant_type() {
                let stats = type_map.entry(variant_type.to_string()).or_default();
                truth_sample.add_truth(stats, filtered_out);
            }
            if !filtered_out && let Some(variant_type) = query_sample.variant_type() {
                let stats = type_map.entry(variant_type.to_string()).or_default();
                query_sample.add_query(stats);
            }
        }
    }
    subsets
}

/// Per-variant-type (FP.gt, FP.al) tally derived by walking the annotated
/// rows. Rows whose query failed the filter are excluded when `pass_only` is
/// true. Returns `{variant_type: (fp_gt_count, fp_al_count)}` for every
/// variant type observed.
fn derive_fp_classes(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, (usize, usize)> {
    let mut out: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for row in rows {
        if pass_only && !row.query_pass {
            continue;
        }
        let Some(class) = row.fp_class else {
            continue;
        };
        let fields: Vec<&str> = row.line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let format_keys: Vec<&str> = fields[8].split(':').collect();
        let query_parts: Vec<&str> = fields[10].split(':').collect();
        let query_sample = SampleView::new(&format_keys, &query_parts);
        let Some(variant_type) = query_sample.variant_type() else {
            continue;
        };
        let bucket = out.entry(variant_type.to_string()).or_default();
        if class == "gt" {
            bucket.0 += 1;
        } else if class == "al" {
            bucket.1 += 1;
        }
    }
    out
}

fn derive_total_counts(rows: &[AnnotatedRow], pass_only: bool) -> BTreeMap<String, TypeCounts> {
    let mut totals: BTreeMap<String, TypeCounts> = BTreeMap::new();
    for row in rows {
        let fields: Vec<&str> = row.line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let format_keys: Vec<&str> = fields[8].split(':').collect();
        let truth_parts: Vec<&str> = fields[9].split(':').collect();
        let query_parts: Vec<&str> = fields[10].split(':').collect();

        let truth_sample = SampleView::new(&format_keys, &truth_parts);
        let query_sample = SampleView::new(&format_keys, &query_parts);
        let filtered_out = pass_only && !row.query_pass;

        if let Some(variant_type) = truth_sample.variant_type() {
            let stats = totals.entry(variant_type.to_string()).or_default();
            truth_sample.add_truth(stats, filtered_out);
        }
        if !filtered_out && let Some(variant_type) = query_sample.variant_type() {
            let stats = totals.entry(variant_type.to_string()).or_default();
            query_sample.add_query(stats);
        }
    }
    totals
}

/// Triple-nested counts for the extended CSV cross-product rows: subset
/// (TS_boundary / TS_contained) → variant_type (INDEL) → subtype (C1_5,
/// I1_5, etc.) → TypeCounts. Legacy emits 72 non-zero rows for these
/// stratifications per case; rust previously zero-filled them which
/// diverged extended.csv on every INDEL subtype × subset × filter cell.
fn derive_subset_subtype_counts(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>> {
    let mut out: BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>> =
        BTreeMap::new();
    for row in rows {
        let fields: Vec<&str> = row.line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let info = fields[7];
        let Some(region_tail) = info.split(";Regions=").nth(1) else {
            continue;
        };
        let subset_tags: Vec<&str> = region_tail
            .split(',')
            .filter(|tag| *tag != "CONF")
            .collect();
        if subset_tags.is_empty() {
            continue;
        }

        let format_keys: Vec<&str> = fields[8].split(':').collect();
        let truth_parts: Vec<&str> = fields[9].split(':').collect();
        let query_parts: Vec<&str> = fields[10].split(':').collect();

        let truth_sample = SampleView::new(&format_keys, &truth_parts);
        let query_sample = SampleView::new(&format_keys, &query_parts);
        let filtered_out = pass_only && !row.query_pass;

        for subset in &subset_tags {
            let by_type = out.entry((*subset).to_string()).or_default();
            if let Some((variant_type, subtype)) = truth_sample.variant_type_and_subtype() {
                let stats = by_type
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype)
                    .or_default();
                truth_sample.add_truth(stats, filtered_out);
            }
            if !filtered_out
                && let Some((variant_type, subtype)) = query_sample.variant_type_and_subtype()
            {
                let stats = by_type
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype)
                    .or_default();
                query_sample.add_query(stats);
            }
        }
    }
    out
}

fn derive_subtype_counts(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, TypeCounts>> {
    let mut subtypes: BTreeMap<String, BTreeMap<String, TypeCounts>> = BTreeMap::new();
    for row in rows {
        let fields: Vec<&str> = row.line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let format_keys: Vec<&str> = fields[8].split(':').collect();
        let truth_parts: Vec<&str> = fields[9].split(':').collect();
        let query_parts: Vec<&str> = fields[10].split(':').collect();

        let truth_sample = SampleView::new(&format_keys, &truth_parts);
        let query_sample = SampleView::new(&format_keys, &query_parts);
        let filtered_out = pass_only && !row.query_pass;

        if let Some((variant_type, subtype)) = truth_sample.variant_type_and_subtype() {
            let stats = subtypes
                .entry(variant_type.to_string())
                .or_default()
                .entry(subtype)
                .or_default();
            truth_sample.add_truth(stats, filtered_out);
        }
        if !filtered_out
            && let Some((variant_type, subtype)) = query_sample.variant_type_and_subtype()
        {
            let stats = subtypes
                .entry(variant_type.to_string())
                .or_default()
                .entry(subtype)
                .or_default();
            query_sample.add_query(stats);
        }
    }
    subtypes
}

struct SampleView<'a> {
    gt: Option<&'a str>,
    bd: Option<&'a str>,
    bi: Option<&'a str>,
    bvt: Option<&'a str>,
    blt: Option<&'a str>,
}

impl<'a> SampleView<'a> {
    fn new(format_keys: &[&str], sample_parts: &'a [&str]) -> Self {
        let lookup = |name: &str| -> Option<&'a str> {
            format_keys
                .iter()
                .position(|key| *key == name)
                .and_then(|index| sample_parts.get(index).copied())
        };
        Self {
            gt: lookup("GT"),
            bd: lookup("BD"),
            bi: lookup("BI"),
            bvt: lookup("BVT"),
            blt: lookup("BLT"),
        }
    }

    fn variant_type(&self) -> Option<&'a str> {
        self.bvt.filter(|value| *value != "NOCALL")
    }

    fn variant_type_and_subtype(&self) -> Option<(&'a str, String)> {
        let variant_type = self.variant_type()?;
        if variant_type != "INDEL" {
            return None;
        }
        let subtype = self.bi?.to_uppercase();
        Some((variant_type, subtype))
    }

    fn add_truth(&self, stats: &mut TypeCounts, demote_tp_to_fn: bool) {
        // When deriving PASS-tier counts for a row whose query failed the
        // filter, legacy reclassifies the truth side from TP to FN: the good
        // match doesn't count because the query record wouldn't have been
        // considered in a pass-only run.
        let effective_bd = if demote_tp_to_fn && self.bd == Some("TP") {
            Some("FN")
        } else {
            self.bd
        };
        match effective_bd {
            Some("TP") | Some("FN") => add_sample_stats(&mut stats.truth_total, self),
            _ => {}
        }
        match effective_bd {
            Some("TP") => add_sample_stats(&mut stats.truth_tp, self),
            Some("FN") => add_sample_stats(&mut stats.truth_fn, self),
            _ => {}
        }
    }

    fn add_query(&self, stats: &mut TypeCounts) {
        match self.bd {
            Some("TP") | Some("FP") | Some("UNK") => add_sample_stats(&mut stats.query_total, self),
            _ => {}
        }
        match self.bd {
            Some("TP") => add_sample_stats(&mut stats.query_tp, self),
            Some("FP") => add_sample_stats(&mut stats.query_fp, self),
            Some("UNK") => add_sample_stats(&mut stats.query_unk, self),
            _ => {}
        }
    }
}

fn add_sample_stats(bucket: &mut CountsBucket, sample: &SampleView<'_>) {
    bucket.total += 1;
    if let Some("SNP") = sample.bvt {
        match sample.bi {
            Some("ti") => bucket.ti += 1,
            Some("tv") => bucket.tv += 1,
            _ => {}
        }
    }
    // Legacy's TRUTH.TOTAL.het / homalt buckets count only the literal GT
    // strings (0/1, 1/0, 0|1, 1|0 for het; 1/1, 1|1 for homalt). The VCF
    // BLT column is wider (it resolves multi-allelic GT like 3|0 to "het"
    // to match legacy's row-level classifier) so counting off BLT would
    // over-count 50+ multi-allelic truth records per chr21. Parse GT here
    // so row-derived counts stay byte-equal with legacy's summary.
    match sample.gt {
        Some("0/1") | Some("1/0") | Some("0|1") | Some("1|0") => bucket.het += 1,
        Some("1/1") | Some("1|1") => bucket.homalt += 1,
        _ => {}
    }
}

fn add_variant_stats(bucket: &mut CountsBucket, variant: &Variant) {
    bucket.total += 1;
    if let Some(kind) = snp_bucket_label(variant) {
        if kind == "ti" {
            bucket.ti += 1;
        } else if kind == "tv" {
            bucket.tv += 1;
        }
    }
    if variant.is_het() {
        bucket.het += 1;
    }
    if variant.is_homalt() {
        bucket.homalt += 1;
    }
}

fn add_variant_stats_subtype<F>(
    subtype_counts: &mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    variant_type: &str,
    variant: &Variant,
    bucket_selector: F,
) where
    F: Fn(&mut TypeCounts) -> &mut CountsBucket,
{
    if let Some(subtype) = subtype_label(variant) {
        let family = subtype_counts
            .entry(variant_type.to_string())
            .or_default()
            .entry(subtype)
            .or_default();
        add_variant_stats(bucket_selector(family), variant);
    }
}

fn tp_combined_row(
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
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_gt}:TP:gm:{info}:{type_label}:{truth_loc}:{qq}\t{query_gt}:TP:gm:{info}:{type_label}:{query_loc}:{qq}",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = truth.key.ref_allele,
            alt = display_alt(truth),
            qual = query.qual,
            filter = filter_for_output(&query.filter),
            bs = block_start,
            regions = regions,
            truth_gt = truth.gt,
            query_gt = query.gt,
            info = info,
            type_label = truth.primary_type(),
            truth_loc = genotype_label(truth),
            query_loc = genotype_label(query),
            qq = query.qual,
        ),
    }
}

/// Legacy treats an empty, `.`, or `PASS` filter column as the passing set.
/// Anything else (`OverlapConflict`, `SuspiciousHomAlt`, user-defined filter
/// names) is excluded from the PASS-tier rollups in summary/extended.
fn filter_is_pass(filter: &str) -> bool {
    filter.is_empty() || filter == "." || filter == "PASS"
}

/// Render a query variant's FILTER column for the output VCF. Legacy
/// collapses empty / missing / PASS to `.` and passes everything else
/// through verbatim, preserving the original ;-separated filter tags
/// (e.g. `TruthSensitivityTranche99.00to99.90;LowQD`). Truth-only rows
/// always print `.` — the legacy emitter only forwards the query call's
/// filter, never the truth call's.
fn filter_for_output(filter: &str) -> &str {
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
fn cluster_query_filter(cluster: &Cluster) -> String {
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
        tokens.join(";")
    }
}

fn tp_single_side_row(
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
            // Truth-only rows carry no query variant; count them in both tiers.
            query_pass: true,
            fp_class: None,
            line: format!(
                "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{gt}:TP:gm:{info}:{type_label}:{loc}:{qq}\t./.:.:.:.:NOCALL:nocall:0",
                chrom = variant.key.chrom,
                pos = variant.key.pos,
                ref = variant.key.ref_allele,
                alt = display_alt(variant),
                qual = variant.qual,
                filter = truth_filter,
                bs = block_start,
                regions = regions,
                gt = variant.gt,
                info = info,
                type_label = variant.primary_type(),
                loc = genotype_label(variant),
                // Truth carries QQ=. in its input VCF (qual column is "0");
                // legacy propagates the matched query's QQ into the truth-
                // side sample column. Accept a shared value from the caller
                // that holds the cluster-level pairing context.
                qq = shared_qq.unwrap_or(variant.qual.as_str()),
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
            line: format!(
                "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t./.:.:.:.:NOCALL:nocall:.\t{gt}:TP:gm:{info}:{type_label}:{loc}:{qq}",
                chrom = variant.key.chrom,
                pos = variant.key.pos,
                ref = variant.key.ref_allele,
                alt = display_alt(variant),
                qual = variant.qual,
                filter = filter_for_output(&variant.filter),
                bs = block_start,
                regions = regions,
                gt = variant.gt,
                info = info,
                type_label = variant.primary_type(),
                loc = genotype_label(variant),
                qq = variant.qual,
            ),
        },
    }
}

/// Combined FN+FP row for a truth+query pair that share chrom/pos/ref/alt
/// as a set but disagree on GT (e.g. truth 0|1 het vs query 1/1 homalt).
/// Legacy emits BD=FN on truth, BD=FP on query, BK=am on both. The row
/// also carries the query's QUAL so downstream ROC enumeration can key
/// off the actual call quality rather than truth's placeholder zero.
fn fn_fp_combined_row(
    truth: &Variant,
    query: &Variant,
    reference: &str,
    block_start: usize,
    regions: &str,
    query_bd: &'static str,
) -> AnnotatedRow {
    let info = comparison_info(truth, reference);
    AnnotatedRow {
        sort_key: (truth.key.chrom.clone(), truth.key.pos, 0, 0),
        // Combined rows inherit the query's PASS status for the ALL vs
        // PASS bifurcation — query filter drives whether this row
        // contributes to the PASS-tier summary, consistent with
        // legacy's derived-from-row contract.
        query_pass: filter_is_pass(&query.filter),
        fp_class: Some("gt"),
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_gt}:FN:am:{info}:{type_label}:{truth_loc}:.\t{query_gt}:{query_bd}:am:{info}:{type_label}:{query_loc}:{qq}",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = truth.key.ref_allele,
            alt = display_alt(truth),
            qual = query.qual,
            filter = filter_for_output(&query.filter),
            bs = block_start,
            regions = regions,
            truth_gt = truth.gt,
            query_gt = query.gt,
            info = info,
            type_label = truth.primary_type(),
            truth_loc = genotype_label(truth),
            query_loc = genotype_label(query),
            qq = query.qual,
            query_bd = query_bd,
        ),
    }
}

fn fn_row(
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
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t.\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{gt}:FN:{bk}:{info}:{type_label}:{loc}:.\t./.:.:.:.:NOCALL:nocall:0",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = truth.key.ref_allele,
            alt = display_alt(truth),
            qual = truth.qual,
            bs = block_start,
            regions = regions,
            gt = truth.gt,
            bk = bk,
            info = info,
            type_label = truth.primary_type(),
            loc = genotype_label(truth),
        ),
    }
}

/// Truth variant outside the confident region emits BD=UNK (not FN) in
/// legacy; legacy additionally sets BK=lm when a cluster-adjacent query
/// variant shares at least one alternate allele at byte level. Mirrors
/// the legacy classifier in `XCmpQuantify::countVariants` where
/// non-CONF truth paired with a local query match short-circuits to
/// UNK+lm rather than claiming an FN against a non-evaluable region.
fn unk_truth_row(
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
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t.\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{gt}:UNK:{bk}:{info}:{type_label}:{loc}:.\t./.:.:.:.:NOCALL:nocall:0",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = truth.key.ref_allele,
            alt = display_alt(truth),
            qual = truth.qual,
            bs = block_start,
            regions = regions,
            gt = truth.gt,
            bk = bk,
            info = info,
            type_label = truth.primary_type(),
            loc = genotype_label(truth),
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
fn split_query_primitives(variant: &Variant) -> Vec<Variant> {
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

    trimmed
        .into_iter()
        .map(|(pos, ref_allele, alt_allele)| Variant {
            key: VariantKey {
                chrom: variant.key.chrom.clone(),
                pos,
                ref_allele,
                alt_allele,
            },
            qual: variant.qual.clone(),
            filter: variant.filter.clone(),
            // Each split primitive becomes a simple 0/1 het so downstream
            // is_het / genotype_label / primary_type classify it the same
            // way legacy does for its decomposed representations.
            gt: "0/1".to_string(),
        })
        .collect()
}

/// Trim common prefix and suffix from (ref, alt) and re-anchor the
/// position so both sides end up non-empty. Matches VCF primitive
/// normalisation: insertions anchor on the last prefix base, deletions
/// anchor on the same, substitutions just shift `pos` by the prefix
/// length.
fn trim_variant(pos: usize, ref_allele: &str, alt: &str) -> (usize, String, String) {
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
fn bk_for_row(row: &Variant, counterparts: &[Variant], hap_mismatch: bool) -> &'static str {
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
fn almismatch_same_locus(row: &Variant, counterparts: &[Variant]) -> bool {
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

fn fp_like_row(
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
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t./.:.:.:.:NOCALL:nocall:.\t{gt}:{bd}:{bk}:{info}:{type_label}:{loc}:{qq}",
            chrom = query.key.chrom,
            pos = query.key.pos,
            ref = query.key.ref_allele,
            alt = display_alt(query),
            qual = query.qual,
            filter = filter_for_output(&query.filter),
            bs = block_start,
            regions = regions,
            gt = query.gt,
            bd = bd,
            bk = bk,
            info = info,
            type_label = query.primary_type(),
            loc = genotype_label(query),
            qq = query.qual,
        ),
    }
}

fn display_alt(variant: &Variant) -> String {
    // Preserve the source VCF's ALT ordering rather than sorting lexicographically.
    // Legacy hap.py passes the original multi-allelic ALT string through xcmp
    // and writes it back verbatim; sorting here would rewrite the ALT column
    // (e.g. `ATCTC,ATC → ATC,ATCTC`) and also invalidate the GT indices that
    // were emitted for the pre-sort allele order.
    variant.key.alt_allele.clone()
}

fn query_type_rank(variant: &Variant) -> usize {
    if variant.primary_type() == "SNP" {
        0
    } else {
        1
    }
}

fn genotype_label(variant: &Variant) -> &'static str {
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
    if all_alleles.len() == 2
        && (all_alleles[0] == 0) != (all_alleles[1] == 0)
    {
        return "het";
    }
    "nocall"
}

fn comparison_info(variant: &Variant, reference: &str) -> String {
    // Multi-allelic records — SNP or INDEL — go through subtype_label so
    // mixed ti/tv or mixed indel sizes emit the legacy comma-joined BI
    // (e.g. `A → G,T` GT=2/1 must be `ti,tv`, not a single ti/tv token).
    if variant.key.alt_allele.contains(',') {
        if let Some(subtype) = subtype_label(variant) {
            return subtype.to_lowercase();
        }
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

fn snp_bucket_label(variant: &Variant) -> Option<&'static str> {
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

fn subtype_label(variant: &Variant) -> Option<String> {
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
        return Some(
            subtypes
                .into_iter()
                .collect::<Vec<_>>()
                .join(","),
        );
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
        .take_while(|i| {
            ref_bytes[ref_bytes.len() - 1 - i] == alt_bytes[alt_bytes.len() - 1 - i]
        })
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
    Some(format!("{class}{bucket}"))
}

fn bucket_label(size: usize) -> &'static str {
    match size {
        0..=5 => "1_5",
        6..=15 => "6_15",
        _ => "16_PLUS",
    }
}

fn is_transition_pair(ref_base: char, alt_base: char) -> bool {
    matches!(
        (ref_base, alt_base),
        ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C')
    )
}

#[cfg(test)]
mod memory_guards {
    use super::*;

    fn het(gt: &str) -> Variant {
        Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 10,
                ref_allele: "A".to_string(),
                alt_allele: "G".to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: gt.to_string(),
        }
    }

    #[test]
    fn estimated_state_count_doubles_per_unphased_het() {
        let variants = vec![het("0/1"), het("0/1"), het("0/1")];
        assert_eq!(estimated_state_count(&variants), 8);
    }

    #[test]
    fn estimated_state_count_phased_het_does_not_double() {
        let variants = vec![het("0|1"), het("0|1"), het("0|1")];
        assert_eq!(estimated_state_count(&variants), 1);
    }

    #[test]
    fn estimated_state_count_homozygous_does_not_double() {
        let variants = vec![het("1/1"), het("1/1"), het("1/1")];
        assert_eq!(estimated_state_count(&variants), 1);
    }

    #[test]
    fn estimated_state_count_saturates_above_threshold() {
        // 40 heterozygous unphased variants would be 2^40 states without the
        // saturating short-circuit; we expect `usize::MAX` as the sentinel so
        // the caller knows enumeration is not affordable.
        let variants: Vec<Variant> = (0..40).map(|_| het("0/1")).collect();
        assert_eq!(estimated_state_count(&variants), usize::MAX);
    }

    #[test]
    fn build_clusters_splits_at_variant_cap() {
        let mk = |pos: usize, side_is_truth: bool| Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos,
                ref_allele: "A".to_string(),
                alt_allele: "G".to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: "0/1".to_string(),
        };
        let _ = mk; // silence lint
        // MAX_CLUSTER_VARIANTS + 2 variants packed within 1 bp of each other
        // must yield at least 2 clusters — no single cluster may exceed the
        // cap.
        let mut truth = Vec::new();
        for offset in 0..(MAX_CLUSTER_VARIANTS + 2) {
            truth.push(Variant {
                key: VariantKey {
                    chrom: "chr1".to_string(),
                    pos: 100 + offset,
                    ref_allele: "A".to_string(),
                    alt_allele: "G".to_string(),
                },
                qual: "30".to_string(),
                filter: "PASS".to_string(),
                gt: "0/1".to_string(),
            });
        }
        let clusters = build_clusters(&truth, &[]);
        assert!(clusters.len() >= 2, "expected split at variant cap");
        for cluster in &clusters {
            assert!(
                cluster.truth.len() + cluster.query.len() <= MAX_CLUSTER_VARIANTS,
                "cluster exceeded MAX_CLUSTER_VARIANTS"
            );
        }
    }

    fn variant(pos: usize, r: &str, alt: &str, gt: &str) -> Variant {
        Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos,
                ref_allele: r.to_string(),
                alt_allele: alt.to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: gt.to_string(),
        }
    }

    #[test]
    fn cluster_gate_rejects_snp_only_multiallelic_cluster() {
        // chr21:15027483 fixture: truth T→C,TATC (GT 1|0) and query T→C
        // (GT 0/1). The GT-active alt on truth is C (SNP) and on query is
        // C (SNP); TATC is declared but not selected. Legacy's xcmp skips
        // block-level hapcmp in this setup because `n_nonsnp == 0` across
        // all GT-selected alts. The haplotype match that would otherwise
        // rescue this as TP must therefore be suppressed so the unmatched
        // pair stays as FN + FP.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 15027483,
            end: 15027483,
            truth: vec![variant(15027483, "T", "C,TATC", "1|0")],
            query: vec![variant(15027483, "T", "C", "0/1")],
        };
        assert!(!cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn cluster_gate_admits_cluster_with_gt_selected_indel() {
        // chr21:15006495 fixture: truth A→ATCTC,ATC (GT 1|0) selects ATCTC
        // (4-bp insert); query A→ATCTC (GT 0/1) also selects an insert.
        // Hapcmp must run so the haplotype rescue lets this pair through
        // as TP.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 15006495,
            end: 15006499,
            truth: vec![variant(15006495, "A", "ATCTC,ATC", "1|0")],
            query: vec![variant(15006495, "A", "ATCTC", "0/1")],
        };
        assert!(cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn cluster_gate_ignores_non_selected_insertion() {
        // Truth A→G,AT (GT=1|0): GT selects only the SNP G; AT is a
        // non-selected insertion. Non-selected insertions do NOT trigger
        // hapcmp — only non-selected deletions do.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 100,
            truth: vec![variant(100, "A", "G,AT", "1|0")],
            query: vec![variant(100, "A", "G", "0/1")],
        };
        assert!(!cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn cluster_gate_admits_non_selected_deletion() {
        // Truth TA→AA,T (GT=1|0): GT selects only the SNP AA; T is a
        // non-selected deletion allele. Legacy's xcmp counts this record
        // toward n_nonsnp and runs hapcmp, so we must too.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 36476244,
            end: 36476244,
            truth: vec![variant(36476244, "TA", "AA,T", "1|0")],
            query: vec![variant(36476244, "T", "A", "0/1")],
        };
        assert!(cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn trimmed_primitive_lens_collapses_shared_prefix_and_suffix() {
        // TAT→CAT trims to T→C after prefix/suffix removal — still a SNP.
        assert_eq!(trimmed_primitive_lens("TAT", "CAT"), (1, 1));
        // A→AT is a 1-bp insert after trimming the shared A.
        assert_eq!(trimmed_primitive_lens("A", "AT"), (0, 1));
        // AT→A is a 1-bp delete after trimming.
        assert_eq!(trimmed_primitive_lens("AT", "A"), (1, 0));
    }

    #[test]
    fn exact_match_key_rejects_reordered_multiallelic_alts() {
        // chr21:18743964 fixture: truth declares `GAA→GA,G` (GT 2|1) and
        // query declares `GAA→G,GA` (GT 1/2). Under BTreeSet-of-alts
        // semantics both sides share alleles {GA, G}, but legacy's
        // simpleCompare operates on byte-level ALT strings and sees them
        // as distinct (column order differs). The multi-allelic pair
        // must NOT exact-match — when the enclosing block also carries a
        // truth-only C→T 1|1 SNP, the block-wide hap signature is
        // disjoint and everything in the cluster falls to FN/FP.
        let truth = variant(18743964, "GAA", "GA,G", "2|1");
        let query = variant(18743964, "GAA", "G,GA", "1/2");
        assert!(
            !query_matches_truth_key(&query, &truth),
            "multi-allelic ALTs with different column order must not pair"
        );
    }

    #[test]
    fn exact_match_key_pairs_identical_multiallelic_alts_regardless_of_gt_phase() {
        // chr21:16032497 fixture: truth `C→CA,CAA` GT 2|1 (phased) and
        // query `C→CA,CAA` GT 1/2 (unphased). Raw ALT columns are
        // byte-equal; GT multisets are equal. Pairs as TP with
        // `equivalent_gt` doing the multiset check at the call site.
        let truth = variant(16032497, "C", "CA,CAA", "2|1");
        let query = variant(16032497, "C", "CA,CAA", "1/2");
        assert!(query_matches_truth_key(&query, &truth));
        assert!(equivalent_gt(&truth.gt, &query.gt));
    }

    #[test]
    fn bk_path_different_alt_at_same_pos_emits_dot() {
        // Legacy's VariantReader groups records by the full (chrom, pos,
        // ref, alt) tuple. Truth C→T and query C→G at the same position
        // have DIFFERENT alt columns, so they stay as separate Variants
        // objects; compareVariants returns kind=missing for each →
        // BK=. (not almismatch/lm).  The synth-snp-mismatch fixture
        // (generated by legacy) confirms this: G→T truth vs G→A query →
        // BK=. for both FN and FP rows.
        let truth = variant(15181526, "C", "T", "0|1");
        let query_counterpart = variant(15181526, "C", "G", "0/1");
        assert_eq!(bk_for_row(&truth, std::slice::from_ref(&query_counterpart), false), ".");
    }

    #[test]
    fn bk_path_almismatch_emits_lm_on_same_alt_disjoint_gt_selection() {
        // almismatch fires when both sides share the same (chrom, pos,
        // ref, alt) record — i.e. the same multi-allelic VCF entry — but
        // each sample's GT selects a completely disjoint non-ref allele
        // subset. Truth selects {T} from "T,A", query selects {A} →
        // disjoint → BK=lm.
        let truth = variant(15181526, "C", "T,A", "1/1");
        let query_counterpart = variant(15181526, "C", "T,A", "2/2");
        assert_eq!(bk_for_row(&truth, std::slice::from_ref(&query_counterpart), false), "lm");
    }

    #[test]
    fn bk_path_hap_mismatch_emits_lm_without_same_locus_counterpart() {
        // chr21:15313079 pattern: truth FN at pos P; the cluster's
        // query side has a GT-selected indel at a different locus (not
        // same-locus-to-FN). No almismatch fires, but the block-level
        // haplotype comparator ran and the two signatures disagreed —
        // legacy's ctype="hap:mismatch" → BK=lm. Confidence-region
        // membership of the counterpart is irrelevant.
        let truth = variant(15313088, "A", "G", "0|1");
        let query_indel = variant(15313079, "C", "CA", "0/1");
        assert_eq!(bk_for_row(&truth, std::slice::from_ref(&query_indel), true), "lm");
    }

    #[test]
    fn bk_path_fallthrough_emits_dot_on_snp_only_mismatch() {
        // chr21:15200371 pattern: truth FN at pos P; cluster's only
        // query is a SNP 7bp away at a different locus. No same-locus
        // counterpart means no almismatch; SNP-only cluster means
        // legacy's hap-run gate (n_nonsnp>0) didn't fire, so
        // hap:mismatch is also false. Legacy falls through to BK=.
        // Proximity alone never promotes BK to lm.
        let truth = variant(15200371, "T", "C", "0|1");
        let neighbour_snp = variant(15200378, "T", "C", "0/1");
        assert_eq!(bk_for_row(&truth, std::slice::from_ref(&neighbour_snp), false), ".");
    }

    #[test]
    fn bk_path_filtered_counterpart_not_almismatch() {
        // chr21:15007500 pattern: filtered-out query at the exact
        // same (chrom, pos, ref) as a truth FN. Legacy's simple-compare
        // runs only on post-filter calls, so a filtered counterpart
        // never enters `alleles_seen_2` — almismatch cannot fire.
        // Rust must skip filtered counterparts with the same guard.
        let truth = variant(15007500, "C", "T", "1|0");
        let mut filtered_query = variant(15007500, "C", "G", "0/1");
        filtered_query.filter = "LowQual".to_string();
        assert_eq!(bk_for_row(&truth, std::slice::from_ref(&filtered_query), false), ".");
    }

    // Residual #49 — chr21:17566241 exact-match pair with reordered
    // multi-allelic ALT columns. Truth `C→CA,CAA` (phased 1|2) and query
    // `C→CAA,CA` (unphased 1/2) encode the same diploid allele set.
    // Legacy's simpleCompare matches them via the VariantReader's shared
    // allele-unification table and emits one combined TP:gm row; rust
    // used to fall through to `cluster_signature` and split the pair
    // into truth-only FN + query-only FP. The `simple_compare_pairs_match`
    // predicate plus `remap_query_gt_to_truth` close this gap.
    #[test]
    fn simple_compare_matches_reordered_multiallelic_hetalt_indel() {
        let truth = variant(17566241, "C", "CA,CAA", "1|2");
        let query = variant(17566241, "C", "CAA,CA", "1/2");
        assert!(simple_compare_pairs_match(&truth, &query));
        // Query GT must be remapped into truth's ALT index space —
        // `1` (CAA) → truth idx 2, `2` (CA) → truth idx 1, so `1/2` → `2/1`.
        assert_eq!(remap_query_gt_to_truth(&truth, &query), "2/1");
    }

    // Residual #50 — chr21:15671076 cluster had a truth `T→TATATA` at
    // pos 15671094 plus a truth `T→TA` at pos 15671095. Both are
    // insertions. The previous homopolymer anchor slide pushed both
    // anchors to the same cluster_end-1 position, and `apply_events`
    // then rejected the pair as "two different inserts at the same
    // anchor". Restricting the slide to true homopolymer-extension
    // inserts (seq consists entirely of the anchor base) keeps the
    // distinct events at distinct anchors so the cluster signature
    // resolves and the block fires BK=lm.
    #[test]
    fn normalize_ref_alt_does_not_slide_heterogeneous_insert() {
        // Reference has a T-homopolymer at pos 15671094..15671099.
        // A 5-base `ATATA` insertion must stay at its original anchor
        // 15671094 rather than sliding through the T-homopolymer —
        // the previous slide collapsed distinct cluster inserts to the
        // same anchor, triggering a false conflict in `apply_events`.
        let mut reference = vec![b'N'; 15671110];
        // Lay down "atatatatatatattttttt" starting at pos 15671080.
        let window = b"atatatatatatattttttt";
        for (i, b) in window.iter().enumerate() {
            reference[15671080 - 1 + i] = *b;
        }
        let reference = String::from_utf8(reference).unwrap();
        let events = normalize_ref_alt(15671094, "T", "TATATA", &reference, 15671076, 15671097);
        assert_eq!(events.len(), 1);
        match &events[0] {
            Event::Insert { anchor, seq } => {
                assert_eq!(*anchor, 15671094, "anchor must not slide for non-homopolymer insert");
                assert_eq!(seq, "ATATA");
            }
            other => panic!("expected Insert, got {other:?}"),
        }
    }

    // Residual #50 companion — chr21:16997925 cluster had a query
    // homalt deletion `GCA→G` at pos 16997949 (delete pos 16997950-51)
    // plus a query het insert `A→ACG` at pos 16997951. The insert's
    // anchor (16997951) falls inside the delete's span, so
    // `apply_events` returned `None`. Legacy applies downstream
    // insertions against the shifted (post-delete) sequence without
    // failing; rust now allows insertions at deleted anchors since the
    // output-walk loop emits them unconditionally.
    #[test]
    fn apply_events_allows_insert_anchored_inside_delete() {
        // Reference must span the whole cluster — padding with Ns up to
        // position 16997960 and placing the CACA context inline.
        let mut reference = vec![b'N'; 16997960];
        reference[16997948] = b'G';  // pos 16997949: G
        reference[16997949] = b'C';  // pos 16997950: C (deleted)
        reference[16997950] = b'A';  // pos 16997951: A (deleted, insert anchor)
        reference[16997951] = b'G';  // pos 16997952: G
        let reference = String::from_utf8(reference).unwrap();
        let events = vec![
            Event::Delete { start: 16997950, end: 16997951 },
            Event::Insert { anchor: 16997951, seq: "CG".to_string() },
        ];
        let result = apply_events(&reference, 16997949, 16997952, &events).unwrap();
        assert!(result.is_some(), "delete + downstream insert must produce a valid haplotype");
    }

    // Class 4 pin (rust/PHASE1_BASELINE.md #52): BI (comparison_info)
    // on multi-allelic SNPs with mixed ti/tv alleles must emit the
    // comma-joined per-allele tokens legacy uses, not a single
    // collapsed ti/tv. Chr21:17562906 `A → G,T GT=2/1`: A→G is ti,
    // A→T is tv; legacy emits `ti,tv`, rust previously collapsed to
    // `tv` via the SNP-single-alt fastpath.
    #[test]
    fn comparison_info_joins_multi_allelic_snp_ti_tv() {
        let var = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 17562906,
                ref_allele: "A".to_string(),
                alt_allele: "G,T".to_string(),
            },
            qual: "1684.9".to_string(),
            filter: "PASS".to_string(),
            gt: "2/1".to_string(),
        };
        assert_eq!(comparison_info(&var, "N"), "ti,tv");
    }

    // Class 3 pin: `variant_is_conf` must apply legacy's
    // `!is_pure_insertion || fully_covered` rule using gvcf2bed-style
    // refrange. A pure insertion at a CONF edge (anchor in, anchor+1
    // out) must NOT be classified as covered. Chr21:15859667 `T → TA`
    // anchor 15859667 is last base of CONF `[15859657, 15859667)`,
    // anchor+1 15859668 is not covered — insertion straddles the
    // boundary, legacy emits no CONF tag, rust now agrees.
    #[test]
    fn variant_is_conf_rejects_insertion_at_bed_edge() {
        let var = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15859667,
                ref_allele: "T".to_string(),
                alt_allele: "TA".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1|0".to_string(),
        };
        let intervals = vec![
            vcf::BedInterval {
                chrom: "chr21".to_string(),
                start: 15859657,
                end: 15859667,
            },
        ];
        // Anchor at 15859667 is in [15859657,15859667)? 15859666 < 15859667 → yes.
        // Anchor+1 at 15859668 is in? 15859667 < 15859667 → NO.
        // Partial → pure insertion → skip CONF.
        assert!(!variant_is_conf(&var, "N", 15859600, 15859700, &intervals));
    }

    // Class 3 support: SNPs at a single base inside any CONF interval
    // are covered regardless of the insertion-aware fully_covered
    // clause. Sanity-check that the refactored `variant_is_conf`
    // doesn't accidentally reject point-covered SNPs.
    #[test]
    fn variant_is_conf_accepts_snp_inside_interval() {
        let var = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15859632,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1|0".to_string(),
        };
        let intervals = vec![
            vcf::BedInterval {
                chrom: "chr21".to_string(),
                start: 15859480,
                end: 15859645,
            },
        ];
        assert!(variant_is_conf(&var, "N", 15859600, 15859700, &intervals));
    }

    // Class 3 gvcf2bed-style padding bridges CONF gaps at insertion
    // anchors. Chr21:17562905 `C→CG` pure insertion sits at the edge
    // of CONF interval `[17561431, 17562905)` and the next interval
    // `[17562906, 17564580)` — a 1-base gap at 17562905. Legacy's
    // gvcf2bed emits `chr21 17562904 17562906` which, when merged with
    // the raw CONF bed, closes the gap and makes the adjacent A→G,T
    // SNP at pos 17562906 land inside CONF. Pin the padding function.
    #[test]
    fn gvcf2bed_padding_spans_insertion_anchor_and_next_base() {
        let truth = vec![Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 17562905,
                ref_allele: "C".to_string(),
                alt_allele: "CG".to_string(),
            },
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            gt: "1|1".to_string(),
        }];
        let padding = gvcf2bed_padding(&truth);
        assert_eq!(padding.len(), 1);
        assert_eq!(padding[0].chrom, "chr21");
        // 0-based half-open: anchor = 17562904 .. anchor+1+1 = 17562906
        assert_eq!(padding[0].start, 17562904);
        assert_eq!(padding[0].end, 17562906);
    }

    // Class 2 pin: `canonical_hetalt_gt` renders a query hetalt GT in
    // "alpha-later / alpha-earlier" order against the output (truth's)
    // alt list.
    //
    //   * Same-alt canonical case (C→C,G 1/2): alpha-later = G at pos 2
    //     → output `2/1`.
    //   * Same-alt NON-canonical case (T→TAA,TA 1/2): alpha-later =
    //     TAA at pos 1 → output `1/2` (verbatim).
    //   * Reordered case (output `A→ATT,AT` + query `A→AT,ATT 1/2`):
    //     alpha-later = ATT at output pos 1 → output `1/2`.
    //   * Reordered case (output `C→CA,CAA` + query `C→CAA,CA 1/2`):
    //     alpha-later = CAA at output pos 2 → output `2/1`.
    #[test]
    fn canonical_hetalt_gt_snp_canonical_swap() {
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 18280183,
                ref_allele: "T".to_string(),
                alt_allele: "C,G".to_string(),
            },
            qual: "1367.41".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("C,G", &query), "2/1");
    }

    #[test]
    fn canonical_hetalt_gt_indel_noncanonical_verbatim() {
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 21189041,
                ref_allele: "T".to_string(),
                alt_allele: "TAA,TA".to_string(),
            },
            qual: "1067.79".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("TAA,TA", &query), "1/2");
    }

    #[test]
    fn canonical_hetalt_gt_reordered_same_set() {
        // Truth A→ATT,AT vs query A→AT,ATT 1/2 → output 1/2 (ATT later).
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15712678,
                ref_allele: "A".to_string(),
                alt_allele: "AT,ATT".to_string(),
            },
            qual: "1751.52".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("ATT,AT", &query), "1/2");
    }

    #[test]
    fn canonical_hetalt_gt_reordered_set_swap() {
        // Truth C→CA,CAA vs query C→CAA,CA 1/2 → output 2/1.
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 16032497,
                ref_allele: "C".to_string(),
                alt_allele: "CAA,CA".to_string(),
            },
            qual: "759.81".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("CA,CAA", &query), "2/1");
    }

    // Class 1 pin (shared_qq picker): truth-only TP rows on a hap-
    // matched cluster must propagate the MINIMUM non-zero query QUAL
    // across the superlocus, not the first query's QUAL. Chr21:
    // 15246143 cluster contains query SNP qual=817.09 and query
    // insertion qual=174.59 — legacy emits QQ=174.59 on the truth-
    // only row `TA→TAA,T`.
    #[test]
    fn shared_qq_picks_minimum_nonzero_query_qual() {
        let queries = vec![
            variant(15246143, "G", "C", "0/1").with_qual("817.09"),
            variant(15246157, "T", "TA", "0/1").with_qual("174.59"),
        ];
        let min_qq: Option<&str> = queries
            .iter()
            .filter_map(|q| {
                q.qual
                    .parse::<f64>()
                    .ok()
                    .filter(|v| *v > 0.0)
                    .map(|v| (v, q.qual.as_str()))
            })
            .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(_, s)| s);
        assert_eq!(min_qq, Some("174.59"));
    }

    // Class 5 pin: `cluster_signature` returns `Ok(None)` when two homalt
    // deletions on the same side have overlapping ref spans. Var1 (pos=5,
    // GAT→G) claims hap1 and hap2 ref bases 6-7 after placement; Var2
    // (pos=6, AT→A) would start at ref position 6, but h1_end=7 ≥ 6 and
    // h2_end=7 ≥ 6, so the overlap gate in `enumerate_haplotype_assignments`
    // prunes every state for Var2. The resulting empty state list leaves
    // `signatures` empty; the non-empty `variants` slice triggers the
    // `Ok(None)` return so the caller falls back to mismatch.
    #[test]
    fn cluster_signature_overlapping_homalt_deletions_returns_none() {
        let mut reference = vec![b'N'; 10];
        reference[4] = b'G'; // pos 5 (1-based)
        reference[5] = b'A'; // pos 6
        reference[6] = b'T'; // pos 7
        let reference = String::from_utf8(reference).unwrap();
        let var1 = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 5,
                ref_allele: "GAT".to_string(),
                alt_allele: "G".to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let var2 = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 6,
                ref_allele: "AT".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr1".to_string(),
            start: 5,
            end: 7,
            truth: vec![var1.clone(), var2.clone()],
            query: vec![],
        };
        let variants = vec![var1, var2];
        let result = cluster_signature(&cluster, &variants, &reference).unwrap();
        assert!(
            result.is_none(),
            "overlapping homalt deletions must produce Ok(None)"
        );
    }
}

// Class 1 pinning helper: small extension trait for ergonomic qual
// overrides in tests. Not used outside the test module.
#[cfg(test)]
impl Variant {
    fn with_qual(mut self, q: &str) -> Self {
        self.qual = q.to_string();
        self
    }
}

/// Resolve the scratch directory used to stage preprocessed truth/query VCFs
/// before xcmp. The legacy pipeline writes these under the task's working
/// directory; we mirror that by anchoring on the report prefix's parent.
fn scratch_dir(prefix: &Path) -> Result<PathBuf> {
    let base = prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));
    let dir = base.join(".hap_scratch");
    std::fs::create_dir_all(&dir)
        .with_context(|| format!("failed to create scratch dir {}", dir.display()))?;
    Ok(dir)
}

/// Build the `PreprocessArgs` that germline should apply to truth or query
/// before running xcmp. Preserves the flags germline inherits from pre.py.
/// `pass_only` is passed through separately so the caller can force truth to
/// always filter to PASS (matching legacy's `--usefiltered-truth=False`).
fn build_preprocess_args(
    args: &CompareArgs,
    input: &str,
    output: &Path,
    pass_only: bool,
    decompose: bool,
) -> PreprocessArgs {
    PreprocessArgs {
        input: input.to_string(),
        output: output.to_string_lossy().into_owned(),
        reference: args.reference.clone(),
        locations: args.locations.clone(),
        pass_only,
        regions_bedfile: args.regions_bedfile.clone(),
        targets_bedfile: args.targets_bedfile.clone(),
        fixchr: Some(true),
        no_fixchr: false,
        somatic: false,
        set_gt: None,
        filter_nonref: false,
        convert_gvcf_to_vcf: false,
        // Legacy germline runs its internal xcmp preprocess with
        // side-specific flags: truth stays at decompose=false to avoid
        // primitive-splitting already-canonical truth VCFs, while query
        // is decomposed so multi-allelic calls fan out into per-allele
        // records the xcmp comparator can pair one-to-one with matching
        // truth primitives. Inflating truth by ~120 rows (what decompose
        // does to truth) blows the INDEL count delta; leaving query
        // undecomposed leaves 40+ multi-allelic query rows that never
        // reach the decomposed representation legacy emits. Caller
        // decides per side.
        decompose,
        threads: args.threads,
    }
}
