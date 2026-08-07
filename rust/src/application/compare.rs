use crate::adapters::vcf::{self, Variant, VariantKey};
use crate::adapters::{fasta, metrics_json, report};
use crate::application::{comparison_io, preprocess, roc_publication};
use crate::cli_compat::cli::{CompareArgs, CompareEngine, PreprocessArgs};
use crate::domain::{AnnotatedRow, CountsBucket, Interval, RawVcfRecord, TypeCounts};
use crate::engines::partial_credit;
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

mod genotype;
mod matching;
mod metrics;
mod output;
mod rows;

use genotype::{equivalent_gt, parse_gt_alleles};
use matching::*;
use metrics::*;
use output::*;
use rows::*;

static SCRATCH_RUN_ID: AtomicU64 = AtomicU64::new(0);

struct ScratchRun {
    path: PathBuf,
    keep: bool,
}

impl ScratchRun {
    fn create(parent: &Path, keep: bool) -> Result<Self> {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create scratch parent {}", parent.display()))?;
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();

        for _ in 0..128 {
            let id = SCRATCH_RUN_ID.fetch_add(1, Ordering::Relaxed);
            let path = parent.join(format!(
                "hap-compare-{}-{timestamp}-{id}",
                std::process::id()
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Ok(Self { path, keep }),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => {
                    return Err(error).with_context(|| {
                        format!("failed to create scratch run directory {}", path.display())
                    });
                }
            }
        }

        bail!(
            "failed to allocate a unique scratch run directory under {}",
            parent.display()
        )
    }

    fn path(&self) -> &Path {
        &self.path
    }

    fn cleanup(mut self) -> Result<()> {
        if self.keep {
            return Ok(());
        }
        fs::remove_dir_all(&self.path).with_context(|| {
            format!(
                "failed to remove scratch run directory {}",
                self.path.display()
            )
        })?;
        self.keep = true;
        Ok(())
    }
}

impl Drop for ScratchRun {
    fn drop(&mut self) {
        if !self.keep {
            let _ = fs::remove_dir_all(&self.path);
        }
    }
}

use crate::adapters::report::suffixed_report_path;

struct ComparisonOutputs<'a> {
    counts: &'a mut BTreeMap<String, TypeCounts>,
    subtype_counts: &'a mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    rows: &'a mut Vec<AnnotatedRow>,
}

#[derive(Clone, Copy, Debug)]
struct ComparisonConfig {
    no_hc: bool,
    max_enum: usize,
    hb_expand: usize,
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

#[cfg(test)]
const CLUSTER_GAP_BP: usize = 50;

/// Matches legacy hap.py's `--xcmp-enumeration-threshold` default. Beyond this
/// number of haplotype-assignment states we stop enumerating and treat the
/// cluster as a mismatch. Without the cap a dense cluster of ~30 het variants
/// explodes to 2^30 × two `Vec<Event>` allocations — hundreds of GB of RAM.
#[cfg(test)]
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

/// Window applied to `cluster_start` when bounding per-primitive left-shift
/// during the multi-allelic fan-out. Matches legacy hap.py's 1 kbp
/// `leftshift_limit` parameter so a slid primitive can never escape further
/// than ~1 kbp upstream of the cluster's anchor.
const SPLIT_LEFT_SHIFT_WINDOW: usize = 1024;

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
    conf_intervals: Vec<Interval>,
    covered_truth: BTreeSet<VariantKey>,
    covered_query: BTreeSet<VariantKey>,
}

impl RegionState {
    fn from_cluster(cluster: &Cluster, reference: &str, conf_bed: Option<&[Interval]>) -> Self {
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
            // Parent-keyed coverage feeds `query_is_conf(parent)` —
            // used by every emit path that hands the original (parent)
            // query record to row_tags / shared_qq filtering. Always
            // register the parent's `query.key` so those lookups stay
            // consistent regardless of the per-primitive fan-out shape.
            let parent_covered =
                variant_is_conf(query, reference, cluster.start, cluster.end, conf_bed);
            if parent_covered {
                state.covered_query.insert(query.key.clone());
            }

            // Per-primitive coverage: legacy's QuantifyRegions::annotate
            // computes Regions tags per-primitive after fan-out (see
            // hap.py:323), so a multi-allelic query whose deletion allele
            // is in CONF but whose insertion allele straddles a CONF edge
            // emits TWO rows with different Regions tags. Without
            // primitive-level tracking, all primitives inherit the
            // parent's "covered" verdict and the cluster-wide
            // any_nonconf bit misses the insertion's non-coverage —
            // collapsing TS_boundary into TS_contained on every record
            // and hiding UNK BD on the insertion itself.
            //
            // Pinning case (Class F — chr21:47906004): the parent
            // `AGAACTAAA→A,AAAA` has an `effective_refrange` reaching
            // into PG_Conf BED at 47906010 even though both fanned-out
            // primitives at 47906001 / 47906005 sit fully upstream of
            // CONF. When a parent record fans into per-primitive rows,
            // its any_conf / any_nonconf vote is replaced by the
            // primitives' votes — otherwise the cluster picks up a
            // spurious TS_boundary tag from the parent's coverage that
            // none of the emitted per-primitive rows actually touches.
            let primitives = split_query_primitives_with_neighbors(
                query,
                reference,
                cluster.start,
                &cluster.query,
                &cluster.truth,
            );
            let fans_out =
                primitives.len() > 1 || (primitives.len() == 1 && primitives[0].key != query.key);
            if !fans_out {
                if parent_covered {
                    state.any_conf = true;
                } else {
                    state.any_nonconf = true;
                }
            }
            for primitive in primitives {
                if primitive.key == query.key {
                    // Same key as parent — already accounted for above.
                    continue;
                }
                let primitive_covered =
                    variant_is_conf(&primitive, reference, cluster.start, cluster.end, conf_bed);
                if primitive_covered {
                    state.any_conf = true;
                    state.covered_query.insert(primitive.key);
                } else {
                    state.any_nonconf = true;
                }
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

fn validate_report_parent(prefix: &Path) -> Result<()> {
    let parent = prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    if !parent.exists() {
        bail!(
            "The output path does not exist. Please specify a valid output path and prefix using -o"
        );
    }
    Ok(())
}

fn comparison_requests_bcf(args: &CompareArgs) -> bool {
    args.bcf || (args.truth.ends_with(".bcf") && args.query.ends_with(".bcf"))
}

/// Preserve the pinned HAP-57 guard, including its allowance for a one-base
/// overlap (`previous_end - 1 == next_start`). Target BEDs intentionally skip
/// this check because legacy accepts the same interval layout with `-T`.
fn validate_regions_bed(path: &Path) -> Result<()> {
    let text = vcf::read_text(path)
        .with_context(|| format!("failed to inspect regions BED {}", path.display()))?;
    let mut previous_chrom: Option<&str> = None;
    let mut previous_end: Option<usize> = None;

    for (line_index, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 3 {
            continue;
        }
        let start = fields[1].parse::<usize>().with_context(|| {
            format!(
                "invalid BED start '{}' in {} at line {}",
                fields[1],
                path.display(),
                line_index + 1
            )
        })?;
        let end = fields[2].parse::<usize>().with_context(|| {
            format!(
                "invalid BED end '{}' in {} at line {}",
                fields[2],
                path.display(),
                line_index + 1
            )
        })?;

        if previous_chrom != Some(fields[0]) {
            previous_end = None;
        }
        if previous_end.is_some_and(|last| last.saturating_sub(1) > start) {
            bail!(
                "The regions bed file (specified using -R) has overlaps, this will not work with xcmp. You can either use -T, or run the file through bedtools merge"
            );
        }
        previous_chrom = Some(fields[0]);
        previous_end = Some(end);
    }
    Ok(())
}

pub(crate) fn run(mut args: CompareArgs) -> Result<()> {
    if args.version {
        println!("{}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    validate_report_parent(Path::new(&args.report_prefix))?;
    let explicit_bcf = args.bcf;
    args.bcf = comparison_requests_bcf(&args);
    if args.reference.is_empty() {
        args.reference = resolve_default_reference()?;
    }
    ensure_aggregate_roc_region(&mut args.roc_regions);
    normalize_engine_preprocessing(&mut args);
    initialize_compare_log(&args)?;
    log_compare_info(&args, "Starting germline comparison")?;
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
    if let Some(path) = args.regions_bedfile.as_deref() {
        validate_regions_bed(Path::new(path))?;
    }
    let regions = args
        .regions_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &contig_set))
        .transpose()?;
    let targets = args
        .targets_bedfile
        .as_ref()
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
    let scratch_parent = scratch_parent(&args, prefix);
    let scratch = ScratchRun::create(&scratch_parent, args.keep_scratch)?;
    let mut preprocessing_args = args.clone();
    // Legacy keeps vcfeval handoff files as VCF even when `--bcf` requests a
    // BCF report. The built-in xcmp and SCMP engines can consume BCF
    // intermediates directly.
    let bcf_intermediates = args.bcf && args.engine != CompareEngine::Vcfeval;
    preprocessing_args.bcf = bcf_intermediates;
    if preprocessing_args.gender == crate::cli_compat::cli::PreprocessGender::Auto {
        preprocessing_args.gender = preprocess::infer_gender(Path::new(&args.truth))?;
    }
    let truth_prep = scratch.path().join(if bcf_intermediates {
        "truth.prep.bcf"
    } else {
        "truth.prep.vcf.gz"
    });
    let query_prep = scratch.path().join(if bcf_intermediates {
        "query.prep.bcf"
    } else {
        "query.prep.vcf.gz"
    });
    log_compare_info(&args, "Preprocessing truth")?;
    preprocess::run(build_preprocess_args(
        &preprocessing_args,
        &args.truth,
        &truth_prep,
        !args.usefiltered_truth,
        args.preprocess_truth,
    ))?;
    if args.locations.is_none() {
        let (_, preprocessed_truth) = vcf::load_raw_vcf(&truth_prep)?;
        if !preprocessed_truth
            .iter()
            .any(|record| contig_set.contains(&record.chrom))
        {
            bail!("Truth and reference have no chromosomes in common!");
        }
    }
    log_compare_info(&args, "Preprocessing query")?;
    preprocess::run(build_preprocess_args(
        &preprocessing_args,
        &args.query,
        &query_prep,
        args.pass_only,
        true,
    ))?;

    if args.engine == CompareEngine::Vcfeval {
        log_compare_info(&args, "Running vcfeval comparison")?;
        return run_vcfeval(
            &args,
            explicit_bcf,
            &truth_prep,
            &query_prep,
            prefix,
            scratch,
        );
    }
    if matches!(
        args.engine,
        CompareEngine::ScmpSomatic | CompareEngine::ScmpDistance
    ) {
        log_compare_info(&args, "Running SCMP comparison")?;
        return run_scmp(
            &args,
            explicit_bcf,
            &truth_prep,
            &query_prep,
            prefix,
            scratch,
        );
    }

    log_compare_info(&args, "Running Rust comparison")?;

    let (truth_headers, truth_raw) = vcf::load_raw_vcf(&truth_prep)?;
    let (query_headers, query_raw) = vcf::load_raw_vcf(&query_prep)?;

    let truth_variant_input =
        materialize_variant_input(&truth_prep, scratch.path(), "truth.variants.vcf.gz")?;
    let query_variant_input =
        materialize_variant_input(&query_prep, scratch.path(), "query.variants.vcf.gz")?;
    let mut truth = vcf::load_variants(
        &truth_variant_input,
        &contig_set,
        false,
        regions.as_deref(),
        targets.as_deref(),
        locations.as_deref(),
    )?;
    // CONF insertion padding is derived from the preprocessed truth stream,
    // including filtered records retained by `--usefiltered-truth`. xcmp
    // excludes those records as calls below, but gvcf2bed sees them first;
    // collisions with PASS records can therefore change the padding lane.
    let truth_for_conf_padding = truth.clone();
    let filtered_truth_keys: BTreeSet<VariantKey> = truth
        .iter()
        .filter(|variant| !variant.is_pass())
        .map(|variant| variant.key.clone())
        .collect();
    // `--usefiltered-truth` preserves filtered records through pre.py, but
    // legacy xcmp still excludes them from haplotype enumeration and output.
    // This distinction matters: letting these records reach Rust comparison
    // creates spurious truth primitives and can turn otherwise query-only
    // calls into TPs. Keep preprocessing byte-compatible, then apply xcmp's
    // PASS-only truth-call contract at the comparison boundary.
    retain_xcmp_truth_calls(&mut truth);
    let query = vcf::load_variants(
        &query_variant_input,
        &contig_set,
        false,
        regions.as_deref(),
        targets.as_deref(),
        locations.as_deref(),
    )?;

    let contigs_in_play = collect_contigs(&truth, &query, locations.as_deref());
    let subset_size = report_subset_size(
        &contig_non_n_lengths,
        &contigs_in_play,
        explicit_bcf,
        args.bcf && !explicit_bcf,
    );
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

    let cluster_gap = match args.engine {
        CompareEngine::ScmpSomatic => 0,
        CompareEngine::ScmpDistance => args.engine_scmp_distance,
        CompareEngine::Xcmp | CompareEngine::Vcfeval => args.window,
    };
    let clusters = build_clusters_with_gap(&truth, &query, cluster_gap);
    // Fold gvcf2bed-style insertion padding (derived from truth) into
    // the raw CONF bed before classification. Legacy hap.py does the
    // same in Python (hap.py:323): `args.strat_regions.append(
    // "CONF_VARS:" + gvcf2bed(vcf1))`, and QuantifyRegions's CONF-label
    // squash folds the CONF_VARS lane back under the "CONF" lane.
    // Merging here bridges the 1-base gaps legacy's CONF beds carry at
    // every interval boundary (e.g. chr21:17562905) whenever an
    // adjacent truth insertion would pad them.
    // gvcf2bed-style padding: produced once, used twice. The merged
    // form (`adjusted_conf_bed`) is the union of raw CONF ∪ padding and
    // drives per-variant `is_conf` classification. The un-merged
    // **per-file BED-length sum** is what legacy reports as
    // `Subset.IS_CONF.Size` (`region_sizes[CONF] += stop-start+1` for
    // every interval loaded, with no cross-file dedup — see
    // QuantifyRegions::load).
    let adjust_conf = args.adjust_conf_regions && !args.no_adjust_conf_regions;
    let gvcf_padding: Option<Vec<Interval>> = conf_bed.as_ref().map(|raw| {
        if adjust_conf {
            gvcf2bed_padding(&truth_for_conf_padding, Some(raw))
        } else {
            Vec::new()
        }
    });
    let adjusted_conf_bed: Option<Vec<Interval>> = conf_bed.as_ref().map(|raw| {
        let mut combined = raw.clone();
        if let Some(pad) = gvcf_padding.as_ref() {
            combined.extend(pad.iter().cloned());
        }
        merge_bed_intervals(&combined)
    });
    // Sum half-open BED lengths per file separately, exactly as legacy's
    // QuantifyRegions::load does. The raw CONF bed is already disjoint
    // (BED files are conventionally non-overlapping), so its sum equals
    // its merged span. The padding output, however, often lies inside
    // the CONF bed — counting it via a cross-file merge would silently
    // drop those overlapping bases (legacy adds them anyway).
    let raw_conf_size = conf_bed
        .as_ref()
        .map(|raw| {
            raw.iter()
                .map(|iv| iv.end.saturating_sub(iv.start))
                .sum::<usize>()
        })
        .unwrap_or(0);
    let padding_size = gvcf_padding
        .as_ref()
        .map(|pad| {
            pad.iter()
                .map(|iv| iv.end.saturating_sub(iv.start))
                .sum::<usize>()
        })
        .unwrap_or(0);
    let conf_size = raw_conf_size + padding_size;
    let mut rows = Vec::new();
    for cluster in clusters {
        process_cluster(
            &cluster,
            &reference_sequences,
            adjusted_conf_bed.as_deref(),
            ComparisonConfig {
                no_hc: args.no_hc || args.engine != CompareEngine::Xcmp,
                max_enum: args.max_enum,
                hb_expand: args.hb_expand,
            },
            &mut counts,
            &mut subtype_counts,
            &mut rows,
        )?;
    }

    sort_comparison_rows(&mut rows, &filtered_truth_keys);
    decorate_output_rows(
        &mut rows,
        &truth_raw,
        &query_raw,
        args.preserve_info,
        args.output_vtc,
        &args.roc,
    )?;

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
    let all_subset_fp = derive_subset_fp_classes(&rows, false);
    let pass_subset_fp = derive_subset_fp_classes(&rows, true);
    let all_subtype_fp = derive_subtype_fp_classes(&rows, false);
    let pass_subtype_fp = derive_subtype_fp_classes(&rows, true);
    let all_subset_subtype_fp = derive_subset_subtype_fp_classes(&rows, false);
    let pass_subset_subtype_fp = derive_subset_subtype_fp_classes(&rows, true);
    report::write_summary(
        &suffixed_report_path(prefix, "summary.csv"),
        &all_counts,
        &pass_counts,
        &all_fp,
        &pass_fp,
    )?;
    let write_counts = args.write_counts && !args.no_write_counts;
    if write_counts {
        report::write_extended(
            &suffixed_report_path(prefix, "extended.csv"),
            &all_counts,
            &pass_counts,
            &all_subtype,
            &pass_subtype,
            subset_size,
            conf_size,
            conf_bed.is_some(),
            &all_subset,
            &pass_subset,
            &all_subset_subtype,
            &pass_subset_subtype,
            &all_fp,
            &pass_fp,
            &all_subset_fp,
            &pass_subset_fp,
            &all_subtype_fp,
            &pass_subtype_fp,
            &all_subset_subtype_fp,
            &pass_subset_subtype_fp,
        )?;
    }
    let vcf_headers = build_vcf_headers(
        &truth_headers,
        &query_headers,
        args.pass_only,
        args.output_vtc,
        args.preserve_info,
        &args.roc,
    );
    let requantify = args.strat_tsv.is_some()
        || !args.strat_regions.is_empty()
        || args.strat_fixchr
        || args.roc != "QUAL"
        || args.roc_filter.is_some()
        || args.roc_regions.iter().any(|region| region != "*")
        || (args.roc_delta - 0.5).abs() > f64::EPSILON
        || args.ci_alpha != 0.0;
    let comparison_vcf = if requantify {
        scratch.path().join("comparison.vcf.gz")
    } else {
        suffixed_report_path(prefix, "vcf.gz")
    };
    let requantify_rows = requantify.then(|| sanitize_requantify_handoff_rows(&rows));
    report::write_vcf(
        &comparison_vcf,
        &vcf_headers,
        requantify_rows.as_deref().unwrap_or(&rows),
    )?;
    let roc_indices = if requantify {
        crate::application::quantify::run_from_compare(
            crate::cli_compat::cli::QuantifyArgs {
                input_vcf: comparison_vcf.display().to_string(),
                report_prefix: args.report_prefix.clone(),
                reference: args.reference.clone(),
                // Rust comparison rows already carry finalized GA4GH BD/BK/BVT
                // sample fields. Re-quantify those decisions while adding user
                // stratifications/ROC controls; XCMP mode would instead expect
                // the legacy pre-quantify INFO/type annotations.
                annotation_type: Some("ga4gh".to_string()),
                fp_bedfile: args.fp_bedfile.clone(),
                strat_tsv: args.strat_tsv.clone(),
                strat_regions: args.strat_regions.clone(),
                strat_fixchr: args.strat_fixchr,
                write_vcf: true,
                write_counts,
                output_vtc: false,
                preserve_info: false,
                // `--adjust-conf-regions` is a no-op when no confidence BED was
                // supplied. Passing the truth VCF through to qfy in that case
                // incorrectly turns the default happy setting into qfy's
                // standalone argument error.
                adjust_conf_regions: (adjust_conf && conf_bed.is_some())
                    .then(|| truth_variant_input.display().to_string()),
                threads: None,
                bcf: false,
                logfile: None,
                verbose: false,
                quiet: false,
                force_interactive: false,
                roc: args.roc.clone(),
                do_roc: !args.no_roc,
                roc_regions: args.roc_regions.clone(),
                roc_filter: args.roc_filter.clone(),
                roc_delta: args.roc_delta,
                ci_alpha: args.ci_alpha,
                no_json: args.no_json,
            },
            crate::application::quantify::CompareQuantifyMode {
                preserve_missing_nocall_bd: args.usefiltered_truth,
                ..Default::default()
            },
        )?
    } else {
        let indices = roc_publication::write_roc_files(prefix, &rows, subset_size, conf_size)?;
        if args.no_roc {
            compact_no_roc_outputs(prefix)?;
        }
        indices
    };
    let commandline = std::env::args().collect::<Vec<_>>().join(" ");
    let run_args = metrics_json::CompareRunArgs {
        truth: &args.truth,
        query: &args.query,
        reference: &args.reference,
        reports_prefix: &args.report_prefix,
        annotation_type: args.annotation_type.as_deref(),
        pass_only: args.pass_only,
        preprocessing_truth: args.preprocess_truth,
        preprocessing_leftshift: args.engine == CompareEngine::ScmpSomatic || !args.no_leftshift,
        preprocessing_decompose: effective_decomposition(&args),
        regions_bedfile: args.regions_bedfile.as_deref(),
        targets_bedfile: args.targets_bedfile.as_deref(),
        fp_bedfile: args.fp_bedfile.as_deref(),
        locations: args.locations.as_deref(),
        threads: args.threads.unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(usize::from)
                .unwrap_or(1)
        }),
        strat_tsv: args.strat_tsv.as_deref(),
        scratch_prefix: args.scratch_prefix.as_deref(),
        keep_scratch: args.keep_scratch,
        bcf: explicit_bcf,
        ci_alpha: args.ci_alpha,
        convert_gvcf_query: args.convert_gvcf_query,
        convert_gvcf_to_vcf: args.convert_gvcf_to_vcf,
        convert_gvcf_truth: args.convert_gvcf_truth,
        do_roc: !args.no_roc,
        engine: args.engine.legacy_name(),
        engine_scmp_distance: args.engine_scmp_distance,
        engine_vcfeval: args.engine_vcfeval.as_deref().unwrap_or("rtg"),
        engine_vcfeval_template: args.engine_vcfeval_template.as_deref(),
        filter_nonref: args.filter_nonref,
        filters_only: args.filters_only.as_deref(),
        fixchr: if args.no_fixchr {
            Some(false)
        } else {
            args.fixchr
        },
        fp_adjust_conf: args.adjust_conf_regions && !args.no_adjust_conf_regions,
        gender: match args.gender {
            crate::cli_compat::cli::PreprocessGender::Male => "male",
            crate::cli_compat::cli::PreprocessGender::Female => "female",
            crate::cli_compat::cli::PreprocessGender::Auto => "auto",
            crate::cli_compat::cli::PreprocessGender::None => "none",
        },
        hb_expand: args.hb_expand,
        logfile: args.logfile.as_deref(),
        max_enum: args.max_enum,
        no_hc: args.no_hc,
        output_vtc: args.output_vtc,
        preprocess_window: args.preprocess_window,
        preprocessing_norm: args.bcftools_norm,
        preserve_info: args.preserve_info,
        quiet: args.quiet,
        roc: &args.roc,
        roc_delta: args.roc_delta,
        roc_filter: args.roc_filter.as_deref(),
        roc_regions: &args.roc_regions,
        somatic: args.somatic,
        somatic_mode: somatic_mode_name(args.set_gt),
        strat_fixchr: args.strat_fixchr,
        strat_regions: &args.strat_regions,
        usefiltered_truth: args.usefiltered_truth,
        verbose: args.verbose,
        window: args.window,
        write_counts,
        write_json: !args.no_json,
        write_vcf: true,
    };
    metrics_json::write_compare_runinfo(
        &suffixed_report_path(prefix, "runinfo.json"),
        &commandline,
        &run_args,
    )?;
    let mut metric_tables = vec![(
        "summary.metrics",
        "summary.metrics",
        suffixed_report_path(prefix, "summary.csv"),
    )];
    if write_counts {
        metric_tables.push((
            "all.metrics",
            "all.metrics",
            suffixed_report_path(prefix, "extended.csv"),
        ));
    }
    for id in &roc_indices.table_order {
        let path = suffixed_report_path(prefix, &format!("{id}.csv.gz"));
        if path.exists() {
            metric_tables.push((id.as_str(), id.as_str(), path));
        }
    }
    let metric_table_refs = metric_tables
        .iter()
        .map(|(id, label, path)| (*id, *label, path.as_path()))
        .collect::<Vec<_>>();
    if !args.no_json {
        metrics_json::write_metrics_gz_for_module_with_indices(
            &suffixed_report_path(prefix, "metrics.json.gz"),
            "hap.py.comparison",
            "hap.py",
            &commandline,
            &metric_table_refs,
            Some(&roc_indices.tables),
        )?;
    }
    publish_bcf_output(&args, prefix)?;
    log_compare_info(&args, "Germline comparison completed successfully")?;
    scratch.cleanup()?;
    Ok(())
}

/// `vcf::load_variants` is intentionally a text-VCF reader. Preserve BCF as
/// the preprocessing format, then materialize a private VCF view only for the
/// xcmp/SCMP variant-selection layer that still consumes textual records.
fn materialize_variant_input(source: &Path, scratch: &Path, name: &str) -> Result<PathBuf> {
    if source.extension().and_then(|value| value.to_str()) != Some("bcf") {
        return Ok(source.to_path_buf());
    }
    let (headers, records) = vcf::load_raw_vcf(source)?;
    let output = scratch.join(name);
    vcf::write_raw_vcf(&output, &headers, &records)?;
    Ok(output)
}

/// Quantification always produces its annotated report as indexed VCF. In
/// legacy `--bcf` mode that report is subsequently published as BCF+CSI and
/// the temporary VCF pair is not part of the final artifact family.
fn publish_bcf_output(args: &CompareArgs, prefix: &Path) -> Result<()> {
    if !args.bcf {
        return Ok(());
    }
    let vcf_path = suffixed_report_path(prefix, "vcf.gz");
    let vcf_index = suffixed_report_path(prefix, "vcf.gz.tbi");
    let bcf_path = suffixed_report_path(prefix, "bcf");
    let (headers, records) = vcf::load_raw_vcf(&vcf_path)
        .with_context(|| format!("failed to prepare BCF report from {}", vcf_path.display()))?;
    vcf::write_raw_vcf(&bcf_path, &headers, &records)?;
    fs::remove_file(&vcf_path)?;
    if vcf_index.exists() {
        fs::remove_file(&vcf_index)?;
    }
    Ok(())
}

fn retain_xcmp_truth_calls(variants: &mut Vec<Variant>) {
    variants.retain(Variant::is_pass);
}

fn sort_comparison_rows(rows: &mut [AnnotatedRow], filtered_truth_keys: &BTreeSet<VariantKey>) {
    rows.sort_by(|left, right| {
        left.sort_key
            .0
            .cmp(&right.sort_key.0)
            .then(left.sort_key.1.cmp(&right.sort_key.1))
            // A query record sharing the exact key of a filtered truth call
            // occupies that merged record's slot in legacy xcmp, even though
            // the truth sample itself is excluded from comparison. It sorts
            // before other records at the locus (PASS truth included).
            .then_with(|| {
                let left_filtered = row_matches_variant_key(left, filtered_truth_keys);
                let right_filtered = row_matches_variant_key(right, filtered_truth_keys);
                right_filtered.cmp(&left_filtered)
            })
            .then(left.sort_key.2.cmp(&right.sort_key.2))
            .then(left.sort_key.3.cmp(&right.sort_key.3))
            .then_with(|| left.record.cmp(&right.record))
    });
}

fn row_matches_variant_key(row: &AnnotatedRow, keys: &BTreeSet<VariantKey>) -> bool {
    let mut fields = row.record.split('\t');
    let (Some(chrom), Some(pos), Some(_id), Some(ref_allele), Some(alt_allele)) = (
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
    ) else {
        return false;
    };
    let Ok(pos) = pos.parse::<usize>() else {
        return false;
    };
    keys.contains(&VariantKey {
        chrom: chrom.to_string(),
        pos,
        ref_allele: ref_allele.to_string(),
        alt_allele: alt_allele.to_string(),
    })
}

fn run_vcfeval(
    args: &CompareArgs,
    explicit_bcf: bool,
    truth_prep: &Path,
    query_prep: &Path,
    prefix: &Path,
    scratch: ScratchRun,
) -> Result<()> {
    let mut strat_regions = args.strat_regions.clone();
    if args.adjust_conf_regions
        && !args.no_adjust_conf_regions
        && let Some(conf_path) = args.fp_bedfile.as_deref()
    {
        let reference = fasta::read_sequences(Path::new(&args.reference))?;
        let contigs = reference.keys().cloned().collect::<BTreeSet<_>>();
        let raw_conf = vcf::load_bed(Path::new(conf_path), &contigs)?;
        let (_, truth_records) = vcf::load_raw_vcf(truth_prep)?;
        let truth = truth_records
            .into_iter()
            .map(|record| -> Result<Variant> {
                let effective_end = record.effective_end_pos(truth_prep)?;
                let span = effective_end.saturating_sub(record.pos).saturating_add(1);
                let mut ref_allele = record.ref_allele;
                if span > ref_allele.len()
                    && let Some(sequence) = reference.get(&record.chrom)
                    && let Some(slice) = sequence
                        .as_bytes()
                        .get(record.pos.saturating_sub(1)..effective_end)
                {
                    ref_allele = String::from_utf8_lossy(slice).to_ascii_uppercase();
                }
                Ok(Variant {
                    key: VariantKey {
                        chrom: record.chrom,
                        pos: record.pos,
                        ref_allele,
                        alt_allele: record.alt_allele,
                    },
                    qual: record.qual,
                    filter: record.filter,
                    gt: String::new(),
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let padding = gvcf2bed_padding(&truth, Some(&raw_conf));
        let padding_path = scratch.path().join("truth.conf-vars.bed");
        let mut output = fs::File::create(&padding_path)
            .with_context(|| format!("failed to create {}", padding_path.display()))?;
        for interval in padding {
            writeln!(
                output,
                "{}\t{}\t{}",
                interval.chrom, interval.start, interval.end
            )
            .with_context(|| format!("failed to write {}", padding_path.display()))?;
        }
        strat_regions.push(format!("CONF_VARS:{}", padding_path.display()));
    }
    if args.engine_vcfeval.is_some() || args.engine_vcfeval_template.is_some() {
        eprintln!(
            "warning: --engine-vcfeval-path and --engine-vcfeval-template are deprecated and ignored; --engine vcfeval is native Rust"
        );
    }
    let vcfeval_vcf = scratch.path().join("vcfeval.comparison.vcf.gz");
    comparison_io::run_vcfeval(
        truth_prep,
        query_prep,
        Path::new(&args.reference),
        crate::engines::vcfeval::Options {
            roc_field: &args.roc,
            loose_match_distance: args.engine_scmp_distance,
        },
        &vcfeval_vcf,
    )?;
    let roc_indices = crate::application::quantify::run_from_compare(
        crate::cli_compat::cli::QuantifyArgs {
            input_vcf: vcfeval_vcf.display().to_string(),
            report_prefix: args.report_prefix.clone(),
            reference: args.reference.clone(),
            annotation_type: Some("ga4gh".to_string()),
            fp_bedfile: args.fp_bedfile.clone(),
            strat_tsv: args.strat_tsv.clone(),
            strat_regions,
            strat_fixchr: args.strat_fixchr,
            write_vcf: true,
            write_counts: args.write_counts && !args.no_write_counts,
            output_vtc: false,
            preserve_info: false,
            adjust_conf_regions: None,
            threads: None,
            bcf: false,
            logfile: None,
            verbose: false,
            quiet: false,
            force_interactive: false,
            roc: args.roc.clone(),
            do_roc: !args.no_roc,
            roc_regions: args.roc_regions.clone(),
            roc_filter: args.roc_filter.clone(),
            roc_delta: args.roc_delta,
            ci_alpha: args.ci_alpha,
            no_json: args.no_json,
        },
        crate::application::quantify::CompareQuantifyMode {
            inherit_same_position_tp_qq: true,
            roc_value_from_qq: true,
            ..Default::default()
        },
    )?;
    if args.preserve_info || args.output_vtc {
        decorate_existing_comparison_vcf(
            &suffixed_report_path(prefix, "vcf.gz"),
            truth_prep,
            query_prep,
            args.preserve_info,
            args.output_vtc,
        )?;
    }
    let commandline = std::env::args().collect::<Vec<_>>().join(" ");
    write_runinfo_for_args(args, explicit_bcf, prefix, &commandline)?;
    rewrite_compare_metrics(args, prefix, &commandline, &roc_indices)?;
    publish_bcf_output(args, prefix)?;
    log_compare_info(args, "Germline comparison completed successfully")?;
    scratch.cleanup()?;
    Ok(())
}

fn run_scmp(
    args: &CompareArgs,
    explicit_bcf: bool,
    truth_prep: &Path,
    query_prep: &Path,
    prefix: &Path,
    scratch: ScratchRun,
) -> Result<()> {
    let comparison_vcf = scratch.path().join("scmp.comparison.vcf.gz");
    let mut strat_regions = args.strat_regions.clone();
    if args.adjust_conf_regions
        && !args.no_adjust_conf_regions
        && let Some(conf_path) = args.fp_bedfile.as_deref()
    {
        let reference = fasta::read_sequences(Path::new(&args.reference))?;
        let contigs = reference.keys().cloned().collect::<BTreeSet<_>>();
        let raw_conf = vcf::load_bed(Path::new(conf_path), &contigs)?;
        let truth_variant_input =
            materialize_variant_input(truth_prep, scratch.path(), "truth.scmp-variants.vcf.gz")?;
        let truth = vcf::load_variants(&truth_variant_input, &contigs, false, None, None, None)?;
        let padding = gvcf2bed_padding(&truth, Some(&raw_conf));
        let padding_path = scratch.path().join("truth.conf-vars.bed");
        let mut output = fs::File::create(&padding_path)
            .with_context(|| format!("failed to create {}", padding_path.display()))?;
        for interval in padding {
            writeln!(
                output,
                "{}\t{}\t{}",
                interval.chrom, interval.start, interval.end
            )
            .with_context(|| format!("failed to write {}", padding_path.display()))?;
        }
        strat_regions.push(format!("CONF_VARS:{}", padding_path.display()));
    }
    let mode = match args.engine {
        CompareEngine::ScmpSomatic => crate::engines::scmp::ScmpMode::Alleles,
        CompareEngine::ScmpDistance => crate::engines::scmp::ScmpMode::Distance {
            max_distance: i64::try_from(args.engine_scmp_distance)
                .context("SCMP match distance exceeds the supported range")?,
        },
        CompareEngine::Xcmp | CompareEngine::Vcfeval => {
            bail!("internal error: run_scmp called for a non-SCMP engine")
        }
    };
    comparison_io::run_scmp(
        truth_prep,
        query_prep,
        Path::new(&args.reference),
        mode,
        &args.roc,
        &comparison_vcf,
    )?;
    let write_counts = args.write_counts && !args.no_write_counts;
    let roc_indices = crate::application::quantify::run_from_compare(
        crate::cli_compat::cli::QuantifyArgs {
            input_vcf: comparison_vcf.display().to_string(),
            report_prefix: args.report_prefix.clone(),
            reference: args.reference.clone(),
            annotation_type: Some("ga4gh".to_string()),
            fp_bedfile: args.fp_bedfile.clone(),
            strat_tsv: args.strat_tsv.clone(),
            strat_regions,
            strat_fixchr: args.strat_fixchr,
            write_vcf: true,
            write_counts,
            output_vtc: args.output_vtc,
            preserve_info: false,
            adjust_conf_regions: None,
            threads: None,
            bcf: false,
            logfile: None,
            verbose: false,
            quiet: false,
            force_interactive: false,
            roc: args.roc.clone(),
            do_roc: !args.no_roc,
            roc_regions: args.roc_regions.clone(),
            roc_filter: args.roc_filter.clone(),
            roc_delta: args.roc_delta,
            ci_alpha: args.ci_alpha,
            no_json: args.no_json,
        },
        crate::application::quantify::CompareQuantifyMode {
            preserve_missing_query_qq: args.engine == CompareEngine::ScmpSomatic,
            ..Default::default()
        },
    )?;
    let commandline = std::env::args().collect::<Vec<_>>().join(" ");
    write_runinfo_for_args(args, explicit_bcf, prefix, &commandline)?;
    rewrite_compare_metrics(args, prefix, &commandline, &roc_indices)?;
    publish_bcf_output(args, prefix)?;
    log_compare_info(args, "Germline comparison completed successfully")?;
    scratch.cleanup()?;
    Ok(())
}

fn ensure_aggregate_roc_region(regions: &mut Vec<String>) {
    if !regions.iter().any(|region| region == "*") {
        regions.insert(0, "*".to_string());
    }
}

fn normalize_engine_preprocessing(args: &mut CompareArgs) {
    match args.engine {
        CompareEngine::ScmpSomatic => {
            if !args.somatic && args.set_gt.is_none() {
                args.somatic = true;
            }
            // Legacy turns partial-credit normalization off for the somatic
            // scmp engine after selecting the synthetic half genotype.
            args.preprocess_truth = false;
            args.leftshift = false;
            args.no_leftshift = true;
            args.decompose = false;
            args.bcftools_norm = false;
        }
        CompareEngine::ScmpDistance => {
            if !args.somatic && args.set_gt.is_none() {
                args.set_gt = Some(crate::cli_compat::cli::SomaticGtMode::First);
            }
            args.decompose = false;
        }
        CompareEngine::Xcmp | CompareEngine::Vcfeval => {}
    }
}

fn somatic_mode_name(mode: Option<crate::cli_compat::cli::SomaticGtMode>) -> Option<&'static str> {
    mode.map(|mode| match mode {
        crate::cli_compat::cli::SomaticGtMode::Half => "half",
        crate::cli_compat::cli::SomaticGtMode::Hemi => "hemi",
        crate::cli_compat::cli::SomaticGtMode::Het => "het",
        crate::cli_compat::cli::SomaticGtMode::Hom => "hom",
        crate::cli_compat::cli::SomaticGtMode::First => "first",
    })
}

fn initialize_compare_log(args: &CompareArgs) -> Result<()> {
    let Some(path) = args.logfile.as_deref() else {
        return Ok(());
    };
    if let Some(parent) = Path::new(path)
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create logfile directory {}", parent.display()))?;
    }
    fs::write(path, "").with_context(|| format!("failed to initialize logfile {path}"))
}

fn log_compare_info(args: &CompareArgs, message: &str) -> Result<()> {
    if !args.verbose || args.quiet {
        return Ok(());
    }
    if let Some(path) = args.logfile.as_deref() {
        let mut file = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(path)
            .with_context(|| format!("failed to open logfile {path}"))?;
        writeln!(file, "INFO {message}")
            .with_context(|| format!("failed to write logfile {path}"))?;
    } else {
        eprintln!("INFO {message}");
    }
    Ok(())
}

fn write_runinfo_for_args(
    args: &CompareArgs,
    explicit_bcf: bool,
    prefix: &Path,
    commandline: &str,
) -> Result<()> {
    let write_counts = args.write_counts && !args.no_write_counts;
    let run_args = metrics_json::CompareRunArgs {
        truth: &args.truth,
        query: &args.query,
        reference: &args.reference,
        reports_prefix: &args.report_prefix,
        annotation_type: args.annotation_type.as_deref(),
        pass_only: args.pass_only,
        preprocessing_truth: args.preprocess_truth,
        preprocessing_leftshift: args.engine == CompareEngine::ScmpSomatic || !args.no_leftshift,
        preprocessing_decompose: effective_decomposition(args),
        regions_bedfile: args.regions_bedfile.as_deref(),
        targets_bedfile: args.targets_bedfile.as_deref(),
        fp_bedfile: args.fp_bedfile.as_deref(),
        locations: args.locations.as_deref(),
        threads: args.threads.unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(usize::from)
                .unwrap_or(1)
        }),
        strat_tsv: args.strat_tsv.as_deref(),
        scratch_prefix: args.scratch_prefix.as_deref(),
        keep_scratch: args.keep_scratch,
        bcf: explicit_bcf,
        ci_alpha: args.ci_alpha,
        convert_gvcf_query: args.convert_gvcf_query,
        convert_gvcf_to_vcf: args.convert_gvcf_to_vcf,
        convert_gvcf_truth: args.convert_gvcf_truth,
        do_roc: !args.no_roc,
        engine: args.engine.legacy_name(),
        engine_scmp_distance: args.engine_scmp_distance,
        engine_vcfeval: args.engine_vcfeval.as_deref().unwrap_or("rtg"),
        engine_vcfeval_template: args.engine_vcfeval_template.as_deref(),
        filter_nonref: args.filter_nonref,
        filters_only: args.filters_only.as_deref(),
        fixchr: if args.no_fixchr {
            Some(false)
        } else {
            args.fixchr
        },
        fp_adjust_conf: args.adjust_conf_regions && !args.no_adjust_conf_regions,
        gender: match args.gender {
            crate::cli_compat::cli::PreprocessGender::Male => "male",
            crate::cli_compat::cli::PreprocessGender::Female => "female",
            crate::cli_compat::cli::PreprocessGender::Auto => "auto",
            crate::cli_compat::cli::PreprocessGender::None => "none",
        },
        hb_expand: args.hb_expand,
        logfile: args.logfile.as_deref(),
        max_enum: args.max_enum,
        no_hc: args.no_hc,
        output_vtc: args.output_vtc,
        preprocess_window: args.preprocess_window,
        preprocessing_norm: args.bcftools_norm,
        preserve_info: args.preserve_info,
        quiet: args.quiet,
        roc: &args.roc,
        roc_delta: args.roc_delta,
        roc_filter: args.roc_filter.as_deref(),
        roc_regions: &args.roc_regions,
        somatic: args.somatic,
        somatic_mode: somatic_mode_name(args.set_gt),
        strat_fixchr: args.strat_fixchr,
        strat_regions: &args.strat_regions,
        usefiltered_truth: args.usefiltered_truth,
        verbose: args.verbose,
        window: args.window,
        write_counts,
        write_json: !args.no_json,
        write_vcf: true,
    };
    metrics_json::write_compare_runinfo(
        &suffixed_report_path(prefix, "runinfo.json"),
        commandline,
        &run_args,
    )
}

fn scratch_parent(args: &CompareArgs, prefix: &Path) -> PathBuf {
    if let Some(parent) = args.scratch_prefix.as_deref() {
        return PathBuf::from(parent);
    }
    let base = prefix
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .map(Path::to_path_buf)
        .unwrap_or_else(|| PathBuf::from("."));
    base.join(".hap_scratch")
}

fn compact_no_roc_outputs(prefix: &Path) -> Result<()> {
    let all_path = suffixed_report_path(prefix, "roc.all.csv.gz");
    let text = vcf::read_text(&all_path)?;
    let rows = text
        .lines()
        .enumerate()
        .filter(|(index, line)| *index == 0 || line.split(',').nth(6).is_some_and(|qq| qq == "*"))
        .map(|(_, line)| line.to_string())
        .collect::<Vec<_>>();
    let file = fs::File::create(&all_path)
        .with_context(|| format!("failed to create {}", all_path.display()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    writeln!(encoder, "{}", rows.join("\n"))?;
    encoder.finish()?;
    for suffix in [
        "roc.Locations.SNP.csv.gz",
        "roc.Locations.SNP.PASS.csv.gz",
        "roc.Locations.INDEL.csv.gz",
        "roc.Locations.INDEL.PASS.csv.gz",
    ] {
        let path = suffixed_report_path(prefix, suffix);
        if path.exists() {
            fs::remove_file(&path)
                .with_context(|| format!("failed to remove {}", path.display()))?;
        }
    }
    Ok(())
}

/// Build the `PreprocessArgs` that germline should apply to truth or query
/// before running xcmp. Preserves the flags germline inherits from pre.py.
/// `pass_only` is passed through separately so the caller can force truth to
/// always filter to PASS (matching legacy's `--usefiltered-truth=False`).
/// `preprocess_enabled` is false for truth by default and true for query,
/// matching hap.py's asymmetric preprocessing policy.
fn build_preprocess_args(
    args: &CompareArgs,
    input: &str,
    output: &Path,
    pass_only: bool,
    preprocess_enabled: bool,
) -> PreprocessArgs {
    let truth_side = input == args.truth;
    let convert_gvcf = args.convert_gvcf_to_vcf
        || if truth_side {
            args.convert_gvcf_truth
        } else {
            args.convert_gvcf_query
        };
    PreprocessArgs {
        input: input.to_string(),
        output: output.to_string_lossy().into_owned(),
        version: false,
        reference: Some(args.reference.clone()),
        locations: args.locations.clone(),
        pass_only,
        filters_only: (!truth_side).then(|| args.filters_only.clone()).flatten(),
        regions_bedfile: args.regions_bedfile.clone(),
        targets_bedfile: args.targets_bedfile.clone(),
        // Germline preprocessing inherits legacy's automatic prefix policy:
        // add `chr` only when the reference uses it and the VCF does not.
        fixchr: args.fixchr,
        no_fixchr: args.no_fixchr,
        somatic: args.somatic,
        set_gt: args.set_gt,
        filter_nonref: args.filter_nonref && (!truth_side || preprocess_enabled),
        convert_gvcf_to_vcf: convert_gvcf,
        bcf: args.bcf,
        bcftools_norm: args.bcftools_norm && (!truth_side || preprocess_enabled),
        leftshift: preprocess_enabled && (args.leftshift || !args.no_leftshift),
        no_leftshift: false,
        decompose: preprocess_enabled && effective_decomposition(args),
        no_decompose: false,
        gender: args.gender,
        window_size: args.preprocess_window as i64,
        threads: args.threads,
        logfile: None,
        verbose: args.verbose,
        quiet: args.quiet,
        force_interactive: args.force_interactive,
    }
}

/// Legacy hap.py disables decomposition for somatic/set-gt preprocessing
/// unless the user opts back in explicitly with `--decompose`.
fn effective_decomposition(args: &CompareArgs) -> bool {
    !args.no_decompose && (!(args.somatic || args.set_gt.is_some()) || args.decompose)
}

#[cfg(test)]
mod test_suite;
