use crate::adapters::vcf::{self, Variant, VariantKey};
use crate::adapters::{fasta, metrics_json, report};
use crate::application::{
    CompareArgs, CompareEngine, PreprocessArgs, PreprocessGender, QuantifyArgs, SomaticGtMode,
    ValidatedCompareArgs, comparison_io, preprocess, roc_publication,
};
#[cfg(test)]
use crate::domain::RawVcfRecord;
use crate::domain::{AnnotatedRow, Interval, TypeCounts};
use crate::output::{OutputTransaction, benchmark_artifacts, stratification_inputs};
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

mod genotype;
mod legacy_graph;
mod matching;
mod metrics;
mod output;
mod rows;
mod spool;

#[cfg(test)]
use genotype::equivalent_gt;
use matching::*;
use metrics::*;
use output::*;
use rows::*;
use spool::*;

static SCRATCH_RUN_ID: AtomicU64 = AtomicU64::new(0);

fn report_phase(name: &str, started: std::time::Instant) {
    if std::env::var_os("HAP_RS_PROFILE").is_some() {
        eprintln!(
            "HAP_RS_PHASE name={name} elapsed_seconds={:.3}",
            started.elapsed().as_secs_f64()
        );
    }
}

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

struct ProcessedCluster {
    cluster: Cluster,
    rows: Vec<AnnotatedRow>,
}

fn process_cluster_work(
    cluster: Cluster,
    reference_sequences: &BTreeMap<String, String>,
    conf_bed: Option<&[Interval]>,
    config: ComparisonConfig,
) -> Result<ProcessedCluster> {
    let mut counts: BTreeMap<String, TypeCounts> = BTreeMap::new();
    let mut subtype_counts: BTreeMap<String, BTreeMap<String, TypeCounts>> = BTreeMap::new();
    for variant in &cluster.truth {
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
    for variant in &cluster.query {
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
    let mut rows = Vec::new();
    process_cluster(
        &cluster,
        reference_sequences,
        conf_bed,
        config,
        &mut counts,
        &mut subtype_counts,
        &mut rows,
    )?;
    Ok(ProcessedCluster { cluster, rows })
}

fn process_clusters_parallel<I, F>(
    mut clusters: I,
    threads: usize,
    reference_sequences: &BTreeMap<String, String>,
    conf_bed: Option<&[Interval]>,
    config: ComparisonConfig,
    mut consume: F,
) -> Result<()>
where
    I: Iterator<Item = Result<Cluster>>,
    F: FnMut(ProcessedCluster) -> Result<()>,
{
    let worker_count = threads.max(1);
    // Keep enough ordered work in flight that one expensive cluster cannot
    // starve the remaining workers.  The old two-jobs-per-worker window
    // counted completed, out-of-order results against its limit; a slow
    // cluster at the front therefore stopped dispatch after only a handful
    // of later clusters completed.  A fixed per-worker window remains
    // bounded (and clusters themselves are capped at MAX_CLUSTER_VARIANTS)
    // while amortising those head-of-line stalls on real GIAB inputs.
    const CLUSTERS_INFLIGHT_PER_WORKER: usize = 64;
    let maximum_inflight = worker_count
        .saturating_mul(CLUSTERS_INFLIGHT_PER_WORKER)
        .max(1);
    if std::env::var_os("HAP_RS_PROFILE").is_some() {
        eprintln!(
            "HAP_RS_PROFILE comparison_workers={worker_count} comparison_max_inflight={maximum_inflight}"
        );
    }
    let queue = (
        std::sync::Mutex::new((VecDeque::<(usize, Cluster)>::new(), false)),
        std::sync::Condvar::new(),
    );
    let (result_sender, result_receiver) = std::sync::mpsc::sync_channel(maximum_inflight);
    let mut first_error = None;
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let queue = &queue;
            let result_sender = result_sender.clone();
            handles.push(scope.spawn(move || {
                loop {
                    let work = {
                        let mut state = queue.0.lock().unwrap_or_else(|error| error.into_inner());
                        while state.0.is_empty() && !state.1 {
                            state = queue
                                .1
                                .wait(state)
                                .unwrap_or_else(|error| error.into_inner());
                        }
                        state.0.pop_front()
                    };
                    let Some((index, cluster)) = work else {
                        break;
                    };
                    let result =
                        process_cluster_work(cluster, reference_sequences, conf_bed, config);
                    if result_sender.send((index, result)).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(result_sender);

        let mut source_finished = false;
        let mut next_index = 0usize;
        let mut next_to_consume = 0usize;
        let mut inflight = 0usize;
        let mut pending = BTreeMap::new();
        while !source_finished || inflight > 0 {
            while first_error.is_none() && !source_finished && inflight < maximum_inflight {
                match clusters.next() {
                    Some(Ok(cluster)) => {
                        let mut state = queue.0.lock().unwrap_or_else(|error| error.into_inner());
                        state.0.push_back((next_index, cluster));
                        next_index += 1;
                        inflight += 1;
                        drop(state);
                        queue.1.notify_one();
                    }
                    Some(Err(error)) => {
                        first_error = Some(error);
                        source_finished = true;
                    }
                    None => source_finished = true,
                }
            }
            if (source_finished || first_error.is_some())
                && !queue.0.lock().unwrap_or_else(|error| error.into_inner()).1
            {
                let mut state = queue.0.lock().unwrap_or_else(|error| error.into_inner());
                state.1 = true;
                drop(state);
                queue.1.notify_all();
                source_finished = true;
            }
            if inflight == 0 {
                break;
            }
            match result_receiver.recv() {
                Ok((index, result)) => {
                    inflight -= 1;
                    pending.insert(index, result);
                    while let Some(result) = pending.remove(&next_to_consume) {
                        if first_error.is_none() {
                            match result {
                                Ok(result) => {
                                    if let Err(error) = consume(result) {
                                        first_error = Some(error);
                                    }
                                }
                                Err(error) => first_error = Some(error),
                            }
                        }
                        next_to_consume += 1;
                    }
                }
                Err(error) => {
                    if first_error.is_none() {
                        first_error = Some(anyhow::anyhow!(
                            "comparison worker result channel closed: {error}"
                        ));
                    }
                    break;
                }
            }
        }
        {
            let mut state = queue.0.lock().unwrap_or_else(|error| error.into_inner());
            state.1 = true;
        }
        queue.1.notify_all();
        for handle in handles {
            if handle.join().is_err() && first_error.is_none() {
                first_error = Some(anyhow::anyhow!("comparison worker panicked"));
            }
        }
    });
    if let Some(error) = first_error {
        return Err(error);
    }
    Ok(())
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
        // `merge_bed_intervals` orders the effective CONF lane by chromosome
        // and coordinate. Retain only the tiny slice that can affect this
        // cluster. The previous implementation cloned the complete BED and
        // linearly searched every interval for every record in every cluster;
        // on a whole-genome query that made confidence annotation quadratic
        // in the number of variants and BED intervals.
        let conf_intervals = conf_bed
            .map(|intervals| cluster_conf_intervals(intervals, cluster))
            .unwrap_or_default();
        let mut state = Self {
            conf_enabled: conf_bed.is_some(),
            conf_intervals,
            ..Self::default()
        };
        if conf_bed.is_none() {
            return state;
        }
        let conf_bed = state.conf_intervals.as_slice();

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

fn cluster_conf_intervals(intervals: &[Interval], cluster: &Cluster) -> Vec<Interval> {
    let chrom_start =
        intervals.partition_point(|interval| interval.chrom.as_str() < cluster.chrom.as_str());
    let chrom_end =
        intervals.partition_point(|interval| interval.chrom.as_str() <= cluster.chrom.as_str());
    let chrom_intervals = &intervals[chrom_start..chrom_end];

    // BED is zero-based half-open while cluster coordinates are one-based
    // inclusive. Include an interval starting at cluster.end so insertion
    // coverage at the base immediately after the anchor remains visible.
    let range_start = chrom_intervals.partition_point(|interval| interval.end < cluster.start);
    let range_end = chrom_intervals.partition_point(|interval| interval.start <= cluster.end);
    chrom_intervals[range_start.min(range_end)..range_end].to_vec()
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

pub(crate) fn run(args: ValidatedCompareArgs) -> Result<()> {
    let mut args = args.into_inner();
    if args.version {
        println!("{}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    validate_report_parent(Path::new(&args.report_prefix))?;
    let explicit_bcf = args.bcf;
    args.bcf = comparison_requests_bcf(&args);
    if args.reference.is_empty() {
        bail!("no reference file found; pass --reference");
    }
    ensure_aggregate_roc_region(&mut args.roc.roc_regions);
    normalize_engine_preprocessing(&mut args);
    let destination_prefix = PathBuf::from(&args.report_prefix);
    let (inputs, mut labels) = compare_inputs(&args)?;
    labels.extend(
        args.roc
            .roc_regions
            .iter()
            .filter(|label| label.as_str() != "*")
            .cloned(),
    );
    let published_logfile = args.logfile.clone();
    let logfile = args.logfile.as_ref().map(PathBuf::from);
    let transaction =
        OutputTransaction::family(&inputs, &destination_prefix, benchmark_artifacts(labels))?
            .with_files(logfile.iter())?;
    if let Some(path) = logfile.as_deref() {
        args.logfile = Some(
            transaction
                .staged_file(path)?
                .to_string_lossy()
                .into_owned(),
        );
    }
    let staged_prefix = transaction.staged_prefix()?.to_path_buf();
    let comparison_started = std::time::Instant::now();
    run_inner(
        args,
        explicit_bcf,
        &staged_prefix,
        published_logfile.as_deref(),
    )
    .map_err(|error| {
        anyhow::anyhow!(
            "failed to produce report generation {}: {error:#}",
            destination_prefix.display()
        )
    })?;
    report_phase("comparison_pipeline", comparison_started);
    let publication_started = std::time::Instant::now();
    transaction.commit()?;
    report_phase("output_publication", publication_started);
    Ok(())
}

fn compare_inputs(args: &CompareArgs) -> Result<(Vec<PathBuf>, Vec<String>)> {
    let (indirect, labels) = stratification_inputs(args.strat_tsv.as_deref(), &args.strat_regions)?;
    let mut inputs = std::iter::once(args.truth.as_str())
        .chain(std::iter::once(args.query.as_str()))
        .chain(std::iter::once(args.reference.as_str()))
        .chain(args.regions_bedfile.as_deref())
        .chain(args.targets_bedfile.as_deref())
        .chain(args.fp_bedfile.as_deref())
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    inputs.extend(indirect);
    Ok((inputs, labels))
}

fn run_inner(
    args: CompareArgs,
    explicit_bcf: bool,
    output_prefix: &Path,
    published_logfile: Option<&str>,
) -> Result<()> {
    initialize_compare_log(&args)?;
    log_compare_info(&args, "Starting germline comparison")?;
    let reference_path = Path::new(&args.reference);
    let prefix = output_prefix;
    let reference_started = std::time::Instant::now();
    let reference_sequences = fasta::read_sequences(reference_path)?;
    report_phase("reference_loading", reference_started);
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
    let bcf_intermediates = args.bcf && args.engine.engine != CompareEngine::Vcfeval;
    preprocessing_args.bcf = bcf_intermediates;
    if preprocessing_args.preprocess.gender == PreprocessGender::Auto {
        preprocessing_args.preprocess.gender = preprocess::infer_gender(Path::new(&args.truth))?;
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
    let truth_preprocess_started = std::time::Instant::now();
    preprocess::run_with_reference(
        build_preprocess_args(
            &preprocessing_args,
            &args.truth,
            &truth_prep,
            !args.usefiltered_truth,
            args.preprocess_truth,
        )
        .validated()?,
        &reference_sequences,
    )?;
    report_phase("truth_preprocessing", truth_preprocess_started);
    if args.locations.is_none() {
        let mut preprocessed_truth = vcf::open_validated_vcf(&truth_prep)?;
        let mut common_contig = false;
        for record in &mut preprocessed_truth {
            if contig_set.contains(&record?.raw().chrom) {
                common_contig = true;
                break;
            }
        }
        if !common_contig {
            bail!("Truth and reference have no chromosomes in common!");
        }
    }
    log_compare_info(&args, "Preprocessing query")?;
    let query_preprocess_started = std::time::Instant::now();
    preprocess::run_with_reference(
        build_preprocess_args(
            &preprocessing_args,
            &args.query,
            &query_prep,
            args.pass_only,
            true,
        )
        .validated()?,
        &reference_sequences,
    )?;
    report_phase("query_preprocessing", query_preprocess_started);

    if args.engine.engine == CompareEngine::Vcfeval {
        log_compare_info(&args, "Running vcfeval comparison")?;
        return run_vcfeval(
            &args,
            explicit_bcf,
            &truth_prep,
            &query_prep,
            prefix,
            scratch,
            published_logfile,
        );
    }
    if matches!(
        args.engine.engine,
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
            published_logfile,
        );
    }

    log_compare_info(&args, "Running Rust comparison")?;

    let truth_headers = vcf::open_validated_vcf(&truth_prep)?.headers().to_vec();
    let query_headers = vcf::open_validated_vcf(&query_prep)?.headers().to_vec();
    let truth_metadata = spool_comparison_contigs(&truth_prep)?;
    // Without an explicit -l selector, legacy builds the comparison
    // chromosome list from truth calls. Query-only chromosomes are not sent
    // to xcmp (for example chrX in an autosome-only Platinum truth set).
    let derived_locations = args.locations.is_none().then(|| {
        truth_metadata
            .keys()
            .cloned()
            .map(vcf::LocationFilter::Contig)
            .collect::<Vec<_>>()
    });
    let comparison_locations = locations
        .as_deref()
        .or_else(|| derived_locations.as_deref());

    let cluster_gap = match args.engine.engine {
        CompareEngine::ScmpSomatic => 0,
        CompareEngine::ScmpDistance => args.engine.engine_scmp_distance,
        CompareEngine::Xcmp | CompareEngine::Vcfeval => args.engine.window,
    };
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
    let gvcf_padding: Option<Vec<Interval>> = if let Some(raw) = conf_bed.as_deref() {
        Some(if adjust_conf {
            gvcf2bed_padding_iter(
                vcf::open_variants(
                    &truth_prep,
                    &contig_set,
                    false,
                    regions.as_deref(),
                    targets.as_deref(),
                    comparison_locations,
                )?,
                Some(raw),
            )?
        } else {
            Vec::new()
        })
    } else {
        None
    };
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

    let truth = vcf::open_variants(
        &truth_prep,
        &contig_set,
        true,
        regions.as_deref(),
        targets.as_deref(),
        comparison_locations,
    )?;
    let query = vcf::open_variants(
        &query_prep,
        &contig_set,
        false,
        regions.as_deref(),
        targets.as_deref(),
        comparison_locations,
    )?;
    let contig_ranks = comparison_contig_ranks(&truth_headers, &query_headers);
    let clusters = StreamingClusters::new(truth, query, cluster_gap, contig_ranks);
    let mut contigs_in_play = comparison_locations
        .map(|locations| {
            locations
                .iter()
                .map(|location| match location {
                    vcf::LocationFilter::Contig(chrom)
                    | vcf::LocationFilter::Range { chrom, .. } => chrom.clone(),
                })
                .collect::<BTreeSet<_>>()
        })
        .unwrap_or_default();
    let mut row_spool = ComparisonRowSpool::new();
    let needs_decoration =
        args.preserve_info || args.output_vtc || !matches!(args.roc.roc.as_str(), "QUAL" | "QQ");
    let query_metadata = if needs_decoration {
        spool_comparison_contigs(&query_prep)?
    } else {
        BTreeMap::new()
    };
    let mut active_metadata: Option<ActiveComparisonMetadata> = None;
    let matching_started = std::time::Instant::now();
    let comparison_threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
    });
    process_clusters_parallel(
        clusters,
        comparison_threads,
        &reference_sequences,
        adjusted_conf_bed.as_deref(),
        ComparisonConfig {
            no_hc: args.no_hc || args.engine.engine != CompareEngine::Xcmp,
            max_enum: args.engine.max_enum,
            hb_expand: args.engine.hb_expand,
        },
        |processed| {
            let cluster = processed.cluster;
            if active_metadata
                .as_ref()
                .is_none_or(|metadata| metadata.chrom != cluster.chrom)
            {
                let truth_cursor = truth_metadata
                    .get(&cluster.chrom)
                    .map(ComparisonMetadataCursor::open)
                    .transpose()?;
                let query_cursor = query_metadata
                    .get(&cluster.chrom)
                    .map(ComparisonMetadataCursor::open)
                    .transpose()?;
                active_metadata = Some(ActiveComparisonMetadata {
                    chrom: cluster.chrom.clone(),
                    truth: truth_cursor,
                    query: query_cursor,
                });
            }
            contigs_in_play.insert(cluster.chrom.clone());
            let cluster_rows = processed.rows;
            let mut emitted_start = cluster.start;
            let mut emitted_end = cluster.end;
            for row in &cluster_rows {
                let record = row.record.raw();
                emitted_start = emitted_start.min(record.pos);
                emitted_end = emitted_end.max(record.end_pos());
            }
            let mut filtered_truth_keys = BTreeSet::new();
            let mut decorations = DecorationIndex::default();
            let metadata = active_metadata
                .as_mut()
                .expect("comparison metadata cursor was initialized");
            if let Some(truth_cursor) = &mut metadata.truth {
                truth_cursor.collect(
                    &cluster.chrom,
                    emitted_start,
                    emitted_end,
                    &mut filtered_truth_keys,
                    &mut decorations,
                    spool::CollectOptions {
                        preserve_info: args.preserve_info,
                        roc_field: &args.roc.roc,
                        collect_filtered: true,
                    },
                )?;
            }
            if needs_decoration && let Some(query_cursor) = &mut metadata.query {
                query_cursor.collect(
                    &cluster.chrom,
                    emitted_start,
                    emitted_end,
                    &mut filtered_truth_keys,
                    &mut decorations,
                    spool::CollectOptions {
                        preserve_info: args.preserve_info,
                        roc_field: &args.roc.roc,
                        collect_filtered: false,
                    },
                )?;
            }
            for mut row in cluster_rows {
                let sort_line = row.record.raw().to_line();
                let filtered_match = row_matches_variant_key(&row, &filtered_truth_keys);
                if needs_decoration {
                    decorate_output_rows_with_index(
                        std::slice::from_mut(&mut row),
                        &decorations,
                        args.preserve_info,
                        args.output_vtc,
                        &args.roc.roc,
                    )?;
                }
                row_spool.push(row, filtered_match, sort_line)?;
            }
            Ok(())
        },
    )?;
    report_phase("matching", matching_started);
    let sorting_started = std::time::Instant::now();
    let row_file = row_spool.finish()?;
    report_phase("external_sorting", sorting_started);
    let subset_size = report_subset_size(
        &contig_non_n_lengths,
        &contigs_in_play,
        explicit_bcf,
        args.bcf && !explicit_bcf,
    );
    let whole_reference_size = contig_non_n_lengths.values().sum();
    if subset_size == 0 {
        bail!("no reference contigs selected for analysis");
    }

    let mut tallies = FoldedComparisonReports::default();
    for row in row_file.rows()? {
        tallies.observe(&row?);
    }
    report::write_summary(
        &suffixed_report_path(prefix, "summary.csv"),
        &tallies.all_counts,
        &tallies.pass_counts,
        &tallies.all_fp,
        &tallies.pass_fp,
    )?;
    let write_counts = args.write_counts && !args.no_write_counts;
    if write_counts {
        report::write_extended(
            &suffixed_report_path(prefix, "extended.csv"),
            &report::ExtendedTables {
                all_counts: &tallies.all_counts,
                pass_counts: &tallies.pass_counts,
                all_subtype: &tallies.all_subtype,
                pass_subtype: &tallies.pass_subtype,
                all_subset: &tallies.all_subset,
                pass_subset: &tallies.pass_subset,
                all_subset_subtype: &tallies.all_subset_subtype,
                pass_subset_subtype: &tallies.pass_subset_subtype,
                all_fp: &tallies.all_fp,
                pass_fp: &tallies.pass_fp,
                all_subset_fp: &tallies.all_subset_fp,
                pass_subset_fp: &tallies.pass_subset_fp,
                all_subtype_fp: &tallies.all_subtype_fp,
                pass_subtype_fp: &tallies.pass_subtype_fp,
                all_subset_subtype_fp: &tallies.all_subset_subtype_fp,
                pass_subset_subtype_fp: &tallies.pass_subset_subtype_fp,
            },
            report::ExtendedSizes {
                subset_size,
                whole_reference_size,
                conf_size,
                has_conf_regions: conf_bed.is_some(),
            },
        )?;
    }
    let vcf_headers = build_vcf_headers(
        &truth_headers,
        &query_headers,
        args.pass_only,
        args.output_vtc,
        args.preserve_info,
        &args.roc.roc,
    );
    let requantify = args.strat_tsv.is_some()
        || !args.strat_regions.is_empty()
        || args.strat_fixchr
        || args.roc.roc != "QUAL"
        || args.roc.roc_filter.is_some()
        || args.roc.roc_regions.iter().any(|region| region != "*")
        || (args.roc.roc_delta - 0.5).abs() > f64::EPSILON
        || args.ci_alpha != 0.0;
    let comparison_vcf = if requantify {
        scratch.path().join("comparison.vcf.gz")
    } else {
        suffixed_report_path(prefix, "vcf.gz")
    };
    vcf::write_validated_vcf_iter(
        &comparison_vcf,
        &vcf_headers,
        row_file.rows()?.map(|row| {
            row.map(|row| {
                if requantify {
                    sanitize_requantify_handoff_row(row)
                } else {
                    row
                }
                .record
                .into_validated()
            })
        }),
    )?;
    let roc_started = std::time::Instant::now();
    let roc_indices = if requantify {
        crate::application::quantify::run_from_compare_path(
            QuantifyArgs {
                input_vcf: comparison_vcf.display().to_string(),
                // The outer comparison transaction owns publication to the
                // user-facing prefix. Quantification must therefore publish
                // into that transaction's staged family, not directly to the
                // destination family recorded in `args`.
                report_prefix: prefix.to_string_lossy().into_owned(),
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
                    .then(|| truth_prep.display().to_string()),
                threads: None,
                bcf: false,
                logfile: None,
                verbose: false,
                quiet: false,
                force_interactive: false,
                roc: args.roc.roc.clone(),
                do_roc: !args.roc.no_roc,
                roc_regions: args.roc.roc_regions.clone(),
                roc_filter: args.roc.roc_filter.clone(),
                roc_delta: args.roc.roc_delta,
                ci_alpha: args.ci_alpha,
                no_json: args.no_json,
            },
            crate::application::quantify::CompareQuantifyMode {
                preserve_missing_nocall_bd: args.usefiltered_truth,
                ..Default::default()
            },
        )?
    } else {
        let roc_options = crate::engines::roc::RocOptions {
            threads: comparison_threads,
            output_rocs: !args.roc.no_roc,
            whole_reference_size: Some(whole_reference_size),
            preserve_raw_table: args.verbose,
            ..Default::default()
        };
        let indices = roc_publication::write_roc_files_with_options_iter(
            prefix,
            row_file.rows()?,
            subset_size,
            conf_size,
            &roc_options,
        )?;
        if args.roc.no_roc {
            compact_no_roc_outputs(prefix)?;
        }
        indices
    };
    report_phase("roc_generation", roc_started);
    let commandline = std::env::args().collect::<Vec<_>>().join(" ");
    write_runinfo_for_args(&args, explicit_bcf, prefix, &commandline, published_logfile)?;
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

#[cfg(test)]
fn retain_xcmp_truth_calls(variants: &mut Vec<Variant>) {
    variants.retain(Variant::is_pass);
}

#[cfg(test)]
fn sort_comparison_rows(rows: &mut [AnnotatedRow], filtered_truth_keys: &BTreeSet<VariantKey>) {
    rows.sort_by(|left, right| {
        left.sort_key
            .chrom
            .cmp(&right.sort_key.chrom)
            .then(left.sort_key.pos.cmp(&right.sort_key.pos))
            // A query record sharing the exact key of a filtered truth call
            // occupies that merged record's slot in legacy xcmp, even though
            // the truth sample itself is excluded from comparison. It sorts
            // before other records at the locus (PASS truth included).
            .then_with(|| {
                let left_filtered = row_matches_variant_key(left, filtered_truth_keys);
                let right_filtered = row_matches_variant_key(right, filtered_truth_keys);
                right_filtered.cmp(&left_filtered)
            })
            .then(left.sort_key.side_rank.cmp(&right.sort_key.side_rank))
            .then(left.sort_key.type_rank.cmp(&right.sort_key.type_rank))
            .then_with(|| {
                (
                    &left.record.ref_allele,
                    &left.record.alt_allele,
                    &left.record.filter,
                    &left.record.info,
                    &left.record.samples,
                )
                    .cmp(&(
                        &right.record.ref_allele,
                        &right.record.alt_allele,
                        &right.record.filter,
                        &right.record.info,
                        &right.record.samples,
                    ))
            })
    });
}

fn row_matches_variant_key(row: &AnnotatedRow, keys: &BTreeSet<VariantKey>) -> bool {
    keys.contains(&VariantKey {
        chrom: row.record.chrom.clone(),
        pos: row.record.pos,
        ref_allele: row.record.ref_allele.clone(),
        alt_allele: row.record.alt_allele.clone(),
    })
}

fn run_vcfeval(
    args: &CompareArgs,
    explicit_bcf: bool,
    truth_prep: &Path,
    query_prep: &Path,
    prefix: &Path,
    scratch: ScratchRun,
    published_logfile: Option<&str>,
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
    let validated = comparison_io::run_vcfeval(
        truth_prep,
        query_prep,
        Path::new(&args.reference),
        crate::engines::vcfeval::Options {
            roc_field: &args.roc.roc,
            loose_match_distance: args.engine.engine_scmp_distance,
        },
    )?;
    let (vcfeval_headers, vcfeval_records) = validated.into_parts();
    let roc_indices = crate::application::quantify::run_from_compare(
        QuantifyArgs {
            input_vcf: String::new(),
            report_prefix: prefix.to_string_lossy().into_owned(),
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
            roc: args.roc.roc.clone(),
            do_roc: !args.roc.no_roc,
            roc_regions: args.roc.roc_regions.clone(),
            roc_filter: args.roc.roc_filter.clone(),
            roc_delta: args.roc.roc_delta,
            ci_alpha: args.ci_alpha,
            no_json: args.no_json,
        },
        vcfeval_headers,
        vcfeval_records,
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
    write_runinfo_for_args(args, explicit_bcf, prefix, &commandline, published_logfile)?;
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
    published_logfile: Option<&str>,
) -> Result<()> {
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
    let mode = match args.engine.engine {
        CompareEngine::ScmpSomatic => crate::engines::scmp::ScmpMode::Alleles,
        CompareEngine::ScmpDistance => crate::engines::scmp::ScmpMode::Distance {
            max_distance: i64::try_from(args.engine.engine_scmp_distance)
                .context("SCMP match distance exceeds the supported range")?,
        },
        CompareEngine::Xcmp | CompareEngine::Vcfeval => {
            bail!("internal error: run_scmp called for a non-SCMP engine")
        }
    };
    let validated = comparison_io::run_scmp(
        truth_prep,
        query_prep,
        Path::new(&args.reference),
        mode,
        &args.roc.roc,
    )?;
    let (scmp_headers, scmp_records) = validated.into_parts();
    let write_counts = args.write_counts && !args.no_write_counts;
    let roc_indices = crate::application::quantify::run_from_compare(
        QuantifyArgs {
            input_vcf: String::new(),
            report_prefix: prefix.to_string_lossy().into_owned(),
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
            roc: args.roc.roc.clone(),
            do_roc: !args.roc.no_roc,
            roc_regions: args.roc.roc_regions.clone(),
            roc_filter: args.roc.roc_filter.clone(),
            roc_delta: args.roc.roc_delta,
            ci_alpha: args.ci_alpha,
            no_json: args.no_json,
        },
        scmp_headers,
        scmp_records,
        crate::application::quantify::CompareQuantifyMode {
            preserve_missing_query_qq: args.engine.engine == CompareEngine::ScmpSomatic,
            ..Default::default()
        },
    )?;
    let commandline = std::env::args().collect::<Vec<_>>().join(" ");
    write_runinfo_for_args(args, explicit_bcf, prefix, &commandline, published_logfile)?;
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
    match args.engine.engine {
        CompareEngine::ScmpSomatic => {
            if !args.somatic && args.preprocess.set_gt.is_none() {
                args.somatic = true;
            }
            // Legacy turns partial-credit normalization off for the somatic
            // scmp engine after selecting the synthetic half genotype.
            args.preprocess_truth = false;
            args.preprocess.leftshift = false;
            args.preprocess.no_leftshift = true;
            args.preprocess.decompose = false;
            args.preprocess.bcftools_norm = false;
        }
        CompareEngine::ScmpDistance => {
            if !args.somatic && args.preprocess.set_gt.is_none() {
                args.preprocess.set_gt = Some(SomaticGtMode::First);
            }
            args.preprocess.decompose = false;
        }
        CompareEngine::Xcmp | CompareEngine::Vcfeval => {}
    }
}

fn somatic_mode_name(mode: Option<SomaticGtMode>) -> Option<&'static str> {
    mode.map(|mode| match mode {
        SomaticGtMode::Half => "half",
        SomaticGtMode::Hemi => "hemi",
        SomaticGtMode::Het => "het",
        SomaticGtMode::Hom => "hom",
        SomaticGtMode::First => "first",
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
    published_logfile: Option<&str>,
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
        preprocessing_leftshift: args.engine.engine == CompareEngine::ScmpSomatic
            || !args.preprocess.no_leftshift,
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
        do_roc: !args.roc.no_roc,
        engine: args.engine.engine.legacy_name(),
        engine_scmp_distance: args.engine.engine_scmp_distance,
        engine_vcfeval: args.engine.engine_vcfeval.as_deref().unwrap_or("rtg"),
        engine_vcfeval_template: args.engine.engine_vcfeval_template.as_deref(),
        filter_nonref: args.preprocess.filter_nonref,
        filters_only: args.filters_only.as_deref(),
        fixchr: if args.preprocess.no_fixchr {
            Some(false)
        } else {
            args.preprocess.fixchr
        },
        fp_adjust_conf: args.adjust_conf_regions && !args.no_adjust_conf_regions,
        gender: match args.preprocess.gender {
            PreprocessGender::Male => "male",
            PreprocessGender::Female => "female",
            PreprocessGender::Auto => "auto",
            PreprocessGender::None => "none",
        },
        hb_expand: args.engine.hb_expand,
        logfile: published_logfile,
        max_enum: args.engine.max_enum,
        no_hc: args.no_hc,
        output_vtc: args.output_vtc,
        preprocess_window: args.preprocess.preprocess_window,
        preprocessing_norm: args.preprocess.bcftools_norm,
        preserve_info: args.preserve_info,
        quiet: args.quiet,
        roc: &args.roc.roc,
        roc_delta: args.roc.roc_delta,
        roc_filter: args.roc.roc_filter.as_deref(),
        roc_regions: &args.roc.roc_regions,
        somatic: args.somatic,
        somatic_mode: somatic_mode_name(args.preprocess.set_gt),
        strat_fixchr: args.strat_fixchr,
        strat_regions: &args.strat_regions,
        usefiltered_truth: args.usefiltered_truth,
        verbose: args.verbose,
        window: args.engine.window,
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
        fixchr: args.preprocess.fixchr,
        no_fixchr: args.preprocess.no_fixchr,
        somatic: args.somatic,
        set_gt: args.preprocess.set_gt,
        filter_nonref: args.preprocess.filter_nonref && (!truth_side || preprocess_enabled),
        convert_gvcf_to_vcf: convert_gvcf,
        bcf: args.bcf,
        bcftools_norm: args.preprocess.bcftools_norm && (!truth_side || preprocess_enabled),
        leftshift: preprocess_enabled
            && (args.preprocess.leftshift || !args.preprocess.no_leftshift),
        no_leftshift: false,
        decompose: preprocess_enabled && effective_decomposition(args),
        no_decompose: false,
        gender: args.preprocess.gender,
        window_size: args.preprocess.preprocess_window as i64,
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
    !args.preprocess.no_decompose
        && (!(args.somatic || args.preprocess.set_gt.is_some()) || args.preprocess.decompose)
}

#[cfg(test)]
mod test_suite;
