use crate::cli::{CompareArgs, CompareEngine, PreprocessArgs};
use crate::fasta;
use crate::metrics_json;
use crate::partial_credit;
use crate::preprocess;
use crate::report::{self, CountsBucket};
use crate::vcf::{self, RawVcfRecord, Variant, VariantKey};
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::borrow::Borrow;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

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
    /// Legacy xcmp annotates every record in a superlocus with the result of
    /// the block-level comparison. Kept out-of-band until final decoration so
    /// ordinary (clean-INFO) output remains unchanged.
    pub xcmp_ctype: Option<&'static str>,
    pub xcmp_hap_match: bool,
}

/// Append one of the legacy report suffixes without treating a dotted report
/// prefix as a filename extension.
pub(crate) fn suffixed_report_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(".");
    path.push(suffix);
    PathBuf::from(path)
}

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

pub fn run(mut args: CompareArgs) -> Result<()> {
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
    if preprocessing_args.gender == crate::cli::PreprocessGender::Auto {
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
        let mut preprocessed_truth = vcf::open_raw_vcf(&truth_prep)?;
        let mut common_contig = false;
        for record in &mut preprocessed_truth {
            if contig_set.contains(&record?.chrom) {
                common_contig = true;
                break;
            }
        }
        if !common_contig {
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

    let truth_reader = vcf::open_raw_vcf(&truth_prep)?;
    let truth_headers = truth_reader.headers().to_vec();
    let query_reader = vcf::open_raw_vcf(&query_prep)?;
    let query_headers = query_reader.headers().to_vec();
    drop(truth_reader);
    drop(query_reader);

    let filtered_truth_keys = vcf::open_variants(
        &truth_prep,
        &contig_set,
        false,
        regions.as_deref(),
        targets.as_deref(),
        locations.as_deref(),
    )?
    .filter_map(|variant| match variant {
        Ok(variant) if !variant.is_pass() => Some(Ok(variant.key)),
        Ok(_) => None,
        Err(error) => Some(Err(error)),
    })
    .collect::<Result<BTreeSet<_>>>()?;
    // CONF insertion padding is derived from the preprocessed truth stream,
    // including filtered records retained by `--usefiltered-truth`. xcmp
    // excludes those records as calls below, but gvcf2bed sees them first;
    // collisions with PASS records can therefore change the padding lane.
    let cluster_gap = match args.engine {
        CompareEngine::ScmpSomatic => 0,
        CompareEngine::ScmpDistance => args.engine_scmp_distance,
        CompareEngine::Xcmp | CompareEngine::Vcfeval => args.window,
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
    let gvcf_padding: Option<Vec<vcf::BedInterval>> = if let Some(raw) = conf_bed.as_deref() {
        Some(if adjust_conf {
            gvcf2bed_padding_iter(
                vcf::open_variants(
                    &truth_prep,
                    &contig_set,
                    false,
                    regions.as_deref(),
                    targets.as_deref(),
                    locations.as_deref(),
                )?,
                Some(raw),
            )?
        } else {
            Vec::new()
        })
    } else {
        None
    };
    let adjusted_conf_bed: Option<Vec<vcf::BedInterval>> = conf_bed.as_ref().map(|raw| {
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
        locations.as_deref(),
    )?;
    let query = vcf::open_variants(
        &query_prep,
        &contig_set,
        false,
        regions.as_deref(),
        targets.as_deref(),
        locations.as_deref(),
    )?;
    let clusters = StreamingClusters::new(truth, query, cluster_gap);
    let mut contigs_in_play = locations
        .as_deref()
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
    let mut counts: BTreeMap<String, TypeCounts> = BTreeMap::new();
    let mut subtype_counts: BTreeMap<String, BTreeMap<String, TypeCounts>> = BTreeMap::new();
    let mut rows = Vec::new();
    for cluster in clusters {
        let cluster = cluster?;
        contigs_in_play.insert(cluster.chrom.clone());
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

    let subset_size = report_subset_size(
        &contig_non_n_lengths,
        &contigs_in_play,
        explicit_bcf,
        args.bcf && !explicit_bcf,
    );
    if subset_size == 0 {
        bail!("no reference contigs selected for analysis");
    }

    sort_comparison_rows(&mut rows, &filtered_truth_keys);
    decorate_output_rows_from_paths(
        &mut rows,
        [&truth_prep, &query_prep],
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
        crate::quantify::run_from_compare(
            crate::cli::QuantifyArgs {
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
                    .then(|| truth_prep.display().to_string()),
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
            crate::quantify::CompareQuantifyMode {
                preserve_missing_nocall_bd: args.usefiltered_truth,
                ..Default::default()
            },
        )?
    } else {
        let indices = crate::roc::write_roc_files(prefix, &rows, subset_size, conf_size)?;
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
            crate::cli::PreprocessGender::Male => "male",
            crate::cli::PreprocessGender::Female => "female",
            crate::cli::PreprocessGender::Auto => "auto",
            crate::cli::PreprocessGender::None => "none",
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
    let reader = vcf::open_raw_vcf(&vcf_path)
        .with_context(|| format!("failed to prepare BCF report from {}", vcf_path.display()))?;
    let headers = reader.headers().to_vec();
    vcf::write_raw_vcf_iter(&bcf_path, &headers, reader)?;
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
            .then_with(|| left.line.cmp(&right.line))
    });
}

fn row_matches_variant_key(row: &AnnotatedRow, keys: &BTreeSet<VariantKey>) -> bool {
    let mut fields = row.line.split('\t');
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
        let truth = vcf::open_raw_vcf(truth_prep)?.map(|record| -> Result<Variant> {
            let record = record?;
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
        });
        let padding = gvcf2bed_padding_iter(truth, Some(&raw_conf))?;
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
    crate::vcfeval::compare_files(
        truth_prep,
        query_prep,
        Path::new(&args.reference),
        crate::vcfeval::Options {
            roc_field: &args.roc,
            loose_match_distance: args.engine_scmp_distance,
        },
        &vcfeval_vcf,
    )?;
    let roc_indices = crate::quantify::run_from_compare(
        crate::cli::QuantifyArgs {
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
        crate::quantify::CompareQuantifyMode {
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
        let truth = vcf::open_variants(truth_prep, &contigs, false, None, None, None)?;
        let padding = gvcf2bed_padding_iter(truth, Some(&raw_conf))?;
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
        CompareEngine::ScmpSomatic => crate::scmp::ScmpMode::Alleles,
        CompareEngine::ScmpDistance => crate::scmp::ScmpMode::Distance {
            max_distance: i64::try_from(args.engine_scmp_distance)
                .context("SCMP match distance exceeds the supported range")?,
        },
        CompareEngine::Xcmp | CompareEngine::Vcfeval => {
            bail!("internal error: run_scmp called for a non-SCMP engine")
        }
    };
    crate::scmp::compare_files(
        truth_prep,
        query_prep,
        Path::new(&args.reference),
        mode,
        &args.roc,
        &comparison_vcf,
    )?;
    let write_counts = args.write_counts && !args.no_write_counts;
    let roc_indices = crate::quantify::run_from_compare(
        crate::cli::QuantifyArgs {
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
        crate::quantify::CompareQuantifyMode {
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
                args.set_gt = Some(crate::cli::SomaticGtMode::First);
            }
            args.decompose = false;
        }
        CompareEngine::Xcmp | CompareEngine::Vcfeval => {}
    }
}

fn somatic_mode_name(mode: Option<crate::cli::SomaticGtMode>) -> Option<&'static str> {
    mode.map(|mode| match mode {
        crate::cli::SomaticGtMode::Half => "half",
        crate::cli::SomaticGtMode::Hemi => "hemi",
        crate::cli::SomaticGtMode::Het => "het",
        crate::cli::SomaticGtMode::Hom => "hom",
        crate::cli::SomaticGtMode::First => "first",
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
            crate::cli::PreprocessGender::Male => "male",
            crate::cli::PreprocessGender::Female => "female",
            crate::cli::PreprocessGender::Auto => "auto",
            crate::cli::PreprocessGender::None => "none",
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

fn rewrite_compare_metrics(
    args: &CompareArgs,
    prefix: &Path,
    commandline: &str,
    roc_indices: &crate::roc::MetricIndices,
) -> Result<()> {
    if args.no_json {
        return Ok(());
    }
    let mut tables = vec![(
        "summary.metrics",
        "summary.metrics",
        suffixed_report_path(prefix, "summary.csv"),
    )];
    if args.write_counts && !args.no_write_counts {
        tables.push((
            "all.metrics",
            "all.metrics",
            suffixed_report_path(prefix, "extended.csv"),
        ));
    }
    for id in &roc_indices.table_order {
        let path = suffixed_report_path(prefix, &format!("{id}.csv.gz"));
        if path.is_file() {
            tables.push((id.as_str(), id.as_str(), path));
        }
    }
    let refs = tables
        .iter()
        .map(|(id, label, path)| (*id, *label, path.as_path()))
        .collect::<Vec<_>>();
    metrics_json::write_metrics_gz_for_module_with_indices(
        &suffixed_report_path(prefix, "metrics.json.gz"),
        "hap.py.comparison",
        "hap.py",
        commandline,
        &refs,
        Some(&roc_indices.tables),
    )
}

type InfoKey = (String, usize, String, String);
type SemanticInfoKey = (String, usize, String, Vec<String>);

fn semantic_info_key(record: &RawVcfRecord) -> SemanticInfoKey {
    let mut alts = record
        .alt_allele
        .split(',')
        .map(str::to_string)
        .collect::<Vec<_>>();
    alts.sort();
    let reference = if alts
        .iter()
        .any(|alt| alt.starts_with('<') && alt.ends_with('>'))
    {
        "*".to_string()
    } else {
        record.ref_allele.clone()
    };
    (record.chrom.clone(), record.pos, reference, alts)
}

#[cfg(test)]
fn decorate_output_rows(
    rows: &mut [AnnotatedRow],
    truth: &[RawVcfRecord],
    query: &[RawVcfRecord],
    preserve_info: bool,
    output_vtc: bool,
    roc_field: &str,
) -> Result<()> {
    if !preserve_info && !output_vtc && matches!(roc_field, "QUAL" | "QQ") {
        return Ok(());
    }
    let mut decorations = DecorationIndex::default();
    for record in truth.iter().chain(query) {
        decorations.observe(record, preserve_info, roc_field);
    }
    decorate_output_rows_with_index(rows, &decorations, preserve_info, output_vtc, roc_field)
}

fn decorate_output_rows_from_paths<const N: usize>(
    rows: &mut [AnnotatedRow],
    paths: [&Path; N],
    preserve_info: bool,
    output_vtc: bool,
    roc_field: &str,
) -> Result<()> {
    if !preserve_info && !output_vtc && matches!(roc_field, "QUAL" | "QQ") {
        return Ok(());
    }
    let mut decorations = DecorationIndex::default();
    for path in paths {
        for record in vcf::open_raw_vcf(path)? {
            decorations.observe(&record?, preserve_info, roc_field);
        }
    }
    decorate_output_rows_with_index(rows, &decorations, preserve_info, output_vtc, roc_field)
}

#[derive(Default)]
struct DecorationIndex {
    preserved: BTreeMap<InfoKey, BTreeSet<String>>,
    semantic_preserved: BTreeMap<SemanticInfoKey, BTreeSet<String>>,
    roc_values: BTreeMap<InfoKey, String>,
}

impl DecorationIndex {
    fn observe(&mut self, record: &RawVcfRecord, preserve_info: bool, roc_field: &str) {
        let key = (
            record.chrom.clone(),
            record.pos,
            record.ref_allele.clone(),
            record.alt_allele.clone(),
        );
        if preserve_info {
            for field in record
                .info
                .split(';')
                .filter(|field| !matches!(*field, "" | "."))
            {
                self.preserved
                    .entry(key.clone())
                    .or_default()
                    .insert(field.to_string());
                self.semantic_preserved
                    .entry(semantic_info_key(record))
                    .or_default()
                    .insert(field.to_string());
            }
        }
        if !matches!(roc_field, "QUAL" | "QQ" | ".") {
            let value = record
                .info
                .split(';')
                .find_map(|entry| {
                    entry
                        .split_once('=')
                        .filter(|(key, _)| *key == roc_field)
                        .map(|(_, value)| value.to_string())
                })
                .or_else(|| record.sample_map(0).get(roc_field).cloned());
            if let Some(value) = value {
                self.roc_values.insert(key, value);
            }
        }
    }
}

fn decorate_output_rows_with_index(
    rows: &mut [AnnotatedRow],
    decorations: &DecorationIndex,
    preserve_info: bool,
    output_vtc: bool,
    roc_field: &str,
) -> Result<()> {
    for row in rows {
        let mut record = RawVcfRecord::from_line(&row.line, Path::new("comparison-output"))?;
        let key = (
            record.chrom.clone(),
            record.pos,
            record.ref_allele.clone(),
            record.alt_allele.clone(),
        );
        let current = info_fields_by_key(&record.info);
        let regions = current.get("Regions").cloned();
        let mut base = BTreeMap::<String, String>::new();
        if preserve_info {
            let fields = decorations.preserved.get(&key).into_iter().chain(
                decorations
                    .semantic_preserved
                    .get(&semantic_info_key(&record)),
            );
            for fields in fields {
                for field in fields {
                    let field_key = field.split_once('=').map_or(field.as_str(), |(key, _)| key);
                    if field_key != "Regions" {
                        base.insert(field_key.to_string(), field.clone());
                    }
                }
            }
        }
        if let Some(bs) = current.get("BS") {
            base.insert("BS".to_string(), bs.clone());
        }
        let truth_fields = record.sample_map(0);
        let query_fields = record.sample_map(1);
        let comparison = legacy_comparison_fields(
            &record,
            &truth_fields,
            &query_fields,
            row.xcmp_ctype.unwrap_or("simple:match"),
        );
        if preserve_info {
            for (name, value) in comparison.preserved_fields(row.xcmp_hap_match) {
                base.insert(name.to_string(), value);
            }
            let iqq = if roc_field == "QUAL" {
                Some(comparison.iqq.as_str())
            } else {
                decorations.roc_values.get(&key).map(String::as_str)
            };
            if let Some(iqq) = iqq {
                base.insert("IQQ".to_string(), format!("IQQ={iqq}"));
            }
        }

        if !matches!(roc_field, "QUAL" | "QQ" | ".") && !decorations.roc_values.contains_key(&key) {
            // When xcmp's literal custom-field lookup misses, quantify reads
            // an absent IQQ: called query samples receive NaN, truth samples
            // remain missing, and no-call queries retain the zero sentinel.
            set_comparison_format_value(&mut record, 0, "QQ", ".");
            let query_called = record
                .sample_map(1)
                .get("BVT")
                .is_some_and(|value| !matches!(value.as_str(), "" | "." | "NOCALL"));
            set_comparison_format_value(
                &mut record,
                1,
                "QQ",
                if query_called { "nan" } else { "." },
            );
        }

        let mut info = base.into_values().collect::<Vec<_>>();
        let append_regions_last = regions.as_deref() == Some("Regions=TS_boundary");
        if !append_regions_last && let Some(regions) = regions.as_ref() {
            info.push(regions.clone());
        }
        if preserve_info {
            info.push(format!("RegionsExtent={}", legacy_regions_extent(&record)));
        }
        if output_vtc {
            let mut xcmp_type = comparison.decision;
            let mut xcmp_kind = comparison.kind.to_string();
            if row.xcmp_hap_match && xcmp_type != "TP" {
                xcmp_kind = format!("hapmatch__{xcmp_type}__{xcmp_kind}");
                xcmp_type = "TP".to_string();
            }
            if !info_list_values(&record.info, "Regions").contains(&"CONF") {
                xcmp_type = "UNK".to_string();
            }
            info.push(format!(
                "XCMP={xcmp_type}:{xcmp_kind}:{}:{}:{}",
                comparison.gtt1, comparison.gtt2, comparison.ctype
            ));
            let vtc = legacy_vtc(&record, &truth_fields, &query_fields);
            if !vtc.is_empty() {
                info.push(format!("VTC={vtc}"));
            }
        }
        if append_regions_last && let Some(regions) = regions {
            info.push(regions);
        }
        record.info = if info.is_empty() {
            ".".to_string()
        } else {
            info.join(";")
        };
        row.line = record.to_line();
    }
    Ok(())
}

fn info_fields_by_key(info: &str) -> BTreeMap<String, String> {
    info.split(';')
        .filter(|field| !matches!(*field, "" | "."))
        .map(|field| {
            let key = field.split_once('=').map_or(field, |(key, _)| key);
            (key.to_string(), field.to_string())
        })
        .collect()
}

/// Internal compare rows have already received synthetic truth-set membership
/// tags. Legacy hands qfy the pre-quantification stream instead, so remove
/// those provisional tags from the private re-quantification handoff and let
/// qfy derive them from the final confidence and stratification inputs.
fn sanitize_requantify_handoff_rows(rows: &[AnnotatedRow]) -> Vec<AnnotatedRow> {
    rows.iter()
        .cloned()
        .map(|mut row| {
            let mut fields = row.line.split('\t').map(str::to_string).collect::<Vec<_>>();
            if let Some(info) = fields.get_mut(7) {
                let entries = info
                    .split(';')
                    .filter_map(|entry| {
                        let Some(regions) = entry.strip_prefix("Regions=") else {
                            return Some(entry.to_string());
                        };
                        let retained = regions
                            .split(',')
                            .filter(|tag| !matches!(*tag, "TS_boundary" | "TS_contained"))
                            .collect::<Vec<_>>();
                        (!retained.is_empty()).then(|| format!("Regions={}", retained.join(",")))
                    })
                    .collect::<Vec<_>>();
                *info = if entries.is_empty() {
                    ".".to_string()
                } else {
                    entries.join(";")
                };
                row.line = fields.join("\t");
            }
            row
        })
        .collect()
}

fn set_comparison_format_value(
    record: &mut RawVcfRecord,
    sample_index: usize,
    key: &str,
    value: &str,
) {
    let Some(index) = record.format_keys().iter().position(|field| *field == key) else {
        return;
    };
    let Some(sample) = record.samples.get_mut(sample_index) else {
        return;
    };
    let mut fields = sample.split(':').map(str::to_string).collect::<Vec<_>>();
    if let Some(field) = fields.get_mut(index) {
        *field = value.to_string();
        *sample = fields.join(":");
    }
}

struct LegacyComparison {
    decision: String,
    kind: String,
    gtt1: String,
    gtt2: String,
    ctype: &'static str,
    iqq: String,
}

impl LegacyComparison {
    fn preserved_fields(&self, hap_match: bool) -> Vec<(&'static str, String)> {
        let mut fields = vec![
            ("ctype", format!("ctype={}", self.ctype)),
            ("kind", format!("kind={}", self.kind)),
            ("type", format!("type={}", self.decision)),
        ];
        if self.gtt1 != "." {
            fields.push(("gtt1", format!("gtt1={}", self.gtt1)));
        }
        if self.gtt2 != "." {
            fields.push(("gtt2", format!("gtt2={}", self.gtt2)));
        }
        if hap_match {
            fields.push(("HapMatch", "HapMatch".to_string()));
        }
        fields
    }
}

fn legacy_comparison_fields(
    record: &RawVcfRecord,
    truth: &BTreeMap<String, String>,
    query: &BTreeMap<String, String>,
    ctype: &'static str,
) -> LegacyComparison {
    let truth_called = sample_is_called(truth);
    let query_called = sample_is_called(query);
    let (decision, kind) = match (truth_called, query_called) {
        (true, false) => ("FN", "missing"),
        (false, true) => ("FP", "missing"),
        (false, false) => ("N", "match"),
        (true, true) => {
            let truth_bd = truth.get("BD").map(String::as_str).unwrap_or(".");
            let query_bd = query.get("BD").map(String::as_str).unwrap_or(".");
            let bk = query
                .get("BK")
                .or_else(|| truth.get("BK"))
                .map(String::as_str)
                .unwrap_or(".");
            if truth_bd == "FN" || query_bd == "FP" || bk == "am" || bk == "lm" {
                ("FP", legacy_mismatch_kind(record, truth, query))
            } else {
                ("TP", "match")
            }
        }
    };
    let iqq = query
        .get("QQ")
        .filter(|value| !matches!(value.as_str(), "" | "."))
        .cloned()
        .unwrap_or_else(|| "0".to_string());
    LegacyComparison {
        decision: decision.to_string(),
        kind: kind.to_string(),
        gtt1: legacy_gt_label(truth, truth_called),
        gtt2: legacy_gt_label(query, query_called),
        ctype,
        iqq,
    }
}

fn sample_is_called(sample: &BTreeMap<String, String>) -> bool {
    sample
        .get("BVT")
        .is_some_and(|value| !matches!(value.as_str(), "" | "." | "NOCALL" | "HOMREF"))
}

fn legacy_gt_label(sample: &BTreeMap<String, String>, called: bool) -> String {
    if !called {
        return ".".to_string();
    }
    sample
        .get("BLT")
        .filter(|value| !matches!(value.as_str(), "" | "." | "nocall"))
        .map(|value| format!("gt_{value}"))
        .unwrap_or_else(|| "gt_unknown".to_string())
}

fn legacy_mismatch_kind<'a>(
    _record: &RawVcfRecord,
    truth: &'a BTreeMap<String, String>,
    query: &'a BTreeMap<String, String>,
) -> &'static str {
    let truth_gt = truth.get("GT").map(String::as_str).unwrap_or(".");
    let query_gt = query.get("GT").map(String::as_str).unwrap_or(".");
    let truth_alleles = parse_gt_alleles(truth_gt)
        .into_iter()
        .collect::<BTreeSet<_>>();
    let query_alleles = parse_gt_alleles(query_gt)
        .into_iter()
        .collect::<BTreeSet<_>>();
    if truth_alleles == query_alleles {
        return "gtmismatch";
    }
    let truth_nonref = truth_alleles
        .into_iter()
        .filter(|allele| *allele > 0)
        .collect::<BTreeSet<_>>();
    let query_nonref = query_alleles
        .into_iter()
        .filter(|allele| *allele > 0)
        .collect::<BTreeSet<_>>();
    if truth_nonref == query_nonref {
        "gtmismatch"
    } else if !truth_nonref.is_disjoint(&query_nonref) {
        "alpartial"
    } else {
        "almismatch"
    }
}

fn legacy_regions_extent(record: &RawVcfRecord) -> String {
    let variant = Variant {
        key: VariantKey {
            chrom: record.chrom.clone(),
            pos: record.pos,
            ref_allele: record.ref_allele.clone(),
            alt_allele: record.alt_allele.clone(),
        },
        qual: record.qual.clone(),
        filter: record.filter.clone(),
        gt: ".".to_string(),
    };
    effective_refrange(&variant)
        .map(|(start, end, _)| format!("{start}-{end}"))
        .unwrap_or_else(|| {
            format!(
                "{}-{}",
                record.pos,
                record.pos + record.ref_allele.len().saturating_sub(1)
            )
        })
}

fn legacy_vtc(
    record: &RawVcfRecord,
    truth: &BTreeMap<String, String>,
    query: &BTreeMap<String, String>,
) -> String {
    let mut types = BTreeMap::<u8, String>::new();
    for sample in [truth, query] {
        if !sample_is_called(sample) {
            types.insert(0x80, "nocall__nc".to_string());
            continue;
        }
        let gt = sample.get("GT").map(String::as_str).unwrap_or(".");
        let gt_alleles = parse_gt_alleles(gt);
        let mut allele_bits = 0u8;
        for allele in gt_alleles.iter().copied().filter(|allele| *allele > 0) {
            let Some(alt) = record.alt_allele.split(',').nth(allele - 1) else {
                continue;
            };
            let bits = allele_edit_bits(&record.ref_allele, alt);
            allele_bits |= bits;
            for bit in [1u8, 2, 4] {
                if bits & bit != 0 {
                    types.insert(bit, format!("nuc__{}", legacy_type_bits(bit)));
                }
            }
            if bits != 0 {
                types.insert(0x10 | bits, format!("al__{}", legacy_type_bits(bits)));
            }
        }
        if allele_bits == 0 {
            continue;
        }
        let location = match sample.get("BLT").map(String::as_str).unwrap_or("") {
            "het" => 0x30,
            "hetalt" => 0x40,
            "homalt" => 0x90,
            "hemi" => 0x50,
            _ => 0xa0,
        };
        let ref_bit = u8::from(gt_alleles.contains(&0)) * 8;
        types.insert(
            location | ref_bit | allele_bits,
            format!(
                "{}__{}",
                sample.get("BLT").map(String::as_str).unwrap_or("unknown"),
                legacy_type_bits(ref_bit | allele_bits)
            ),
        );
    }
    types.into_values().collect::<Vec<_>>().join(",")
}

fn allele_edit_bits(reference: &str, alternate: &str) -> u8 {
    if alternate.starts_with('<') {
        return if alternate.starts_with("<DEL") { 4 } else { 2 };
    }
    let ref_bytes = reference.as_bytes();
    let alt_bytes = alternate.as_bytes();
    let prefix = ref_bytes
        .iter()
        .zip(alt_bytes)
        .take_while(|(left, right)| left == right)
        .count();
    let suffix_limit = (ref_bytes.len() - prefix).min(alt_bytes.len() - prefix);
    let suffix = (0..suffix_limit)
        .take_while(|offset| {
            ref_bytes[ref_bytes.len() - 1 - offset] == alt_bytes[alt_bytes.len() - 1 - offset]
        })
        .count();
    let ref_remaining = ref_bytes.len() - prefix - suffix;
    let alt_remaining = alt_bytes.len() - prefix - suffix;
    match (ref_remaining, alt_remaining) {
        (0, 0) => 0,
        (0, _) => 2,
        (_, 0) => 4,
        (left, right) if left == right => 1,
        (left, right) if left < right => 1 | 2,
        _ => 1 | 4,
    }
}

fn legacy_type_bits(bits: u8) -> &'static str {
    const NAMES: [&str; 16] = [
        "nc", "s", "i", "si", "d", "sd", "id", "sid", "r", "rs", "ri", "rsi", "rd", "rsd", "rid",
        "rsid",
    ];
    NAMES[usize::from(bits & 0x0f)]
}

fn decorate_existing_comparison_vcf(
    output_path: &Path,
    truth_path: &Path,
    query_path: &Path,
    preserve_info: bool,
    output_vtc: bool,
) -> Result<()> {
    let reader = vcf::open_raw_vcf(output_path)?;
    let mut headers = reader.headers().to_vec();
    if output_vtc {
        let mut chrom_index = headers
            .iter()
            .position(|line| line.starts_with("#CHROM"))
            .unwrap_or(headers.len());
        for declaration in [
            "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">",
            "##INFO=<ID=XCMP,Number=.,Type=String,Description=\"XCMP extra information.\">",
        ] {
            let identity = preprocess::structured_header_identity(declaration);
            let present = headers
                .iter()
                .any(|line| preprocess::structured_header_identity(line) == identity);
            if !present {
                headers.insert(chrom_index, declaration.to_string());
                chrom_index += 1;
            }
        }
    }
    let mut decorations = DecorationIndex::default();
    for path in [truth_path, query_path] {
        for record in vcf::open_raw_vcf(path)? {
            decorations.observe(&record?, preserve_info, "QUAL");
        }
    }
    let decorated = reader.map(|record| {
        let record = record?;
        let mut row = AnnotatedRow {
            sort_key: (record.chrom.clone(), record.pos, 0, 0),
            line: record.to_line(),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        };
        decorate_output_rows_with_index(
            std::slice::from_mut(&mut row),
            &decorations,
            preserve_info,
            output_vtc,
            "QUAL",
        )?;
        RawVcfRecord::from_line(&row.line, output_path)
    });
    vcf::write_raw_vcf_iter(output_path, &headers, decorated)
}

fn resolve_default_reference() -> Result<String> {
    for candidate in [
        std::env::var_os("HG19"),
        std::env::var_os("HGREF"),
        Some("/opt/hap.py-data/hg19.fa".into()),
    ]
    .into_iter()
    .flatten()
    {
        let path = PathBuf::from(candidate);
        if path.is_file() {
            return Ok(path.display().to_string());
        }
    }
    bail!("no reference file found; pass --reference or set HG19/HGREF")
}

fn build_vcf_headers(
    truth_headers: &[String],
    query_headers: &[String],
    apply_filters_query: bool,
    output_vtc: bool,
    preserve_info: bool,
    roc_field: &str,
) -> Vec<String> {
    let mut supplied: Vec<String> = truth_headers
        .iter()
        .chain(query_headers)
        .filter(|line| line.starts_with("##"))
        .cloned()
        .collect();
    supplied.extend([
        "##INFO=<ID=gtt1,Number=1,Type=String,Description=\"GT of truth call\">".to_string(),
        "##INFO=<ID=gtt2,Number=1,Type=String,Description=\"GT of query call\">".to_string(),
        "##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">".to_string(),
        "##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">".to_string(),
        "##INFO=<ID=ctype,Number=1,Type=String,Description=\"Type of comparison performed\">".to_string(),
        "##INFO=<ID=HapMatch,Number=0,Type=Flag,Description=\"Variant is in matching haplotype block\">".to_string(),
        "##INFO=<ID=BS,Number=1,Type=Integer,Description=\"Start position of the benchmarking superlocus on current chromosome\">".to_string(),
        format!("##INFO=<ID=IQQ,Number=1,Type=Float,Description=\"Quality value for query variants ({roc_field}).\">")
    ]);
    if apply_filters_query {
        supplied.push(
            "##INFO=<ID=Q_FILTERED,Number=0,Type=Flag,Description=\"Filtered call in query\">"
                .to_string(),
        );
    }
    supplied
        .push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string());

    // xcmp first writes through VariantWriter (base + sorted merged input
    // headers). Quantify then appends its region/FORMAT declarations in this
    // exact order before writing the final two-sample VCF.
    let mut headers = preprocess::canonicalize_legacy_headers(&supplied);
    let chrom = headers.pop().expect("comparison header has #CHROM line");
    let mut quantified_headers =
        vec!["##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">"];
    if preserve_info {
        quantified_headers.push(
            "##INFO=<ID=RegionsExtent,Number=.,Type=String,Description=\"Trimmed reference coordinates matched to regions for this record.\">",
        );
    }
    quantified_headers.extend([
        "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">",
        "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">",
        "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">",
        "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">",
        "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"High-level variant type (SNP|INDEL).\">",
        "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"High-level location type (het|homref|hetalt|homalt|nocall).\">",
    ]);
    if output_vtc {
        quantified_headers.extend([
            "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">",
            "##INFO=<ID=XCMP,Number=.,Type=String,Description=\"XCMP extra information.\">",
        ]);
    }
    for line in quantified_headers {
        let identity = preprocess::structured_header_identity(line);
        let already_present = identity.as_ref().is_some_and(|wanted| {
            headers.iter().any(|existing| {
                preprocess::structured_header_identity(existing).as_ref() == Some(wanted)
            })
        });
        if !already_present {
            headers.push(line.to_string());
        }
    }
    headers.push(chrom);
    headers
}

#[cfg(test)]
fn build_clusters(truth: &[Variant], query: &[Variant]) -> Vec<Cluster> {
    build_clusters_with_gap(truth, query, CLUSTER_GAP_BP)
}

#[cfg(test)]
fn build_clusters_with_gap(
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

struct StreamingClusters<T, Q>
where
    T: Iterator<Item = Result<Variant>>,
    Q: Iterator<Item = Result<Variant>>,
{
    truth: std::iter::Peekable<T>,
    query: std::iter::Peekable<Q>,
    pending: Option<Entry>,
    cluster_gap: usize,
}

impl<T, Q> StreamingClusters<T, Q>
where
    T: Iterator<Item = Result<Variant>>,
    Q: Iterator<Item = Result<Variant>>,
{
    fn new(truth: T, query: Q, cluster_gap: usize) -> Self {
        Self {
            truth: truth.peekable(),
            query: query.peekable(),
            pending: None,
            cluster_gap,
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
                variant_stream_key(truth) <= variant_stream_key(query)
            }
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => return Ok(None),
            (Some(_), Some(_)) => unreachable!("errors handled before ordering"),
        };
        let variant = if take_truth {
            self.truth.next().transpose()?.map(|variant| Entry {
                side: Side::Truth,
                variant,
            })
        } else {
            self.query.next().transpose()?.map(|variant| Entry {
                side: Side::Query,
                variant,
            })
        };
        Ok(variant)
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
                && cluster.truth.len() + cluster.query.len() < MAX_CLUSTER_VARIANTS
            {
                cluster.end = cluster.end.max(entry.variant.end_pos());
                push_cluster_entry(&mut cluster, entry);
            } else {
                self.pending = Some(entry);
                return Some(Ok(cluster));
            }
        }
    }
}

fn variant_stream_key(variant: &Variant) -> (&str, usize, usize) {
    (&variant.key.chrom, variant.key.pos, variant.end_pos())
}

fn push_cluster_entry(cluster: &mut Cluster, entry: Entry) {
    match entry.side {
        Side::Truth => cluster.truth.push(entry.variant),
        Side::Query => cluster.query.push(entry.variant),
    }
}

fn process_cluster(
    cluster: &Cluster,
    reference_sequences: &BTreeMap<String, String>,
    conf_bed: Option<&[vcf::BedInterval]>,
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
            if row.line.contains(":UNK:lm:") {
                row.line = row.line.replace(":UNK:lm:", ":UNK:.:");
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
                && row.line.contains(":FN:am:")
                && row.line.contains(":FP:am:")
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

fn identical_gt_exact_indel_keys(cluster: &Cluster) -> BTreeSet<VariantKey> {
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

fn degrade_identical_exact_unk_rows(rows: &mut [AnnotatedRow], keys: &BTreeSet<VariantKey>) {
    for row in rows {
        if row_matches_variant_key(row, keys) && row.line.contains(":UNK:lm:") {
            row.line = row.line.replace(":UNK:lm:", ":UNK:.:");
        }
    }
}

fn legacy_preprocessed_snp_first_positions(cluster: &Cluster) -> BTreeSet<usize> {
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

fn annotated_row_is_snp(row: &AnnotatedRow) -> bool {
    let fields = row.line.split('\t').collect::<Vec<_>>();
    fields.get(3).is_some_and(|reference| reference.len() == 1)
        && fields
            .get(4)
            .is_some_and(|alternate| alternate.split(',').all(|allele| allele.len() == 1))
}

fn legacy_repetitive_indel_hap_promotions(
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

fn same_single_alt_indel_edit(left: &Variant, right: &Variant) -> bool {
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

fn set_xcmp_context(rows: &mut [AnnotatedRow], ctype: &'static str, hap_match: bool) {
    for row in rows {
        row.xcmp_ctype = Some(ctype);
        row.xcmp_hap_match = hap_match;
    }
}

fn annotated_row_has_unreconciled_allele(row: &AnnotatedRow) -> bool {
    let fields = row.line.split('\t').collect::<Vec<_>>();
    if fields.len() < 11 {
        return false;
    }
    let format = fields[8].split(':').collect::<Vec<_>>();
    let Some(bk_index) = format.iter().position(|key| *key == "BK") else {
        return false;
    };
    for sample in [&fields[9], &fields[10]] {
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
fn cluster_signature(
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

fn cluster_signature_with_limit(
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
fn compute_paired_bk(truth: &Variant, query: &Variant) -> &'static str {
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
fn compute_class_c_relaxation_positions(query: &[Variant], truth: &[Variant]) -> BTreeSet<usize> {
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
fn estimated_state_count(variants: &[Variant]) -> usize {
    estimated_state_count_with_limit(variants, XCMP_ENUMERATION_THRESHOLD)
}

fn estimated_state_count_with_limit(variants: &[Variant], max_enum: usize) -> usize {
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
fn simple_compare_pairs_match(
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
fn query_primitive_splits(
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
fn compute_shift_target(
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
fn try_split_same_anchor_via_shift(
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
fn truth_subset_match(
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
fn remap_query_gt_subset(truth: &Variant, query: &Variant) -> String {
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

fn is_multi_allelic(variant: &Variant) -> bool {
    variant.key.alt_allele.contains(',')
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

fn has_nonconf_split_sibling(
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

fn matched_query_unk_bk(
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
fn fp_class_from_bk(bk: &str) -> Option<&'static str> {
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

/// Same physical VCF locus with the same declared alleles, allowing ALT
/// columns to use different index orders. Legacy's shared allele table uses
/// this equivalence when it emits a combined FN/FP genotype-mismatch row;
/// the query GT is remapped into truth's ALT order before serialization.
fn query_matches_truth_allele_set(query: &Variant, truth: &Variant) -> bool {
    query.key.chrom == truth.key.chrom
        && query.key.pos == truth.key.pos
        && query.key.ref_allele == truth.key.ref_allele
        && query.key.alt_allele.split(',').collect::<BTreeSet<_>>()
            == truth.key.alt_allele.split(',').collect::<BTreeSet<_>>()
}

fn is_distinct_hetalt(gt: &str) -> bool {
    let alleles = parse_gt_alleles(gt);
    alleles.len() == 2 && alleles[0] > 0 && alleles[1] > 0 && alleles[0] != alleles[1]
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
///
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
#[cfg(test)]
fn gvcf2bed_padding(
    truth: &[Variant],
    target_bed: Option<&[vcf::BedInterval]>,
) -> Vec<vcf::BedInterval> {
    gvcf2bed_padding_iter(truth.iter().map(Ok::<_, anyhow::Error>), target_bed)
        .expect("in-memory variants are infallible")
}

fn gvcf2bed_padding_iter<I, V>(
    truth: I,
    target_bed: Option<&[vcf::BedInterval]>,
) -> Result<Vec<vcf::BedInterval>>
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
    Ok(out)
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

fn report_subset_size(
    contig_non_n_lengths: &BTreeMap<String, usize>,
    contigs_in_play: &BTreeSet<String>,
    explicit_bcf: bool,
    implicit_bcf: bool,
) -> usize {
    if explicit_bcf {
        // The legacy BCF path initializes its aggregate region from the
        // complete FASTA dictionary, while the VCF path restricts it to the
        // contigs participating in the comparison. Preserve that observable
        // reporter quirk even though output encoding does not change calls.
        contig_non_n_lengths.values().sum()
    } else {
        let mut selected = contigs_in_play.clone();
        if implicit_bcf {
            // Paired BCF inputs select BCF intermediates and reports, but do
            // not set argparse's explicit `bcf` value. The pinned wrapper's
            // default chromosome discovery consequently retains both aliases
            // when the FASTA declares, for example, `1` and `chr1`.
            for contig in contigs_in_play {
                let alias = contig
                    .strip_prefix("chr")
                    .map(str::to_string)
                    .unwrap_or_else(|| format!("chr{contig}"));
                if contig_non_n_lengths.contains_key(&alias) {
                    selected.insert(alias);
                }
            }
        }
        selected
            .iter()
            .filter_map(|contig| contig_non_n_lengths.get(contig))
            .sum()
    }
}

/// Return the comma-delimited values for one exact INFO key.
///
/// `Regions` is normally the final field, but `--preserve-info` appends the
/// source annotations after it. Parsing the remainder of the INFO string as
/// region names consequently folds the next `;KEY=value` into the last tag.
fn info_list_values<'a>(info: &'a str, key: &str) -> Vec<&'a str> {
    info.split(';')
        .find_map(|entry| entry.split_once('=').filter(|(name, _)| *name == key))
        .map(|(_, value)| value.split(',').filter(|value| !value.is_empty()).collect())
        .unwrap_or_default()
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
        let subset_tags = info_list_values(fields[7], "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
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
fn derive_fp_classes(rows: &[AnnotatedRow], pass_only: bool) -> BTreeMap<String, (usize, usize)> {
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

/// Per-subset variant of `derive_fp_classes`. Each row's `Regions=` tail
/// contributes to every named subset it carries (CONF is filtered out
/// the same way `derive_subset_counts` does so the keys align with the
/// subset rows in extended.csv).
fn derive_subset_fp_classes(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, (usize, usize)>> {
    let mut out: BTreeMap<String, BTreeMap<String, (usize, usize)>> = BTreeMap::new();
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
        let subset_tags = info_list_values(fields[7], "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
        if subset_tags.is_empty() {
            continue;
        }
        let format_keys: Vec<&str> = fields[8].split(':').collect();
        let query_parts: Vec<&str> = fields[10].split(':').collect();
        let query_sample = SampleView::new(&format_keys, &query_parts);
        let Some(variant_type) = query_sample.variant_type() else {
            continue;
        };
        for subset in subset_tags {
            let bucket = out
                .entry(subset.to_string())
                .or_default()
                .entry(variant_type.to_string())
                .or_default();
            if class == "gt" {
                bucket.0 += 1;
            } else if class == "al" {
                bucket.1 += 1;
            }
        }
    }
    out
}

/// Per-(variant_type, subtype) FP class tally for INDEL subtype rows in
/// extended.csv. Legacy emits FP.gt / FP.al at every (INDEL, subtype, *,
/// filter) row and at every (INDEL, subtype, TS_boundary|TS_contained,
/// filter) row. Each FP query row contributes to every indel-class token
/// in its multi-allelic BI (e.g. a hetalt FP with BI `i1_5,i6_15` adds
/// one to both I1_5 and I6_15 — same fanout rule the truth/query stats
/// use).
type SubtypeFpClasses = BTreeMap<String, BTreeMap<String, (usize, usize)>>;
type SubsetSubtypeFpClasses = BTreeMap<String, SubtypeFpClasses>;

fn derive_subtype_fp_classes(rows: &[AnnotatedRow], pass_only: bool) -> SubtypeFpClasses {
    let mut out = SubtypeFpClasses::new();
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
        let Some((variant_type, subtypes)) = query_sample.variant_type_and_subtypes() else {
            continue;
        };
        for subtype in subtypes {
            let bucket = out
                .entry(variant_type.to_string())
                .or_default()
                .entry(subtype)
                .or_default();
            if class == "gt" {
                bucket.0 += 1;
            } else if class == "al" {
                bucket.1 += 1;
            }
        }
    }
    out
}

/// Per-(subset, variant_type, subtype) FP class tally — the cross-product
/// counterpart of `derive_subtype_fp_classes`, plumbed into the
/// (INDEL, subtype, TS_*, filter) extended-csv rows.
fn derive_subset_subtype_fp_classes(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> SubsetSubtypeFpClasses {
    let mut out = SubsetSubtypeFpClasses::new();
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
        let subset_tags = info_list_values(fields[7], "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
        if subset_tags.is_empty() {
            continue;
        }
        let format_keys: Vec<&str> = fields[8].split(':').collect();
        let query_parts: Vec<&str> = fields[10].split(':').collect();
        let query_sample = SampleView::new(&format_keys, &query_parts);
        let Some((variant_type, subtypes)) = query_sample.variant_type_and_subtypes() else {
            continue;
        };
        for subset in subset_tags {
            for subtype in &subtypes {
                let bucket = out
                    .entry(subset.to_string())
                    .or_default()
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default();
                if class == "gt" {
                    bucket.0 += 1;
                } else if class == "al" {
                    bucket.1 += 1;
                }
            }
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
    let mut out: BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>> = BTreeMap::new();
    for row in rows {
        let fields: Vec<&str> = row.line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }
        let subset_tags = info_list_values(fields[7], "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
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
            if let Some((variant_type, subtypes)) = truth_sample.variant_type_and_subtypes() {
                for subtype in subtypes {
                    let stats = by_type
                        .entry(variant_type.to_string())
                        .or_default()
                        .entry(subtype)
                        .or_default();
                    truth_sample.add_truth(stats, filtered_out);
                }
            }
            if !filtered_out
                && let Some((variant_type, subtypes)) = query_sample.variant_type_and_subtypes()
            {
                for subtype in subtypes {
                    let stats = by_type
                        .entry(variant_type.to_string())
                        .or_default()
                        .entry(subtype)
                        .or_default();
                    query_sample.add_query(stats);
                }
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

        if let Some((variant_type, sub_list)) = truth_sample.variant_type_and_subtypes() {
            for subtype in sub_list {
                let stats = subtypes
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype)
                    .or_default();
                truth_sample.add_truth(stats, filtered_out);
            }
        }
        if !filtered_out
            && let Some((variant_type, sub_list)) = query_sample.variant_type_and_subtypes()
        {
            for subtype in sub_list {
                let stats = subtypes
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype)
                    .or_default();
                query_sample.add_query(stats);
            }
        }
    }
    subtypes
}

struct SampleView<'a> {
    gt: Option<&'a str>,
    bd: Option<&'a str>,
    bi: Option<&'a str>,
    bvt: Option<&'a str>,
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
        }
    }

    fn variant_type(&self) -> Option<&'a str> {
        self.bvt.filter(|value| *value != "NOCALL")
    }

    /// Per-row subtypes for INDEL aggregation. Multi-allelic INDELs emit
    /// comma-joined BI strings (e.g. `i1_5,i6_15` for a 1|2 hetalt with one
    /// primitive in each size bucket). Legacy quantify fans these out per
    /// token — each indel-class primitive contributes to its own subtype
    /// bucket. The `ti`/`tv` tokens that decorate complex INDELs (BI like
    /// `c6_15,tv`) describe the SNP-side of a single complex primitive and
    /// must NOT spawn an extra bucket here, so we drop them.
    fn variant_type_and_subtypes(&self) -> Option<(&'a str, Vec<String>)> {
        let variant_type = self.variant_type()?;
        if variant_type != "INDEL" {
            return None;
        }
        let bi = self.bi?;
        let subtypes: Vec<String> = bi
            .split(',')
            .filter(|tok| !matches!(*tok, "ti" | "tv"))
            .map(|tok| tok.to_uppercase())
            .collect();
        if subtypes.is_empty() {
            return None;
        }
        Some((variant_type, subtypes))
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
        // BI on multi-allelic hetalt SNPs is comma-separated (e.g.
        // `ti,tv` for GT=1|2 with one transition and one transversion
        // active alt). Legacy fans these out per primitive — one ti
        // and one tv contribution. Iterate the comma list and count
        // each tag once so single-allelic rows still increment by 1.
        if let Some(bi) = sample.bi {
            for tag in bi.split(',') {
                match tag {
                    "ti" => bucket.ti += 1,
                    "tv" => bucket.tv += 1,
                    _ => {}
                }
            }
        }
    }
    // Legacy's summary counts:
    // * het = exactly one allele is the reference index 0 (covers 0/1,
    //   1/0, 0|1, 1|0 AND 0|2, 2|0, 0|3, 3|0, …) — anything heterozygous
    //   with the reference base.
    // * homalt = both alleles equal AND non-zero (1/1, 1|1, 2|2, 3|3, …).
    // Hetalt (1|2, 2|1, …) lands in NEITHER bucket. Earlier the
    // classifier matched only literal `1/1`/`0/1` to mirror QUERY-side
    // counts (queries are split into per-primitive `0/1`/`1/1` rows),
    // but TRUTH-side rows preserve their original multi-allelic GT
    // through bcftools merge so the literal-only rule under-counted
    // every truth-only multi-allelic record.
    if let Some(gt) = sample.gt {
        let alleles: Vec<usize> = gt
            .split(['/', '|'])
            .map(|part| part.parse::<usize>().unwrap_or(0))
            .collect();
        if alleles.len() == 2 {
            let zero_count = alleles.iter().filter(|a| **a == 0).count();
            if zero_count == 1 {
                bucket.het += 1;
            } else if alleles[0] != 0 && alleles[0] == alleles[1] {
                bucket.homalt += 1;
            }
        }
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

/// Legacy `VariantWriter` assigns the record-level QUAL as the maximum QUAL
/// across every sample call on the merged record. Keep per-sample QQ sourced
/// from the query, but select the combined VCF column independently.
fn combined_record_qual<'a>(truth: &'a Variant, query: &'a Variant) -> &'a str {
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
        xcmp_ctype: None,
        xcmp_hap_match: false,
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_gt}:TP:gm:{info}:{type_label}:{truth_loc}:{qq}\t{query_gt}:TP:gm:{info}:{type_label}:{query_loc}:{qq}",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = display_ref(truth, reference),
            alt = display_alt(truth),
            qual = combined_record_qual(truth, query),
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

/// Combined UNK+UNK row for a same-key truth+query pair where the locus
/// falls outside the confident region. Legacy emits BD=UNK on both
/// samples with BK=lm — the locus-match heuristic fires because the two
/// sides share the variant exactly, and the unconditional non-CONF →
/// UNK rewrite trumps the gm verdict from xcmp. Truth-side QQ stays `.`
/// (truth's input qual is always "0") while query-side carries the
/// query's own qual.
fn unk_combined_row(
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
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_gt}:UNK:lm:{info}:{type_label}:{truth_loc}:.\t{query_gt}:UNK:lm:{info}:{type_label}:{query_loc}:{qq}",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = display_ref(truth, reference),
            alt = display_alt(truth),
            qual = combined_record_qual(truth, query),
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
            line: format!(
                "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{gt}:TP:gm:{info}:{type_label}:{loc}:{qq}\t./.:.:.:.:NOCALL:nocall:0",
                chrom = variant.key.chrom,
                pos = variant.key.pos,
                ref = display_ref(variant, reference),
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
            xcmp_ctype: None,
            xcmp_hap_match: false,
            line: format!(
                "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t./.:.:.:.:NOCALL:nocall:.\t{gt}:TP:gm:{info}:{type_label}:{loc}:{qq}",
                chrom = variant.key.chrom,
                pos = variant.key.pos,
                ref = display_ref(variant, reference),
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
/// Legacy emits BD=FN on truth, BD=FP on query, BK=am on both. Record QUAL
/// is the maximum call QUAL, while query QQ retains the query's own score
/// for downstream ROC enumeration.
#[allow(clippy::too_many_arguments)] // Mirrors the two-sample legacy VCF row contract.
fn fn_fp_combined_row(
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
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_gt}:{truth_bd}:{bk}:{truth_info}:{type_label}:{truth_loc}:.\t{query_gt}:{query_bd}:{bk}:{query_info}:{type_label}:{query_loc}:{qq}",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = display_ref(truth, reference),
            alt = display_alt(truth),
            qual = combined_record_qual(truth, query),
            filter = filter_for_output(&query.filter),
            bs = block_start,
            regions = regions,
            truth_gt = truth.gt,
            query_gt = query.gt,
            truth_info = truth_info,
            query_info = query_info,
            type_label = truth.primary_type(),
            truth_loc = genotype_label(truth),
            query_loc = genotype_label(query),
            qq = query.qual,
            truth_bd = truth_bd,
            query_bd = query_bd,
            bk = bk,
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
        xcmp_ctype: None,
        xcmp_hap_match: false,
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t.\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{gt}:FN:{bk}:{info}:{type_label}:{loc}:.\t./.:.:.:.:NOCALL:nocall:0",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = display_ref(truth, reference),
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
        xcmp_ctype: None,
        xcmp_hap_match: false,
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t.\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t{gt}:UNK:{bk}:{info}:{type_label}:{loc}:.\t./.:.:.:.:NOCALL:nocall:0",
            chrom = truth.key.chrom,
            pos = truth.key.pos,
            ref = display_ref(truth, reference),
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
fn split_query_primitives_with_neighbors(
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
        xcmp_ctype: None,
        xcmp_hap_match: false,
        line: format!(
            "{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter}\tBS={bs}{regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t./.:.:.:.:NOCALL:nocall:.\t{gt}:{bd}:{bk}:{info}:{type_label}:{loc}:{qq}",
            chrom = query.key.chrom,
            pos = query.key.pos,
            ref = display_ref(query, reference),
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

fn display_ref(variant: &Variant, reference: &str) -> String {
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
    if all_alleles.len() == 2 && (all_alleles[0] == 0) != (all_alleles[1] == 0) {
        return "het";
    }
    "nocall"
}

fn comparison_info(variant: &Variant, reference: &str) -> String {
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

/// Resolve the parent under which a unique per-invocation scratch directory is
/// created. `--scratch-prefix` is a caller-selected parent; without it, keep
/// scratch beside the reports as legacy hap.py does.
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
mod scratch_tests {
    use super::*;
    use std::thread;

    #[test]
    fn comparison_headers_merge_inputs_and_append_legacy_annotations() {
        let truth = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=TRUTH_ONLY,Number=1,Type=String,Description=\"truth\">".to_string(),
            "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"wrong\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT".to_string(),
        ];
        let query = vec![
            "##INFO=<ID=QUERY_ONLY,Number=1,Type=String,Description=\"query\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQ".to_string(),
        ];
        let headers = build_vcf_headers(&truth, &query, true, false, false, "QUAL");

        assert!(headers.iter().any(|line| line.contains("ID=TRUTH_ONLY,")));
        assert!(headers.iter().any(|line| line.contains("ID=QUERY_ONLY,")));
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=Q_FILTERED,"))
        );
        assert!(headers.iter().any(|line| line == "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">"));
        assert!(headers.iter().any(|line| line == "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">"));
        assert!(
            !headers
                .iter()
                .any(|line| line.contains("Description=\"wrong\""))
        );
        assert_eq!(
            headers.last().unwrap(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY"
        );
    }

    #[test]
    fn comparison_headers_only_declare_filtered_calls_when_enabled() {
        let headers = build_vcf_headers(&[], &[], false, false, false, "QUAL");
        assert!(!headers.iter().any(|line| line.contains("ID=Q_FILTERED,")));
    }

    fn fixture_path(file: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/synth-snp-match")
            .join(file)
    }

    fn test_root(label: &str) -> PathBuf {
        let id = SCRATCH_RUN_ID.fetch_add(1, Ordering::Relaxed);
        let root =
            std::env::temp_dir().join(format!("hap-compare-{label}-{}-{id}", std::process::id()));
        fs::create_dir_all(&root).expect("create test root");
        root
    }

    fn args(report_prefix: PathBuf, scratch_parent: &Path, keep_scratch: bool) -> CompareArgs {
        let mut args = CompareArgs::with_paths(
            fixture_path("truth.vcf").display().to_string(),
            fixture_path("query.vcf").display().to_string(),
            fixture_path("ref.fa").display().to_string(),
            report_prefix.display().to_string(),
        );
        args.scratch_prefix = Some(scratch_parent.display().to_string());
        args.keep_scratch = keep_scratch;
        args
    }

    #[test]
    fn preprocessing_defaults_are_asymmetric_between_truth_and_query() {
        let root = test_root("preprocessing-defaults");
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.bcftools_norm = true;
        let truth = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth.vcf.gz"),
            true,
            options.preprocess_truth,
        );
        let query = build_preprocess_args(
            &options,
            &options.query,
            &root.join("query.vcf.gz"),
            false,
            true,
        );

        assert!(!truth.leftshift);
        assert!(!truth.decompose);
        assert!(!truth.bcftools_norm);
        assert!(query.leftshift);
        assert!(query.decompose);
        assert!(query.bcftools_norm);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn custom_roc_regions_retain_the_implicit_aggregate_region() {
        let mut regions = vec!["CONF".to_string()];
        ensure_aggregate_roc_region(&mut regions);
        ensure_aggregate_roc_region(&mut regions);
        assert_eq!(regions, ["*", "CONF"]);
    }

    #[test]
    fn scmp_engines_apply_legacy_preprocessing_defaults() {
        let root = test_root("scmp-policy");
        let mut somatic = args(root.join("somatic"), &root, false);
        somatic.engine = CompareEngine::ScmpSomatic;
        somatic.preprocess_truth = true;
        somatic.bcftools_norm = true;
        normalize_engine_preprocessing(&mut somatic);
        assert!(somatic.somatic);
        assert_eq!(somatic.set_gt, None);
        assert!(!somatic.preprocess_truth);
        assert!(somatic.no_leftshift);
        assert!(!somatic.bcftools_norm);
        assert!(!effective_decomposition(&somatic));

        let mut distance = args(root.join("distance"), &root, false);
        distance.engine = CompareEngine::ScmpDistance;
        normalize_engine_preprocessing(&mut distance);
        assert_eq!(distance.set_gt, Some(crate::cli::SomaticGtMode::First));
        assert!(!effective_decomposition(&distance));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn comparison_log_honors_verbose_and_quiet_levels() {
        let root = test_root("comparison-log");
        let logfile = root.join("comparison.log");
        let mut options = args(root.join("result"), &root, false);
        options.logfile = Some(logfile.display().to_string());
        options.verbose = true;
        initialize_compare_log(&options).unwrap();
        log_compare_info(&options, "comparison stage").unwrap();
        assert_eq!(
            fs::read_to_string(&logfile).unwrap(),
            "INFO comparison stage\n"
        );

        options.quiet = true;
        log_compare_info(&options, "suppressed").unwrap();
        assert!(!fs::read_to_string(&logfile).unwrap().contains("suppressed"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn preprocess_truth_and_negative_switches_control_both_sides() {
        let root = test_root("preprocessing-overrides");
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.preprocess_truth = true;

        let truth = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth-enabled.vcf.gz"),
            true,
            options.preprocess_truth,
        );
        assert!(truth.leftshift);
        assert!(truth.decompose);

        options.no_leftshift = true;
        options.no_decompose = true;
        let truth_disabled = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth-disabled.vcf.gz"),
            true,
            options.preprocess_truth,
        );
        let query_disabled = build_preprocess_args(
            &options,
            &options.query,
            &root.join("query-disabled.vcf.gz"),
            false,
            true,
        );
        assert!(!truth_disabled.leftshift);
        assert!(!truth_disabled.decompose);
        assert!(!query_disabled.leftshift);
        assert!(!query_disabled.decompose);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn remaining_germline_preprocess_controls_propagate_per_side() {
        let root = test_root("preprocessing-controls");
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.usefiltered_truth = true;
        options.filters_only = Some("LowQual,q10".to_string());
        options.convert_gvcf_truth = true;
        options.convert_gvcf_query = false;
        options.filter_nonref = true;
        options.preprocess_truth = true;
        options.bcftools_norm = true;
        options.fixchr = Some(true);
        options.gender = crate::cli::PreprocessGender::Male;
        options.preprocess_window = 4096;

        let truth = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth.bcf"),
            !options.usefiltered_truth,
            options.preprocess_truth,
        );
        let query = build_preprocess_args(
            &options,
            &options.query,
            &root.join("query.bcf"),
            options.pass_only,
            true,
        );
        assert!(!truth.pass_only);
        assert_eq!(truth.filters_only, None);
        assert!(truth.convert_gvcf_to_vcf);
        assert!(truth.filter_nonref);
        assert_eq!(query.filters_only.as_deref(), Some("LowQual,q10"));
        assert!(!query.convert_gvcf_to_vcf);
        assert!(query.filter_nonref);
        for side in [&truth, &query] {
            assert!(side.bcftools_norm);
            assert_eq!(side.fixchr, Some(true));
            assert_eq!(side.gender, crate::cli::PreprocessGender::Male);
            assert_eq!(side.window_size, 4096);
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn no_roc_no_counts_no_json_emits_the_legacy_artifact_shape() {
        let root = test_root("minimal-artifacts");
        let prefix = root.join("result");
        let mut options = args(prefix.clone(), &root.join("scratch"), false);
        options.no_roc = true;
        options.no_write_counts = true;
        options.no_json = true;
        run(options).unwrap();

        for suffix in [
            "summary.csv",
            "vcf.gz",
            "vcf.gz.tbi",
            "roc.all.csv.gz",
            "runinfo.json",
        ] {
            assert!(suffixed_report_path(&prefix, suffix).is_file(), "{suffix}");
        }
        for suffix in [
            "extended.csv",
            "metrics.json.gz",
            "roc.Locations.SNP.csv.gz",
            "roc.Locations.SNP.PASS.csv.gz",
            "roc.Locations.INDEL.csv.gz",
            "roc.Locations.INDEL.PASS.csv.gz",
        ] {
            assert!(!suffixed_report_path(&prefix, suffix).exists(), "{suffix}");
        }
        let roc = vcf::read_text(&suffixed_report_path(&prefix, "roc.all.csv.gz")).unwrap();
        for line in roc.lines().skip(1) {
            let cells = line.split(',').collect::<Vec<_>>();
            assert_eq!(cells[6], "*");
            if cells[0] == "INDEL" {
                for block_start in [16usize, 23, 30, 37, 44, 51, 58] {
                    assert_eq!(cells[block_start + 1], ".");
                    assert_eq!(cells[block_start + 2], ".");
                }
            }
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bcf_mode_keeps_bcf_intermediates_and_publishes_only_the_bcf_report_pair() {
        let root = test_root("bcf-artifacts");
        let prefix = root.join("result");
        let scratch_parent = root.join("scratch");
        let mut options = args(prefix.clone(), &scratch_parent, true);
        options.bcf = true;
        run(options).unwrap();

        let bcf_report = suffixed_report_path(&prefix, "bcf");
        assert!(bcf_report.is_file());
        assert!(suffixed_report_path(&prefix, "bcf.csi").is_file());
        assert!(!suffixed_report_path(&prefix, "vcf.gz").exists());
        assert!(!suffixed_report_path(&prefix, "vcf.gz.tbi").exists());
        let (_, records) = vcf::load_raw_vcf(&bcf_report).unwrap();
        assert!(!records.is_empty());

        let scratch_runs = fs::read_dir(&scratch_parent)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect::<Vec<_>>();
        assert_eq!(scratch_runs.len(), 1);
        let scratch = &scratch_runs[0];
        for name in [
            "truth.prep.bcf",
            "truth.prep.bcf.csi",
            "query.prep.bcf",
            "query.prep.bcf.csi",
        ] {
            assert!(scratch.join(name).is_file(), "{name}");
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn paired_bcf_inputs_implicitly_enable_bcf_reports_and_intermediates() {
        let root = test_root("implicit-bcf-artifacts");
        let truth_bcf = root.join("truth.bcf");
        let query_bcf = root.join("query.bcf");
        for (source, destination) in [
            (fixture_path("truth.vcf"), &truth_bcf),
            (fixture_path("query.vcf"), &query_bcf),
        ] {
            let (headers, records) = vcf::load_raw_vcf(&source).unwrap();
            vcf::write_raw_vcf(destination, &headers, &records).unwrap();
        }

        let prefix = root.join("result");
        let scratch_parent = root.join("scratch");
        let mut options = args(prefix.clone(), &scratch_parent, true);
        options.truth = truth_bcf.display().to_string();
        options.query = query_bcf.display().to_string();
        assert!(!options.bcf, "the CLI flag is intentionally absent");
        run(options).unwrap();

        assert!(suffixed_report_path(&prefix, "bcf").is_file());
        assert!(suffixed_report_path(&prefix, "bcf.csi").is_file());
        assert!(!suffixed_report_path(&prefix, "vcf.gz").exists());
        let runinfo = fs::read_to_string(suffixed_report_path(&prefix, "runinfo.json")).unwrap();
        assert!(
            runinfo.contains("\"bcf\":false"),
            "implicit output selection must not rewrite the explicit CLI flag: {runinfo}"
        );
        let scratch_runs = child_directories(&scratch_parent);
        assert_eq!(scratch_runs.len(), 1);
        assert!(scratch_runs[0].join("truth.prep.bcf").is_file());
        assert!(scratch_runs[0].join("query.prep.bcf").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn missing_report_parent_is_rejected_without_creating_artifacts() {
        let root = test_root("missing-report-parent");
        let missing_parent = root.join("missing");
        let scratch_parent = root.join("scratch");
        let options = args(missing_parent.join("result"), &scratch_parent, false);

        let error = run(options).expect_err("missing report parents must be rejected");
        assert!(error.to_string().contains("output path does not exist"));
        assert!(!missing_parent.exists());
        assert!(!scratch_parent.exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn regions_reject_legacy_overlaps_and_ordering_but_targets_do_not() {
        let root = test_root("region-overlap-check");
        for (label, contents) in [
            ("overlap", "chr1\t0\t10\nchr1\t5\t12\n"),
            ("out-of-order", "chr1\t10\t12\nchr1\t0\t6\n"),
        ] {
            let bed = root.join(format!("{label}.bed"));
            fs::write(&bed, contents).unwrap();
            let mut options = args(
                root.join(format!("{label}-result")),
                &root.join(format!("{label}-scratch")),
                false,
            );
            options.regions_bedfile = Some(bed.display().to_string());
            let error = run(options).expect_err("invalid -R BED must fail before comparison");
            assert!(
                error
                    .to_string()
                    .contains("The regions bed file (specified using -R) has overlaps")
            );
        }

        let targets = root.join("targets.bed");
        fs::write(&targets, "chr1\t10\t12\nchr1\t0\t12\n").unwrap();
        let target_prefix = root.join("target-result");
        let mut options = args(target_prefix.clone(), &root.join("target-scratch"), false);
        options.targets_bedfile = Some(targets.display().to_string());
        run(options).expect("-T keeps accepting overlapping or out-of-order intervals");
        assert!(suffixed_report_path(&target_prefix, "summary.csv").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn empty_truth_with_default_locations_fails_before_query_comparison() {
        let root = test_root("empty-truth-default-locations");
        let truth = root.join("truth.vcf");
        fs::write(
            &truth,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=16>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n",
            ),
        )
        .unwrap();
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.truth = truth.display().to_string();

        let error = run(options).expect_err("legacy derives default contigs from truth calls");
        assert!(
            error
                .to_string()
                .contains("Truth and reference have no chromosomes in common")
        );

        let explicit_prefix = root.join("explicit-result");
        let mut explicit = args(
            explicit_prefix.clone(),
            &root.join("explicit-scratch"),
            false,
        );
        explicit.truth = truth.display().to_string();
        explicit.locations = Some("chr1".to_string());
        run(explicit).expect("an explicit contig bypasses legacy default-contig discovery");
        assert!(suffixed_report_path(&explicit_prefix, "summary.csv").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bcf_report_size_spans_the_full_reference() {
        let contig_lengths = BTreeMap::from([("chr1".to_string(), 100), ("chrX".to_string(), 40)]);
        let contigs_in_play = BTreeSet::from(["chr1".to_string()]);

        assert_eq!(
            report_subset_size(&contig_lengths, &contigs_in_play, false, false),
            100,
            "ordinary reports retain their active-contig size"
        );
        assert_eq!(
            report_subset_size(&contig_lengths, &contigs_in_play, true, false),
            140,
            "legacy BCF reports size the complete reference dictionary"
        );

        let aliased_lengths = BTreeMap::from([
            ("1".to_string(), 100),
            ("chr1".to_string(), 100),
            ("chrX".to_string(), 40),
        ]);
        assert_eq!(
            report_subset_size(&aliased_lengths, &contigs_in_play, false, true),
            200,
            "implicit BCF reports retain both declared chromosome aliases"
        );
    }

    #[test]
    fn scmp_bcf_mode_materializes_confidence_padding_without_exposing_vcf() {
        let root = test_root("scmp-bcf-artifacts");
        let prefix = root.join("result");
        let confidence = root.join("confident.bed");
        fs::write(&confidence, "chr1\t0\t16\n").unwrap();
        let mut options = args(prefix.clone(), &root.join("scratch"), false);
        options.engine = CompareEngine::ScmpDistance;
        options.bcf = true;
        options.fp_bedfile = Some(confidence.display().to_string());
        run(options).unwrap();

        assert!(suffixed_report_path(&prefix, "bcf").is_file());
        assert!(suffixed_report_path(&prefix, "bcf.csi").is_file());
        assert!(!suffixed_report_path(&prefix, "vcf.gz").exists());
        assert!(!suffixed_report_path(&prefix, "vcf.gz.tbi").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn scmp_output_vtc_is_forwarded_to_ga4gh_quantification() {
        let root = test_root("scmp-output-vtc");
        let prefix = root.join("result");
        let mut options = args(prefix.clone(), &root.join("scratch"), false);
        options.engine = CompareEngine::ScmpDistance;
        options.output_vtc = true;
        run(options).unwrap();

        let (headers, records) =
            vcf::load_raw_vcf(&suffixed_report_path(&prefix, "vcf.gz")).unwrap();
        assert!(headers.iter().any(|line| line.contains("##INFO=<ID=VTC,")));
        assert!(records.iter().any(|record| record.info.contains("VTC=")));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn subset_derivation_stops_regions_at_the_next_info_field() {
        let rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            line: concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\t",
                "BS=7;Regions=CONF,TS_boundary,TS_contained;AF=0.5;VTC=nuc__s\t",
                "GT:BD:BK:BVT:BLT:BI\t",
                "0/1:TP:gm:SNP:het:ti\t0/1:TP:gm:SNP:het:ti"
            )
            .to_string(),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];

        let counts = derive_subset_counts(&rows, false);
        assert_eq!(
            counts.keys().cloned().collect::<Vec<_>>(),
            vec!["TS_boundary", "TS_contained"]
        );
        for subset in ["TS_boundary", "TS_contained"] {
            let snp = &counts[subset]["SNP"];
            assert_eq!(snp.truth_total.total, 1);
            assert_eq!(snp.query_total.total, 1);
        }
    }

    #[test]
    fn requantify_handoff_drops_only_provisional_truth_set_membership() {
        let rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            line: concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\t",
                "BS=7;Regions=CONF,TS_boundary,EXTRA,TS_contained;RegionsExtent=7-7\t",
                "GT:BD\t0/1:TP\t0/1:TP"
            )
            .to_string(),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];

        let sanitized = sanitize_requantify_handoff_rows(&rows);

        assert!(
            rows[0]
                .line
                .contains("Regions=CONF,TS_boundary,EXTRA,TS_contained")
        );
        assert!(
            sanitized[0]
                .line
                .contains("BS=7;Regions=CONF,EXTRA;RegionsExtent=7-7")
        );
        assert!(!sanitized[0].line.contains("TS_boundary"));
        assert!(!sanitized[0].line.contains("TS_contained"));
    }

    #[test]
    fn vcfeval_ignores_deprecated_external_runtime_flags() {
        let root = test_root("vcfeval-deprecated-flags");
        let mut options = CompareArgs::with_paths(
            fixture_path("truth.vcf").display().to_string(),
            fixture_path("query.vcf").display().to_string(),
            fixture_path("ref.fa").display().to_string(),
            root.join("result").display().to_string(),
        );
        options.scratch_prefix = Some(root.join("scratch").display().to_string());
        options.engine = CompareEngine::Vcfeval;
        options.engine_vcfeval = Some("definitely-absent-rtg-for-test".to_string());
        options.engine_vcfeval_template = Some(root.join("absent.sdf").display().to_string());
        run(options).unwrap();
        assert!(root.join("result.summary.csv").is_file());
        let (_, records) = vcf::load_raw_vcf(&root.join("result.vcf.gz")).unwrap();
        assert_eq!(records.len(), 1);
        assert!(
            records[0]
                .samples
                .iter()
                .all(|sample| sample.contains(":TP:gm"))
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn native_vcfeval_with_paired_bcf_inputs_publishes_bcf_and_preserves_info() {
        let root = test_root("vcfeval-contract");
        let truth_bcf = root.join("truth.bcf");
        let query_bcf = root.join("query.bcf");
        for (source, destination) in [
            (fixture_path("truth.vcf"), &truth_bcf),
            (fixture_path("query.vcf"), &query_bcf),
        ] {
            let (headers, records) = vcf::load_raw_vcf(&source).unwrap();
            vcf::write_raw_vcf(destination, &headers, &records).unwrap();
        }

        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.truth = truth_bcf.display().to_string();
        options.query = query_bcf.display().to_string();
        options.engine = CompareEngine::Vcfeval;
        options.output_vtc = true;
        run(options).unwrap();

        assert!(!root.join("result.vcf.gz").exists());
        let (headers, records) = vcf::load_raw_vcf(&root.join("result.bcf")).unwrap();
        assert!(headers.iter().any(|line| line.contains("ID=VTC,")));
        assert!(headers.iter().any(|line| line.contains("ID=XCMP,")));
        assert!(
            records[0].info.contains("VTC=nuc__s,al__s,homalt__s"),
            "{}",
            records[0].info
        );
        assert!(records[0].info.contains("XCMP="), "{}", records[0].info);

        let preserve_prefix = root.join("preserve-result");
        let preserve_scratch = root.join("preserve-scratch");
        let mut preserve = args(preserve_prefix.clone(), &preserve_scratch, false);
        preserve.engine = CompareEngine::Vcfeval;
        preserve.preserve_info = true;
        run(preserve).unwrap();
        assert!(suffixed_report_path(&preserve_prefix, "runinfo.json").is_file());
        for suffix in ["summary.csv", "extended.csv", "vcf.gz", "metrics.json.gz"] {
            assert!(suffixed_report_path(&preserve_prefix, suffix).is_file());
        }
        let (_, preserved_records) =
            vcf::load_raw_vcf(&suffixed_report_path(&preserve_prefix, "vcf.gz")).unwrap();
        assert_eq!(preserved_records.len(), 1);
        assert!(child_directories(&preserve_scratch).is_empty());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn scmp_engines_use_non_haplotype_comparison_semantics() {
        let root = test_root("scmp-engines");
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/synth-homopolymer-insertion");
        let engine_args = |prefix: PathBuf| {
            let mut options = CompareArgs::with_paths(
                fixture.join("truth.vcf").display().to_string(),
                fixture.join("query.vcf").display().to_string(),
                fixture.join("ref.fa").display().to_string(),
                prefix.display().to_string(),
            );
            options.scratch_prefix = Some(root.join("scratch").display().to_string());
            options
        };
        let mut somatic = engine_args(root.join("somatic"));
        somatic.engine = CompareEngine::ScmpSomatic;
        run(somatic).unwrap();
        let somatic_summary = fs::read_to_string(root.join("somatic.summary.csv")).unwrap();

        let mut distance = engine_args(root.join("distance"));
        distance.engine = CompareEngine::ScmpDistance;
        distance.engine_scmp_distance = 30;
        run(distance).unwrap();
        let distance_summary = fs::read_to_string(root.join("distance.summary.csv")).unwrap();
        // Legacy AlleleMatcher's RefVar constructor uses ALT length for the
        // reference end. Consequently these two ordinary VCF-equivalent
        // homopolymer insertions do not hash alike in scmp-somatic, while
        // distance mode still pairs their overlapping intervals.
        assert_ne!(somatic_summary, distance_summary);

        run(engine_args(root.join("xcmp"))).unwrap();
        let xcmp_summary = fs::read_to_string(root.join("xcmp.summary.csv")).unwrap();
        assert_ne!(somatic_summary, xcmp_summary);
        assert!(distance_summary.contains("INDEL,ALL,1,1,0,1,0,0"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn preserve_info_and_vtc_decorate_rows_and_headers() {
        let source = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tC\t30\tPASS\tSCORE=9\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        let mut rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            line: concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\tBS=7;Regions=CONF\t",
                "GT:BD:BK:BVT:BLT:QQ\t",
                "0/1:TP:gm:SNP:het:30\t0/1:TP:gm:SNP:het:30"
            )
            .to_string(),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];
        decorate_output_rows(&mut rows, &[source], &[], true, true, "QUAL").unwrap();
        assert!(
            rows[0].line.contains(concat!(
                "BS=7;IQQ=30;SCORE=9;ctype=simple:match;gtt1=gt_het;",
                "gtt2=gt_het;kind=match;type=TP;Regions=CONF;RegionsExtent=7-7;",
                "XCMP=TP:match:gt_het:gt_het:simple:match;",
                "VTC=nuc__s,al__s,het__rs"
            )),
            "{}",
            rows[0].line
        );
        let headers = build_vcf_headers(&[], &[], false, true, true, "QUAL");
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=VTC,"))
        );
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=XCMP,"))
        );
        let regions = headers
            .iter()
            .position(|line| line.starts_with("##INFO=<ID=Regions,"))
            .unwrap();
        let extent = headers
            .iter()
            .position(|line| line.starts_with("##INFO=<ID=RegionsExtent,"))
            .unwrap();
        let vtc = headers
            .iter()
            .position(|line| line.starts_with("##INFO=<ID=VTC,"))
            .unwrap();
        assert_eq!(extent, regions + 1);
        assert!(vtc > extent);
    }

    #[test]
    fn prefixed_custom_roc_field_preserves_xcmp_literal_lookup_miss() {
        let source = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tC\t30\tPASS\tSCORE=9\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        let mut rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            line: concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\tBS=7;Regions=CONF\t",
                "GT:BD:BK:BVT:BLT:QQ\t",
                "0/1:TP:gm:SNP:het:30\t0/1:TP:gm:SNP:het:30"
            )
            .to_string(),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];

        decorate_output_rows(&mut rows, &[source], &[], true, false, "INFO.SCORE").unwrap();

        let record = RawVcfRecord::from_line(&rows[0].line, Path::new("output.vcf")).unwrap();
        assert!(record.info.contains("SCORE=9"));
        assert!(!record.info.contains("INFO.SCORE="));
        assert!(!record.info.contains("IQQ="));
        assert_eq!(
            record.sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("nan")
        );

        let headers = build_vcf_headers(&[], &[], false, false, true, "INFO.SCORE");
        assert!(headers.iter().any(|line| {
            line == "##INFO=<ID=IQQ,Number=1,Type=Float,Description=\"Quality value for query variants (INFO.SCORE).\">"
        }));
    }

    #[test]
    fn metadata_helpers_preserve_legacy_multiallelic_contracts() {
        let record = RawVcfRecord::from_line(
            concat!(
                "chr21\t19323424\t.\tCGTGT\tC,CGTGTGTGTGT,CGT\t0\t.\t.\t",
                "GT:BVT:BLT\t1/2:INDEL:hetalt\t./.:NOCALL:nocall"
            ),
            Path::new("metadata.vcf"),
        )
        .unwrap();
        assert_eq!(legacy_regions_extent(&record), "19323424-19323428");

        let mixed = RawVcfRecord::from_line(
            concat!(
                "chr21\t10\t.\tAGT\tA,AGTGT\t0\t.\t.\t",
                "GT:BVT:BLT\t1/2:INDEL:hetalt\t./.:NOCALL:nocall"
            ),
            Path::new("metadata.vcf"),
        )
        .unwrap();
        assert_eq!(
            legacy_vtc(&mixed, &mixed.sample_map(0), &mixed.sample_map(1)),
            "nuc__i,nuc__d,al__i,al__d,hetalt__id,nocall__nc"
        );
    }

    #[test]
    fn semantic_preserve_key_reorders_alts_without_colliding_deletions() {
        let first = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tAT,ATT\t0\t.\t.\tGT\t1/2",
            Path::new("source.vcf"),
        )
        .unwrap();
        let reordered = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tATT,AT\t0\t.\t.\tGT\t1/2",
            Path::new("source.vcf"),
        )
        .unwrap();
        let other_ref = RawVcfRecord::from_line(
            "chr1\t7\t.\tAA\tAT,ATT\t0\t.\t.\tGT\t1/2",
            Path::new("source.vcf"),
        )
        .unwrap();
        assert_eq!(semantic_info_key(&first), semantic_info_key(&reordered));
        assert_ne!(semantic_info_key(&first), semantic_info_key(&other_ref));

        let symbolic_n = RawVcfRecord::from_line(
            "chr1\t7\t.\tN\t<DEL>\t0\t.\t.\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        let symbolic_a = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\t<DEL>\t0\t.\t.\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        assert_eq!(
            semantic_info_key(&symbolic_n),
            semantic_info_key(&symbolic_a)
        );
    }

    fn child_directories(parent: &Path) -> Vec<PathBuf> {
        fs::read_dir(parent)
            .expect("read scratch parent")
            .filter_map(|entry| entry.ok())
            .map(|entry| entry.path())
            .filter(|path| path.is_dir())
            .collect()
    }

    #[test]
    fn concurrent_runs_share_parent_without_colliding_and_cleanup() {
        let root = test_root("concurrent");
        let scratch_parent = root.join("scratch");
        fs::create_dir_all(root.join("first")).unwrap();
        fs::create_dir_all(root.join("second")).unwrap();
        let first = args(root.join("first/result"), &scratch_parent, false);
        let second = args(root.join("second/result"), &scratch_parent, false);

        let first_run = thread::spawn(move || run(first));
        let second_run = thread::spawn(move || run(second));
        first_run.join().expect("first thread panicked").unwrap();
        second_run.join().expect("second thread panicked").unwrap();

        assert_eq!(
            fs::read(root.join("first/result.summary.csv")).unwrap(),
            fs::read(root.join("second/result.summary.csv")).unwrap()
        );
        assert!(
            child_directories(&scratch_parent).is_empty(),
            "successful invocations must delete only their own run directories"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn keep_scratch_retains_run_but_errors_cleanup_by_default() {
        let root = test_root("lifecycle");
        let kept_parent = root.join("kept");
        fs::create_dir(root.join("kept-output")).unwrap();
        run(args(root.join("kept-output/result"), &kept_parent, true)).unwrap();

        let kept = child_directories(&kept_parent);
        assert_eq!(kept.len(), 1, "--keep-scratch retains the unique run");
        assert!(kept[0].join("truth.prep.vcf.gz").is_file());
        assert!(kept[0].join("truth.prep.vcf.gz.tbi").is_file());
        assert!(kept[0].join("query.prep.vcf.gz").is_file());
        assert!(kept[0].join("query.prep.vcf.gz.tbi").is_file());

        let error_parent = root.join("error");
        fs::create_dir(root.join("error-output")).unwrap();
        let mut failing = args(root.join("error-output/result"), &error_parent, false);
        failing.truth = root.join("missing.vcf").display().to_string();
        assert!(run(failing).is_err());
        assert!(
            child_directories(&error_parent).is_empty(),
            "error paths must delete their invocation directory"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn report_suffixes_preserve_dotted_prefixes() {
        let prefix = Path::new("reports/sample.v1");
        assert_eq!(
            suffixed_report_path(prefix, "summary.csv"),
            Path::new("reports/sample.v1.summary.csv")
        );
    }

    #[test]
    fn stratification_tsv_requantifies_reports_and_vcf() {
        let root = test_root("stratification");
        let bed = root.join("focus.bed");
        let tsv = root.join("regions.tsv");
        fs::write(&bed, "chr1\t4\t5\n").unwrap();
        fs::write(&tsv, "FOCUS\tfocus.bed\n").unwrap();
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.strat_tsv = Some(tsv.display().to_string());
        run(options).unwrap();

        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        assert!(extended.lines().any(|line| line.contains(",FOCUS,")));
        let (_, records) = vcf::load_raw_vcf(&root.join("result.vcf.gz")).unwrap();
        assert!(
            records
                .iter()
                .any(|record| record.info.contains("Regions=FOCUS"))
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn explicit_scratch_cleanup_propagates_removal_errors() {
        let root = test_root("cleanup-error");
        let scratch = ScratchRun::create(&root, false).unwrap();
        let scratch_path = scratch.path().to_path_buf();
        fs::remove_dir_all(&scratch_path).unwrap();
        fs::write(&scratch_path, "not a directory").unwrap();

        let error = scratch.cleanup().unwrap_err();
        assert!(error.to_string().contains("failed to remove scratch run"));

        fs::remove_file(scratch_path).unwrap();
        fs::remove_dir_all(root).unwrap();
    }
}

#[cfg(test)]
mod memory_guards {
    use super::*;

    // Class 1 pinning helper for ergonomic qual overrides in tests.
    impl Variant {
        fn with_qual(mut self, q: &str) -> Self {
            self.qual = q.to_string();
            self
        }
    }

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
    fn custom_enumeration_threshold_above_default_is_honored() {
        let variants = (0..15).map(|_| het("0/1")).collect::<Vec<_>>();
        assert_eq!(estimated_state_count(&variants), usize::MAX);
        assert_eq!(estimated_state_count_with_limit(&variants, 32_768), 32_768);
    }

    #[test]
    fn build_clusters_splits_at_variant_cap() {
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

    #[test]
    fn streaming_clusters_match_the_in_memory_cluster_boundaries() -> Result<()> {
        let make = |chrom: &str, pos, reference: &str, alternate: &str| Variant {
            key: VariantKey {
                chrom: chrom.to_string(),
                pos,
                ref_allele: reference.to_string(),
                alt_allele: alternate.to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: "0/1".to_string(),
        };
        let truth = vec![
            make("chr1", 10, "A", "C"),
            make("chr1", 25, "A", "G"),
            make("chr2", 3, "T", "TA"),
        ];
        let query = vec![make("chr1", 11, "A", "T"), make("chr1", 40, "C", "G")];
        let expected = build_clusters_with_gap(&truth, &query, 2);
        let actual = StreamingClusters::new(
            truth.into_iter().map(Ok::<_, anyhow::Error>),
            query.into_iter().map(Ok::<_, anyhow::Error>),
            2,
        )
        .collect::<Result<Vec<_>>>()?;
        assert_eq!(actual.len(), expected.len());
        for (actual, expected) in actual.iter().zip(expected) {
            assert_eq!(actual.chrom, expected.chrom);
            assert_eq!((actual.start, actual.end), (expected.start, expected.end));
            assert_eq!(actual.truth.len(), expected.truth.len());
            assert_eq!(actual.query.len(), expected.query.len());
        }
        Ok(())
    }

    /// Class D pin: chr21:38861935 INDEL hetalt-vs-homalt-of-shared-allele.
    /// Truth `T→TAA,TA 1|1` selects {TAA}, query `T→TAA,TA 1/2` selects
    /// {TAA, TA}. Selected sets differ but overlap on TAA. INDEL → BK=lm.
    #[test]
    fn class_d_indel_overlap_emits_lm() {
        let truth = variant(38861935, "T", "TAA,TA", "1|1");
        let query = variant(38861935, "T", "TAA,TA", "1/2");
        assert_eq!(compute_paired_bk(&truth, &query), "lm");
    }

    /// Class D pin: chr21:9922359 SNP hetalt-vs-het-overlapping. Truth
    /// `T→A,C 1|0` selects {A}, query `T→A,C 1/2` selects {A, C}. SNP →
    /// BK=`.` (legacy quirk).
    #[test]
    fn class_d_snp_overlap_emits_dot() {
        let truth = variant(9922359, "T", "A,C", "1|0");
        let query = variant(9922359, "T", "A,C", "1/2");
        assert_eq!(compute_paired_bk(&truth, &query), ".");
    }

    /// Class D pin: same selected set, different multiset (truth het,
    /// query homalt of same allele) → BK=`am`.
    #[test]
    fn class_d_same_set_diff_multiset_emits_am() {
        let truth = variant(100, "T", "A", "0|1");
        let query = variant(100, "T", "A", "1/1");
        assert_eq!(compute_paired_bk(&truth, &query), "am");
    }

    /// Class C pin (post-#79): chr21:30374431-435 cluster signatures.
    /// Truth has multi-position multi-allelic (insert + multi-allelic
    /// G→GT,T) and query has overlapping insert + subst single-allelic
    /// records. The narrow relaxation (truth's alts cover both query
    /// alleles at the conflict pos) lets the query enumeration produce
    /// the truth-matching haplotype pair.
    #[test]
    fn class_c_cluster_signatures_match_after_relaxation() {
        let mut reference = vec![b'N'; 30374450];
        let window = b"ggccTAATTTGTTTTTTTTTT";
        for (i, b) in window.iter().enumerate() {
            reference[30374425 - 1 + i] = *b;
        }
        let reference = String::from_utf8(reference).unwrap();
        let truth = vec![
            variant(30374431, "A", "AT", "1|0"),
            variant(30374435, "G", "GT,T", "2|1"),
        ];
        let query = vec![
            variant(30374435, "G", "GT", "1/1"),
            variant(30374435, "G", "T", "0/1"),
        ];
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 30374431,
            end: 30374435,
            truth: truth.clone(),
            query: query.clone(),
        };
        let relax = compute_class_c_relaxation_positions(&query, &truth);
        assert_eq!(relax, BTreeSet::from([30374435usize]));
        let truth_sig =
            cluster_signature(&cluster, &truth, &reference, None, &BTreeSet::new()).unwrap();
        let query_sig = cluster_signature(&cluster, &query, &reference, None, &relax).unwrap();
        assert!(truth_sig.is_some() && query_sig.is_some());
        assert!(
            truth_sig
                .as_ref()
                .unwrap()
                .intersection(query_sig.as_ref().unwrap())
                .next()
                .is_some()
        );
    }

    /// Class C negative pin: chr21:16328989 — truth `G→GA 1|1` (homalt
    /// insert) vs query `G→GA 1/1 + G→A 0/1`. Truth's alts {GA} do not
    /// cover the SNP `A` so the relaxation must NOT fire; legacy keeps
    /// the strict drain semantics → BK=`.` on the FP query record.
    #[test]
    fn class_c_relaxation_skips_insert_only_truth_counterpart() {
        let truth = vec![variant(16328989, "G", "GA", "1|1")];
        let query = vec![
            variant(16328989, "G", "GA", "1/1"),
            variant(16328989, "G", "A", "0/1"),
        ];
        let relax = compute_class_c_relaxation_positions(&query, &truth);
        assert!(
            relax.is_empty(),
            "truth must include both Insert and Subst alleles for relaxation"
        );
    }

    /// Debug helper retained (and pinned via the test above).
    /// Reference at chr21:30374425-30374445 = "ggccTAATTTGTTTTTTTTTT".
    #[test]
    #[ignore]
    fn debug_class_c_cluster_signatures() {
        // Build a synthetic reference exposing the chr21:30374431-435 window.
        // We need positions 30374431-435 to be "ATTTG" and 30374436-445 = "TTTTTTTTTT".
        let mut reference = vec![b'N'; 30374450];
        let window = b"ggccTAATTTGTTTTTTTTTT"; // 30374425-30374445
        for (i, b) in window.iter().enumerate() {
            reference[30374425 - 1 + i] = *b;
        }
        let reference = String::from_utf8(reference).unwrap();
        let truth = vec![
            variant(30374431, "A", "AT", "1|0"),
            variant(30374435, "G", "GT,T", "2|1"),
        ];
        let query = vec![
            variant(30374435, "G", "GT", "1/1"),
            variant(30374435, "G", "T", "0/1"),
        ];
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 30374431,
            end: 30374435,
            truth: truth.clone(),
            query: query.clone(),
        };
        let relax = compute_class_c_relaxation_positions(&query, &truth);
        let truth_sig =
            cluster_signature(&cluster, &truth, &reference, None, &BTreeSet::new()).unwrap();
        let query_sig = cluster_signature(&cluster, &query, &reference, None, &relax).unwrap();
        eprintln!("relax_positions = {:?}", relax);
        eprintln!("truth_sig = {:?}", truth_sig);
        eprintln!("query_sig = {:?}", query_sig);
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
    fn combined_tp_uses_max_call_qual_while_qq_stays_query_sourced() {
        let truth = variant(25, "A", "G", "1/1").with_qual("60");
        let query = variant(25, "A", "G", "1/1").with_qual("55");
        let row = tp_combined_row(&truth, &query, "A", 25, "");
        let fields = row.line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[5], "60");
        assert!(fields[9].ends_with(":55"));
        assert!(fields[10].ends_with(":55"));

        let higher_query = query.with_qual("65");
        let row = tp_combined_row(&truth, &higher_query, "A", 25, "");
        assert_eq!(row.line.split('\t').nth(5), Some("65"));
    }

    #[test]
    fn exact_only_unphased_indel_block_reaches_legacy_hap_match_verdict() {
        let identical = Cluster {
            chrom: "chr21".to_string(),
            start: 15006495,
            end: 15006495,
            truth: vec![variant(15006495, "A", "ATCTC", "0/1")],
            query: vec![variant(15006495, "A", "ATCTC", "0/1")],
        };
        assert_eq!(
            identical_gt_exact_indel_keys(&identical),
            BTreeSet::from([identical.truth[0].key.clone()])
        );

        let mut phased_truth = identical.clone();
        phased_truth.truth[0].gt = "1|0".to_string();
        assert!(
            identical_gt_exact_indel_keys(&phased_truth).is_empty(),
            "the standard phased-truth fixture must retain its legacy lm verdict"
        );
    }

    #[test]
    fn repetitive_indel_block_promotes_legacy_hap_matched_gt_mismatch() {
        let earlier_truth = variant(15181523, "A", "AT", "0/1");
        let paired_truth = variant(15181526, "A", "AT", "0/1");
        let shared_truth_snp = variant(15181526, "A", "T", "0/1");
        let paired_query = variant(15181526, "A", "AT", "1/1");
        let shared_query_snp = variant(15181526, "A", "T", "0/1");
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 15181523,
            end: 15181526,
            truth: vec![
                earlier_truth.clone(),
                paired_truth.clone(),
                shared_truth_snp.clone(),
            ],
            query: vec![paired_query.clone(), shared_query_snp.clone()],
        };
        let region_state = RegionState {
            conf_enabled: true,
            any_conf: true,
            any_nonconf: true,
            covered_truth: BTreeSet::from([paired_truth.key.clone(), shared_truth_snp.key.clone()]),
            covered_query: BTreeSet::from([paired_query.key.clone(), shared_query_snp.key.clone()]),
            ..RegionState::default()
        };

        assert_eq!(
            legacy_repetitive_indel_hap_promotions(&cluster, &region_state),
            BTreeSet::from([paired_truth.key])
        );
        assert_eq!(
            legacy_preprocessed_snp_first_positions(&cluster),
            BTreeSet::from([15181526])
        );

        let all_conf = RegionState {
            covered_truth: cluster.truth.iter().map(|v| v.key.clone()).collect(),
            ..region_state.clone()
        };
        assert!(
            legacy_repetitive_indel_hap_promotions(&cluster, &all_conf).is_empty(),
            "promotion requires the balancing truth copy outside CONF"
        );

        let mut phased = cluster.clone();
        phased.truth[1].gt = "0|1".to_string();
        assert!(legacy_preprocessed_snp_first_positions(&phased).is_empty());
    }

    #[test]
    fn xcmp_excludes_filtered_truth_after_preprocessing() {
        let pass = variant(100, "A", "G", "1|1");
        let mut filtered = variant(101, "C", "T", "1|1");
        filtered.filter = "OverlapConflict".to_string();
        let mut variants = vec![pass.clone(), filtered];

        retain_xcmp_truth_calls(&mut variants);

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].key, pass.key);
    }

    #[test]
    fn filtered_truth_counterpart_sorts_first_at_shared_locus() {
        let row = |alt: &str, side_rank| AnnotatedRow {
            sort_key: ("chr21".to_string(), 15576177, side_rank, 0),
            line: format!("chr21\t15576177\t.\tG\t{alt}\t0\t.\tBS=15576177\tGT\t./.\t0/1"),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        };
        let mut rows = vec![row("GAAAGAA", 0), row("A", 2)];
        let keys = BTreeSet::from([VariantKey {
            chrom: "chr21".to_string(),
            pos: 15576177,
            ref_allele: "G".to_string(),
            alt_allele: "A".to_string(),
        }]);

        sort_comparison_rows(&mut rows, &keys);

        assert!(rows[0].line.contains("\tG\tA\t"));
    }

    #[test]
    fn adjacent_decomposed_indels_are_split_siblings() {
        let deletion = variant(18757292, "AT", "A", "0/1").with_qual("1110.88");
        let insertion = variant(18757293, "T", "TT", "0/1").with_qual("1110.88");
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: deletion.key.pos,
            end: insertion.end_pos(),
            truth: Vec::new(),
            query: vec![deletion.clone(), insertion.clone()],
        };
        let regions = RegionState::from_cluster(&cluster, "ATTT", Some(&[]));

        assert!(has_nonconf_split_sibling(&deletion, &cluster, &regions));
        assert!(has_nonconf_split_sibling(&insertion, &cluster, &regions));
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
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&query_counterpart), false),
            "."
        );
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
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&query_counterpart), false),
            "lm"
        );
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
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&query_indel), true),
            "lm"
        );
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
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&neighbour_snp), false),
            "."
        );
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
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&filtered_query), false),
            "."
        );
    }

    // Residual #49 — chr21:17566241 exact-match pair with reordered
    // multi-allelic ALT columns. Truth `C→CA,CAA` (phased 1|2) and query
    // `C→CAA,CA` (unphased 1/2) encode the same diploid allele set.
    // Legacy's simpleCompare matches them via the VariantReader's shared
    // allele-unification table and emits one combined TP:gm row; rust
    // used to fall through to `cluster_signature` and split the pair
    // into truth-only FN + query-only FP. The `simple_compare_pairs_match`
    // predicate plus `canonical_hetalt_gt` close this gap.
    #[test]
    fn simple_compare_matches_reordered_multiallelic_hetalt_indel() {
        let truth = variant(17566241, "C", "CA,CAA", "1|2");
        let query = variant(17566241, "C", "CAA,CA", "1/2");
        let reference = "N".repeat(17566250);
        assert!(simple_compare_pairs_match(
            &truth,
            &query,
            &reference,
            17566241,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
        // Query GT must be remapped into truth's ALT index space —
        // `1` (CAA) → truth idx 2, `2` (CA) → truth idx 1, so `1/2` → `2/1`.
        assert_eq!(canonical_hetalt_gt(&truth.key.alt_allele, &query), "2/1");
    }

    #[test]
    fn reordered_multiallelic_genotype_mismatch_uses_truth_allele_indices() {
        // happy:chr21 at 38861935. Legacy unifies the reordered ALT columns
        // into truth order and emits one combined FN/FP row. Query `2/1`
        // against `TA,TAA` therefore becomes `1/2` against `TAA,TA`.
        let truth = variant(38_861_935, "T", "TAA,TA", "1|1");
        let query = variant(38_861_935, "T", "TA,TAA", "2/1");

        assert!(!query_matches_truth_key(&query, &truth));
        assert!(query_matches_truth_allele_set(&query, &truth));
        assert_ne!(
            selected_alt_sequences(&truth),
            selected_alt_sequences(&query)
        );
        assert_eq!(canonical_hetalt_gt(&truth.key.alt_allele, &query), "1/2");
        assert_eq!(compute_paired_bk(&truth, &query), "lm");
    }

    // Class A pin (post-#79): chr21:27249918 truth-subset match. Truth
    // `CTAAATAAA→C` GT 1|0 selects {C}; query `CTAAATAAA→C,CTAAA` GT 1/2
    // selects {C, CTAAA}. Truth's {C} ⊊ query's selected, the C allele
    // matches between sides, and the multi-allelic query primitive-splits
    // (CTAAA trims to ATAAA→A at pos+4, distinct from the C primitive at
    // pos). Legacy emits a combined TP/gm row at truth's representation
    // plus a residual FP at the orphan primitive.
    #[test]
    fn truth_subset_match_fires_for_chr21_27249918_shape() {
        let truth = variant(27249918, "CTAAATAAA", "C", "1|0");
        let query = variant(27249918, "CTAAATAAA", "C,CTAAA", "1/2");
        let reference = "N".repeat(27249930);
        assert!(truth_subset_match(
            &truth,
            &query,
            &reference,
            27249918,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
        // Query GT remap into truth's index space: query allele 1 (C)
        // matches truth idx 1; query allele 2 (CTAAA) drops to ref 0.
        // Unphased canonicalisation places the smaller index first.
        assert_eq!(remap_query_gt_subset(&truth, &query), "0/1");
    }

    #[test]
    fn truth_subset_match_rejects_homalt_truth_against_hetalt_query() {
        // chr21:10716541 shape: truth `C→G` GT 1|1 (homalt, multiset
        // [G×2]) vs query `C→A,G` GT 2/1 (hetalt, multiset [G×1, A×1]).
        // Set-subset would match {G} ⊆ {A,G} but the multiset check
        // rejects: truth needs G twice, query has it once.
        let truth = variant(10716541, "C", "G", "1|1");
        let query = variant(10716541, "C", "A,G", "2/1");
        let reference = "N".repeat(10716550);
        assert!(!truth_subset_match(
            &truth,
            &query,
            &reference,
            10716541,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
    }

    #[test]
    fn truth_subset_match_rejects_same_anchor_multiallelic_query() {
        // chr21:40096658 shape: truth `T→TAGATAGAG` GT 1|0 vs query
        // `T→TAGATAGAG,TAGATAGAT` GT 1/2. Both query alts trim to the
        // same anchor (T at pos), so `query_primitive_splits` is false
        // and legacy keeps two separate rows rather than emitting a
        // combined TP. The gate must reject this case.
        let truth = variant(40096658, "T", "TAGATAGAG", "1|0");
        let query = variant(40096658, "T", "TAGATAGAG,TAGATAGAT", "1/2");
        let reference = "N".repeat(40096670);
        assert!(!truth_subset_match(
            &truth,
            &query,
            &reference,
            40096658,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
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
                assert_eq!(
                    *anchor, 15671094,
                    "anchor must not slide for non-homopolymer insert"
                );
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
        reference[16997948] = b'G'; // pos 16997949: G
        reference[16997949] = b'C'; // pos 16997950: C (deleted)
        reference[16997950] = b'A'; // pos 16997951: A (deleted, insert anchor)
        reference[16997951] = b'G'; // pos 16997952: G
        let reference = String::from_utf8(reference).unwrap();
        let events = vec![
            Event::Delete {
                start: 16997950,
                end: 16997951,
            },
            Event::Insert {
                anchor: 16997951,
                seq: "CG".to_string(),
            },
        ];
        let result = apply_events(&reference, 16997949, 16997952, &events).unwrap();
        assert!(
            result.is_some(),
            "delete + downstream insert must produce a valid haplotype"
        );
    }

    // Class 4 pin: BI (comparison_info)
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

    #[test]
    fn symbolic_output_ref_uses_the_reference_base() {
        let variant = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 2,
                ref_allele: "N".to_string(),
                alt_allele: "<DEL>".to_string(),
            },
            qual: "0".to_string(),
            filter: ".".to_string(),
            gt: "1|0".to_string(),
        };
        assert_eq!(display_ref(&variant, "aTg"), "T");
    }

    #[test]
    fn fully_nonconf_matched_fanout_uses_local_match_bk() {
        assert_eq!(matched_query_unk_bk(true, false, "."), "lm");
        assert_eq!(matched_query_unk_bk(true, true, "."), ".");
        assert_eq!(matched_query_unk_bk(false, false, "."), ".");
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
        let intervals = vec![vcf::BedInterval {
            chrom: "chr21".to_string(),
            start: 15859657,
            end: 15859667,
        }];
        // Anchor at 15859667 is in [15859657,15859667)? 15859666 < 15859667 → yes.
        // Anchor+1 at 15859668 is in? 15859667 < 15859667 → NO.
        // Partial → pure insertion → skip CONF.
        assert!(!variant_is_conf(&var, "N", 15859600, 15859700, &intervals));
    }

    // Per-primitive CONF coverage: a multi-allelic
    // query whose deletion primitive sits inside CONF but whose insertion
    // primitive straddles a CONF edge must produce a RegionState where
    // `covered_query` contains the deletion primitive's key but NOT the
    // insertion primitive's, and `any_nonconf` is true. Without this,
    // the cluster collapses to TS_contained on every record and the
    // insertion primitive incorrectly carries a CONF tag.
    //
    // Mirrors chr21:48036437 — query `AGTGTGT → AGTGTGTGT,A` GT=1/2 at
    // pos 37002776 splits into a deletion at pos 37002776 (in CONF) and
    // an insertion T→TGT at pos 37002782 (straddles a CONF gap).
    #[test]
    fn region_state_marks_multi_allelic_primitives_separately() {
        let parent = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 100,
                ref_allele: "AGT".to_string(),
                alt_allele: "AGTGT,A".to_string(),
            },
            qual: "100".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 102,
            truth: vec![],
            query: vec![parent.clone()],
        };
        // CONF covers 1-based positions 100..101, leaving 102 uncovered.
        // The deletion primitive's effective ref range
        // (refstart=101, refend=102, is_pure_insertion=false) overlaps
        // CONF at position 101 → covered (subst/del has_overlap path).
        // The insertion primitive after right-anchor canonicalisation
        // emits at pos 102 with anchor T → bracket [101, 102]. Position
        // 102 is NOT in CONF → fully_covered=false → primitive must NOT
        // be marked covered.
        let intervals = vec![vcf::BedInterval {
            chrom: "chr21".to_string(),
            start: 99, // 0-based half-open → covers 1-based 100..101
            end: 101,
        }];
        let state = RegionState::from_cluster(&cluster, "N", Some(&intervals));
        assert!(state.any_conf, "deletion primitive must register coverage");
        // The insertion primitive's anchor falls at the CONF edge with
        // anchor+1 outside coverage — primitive must NOT be in
        // covered_query, and any_nonconf must be set.
        let primitives = split_query_primitives_with_neighbors(&parent, "N", 100, &[], &[]);
        let insertion_primitive = primitives
            .iter()
            .find(|p| p.key.alt_allele.len() > p.key.ref_allele.len())
            .expect("split_query_primitives must produce one insertion");
        assert!(
            !state.covered_query.contains(&insertion_primitive.key),
            "insertion primitive at CONF edge must NOT be in covered_query"
        );
        assert!(
            state.any_nonconf,
            "presence of an uncovered insertion primitive must mark cluster as non-CONF"
        );
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
        let intervals = vec![vcf::BedInterval {
            chrom: "chr21".to_string(),
            start: 15859480,
            end: 15859645,
        }];
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
        let padding = gvcf2bed_padding(&truth, None);
        assert_eq!(padding.len(), 1);
        assert_eq!(padding[0].chrom, "chr21");
        // 0-based half-open: anchor = 17562904 .. anchor+1+1 = 17562906
        assert_eq!(padding[0].start, 17562904);
        assert_eq!(padding[0].end, 17562906);
    }

    // gvcf2bed `-T <bed>` filter: legacy uses `bcf_sr_set_targets(.., 1, 0)`
    // which gates on the **start position only** — not the full ref span.
    // A record whose 1-based pos (→ 0-based start) lies outside every
    // raw CONF interval is dropped before emission. This is what
    // `IS_CONF.Size` parity hinges on (legacy sums per-file BED lengths
    // without cross-file dedup, so dropping out-of-target records keeps
    // the padding budget honest).
    #[test]
    fn gvcf2bed_padding_target_filter_excludes_out_of_target_record() {
        let truth = vec![
            // pos 100 → pos_0b 99, INSIDE conf [50, 150)
            Variant {
                key: VariantKey {
                    chrom: "chr1".to_string(),
                    pos: 100,
                    ref_allele: "A".to_string(),
                    alt_allele: "G".to_string(),
                },
                qual: ".".to_string(),
                filter: "PASS".to_string(),
                gt: "0/1".to_string(),
            },
            // pos 200 → pos_0b 199, OUTSIDE conf — should be dropped
            Variant {
                key: VariantKey {
                    chrom: "chr1".to_string(),
                    pos: 200,
                    ref_allele: "A".to_string(),
                    alt_allele: "G".to_string(),
                },
                qual: ".".to_string(),
                filter: "PASS".to_string(),
                gt: "0/1".to_string(),
            },
        ];
        let conf = vec![vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 50,
            end: 150,
        }];
        let padding = gvcf2bed_padding(&truth, Some(&conf));
        assert_eq!(padding.len(), 1, "only in-target record should emit");
        assert_eq!(padding[0].start, 99);
        assert_eq!(padding[0].end, 100);
    }

    // Legacy gvcf2bed emits a BED line **per record** unconditionally,
    // even when every alt is symbolic (`<DEL>`, `<NON_REF>`, etc.). The
    // alt loop's `break` on the first non-NUC alt leaves
    // `nuc_alleles=false`, so refstart/refend stay at the raw
    // [pos, pos+reflen-1] from getLocation — and emission proceeds. Our
    // truth fixtures contain ~14 such records on chr21 (`<DEL>` calls);
    // skipping them under-counts IS_CONF.Size by ~14 bp.
    #[test]
    fn gvcf2bed_padding_emits_symbolic_only_record_with_raw_ref_span() {
        let truth = vec![Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15847471, // 1-based — pos_0b = 15847470
                ref_allele: "N".to_string(),
                alt_allele: "<DEL>".to_string(),
            },
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            gt: "0/1".to_string(),
        }];
        let padding = gvcf2bed_padding(&truth, None);
        assert_eq!(padding.len(), 1, "symbolic-only record must still emit");
        // Raw refrange [pos_0b, pos_0b + reflen - 1] = [15847470, 15847470].
        // Half-open BED: [15847470, 15847471) → 1 bp.
        assert_eq!(padding[0].start, 15847470);
        assert_eq!(padding[0].end, 15847471);
    }

    #[test]
    fn gvcf2bed_padding_preserves_preprocessed_truth_spans() {
        let truth = [
            (10, "C", "A"),
            (20, "T", "A"),
            (29, "ACGTACG", "A"),
            (40, "T", "."),
            (50, "CG", "C"),
        ]
        .into_iter()
        .map(|(pos, reference, alternate)| Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos,
                ref_allele: reference.to_string(),
                alt_allele: alternate.to_string(),
            },
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            gt: String::new(),
        })
        .collect::<Vec<_>>();
        let confidence = [vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 0,
            end: 120,
        }];
        let padding = gvcf2bed_padding(&truth, Some(&confidence));
        let size = padding
            .iter()
            .map(|interval| interval.end - interval.start)
            .sum::<usize>();

        assert_eq!(size, 10);
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
        let queries = [
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
        let result =
            cluster_signature(&cluster, &variants, &reference, None, &BTreeSet::new()).unwrap();
        assert!(
            result.is_none(),
            "overlapping homalt deletions must produce Ok(None)"
        );
    }

    /// Pin Class G: when `cluster_query_filter` aggregates filter tokens
    /// across multiple query records (truth-side TP-row stamping path),
    /// the joined string must be byte-wise sorted to match legacy's
    /// bcftools-merged ordering. chr21:40875336 cluster sources are:
    ///   * pos 40875343 T→A: `TruthSensitivityTranche99.90to100.00;LowGQX`
    ///   * pos 40875344 T→A: `TruthSensitivityTranche99.00to99.90`
    ///   * pos 40875347 A→G: `TruthSensitivityTranche99.90to100.00`
    ///
    /// Source-order union is `T99.90to100.00;LowGQX;T99.00to99.90` (rust
    /// pre-fix). Legacy emits `LowGQX;T99.00to99.90;T99.90to100.00` —
    /// alphabetic.
    #[test]
    fn cluster_query_filter_sorts_aggregated_tokens() {
        let mk = |pos: usize, filter: &str| Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos,
                ref_allele: "T".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: filter.to_string(),
            gt: "0/1".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 40875336,
            end: 40875347,
            truth: vec![],
            query: vec![
                mk(40875343, "TruthSensitivityTranche99.90to100.00;LowGQX"),
                mk(40875344, "TruthSensitivityTranche99.00to99.90"),
                mk(40875347, "TruthSensitivityTranche99.90to100.00"),
            ],
        };
        let got = cluster_query_filter(&cluster);
        assert_eq!(
            got,
            "LowGQX;TruthSensitivityTranche99.00to99.90;TruthSensitivityTranche99.90to100.00"
        );
    }

    /// Empty cluster (no PASS-bearing query) must yield ".".
    #[test]
    fn cluster_query_filter_empty_returns_dot() {
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 0,
            end: 0,
            truth: vec![],
            query: vec![],
        };
        assert_eq!(cluster_query_filter(&cluster), ".");
    }

    /// Pin Class E: chr21:44049606 cluster (chr21_passonly shape) — query
    /// has a homalt deletion AATGATAGATAG→A at 44049606 covering positions
    /// 44049607..44049617, plus a 1/2 multi-allelic at 44049615 whose
    /// second alt registers as an "insert" in `query_insert_conflict_…`.
    /// The deletion blocks the insert anchor on both haplotypes, draining
    /// query enumeration. Truth's only record (TGATA→T at 44049663) is
    /// far outside the deletion's claimed range, so legacy treats this as
    /// an internal query conflict — `BK=.`. Pre-fix the loose
    /// `!truth_remaining.is_empty()` disjunct made this fire `Some(false)`
    /// (→ BK=lm); the tightened gate must return `None`.
    #[test]
    fn deletion_covers_insert_no_proximate_truth_returns_none() {
        let q_del = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 44049606,
                ref_allele: "AATGATAGATAG".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let q_multi = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 44049615,
                ref_allele: "TAGATGATAGAT".to_string(),
                alt_allele: "T,TAGACAGATGATAGAT".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let truth_far = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 44049663,
                ref_allele: "TGATA".to_string(),
                alt_allele: "T".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "0|1".to_string(),
        };
        let query = vec![q_del, q_multi];
        let truth = vec![truth_far.clone()];
        // truth_remaining still contains the unmatched 44049663 record.
        let result = query_insert_conflict_has_truth_counterpart(&query, &truth, &truth);
        assert_eq!(
            result, None,
            "deletion-covers-insert with truth outside the deletion range \
             must return None so hap_mismatch stays false"
        );
    }

    /// Positive pin for Class E gate: when truth has a variant at the
    /// blocked insert anchor, `Some(false)` must still fire so the BK=lm
    /// branch keeps working for genuinely truth-anchored mismatches.
    #[test]
    fn deletion_covers_insert_truth_at_anchor_returns_some_false() {
        let q_del = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 100,
                ref_allele: "ATGATGATGAT".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        // Multi-allelic with a 16-base alt → registers as `insert` at pos 105.
        let q_multi = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 105,
                ref_allele: "GATGAT".to_string(),
                alt_allele: "G,GATCATGAT".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let truth_at_anchor = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 105,
                ref_allele: "G".to_string(),
                alt_allele: "GATC".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "0|1".to_string(),
        };
        let query = vec![q_del, q_multi];
        let truth = vec![truth_at_anchor];
        let result = query_insert_conflict_has_truth_counterpart(&query, &truth, &truth);
        assert_eq!(
            result,
            Some(false),
            "truth at the blocked anchor must keep the BK=lm path firing"
        );
    }

    /// Positive pin for Class E gate: when truth_remaining contains a
    /// variant that overlaps the blocking deletion's claimed range, the
    /// drain represents a genuine mismatch — `Some(false)` must fire.
    #[test]
    fn deletion_covers_insert_truth_in_del_range_returns_some_false() {
        let q_del = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 100,
                ref_allele: "ATGATGATGAT".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let q_multi = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 105,
                ref_allele: "GATGAT".to_string(),
                alt_allele: "G,GATCATGAT".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        // Truth variant at pos 107 — inside the deletion's [101, 110] range.
        let truth_in_range = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 107,
                ref_allele: "T".to_string(),
                alt_allele: "C".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "0|1".to_string(),
        };
        let query = vec![q_del, q_multi];
        let truth = vec![truth_in_range];
        let result = query_insert_conflict_has_truth_counterpart(&query, &truth, &truth);
        assert_eq!(
            result,
            Some(false),
            "unmatched truth inside the blocking deletion range must keep \
             BK=lm firing"
        );
    }

    /// Class B (chr21:21690513). Reproduces the chr21 reference window
    /// `ttatatatatatatatatatatacacacacacacacatacatacatacata` at synthetic
    /// 1-based positions 1..51 (pos 1 ↔ chr21:21690480). Anchor at pos 34
    /// corresponds to chr21:21690513 (`C`). The CA-microsat upstream lets
    /// `partial_credit::left_shift` canonicalize the CACAC primitive at
    /// pos 22 (`T→TACAC`) — distinct from CACAT's stayed-put pos 34.
    /// Truth declares `C→CACAT` at pos 34 only, so the shifted CACAC has
    /// no truth representation → fan out.
    #[test]
    fn class_b_same_anchor_insertion_fans_out_when_truth_at_original_only() {
        let reference = b"ttatatatatatatatatatatacacacacacacacatacatacatacata";
        let trimmed = vec![
            (34, "C".to_string(), "CACAC".to_string()),
            (34, "C".to_string(), "CACAT".to_string()),
        ];
        let cluster_truth = vec![variant(34, "C", "CACAT", "0|1")];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            34,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &[],
        );
        let shifted = result.expect("must fan out — truth at original, orphan at shifted");
        assert_eq!(shifted.len(), 2);
        let cacac = shifted
            .iter()
            .find(|(_, _, alt)| alt.ends_with('C'))
            .expect("CACAC primitive present");
        let cacat = shifted
            .iter()
            .find(|(_, _, alt)| alt.ends_with('T'))
            .expect("CACAT primitive present");
        assert_eq!(
            cacac.0, 22,
            "CACAC must slide through CA-microsat to pos 22"
        );
        assert_eq!(cacac.1, "T", "ref byte at the shifted anchor (pos 22) is T");
        assert_eq!(cacac.2, "TACAC", "alt rotates to T-prefixed canonical form");
        assert_eq!(
            cacat.0, 34,
            "CACAT must stay at the original anchor (pos 34)"
        );
        assert_eq!(cacat.1, "C");
        assert_eq!(cacat.2, "CACAT");
    }

    /// Class B negative — chr21:40096658 shape. The TAGATAGAT primitive
    /// canonicalizes via `left_shift` to a position where truth ALREADY
    /// declares that allele as part of a multi-allelic. The fan-out
    /// gate must suppress the split so the block-level haplotype matcher
    /// can reconcile the multi-allelic record. Reproduces the real chr21
    /// AGAT-microsat upstream of pos 40096658.
    #[test]
    fn class_b_same_anchor_insertion_blocked_when_truth_at_shifted_anchor() {
        // Synthetic positions 1..50: pos 1 ↔ chr21:40096640. Pos 19 ↔
        // chr21:40096658 (anchor `T`). Pos 11 ↔ chr21:40096650 (anchor
        // `C`) where TAGATAGAT shifts to `CAGATAGAT`.
        let reference = b"ctgaagagttcagatagatagatagatagatagatagatagatagacaga";
        let trimmed = vec![
            (19, "T".to_string(), "TAGATAGAG".to_string()),
            (19, "T".to_string(), "TAGATAGAT".to_string()),
        ];
        // Truth declares TAGATAGAG at the original anchor AND
        // CAGATAGAT,CAGATAGATAGAT at the shifted anchor.
        let cluster_truth = vec![
            variant(11, "C", "CAGATAGAT,CAGATAGATAGAT", "0|1"),
            variant(19, "T", "TAGATAGAG", "1|0"),
        ];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            19,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &[],
        );
        assert!(
            result.is_none(),
            "fan-out must be blocked when truth declares the shifted alt"
        );
    }

    /// Class B negative — chr21:32767041 shape. No truth records exist in
    /// the cluster. Even though the primitives shift apart (TCTCTCT slides
    /// to a CTCT-microsat anchor, TCTCACA stays put), the fan-out must be
    /// suppressed because no query alt has truth representation at the
    /// original anchor — there's no truth_subset_match-style emit shape
    /// to license the split.
    #[test]
    fn class_b_same_anchor_insertion_blocked_without_truth_at_original() {
        let reference = b"ttatatatatatatatatatatacacacacacacacatacatacatacata";
        let trimmed = vec![
            (34, "C".to_string(), "CACAC".to_string()),
            (34, "C".to_string(), "CACAT".to_string()),
        ];
        let cluster_truth: Vec<Variant> = vec![];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            34,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &[],
        );
        assert!(
            result.is_none(),
            "fan-out must be blocked when no query alt has truth at the original anchor"
        );
    }

    /// Class B (chr21:21690513 chr21 case). When a neighboring cluster
    /// query record sits at the natural slide target, the slide must
    /// stop one position above so the shifted primitive doesn't share
    /// an anchor with the existing record. Pos 22 holds a SNP neighbor;
    /// the CACAC primitive must canonicalize at pos 23 (`A→ACACA`) — the
    /// same legacy verdict reproduced in the chr21 (no --pass-only) case.
    #[test]
    fn class_b_same_anchor_neighbor_floor_clamps_slide_target() {
        let reference = b"ttatatatatatatatatatatacacacacacacacatacatacatacata";
        let trimmed = vec![
            (34, "C".to_string(), "CACAC".to_string()),
            (34, "C".to_string(), "CACAT".to_string()),
        ];
        let cluster_truth = vec![variant(34, "C", "CACAT", "0|1")];
        // Neighboring SNP at pos 22 (`T→C`) — sliding CACAC onto pos 22
        // would clobber that record's anchor.
        let cluster_neighbors = vec![variant(22, "T", "C", "0/1")];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            34,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &cluster_neighbors,
        );
        let shifted = result.expect("must still fan out — neighbor only clamps slide depth");
        // The CACAC primitive rotates to (23, A, ACACA) — the slide
        // reduces start by 1 per iteration, alternating the alt's last
        // base. With pos_min clamped to 22 (after pure-insertion bump
        // becomes 23), the slide stops with start = 23 holding 'A' anchor.
        let stayer = shifted
            .iter()
            .find(|(p, _, _)| *p == 34)
            .expect("CACAT primitive must stay at pos 34");
        let shifter = shifted
            .iter()
            .find(|(p, _, _)| *p != 34)
            .expect("shifted primitive must land at a distinct anchor");
        assert_eq!(stayer.1, "C");
        assert_eq!(stayer.2, "CACAT");
        assert_eq!(
            shifter.0, 23,
            "CACAC slide must stop at pos 23 (one above neighbor at pos 22)"
        );
        assert_eq!(shifter.1, "A", "ref byte at pos 23 is A");
        assert_eq!(
            shifter.2, "ACACA",
            "alt rotates to A-prefixed canonical form"
        );
    }

    /// Class F (chr21:47906004). A multi-allelic deletion's parent record
    /// has an `effective_refrange` that reaches into a CONF interval, but
    /// neither fanned-out primitive's range (after the per-primitive
    /// left-shift) touches CONF. Legacy emits no Regions tag on the
    /// per-primitive rows because it operates on post-fan-out records
    /// only — the parent-path's any_conf vote must be suppressed for
    /// fanned-out multi-allelics. Without this gate the cluster picks
    /// up a spurious TS_boundary tag from the parent.
    ///
    /// Synthetic layout (mirrors chr21:47906xxx at smaller positions):
    /// reference `aaaaaaaaaaaaaaaaaaaaagaactaaagt` covers 1-based positions
    /// 1..31. The variant `AGAACTAAA→A,AAAA` at pos 21 produces
    /// per-primitive rows at (17, AAAAGAACT, A) and (21, AGAACT, A) after
    /// the deletion-only slide (range 18..25 and 22..26). The parent's
    /// effective_refrange is 21..28 (alt A reaches further right than
    /// either fanned-out primitive does).
    #[test]
    fn class_f_region_state_skips_parent_path_for_fanned_out_multiallelic() {
        let reference = "aaaaaaaaaaaaaaaaaaaaagaactaaagt".to_string();
        let parent = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 21,
                ref_allele: "AGAACTAAA".to_string(),
                alt_allele: "A,AAAA".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 21,
            end: 30,
            truth: vec![],
            query: vec![parent],
        };
        // CONF covers 1-based positions 27..30 — the parent's effective
        // range reaches into 27..28 but the fanned-out primitives' ranges
        // (post-slide 18..25 and 22..26) both stop at or before pos 26.
        let intervals = vec![vcf::BedInterval {
            chrom: "chr21".to_string(),
            start: 26,
            end: 30,
        }];
        let state = RegionState::from_cluster(&cluster, &reference, Some(&intervals));
        assert!(
            !state.any_conf,
            "fanned-out multi-allelic with all primitives outside CONF must NOT \
             register any_conf via the parent path (got any_conf=true)"
        );
        assert!(
            state.any_nonconf,
            "primitives outside CONF must vote any_nonconf"
        );
    }
}
