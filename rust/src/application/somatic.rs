use crate::adapters::report::suffixed_report_path;
use crate::application::{SomaticArgs, ValidatedSomaticArgs};
use crate::domain::{Interval, RawVcfRecord};
use crate::{
    adapters::{fasta, vcf},
    application::ftx,
    engines::strelka,
    output::OutputTransaction,
};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

mod allele_frequency;
mod features;
mod metrics;
mod normalization;
mod reports;

#[cfg(test)]
mod test_suite;

use allele_frequency::{format_af_interval, parse_af_bins};
use features::*;
use metrics::*;
use normalization::*;
use reports::*;

/// Version string emitted in the `sompyversion` column and JSON metadata.
/// Matches legacy som.py when `Haplo.version.__version__` is empty — the
/// pinned reference container resolves this to `"som.py-"` (trailing hyphen with
/// nothing after). We reproduce the exact literal so stats.csv byte-matches.
const SOM_VERSION: &str = "som.py-";
const MAX_AF_BINS: usize = 100;
const SOMATIC_ROC_CHUNK: usize = 16_384;
const SOMATIC_ROC_MERGE_FAN_IN: usize = 32;
const STATS_TYPE_ROWS: [(usize, &str); 4] =
    [(0, "indels"), (1, "SNVs"), (6, "MNPs"), (7, "others")];
static SOMATIC_SCRATCH_RUN_ID: AtomicU64 = AtomicU64::new(0);

struct SomaticScratch {
    path: PathBuf,
    keep: bool,
}

impl SomaticScratch {
    fn create(prefix: Option<&str>, keep_requested: bool) -> Result<Self> {
        if let Some(prefix) = prefix {
            let path = PathBuf::from(prefix);
            fs::create_dir_all(&path).with_context(|| {
                format!(
                    "failed to create somatic scratch directory {}",
                    path.display()
                )
            })?;
            // Legacy som.py always retains an explicit scratch prefix so it
            // can be consumed by a later --continue invocation.
            return Ok(Self { path, keep: true });
        }

        let parent = std::env::temp_dir().join("hap-somatic");
        fs::create_dir_all(&parent)
            .with_context(|| format!("failed to create scratch parent {}", parent.display()))?;
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap_or_default()
            .as_nanos();
        for _ in 0..128 {
            let id = SOMATIC_SCRATCH_RUN_ID.fetch_add(1, Ordering::Relaxed);
            let path = parent.join(format!("run-{}-{timestamp}-{id}", std::process::id()));
            match fs::create_dir(&path) {
                Ok(()) => {
                    return Ok(Self {
                        path,
                        keep: keep_requested,
                    });
                }
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => {
                    return Err(error).with_context(|| {
                        format!(
                            "failed to create somatic scratch directory {}",
                            path.display()
                        )
                    });
                }
            }
        }
        bail!(
            "failed to allocate a unique somatic scratch directory under {}",
            parent.display()
        )
    }

    fn cleanup(mut self) -> Result<()> {
        if self.keep {
            return Ok(());
        }
        fs::remove_dir_all(&self.path).with_context(|| {
            format!(
                "failed to remove somatic scratch directory {}",
                self.path.display()
            )
        })?;
        self.keep = true;
        Ok(())
    }
}

impl Drop for SomaticScratch {
    fn drop(&mut self) {
        if !self.keep {
            let _ = fs::remove_dir_all(&self.path);
        }
    }
}

struct SomaticOperationalControls {
    scratch: SomaticScratch,
    logfile: Option<fs::File>,
    verbose: bool,
    quiet: bool,
}

impl SomaticOperationalControls {
    fn prepare(args: &SomaticArgs) -> Result<Self> {
        let logfile = args
            .logfile
            .as_deref()
            .map(|logfile| {
                let path = Path::new(logfile);
                if let Some(parent) = path.parent()
                    && !parent.as_os_str().is_empty()
                {
                    fs::create_dir_all(parent).with_context(|| {
                        format!("failed to create log directory {}", parent.display())
                    })?;
                }
                fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(path)
                    .with_context(|| format!("failed to open somatic log file {}", path.display()))
            })
            .transpose()?;
        let scratch = SomaticScratch::create(args.scratch_prefix.as_deref(), args.keep_scratch)?;
        Ok(Self {
            scratch,
            logfile,
            verbose: args.verbose,
            quiet: args.quiet,
        })
    }

    fn info(&mut self, message: &str) -> Result<()> {
        if !self.verbose {
            return Ok(());
        }
        let line = format!("INFO     {message}\n");
        if let Some(logfile) = self.logfile.as_mut() {
            logfile
                .write_all(line.as_bytes())
                .context("failed to write somatic log")?;
            logfile.flush().context("failed to flush somatic log")?;
        } else {
            eprint!("{line}");
        }
        Ok(())
    }

    fn print_default_summary(&self, lines: &[String]) {
        if self.should_print_summary() {
            println!("\n{}", lines.join("\n"));
        }
    }

    fn should_print_summary(&self) -> bool {
        !self.quiet && !self.verbose
    }
}

#[derive(Clone, Copy, Debug)]
struct SomaticRocConfig {
    feature_table: &'static str,
    score: &'static str,
    filter_name: Option<&'static str>,
    zero_score_unless_nt_ref: bool,
}

fn somatic_roc_config(name: &str) -> Option<SomaticRocConfig> {
    let config = match name {
        "strelka.snv.qss" => SomaticRocConfig {
            feature_table: "hcc.strelka.snv",
            score: "QSS_NT",
            filter_name: Some("QSS_ref"),
            zero_score_unless_nt_ref: true,
        },
        "strelka.snv.vqsr" => SomaticRocConfig {
            feature_table: "hcc.strelka.snv",
            score: "VQSR",
            filter_name: Some("LowQscore"),
            zero_score_unless_nt_ref: true,
        },
        "strelka.snv" => SomaticRocConfig {
            feature_table: "hcc.strelka.snv",
            score: "EVS",
            filter_name: Some("LowEVS"),
            zero_score_unless_nt_ref: true,
        },
        "strelka.indel" => SomaticRocConfig {
            feature_table: "hcc.strelka.indel",
            score: "QSI_NT",
            filter_name: Some("QSI_ref"),
            zero_score_unless_nt_ref: true,
        },
        "strelka.indel.evs" => SomaticRocConfig {
            feature_table: "hcc.strelka.indel",
            score: "EVS",
            filter_name: Some("LowEVS"),
            zero_score_unless_nt_ref: false,
        },
        "varscan2.snv" => SomaticRocConfig {
            feature_table: "hcc.varscan2.snv",
            score: "SSC",
            filter_name: None,
            zero_score_unless_nt_ref: false,
        },
        "varscan2.indel" => SomaticRocConfig {
            feature_table: "hcc.varscan2.indel",
            score: "SSC",
            filter_name: None,
            zero_score_unless_nt_ref: false,
        },
        "mutect.snv" => SomaticRocConfig {
            feature_table: "hcc.mutect.snv",
            score: "TLOD",
            filter_name: Some("t_lod_fstar"),
            zero_score_unless_nt_ref: false,
        },
        "mutect.indel" => SomaticRocConfig {
            feature_table: "hcc.mutect.indel",
            score: "TLOD",
            filter_name: Some("t_lod_fstar"),
            zero_score_unless_nt_ref: false,
        },
        _ => return None,
    };
    Some(config)
}

#[derive(Clone, Copy, Debug, Default)]
struct SomaticCounts {
    truth_total: usize,
    query_total: usize,
    tp: usize,
    fp: usize,
    fn_count: usize,
    unk: usize,
    ambi: usize,
}

#[derive(Clone, Debug)]
struct FilteredRawRecord {
    record: RawVcfRecord,
}

#[derive(Clone, Copy, Debug, Default)]
struct FilteredCounts {
    tp: usize,
    fp: usize,
    unk: usize,
    ambi: usize,
}

#[derive(Clone, Copy, Debug)]
struct StatsRowContext<'a> {
    fp_region_size: usize,
    ci_alpha: f64,
    filtered: Option<FilteredCounts>,
    include_filtered_columns: bool,
    commandline: &'a str,
}

#[derive(Clone, Debug)]
struct AmbiguousInterval {
    interval: Interval,
    label: String,
    details: Vec<String>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum QueryClass {
    Fp,
    Unk,
    Ambi,
}

pub(crate) fn run(args: ValidatedSomaticArgs) -> Result<()> {
    let mut args = args.into_inner();
    let destination_prefix = PathBuf::from(&args.output);
    let inputs = somatic_inputs(&args);
    let logfile = args.logfile.as_ref().map(PathBuf::from);
    let artifacts =
        somatic_artifacts(&args).context("failed to plan somatic AF-bin ROC artifacts")?;
    let transaction = OutputTransaction::family(&inputs, &destination_prefix, artifacts)?
        .with_files(logfile.iter())?;
    if let Some(path) = logfile.as_deref() {
        args.logfile = Some(
            transaction
                .staged_file(path)?
                .to_string_lossy()
                .into_owned(),
        );
    }
    args.output = transaction.staged_prefix()?.to_string_lossy().into_owned();
    run_inner(args).map_err(|error| {
        anyhow::anyhow!(
            "failed to produce somatic report generation {}: {error:#}",
            destination_prefix.display()
        )
    })?;
    transaction.commit()
}

fn somatic_artifacts(args: &SomaticArgs) -> Result<Vec<String>> {
    let mut artifacts = [
        "ambiclasses.csv",
        "ambireasons.csv",
        "features.csv",
        "roc.csv",
        "stats.csv",
        "metrics.json",
        "summary.csv",
        "extended.csv",
    ]
    .into_iter()
    .map(str::to_string)
    .collect::<Vec<_>>();
    if args.af_strat {
        for (start, end) in parse_af_bins(&args.af_strat_binsize)? {
            let interval = format_af_interval(start, end);
            for prefix in ["records", "SNVs", "indels"] {
                artifacts.push(format!("{prefix}.{interval}.roc.csv"));
            }
        }
    }
    Ok(artifacts)
}

fn somatic_inputs(args: &SomaticArgs) -> Vec<PathBuf> {
    std::iter::once(args.truth.as_str())
        .chain(std::iter::once(args.query.as_str()))
        .chain(std::iter::once(args.reference.as_str()))
        .chain(args.regions_bedfile.as_deref())
        .chain(args.targets_bedfile.as_deref())
        .chain(args.fp_bedfile.as_deref())
        .chain(args.ambiguous_beds.iter().map(String::as_str))
        .chain(args.bams.iter().map(String::as_str))
        .map(PathBuf::from)
        .collect()
}

fn run_inner(mut args: SomaticArgs) -> Result<()> {
    if let Some(config) = args.roc.as_deref().and_then(somatic_roc_config) {
        args.feature_table = Some(config.feature_table.to_string());
    }
    validate_args(&args)?;
    let af_bins = if args.af_strat {
        parse_af_bins(&args.af_strat_binsize)?
    } else {
        Vec::new()
    };
    let mut controls = SomaticOperationalControls::prepare(&args)?;
    controls.info(&format!(
        "Scratch path is {}",
        controls.scratch.path.display()
    ))?;
    let count_unk = resolve_toggle(args.count_unk, args.no_count_unk);
    let ambi_fp = resolve_toggle(args.ambi_fp, args.no_ambi_fp);
    let fixchr_truth = args.fixchr_truth.unwrap_or(true) && !args.no_fixchr_truth;
    let fixchr_query = args.fixchr_query.unwrap_or(true) && !args.no_fixchr_query;
    let ci_alpha = 1.0 - args.ci_level;
    let (normalize_truth, normalize_query) = selected_normalizations(&args);

    // Legacy som.py opens the reference only for bcftools normalization or
    // when it must derive an automatic reference-sized FP denominator. Plain
    // allele comparison with an explicit/FP-BED denominator must therefore
    // remain usable even when the default hg19 path is absent.
    let reference_sequences = if normalize_truth || normalize_query {
        Some(fasta::read_sequences(Path::new(&args.reference))?)
    } else {
        None
    };
    let mut reference_lengths: BTreeMap<String, usize> = reference_sequences
        .as_ref()
        .map(|sequences| {
            sequences
                .iter()
                .map(|(name, sequence)| (name.clone(), sequence.len()))
                .collect()
        })
        .unwrap_or_default();
    let reference_contigs: BTreeSet<String> = reference_sequences
        .as_ref()
        .map(|sequences| sequences.keys().cloned().collect())
        .unwrap_or_default();
    let literal_contigs = BTreeSet::new();
    let regions = args
        .regions_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &literal_contigs))
        .transpose()?;
    let targets = args
        .targets_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &literal_contigs))
        .transpose()?;
    let fp_regions = args
        .fp_bedfile
        .as_ref()
        .map(|path| load_classification_bed(Path::new(path), &reference_contigs, fixchr_truth))
        .transpose()?;
    let ambiguous_regions =
        load_ambiguous_beds(&args.ambiguous_beds, &reference_contigs, fixchr_truth)?;
    let mut explanation_regions = ambiguous_regions.clone();
    if args.explain_ambiguous
        && !args.ambiguous_beds.is_empty()
        && let Some(path) = args.fp_bedfile.as_deref()
    {
        explanation_regions.extend(load_fp_explanation_bed(
            Path::new(path),
            &reference_contigs,
            fixchr_truth,
        )?);
    }
    let locations = args
        .location
        .as_deref()
        .map(|text| vcf::parse_locations(text, &literal_contigs))
        .transpose()?;

    controls.info("Normalizing/reading inputs")?;
    let truth_cache = controls.scratch.path.join("normalized_truth.vcf.gz");
    let query_cache = controls.scratch.path.join("normalized_query.vcf.gz");
    let reuse_truth = args.cont && truth_cache.exists();
    let reuse_query = args.cont && query_cache.exists();
    if reuse_truth {
        controls.info(&format!("Continuing from {}", truth_cache.display()))?;
    }
    if reuse_query {
        controls.info(&format!("Continuing from {}", query_cache.display()))?;
    }
    let truth_source = Path::new(&args.truth);
    let query_source = Path::new(&args.query);
    let truth_headers = vcf::open_validated_vcf(if reuse_truth {
        &truth_cache
    } else {
        truth_source
    })?
    .headers()
    .to_vec();
    let query_headers = vcf::open_validated_vcf(if reuse_query {
        &query_cache
    } else {
        query_source
    })?
    .headers()
    .to_vec();
    let bam_depths = ftx::bam_normalization_depths(&args.bams)?;
    if !reuse_truth {
        prepare_somatic_cache(
            truth_source,
            &truth_cache,
            &truth_headers,
            normalize_truth,
            reference_sequences.as_ref(),
        )?;
    }
    if !reuse_query {
        prepare_somatic_cache(
            query_source,
            &query_cache,
            &query_headers,
            normalize_query,
            reference_sequences.as_ref(),
        )?;
    }
    let truth_spools = spool_filtered_contigs(
        &truth_cache,
        Path::new(&args.truth),
        &RawFilterOptions {
            reference_contigs: &reference_contigs,
            fixchr: fixchr_truth,
            pass_only: true,
            regions: regions.as_deref(),
            targets: targets.as_deref(),
            locations: locations.as_deref(),
        },
    )?;
    let query_spools = spool_filtered_contigs(
        &query_cache,
        Path::new(&args.query),
        &RawFilterOptions {
            reference_contigs: &reference_contigs,
            fixchr: fixchr_query,
            pass_only: !args.include_nonpass,
            regions: regions.as_deref(),
            targets: targets.as_deref(),
            locations: locations.as_deref(),
        },
    )?;
    let query_depths = if bam_depths.is_empty() {
        strelka::parse_depths(&query_headers)
    } else {
        bam_depths.clone()
    };

    let mut by_type: BTreeMap<&'static str, SomaticCounts> = BTreeMap::new();
    let mut filtered_by_type: BTreeMap<&'static str, FilteredCounts> = BTreeMap::new();
    let mut filtered_records = FilteredCounts::default();
    let mut record_counts = SomaticCounts {
        truth_total: truth_spools.iter().map(|spool| spool.count).sum(),
        query_total: query_spools.iter().map(|spool| spool.count).sum(),
        ..SomaticCounts::default()
    };
    let feature_table_name = args.feature_table.as_deref().unwrap_or("");
    let use_strelka_hcc_indel = feature_table_name == "hcc.strelka.indel";
    let use_caller_feature_table =
        !feature_table_name.is_empty() && feature_table_name != "generic" && !use_strelka_hcc_indel;
    let use_generic_feature_table = feature_table_name == "generic";
    let mut feature_header = if use_strelka_hcc_indel {
        Some(STRELKA_HCC_INDEL_HEADER.to_string())
    } else if use_generic_feature_table {
        Some(
            ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,FILTER,FILTER.truth,QUAL,QUAL.truth"
                .to_string(),
        )
    } else {
        None
    };
    let mut feature_rows = FeatureRowSpools::new()?;
    let mut ambiguous_classes = BTreeMap::new();
    let mut ambiguous_reasons = BTreeMap::new();
    let truth_spool_index = truth_spools
        .iter()
        .enumerate()
        .map(|(index, spool)| (spool.chrom.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let query_spool_index = query_spools
        .iter()
        .enumerate()
        .map(|(index, spool)| (spool.chrom.clone(), index))
        .collect::<BTreeMap<_, _>>();
    let mut ordered_contigs = truth_spools
        .iter()
        .map(|spool| spool.chrom.clone())
        .collect::<Vec<_>>();
    ordered_contigs.extend(
        query_spools
            .iter()
            .map(|spool| spool.chrom.clone())
            .filter(|chrom| !truth_spool_index.contains_key(chrom)),
    );
    let truth_contigs = truth_spool_index.keys().cloned().collect::<BTreeSet<_>>();
    for chrom in ordered_contigs {
        let truth_raw_filtered = truth_spool_index
            .get(&chrom)
            .map(|index| truth_spools[*index].load())
            .transpose()?
            .unwrap_or_default();
        let query_raw_filtered = query_spool_index
            .get(&chrom)
            .map(|index| query_spools[*index].load())
            .transpose()?
            .unwrap_or_default();
        let (truth_matches, query_matches) =
            pair_exact_records(&truth_raw_filtered, &query_raw_filtered);
        let matched = truth_matches.iter().filter(|entry| entry.is_some()).count();
        record_counts.tp += matched;
        record_counts.fn_count += truth_raw_filtered.len().saturating_sub(matched);
        filtered_records.tp += query_raw_filtered
            .iter()
            .zip(&query_matches)
            .filter(|(record, matched)| matched.is_some() && !record.record.is_pass())
            .count();
        let mut caller_tp_truth = Vec::new();
        let mut caller_tp_query = Vec::new();
        let mut caller_fn_truth = Vec::new();
        let mut caller_fp_query = Vec::new();
        let mut caller_ambi_query = Vec::new();
        let mut caller_unk_query = Vec::new();
        for (truth_index, truth_record) in truth_raw_filtered.iter().enumerate() {
            let Some(label) = raw_type_label(&truth_record.record) else {
                continue;
            };
            by_type.entry(label).or_default().truth_total += 1;
            if let Some(query_index) = truth_matches[truth_index] {
                by_type.entry(label).or_default().tp += 1;
                let query_record = &query_raw_filtered[query_index];
                if !query_record.record.is_pass() {
                    filtered_by_type.entry(label).or_default().tp += 1;
                }
                if use_strelka_hcc_indel {
                    feature_rows.push(
                        0,
                        render_strelka_hcc_indel_tp_row(
                            0,
                            &truth_record.record,
                            &query_record.record,
                            &query_depths,
                        ),
                    )?;
                } else if use_generic_feature_table {
                    feature_rows.push(
                        0,
                        render_generic_tp_row(0, &truth_record.record, &query_record.record),
                    )?;
                } else if use_caller_feature_table {
                    caller_tp_truth.push(truth_record.record.clone());
                    caller_tp_query.push(query_record.record.clone());
                }
            } else {
                by_type.entry(label).or_default().fn_count += 1;
                if use_strelka_hcc_indel {
                    feature_rows
                        .push(2, render_strelka_hcc_indel_fn_row(0, &truth_record.record))?;
                } else if use_generic_feature_table {
                    feature_rows.push(2, render_generic_fn_row(0, &truth_record.record))?;
                } else if use_caller_feature_table {
                    caller_fn_truth.push(truth_record.record.clone());
                }
            }
        }

        for (query_index, query_record) in query_raw_filtered.iter().enumerate() {
            let type_label = raw_type_label(&query_record.record);
            if let Some(label) = type_label {
                by_type.entry(label).or_default().query_total += 1;
            }
            if query_matches[query_index].is_some() {
                continue;
            }
            // som.py classifies FP/AMBI intervals from POS through POS+len(REF),
            // even for symbolic and gVCF records. INFO/END and FORMAT/LEN are used
            // by bcftools -R preprocessing, not by this classification step.
            let class = classify_query(
                &query_record.record.chrom,
                query_record.record.pos,
                query_record.record.end_pos(),
                fp_regions.as_deref().unwrap_or(&[]),
                &ambiguous_regions,
                count_unk,
                ambi_fp,
            );
            if args.explain_ambiguous {
                record_ambiguous_explanation(
                    &query_record.record.chrom,
                    query_record.record.pos,
                    query_record.record.end_pos(),
                    &explanation_regions,
                    ambi_fp,
                    &mut ambiguous_classes,
                    &mut ambiguous_reasons,
                );
            }
            if let Some(label) = type_label {
                let row = by_type.entry(label).or_default();
                match class {
                    QueryClass::Fp => row.fp += 1,
                    QueryClass::Unk => row.unk += 1,
                    QueryClass::Ambi => row.ambi += 1,
                }
                if !query_record.record.is_pass() {
                    let filtered = filtered_by_type.entry(label).or_default();
                    match class {
                        QueryClass::Fp => filtered.fp += 1,
                        QueryClass::Unk => filtered.unk += 1,
                        QueryClass::Ambi => filtered.ambi += 1,
                    }
                }
            }
            let tag = match class {
                QueryClass::Fp => {
                    record_counts.fp += 1;
                    "FP"
                }
                QueryClass::Unk => {
                    record_counts.unk += 1;
                    "UNK"
                }
                QueryClass::Ambi => {
                    record_counts.ambi += 1;
                    "AMBI"
                }
            };
            if !query_record.record.is_pass() {
                match class {
                    QueryClass::Fp => filtered_records.fp += 1,
                    QueryClass::Unk => filtered_records.unk += 1,
                    QueryClass::Ambi => filtered_records.ambi += 1,
                }
            }
            if args.feature_table.is_some() {
                if use_strelka_hcc_indel {
                    let row = render_strelka_hcc_indel_query_row(
                        0,
                        &query_record.record,
                        tag,
                        &query_depths,
                    );
                    match class {
                        QueryClass::Fp => feature_rows.push(1, row)?,
                        QueryClass::Unk => feature_rows.push(4, row)?,
                        QueryClass::Ambi => feature_rows.push(3, row)?,
                    }
                } else if use_generic_feature_table {
                    let row = render_generic_query_row(0, &query_record.record, tag);
                    match class {
                        QueryClass::Fp => feature_rows.push(1, row)?,
                        QueryClass::Unk => feature_rows.push(4, row)?,
                        QueryClass::Ambi => feature_rows.push(3, row)?,
                    }
                } else if use_caller_feature_table {
                    match class {
                        QueryClass::Fp => caller_fp_query.push(query_record.record.clone()),
                        QueryClass::Unk => caller_unk_query.push(query_record.record.clone()),
                        QueryClass::Ambi => caller_ambi_query.push(query_record.record.clone()),
                    }
                }
            }
        }
        if use_caller_feature_table {
            let caller_table = build_caller_feature_table(
                feature_table_name,
                &truth_headers,
                &query_headers,
                (!bam_depths.is_empty()).then_some(&bam_depths),
                !args.no_order_check,
                &CallerRecordGroups {
                    tp_truth: &caller_tp_truth,
                    tp_query: &caller_tp_query,
                    fn_truth: &caller_fn_truth,
                    fp_query: &caller_fp_query,
                    ambi_query: &caller_ambi_query,
                    unk_query: &caller_unk_query,
                },
            )?;
            if let Some(existing) = feature_header.as_deref()
                && existing != caller_table.header
            {
                bail!("caller feature header changed between somatic contigs");
            }
            feature_header = Some(caller_table.header);
            feature_rows.extend(0, caller_table.tp)?;
            feature_rows.extend(1, caller_table.fp)?;
            feature_rows.extend(2, caller_table.fn_rows)?;
            feature_rows.extend(3, caller_table.ambi)?;
            feature_rows.extend(4, caller_table.unk)?;
        }
    }
    let ordered_feature_rows = feature_header
        .as_ref()
        .map(|_| feature_rows.renumber())
        .transpose()?;

    // som.py writes feature and ambiguity detail artifacts before it derives
    // the FP denominator. Preserve that order because the legacy range/FP
    // bug below exits after these files have already been created.
    if !ambiguous_classes.is_empty() {
        write_legacy_count_table(
            &suffixed_report_path(Path::new(&args.output), "ambiclasses.csv"),
            "class",
            &ambiguous_classes,
        )?;
    }
    if !ambiguous_reasons.is_empty() {
        write_legacy_count_table(
            &suffixed_report_path(Path::new(&args.output), "ambireasons.csv"),
            "reason",
            &ambiguous_reasons,
        )?;
    }
    if let Some(header) = feature_header.as_deref() {
        let ordered_rows = ordered_feature_rows
            .as_ref()
            .context("missing ordered somatic feature spool")?;
        let features_path = suffixed_report_path(Path::new(&args.output), "features.csv");
        let mut features = BufWriter::new(File::create(&features_path)?);
        writeln!(features, "{header}")?;
        std::io::copy(&mut File::open(ordered_rows.path())?, &mut features)?;
        features
            .flush()
            .with_context(|| format!("failed to write {}", features_path.display()))?;
        if let Some(roc_name) = args.roc.as_deref() {
            write_somatic_roc(
                &suffixed_report_path(Path::new(&args.output), "roc.csv"),
                header,
                ordered_rows.path(),
                roc_name,
            )?;
            if args.af_strat {
                for &(start, end) in &af_bins {
                    let Some(rows) = feature_rows_for_af_roc(
                        header,
                        ordered_rows.path(),
                        start,
                        end,
                        &args.af_strat_truth,
                        &args.af_strat_query,
                    )?
                    else {
                        continue;
                    };
                    for prefix in ["records", "SNVs", "indels"] {
                        let path = PathBuf::from(format!(
                            "{}.{}.{}.roc.csv",
                            args.output,
                            prefix,
                            format_af_interval(start, end)
                        ));
                        write_somatic_roc(&path, header, rows.path(), roc_name)?;
                    }
                }
            }
        }
    }

    if fp_region_size_requires_reference(
        args.fp_region_size.as_deref(),
        fp_regions.as_deref().unwrap_or(&[]),
        &ambiguous_regions,
    ) && reference_lengths.is_empty()
    {
        reference_lengths = fasta::contig_lengths(Path::new(&args.reference))?;
    }
    validate_legacy_fp_location_denominator(
        args.fp_region_size.as_deref(),
        args.location.as_deref(),
        has_automatic_fp_bases(fp_regions.as_deref().unwrap_or(&[]), &ambiguous_regions),
    )?;
    let fp_region_size = calculate_fp_region_size_for_contigs(
        args.fp_region_size.as_deref(),
        fp_regions.as_deref().unwrap_or(&[]),
        &ambiguous_regions,
        locations.as_deref(),
        &reference_lengths,
        &truth_contigs,
    );

    let commandline = somatic_commandline(&args);

    let mut lines = Vec::new();
    let use_af_column_order = args.af_strat && !af_bins.is_empty();
    lines.push(if use_af_column_order {
        stats_header_af(args.count_filtered_fn)
    } else {
        stats_header(args.count_filtered_fn)
    });
    // These indexes are the stable pandas indexes produced by the pinned
    // Python-2 reference after its sequence of bcftools-stat merges. They are
    // serialized by DataFrame.to_csv and are therefore part of the contract.
    if args.af_strat {
        for (index, label) in [
            (0, "indels"),
            (1, "SNVs"),
            (2, "no-ALTs"),
            (5, "records"),
            (6, "MNPs"),
            (7, "others"),
        ] {
            let row = if label == "records" {
                record_counts
            } else {
                by_type.get(label).copied().unwrap_or_default()
            };
            let filtered = if label == "records" {
                args.count_filtered_fn.then_some(filtered_records)
            } else {
                filtered_counts_for_type(
                    args.count_filtered_fn,
                    args.feature_table.as_deref(),
                    label,
                    &filtered_by_type,
                )
            };
            let context = StatsRowContext {
                fp_region_size,
                ci_alpha,
                filtered,
                include_filtered_columns: args.count_filtered_fn,
                commandline: &commandline,
            };
            lines.push(if use_af_column_order {
                render_row_af(index, label, row, &context)
            } else if args.af_strat {
                render_row_af_without_bins(index, label, row, &context)
            } else {
                render_row(index, label, row, &context)
            });
        }
    } else {
        for (index, label) in STATS_TYPE_ROWS {
            let Some(row) = by_type.get(label).filter(|row| row.truth_total > 0) else {
                continue;
            };
            lines.push(render_row(
                index,
                label,
                *row,
                &StatsRowContext {
                    fp_region_size,
                    ci_alpha,
                    filtered: filtered_counts_for_type(
                        args.count_filtered_fn,
                        args.feature_table.as_deref(),
                        label,
                        &filtered_by_type,
                    ),
                    include_filtered_columns: args.count_filtered_fn,
                    commandline: &commandline,
                },
            ));
        }
        if record_counts.truth_total > 0 {
            lines.push(render_row(
                5,
                "records",
                record_counts,
                &StatsRowContext {
                    fp_region_size,
                    ci_alpha,
                    filtered: args.count_filtered_fn.then_some(filtered_records),
                    include_filtered_columns: args.count_filtered_fn,
                    commandline: &commandline,
                },
            ));
        }
    }
    if args.af_strat
        && let Some(header) = feature_header.as_deref()
    {
        for prefix in ["records", "SNVs", "indels"] {
            let af_counts = calculate_af_stats(
                header,
                ordered_feature_rows
                    .as_ref()
                    .context("missing ordered somatic feature spool")?
                    .path(),
                &args.af_strat_binsize,
                &args.af_strat_truth,
                &args.af_strat_query,
            )?;
            for (start, end, counts, filtered) in &af_counts {
                let label = format!("{prefix}.{}", format_af_interval(*start, *end));
                lines.push(render_row_af(
                    0,
                    &label,
                    *counts,
                    &StatsRowContext {
                        fp_region_size,
                        ci_alpha,
                        filtered: args.count_filtered_fn.then_some(*filtered),
                        include_filtered_columns: args.count_filtered_fn,
                        commandline: &commandline,
                    },
                ));
            }
        }
    }

    let output = suffixed_report_path(Path::new(&args.output), "stats.csv");
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }
    fs::write(&output, format!("{}\n", lines.join("\n")))
        .with_context(|| format!("failed to write {}", output.display()))?;
    controls.print_default_summary(&lines);

    let metrics_json = suffixed_report_path(Path::new(&args.output), "metrics.json");
    let explanation_enabled = args.explain_ambiguous && !args.ambiguous_beds.is_empty();
    write_legacy_metrics_json(
        &metrics_json,
        &commandline,
        &output,
        ci_alpha,
        explanation_enabled.then_some(&ambiguous_classes),
        explanation_enabled.then_some(&ambiguous_reasons),
    )?;

    if let Some(header) = feature_header.as_deref() {
        let ordered_rows = ordered_feature_rows
            .as_ref()
            .context("missing ordered somatic feature spool")?;
        if args.happy_stats {
            write_happy_style_summary(
                &suffixed_report_path(Path::new(&args.output), "summary.csv"),
                header,
                ordered_rows.path(),
                feature_table_name,
            )?;
            if args.af_strat {
                write_happy_style_extended(
                    &suffixed_report_path(Path::new(&args.output), "extended.csv"),
                    header,
                    ordered_rows.path(),
                    feature_table_name,
                    &args.af_strat_binsize,
                    &args.af_strat_truth,
                    &args.af_strat_query,
                )?;
            }
        }
    }
    controls.info("Somatic comparison complete")?;
    controls.scratch.cleanup()
}
