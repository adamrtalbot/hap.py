use crate::cli::SomaticArgs;
use crate::compare::suffixed_report_path;
use crate::{fasta, ftx, strelka, vcf};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

/// Version string emitted in the `sompyversion` column and JSON metadata.
/// Matches legacy som.py when `Haplo.version.__version__` is empty — the
/// pinned oracle container resolves this to `"som.py-"` (trailing hyphen with
/// nothing after). We reproduce the exact literal so stats.csv byte-matches.
const SOM_VERSION: &str = "som.py-";
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
    key: vcf::VariantKey,
    record: vcf::RawVcfRecord,
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
    interval: vcf::BedInterval,
    label: String,
    details: Vec<String>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum QueryClass {
    Fp,
    Unk,
    Ambi,
}

pub fn run(mut args: SomaticArgs) -> Result<()> {
    if let Some(config) = args.roc.as_deref().and_then(somatic_roc_config) {
        args.feature_table = Some(config.feature_table.to_string());
    }
    validate_args(&args)?;
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
    let mut reference_sequences = if normalize_truth || normalize_query {
        Some(fasta::read_sequences(Path::new(&args.reference))?)
    } else {
        None
    };
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
    let (truth_headers, mut truth_raw) = vcf::load_raw_vcf(if reuse_truth {
        &truth_cache
    } else {
        Path::new(&args.truth)
    })?;
    let (query_headers, mut query_raw) = vcf::load_raw_vcf(if reuse_query {
        &query_cache
    } else {
        Path::new(&args.query)
    })?;
    let bam_depths = ftx::bam_normalization_depths(&args.bams)?;
    if normalize_truth {
        truth_raw = normalize_somatic_records(
            truth_raw,
            reference_sequences
                .as_ref()
                .expect("normalization loaded the reference"),
        );
    }
    if normalize_query {
        query_raw = normalize_somatic_records(
            query_raw,
            reference_sequences
                .as_ref()
                .expect("normalization loaded the reference"),
        );
    }
    if !reuse_truth {
        vcf::write_raw_vcf(&truth_cache, &truth_headers, &truth_raw)?;
    }
    if !reuse_query {
        vcf::write_raw_vcf(&query_cache, &query_headers, &query_raw)?;
    }
    let truth_raw_filtered = filter_raw_records(
        truth_raw,
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
    let query_raw_filtered = filter_raw_records(
        query_raw,
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
    let (truth_matches, query_matches) =
        pair_exact_records(&truth_raw_filtered, &query_raw_filtered);
    let query_depths = if bam_depths.is_empty() {
        strelka::parse_depths(&query_headers)
    } else {
        bam_depths.clone()
    };

    let mut by_type: BTreeMap<&'static str, SomaticCounts> = BTreeMap::new();
    let mut filtered_by_type: BTreeMap<&'static str, FilteredCounts> = BTreeMap::new();
    let mut filtered_records = FilteredCounts::default();
    let mut record_counts = SomaticCounts {
        truth_total: truth_raw_filtered.len(),
        query_total: query_raw_filtered.len(),
        ..SomaticCounts::default()
    };
    record_counts.tp = truth_matches.iter().filter(|entry| entry.is_some()).count();
    record_counts.fn_count = truth_raw_filtered.len().saturating_sub(record_counts.tp);
    filtered_records.tp = query_raw_filtered
        .iter()
        .zip(&query_matches)
        .filter(|(record, matched)| matched.is_some() && !record.record.is_pass())
        .count();
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
    let mut tp_rows: Vec<String> = Vec::new();
    let mut fp_rows: Vec<String> = Vec::new();
    let mut fn_rows: Vec<String> = Vec::new();
    let mut ambi_rows: Vec<String> = Vec::new();
    let mut unk_rows: Vec<String> = Vec::new();
    let mut caller_tp_truth = Vec::new();
    let mut caller_tp_query = Vec::new();
    let mut caller_fn_truth = Vec::new();
    let mut caller_fp_query = Vec::new();
    let mut caller_ambi_query = Vec::new();
    let mut caller_unk_query = Vec::new();
    let mut ambiguous_classes = BTreeMap::new();
    let mut ambiguous_reasons = BTreeMap::new();
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
                tp_rows.push(render_strelka_hcc_indel_tp_row(
                    0,
                    &truth_record.record,
                    &query_record.record,
                    &query_depths,
                ));
            } else if use_generic_feature_table {
                tp_rows.push(render_generic_tp_row(
                    0,
                    &truth_record.record,
                    &query_record.record,
                ));
            } else if use_caller_feature_table {
                caller_tp_truth.push(truth_record.record.clone());
                caller_tp_query.push(query_record.record.clone());
            }
        } else {
            by_type.entry(label).or_default().fn_count += 1;
            if use_strelka_hcc_indel {
                fn_rows.push(render_strelka_hcc_indel_fn_row(0, &truth_record.record));
            } else if use_generic_feature_table {
                fn_rows.push(render_generic_fn_row(0, &truth_record.record));
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
            &query_record.key.chrom,
            query_record.key.pos,
            query_record.record.end_pos(),
            fp_regions.as_deref().unwrap_or(&[]),
            &ambiguous_regions,
            count_unk,
            ambi_fp,
        );
        if args.explain_ambiguous {
            record_ambiguous_explanation(
                &query_record.key.chrom,
                query_record.key.pos,
                query_record.record.end_pos(),
                &ambiguous_regions,
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
                let row =
                    render_strelka_hcc_indel_query_row(0, &query_record.record, tag, &query_depths);
                match class {
                    QueryClass::Fp => fp_rows.push(row),
                    QueryClass::Unk => unk_rows.push(row),
                    QueryClass::Ambi => ambi_rows.push(row),
                }
            } else if use_generic_feature_table {
                let row = render_generic_query_row(0, &query_record.record, tag);
                match class {
                    QueryClass::Fp => fp_rows.push(row),
                    QueryClass::Unk => unk_rows.push(row),
                    QueryClass::Ambi => ambi_rows.push(row),
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
        feature_header = Some(caller_table.header);
        tp_rows = caller_table.tp;
        fp_rows = caller_table.fp;
        fn_rows = caller_table.fn_rows;
        ambi_rows = caller_table.ambi;
        unk_rows = caller_table.unk;
    }
    let ordered_feature_rows = feature_header
        .as_ref()
        .map(|_| renumber_feature_rows(&[tp_rows, fp_rows, fn_rows, ambi_rows, unk_rows]));

    if fp_region_size_requires_reference(
        args.fp_region_size.as_deref(),
        fp_regions.as_deref().unwrap_or(&[]),
    ) && reference_sequences.is_none()
    {
        reference_sequences = Some(fasta::read_sequences(Path::new(&args.reference))?);
    }
    let empty_reference = BTreeMap::new();
    let fp_region_size = calculate_fp_region_size(
        args.fp_region_size.as_deref(),
        fp_regions.as_deref().unwrap_or(&[]),
        locations.as_deref(),
        reference_sequences.as_ref().unwrap_or(&empty_reference),
        &truth_raw_filtered,
        &query_raw_filtered,
    );

    let commandline = legacy_som_commandline(&args);

    let mut lines = Vec::new();
    lines.push(if args.af_strat {
        stats_header_af(args.count_filtered_fn)
    } else {
        stats_header(args.count_filtered_fn)
    });
    // These indexes are the stable pandas indexes produced by the pinned
    // Python-2 oracle after its sequence of bcftools-stat merges. They are
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
            lines.push(render_row_af(
                index,
                label,
                row,
                &StatsRowContext {
                    fp_region_size,
                    ci_alpha,
                    filtered,
                    include_filtered_columns: args.count_filtered_fn,
                    commandline: &commandline,
                },
            ));
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
        let af_counts = calculate_af_stats(
            header,
            ordered_feature_rows.as_deref().unwrap_or_default(),
            &args.af_strat_binsize,
            &args.af_strat_truth,
            &args.af_strat_query,
        )?;
        for prefix in ["records", "SNVs", "indels"] {
            for (start, end, counts, filtered) in &af_counts {
                let label = format!("{prefix}.{start:.6}-{end:.6}");
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
        explanation_enabled.then_some(&ambiguous_classes),
        explanation_enabled.then_some(&ambiguous_reasons),
    )?;

    if let Some(header) = feature_header {
        let features = suffixed_report_path(Path::new(&args.output), "features.csv");
        let ordered_rows = ordered_feature_rows.as_deref().unwrap_or_default();
        fs::write(
            &features,
            format!("{header}\n{}\n", ordered_rows.join("\n")),
        )
        .with_context(|| format!("failed to write {}", features.display()))?;
        if let Some(roc_name) = args.roc.as_deref() {
            write_somatic_roc(
                &suffixed_report_path(Path::new(&args.output), "roc.csv"),
                &header,
                ordered_rows,
                roc_name,
            )?;
            if args.af_strat {
                for (start, end) in parse_af_bins(&args.af_strat_binsize) {
                    let rows = feature_rows_for_af_roc(
                        &header,
                        ordered_rows,
                        start,
                        end,
                        &args.af_strat_truth,
                        &args.af_strat_query,
                    )?;
                    if rows.is_empty() {
                        continue;
                    }
                    for prefix in ["records", "SNVs", "indels"] {
                        let path = PathBuf::from(format!(
                            "{}.{}.{start:.6}-{end:.6}.roc.csv",
                            args.output, prefix
                        ));
                        write_somatic_roc(&path, &header, &rows, roc_name)?;
                    }
                }
            }
        }
        if args.happy_stats {
            write_happy_style_summary(
                &suffixed_report_path(Path::new(&args.output), "summary.csv"),
                &header,
                ordered_rows,
                feature_table_name,
            )?;
            if args.af_strat {
                write_happy_style_extended(
                    &suffixed_report_path(Path::new(&args.output), "extended.csv"),
                    &header,
                    ordered_rows,
                    feature_table_name,
                    &args.af_strat_binsize,
                    &args.af_strat_truth,
                    &args.af_strat_query,
                )?;
            }
        }
    }
    if args.explain_ambiguous && !args.ambiguous_beds.is_empty() {
        write_legacy_count_table(
            &suffixed_report_path(Path::new(&args.output), "ambiclasses.csv"),
            "class",
            &ambiguous_classes,
        )?;
        write_legacy_count_table(
            &suffixed_report_path(Path::new(&args.output), "ambireasons.csv"),
            "reason",
            &ambiguous_reasons,
        )?;
    }
    controls.info("Somatic comparison complete")?;
    controls.scratch.cleanup()
}

fn renumber_feature_rows(groups: &[Vec<String>]) -> Vec<String> {
    let mut out = Vec::new();
    for group in groups {
        for (idx, row) in group.iter().enumerate() {
            let tail = row.split_once(',').map(|(_, tail)| tail).unwrap_or(row);
            out.push(format!("{idx},{tail}"));
        }
    }
    out
}

struct CallerFeatureTable {
    header: String,
    tp: Vec<String>,
    fp: Vec<String>,
    fn_rows: Vec<String>,
    ambi: Vec<String>,
    unk: Vec<String>,
}

struct CallerRecordGroups<'a> {
    tp_truth: &'a [vcf::RawVcfRecord],
    tp_query: &'a [vcf::RawVcfRecord],
    fn_truth: &'a [vcf::RawVcfRecord],
    fp_query: &'a [vcf::RawVcfRecord],
    ambi_query: &'a [vcf::RawVcfRecord],
    unk_query: &'a [vcf::RawVcfRecord],
}

struct ParsedFeatureTable {
    columns: BTreeSet<String>,
    rows: Vec<BTreeMap<String, String>>,
}

impl ParsedFeatureTable {
    fn from_lines(lines: Vec<String>) -> Result<Self> {
        let mut lines = lines.into_iter();
        let header = parse_csv_line(&lines.next().context("feature table has no header")?);
        let named_columns = header.iter().skip(1).cloned().collect::<Vec<_>>();
        let columns = named_columns.iter().cloned().collect();
        let rows = lines
            .map(|line| {
                let cells = parse_csv_line(&line);
                named_columns
                    .iter()
                    .cloned()
                    .zip(cells.into_iter().skip(1))
                    .collect()
            })
            .collect();
        Ok(Self { columns, rows })
    }
}

fn build_caller_feature_table(
    feature: &str,
    truth_headers: &[String],
    query_headers: &[String],
    depths: Option<&BTreeMap<String, f64>>,
    check_order: bool,
    groups: &CallerRecordGroups<'_>,
) -> Result<CallerFeatureTable> {
    let truth_tp = ParsedFeatureTable::from_lines(ftx::emit_feature_table_with_depths(
        feature,
        groups.tp_truth,
        truth_headers,
        "TP",
        depths,
    )?)?;
    let query_tp = ParsedFeatureTable::from_lines(ftx::emit_feature_table_with_depths(
        feature,
        groups.tp_query,
        query_headers,
        "TP_r",
        depths,
    )?)?;
    if truth_tp.rows.len() != query_tp.rows.len() {
        bail!(
            "cannot merge TP features: truth and query lengths differ ({} != {})",
            truth_tp.rows.len(),
            query_tp.rows.len()
        );
    }
    if check_order {
        for (index, (truth, query)) in truth_tp.rows.iter().zip(&query_tp.rows).enumerate() {
            for column in ["CHROM", "POS"] {
                if truth.get(column) != query.get(column) {
                    bail!(
                        "cannot merge TP features: inputs are out of order at row {index} ({column})"
                    );
                }
            }
        }
    }

    let (tp_columns, tp) = merge_caller_tp_tables(&truth_tp, &query_tp);
    let truth_fn = ParsedFeatureTable::from_lines(ftx::emit_feature_table_with_depths(
        feature,
        groups.fn_truth,
        truth_headers,
        "FN",
        depths,
    )?)?;
    let query_fp = ParsedFeatureTable::from_lines(ftx::emit_feature_table_with_depths(
        feature,
        groups.fp_query,
        query_headers,
        "FP",
        depths,
    )?)?;
    let query_ambi = ParsedFeatureTable::from_lines(ftx::emit_feature_table_with_depths(
        feature,
        groups.ambi_query,
        query_headers,
        "AMBI",
        depths,
    )?)?;
    let query_unk = ParsedFeatureTable::from_lines(ftx::emit_feature_table_with_depths(
        feature,
        groups.unk_query,
        query_headers,
        "UNK",
        depths,
    )?)?;
    let (fn_columns, fn_rows) = suffix_truth_columns(truth_fn, &tp_columns);

    let mut all_columns = tp_columns;
    all_columns.extend(fn_columns);
    all_columns.extend(query_fp.columns.iter().cloned());
    all_columns.extend(query_ambi.columns.iter().cloned());
    all_columns.extend(query_unk.columns.iter().cloned());
    let columns = ordered_somatic_feature_columns(&all_columns);
    let header = format!(",{}", columns.join(","));

    Ok(CallerFeatureTable {
        header,
        tp: render_feature_group(&tp, &columns),
        fp: render_feature_group(&query_fp.rows, &columns),
        fn_rows: render_feature_group(&fn_rows, &columns),
        ambi: render_feature_group(&query_ambi.rows, &columns),
        unk: render_feature_group(&query_unk.rows, &columns),
    })
}

fn merge_caller_tp_tables(
    truth: &ParsedFeatureTable,
    query: &ParsedFeatureTable,
) -> (BTreeSet<String>, Vec<BTreeMap<String, String>>) {
    let shared_keys = ["CHROM", "POS", "tag"];
    let union = truth
        .columns
        .union(&query.columns)
        .cloned()
        .collect::<BTreeSet<_>>();
    let mut columns = BTreeSet::new();
    for column in &union {
        if shared_keys.contains(&column.as_str()) {
            columns.insert(column.clone());
        } else if truth.columns.contains(column) && query.columns.contains(column) {
            columns.insert(column.clone());
            columns.insert(format!("{column}.truth"));
        } else {
            columns.insert(column.clone());
        }
    }

    let rows = truth
        .rows
        .iter()
        .zip(&query.rows)
        .map(|(truth_row, query_row)| {
            let mut row = BTreeMap::new();
            for column in &union {
                if shared_keys.contains(&column.as_str()) {
                    if let Some(value) = truth_row.get(column).or_else(|| query_row.get(column)) {
                        row.insert(column.clone(), value.clone());
                    }
                } else if let (Some(truth_value), Some(query_value)) =
                    (truth_row.get(column), query_row.get(column))
                {
                    row.insert(column.clone(), query_value.clone());
                    row.insert(format!("{column}.truth"), truth_value.clone());
                } else if let Some(value) = query_row.get(column).or_else(|| truth_row.get(column))
                {
                    row.insert(column.clone(), value.clone());
                }
            }
            row
        })
        .collect();
    (columns, rows)
}

fn suffix_truth_columns(
    table: ParsedFeatureTable,
    tp_columns: &BTreeSet<String>,
) -> (BTreeSet<String>, Vec<BTreeMap<String, String>>) {
    let rename = |column: &str| {
        let truth_name = format!("{column}.truth");
        if tp_columns.contains(&truth_name) {
            truth_name
        } else {
            column.to_string()
        }
    };
    let columns = table.columns.iter().map(|column| rename(column)).collect();
    let rows = table
        .rows
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|(column, value)| (rename(&column), value))
                .collect()
        })
        .collect();
    (columns, rows)
}

fn ordered_somatic_feature_columns(columns: &BTreeSet<String>) -> Vec<String> {
    let first = [
        "CHROM",
        "POS",
        "tag",
        "REF",
        "REF.truth",
        "ALT",
        "ALT.truth",
    ];
    first
        .iter()
        .filter(|column| columns.contains(**column))
        .map(|column| (*column).to_string())
        .chain(
            columns
                .iter()
                .filter(|column| !first.contains(&column.as_str()))
                .cloned(),
        )
        .collect()
}

fn render_feature_group(rows: &[BTreeMap<String, String>], columns: &[String]) -> Vec<String> {
    rows.iter()
        .enumerate()
        .map(|(index, row)| {
            csv_join(
                std::iter::once(index.to_string()).chain(columns.iter().map(|column| {
                    format_somatic_feature_cell(column, row.get(column).map_or("", String::as_str))
                })),
            )
        })
        .collect()
}

fn format_somatic_feature_cell(column: &str, value: &str) -> String {
    const TEXT_COLUMNS: &[&str] = &[
        "CHROM",
        "tag",
        "REF",
        "REF.truth",
        "ALT",
        "ALT.truth",
        "FILTER",
        "NT",
        "SGT",
        "S.2.GT",
    ];
    if value.is_empty() || column == "POS" || TEXT_COLUMNS.contains(&column) {
        value.to_string()
    } else if let Ok(number) = value.parse::<f64>() {
        format!("{number:.8}")
    } else {
        value.to_string()
    }
}

fn render_generic_tp_row(
    index: usize,
    truth: &vcf::RawVcfRecord,
    query: &vcf::RawVcfRecord,
) -> String {
    csv_join([
        index.to_string(),
        query.chrom.clone(),
        query.pos.to_string(),
        "TP".to_string(),
        query.ref_allele.clone(),
        truth.ref_allele.clone(),
        query.alt_allele.clone(),
        truth.alt_allele.clone(),
        normalized_filter(&query.filter),
        normalized_filter(&truth.filter),
        format_feature_number_or_text(&normalized_qual(&query.qual)),
        format_feature_number_or_text(&normalized_qual(&truth.qual)),
    ])
}

fn render_generic_fn_row(index: usize, truth: &vcf::RawVcfRecord) -> String {
    csv_join([
        index.to_string(),
        truth.chrom.clone(),
        truth.pos.to_string(),
        "FN".to_string(),
        String::new(),
        truth.ref_allele.clone(),
        String::new(),
        truth.alt_allele.clone(),
        String::new(),
        normalized_filter(&truth.filter),
        String::new(),
        format_feature_number_or_text(&normalized_qual(&truth.qual)),
    ])
}

fn render_generic_query_row(index: usize, query: &vcf::RawVcfRecord, tag: &str) -> String {
    csv_join([
        index.to_string(),
        query.chrom.clone(),
        query.pos.to_string(),
        tag.to_string(),
        query.ref_allele.clone(),
        String::new(),
        query.alt_allele.clone(),
        String::new(),
        normalized_filter(&query.filter),
        String::new(),
        format_feature_number_or_text(&normalized_qual(&query.qual)),
        String::new(),
    ])
}

const STRELKA_HCC_INDEL_HEADER: &str = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,EVS,FILTER,I.DP_normal,I.DP_tumor,I.T_ALT_RATE,I.count,I.tag,IC,IHP,INDELTYPE,LENGTH,MQ,MQ0,NT,NT_REF,N_AF,N_BCN,N_DP,N_DP_RATE,N_FDP,QSI_NT,QUAL,RC,RU,RU_LEN,S.1.VT,SGT,T_AF,T_BCN,T_DP,T_DP_RATE,T_FDP";

fn render_strelka_hcc_indel_tp_row(
    index: usize,
    truth: &vcf::RawVcfRecord,
    query: &vcf::RawVcfRecord,
    avg_depth: &BTreeMap<String, f64>,
) -> String {
    let query_row = strelka_indel_row(query, avg_depth);
    let truth_row = generic_truth_row(truth);
    merge_feature_rows(index, "TP", &query_row, &truth_row)
}

fn render_strelka_hcc_indel_fn_row(index: usize, truth: &vcf::RawVcfRecord) -> String {
    let query_row = blank_strelka_query_row(truth);
    let truth_row = generic_truth_row(truth);
    merge_feature_rows(index, "FN", &query_row, &truth_row)
}

fn render_strelka_hcc_indel_query_row(
    index: usize,
    query: &vcf::RawVcfRecord,
    tag: &str,
    avg_depth: &BTreeMap<String, f64>,
) -> String {
    let query_row = strelka_indel_row(query, avg_depth);
    merge_feature_rows(index, tag, &query_row, &BTreeMap::new())
}

fn merge_feature_rows(
    index: usize,
    tag: &str,
    query_row: &BTreeMap<String, String>,
    truth_row: &BTreeMap<String, String>,
) -> String {
    let cols = vec![
        index.to_string(),
        query_row.get("CHROM").cloned().unwrap_or_default(),
        query_row.get("POS").cloned().unwrap_or_default(),
        tag.to_string(),
        query_row.get("REF").cloned().unwrap_or_default(),
        truth_row.get("REF").cloned().unwrap_or_default(),
        query_row.get("ALT").cloned().unwrap_or_default(),
        truth_row.get("ALT").cloned().unwrap_or_default(),
        query_row.get("EVS").cloned().unwrap_or_default(),
        query_row.get("FILTER").cloned().unwrap_or_default(),
        truth_row.get("I.DP_normal").cloned().unwrap_or_default(),
        truth_row.get("I.DP_tumor").cloned().unwrap_or_default(),
        truth_row.get("I.T_ALT_RATE").cloned().unwrap_or_default(),
        truth_row.get("I.count").cloned().unwrap_or_default(),
        truth_row.get("I.tag").cloned().unwrap_or_default(),
        query_row.get("IC").cloned().unwrap_or_default(),
        query_row.get("IHP").cloned().unwrap_or_default(),
        query_row.get("INDELTYPE").cloned().unwrap_or_default(),
        query_row.get("LENGTH").cloned().unwrap_or_default(),
        query_row.get("MQ").cloned().unwrap_or_default(),
        query_row.get("MQ0").cloned().unwrap_or_default(),
        query_row.get("NT").cloned().unwrap_or_default(),
        query_row.get("NT_REF").cloned().unwrap_or_default(),
        query_row.get("N_AF").cloned().unwrap_or_default(),
        query_row.get("N_BCN").cloned().unwrap_or_default(),
        query_row.get("N_DP").cloned().unwrap_or_default(),
        query_row.get("N_DP_RATE").cloned().unwrap_or_default(),
        query_row.get("N_FDP").cloned().unwrap_or_default(),
        query_row.get("QSI_NT").cloned().unwrap_or_default(),
        truth_row.get("QUAL").cloned().unwrap_or_default(),
        query_row.get("RC").cloned().unwrap_or_default(),
        query_row.get("RU").cloned().unwrap_or_default(),
        query_row.get("RU_LEN").cloned().unwrap_or_default(),
        truth_row.get("S.1.VT").cloned().unwrap_or_default(),
        query_row.get("SGT").cloned().unwrap_or_default(),
        query_row.get("T_AF").cloned().unwrap_or_default(),
        query_row.get("T_BCN").cloned().unwrap_or_default(),
        query_row.get("T_DP").cloned().unwrap_or_default(),
        query_row.get("T_DP_RATE").cloned().unwrap_or_default(),
        query_row.get("T_FDP").cloned().unwrap_or_default(),
    ];
    csv_join(cols)
}

fn csv_join<I>(cols: I) -> String
where
    I: IntoIterator<Item = String>,
{
    cols.into_iter()
        .map(|value| csv_escape(&value))
        .collect::<Vec<_>>()
        .join(",")
}

fn csv_escape(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn blank_strelka_query_row(truth: &vcf::RawVcfRecord) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("CHROM".to_string(), truth.chrom.clone());
    row.insert("POS".to_string(), truth.pos.to_string());
    row.insert("REF".to_string(), String::new());
    row.insert("ALT".to_string(), String::new());
    for key in [
        "EVS",
        "FILTER",
        "IC",
        "IHP",
        "INDELTYPE",
        "LENGTH",
        "MQ",
        "MQ0",
        "NT",
        "NT_REF",
        "N_AF",
        "N_BCN",
        "N_DP",
        "N_DP_RATE",
        "N_FDP",
        "QSI_NT",
        "RC",
        "RU",
        "RU_LEN",
        "SGT",
        "T_AF",
        "T_BCN",
        "T_DP",
        "T_DP_RATE",
        "T_FDP",
    ] {
        row.insert(key.to_string(), String::new());
    }
    row
}

fn generic_truth_row(record: &vcf::RawVcfRecord) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("REF".to_string(), record.ref_allele.clone());
    row.insert("ALT".to_string(), record.alt_allele.clone());
    row.insert(
        "QUAL".to_string(),
        format_feature_number_or_text(&normalized_qual(&record.qual)),
    );
    let sample0 = record.sample_map(0);
    for key in [
        "I.T_ALT_RATE",
        "I.DP_normal",
        "I.DP_tumor",
        "I.tag",
        "I.count",
    ] {
        row.insert(
            key.to_string(),
            format_feature_number_or_text(
                &strelka::info_value(&record.info, key.trim_start_matches("I."))
                    .unwrap_or_default(),
            ),
        );
    }
    row.insert(
        "S.1.VT".to_string(),
        sample0.get("VT").cloned().unwrap_or_default(),
    );
    row
}

fn format_feature_number_or_text(value: &str) -> String {
    value
        .parse::<f64>()
        .map(|number| format!("{number:.8}"))
        .unwrap_or_else(|_| value.to_string())
}

fn strelka_indel_row(
    record: &vcf::RawVcfRecord,
    avg_depth: &BTreeMap<String, f64>,
) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    let n = record.sample_map(0);
    let t = record.sample_map(1);
    let ref_len = record.ref_allele.len();
    let alt_lens: Vec<usize> = record.alt_allele.split(',').map(str::len).collect();
    let max_len = alt_lens
        .iter()
        .copied()
        .max()
        .unwrap_or(ref_len)
        .max(ref_len);
    let min_len = alt_lens
        .iter()
        .copied()
        .min()
        .unwrap_or(ref_len)
        .min(ref_len);
    let mut indel_type = 0;
    for alt in record.alt_allele.split(',') {
        if alt.len() > ref_len {
            indel_type |= 1;
        } else {
            indel_type |= 2;
        }
    }
    let nt = strelka::info_value(&record.info, "NT").unwrap_or_default();
    let n_dp = strelka::parse_first_number(n.get("DP")).unwrap_or(0.0);
    let t_dp = strelka::parse_first_number(t.get("DP")).unwrap_or(0.0);
    let norm = avg_depth.get(&record.chrom).copied().unwrap_or(0.0);

    row.insert("CHROM".to_string(), record.chrom.clone());
    row.insert("POS".to_string(), record.pos.to_string());
    row.insert("REF".to_string(), record.ref_allele.clone());
    row.insert("ALT".to_string(), record.alt_allele.clone());
    row.insert(
        "EVS".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "SomaticEVS")
                .or_else(|| strelka::info_float(&record.info, "EVS"))
                .unwrap_or(-1.0)
        ),
    );
    row.insert("FILTER".to_string(), normalized_filter(&record.filter));
    row.insert(
        "IC".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "IC").unwrap_or(0.0)
        ),
    );
    row.insert(
        "IHP".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "IHP").unwrap_or(0.0)
        ),
    );
    row.insert("INDELTYPE".to_string(), format!("{:.8}", indel_type as f64));
    row.insert(
        "LENGTH".to_string(),
        format!("{:.8}", (max_len - min_len) as f64),
    );
    row.insert(
        "MQ".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "MQ").unwrap_or(0.0)
        ),
    );
    row.insert(
        "MQ0".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "MQ0").unwrap_or(0.0)
        ),
    );
    row.insert("NT".to_string(), nt.clone());
    row.insert(
        "NT_REF".to_string(),
        format!("{:.8}", if nt == "ref" { 1.0 } else { 0.0 }),
    );
    row.insert(
        "N_AF".to_string(),
        format!("{:.8}", strelka::af_from_tir_tar(&n, "TIR", "TAR")),
    );
    row.insert(
        "N_BCN".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(n.get("BCN50")).unwrap_or(0.0)
        ),
    );
    row.insert("N_DP".to_string(), format!("{:.8}", n_dp));
    row.insert(
        "N_DP_RATE".to_string(),
        format!("{:.8}", if norm > 0.0 { n_dp / norm } else { 0.0 }),
    );
    row.insert(
        "N_FDP".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(n.get("FDP50")).unwrap_or(0.0)
        ),
    );
    row.insert(
        "QSI_NT".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "QSI_NT").unwrap_or(0.0)
        ),
    );
    row.insert(
        "RC".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "RC").unwrap_or(0.0)
        ),
    );
    row.insert(
        "RU".to_string(),
        strelka::info_value(&record.info, "RU").unwrap_or_default(),
    );
    row.insert(
        "RU_LEN".to_string(),
        format!(
            "{:.8}",
            strelka::info_value(&record.info, "RU")
                .map(|v| v.len() as f64)
                .unwrap_or(0.0)
        ),
    );
    row.insert(
        "SGT".to_string(),
        strelka::info_value(&record.info, "SGT").unwrap_or_default(),
    );
    row.insert(
        "T_AF".to_string(),
        format!("{:.8}", strelka::af_from_tir_tar(&t, "TIR", "TAR")),
    );
    row.insert(
        "T_BCN".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(t.get("BCN50")).unwrap_or(0.0)
        ),
    );
    row.insert("T_DP".to_string(), format!("{:.8}", t_dp));
    row.insert(
        "T_DP_RATE".to_string(),
        format!("{:.8}", if norm > 0.0 { t_dp / norm } else { 0.0 }),
    );
    row.insert(
        "T_FDP".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(t.get("FDP50")).unwrap_or(0.0)
        ),
    );
    row
}

fn normalized_filter(filter: &str) -> String {
    if filter == "PASS" || filter == "." {
        String::new()
    } else {
        filter.to_string()
    }
}

fn normalized_qual(qual: &str) -> String {
    if qual == "." {
        String::new()
    } else {
        qual.to_string()
    }
}

fn write_simple_table(path: &Path, body: &str) -> Result<()> {
    fs::write(path, body).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn write_happy_style_summary(
    path: &Path,
    feature_header: &str,
    feature_rows: &[String],
    feature_table: &str,
) -> Result<()> {
    let mut lines = vec![
        ",Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio".to_string()
    ];
    let headers = parse_csv_line(feature_header);
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = csv_column_index(&headers, "FILTER")?;
    let ref_index = csv_column_index(&headers, "REF")?;
    let alt_index = csv_column_index(&headers, "ALT")?;
    let truth_ref_index = csv_column_index(&headers, "REF.truth")?;
    let rows = feature_rows
        .iter()
        .map(|row| parse_csv_line(row))
        .collect::<Vec<_>>();
    let happy_type = match feature_table.rsplit('.').next() {
        Some("snv") => "SNP",
        Some("indel") => "INDEL",
        _ => "NA",
    };
    let truth_tags = rows
        .iter()
        .filter(|row| nonempty_csv_field(row, truth_ref_index))
        .filter_map(|row| row.get(tag_index).map(String::as_str))
        .collect::<Vec<_>>();
    let truth_total = truth_tags.len();
    let truth_tp = truth_tags.iter().filter(|tag| **tag == "TP").count();
    let truth_fn = truth_tags.iter().filter(|tag| **tag == "FN").count();

    for filter in ["PASS", "ALL"] {
        let query_tags = rows
            .iter()
            .filter(|row| {
                nonempty_csv_field(row, ref_index)
                    && nonempty_csv_field(row, alt_index)
                    && (filter == "ALL"
                        || row.get(filter_index).is_none_or(|value| {
                            value.is_empty() || value == "." || value == "PASS"
                        }))
            })
            .filter_map(|row| row.get(tag_index).map(String::as_str))
            .collect::<Vec<_>>();
        let query_fp = query_tags.iter().filter(|tag| **tag == "FP").count();
        let query_unk = query_tags
            .iter()
            .filter(|tag| **tag == "UNK" || **tag == "AMBI")
            .count();
        let query_total = truth_tp + query_fp + query_unk;
        let recall = rounded_metric(truth_tp, truth_total);
        let precision = rounded_metric(truth_tp, truth_tp + query_fp);
        let frac_na = rounded_metric(query_unk, query_total);
        let f1 = match (recall, precision) {
            (Some(recall), Some(precision)) if recall + precision > 0.0 => {
                Some(round_four(2.0 * recall * precision / (recall + precision)))
            }
            _ => None,
        };
        lines.push(csv_join([
            "0".to_string(),
            happy_type.to_string(),
            filter.to_string(),
            truth_total.to_string(),
            truth_tp.to_string(),
            truth_fn.to_string(),
            query_total.to_string(),
            query_fp.to_string(),
            query_unk.to_string(),
            "NA".to_string(),
            render_summary_metric(recall),
            render_summary_metric(precision),
            render_summary_metric(frac_na),
            render_summary_metric(f1),
            "NA".to_string(),
            "NA".to_string(),
            "NA".to_string(),
            "NA".to_string(),
        ]));
    }
    write_simple_table(path, &format!("{}\n", lines.join("\n")))
}

fn render_summary_metric(value: Option<f64>) -> String {
    value.map_or_else(String::new, py_float)
}

fn write_happy_style_extended(
    path: &Path,
    feature_header: &str,
    feature_rows: &[String],
    feature_table: &str,
    bin_sizes: &str,
    truth_af_field: &str,
    query_af_field: &str,
) -> Result<()> {
    let mut lines = vec![
        ",Type,Subtype,Subset,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio".to_string()
    ];
    let headers = parse_csv_line(feature_header);
    let truth_af_index = headers
        .iter()
        .position(|header| header == truth_af_field)
        .with_context(|| {
            format!("truth AF feature '{truth_af_field}' is not in the feature table")
        })?;
    let query_af_index = headers
        .iter()
        .position(|header| header == query_af_field)
        .with_context(|| {
            format!("query AF feature '{query_af_field}' is not in the feature table")
        })?;
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = csv_column_index(&headers, "FILTER")?;
    let ref_index = csv_column_index(&headers, "REF")?;
    let alt_index = csv_column_index(&headers, "ALT")?;
    let truth_ref_index = csv_column_index(&headers, "REF.truth")?;
    let rows = feature_rows
        .iter()
        .map(|row| parse_csv_line(row))
        .collect::<Vec<_>>();
    let happy_type = match feature_table.rsplit('.').next() {
        Some("snv") => "SNP",
        Some("indel") => "INDEL",
        _ => "NA",
    };

    let mut row_index = 0usize;
    for (start, end) in parse_af_bins(bin_sizes) {
        let inclusive_last = end >= 1.0;
        let subset = if inclusive_last {
            format!("[{start:.2},1.00]")
        } else {
            format!("[{start:.2},{end:.2})")
        };
        let in_bin = |value: &str| {
            value.parse::<f64>().is_ok_and(|value| {
                value >= start && (value < end || (inclusive_last && value <= 1.0))
            })
        };
        let truth_rows = rows.iter().filter(|row| {
            nonempty_csv_field(row, truth_ref_index)
                && row.get(truth_af_index).is_some_and(|value| in_bin(value))
        });
        let truth_tags = truth_rows
            .filter_map(|row| row.get(tag_index).map(String::as_str))
            .collect::<Vec<_>>();
        let truth_total = truth_tags.len();
        let truth_tp = truth_tags.iter().filter(|tag| **tag == "TP").count();
        let truth_fn = truth_tags.iter().filter(|tag| **tag == "FN").count();

        for filter in ["PASS", "ALL"] {
            let query_rows = rows.iter().filter(|row| {
                nonempty_csv_field(row, ref_index)
                    && nonempty_csv_field(row, alt_index)
                    && row.get(query_af_index).is_some_and(|value| in_bin(value))
                    && (filter == "ALL"
                        || row.get(filter_index).is_none_or(|value| {
                            value.is_empty() || value == "." || value == "PASS"
                        }))
            });
            let query_tags = query_rows
                .filter_map(|row| row.get(tag_index).map(String::as_str))
                .collect::<Vec<_>>();
            let query_fp = query_tags.iter().filter(|tag| **tag == "FP").count();
            let query_unk = query_tags
                .iter()
                .filter(|tag| **tag == "UNK" || **tag == "AMBI")
                .count();
            let query_total = truth_tp + query_fp + query_unk;
            let recall = rounded_metric(truth_tp, truth_total);
            let precision = rounded_metric(truth_tp, truth_tp + query_fp);
            let frac_na = rounded_metric(query_unk, query_total);
            let f1 = match (recall, precision) {
                (Some(recall), Some(precision)) if recall + precision > 0.0 => {
                    Some(round_four(2.0 * recall * precision / (recall + precision)))
                }
                _ => None,
            };
            lines.push(csv_join([
                row_index.to_string(),
                happy_type.to_string(),
                "*".to_string(),
                subset.clone(),
                filter.to_string(),
                truth_total.to_string(),
                truth_tp.to_string(),
                truth_fn.to_string(),
                query_total.to_string(),
                query_fp.to_string(),
                query_unk.to_string(),
                "NA".to_string(),
                render_summary_metric(recall),
                render_summary_metric(precision),
                render_summary_metric(frac_na),
                render_summary_metric(f1),
                "NA".to_string(),
                "NA".to_string(),
                "NA".to_string(),
                "NA".to_string(),
            ]));
            row_index += 1;
        }
    }
    write_simple_table(path, &format!("{}\n", lines.join("\n")))
}

fn csv_column_index(headers: &[String], name: &str) -> Result<usize> {
    headers
        .iter()
        .position(|header| header == name)
        .with_context(|| format!("feature table is missing required column '{name}'"))
}

fn nonempty_csv_field(row: &[String], index: usize) -> bool {
    row.get(index)
        .is_some_and(|value| !value.is_empty() && value != ".")
}

fn parse_csv_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => fields.push(std::mem::take(&mut field)),
            _ => field.push(ch),
        }
    }
    fields.push(field);
    fields
}

fn rounded_metric(numerator: usize, denominator: usize) -> Option<f64> {
    (denominator != 0).then(|| round_four(numerator as f64 / denominator as f64))
}

fn round_four(value: f64) -> f64 {
    (value * 10_000.0).round() / 10_000.0
}

fn parse_af_bins(raw: &str) -> Vec<(f64, f64)> {
    let mut bins = raw
        .split(',')
        .filter_map(|part| part.parse::<f64>().ok())
        .collect::<Vec<_>>();
    if bins.is_empty() {
        bins.push(0.2);
    }
    let mut out = Vec::new();
    let mut start = 0.0;
    let mut idx = 0usize;
    while start < 1.0 {
        let mut end = start + bins[idx];
        if end >= 1.0 {
            end = 1.0;
        }
        out.push((start, end));
        if end >= 1.0 {
            break;
        }
        start = end;
        idx = (idx + 1) % bins.len();
    }
    out
}

fn calculate_af_stats(
    feature_header: &str,
    feature_rows: &[String],
    bin_sizes: &str,
    truth_af_field: &str,
    query_af_field: &str,
) -> Result<Vec<(f64, f64, SomaticCounts, FilteredCounts)>> {
    let headers = parse_csv_line(feature_header);
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = csv_column_index(&headers, "FILTER")?;
    let truth_af_index = csv_column_index(&headers, truth_af_field)?;
    let query_af_index = csv_column_index(&headers, query_af_field)?;
    let rows = feature_rows
        .iter()
        .map(|row| parse_csv_line(row))
        .collect::<Vec<_>>();

    let mut output = Vec::new();
    for (start, end) in parse_af_bins(bin_sizes) {
        let in_bin = |row: &[String], index: usize| {
            row.get(index)
                .and_then(|value| value.parse::<f64>().ok())
                .is_some_and(|value| value >= start && value <= 1.0 && value < end)
                || (end >= 1.0
                    && row
                        .get(index)
                        .and_then(|value| value.parse::<f64>().ok())
                        .is_some_and(|value| value == 1.0))
        };

        let mut counts = SomaticCounts::default();
        let mut filtered = FilteredCounts::default();
        for row in &rows {
            let tag = row.get(tag_index).map(String::as_str).unwrap_or_default();
            let filtered_call = row.get(filter_index).is_some_and(|value| !value.is_empty());
            match tag {
                "TP" if in_bin(row, truth_af_index) => {
                    counts.tp += 1;
                    if filtered_call {
                        filtered.tp += 1;
                    }
                }
                "FN" if in_bin(row, truth_af_index) => counts.fn_count += 1,
                "FP" if in_bin(row, query_af_index) => {
                    counts.fp += 1;
                    if filtered_call {
                        filtered.fp += 1;
                    }
                }
                "UNK" if in_bin(row, query_af_index) => {
                    counts.unk += 1;
                    if filtered_call {
                        filtered.unk += 1;
                    }
                }
                "AMBI" if in_bin(row, query_af_index) => {
                    counts.ambi += 1;
                    if filtered_call {
                        filtered.ambi += 1;
                    }
                }
                _ => {}
            }
        }
        counts.truth_total = counts.tp + counts.fn_count;
        counts.query_total = counts.tp + counts.fp + counts.unk + counts.ambi;
        output.push((start, end, counts, filtered));
    }
    Ok(output)
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum SomaticRocTag {
    Tp,
    Fp,
    Fn,
}

fn write_somatic_roc(
    path: &Path,
    feature_header: &str,
    feature_rows: &[String],
    roc_name: &str,
) -> Result<()> {
    let config = somatic_roc_config(roc_name)
        .with_context(|| format!("unsupported somatic ROC mode '{roc_name}'"))?;
    let headers = parse_csv_line(feature_header);
    let score_index = csv_column_index(&headers, config.score)?;
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = headers.iter().position(|field| field == "FILTER");
    let nt_index = headers.iter().position(|field| field == "NT");
    let mut observations = Vec::new();

    for line in feature_rows {
        let fields = parse_csv_line(line);
        let Some(tag) = fields.get(tag_index).and_then(|value| {
            let lower = value.to_ascii_lowercase();
            if lower.starts_with("tp") {
                Some(SomaticRocTag::Tp)
            } else if lower.starts_with("fp") {
                Some(SomaticRocTag::Fp)
            } else if lower.starts_with("fn") {
                Some(SomaticRocTag::Fn)
            } else {
                None
            }
        }) else {
            continue;
        };
        let mut score = fields
            .get(score_index)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
        if config.zero_score_unless_nt_ref
            && nt_index
                .and_then(|index| fields.get(index))
                .is_some_and(|value| value != "ref")
        {
            score = 0.0;
        }
        if let (Some(index), Some(filter_name)) = (filter_index, config.filter_name) {
            let mut filters = fields
                .get(index)
                .map(String::as_str)
                .unwrap_or_default()
                .split([';', ','])
                .filter(|value| !value.is_empty() && *value != "." && *value != "PASS")
                .collect::<Vec<_>>();
            filters.retain(|value| *value != filter_name);
            if !filters.is_empty() {
                score = f64::MIN_POSITIVE;
            }
        }
        observations.push((score, tag));
    }
    observations.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.cmp(&right.1))
    });

    let mut tp = observations
        .iter()
        .filter(|(_, tag)| *tag == SomaticRocTag::Tp)
        .count();
    let mut fp = observations
        .iter()
        .filter(|(_, tag)| *tag == SomaticRocTag::Fp)
        .count();
    let mut fn_count = observations
        .iter()
        .filter(|(_, tag)| *tag == SomaticRocTag::Fn)
        .count();
    let mut roc_rows = Vec::new();
    let mut previous = None;
    for (score, tag) in observations {
        if previous != Some(score) {
            let precision = if tp + fp == 0 {
                1.0
            } else {
                tp as f64 / (tp + fp) as f64
            };
            let recall = if tp + fn_count == 0 {
                0.0
            } else {
                tp as f64 / (tp + fn_count) as f64
            };
            roc_rows.push((
                cpp_default_six(score),
                tp,
                fp,
                fn_count,
                cpp_default_six(precision),
                cpp_default_six(recall),
            ));
            previous = Some(score);
        }
        match tag {
            SomaticRocTag::Tp => {
                tp = tp.saturating_sub(1);
                fn_count += 1;
            }
            SomaticRocTag::Fp => fp = fp.saturating_sub(1),
            SomaticRocTag::Fn => {}
        }
    }
    let integer_scores = roc_rows
        .iter()
        .all(|(score, ..)| score.parse::<i64>().is_ok());
    let mut lines = vec![format!(",{},tp,fp,fn,precision,recall", config.score)];
    for (row_index, (score, tp, fp, fn_count, precision, recall)) in
        roc_rows.into_iter().enumerate()
    {
        let score = if integer_scores {
            score
        } else {
            score
                .parse::<f64>()
                .map(|value| format!("{value:.8}"))
                .unwrap_or(score)
        };
        let precision = precision
            .parse::<f64>()
            .map(|value| format!("{value:.8}"))
            .unwrap_or(precision);
        let recall = recall
            .parse::<f64>()
            .map(|value| format!("{value:.8}"))
            .unwrap_or(recall);
        lines.push(format!(
            "{row_index},{score},{tp},{fp},{fn_count},{precision},{recall}"
        ));
    }
    write_simple_table(path, &format!("{}\n", lines.join("\n")))
}

fn cpp_default_six(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    let exponent = value.abs().log10().floor() as i32;
    if !(-4..6).contains(&exponent) {
        let scientific = format!("{value:.5e}");
        let Some((mantissa, exponent)) = scientific.split_once('e') else {
            return scientific;
        };
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        return format!("{mantissa}e{exponent}");
    }
    let decimals = (5 - exponent).max(0) as usize;
    let fixed = format!("{value:.decimals$}");
    if fixed.contains('.') {
        fixed
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    } else {
        fixed
    }
}

fn feature_rows_for_af_roc(
    feature_header: &str,
    feature_rows: &[String],
    start: f64,
    end: f64,
    truth_af_field: &str,
    query_af_field: &str,
) -> Result<Vec<String>> {
    let headers = parse_csv_line(feature_header);
    let tag_index = csv_column_index(&headers, "tag")?;
    let truth_index = csv_column_index(&headers, truth_af_field)?;
    let query_index = csv_column_index(&headers, query_af_field)?;
    let in_bin = |fields: &[String], index: usize| {
        fields
            .get(index)
            .and_then(|value| value.parse::<f64>().ok())
            .is_some_and(|value| value >= start && (value < end || (end >= 1.0 && value == 1.0)))
    };
    Ok(feature_rows
        .iter()
        .filter(|line| {
            let fields = parse_csv_line(line);
            match fields.get(tag_index).map(String::as_str) {
                Some("TP" | "FN") => in_bin(&fields, truth_index),
                Some("FP") => in_bin(&fields, query_index),
                _ => false,
            }
        })
        .cloned()
        .collect())
}

fn raw_type_label(record: &vcf::RawVcfRecord) -> Option<&'static str> {
    let alternates = record.alt_allele.split(',').collect::<Vec<_>>();
    if alternates.iter().all(|alternate| *alternate == ".") {
        return None;
    }
    if alternates.iter().any(|alternate| {
        alternate == &"*" || alternate.starts_with('<') || alternate.contains(['[', ']'])
    }) {
        return Some("others");
    }

    let reference_len = record.ref_allele.len();
    if reference_len == 1 && alternates.iter().all(|alternate| alternate.len() == 1) {
        Some("SNVs")
    } else if reference_len > 1
        && alternates
            .iter()
            .all(|alternate| alternate.len() == reference_len)
    {
        Some("MNPs")
    } else {
        Some("indels")
    }
}

fn contigs_in_play(truth: &[FilteredRawRecord], query: &[FilteredRawRecord]) -> BTreeSet<String> {
    truth
        .iter()
        .chain(query.iter())
        .map(|record| record.key.chrom.clone())
        .collect()
}

fn fp_region_size_requires_reference(
    requested: Option<&str>,
    fp_regions: &[vcf::BedInterval],
) -> bool {
    requested
        .and_then(|value| value.parse::<usize>().ok())
        .is_none()
        && fp_regions.is_empty()
}

fn calculate_fp_region_size(
    requested: Option<&str>,
    fp_regions: &[vcf::BedInterval],
    locations: Option<&[vcf::LocationFilter]>,
    reference_sequences: &BTreeMap<String, String>,
    truth: &[FilteredRawRecord],
    query: &[FilteredRawRecord],
) -> usize {
    if let Some(size) = requested.and_then(|value| value.parse::<usize>().ok()) {
        return size;
    }

    if !fp_regions.is_empty() {
        return fp_regions
            .iter()
            .map(|interval| match locations {
                None => interval.end.saturating_sub(interval.start),
                Some(locations) => locations
                    .iter()
                    .map(|location| match location {
                        vcf::LocationFilter::Contig(chrom) if chrom == &interval.chrom => {
                            interval.end.saturating_sub(interval.start)
                        }
                        vcf::LocationFilter::Range { chrom, start, end }
                            if chrom == &interval.chrom =>
                        {
                            let location_start = start.saturating_sub(1);
                            interval
                                .end
                                .min(*end)
                                .saturating_sub(interval.start.max(location_start))
                        }
                        _ => 0,
                    })
                    .sum(),
            })
            .sum();
    }

    if let Some(locations) = locations {
        return locations
            .iter()
            .map(|location| match location {
                vcf::LocationFilter::Contig(chrom) => reference_sequences
                    .get(chrom)
                    .map_or(0, |sequence| sequence.len()),
                vcf::LocationFilter::Range { chrom, start, end } => {
                    reference_sequences.get(chrom).map_or(0, |sequence| {
                        end.min(&sequence.len())
                            .saturating_sub(start.saturating_sub(1))
                    })
                }
            })
            .sum();
    }

    contigs_in_play(truth, query)
        .into_iter()
        .filter_map(|contig| {
            reference_sequences
                .get(&contig)
                .map(|sequence| sequence.len())
        })
        .sum()
}

fn pair_exact_records(
    truth: &[FilteredRawRecord],
    query: &[FilteredRawRecord],
) -> (Vec<Option<usize>>, Vec<Option<usize>>) {
    let mut query_by_key: BTreeMap<vcf::VariantKey, VecDeque<usize>> = BTreeMap::new();
    for (index, record) in query.iter().enumerate() {
        query_by_key
            .entry(record.key.clone())
            .or_default()
            .push_back(index);
    }

    let mut truth_matches = vec![None; truth.len()];
    let mut query_matches = vec![None; query.len()];
    for (truth_index, truth_record) in truth.iter().enumerate() {
        let Some(query_index) = query_by_key
            .get_mut(&truth_record.key)
            .and_then(VecDeque::pop_front)
        else {
            continue;
        };
        truth_matches[truth_index] = Some(query_index);
        query_matches[query_index] = Some(truth_index);
    }
    (truth_matches, query_matches)
}

fn render_row(
    index: usize,
    label: &str,
    counts: SomaticCounts,
    context: &StatsRowContext<'_>,
) -> String {
    let (recall, recall_lower, recall_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fn_count, context.ci_alpha);
    let (precision, precision_lower, precision_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fp, context.ci_alpha);
    let fp_rate = if context.fp_region_size == 0 {
        0.0
    } else {
        1_000_000.0 * counts.fp as f64 / context.fp_region_size as f64
    };
    let version = SOM_VERSION;
    let mut columns = vec![
        index.to_string(),
        label.to_string(),
        counts.truth_total.to_string(),
        counts.query_total.to_string(),
        counts.tp.to_string(),
        counts.fp.to_string(),
        counts.fn_count.to_string(),
        counts.unk.to_string(),
        counts.ambi.to_string(),
    ];
    if context.include_filtered_columns {
        if let Some(filtered) = context.filtered {
            columns.extend([
                py_float(filtered.fp as f64),
                py_float(filtered.tp as f64),
                py_float(filtered.unk as f64),
                py_float(filtered.ambi as f64),
            ]);
        } else {
            columns.extend([String::new(), String::new(), String::new(), String::new()]);
        }
    }
    columns.extend([
        py_float(recall),
        py_float(recall_lower),
        py_float(recall_upper),
        py_float(ratio(counts.tp, counts.truth_total)),
        py_float(precision),
        py_float(precision_lower),
        py_float(precision_upper),
        py_float(ratio(counts.unk, counts.query_total)),
        py_float(ratio(counts.ambi, counts.query_total)),
        context.fp_region_size.to_string(),
        py_float(fp_rate),
    ]);
    if context.include_filtered_columns {
        if let Some(filtered) = context.filtered {
            let unfiltered_tp = counts.tp.saturating_sub(filtered.tp);
            let unfiltered_fp = counts.fp.saturating_sub(filtered.fp);
            columns.extend([
                py_float(ratio(unfiltered_tp, counts.tp + counts.fn_count)),
                py_float(ratio(unfiltered_tp, unfiltered_tp + unfiltered_fp)),
                py_float(if context.fp_region_size == 0 {
                    0.0
                } else {
                    1_000_000.0 * unfiltered_fp as f64 / context.fp_region_size as f64
                }),
                py_float(ratio(
                    counts.unk.saturating_sub(filtered.unk),
                    counts.query_total,
                )),
                py_float(ratio(
                    counts.ambi.saturating_sub(filtered.ambi),
                    counts.query_total,
                )),
            ]);
        } else {
            columns.extend([
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
            ]);
        }
    }
    columns.push(version.to_string());
    columns.push(context.commandline.to_string());
    columns.join(",")
}

/// Render the column order produced by the legacy pandas concat used when
/// allele-frequency stratification is enabled.  Adding the bin rows causes
/// pandas 0.x to alphabetize the count columns while leaving subsequently
/// calculated metrics in assignment order; that accidental ordering is part
/// of the byte-for-byte CSV contract.
fn render_row_af(
    index: usize,
    label: &str,
    counts: SomaticCounts,
    context: &StatsRowContext<'_>,
) -> String {
    let (recall, recall_lower, recall_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fn_count, context.ci_alpha);
    let (precision, precision_lower, precision_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fp, context.ci_alpha);
    let fp_rate = if context.fp_region_size == 0 {
        0.0
    } else {
        1_000_000.0 * counts.fp as f64 / context.fp_region_size as f64
    };
    let mut columns = vec![index.to_string(), counts.ambi.to_string()];
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.ambi as f64)),
        );
    }
    columns.extend([counts.fn_count.to_string(), counts.fp.to_string()]);
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.fp as f64)),
        );
    }
    columns.extend([
        counts.query_total.to_string(),
        counts.truth_total.to_string(),
        counts.tp.to_string(),
    ]);
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.tp as f64)),
        );
    }
    columns.extend([label.to_string(), counts.unk.to_string()]);
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.unk as f64)),
        );
    }
    columns.extend([
        py_float(recall),
        py_float(recall_lower),
        py_float(recall_upper),
        ratio_or_blank(counts.tp, counts.truth_total),
        py_float(precision),
        py_float(precision_lower),
        py_float(precision_upper),
        ratio_or_blank(counts.unk, counts.query_total),
        ratio_or_blank(counts.ambi, counts.query_total),
        context.fp_region_size.to_string(),
        py_float(fp_rate),
    ]);
    if context.include_filtered_columns {
        if let Some(filtered) = context.filtered {
            let unfiltered_tp = counts.tp.saturating_sub(filtered.tp);
            let unfiltered_fp = counts.fp.saturating_sub(filtered.fp);
            columns.extend([
                ratio_or_blank(unfiltered_tp, counts.tp + counts.fn_count),
                ratio_or_blank(unfiltered_tp, unfiltered_tp + unfiltered_fp),
                py_float(if context.fp_region_size == 0 {
                    0.0
                } else {
                    1_000_000.0 * unfiltered_fp as f64 / context.fp_region_size as f64
                }),
                ratio_or_blank(counts.unk.saturating_sub(filtered.unk), counts.query_total),
                ratio_or_blank(
                    counts.ambi.saturating_sub(filtered.ambi),
                    counts.query_total,
                ),
            ]);
        } else {
            columns.extend(std::iter::repeat_n(String::new(), 5));
        }
    }
    columns.push(SOM_VERSION.to_string());
    columns.push(context.commandline.to_string());
    columns.join(",")
}

fn ratio_or_blank(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        String::new()
    } else {
        py_float(numerator as f64 / denominator as f64)
    }
}

fn stats_header(count_filtered_fn: bool) -> String {
    let mut columns = vec![
        "",
        "type",
        "total.truth",
        "total.query",
        "tp",
        "fp",
        "fn",
        "unk",
        "ambi",
    ];
    if count_filtered_fn {
        columns.extend([
            "fp.filtered",
            "tp.filtered",
            "unk.filtered",
            "ambi.filtered",
        ]);
    }
    columns.extend([
        "recall",
        "recall_lower",
        "recall_upper",
        "recall2",
        "precision",
        "precision_lower",
        "precision_upper",
        "na",
        "ambiguous",
        "fp.region.size",
        "fp.rate",
    ]);
    if count_filtered_fn {
        columns.extend([
            "recall.filtered",
            "precision.filtered",
            "fp.rate.filtered",
            "na.filtered",
            "ambiguous.filtered",
        ]);
    }
    columns.extend(["sompyversion", "sompycmd"]);
    columns.join(",")
}

fn stats_header_af(count_filtered_fn: bool) -> String {
    let mut columns = vec![""];
    columns.push("ambi");
    if count_filtered_fn {
        columns.push("ambi.filtered");
    }
    columns.extend(["fn", "fp"]);
    if count_filtered_fn {
        columns.push("fp.filtered");
    }
    columns.extend(["total.query", "total.truth", "tp"]);
    if count_filtered_fn {
        columns.push("tp.filtered");
    }
    columns.extend(["type", "unk"]);
    if count_filtered_fn {
        columns.push("unk.filtered");
    }
    columns.extend([
        "recall",
        "recall_lower",
        "recall_upper",
        "recall2",
        "precision",
        "precision_lower",
        "precision_upper",
        "na",
        "ambiguous",
        "fp.region.size",
        "fp.rate",
    ]);
    if count_filtered_fn {
        columns.extend([
            "recall.filtered",
            "precision.filtered",
            "fp.rate.filtered",
            "na.filtered",
            "ambiguous.filtered",
        ]);
    }
    columns.extend(["sompyversion", "sompycmd"]);
    columns.join(",")
}

fn filtered_counts_for_type(
    enabled: bool,
    feature_table: Option<&str>,
    label: &str,
    filtered_by_type: &BTreeMap<&'static str, FilteredCounts>,
) -> Option<FilteredCounts> {
    if !enabled {
        return None;
    }
    if feature_table == Some("generic") {
        return Some(filtered_by_type.get(label).copied().unwrap_or_default());
    }
    let selected_type = match feature_table.and_then(|name| name.rsplit('.').next()) {
        Some("snv") => "SNVs",
        Some("indel") => "indels",
        _ => return None,
    };
    if label == selected_type {
        Some(
            filtered_by_type
                .get(selected_type)
                .copied()
                .unwrap_or_default(),
        )
    } else {
        None
    }
}

fn validate_args(args: &SomaticArgs) -> Result<()> {
    if args.af_strat && args.feature_table.is_none() {
        bail!("--bin-afs requires --feature-table");
    }
    if !(0.0 < args.ci_level && args.ci_level < 1.0) {
        bail!("confidence interval level must be > 0.0 and < 1.0");
    }
    validate_af_bins(&args.af_strat_binsize)?;
    if args.af_strat && (args.af_strat_truth.is_empty() || args.af_strat_query.is_empty()) {
        bail!("AF truth and query feature names must not be empty");
    }
    if args.count_filtered_fn && (!args.include_nonpass || args.feature_table.is_none()) {
        bail!("--count-filtered-fn requires -P/--include-nonpass and --feature-table");
    }
    if args.happy_stats && (!args.include_nonpass || args.feature_table.is_none()) {
        bail!("--happy-stats requires -P/--include-nonpass and --feature-table");
    }
    Ok(())
}

fn resolve_toggle(enabled: bool, explicitly_disabled: bool) -> bool {
    enabled && !explicitly_disabled
}

fn selected_normalizations(args: &SomaticArgs) -> (bool, bool) {
    (
        args.normalize_truth || args.normalize_all,
        args.normalize_query || args.normalize_all,
    )
}

fn validate_af_bins(raw: &str) -> Result<()> {
    let bins = raw.split(',').collect::<Vec<_>>();
    if bins.is_empty() || bins.iter().any(|bin| bin.trim().is_empty()) {
        bail!("AF bin size list must not be empty");
    }
    for bin in bins {
        let value = bin
            .parse::<f64>()
            .with_context(|| format!("failed to parse AF bin size '{bin}'"))?;
        if !value.is_finite() || value <= 0.0 {
            bail!("AF bin sizes must be finite and greater than zero");
        }
    }
    Ok(())
}

fn classification_bed_chrom(
    chrom: &str,
    reference_contigs: &BTreeSet<String>,
    fixchr_truth: bool,
) -> String {
    if fixchr_truth {
        somatic_chrom(chrom, reference_contigs, true)
    } else {
        chrom.to_string()
    }
}

fn load_classification_bed(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    fixchr_truth: bool,
) -> Result<Vec<vcf::BedInterval>> {
    let mut intervals = vcf::load_bed(path, &BTreeSet::new())?;
    for interval in &mut intervals {
        interval.chrom = classification_bed_chrom(&interval.chrom, reference_contigs, fixchr_truth);
    }
    Ok(intervals)
}

fn load_ambiguous_beds(
    paths: &[String],
    reference_contigs: &BTreeSet<String>,
    fixchr_truth: bool,
) -> Result<Vec<AmbiguousInterval>> {
    let mut intervals = Vec::new();
    for path in paths {
        let path = Path::new(path);
        let text = vcf::read_text(path)
            .with_context(|| format!("failed to read ambiguous BED {}", path.display()))?;
        for (line_index, line) in text.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let fields = line.split(['\t', ',', ';']).collect::<Vec<_>>();
            if fields.len() < 5 {
                bail!(
                    "ambiguous BED line {} has fewer than 5 columns in {}",
                    line_index + 1,
                    path.display()
                );
            }
            let start = fields[1].parse::<usize>().with_context(|| {
                format!(
                    "invalid BED start '{}' on line {} in {}",
                    fields[1],
                    line_index + 1,
                    path.display()
                )
            })?;
            let end = fields[2].parse::<usize>().with_context(|| {
                format!(
                    "invalid BED end '{}' on line {} in {}",
                    fields[2],
                    line_index + 1,
                    path.display()
                )
            })?;
            intervals.push(AmbiguousInterval {
                interval: vcf::BedInterval {
                    chrom: classification_bed_chrom(fields[0], reference_contigs, fixchr_truth),
                    start,
                    end,
                },
                // som.py's implementation selects xe[4], i.e. the fifth BED
                // field despite the historical help text calling it field 4.
                label: fields[4].to_string(),
                // BedIntervalTree stores the derived label followed by every
                // original BED field after CHROM/start/end. som.py then uses
                // value[1] for replicate count and value[3..] for reasons.
                details: fields[3..]
                    .iter()
                    .map(|value| (*value).to_string())
                    .collect(),
            });
        }
    }
    Ok(intervals)
}

fn record_ambiguous_explanation(
    chrom: &str,
    start: usize,
    end: usize,
    ambiguous_regions: &[AmbiguousInterval],
    ambi_fp: bool,
    classes: &mut BTreeMap<String, usize>,
    reasons: &mut BTreeMap<String, usize>,
) {
    let mut classes_this_position = BTreeSet::new();
    for entry in ambiguous_regions
        .iter()
        .filter(|entry| entry.interval.overlaps(chrom, start, end))
    {
        let reason = match entry.label.as_str() {
            "fp" if ambi_fp => "FP",
            "fp" => "ambi-fp",
            "unk" => "ambi-unk",
            label => label,
        };
        classes_this_position.insert(reason.to_string());
        let replicate = entry.details.first().map(String::as_str).unwrap_or("*");
        *reasons
            .entry(format!("{reason}: rep. count {replicate}"))
            .or_default() += 1;
        for detail in entry.details.iter().skip(2) {
            *reasons.entry(format!("{reason}: {detail}")).or_default() += 1;
        }
    }
    for class in classes_this_position {
        *classes.entry(class).or_default() += 1;
    }
}

fn write_legacy_count_table(
    path: &Path,
    label: &str,
    counts: &BTreeMap<String, usize>,
) -> Result<()> {
    if counts.is_empty() {
        return Ok(());
    }
    let legacy_indices = python2_counter_indices(counts);
    let mut lines = vec![format!(",{label},count")];
    for (value, count) in counts {
        lines.push(csv_join([
            legacy_indices[value].to_string(),
            value.to_string(),
            count.to_string(),
        ]));
    }
    write_simple_table(path, &format!("{}\n", lines.join("\n")))
}

fn classify_query(
    chrom: &str,
    start: usize,
    end: usize,
    fp_regions: &[vcf::BedInterval],
    ambiguous_regions: &[AmbiguousInterval],
    count_unk: bool,
    ambi_fp: bool,
) -> QueryClass {
    if fp_regions
        .iter()
        .any(|interval| interval.overlaps(chrom, start, end))
    {
        return QueryClass::Fp;
    }

    let overlapping = ambiguous_regions
        .iter()
        .filter(|entry| entry.interval.overlaps(chrom, start, end))
        .collect::<Vec<_>>();
    if overlapping
        .iter()
        .any(|entry| entry.label == "FP" || (ambi_fp && entry.label == "fp"))
    {
        return QueryClass::Fp;
    }
    if !overlapping.is_empty() {
        return QueryClass::Ambi;
    }
    if count_unk {
        QueryClass::Unk
    } else {
        QueryClass::Fp
    }
}

fn ratio(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

fn py_float(value: f64) -> String {
    let text = value.to_string();
    let clean = if text == "-0" { "0" } else { text.as_str() };
    if clean.contains(['.', 'e', 'E']) {
        clean.to_string()
    } else {
        format!("{clean}.0")
    }
}

fn normalize_somatic_records(
    records: Vec<vcf::RawVcfRecord>,
    reference_sequences: &BTreeMap<String, String>,
) -> Vec<vcf::RawVcfRecord> {
    let reference_contigs = reference_sequences.keys().cloned().collect::<BTreeSet<_>>();
    let mut seen = BTreeSet::new();
    let mut output = Vec::new();
    for mut record in records {
        let lookup_chrom = vcf::normalize_chrom(&record.chrom, &reference_contigs);
        let Some(reference) = reference_sequences.get(&lookup_chrom) else {
            continue;
        };
        let reference = reference.as_bytes();
        let start = record.pos.saturating_sub(1);
        let end = start.saturating_add(record.ref_allele.len());
        if end > reference.len()
            || !reference[start..end].eq_ignore_ascii_case(record.ref_allele.as_bytes())
        {
            // bcftools norm `-c x` excludes reference-mismatch records.
            continue;
        }
        if !record.alt_allele.contains(',')
            && !record.alt_allele.starts_with('<')
            && !record.alt_allele.contains(['[', ']'])
            && record.alt_allele != "*"
            && record.alt_allele != "."
        {
            normalize_somatic_alleles(&mut record, reference);
        }
        let key = (
            record.chrom.clone(),
            record.pos,
            record.ref_allele.clone(),
            record.alt_allele.clone(),
        );
        if seen.insert(key) {
            output.push(record);
        }
    }
    output
}

fn normalize_somatic_alleles(record: &mut vcf::RawVcfRecord, reference: &[u8]) {
    let mut reference_allele = record.ref_allele.as_bytes().to_vec();
    let mut alternate = record.alt_allele.as_bytes().to_vec();
    while reference_allele.len() > 1
        && alternate.len() > 1
        && reference_allele.last().map(u8::to_ascii_uppercase)
            == alternate.last().map(u8::to_ascii_uppercase)
    {
        reference_allele.pop();
        alternate.pop();
    }
    while reference_allele.len() > 1
        && alternate.len() > 1
        && reference_allele.first().map(u8::to_ascii_uppercase)
            == alternate.first().map(u8::to_ascii_uppercase)
    {
        reference_allele.remove(0);
        alternate.remove(0);
        record.pos += 1;
    }
    while reference_allele.len() != alternate.len() && record.pos > 1 {
        let previous = reference[record.pos - 2].to_ascii_uppercase();
        let longer_last = if reference_allele.len() > alternate.len() {
            reference_allele.last()
        } else {
            alternate.last()
        }
        .copied()
        .map(|base| base.to_ascii_uppercase());
        if longer_last != Some(previous) {
            break;
        }
        reference_allele.pop();
        alternate.pop();
        reference_allele.insert(0, previous);
        alternate.insert(0, previous);
        record.pos -= 1;
    }
    record.ref_allele = String::from_utf8(reference_allele).unwrap_or_default();
    record.alt_allele = String::from_utf8(alternate).unwrap_or_default();
}

struct RawFilterOptions<'a> {
    reference_contigs: &'a BTreeSet<String>,
    fixchr: bool,
    pass_only: bool,
    regions: Option<&'a [vcf::BedInterval]>,
    targets: Option<&'a [vcf::BedInterval]>,
    locations: Option<&'a [vcf::LocationFilter]>,
}

fn filter_raw_records(
    records: Vec<vcf::RawVcfRecord>,
    path: &Path,
    options: &RawFilterOptions<'_>,
) -> Result<Vec<FilteredRawRecord>> {
    records
        .into_iter()
        .map(|record| {
            let key = vcf::VariantKey {
                chrom: somatic_chrom(&record.chrom, options.reference_contigs, options.fixchr),
                pos: record.pos,
                ref_allele: record.ref_allele.clone(),
                alt_allele: record.alt_allele.clone(),
            };
            if options.pass_only && !record.is_pass() {
                return Ok(None);
            }
            if calls_terminal_non_ref(&record) {
                return Ok(None);
            }
            let effective_end = record.effective_end_pos(path)?;
            if !vcf::matches_interval_filters(
                &key.chrom,
                key.pos,
                effective_end,
                options.regions,
                options.targets,
                options.locations,
            ) {
                return Ok(None);
            }
            let mut normalized = record;
            normalized.chrom = key.chrom.clone();
            Ok(Some(FilteredRawRecord {
                key,
                record: normalized,
            }))
        })
        .filter_map(|result| result.transpose())
        .collect()
}

fn somatic_chrom(chrom: &str, _reference_contigs: &BTreeSet<String>, fixchr: bool) -> String {
    if !fixchr {
        return chrom.to_string();
    }
    if chrom == "MT" || chrom == "chrMT" {
        return "chrM".to_string();
    }
    if chrom.starts_with([
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y', 'M',
    ]) {
        return format!("chr{chrom}");
    }
    // Non-canonical contigs are not rewritten by the legacy Perl expression.
    chrom.to_string()
}

fn calls_terminal_non_ref(record: &vcf::RawVcfRecord) -> bool {
    let alternates = record.alt_allele.split(',').collect::<Vec<_>>();
    if alternates.last() != Some(&"<NON_REF>") {
        return false;
    }
    let non_ref_index = alternates.len();
    record.samples.iter().any(|sample| {
        // The legacy streaming filter assumes GT is the first FORMAT value.
        sample
            .split(':')
            .next()
            .unwrap_or_default()
            .split(['/', '|'])
            .filter_map(|allele| allele.parse::<usize>().ok())
            .any(|allele| allele == non_ref_index)
    })
}

fn legacy_som_commandline(args: &SomaticArgs) -> String {
    let process_args = std::env::args().collect::<Vec<_>>();
    if let Some(subcommand) = process_args.iter().position(|value| value == "somatic") {
        let mut parts = vec!["/opt/conda/bin/som.py".to_string()];
        parts.extend(process_args.into_iter().skip(subcommand + 1));
        return parts.join(" ");
    }

    let mut parts = vec!["/opt/conda/bin/som.py".to_string()];

    if args.include_nonpass {
        parts.push("-P".to_string());
    }
    if args.count_unk {
        parts.push("--count-unk".to_string());
    }
    if args.happy_stats {
        parts.push("--happy-stats".to_string());
    }
    if args.af_strat {
        parts.push("--bin-afs".to_string());
    }

    parts.push(args.truth.clone());
    parts.push(args.query.clone());
    parts.push("-o".to_string());
    parts.push(args.output.clone());
    parts.push("--reference".to_string());
    parts.push(args.reference.clone());
    if let Some(fp_bed) = &args.fp_bedfile {
        parts.push("--false-positives".to_string());
        parts.push(fp_bed.clone());
    }

    if let Some(features) = &args.feature_table {
        parts.push("--feature-table".to_string());
        parts.push(features.clone());
    }

    parts.join(" ")
}

fn write_legacy_metrics_json(
    path: &Path,
    commandline: &str,
    stats_csv: &Path,
    ambiguous_classes: Option<&BTreeMap<String, usize>>,
    ambiguous_reasons: Option<&BTreeMap<String, usize>>,
) -> Result<()> {
    let text = fs::read_to_string(stats_csv)
        .with_context(|| format!("failed to read {}", stats_csv.display()))?;
    let mut lines = text.lines();
    let headers: Vec<String> = lines
        .next()
        .unwrap_or_default()
        .split(',')
        .map(str::to_string)
        .collect();
    let rows: Vec<Vec<String>> = lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split(',').map(str::to_string).collect())
        .collect();
    let metric_headers = headers[..headers.len().saturating_sub(2)].to_vec();
    let metric_rows: Vec<Vec<String>> = rows
        .iter()
        .map(|row| row[..row.len().saturating_sub(2)].to_vec())
        .collect();

    let mut body = String::new();
    body.push_str("{\"runInfo\": [{\"value\": ");
    body.push_str(&json_string(commandline));
    body.push_str(", \"key\": \"commandline\"}], \"metrics\": [");
    for (id, column, counts) in [
        ("ambiclasses", "class", ambiguous_classes),
        ("ambireasons", "reason", ambiguous_reasons),
    ] {
        if let Some(counts) = counts.filter(|counts| !counts.is_empty()) {
            body.push_str(&count_metric_json(id, column, counts));
            body.push_str(", ");
        }
    }
    body.push_str("{\"data\": [");
    body.push_str(&column_json(
        "types",
        "types",
        "string",
        &metric_rows.iter().map(|r| r[0].clone()).collect::<Vec<_>>(),
        true,
    ));
    for (i, header) in metric_headers.iter().enumerate().skip(1) {
        let values = metric_rows
            .iter()
            .map(|r| r.get(i).cloned().unwrap_or_default())
            .collect::<Vec<_>>();
        let kind = infer_type(&values);
        body.push_str(", ");
        body.push_str(&column_json(header, header, kind, &values, false));
    }
    body.push_str(
        "], \"properties\": [], \"type\": \"Table\", \"id\": \"result\", \"label\": \"result\"}], ",
    );
    body.push_str("\"version\": \"\", \"sampleInfo\": [], \"name\": \"som.py.comparison\", \"parameters\": [], ");
    body.push_str("\"timestamp\": ");
    body.push_str(&json_string(&iso_timestamp_now()));
    body.push_str(", \"metadata\": {\"required\": {\"version\": \"\", \"id\": \"haplotypes\", \"module\": \"som.py\", \"description\": ");
    body.push_str(&json_string(&format!(
        "/opt/conda/bin/som.py generated this JSON file via command line {}",
        commandline
    )));
    body.push_str("}}}");
    fs::write(path, body).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn count_metric_json(id: &str, column: &str, counts: &BTreeMap<String, usize>) -> String {
    let legacy_indices = python2_counter_indices(counts);
    let indices = counts
        .keys()
        .map(|key| legacy_indices[key].to_string())
        .collect::<Vec<_>>();
    let labels = counts.keys().cloned().collect::<Vec<_>>();
    let values = counts
        .values()
        .map(|count| count.to_string())
        .collect::<Vec<_>>();
    format!(
        "{{\"data\": [{}, {}, {}], \"properties\": [], \"type\": \"Table\", \"id\": {}, \"label\": {}}}",
        column_json("types", "types", "string", &indices, true),
        column_json(column, column, "string", &labels, false),
        column_json("count", "count", "int64", &values, false),
        json_string(id),
        json_string(id),
    )
}

fn python2_counter_indices(counts: &BTreeMap<String, usize>) -> BTreeMap<String, usize> {
    fn hash(value: &str) -> u64 {
        let bytes = value.as_bytes();
        let mut hash = bytes.first().copied().unwrap_or_default() as u64 * 128;
        for byte in bytes {
            hash = hash.wrapping_mul(1_000_003) ^ u64::from(*byte);
        }
        hash ^= bytes.len() as u64;
        if hash == u64::MAX { u64::MAX - 1 } else { hash }
    }

    fn insert(table: &mut [Option<String>], key: String) {
        let mask = table.len() - 1;
        let hashed = hash(&key);
        let mut slot = hashed as usize & mask;
        let mut perturb = hashed;
        while table[slot].is_some() {
            slot = slot
                .wrapping_mul(5)
                .wrapping_add(perturb as usize)
                .wrapping_add(1)
                & mask;
            perturb >>= 5;
        }
        table[slot] = Some(key);
    }

    let mut table = vec![None; 8];
    let mut used = 0usize;
    for key in counts.keys() {
        insert(&mut table, key.clone());
        used += 1;
        if used * 3 >= table.len() * 2 {
            let minimum = if used > 50_000 { used * 2 } else { used * 4 };
            let mut size = 8usize;
            while size <= minimum {
                size *= 2;
            }
            let old = std::mem::replace(&mut table, vec![None; size]);
            for key in old.into_iter().flatten() {
                insert(&mut table, key);
            }
        }
    }

    table
        .into_iter()
        .flatten()
        .enumerate()
        .map(|(index, key)| (key, index))
        .collect()
}

fn infer_type(values: &[String]) -> &'static str {
    let is_missing = |value: &str| value.is_empty() || value == ".";
    let populated = values
        .iter()
        .filter(|value| !is_missing(value))
        .collect::<Vec<_>>();
    let has_missing = populated.len() != values.len();
    if !populated.is_empty() && populated.iter().all(|value| value.parse::<i64>().is_ok()) {
        // pandas promotes integer columns containing NaN to float64.
        if has_missing { "double" } else { "int64" }
    } else if !populated.is_empty() && populated.iter().all(|value| value.parse::<f64>().is_ok()) {
        "double"
    } else {
        "string"
    }
}

fn column_json(
    id: &str,
    label: &str,
    kind: &str,
    values: &[String],
    numeric_strings: bool,
) -> String {
    let rendered = values
        .iter()
        .map(|value| match kind {
            "int64" | "double" if value.is_empty() || value == "." => "null".to_string(),
            "int64" | "double" => value.to_string(),
            _ if numeric_strings && value.parse::<i64>().is_ok() => value.to_string(),
            _ => json_string(value),
        })
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{{\"values\": [{rendered}], \"type\": {}, \"id\": {}, \"label\": {}}}",
        json_string(kind),
        json_string(id),
        json_string(label)
    )
}

fn json_string(value: &str) -> String {
    let escaped = value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n");
    format!("\"{escaped}\"")
}

/// Current UTC time formatted the same way Python's `datetime.datetime.now().isoformat()`
/// renders it (microsecond precision, no timezone suffix) — matches legacy `som.py` JSON.
fn iso_timestamp_now() -> String {
    let duration = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    let secs = duration.as_secs() as i64;
    let micros = duration.subsec_micros();

    let days = secs.div_euclid(86_400);
    let time_of_day = secs.rem_euclid(86_400);
    let hour = time_of_day / 3_600;
    let minute = (time_of_day / 60) % 60;
    let second = time_of_day % 60;

    // Civil-from-days algorithm (Howard Hinnant, public domain).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if m <= 2 { y + 1 } else { y };

    format!("{year:04}-{m:02}-{d:02}T{hour:02}:{minute:02}:{second:02}.{micros:06}")
}

fn jeffreys_ci(x: usize, n: usize, alpha: f64) -> (f64, f64, f64) {
    if n == 0 {
        return (0.0, 0.0, 1.0);
    }
    let p = x as f64 / n as f64;
    let lower = if x == n {
        (alpha / 2.0).powf(1.0 / n as f64)
    } else if x <= 1 {
        0.0
    } else {
        beta_ppf(alpha / 2.0, x as f64 + 0.5, (n - x) as f64 + 0.5)
    };
    let upper = if x == 0 {
        1.0 - (alpha / 2.0).powf(1.0 / n as f64)
    } else if x >= n.saturating_sub(1) {
        1.0
    } else {
        beta_isf(alpha / 2.0, x as f64 + 0.5, (n - x) as f64 + 0.5)
    };
    (p, lower, upper)
}

fn beta_isf(q: f64, a: f64, b: f64) -> f64 {
    beta_ppf(1.0 - q, a, b)
}

fn beta_ppf(p: f64, a: f64, b: f64) -> f64 {
    // Bit-exact port of scipy 1.2.1's beta.ppf via cephes_incbi.
    // Used for `Tools/ci.py::jeffreys` confidence intervals so sompy
    // recall_lower / precision_lower / etc. match legacy hap.py output.
    crate::cephes::incbi(a, b, p)
}

#[allow(dead_code)]
fn beta_ppf_legacy_binsearch(p: f64, a: f64, b: f64) -> f64 {
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }
    let mut lo = 0.0f64;
    let mut hi = 1.0f64;
    for _ in 0..256 {
        let mid = (lo + hi) / 2.0;
        let cdf = regularized_beta(a, b, mid);
        if cdf > p {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

fn regularized_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    let symm_transform = x >= (a + 1.0) / (a + b + 2.0);
    let eps = 1.110_223_024_625_156_5e-16;
    let fpmin = f64::MIN_POSITIVE / eps;

    let (mut a_work, mut b_work, mut x_work) = (a, b, x);
    if symm_transform {
        x_work = 1.0 - x_work;
        (a_work, b_work) = (b_work, a_work);
    }

    let qab = a_work + b_work;
    let qap = a_work + 1.0;
    let qam = a_work - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x_work / qap;

    if d.abs() < fpmin {
        d = fpmin;
    }
    d = 1.0 / d;
    let mut h = d;

    for m in 1..141 {
        let m = f64::from(m);
        let m2 = m * 2.0;
        let mut aa = m * (b_work - m) * x_work / ((qam + m2) * (a_work + m2));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        h *= d * c;

        aa = -(a_work + m) * (qab + m) * x_work / ((a_work + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = 1.0 + aa / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;

        if (del - 1.0).abs() <= eps {
            return if symm_transform {
                1.0 - bt * h / a_work
            } else {
                bt * h / a_work
            };
        }
    }

    if symm_transform {
        1.0 - bt * h / a_work
    } else {
        bt * h / a_work
    }
}

#[allow(clippy::excessive_precision)] // Preserve the source approximation coefficients exactly.
fn ln_gamma(z: f64) -> f64 {
    const GAMMA_R: f64 = 10.900_511;
    const LN_PI: f64 = 1.144_729_885_849_400_174_143_427_351_353_058_711_647_294_812_915_3;
    const LN_2_SQRT_E_OVER_PI: f64 =
        0.620_782_237_635_245_222_345_518_445_781_647_212_251_852_727_902_597_8;
    const GAMMA_DK: [f64; 11] = [
        2.485_740_891_387_535_655_46e-5,
        1.051_423_785_817_219_742_10,
        -3.456_870_972_220_162_354_69,
        4.512_277_094_668_948_237_00,
        -2.982_852_253_235_766_557_21,
        1.056_397_115_771_267_130_77,
        -1.954_287_731_916_458_695_83e-1,
        1.709_705_434_044_412_243_07e-2,
        -5.719_261_174_043_057_812_83e-4,
        4.633_994_733_599_056_367_08e-6,
        -2.719_949_084_886_077_039_10e-9,
    ];

    if z < 0.5 {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| s + t.1 / (t.0 as f64 - z));
        LN_PI
            - (std::f64::consts::PI * z).sin().ln()
            - s.ln()
            - LN_2_SQRT_E_OVER_PI
            - (0.5 - z) * ((0.5 - z + GAMMA_R) / std::f64::consts::E).ln()
    } else {
        let s = GAMMA_DK
            .iter()
            .enumerate()
            .skip(1)
            .fold(GAMMA_DK[0], |s, t| s + t.1 / (z + t.0 as f64 - 1.0));
        s.ln() + LN_2_SQRT_E_OVER_PI + (z - 0.5) * ((z - 0.5 + GAMMA_R) / std::f64::consts::E).ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Command};
    use clap::Parser;

    fn parsed_somatic(extra: &[&str]) -> SomaticArgs {
        let mut argv = vec![
            "hap",
            "somatic",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
        ];
        argv.extend_from_slice(extra);
        let cli = Cli::try_parse_from(argv).expect("test arguments should parse");
        let Command::Somatic(args) = cli.command else {
            panic!("somatic command expected");
        };
        args
    }

    fn interval(start: usize, end: usize, label: &str) -> AmbiguousInterval {
        AmbiguousInterval {
            interval: vcf::BedInterval {
                chrom: "chr1".to_string(),
                start,
                end,
            },
            label: label.to_string(),
            details: Vec::new(),
        }
    }

    fn raw_record(line: &str) -> vcf::RawVcfRecord {
        vcf::RawVcfRecord::from_line(line, Path::new("test.vcf")).expect("valid test VCF record")
    }

    fn unique_test_dir(label: &str) -> PathBuf {
        let id = SOMATIC_SCRATCH_RUN_ID.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "hap-somatic-test-{label}-{}-{id}",
            std::process::id()
        ))
    }

    fn write_test_vcf(path: &Path, position: usize) {
        fs::write(
            path,
            format!(
                "##fileformat=VCFv4.1\n##contig=<ID=chr1,length=20>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t{position}\t.\tA\tC\t.\tPASS\t.\n"
            ),
        )
        .expect("write test VCF");
    }

    #[test]
    fn count_unk_without_fp_regions_classifies_unmatched_calls_as_unknown() {
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &[], true, false),
            QueryClass::Unk
        );
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &[], false, false),
            QueryClass::Fp
        );
    }

    #[test]
    fn explicit_negative_toggle_wins() {
        assert!(resolve_toggle(true, false));
        assert!(!resolve_toggle(true, true));
        assert!(!resolve_toggle(false, false));
    }

    #[test]
    fn paired_somatic_toggles_use_last_token_wins_precedence() {
        let args = parsed_somatic(&["--count-unk", "--no-count-unk"]);
        assert!(!resolve_toggle(args.count_unk, args.no_count_unk));
        let args = parsed_somatic(&["--no-count-unk", "--count-unk"]);
        assert!(resolve_toggle(args.count_unk, args.no_count_unk));

        let args = parsed_somatic(&["--ambi-fp", "--no-ambi-fp"]);
        assert!(!resolve_toggle(args.ambi_fp, args.no_ambi_fp));
        let args = parsed_somatic(&["--no-ambi-fp", "--ambi-fp"]);
        assert!(resolve_toggle(args.ambi_fp, args.no_ambi_fp));
    }

    #[test]
    fn fixchr_pairs_use_last_token_wins_precedence() {
        let args = parsed_somatic(&["--no-fixchr-truth", "--fixchr-truth"]);
        assert!(args.fixchr_truth.unwrap_or(true) && !args.no_fixchr_truth);
        let args = parsed_somatic(&["--fixchr-truth", "--no-fixchr-truth"]);
        assert!(!args.fixchr_truth.unwrap_or(true) || args.no_fixchr_truth);

        let args = parsed_somatic(&["--no-fixchr-query", "--fixchr-query"]);
        assert!(args.fixchr_query.unwrap_or(true) && !args.no_fixchr_query);
        let args = parsed_somatic(&["--fixchr-query", "--no-fixchr-query"]);
        assert!(!args.fixchr_query.unwrap_or(true) || args.no_fixchr_query);
    }

    #[test]
    fn normalize_truth_query_and_all_select_the_legacy_inputs() {
        let truth = parsed_somatic(&["--normalize-truth"]);
        assert_eq!(selected_normalizations(&truth), (true, false));
        let query = parsed_somatic(&["--normalize-query"]);
        assert_eq!(selected_normalizations(&query), (false, true));
        let mut all = parsed_somatic(&["--normalize-all"]);
        all.count_filtered_fn = true;
        assert_eq!(
            selected_normalizations(&all),
            (true, true),
            "-FN reporting must remain independent of -N normalization"
        );
    }

    #[test]
    fn governed_somatic_defaults_remain_supported() {
        validate_args(&parsed_somatic(&[])).expect("default comparison must remain supported");
        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "hcc.strelka.indel",
        ]))
        .expect("nf-test feature-table comparison must remain supported");
        validate_args(&parsed_somatic(&["--feature-table", "admix.strelka.snv"]))
            .expect("caller-specific SNV comparison must remain supported");
    }

    #[test]
    fn transformation_and_reporting_controls_are_accepted() {
        let cases: &[&[&str]] = &[
            &["--normalize-truth"],
            &["--normalize-query"],
            &["--normalize-all"],
            &["--no-fixchr-truth"],
            &["--no-fixchr-query"],
            &["--roc", "strelka.snv"],
        ];
        for case in cases {
            validate_args(&parsed_somatic(case))
                .unwrap_or_else(|error| panic!("legacy control {case:?} was rejected: {error}"));
        }
    }

    #[test]
    fn scratch_lifecycle_logging_and_verbosity_are_operational() {
        let root = unique_test_dir("controls");
        let scratch = root.join("explicit-scratch");
        let logfile = root.join("som.log");
        fs::create_dir_all(&root).expect("create operational test root");

        let mut args = parsed_somatic(&[]);
        args.scratch_prefix = Some(scratch.display().to_string());
        args.logfile = Some(logfile.display().to_string());
        args.verbose = true;
        let mut controls = SomaticOperationalControls::prepare(&args).expect("prepare controls");
        assert!(controls.scratch.path.is_dir());
        assert!(!controls.should_print_summary());
        controls.info("operational log marker").expect("write log");
        controls.scratch.cleanup().expect("retain explicit scratch");
        assert!(scratch.is_dir());
        assert!(
            fs::read_to_string(&logfile)
                .expect("read logfile")
                .contains("operational log marker")
        );

        let mut quiet_args = parsed_somatic(&[]);
        quiet_args.quiet = true;
        let quiet = SomaticOperationalControls::prepare(&quiet_args).expect("prepare quiet mode");
        assert!(!quiet.should_print_summary());
        let quiet_path = quiet.scratch.path.clone();
        quiet.scratch.cleanup().expect("clean quiet scratch");
        assert!(!quiet_path.exists());

        let mut keep_args = parsed_somatic(&[]);
        keep_args.keep_scratch = true;
        let keep = SomaticOperationalControls::prepare(&keep_args).expect("prepare kept scratch");
        let keep_path = keep.scratch.path.clone();
        keep.scratch.cleanup().expect("keep scratch");
        assert!(keep_path.is_dir());

        fs::remove_dir_all(&keep_path).expect("remove retained generated scratch");
        fs::remove_dir_all(&root).expect("remove operational test root");
    }

    #[test]
    fn continue_reuses_cached_normalized_inputs() {
        let root = unique_test_dir("continue");
        let scratch = root.join("scratch");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let reference = root.join("reference.fa");
        fs::create_dir_all(&root).expect("create continue test root");
        fs::write(&reference, ">chr1\nAAAAAAAAAAAAAAAAAAAA\n").expect("write reference");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 7);

        let comparison = |output: &Path, cont: bool| {
            let mut args = parsed_somatic(&[]);
            args.truth = truth.display().to_string();
            args.query = query.display().to_string();
            args.reference = reference.display().to_string();
            args.output = output.display().to_string();
            args.scratch_prefix = Some(scratch.display().to_string());
            args.cont = cont;
            args.quiet = true;
            run(args).expect("run somatic comparison");
        };

        comparison(&root.join("first"), false);
        assert!(scratch.join("normalized_truth.vcf.gz").is_file());
        assert!(scratch.join("normalized_query.vcf.gz").is_file());

        write_test_vcf(&query, 8);
        comparison(&root.join("continued"), true);
        let stats =
            fs::read_to_string(root.join("continued.stats.csv")).expect("read continued stats");
        let snv = stats
            .lines()
            .find(|line| line.starts_with("1,SNVs,"))
            .expect("continued SNV row");
        assert!(snv.starts_with("1,SNVs,1,1,1,0,0,0,0,"));

        fs::remove_dir_all(&root).expect("remove continue test root");
    }

    #[test]
    fn explain_ambiguous_without_features_writes_csv_and_metrics_tables() {
        let root = unique_test_dir("explain-no-features");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let ambiguous = root.join("ambiguous.bed");
        fs::create_dir_all(&root).expect("create explanation test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 8);
        fs::write(
            &ambiguous,
            "chr1\t7\t8\tignored\tunk\t2\tignored\tlow-vaf\n",
        )
        .expect("write ambiguous BED");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.ambiguous_beds = vec![ambiguous.display().to_string()];
        args.explain_ambiguous = true;
        args.fp_region_size = Some("10".to_string());
        args.quiet = true;
        run(args).expect("explanation run must not require a reference or feature table");

        assert!(!root.join("result.features.csv").exists());
        assert!(
            fs::read_to_string(root.join("result.ambiclasses.csv"))
                .expect("read ambiguity classes")
                .contains("ambi-unk,1")
        );
        assert!(
            fs::read_to_string(root.join("result.ambireasons.csv"))
                .expect("read ambiguity reasons")
                .contains("ambi-unk: low-vaf,1")
        );
        let metrics =
            fs::read_to_string(root.join("result.metrics.json")).expect("read explanation metrics");
        assert!(metrics.contains("\"id\": \"ambiclasses\""));
        assert!(metrics.contains("\"id\": \"ambireasons\""));
        assert!(metrics.contains("ambi-unk: low-vaf"));

        fs::remove_dir_all(&root).expect("remove explanation test root");
    }

    #[test]
    fn usable_fp_bed_avoids_loading_a_missing_reference() {
        let root = unique_test_dir("lazy-reference-fp-bed");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let fp = root.join("fp.bed");
        fs::create_dir_all(&root).expect("create lazy-reference test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 8);
        fs::write(&fp, "chr1\t0\t20\n").expect("write FP BED");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.fp_bedfile = Some(fp.display().to_string());
        args.quiet = true;
        run(args).expect("usable FP BED must avoid loading the missing reference");

        let stats = fs::read_to_string(root.join("result.stats.csv")).expect("read stats");
        let snv = stats
            .lines()
            .find(|line| line.starts_with("1,SNVs,"))
            .expect("SNV row");
        assert!(snv.contains(",20,50000.0,"));

        fs::remove_dir_all(&root).expect("remove lazy-reference test root");
    }

    #[test]
    fn normalize_all_with_filtered_fn_still_requires_and_uses_reference() {
        let root = unique_test_dir("normalize-all-filtered-fn");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        fs::create_dir_all(&root).expect("create normalization test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 7);

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.normalize_all = true;
        args.count_filtered_fn = true;
        args.include_nonpass = true;
        args.feature_table = Some("generic".to_string());
        args.fp_region_size = Some("10".to_string());
        args.quiet = true;
        let error = run(args).expect_err("normalization must load the missing reference");
        assert!(error.to_string().contains("failed to read"));

        fs::remove_dir_all(&root).expect("remove normalization test root");
    }

    #[test]
    fn af_controls_match_legacy_validation_independently_of_happy_stats() {
        validate_args(&parsed_somatic(&[
            "--feature-table",
            "hcc.strelka.indel",
            "--bin-afs",
        ]))
        .expect("AF stats rows do not require --happy-stats");

        validate_args(&parsed_somatic(&[
            "--af-binsize",
            "0.1",
            "--af-truth",
            "TRUTH_AF",
            "--af-query",
            "QUERY_AF",
        ]))
        .expect("AF controls are inert unless --bin-afs is selected");

        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "generic",
            "--happy-stats",
            "--bin-afs",
        ]))
        .expect("legacy validates selected AF columns after extracting the feature table");

        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "hcc.strelka.indel",
            "--happy-stats",
            "--bin-afs",
            "--af-query",
            "MISSING_AF",
        ]))
        .expect("legacy validates missing AF selectors after input I/O");

        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "hcc.strelka.indel",
            "--happy-stats",
            "--bin-afs",
            "--af-binsize",
            "0.25",
            "--af-truth",
            "I.T_ALT_RATE",
            "--af-query",
            "T_AF",
        ]))
        .expect("implemented AF extended-summary controls should remain supported");
    }

    #[test]
    fn ambiguous_fp_toggle_matches_legacy_label_semantics() {
        let ambiguous = vec![
            interval(10, 20, "fp"),
            interval(30, 40, "Fp"),
            interval(50, 60, "FP"),
        ];
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &ambiguous, true, false),
            QueryClass::Ambi
        );
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &ambiguous, true, true),
            QueryClass::Fp
        );
        assert_eq!(
            classify_query("chr1", 31, 31, &[], &ambiguous, true, true),
            QueryClass::Ambi,
            "legacy labels are case-sensitive"
        );
        assert_eq!(
            classify_query("chr1", 51, 51, &[], &ambiguous, true, false),
            QueryClass::Fp,
            "uppercase FP is unconditional"
        );
    }

    #[test]
    fn ambiguous_bed_uses_the_fifth_column_as_the_label() {
        let mut bed = tempfile::NamedTempFile::new().expect("temporary BED");
        writeln!(bed, "chr1\t0\t10\tignored\tFP\textra").expect("write BED");
        let contigs = BTreeSet::from(["chr1".to_string()]);
        let intervals =
            load_ambiguous_beds(&[bed.path().to_string_lossy().into_owned()], &contigs, true)
                .expect("load ambiguous BED");
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].label, "FP");
        assert_eq!(intervals[0].details, ["ignored", "FP", "extra"]);
    }

    #[test]
    fn classification_beds_follow_the_truth_fixchr_switch() {
        let root = unique_test_dir("classification-bed");
        let fp = root.join("fp.bed");
        let ambiguous = root.join("ambiguous.bed");
        fs::create_dir_all(&root).expect("create BED test root");
        fs::write(&fp, "1\t0\t10\n").expect("write FP BED");
        fs::write(&ambiguous, "1\t0\t10\tignored\tfp\n").expect("write ambiguous BED");
        let contigs = BTreeSet::from(["chr1".to_string()]);

        assert_eq!(
            load_classification_bed(&fp, &contigs, true).unwrap()[0].chrom,
            "chr1"
        );
        assert_eq!(
            load_classification_bed(&fp, &contigs, false).unwrap()[0].chrom,
            "1"
        );
        assert_eq!(
            load_ambiguous_beds(&[ambiguous.display().to_string()], &contigs, true).unwrap()[0]
                .interval
                .chrom,
            "chr1"
        );
        assert_eq!(
            load_ambiguous_beds(&[ambiguous.display().to_string()], &contigs, false).unwrap()[0]
                .interval
                .chrom,
            "1"
        );

        fs::remove_dir_all(&root).expect("remove BED test root");
    }

    #[test]
    fn ambiguous_explanation_records_classes_and_reasons() {
        let ambiguous = vec![AmbiguousInterval {
            interval: vcf::BedInterval {
                chrom: "chr1".to_string(),
                start: 10,
                end: 20,
            },
            label: "unk".to_string(),
            details: vec![
                "2".to_string(),
                "ignored".to_string(),
                "low-vaf".to_string(),
            ],
        }];
        let mut classes = BTreeMap::new();
        let mut reasons = BTreeMap::new();
        record_ambiguous_explanation(
            "chr1",
            11,
            11,
            &ambiguous,
            false,
            &mut classes,
            &mut reasons,
        );
        assert_eq!(classes.get("ambi-unk"), Some(&1));
        assert_eq!(reasons.get("ambi-unk: rep. count 2"), Some(&1));
        assert_eq!(reasons.get("ambi-unk: low-vaf"), Some(&1));
    }

    #[test]
    fn ambiguity_reason_csv_preserves_python2_counter_indices() {
        let counts = BTreeMap::from([
            ("ambi-unk: 2".to_string(), 1),
            ("ambi-unk: ignored".to_string(), 1),
            ("ambi-unk: low-vaf".to_string(), 1),
            ("ambi-unk: rep. count ignored".to_string(), 1),
        ]);
        let output = tempfile::NamedTempFile::new().expect("temporary ambiguity CSV");
        write_legacy_count_table(output.path(), "reason", &counts)
            .expect("write ambiguity reasons");
        assert_eq!(
            fs::read_to_string(output.path()).expect("read ambiguity reasons"),
            concat!(
                ",reason,count\n",
                "2,ambi-unk: 2,1\n",
                "1,ambi-unk: ignored,1\n",
                "3,ambi-unk: low-vaf,1\n",
                "0,ambi-unk: rep. count ignored,1\n",
            )
        );
        let metric = count_metric_json("ambireasons", "reason", &counts);
        assert!(metric.contains("\"values\": [2, 1, 3, 0]"));
    }

    #[test]
    fn explicit_fp_regions_take_priority_over_ambiguous_regions() {
        let fp = vec![vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 10,
            end: 20,
        }];
        let ambiguous = vec![interval(10, 20, "unk")];
        assert_eq!(
            classify_query("chr1", 11, 11, &fp, &ambiguous, true, false),
            QueryClass::Fp
        );
    }

    #[test]
    fn raw_filtering_uses_symbolic_end_for_regions_and_start_for_targets() {
        let path = Path::new("query.vcf");
        let record =
            vcf::RawVcfRecord::from_line("chr1\t3\t.\tC\t<DEL>\t.\tPASS\tEND=5\tGT\t0/1", path)
                .unwrap();
        let contigs = BTreeSet::from(["chr1".to_string()]);
        let boundary = vec![vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 4,
            end: 5,
        }];

        let by_region = filter_raw_records(
            vec![record.clone()],
            path,
            &RawFilterOptions {
                reference_contigs: &contigs,
                fixchr: true,
                pass_only: false,
                regions: Some(&boundary),
                targets: None,
                locations: None,
            },
        )
        .unwrap();
        let by_target = filter_raw_records(
            vec![record.clone()],
            path,
            &RawFilterOptions {
                reference_contigs: &contigs,
                fixchr: true,
                pass_only: false,
                regions: None,
                targets: Some(&boundary),
                locations: None,
            },
        )
        .unwrap();
        let by_location = filter_raw_records(
            vec![record],
            path,
            &RawFilterOptions {
                reference_contigs: &contigs,
                fixchr: true,
                pass_only: false,
                regions: None,
                targets: None,
                locations: Some(&[vcf::LocationFilter::Range {
                    chrom: "chr1".to_string(),
                    start: 3,
                    end: 3,
                }]),
            },
        )
        .unwrap();

        assert_eq!(by_region.len(), 1, "-R must use INFO/END");
        assert!(by_target.is_empty(), "-T must use POS");
        assert_eq!(by_location.len(), 1, "-l must use POS");
    }

    #[test]
    fn fp_region_size_is_restricted_by_location() {
        let fp_regions = vec![vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 0,
            end: 50,
        }];
        let references = BTreeMap::from([("chr1".to_string(), "A".repeat(100))]);
        let range = [vcf::LocationFilter::Range {
            chrom: "chr1".to_string(),
            start: 21,
            end: 30,
        }];
        let contig = [vcf::LocationFilter::Contig("chr1".to_string())];

        assert_eq!(
            calculate_fp_region_size(None, &fp_regions, Some(&range), &references, &[], &[]),
            10
        );
        assert_eq!(
            calculate_fp_region_size(None, &fp_regions, Some(&contig), &references, &[], &[]),
            50
        );
        assert_eq!(
            calculate_fp_region_size(None, &[], Some(&range), &references, &[], &[]),
            10
        );
        assert_eq!(
            calculate_fp_region_size(Some("7"), &fp_regions, Some(&range), &references, &[], &[],),
            7
        );
        assert!(!fp_region_size_requires_reference(Some("7"), &[]));
        assert!(!fp_region_size_requires_reference(None, &fp_regions));
        assert!(fp_region_size_requires_reference(None, &[]));
        assert!(fp_region_size_requires_reference(Some("auto"), &[]));
    }

    #[test]
    fn fp_classification_uses_the_reference_span_not_symbolic_end() {
        let record = raw_record("chr1\t10\t.\tA\t<DEL>\t.\tPASS\tEND=100");
        let fp = vec![vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 49,
            end: 50,
        }];
        assert_eq!(
            classify_query(
                &record.chrom,
                record.pos,
                record.end_pos(),
                &fp,
                &[],
                true,
                false,
            ),
            QueryClass::Unk
        );
        assert_eq!(
            classify_query(
                &record.chrom,
                record.pos,
                record.effective_end_pos(Path::new("test.vcf")).unwrap(),
                &fp,
                &[],
                true,
                false,
            ),
            QueryClass::Fp,
            "the symbolic span would have produced the wrong class"
        );
    }

    #[test]
    fn terminal_non_ref_is_removed_only_when_its_allele_is_called() {
        assert!(calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\tC,<NON_REF>\t.\tPASS\t.\tGT\t0/2"
        )));
        assert!(!calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\tC,<NON_REF>\t.\tPASS\t.\tGT\t0/1"
        )));
        assert!(!calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\t<NON_REF>,C\t.\tPASS\t.\tGT\t0/1"
        )));
        assert!(!calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\tC,<NON_REF>\t.\tPASS\t."
        )));
    }

    #[test]
    fn sites_only_records_survive_raw_filtering() {
        let record = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\t.");
        let filtered = filter_raw_records(
            vec![record],
            Path::new("test.vcf"),
            &RawFilterOptions {
                reference_contigs: &BTreeSet::from(["chr1".to_string()]),
                fixchr: true,
                pass_only: true,
                regions: None,
                targets: None,
                locations: None,
            },
        )
        .expect("filter sites-only VCF");
        assert_eq!(filtered.len(), 1);
        assert!(filtered[0].record.samples.is_empty());
    }

    #[test]
    fn exact_pairing_preserves_duplicate_occurrences() {
        let record = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\t.");
        let filtered = |record: vcf::RawVcfRecord| FilteredRawRecord {
            key: vcf::VariantKey {
                chrom: record.chrom.clone(),
                pos: record.pos,
                ref_allele: record.ref_allele.clone(),
                alt_allele: record.alt_allele.clone(),
            },
            record,
        };
        let truth = vec![filtered(record.clone()), filtered(record.clone())];
        let query = vec![filtered(record.clone()), filtered(record)];
        let (truth_matches, query_matches) = pair_exact_records(&truth, &query);
        assert_eq!(truth_matches, vec![Some(0), Some(1)]);
        assert_eq!(query_matches, vec![Some(0), Some(1)]);
    }

    #[test]
    fn caller_feature_merge_suffixes_truth_and_preserves_numeric_csv_types() {
        let truth = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\teditDistance=1\tGT\t0/0\t0/1");
        let query = raw_record(
            "chr1\t10\t.\tA\tC\t60\tPASS\tNT=ref;QSS_NT=10;VQSR=2;SomaticEVS=3;MQ=50;MQ0=1;SNVSB=0.2;ReadPosRankSum=0.3\tGT:DP:FDP:SDP:AU:CU:GU:TU\t0/0:20:2:1:10,0:0,0:0,0:0,0\t0/1:30:3:2:15,0:5,0:0,0:0,0",
        );
        let table = build_caller_feature_table(
            "admix.strelka.snv",
            &[],
            &[],
            None,
            true,
            &CallerRecordGroups {
                tp_truth: &[truth],
                tp_query: &[query],
                fn_truth: &[],
                fp_query: &[],
                ambi_query: &[],
                unk_query: &[],
            },
        )
        .expect("caller feature table");
        let headers = parse_csv_line(&table.header);
        let row = parse_csv_line(&table.tp[0]);
        let field = |name: &str| {
            let index = headers.iter().position(|header| header == name).unwrap();
            row[index].as_str()
        };
        assert_eq!(
            &headers[..8],
            [
                "",
                "CHROM",
                "POS",
                "tag",
                "REF",
                "REF.truth",
                "ALT",
                "ALT.truth"
            ]
        );
        assert_eq!(field("tag"), "TP");
        assert_eq!(field("REF.truth"), "A");
        assert_eq!(field("I.editDistance"), "1.00000000");
        assert_eq!(field("QSS_NT"), "10.00000000");
        assert_eq!(field("S.2.GT"), "0/1");
        assert_eq!(field("POS"), "10");
    }

    #[test]
    fn no_order_check_controls_caller_tp_order_validation() {
        let truth = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1");
        let query = raw_record(
            "chr1\t11\t.\tA\tC\t60\tPASS\tNT=ref;QSS_NT=10\tGT:DP:FDP:SDP:AU:CU:GU:TU\t0/0:20:0:0:20,0:0,0:0,0:0,0\t0/1:30:0:0:15,0:15,0:0,0:0,0",
        );
        let groups = CallerRecordGroups {
            tp_truth: std::slice::from_ref(&truth),
            tp_query: std::slice::from_ref(&query),
            fn_truth: &[],
            fp_query: &[],
            ambi_query: &[],
            unk_query: &[],
        };

        let error =
            match build_caller_feature_table("hcc.strelka.snv", &[], &[], None, true, &groups) {
                Ok(_) => panic!("default order check must reject mismatched TP rows"),
                Err(error) => error,
            };
        assert!(error.to_string().contains("out of order"));
        build_caller_feature_table("hcc.strelka.snv", &[], &[], None, false, &groups)
            .expect("--no-order-check must bypass the developer safety check");
    }

    #[cfg(feature = "verification")]
    #[test]
    fn bam_depths_flow_into_somatic_caller_features() {
        let fixture_dir =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("verification/assets/fixtures/ftx-bam");
        let vcf_path = fixture_dir.join("input.vcf");
        let (headers, records) = vcf::load_raw_vcf(&vcf_path).expect("load BAM parity VCF");
        let depths =
            ftx::bam_normalization_depths(&[fixture_dir.join("reads.bam").display().to_string()])
                .expect("scan BAM normalization depths");

        let table = build_caller_feature_table(
            "hcc.strelka.snv",
            &headers,
            &headers,
            Some(&depths),
            true,
            &CallerRecordGroups {
                tp_truth: &records,
                tp_query: &records,
                fn_truth: &[],
                fp_query: &[],
                ambi_query: &[],
                unk_query: &[],
            },
        )
        .expect("caller feature table with BAM depths");
        let columns = parse_csv_line(&table.header);
        let row = parse_csv_line(&table.tp[0]);
        let field = |name: &str| {
            let index = columns.iter().position(|column| column == name).unwrap();
            row[index].as_str()
        };

        assert!((depths["chr1"] - 0.9).abs() < f64::EPSILON);
        assert_eq!(field("N_DP_RATE"), "10.00000000");
        assert_eq!(field("T_DP_RATE"), "20.00000000");
    }

    #[test]
    fn raw_stats_types_match_legacy_labels_and_indexes() {
        assert_eq!(STATS_TYPE_ROWS[0], (0, "indels"));
        assert_eq!(STATS_TYPE_ROWS[1], (1, "SNVs"));
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\tC\t.\tPASS\t.")),
            Some("SNVs")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\tAC\t.\tPASS\t.")),
            Some("indels")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tAC\tGT\t.\tPASS\t.")),
            Some("MNPs")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\t<DEL>\t.\tPASS\t.")),
            Some("others")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\t.\t.\tPASS\t.")),
            None
        );
    }

    #[test]
    fn confidence_level_changes_interval_without_changing_point_estimate() {
        let counts = SomaticCounts {
            truth_total: 10,
            query_total: 10,
            tp: 8,
            fp: 2,
            fn_count: 2,
            ..SomaticCounts::default()
        };
        let ci_95 = render_row(
            1,
            "SNVs",
            counts,
            &StatsRowContext {
                fp_region_size: 100,
                ci_alpha: 0.05,
                filtered: None,
                include_filtered_columns: false,
                commandline: "som.py",
            },
        );
        let ci_80 = render_row(
            1,
            "SNVs",
            counts,
            &StatsRowContext {
                fp_region_size: 100,
                ci_alpha: 0.20,
                filtered: None,
                include_filtered_columns: false,
                commandline: "som.py",
            },
        );
        let ci_95 = ci_95.split(',').collect::<Vec<_>>();
        let ci_80 = ci_80.split(',').collect::<Vec<_>>();
        assert_eq!(ci_95[9], ci_80[9]);
        assert_ne!(ci_95[10], ci_80[10]);
        assert_ne!(ci_95[11], ci_80[11]);
    }

    #[test]
    fn filtered_fn_columns_align_with_their_header() {
        let counts = SomaticCounts {
            truth_total: 10,
            query_total: 11,
            tp: 8,
            fp: 3,
            fn_count: 2,
            unk: 0,
            ambi: 0,
        };
        let filtered = FilteredCounts {
            tp: 1,
            fp: 2,
            ..FilteredCounts::default()
        };
        let header = stats_header(true).split(',').count();
        let row = render_row(
            1,
            "SNVs",
            counts,
            &StatsRowContext {
                fp_region_size: 100,
                ci_alpha: 0.05,
                filtered: Some(filtered),
                include_filtered_columns: true,
                commandline: "som.py",
            },
        );
        assert_eq!(header, row.split(',').count());
        assert!(row.contains(",2.0,1.0,0.0,0.0,"));
    }

    #[test]
    fn generic_filtered_fn_counts_default_each_reported_variant_type_to_zero() {
        let filtered = filtered_counts_for_type(true, Some("generic"), "indels", &BTreeMap::new())
            .expect("generic feature tables report filtered counts for every variant type");

        assert_eq!(
            (filtered.fp, filtered.tp, filtered.unk, filtered.ambi),
            (0, 0, 0, 0)
        );

        let row = render_row(
            0,
            "indels",
            SomaticCounts {
                truth_total: 1,
                query_total: 1,
                tp: 1,
                ..SomaticCounts::default()
            },
            &StatsRowContext {
                fp_region_size: 10,
                ci_alpha: 0.05,
                filtered: Some(filtered),
                include_filtered_columns: true,
                commandline: "som.py",
            },
        );
        assert_eq!(
            row,
            concat!(
                "0,indels,1,1,1,0,0,0,0,0.0,0.0,0.0,0.0,",
                "1.0,0.025,1.0,1.0,1.0,",
                "0.025,1.0,0.0,0.0,10,0.0,",
                "1.0,1.0,0.0,0.0,0.0,som.py-,som.py"
            )
        );

        let filtered_values = vec!["0.0".to_string(), "0.0".to_string()];
        assert_eq!(infer_type(&filtered_values), "double");
        assert_eq!(
            column_json(
                "fp.filtered",
                "fp.filtered",
                infer_type(&filtered_values),
                &filtered_values,
                false,
            ),
            "{\"values\": [0.0, 0.0], \"type\": \"double\", \"id\": \"fp.filtered\", \"label\": \"fp.filtered\"}"
        );
    }

    #[test]
    fn af_extended_uses_the_selected_truth_and_query_fields() {
        let output = tempfile::NamedTempFile::new().expect("temporary output");
        let header = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,FILTER,TRUTH_AF,QUERY_AF";
        let rows = vec![
            "0,chr1,10,TP,A,A,C,C,,0.1,0.9".to_string(),
            "0,chr1,20,FN,,A,,C,,0.8,".to_string(),
            "0,chr1,30,FP,A,,C,,,0.7,0.1".to_string(),
        ];
        write_happy_style_extended(
            output.path(),
            header,
            &rows,
            "hcc.strelka.snv",
            "0.5",
            "TRUTH_AF",
            "QUERY_AF",
        )
        .expect("extended output");
        let text = fs::read_to_string(output.path()).expect("read extended output");
        assert!(text.contains("SNP,*,\"[0.00,0.50)\",PASS,1,1,0,2,1,0"));
        assert!(text.contains("SNP,*,\"[0.50,1.00]\",PASS,1,0,1,0,0,0"));
    }

    #[test]
    fn af_stats_use_truth_af_for_tp_fn_and_query_af_for_other_tags() {
        let header = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,FILTER,TRUTH_AF,QUERY_AF";
        let rows = vec![
            "0,chr1,10,TP,A,A,C,C,,0.1,0.9".to_string(),
            "0,chr1,20,FN,,A,,C,,0.8,".to_string(),
            "0,chr1,30,FP,A,,C,,LowQual,,0.1".to_string(),
            "0,chr1,40,UNK,A,,C,,,,0.7".to_string(),
            "0,chr1,50,AMBI,A,,C,,,,1.0".to_string(),
        ];
        let bins =
            calculate_af_stats(header, &rows, "0.5", "TRUTH_AF", "QUERY_AF").expect("AF stats");
        assert_eq!(bins.len(), 2);
        assert_eq!(bins[0].2.truth_total, 1);
        assert_eq!(bins[0].2.query_total, 2);
        assert_eq!(bins[0].2.tp, 1);
        assert_eq!(bins[0].2.fp, 1);
        assert_eq!(bins[0].3.fp, 1);
        assert_eq!(bins[1].2.truth_total, 1);
        assert_eq!(bins[1].2.fn_count, 1);
        assert_eq!(bins[1].2.query_total, 2);
        assert_eq!(bins[1].2.unk, 1);
        assert_eq!(bins[1].2.ambi, 1);
    }

    #[test]
    fn happy_summary_matches_the_legacy_dataframe_shape_and_counts() {
        let output = tempfile::NamedTempFile::new().expect("temporary output");
        let header = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,FILTER";
        let rows = vec![
            "0,chr1,10,TP,A,A,C,C,".to_string(),
            "0,chr1,20,FN,,A,,C,".to_string(),
            "0,chr1,30,FP,A,,C,,LowQual".to_string(),
            "0,chr1,40,AMBI,A,,C,,".to_string(),
            "0,chr1,50,UNK,A,,C,,".to_string(),
        ];
        write_happy_style_summary(output.path(), header, &rows, "hcc.strelka.indel")
            .expect("happy summary");
        let text = fs::read_to_string(output.path()).expect("read happy summary");
        let expected = concat!(
            ",Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio\n",
            "0,INDEL,PASS,2,1,1,3,0,2,NA,0.5,1.0,0.6667,0.6667,NA,NA,NA,NA\n",
            "0,INDEL,ALL,2,1,1,4,1,2,NA,0.5,0.5,0.5,0.5,NA,NA,NA,NA\n"
        );
        assert_eq!(text, expected);
    }

    #[test]
    fn metrics_json_preserves_pandas_nullable_numeric_columns() {
        let values = vec!["0".to_string(), String::new(), ".".to_string()];
        assert_eq!(infer_type(&values), "double");
        assert_eq!(
            column_json("recall2", "recall2", "double", &values, false),
            "{\"values\": [0, null, null], \"type\": \"double\", \"id\": \"recall2\", \"label\": \"recall2\"}"
        );

        assert_eq!(infer_type(&["1".to_string(), "2".to_string()]), "int64");
        assert_eq!(infer_type(&["1.5".to_string(), String::new()]), "double");
        assert_eq!(infer_type(&[String::new(), ".".to_string()]), "string");
    }

    #[test]
    fn normalization_left_aligns_trims_and_deduplicates_like_bcftools_norm() {
        let reference = BTreeMap::from([("chr1".to_string(), "AAAAAC".to_string())]);
        let deletion = raw_record("chr1\t3\t.\tAA\tA\t.\tPASS\t.");
        let duplicate = deletion.clone();
        let normalized = normalize_somatic_records(vec![deletion, duplicate], &reference);
        assert_eq!(normalized.len(), 1);
        assert_eq!(normalized[0].pos, 1);
        assert_eq!(normalized[0].ref_allele, "AA");
        assert_eq!(normalized[0].alt_allele, "A");

        let mismatch = raw_record("chr1\t1\t.\tC\tT\t.\tPASS\t.");
        assert!(normalize_somatic_records(vec![mismatch], &reference).is_empty());
    }

    #[test]
    fn chromosome_rewrite_switch_matches_the_legacy_prefix_expression() {
        let prefixed = BTreeSet::from(["chr1".to_string(), "chrM".to_string()]);
        assert_eq!(somatic_chrom("1", &prefixed, true), "chr1");
        assert_eq!(somatic_chrom("MT", &prefixed, true), "chrM");
        assert_eq!(somatic_chrom("1", &prefixed, false), "1");
        assert_eq!(somatic_chrom("GL0001", &prefixed, true), "GL0001");
    }

    #[test]
    fn caller_specific_roc_matches_the_legacy_pandas_csv_shape() {
        let output = tempfile::NamedTempFile::new().expect("temporary ROC");
        let header = ",CHROM,POS,tag,EVS,FILTER,NT";
        let rows = vec![
            "0,chr1,10,TP,10.00000000,,ref".to_string(),
            "0,chr1,20,FP,5.00000000,LowEVS,ref".to_string(),
            "0,chr1,30,FN,,,".to_string(),
        ];
        write_somatic_roc(output.path(), header, &rows, "strelka.snv").expect("write somatic ROC");
        assert_eq!(
            fs::read_to_string(output.path()).expect("read ROC"),
            concat!(
                ",EVS,tp,fp,fn,precision,recall\n",
                "0,0,1,1,1,0.50000000,0.50000000\n",
                "1,5,1,1,1,0.50000000,0.50000000\n",
                "2,10,1,0,1,1.00000000,0.50000000\n"
            )
        );
    }

    #[test]
    fn generic_feature_qual_is_serialized_as_pandas_float() {
        let truth = raw_record("chr1\t1\t.\tA\tC\t60\tPASS\t.");
        let row = render_generic_tp_row(0, &truth, &truth);
        assert!(row.ends_with(",60.00000000,60.00000000"));
        assert_eq!(cpp_default_six(0.872_429_31), "0.872429");
        assert_eq!(cpp_default_six(-1.0), "-1");
    }
}
