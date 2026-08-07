use crate::adapters::metrics_json;
use crate::adapters::report::{self, suffixed_report_path};
use crate::adapters::vcf::{self};
use crate::application::roc_publication;
use crate::cli_compat::cli::QuantifyArgs;
use crate::domain::{AnnotatedRow, CountsBucket, Interval, RawVcfRecord, TypeCounts};
use crate::engines::roc;
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};

mod annotations;
mod counting;
mod regions;
mod reporting;
mod stratification;

use annotations::*;
use counting::*;
use regions::{region_intersection_size, region_size};
use reporting::*;
use stratification::*;

const INDEL_SUBTYPES: [&str; 9] = [
    "C16_PLUS", "C1_5", "C6_15", "D16_PLUS", "D1_5", "D6_15", "I16_PLUS", "I1_5", "I6_15",
];

type RegionMap = BTreeMap<String, Vec<Interval>>;
type RegionLevels = BTreeMap<String, usize>;
type LoadedRegions = (Option<Vec<Interval>>, RegionMap, RegionLevels);

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct BenchmarkSamples {
    truth: Option<usize>,
    query: Option<usize>,
}

impl BenchmarkSamples {
    #[cfg(test)]
    const POSITIONAL: Self = Self {
        truth: Some(0),
        query: Some(1),
    };

    fn has_both(self) -> bool {
        self.truth.is_some() && self.query.is_some()
    }
}

#[derive(Clone, Debug)]
struct ClassifiedVariant {
    variant_type: String,
    subtypes: Vec<String>,
    ti: usize,
    tv: usize,
    het: bool,
    homalt: bool,
    status: String,
    passes_filter: bool,
    subsets: Vec<String>,
    fp_class: Option<&'static str>,
}

#[derive(Clone, Debug, Default)]
struct QuantifyTypeCounts {
    counts: TypeCounts,
    fp_gt: usize,
    fp_al: usize,
}

impl Deref for QuantifyTypeCounts {
    type Target = TypeCounts;

    fn deref(&self) -> &Self::Target {
        &self.counts
    }
}

impl DerefMut for QuantifyTypeCounts {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.counts
    }
}

#[derive(Clone, Debug, Default)]
struct QuantifyCountMaps {
    by_type: BTreeMap<String, QuantifyTypeCounts>,
    by_subset_type: BTreeMap<String, BTreeMap<String, QuantifyTypeCounts>>,
    by_subtype: BTreeMap<String, BTreeMap<String, QuantifyTypeCounts>>,
    by_subset_subtype: BTreeMap<String, BTreeMap<String, BTreeMap<String, QuantifyTypeCounts>>>,
}

struct ExtendedTableOptions<'a> {
    subset_size: usize,
    whole_reference_size: usize,
    stratification_sizes: &'a BTreeMap<String, usize>,
    stratification_confidence_sizes: &'a BTreeMap<String, usize>,
    stratification_levels: &'a RegionLevels,
    confidence_size: Option<usize>,
    qq_field: &'a str,
    ci_alpha: f64,
}

pub(crate) fn run(args: QuantifyArgs) -> Result<()> {
    run_with_metric_indices(args).map(|_| ())
}

/// Run qfy and return the legacy table row indices used by metrics JSON.
///
/// `hap.py` invokes qfy internally and then rewrites the JSON container with
/// hap.py metadata. Returning the indices keeps that rewrite from replacing
/// qfy's unordered-table row identifiers with synthetic `0..N` values.
pub(crate) fn run_with_metric_indices(args: QuantifyArgs) -> Result<roc::MetricIndices> {
    run_with_metric_indices_mode(args, CompareQuantifyMode::default())
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct CompareQuantifyMode {
    pub preserve_missing_query_qq: bool,
    pub preserve_missing_nocall_bd: bool,
    pub inherit_same_position_tp_qq: bool,
    pub roc_value_from_qq: bool,
}

pub(crate) fn run_from_compare(
    args: QuantifyArgs,
    mode: CompareQuantifyMode,
) -> Result<roc::MetricIndices> {
    run_with_metric_indices_mode(args, mode)
}

fn run_with_metric_indices_mode(
    args: QuantifyArgs,
    mode: CompareQuantifyMode,
) -> Result<roc::MetricIndices> {
    validate_options(&args)?;
    if let Some(logfile) = args.logfile.as_deref() {
        OpenOptions::new()
            .create(true)
            .append(true)
            .open(logfile)
            .with_context(|| format!("failed to open qfy logfile {logfile}"))?;
    }
    let annotation_type = args.annotation_type.as_deref().unwrap_or("xcmp");
    let prefix = Path::new(&args.report_prefix);
    let output_vcf = suffixed_report_path(prefix, if args.bcf { "bcf" } else { "vcf.gz" });
    if args.write_vcf && paths_refer_to_same_file(Path::new(&args.input_vcf), &output_vcf) {
        bail!(
            "cannot overwrite input VCF: {} would be overwritten by output {}",
            args.input_vcf,
            output_vcf.display()
        );
    }
    require_quantifier_index(Path::new(&args.input_vcf))?;
    let (mut headers, mut records) = vcf::load_raw_vcf(Path::new(&args.input_vcf))?;
    if annotation_type == "ga4gh" {
        validate_ga4gh_qq_fields(&headers, &records)?;
    }
    let benchmark_samples = benchmark_sample_indices(&headers);
    let do_roc = args.do_roc && benchmark_samples.has_both();
    let reference = crate::adapters::fasta::read_sequences(Path::new(&args.reference))?;
    let reference_contigs = reference.keys().cloned().collect::<BTreeSet<_>>();
    let (confidence, stratifications, stratification_levels) =
        load_regions(&args, &reference_contigs)?;
    for record in &mut records {
        if args.preserve_info {
            let extent = legacy_regions_extent(record);
            set_info_value(&mut record.info, "RegionsExtent", &extent);
        }
        annotate_regions_for_samples(
            record,
            confidence.as_deref(),
            &stratifications,
            annotation_type == "ga4gh",
            benchmark_samples,
            mode.preserve_missing_nocall_bd,
        );
        if annotation_type == "xcmp" {
            reannotate_xcmp_record_for_samples(
                record,
                confidence.is_some(),
                &args.roc,
                benchmark_samples,
            );
        } else {
            reannotate_ga4gh_record(record);
        }
    }
    propagate_superlocus_annotations_for_samples(
        &mut records,
        annotation_type,
        benchmark_samples,
        mode.preserve_missing_query_qq,
        mode.inherit_same_position_tp_qq,
    );
    for record in &mut records {
        decorate_quantified_record_for_samples(
            record,
            annotation_type,
            args.preserve_info,
            args.output_vtc,
            confidence.is_some(),
            benchmark_samples,
        );
        if annotation_type == "ga4gh" {
            normalize_integer_like_format_values(record, "QQ");
        }
    }
    let subset_size = contigs_in_input(&records)
        .into_iter()
        .filter_map(|contig| {
            reference
                .get(&contig)
                .map(|sequence| n_trimmed_length(sequence))
        })
        .sum::<usize>();
    let whole_reference_size = reference
        .values()
        .map(|sequence| n_trimmed_length(sequence))
        .sum::<usize>();

    let mut all_counts = QuantifyCountMaps::default();
    let mut pass_counts = QuantifyCountMaps::default();
    let mut subsets_present = BTreeMap::<String, BTreeSet<String>>::new();
    let stratification_sizes = stratifications
        .iter()
        .map(|(name, intervals)| (name.clone(), region_size(intervals)))
        .collect::<BTreeMap<_, _>>();
    let confidence_size = confidence.as_deref().map(region_size);
    let stratification_confidence_sizes = confidence
        .as_deref()
        .map(|confidence| {
            stratifications
                .iter()
                .map(|(name, intervals)| {
                    (
                        name.clone(),
                        region_intersection_size(intervals, confidence),
                    )
                })
                .collect::<BTreeMap<_, _>>()
        })
        .unwrap_or_default();

    for record in &records {
        let truth = benchmark_samples
            .truth
            .and_then(|sample_index| classify_side(record, sample_index));
        let query = benchmark_samples
            .query
            .and_then(|sample_index| classify_side(record, sample_index));

        if let Some(classified) = truth {
            record_truth(&mut all_counts, &classified);
            if classified.passes_filter {
                record_truth_filtered(&mut pass_counts, &classified);
            } else {
                record_truth_total_only(&mut pass_counts, &classified);
            }
            register_subsets(&mut subsets_present, &classified);
        }

        if let Some(classified) = query {
            record_query(&mut all_counts, &classified);
            if classified.passes_filter {
                record_query(&mut pass_counts, &classified);
            }
            register_subsets(&mut subsets_present, &classified);
        }
    }
    // PASS metrics keep the complete truth denominator. Any truth call that
    // ceased to be a TP because its paired query record was filtered becomes
    // a PASS-level FN. Legacy qfy derives this lane algebraically.
    derive_pass_truth_false_negatives(&mut pass_counts);
    let reported_types = all_counts
        .by_type
        .keys()
        .chain(pass_counts.by_type.keys())
        .cloned()
        .collect::<BTreeSet<_>>();
    for variant_type in reported_types {
        subsets_present
            .entry(variant_type)
            .or_default()
            .extend(stratifications.keys().cloned());
    }

    if let Some(parent) = prefix.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }

    // Legacy qfy always emits the compact summary. `--no-write-counts`
    // disables only the extended count table.
    write_quantify_summary(
        &suffixed_report_path(prefix, "summary.csv"),
        &all_counts.by_type,
        &pass_counts.by_type,
    )?;
    if args.write_counts {
        write_quantify_extended(
            &suffixed_report_path(prefix, "extended.csv"),
            &all_counts,
            &pass_counts,
            &subsets_present,
            &ExtendedTableOptions {
                subset_size,
                whole_reference_size,
                stratification_sizes: &stratification_sizes,
                stratification_confidence_sizes: &stratification_confidence_sizes,
                stratification_levels: &stratification_levels,
                confidence_size,
                qq_field: &args.roc,
                ci_alpha: args.ci_alpha,
            },
        )?;
    }
    if args.write_vcf {
        ensure_info_header(
            &mut headers,
            "Regions",
            "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">",
        );
        if args.preserve_info {
            ensure_info_header(
                &mut headers,
                "RegionsExtent",
                "##INFO=<ID=RegionsExtent,Number=.,Type=String,Description=\"Trimmed reference coordinates matched to regions for this record.\">",
            );
        }
        if annotation_type == "ga4gh" {
            ensure_ga4gh_headers(&mut headers);
            canonicalize_ga4gh_header_order(&mut headers);
        }
        if args.output_vtc {
            ensure_info_header(
                &mut headers,
                "VTC",
                "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">",
            );
            if annotation_type == "xcmp" {
                ensure_info_header(
                    &mut headers,
                    "XCMP",
                    "##INFO=<ID=XCMP,Number=.,Type=String,Description=\"XCMP extra information.\">",
                );
            }
        }
        vcf::write_raw_vcf(&output_vcf, &headers, &records)?;
    }
    let mut rows = records
        .iter()
        .enumerate()
        .filter_map(|(index, record)| {
            Some(AnnotatedRow {
                sort_key: (record.chrom.clone(), record.pos, index, 0),
                record: roc_record(record, benchmark_samples)?,
                query_pass: record.is_pass(),
                fp_class: benchmark_samples
                    .query
                    .and_then(|sample_index| query_fp_class_for_sample(record, sample_index)),
                xcmp_ctype: None,
                xcmp_hap_match: false,
            })
        })
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| left.sort_key.cmp(&right.sort_key));
    let roc_options = roc::RocOptions {
        qq_field: args.roc.clone(),
        score_field: mode.roc_value_from_qq.then(|| "QQ".to_string()),
        ignored_filters: args
            .roc_filter
            .as_deref()
            .into_iter()
            .flat_map(|filters| filters.split(|ch: char| ch.is_whitespace() || ";,".contains(ch)))
            .filter(|filter| !filter.is_empty())
            .map(str::to_string)
            .collect(),
        roc_regions: args.roc_regions.iter().cloned().collect(),
        delta: legacy_roc_delta(args.roc_delta),
        ci_alpha: args.ci_alpha,
        preserve_raw_table: args.verbose,
        output_rocs: do_roc,
        whole_reference_size: Some(whole_reference_size),
        subset_sizes: stratification_sizes.clone(),
        subset_confidence_sizes: stratification_confidence_sizes.clone(),
    };
    let roc_indices = roc_publication::write_roc_files_with_options(
        prefix,
        &rows,
        subset_size,
        confidence_size.unwrap_or(0),
        &roc_options,
    )?;
    apply_stratification_levels(prefix, &stratification_levels, args.verbose)?;
    if !do_roc {
        compact_no_roc_outputs(prefix)?;
    }
    if !args.no_json {
        write_metrics_json(prefix, args.write_counts, &roc_indices)?;
    }
    Ok(roc_indices)
}

fn paths_refer_to_same_file(left: &Path, right: &Path) -> bool {
    if left == right {
        return true;
    }
    match (fs::canonicalize(left), fs::canonicalize(right)) {
        (Ok(left), Ok(right)) => left == right,
        _ => false,
    }
}

fn require_quantifier_index(path: &Path) -> Result<()> {
    if !path.exists() {
        return Ok(());
    }
    let extension = path.extension().and_then(|extension| extension.to_str());
    if !matches!(extension, Some("gz" | "bgz" | "bgzf" | "bcf")) {
        bail!(
            "quantifier input {} must be compressed and indexed",
            path.display()
        );
    }

    let index_suffixes = if extension == Some("bcf") {
        &["csi"][..]
    } else {
        &["tbi", "csi"][..]
    };
    let has_index = index_suffixes
        .iter()
        .map(|suffix| PathBuf::from(format!("{}.{suffix}", path.display())))
        .any(|index| index.is_file());
    if !has_index {
        bail!(
            "quantifier input {} requires a companion .tbi or .csi index",
            path.display()
        );
    }
    Ok(())
}

fn benchmark_sample_indices(headers: &[String]) -> BenchmarkSamples {
    let sample_names = headers
        .iter()
        .rev()
        .find(|header| header.starts_with("#CHROM"))
        .map(|header| header.split('\t').skip(9).collect::<Vec<_>>())
        .unwrap_or_default();
    BenchmarkSamples {
        truth: sample_names.iter().rposition(|name| *name == "TRUTH"),
        query: sample_names.iter().rposition(|name| *name == "QUERY"),
    }
}

fn roc_record(record: &RawVcfRecord, samples: BenchmarkSamples) -> Option<RawVcfRecord> {
    let mut normalized = record.clone();
    normalized.samples = vec![
        record.samples.get(samples.truth?)?.clone(),
        record.samples.get(samples.query?)?.clone(),
    ];
    Some(normalized)
}

fn legacy_roc_delta(delta: f64) -> f64 {
    if delta == 0.0 { 0.1 } else { delta }
}

fn validate_options(args: &QuantifyArgs) -> Result<()> {
    match args.annotation_type.as_deref().unwrap_or("xcmp") {
        "xcmp" | "ga4gh" => {}
        other => bail!("unsupported qfy annotation type: {other}"),
    }
    if args.roc.is_empty() {
        bail!("qfy --roc field cannot be empty");
    }
    if args.roc_regions.iter().any(|region| region.is_empty()) {
        bail!("qfy --roc-regions entries cannot be empty");
    }
    if !args.roc_delta.is_finite() || args.roc_delta < 0.0 {
        bail!("qfy --roc-delta must be finite and nonnegative");
    }
    if !args.ci_alpha.is_finite()
        || (args.ci_alpha != 0.0 && !(0.0 < args.ci_alpha && args.ci_alpha < 1.0))
    {
        bail!("qfy --ci-alpha must be 0 (disabled) or strictly between 0 and 1");
    }
    Ok(())
}

#[cfg(test)]
mod test_suite;
