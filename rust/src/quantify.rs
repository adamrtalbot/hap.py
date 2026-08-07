use crate::cli::QuantifyArgs;
use crate::compare::{AnnotatedRow, TypeCounts, suffixed_report_path};
use crate::metrics_json;
use crate::report::{self, CountsBucket};
use crate::roc;
use crate::vcf::{self, RawVcfRecord};
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::ops::{Deref, DerefMut};
use std::path::{Path, PathBuf};

const INDEL_SUBTYPES: [&str; 9] = [
    "C16_PLUS", "C1_5", "C6_15", "D16_PLUS", "D1_5", "D6_15", "I16_PLUS", "I1_5", "I6_15",
];

type RegionMap = BTreeMap<String, Vec<vcf::BedInterval>>;
type RegionLevels = BTreeMap<String, usize>;
type LoadedRegions = (Option<Vec<vcf::BedInterval>>, RegionMap, RegionLevels);

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

pub fn run(args: QuantifyArgs) -> Result<()> {
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
    let reference = crate::fasta::read_sequences(Path::new(&args.reference))?;
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
                line: roc_record_line(record, benchmark_samples)?,
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
    let roc_indices = roc::write_roc_files_with_options(
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

fn roc_record_line(record: &RawVcfRecord, samples: BenchmarkSamples) -> Option<String> {
    let mut normalized = record.clone();
    normalized.samples = vec![
        record.samples.get(samples.truth?)?.clone(),
        record.samples.get(samples.query?)?.clone(),
    ];
    Some(normalized.to_line())
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

fn load_regions(
    args: &QuantifyArgs,
    reference_contigs: &BTreeSet<String>,
) -> Result<LoadedRegions> {
    let mut confidence = args
        .fp_bedfile
        .as_deref()
        .map(|path| load_region_bed(Path::new(path), reference_contigs, args.strat_fixchr))
        .transpose()
        .with_context(|| "failed to load qfy confidence regions")?;
    if let Some(truth_vcf) = args.adjust_conf_regions.as_deref() {
        let raw_confidence = confidence.as_deref().ok_or_else(|| {
            anyhow::anyhow!("qfy --adjust-conf-regions requires --false-positives")
        })?;
        let (_, truth_records) = vcf::load_raw_vcf(Path::new(truth_vcf)).with_context(|| {
            format!("failed to load --adjust-conf-regions truth VCF {truth_vcf}")
        })?;
        let padding = truth_confidence_padding(&truth_records, raw_confidence);
        confidence.get_or_insert_default().extend(padding);
    }
    let mut paths = BTreeMap::<String, PathBuf>::new();

    if let Some(tsv_path) = args.strat_tsv.as_deref() {
        let tsv_path = Path::new(tsv_path);
        let text = fs::read_to_string(tsv_path)
            .with_context(|| format!("failed to read stratification TSV {}", tsv_path.display()))?;
        for (line_index, line) in text.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let (name, raw_path) = line.split_once('\t').ok_or_else(|| {
                anyhow::anyhow!(
                    "stratification TSV line {} has no region file",
                    line_index + 1
                )
            })?;
            let path = resolve_stratification_path(raw_path.trim(), tsv_path);
            insert_region_path(&mut paths, name.trim(), path)?;
        }
    }

    for spec in &args.strat_regions {
        let (name, path) = spec.split_once(':').ok_or_else(|| {
            anyhow::anyhow!("invalid --stratification-region '{spec}'; expected NAME:BED")
        })?;
        insert_region_path(&mut paths, name.trim(), PathBuf::from(path.trim()))?;
    }

    let mut regions = BTreeMap::<String, Vec<vcf::BedInterval>>::new();
    let mut levels = BTreeMap::new();
    for (raw_name, path) in paths {
        let (name, fixed_label) = dynamic_region_label(&raw_name);
        // QuantifyRegions::load collapses every lane whose label begins
        // with CONF into the reserved confidence lane. hap.py relies on
        // that behavior for its generated CONF_VARS truth padding.
        if name.starts_with("CONF") {
            let intervals = load_region_bed(&path, reference_contigs, args.strat_fixchr)
                .with_context(|| format!("failed to load stratification region {name}"))?;
            confidence.get_or_insert_default().extend(intervals);
        } else {
            let (loaded, loaded_levels) = load_stratification_bed(
                &path,
                &name,
                fixed_label,
                reference_contigs,
                args.strat_fixchr,
            )
            .with_context(|| format!("failed to load stratification region {name}"))?;
            for (label, intervals) in loaded {
                regions.entry(label).or_default().extend(intervals);
            }
            levels.extend(loaded_levels);
        }
    }
    Ok((confidence, regions, levels))
}

fn dynamic_region_label(raw_name: &str) -> (String, bool) {
    if let Some(name) = raw_name.strip_prefix('=') {
        (name.to_string(), true)
    } else {
        (raw_name.to_string(), raw_name.starts_with("CONF"))
    }
}

fn truth_confidence_padding(
    records: &[RawVcfRecord],
    targets: &[vcf::BedInterval],
) -> Vec<vcf::BedInterval> {
    let mut records = records.iter().collect::<Vec<_>>();
    records.sort_by(|left, right| left.chrom.cmp(&right.chrom).then(left.pos.cmp(&right.pos)));
    let mut output = Vec::new();
    let mut active: Option<vcf::BedInterval> = None;
    for record in records {
        if !targets
            .iter()
            .any(|target| target.matches(&record.chrom, record.pos))
        {
            continue;
        }
        let (start, end) = effective_reference_range(record)
            .map(|(start, end, _)| (start.saturating_sub(1), end))
            .unwrap_or((
                record.pos.saturating_sub(1),
                record.pos.saturating_sub(1) + record.ref_allele.len(),
            ));
        let overlaps_active = active.as_ref().is_some_and(|interval| {
            interval.chrom == record.chrom && start < interval.end && end > interval.start
        });
        if overlaps_active {
            if let Some(interval) = active.as_mut() {
                interval.start = interval.start.min(start);
                // gvcf2bed intentionally extends with the next record's start,
                // not its end; qfy relies on this historical under-extension.
                interval.end = interval.end.max(start + 1);
            }
        } else {
            if let Some(interval) = active.take() {
                output.push(interval);
            }
            active = Some(vcf::BedInterval {
                chrom: record.chrom.clone(),
                start,
                end,
            });
        }
    }
    if let Some(interval) = active {
        output.push(interval);
    }
    output
}

fn resolve_stratification_path(raw_path: &str, tsv_path: &Path) -> PathBuf {
    let path = PathBuf::from(raw_path);
    if path.exists() || path.is_absolute() {
        path
    } else {
        tsv_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(path)
    }
}

fn insert_region_path(
    paths: &mut BTreeMap<String, PathBuf>,
    name: &str,
    path: PathBuf,
) -> Result<()> {
    if name.is_empty() || name == "=" {
        bail!("stratification region name cannot be empty");
    }
    if name == "CONF" {
        bail!("stratification region name CONF is reserved for --false-positives");
    }
    if path.as_os_str().is_empty() {
        bail!("no file for stratification region {name}");
    }
    if paths.insert(name.to_string(), path).is_some() {
        bail!("duplicate stratification region ID: {name}");
    }
    Ok(())
}

fn load_region_bed(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    fixchr: bool,
) -> Result<Vec<vcf::BedInterval>> {
    Ok(load_region_bed_rows(path, reference_contigs, fixchr)?
        .into_iter()
        .map(|(interval, _)| interval)
        .collect())
}

fn load_stratification_bed(
    path: &Path,
    parent_label: &str,
    fixed_label: bool,
    reference_contigs: &BTreeSet<String>,
    fixchr: bool,
) -> Result<(RegionMap, RegionLevels)> {
    let mut regions = RegionMap::new();
    let mut levels = RegionLevels::from([(parent_label.to_string(), 0)]);
    for (interval, dynamic_label) in load_region_bed_rows(path, reference_contigs, fixchr)? {
        regions
            .entry(parent_label.to_string())
            .or_default()
            .push(interval.clone());
        if !fixed_label && let Some(dynamic_label) = dynamic_label {
            // Pinned QuantifyRegions adds the complete fourth-column value
            // as one level-1 child. Its apparent trailing-number hierarchy
            // branch starts at string[size] and is therefore unreachable:
            // `coding_1` is EXTRA_coding_1, never EXTRA_coding + child.
            let label = format!("{parent_label}_{dynamic_label}");
            regions.entry(label.clone()).or_default().push(interval);
            levels.insert(label, 1);
        }
    }
    // Empty BEDs still register their parent lane in legacy QuantifyRegions.
    regions.entry(parent_label.to_string()).or_default();
    Ok((regions, levels))
}

fn load_region_bed_rows(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    fixchr: bool,
) -> Result<Vec<(vcf::BedInterval, Option<String>)>> {
    let text =
        vcf::read_text(path).with_context(|| format!("failed to read BED {}", path.display()))?;
    let mut intervals = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 3 {
            bail!(
                "BED line {} has fewer than 3 columns in {}",
                line_index + 1,
                path.display()
            );
        }
        let chrom = if fixchr {
            vcf::normalize_chrom(fields[0], reference_contigs)
        } else {
            fields[0].to_string()
        };
        let start = fields[1]
            .parse::<usize>()
            .with_context(|| format!("invalid BED start '{}' in {}", fields[1], path.display()))?;
        let end = fields[2]
            .parse::<usize>()
            .with_context(|| format!("invalid BED end '{}' in {}", fields[2], path.display()))?;
        if end < start {
            bail!(
                "BED end precedes start on line {} in {}",
                line_index + 1,
                path.display()
            );
        }
        intervals.push((
            vcf::BedInterval { chrom, start, end },
            fields.get(3).map(|label| (*label).to_string()),
        ));
    }
    Ok(intervals)
}

#[cfg(test)]
fn annotate_regions(
    record: &mut RawVcfRecord,
    confidence: Option<&[vcf::BedInterval]>,
    stratifications: &RegionMap,
    rewrite_ga4gh_decisions: bool,
) {
    annotate_regions_for_samples(
        record,
        confidence,
        stratifications,
        rewrite_ga4gh_decisions,
        BenchmarkSamples::POSITIONAL,
        false,
    );
}

fn annotate_regions_for_samples(
    record: &mut RawVcfRecord,
    confidence: Option<&[vcf::BedInterval]>,
    stratifications: &RegionMap,
    rewrite_ga4gh_decisions: bool,
    samples: BenchmarkSamples,
    preserve_missing_nocall_bd: bool,
) {
    let effective_range = effective_reference_range(record);
    let record_chrom = record.chrom.clone();
    let record_pos = record.pos;
    let record_end = record.end_pos();
    let overlaps = |intervals: &[vcf::BedInterval]| match effective_range {
        Some((start, end, pure_insertion)) => {
            if pure_insertion {
                // Insertions occupy the gap between two reference anchors.
                // Legacy only marks them confident when both anchors are
                // covered (possibly by different CONF input lanes).
                (start..=end).all(|pos| {
                    intervals
                        .iter()
                        .any(|interval| interval.matches(&record_chrom, pos))
                })
            } else {
                intervals
                    .iter()
                    .any(|interval| interval.overlaps(&record_chrom, start, end))
            }
        }
        None => intervals
            .iter()
            .any(|interval| interval.overlaps(&record_chrom, record_pos, record_end)),
    };
    let fully_covered = |intervals: &[vcf::BedInterval]| {
        let (start, end) = effective_range
            .map(|(start, end, _)| (start, end))
            .unwrap_or((record_pos, record_end));
        (start..=end).all(|pos| {
            intervals
                .iter()
                .any(|interval| interval.matches(&record_chrom, pos))
        })
    };
    let mut additions = Vec::new();
    let mut crosses_stratification_boundary = false;
    for (name, intervals) in stratifications {
        if overlaps(intervals) {
            additions.push(name.clone());
            crosses_stratification_boundary |= !fully_covered(intervals);
        }
    }
    if crosses_stratification_boundary {
        additions.push("TS_boundary".to_string());
    }

    if let Some(confidence) = confidence {
        remove_region_tag(&mut record.info, "CONF");
        remove_region_tag(&mut record.info, "TS_boundary");
        if overlaps(confidence) {
            additions.insert(0, "CONF".to_string());
        } else if rewrite_ga4gh_decisions {
            // GA4GHQuantify's `count_unk` rule rewrites both samples outside
            // CONF. Truth-side UNK is omitted from truth counts; query-side
            // UNK supplies QUERY.UNK and ROC unknown counts.
            if let Some(truth) = samples.truth
                && !preserve_missing_nocall_decision(record, truth, preserve_missing_nocall_bd)
            {
                replace_existing_decision(record, truth, "UNK");
            }
            if let Some(query) = samples.query
                && !preserve_missing_nocall_decision(record, query, preserve_missing_nocall_bd)
            {
                replace_existing_decision(record, query, "UNK");
            }
        }
    }
    if !additions.is_empty() {
        merge_region_tags(&mut record.info, &additions);
    }
    if has_region(&record.info, "CONF") {
        move_region_to_front(&mut record.info, "CONF");
    }
}

/// `BlockQuantify::count` finalizes region flags per benchmarking superlocus.
/// A block crossing the recomputed CONF boundary marks every member as
/// `TS_boundary`; a wholly confident block marks every member `TS_contained`.
/// Records without a non-negative BS value form singleton blocks.
#[cfg(test)]
fn propagate_superlocus_annotations(records: &mut [RawVcfRecord], annotation_type: &str) {
    propagate_superlocus_annotations_for_samples(
        records,
        annotation_type,
        BenchmarkSamples::POSITIONAL,
        false,
        false,
    );
}

fn propagate_superlocus_annotations_for_samples(
    records: &mut [RawVcfRecord],
    annotation_type: &str,
    samples: BenchmarkSamples,
    preserve_missing_query_qq: bool,
    inherit_same_position_tp_qq: bool,
) {
    let mut start = 0usize;
    while start < records.len() {
        let chrom = records[start].chrom.clone();
        let bs = benchmark_superlocus(&records[start].info);
        let mut end = start + 1;
        if bs.is_some() {
            while end < records.len()
                && records[end].chrom == chrom
                && benchmark_superlocus(&records[end].info) == bs
            {
                end += 1;
            }
        }

        let has_confident = records[start..end]
            .iter()
            .any(|record| has_region(&record.info, "CONF"));
        let has_non_confident = records[start..end].iter().any(|record| {
            !has_region(&record.info, "CONF") || has_region(&record.info, "TS_boundary")
        });
        let region = if has_confident && has_non_confident {
            Some("TS_boundary")
        } else if has_confident {
            Some("TS_contained")
        } else {
            None
        };
        if let Some(region) = region {
            for record in &mut records[start..end] {
                merge_region_tags(&mut record.info, &[region.to_string()]);
            }
        }

        // XCMP has no truth TP in the compatibility lane, so the legacy
        // post-pass clears truth QQ rather than retaining record QUAL.
        if annotation_type == "xcmp" {
            if let Some(truth) = samples.truth {
                for record in &mut records[start..end] {
                    if record.sample_map(truth).get("BD").map(String::as_str) != Some("TP") {
                        set_format_value(record, truth, "QQ", ".");
                    }
                }
            }
        } else {
            propagate_ga4gh_superlocus_for_samples(
                &mut records[start..end],
                samples,
                preserve_missing_query_qq,
                inherit_same_position_tp_qq,
            );
        }
        start = end;
    }
}

fn benchmark_superlocus(info: &str) -> Option<i64> {
    info_value(info, "BS")?
        .split(',')
        .next()?
        .parse::<i64>()
        .ok()
        .filter(|value| *value >= 0)
}

fn replace_existing_decision(record: &mut RawVcfRecord, sample_index: usize, decision: &str) {
    if record.sample_map(sample_index).contains_key("BD") {
        set_format_value(record, sample_index, "BD", decision);
    }
}

fn preserve_missing_nocall_decision(
    record: &RawVcfRecord,
    sample_index: usize,
    enabled: bool,
) -> bool {
    if !enabled {
        return false;
    }
    let sample = record.sample_map(sample_index);
    sample.get("BVT").map(String::as_str) == Some("NOCALL")
        && sample
            .get("BD")
            .is_none_or(|decision| decision.is_empty() || decision == ".")
}

#[derive(Debug, Eq, PartialEq)]
struct Ga4ghAnnotation {
    bi: String,
    bvt: &'static str,
    blt: &'static str,
}

/// Derive the fields that legacy `GA4GHQuantify` computes from each sample's
/// selected genotype alleles. RTG's GA4GH intermediate supplies BD/BK/QQ but
/// does not supply these count-oriented annotations.
fn reannotate_ga4gh_record(record: &mut RawVcfRecord) {
    ensure_format_fields(record, &["BI", "BVT", "BLT", "QQ"]);
    let has_overwide_genotype = (0..record.samples.len()).any(|sample_index| {
        record
            .sample_map(sample_index)
            .get("GT")
            .is_some_and(|gt| gt.split(['/', '|']).count() > 2)
    });
    for sample_index in 0..record.samples.len() {
        let gt = record
            .sample_map(sample_index)
            .get("GT")
            .cloned()
            .unwrap_or_else(|| "./.".to_string());
        if has_overwide_genotype {
            set_format_value(record, sample_index, "BI", ".");
            set_format_value(record, sample_index, "BVT", "UNK");
            set_format_value(record, sample_index, "BLT", "ambi");
            let qq = if gt.split(['/', '|']).count() > 2 {
                "."
            } else {
                "0"
            };
            set_format_value(record, sample_index, "QQ", qq);
            continue;
        }
        let annotation = ga4gh_annotation(record, &gt);
        set_format_value(record, sample_index, "BI", &annotation.bi);
        set_format_value(record, sample_index, "BVT", annotation.bvt);
        set_format_value(record, sample_index, "BLT", annotation.blt);
    }
}

fn ga4gh_annotation(record: &RawVcfRecord, gt: &str) -> Ga4ghAnnotation {
    let alleles = gt
        .split(['/', '|'])
        .filter(|part| !part.is_empty())
        .map(|part| part.parse::<usize>().ok())
        .collect::<Vec<_>>();
    let all_missing = alleles.is_empty() || alleles.iter().all(Option::is_none);
    let homref = !alleles.is_empty() && alleles.iter().all(|allele| *allele == Some(0));
    let blt = ga4gh_location_type(&alleles);

    let alts = record.alt_allele.split(',').collect::<Vec<_>>();
    let mut selected = alleles
        .iter()
        .filter_map(|allele| allele.as_ref().copied())
        .filter(|allele| *allele > 0)
        .collect::<BTreeSet<_>>();
    // VariantStatistics counts a homozygous alternate allele once. Using a
    // set is equivalent for BI/BVT and also de-duplicates repeated indexes in
    // polyploid input.
    let mut type_bits = 0u8;
    let mut extras = BTreeSet::<String>::new();
    let mut invalid_allele = false;
    for allele in std::mem::take(&mut selected) {
        let Some(alt) = alts.get(allele - 1) else {
            invalid_allele = true;
            continue;
        };
        let (bits, allele_extras) = ga4gh_allele_statistics(record, alt);
        type_bits |= bits;
        extras.extend(allele_extras);
    }

    const SNP: u8 = 1;
    const INS: u8 = 2;
    const DEL: u8 = 4;
    let bvt = if all_missing {
        "NOCALL"
    } else if homref {
        "HOMREF"
    } else if invalid_allele {
        "UNK"
    } else if type_bits == SNP {
        "SNP"
    } else if type_bits & (SNP | INS | DEL) != 0 {
        "INDEL"
    } else {
        "UNK"
    };
    Ga4ghAnnotation {
        bi: if extras.is_empty() {
            ".".to_string()
        } else {
            extras.into_iter().collect::<Vec<_>>().join(",")
        },
        bvt,
        blt,
    }
}

fn ga4gh_location_type(alleles: &[Option<usize>]) -> &'static str {
    if alleles.len() > 2 {
        return "ambi";
    }
    if !alleles.is_empty() && alleles.iter().all(|allele| *allele == Some(0)) {
        return "homref";
    }
    if alleles.len() == 2 {
        match (alleles[0], alleles[1]) {
            (Some(0), Some(right)) | (Some(right), Some(0)) if right > 0 => "het",
            (Some(left), Some(right)) if left > 0 && left == right => "homalt",
            (Some(left), Some(right)) if left > 0 && right > 0 => "hetalt",
            (Some(_), None) | (None, Some(_)) => "halfcall",
            _ => "nocall",
        }
    } else if alleles.iter().all(Option::is_none) || alleles.is_empty() {
        "nocall"
    } else if alleles.len() == 1 {
        // GA4GH's legacy VariantStatistics categorises a one-allele call in
        // the same partial-call bucket as `./1`; it does not emit `hemi`.
        "halfcall"
    } else {
        "unknown"
    }
}

/// Return VariantStatistics' low type bits plus its lexically sorted BI set.
fn ga4gh_allele_statistics(record: &RawVcfRecord, alt: &str) -> (u8, BTreeSet<String>) {
    const SNP: u8 = 1;
    const INS: u8 = 2;
    const DEL: u8 = 4;

    if alt.starts_with('<') && alt.ends_with('>') {
        let tag = alt[1..alt.len() - 1].to_ascii_uppercase();
        let size = record.ref_allele.len().max(1);
        let (bits, token) = if tag.starts_with("DEL") {
            (DEL, ga4gh_size_token('d', size))
        } else if tag.starts_with("INS") || tag.starts_with("DUP") {
            (INS, ga4gh_size_token('i', size))
        } else {
            (INS | DEL, ga4gh_size_token('c', size))
        };
        return (bits, BTreeSet::from([token]));
    }

    let primitives =
        crate::align::realign_ref_var(record.pos, record.ref_allele.as_bytes(), alt.as_bytes());
    let mut total_snp = 0usize;
    let mut total_ins = 0usize;
    let mut total_del = 0usize;
    let mut ti = 0usize;
    let mut tv = 0usize;
    for primitive in primitives {
        if primitive.end < primitive.start {
            total_ins += primitive.alt.len();
        } else if primitive.alt.is_empty() {
            total_del += primitive.end - primitive.start + 1;
        } else if primitive.end == primitive.start && primitive.alt.len() == 1 {
            total_snp += 1;
            let ref_offset = primitive.start.saturating_sub(record.pos);
            let ref_base = record.ref_allele.as_bytes().get(ref_offset).copied();
            let alt_base = primitive.alt.as_bytes().first().copied();
            if ref_base
                .zip(alt_base)
                .is_some_and(|(reference, alternate)| {
                    is_transition_pair(reference as char, alternate as char)
                })
            {
                ti += 1;
            } else {
                tv += 1;
            }
        } else {
            // The shared aligner normally decomposes every concrete allele to
            // SNP/INS/DEL primitives. Retain an INDEL classification if a
            // future representation reaches this fallback.
            total_del += primitive.end - primitive.start + 1;
            total_ins += primitive.alt.len();
        }
    }

    let mut extras = BTreeSet::new();
    if ti > 0 {
        extras.insert("ti".to_string());
    }
    if tv > 0 {
        extras.insert("tv".to_string());
    }
    if total_snp == 0 {
        if total_ins > 0 {
            extras.insert(ga4gh_size_token('i', total_ins));
        }
        if total_del > 0 {
            extras.insert(ga4gh_size_token('d', total_del));
        }
    } else if total_ins + total_del > 0 {
        extras.insert(ga4gh_size_token('c', total_ins + total_del));
    }

    let mut bits = 0u8;
    if total_snp > 0 {
        bits |= SNP;
    }
    if total_ins > 0 {
        bits |= INS;
    }
    if total_del > 0 {
        bits |= DEL;
    }
    (bits, extras)
}

fn ga4gh_size_token(prefix: char, size: usize) -> String {
    let bucket = match size {
        0..=5 => "1_5",
        6..=15 => "6_15",
        _ => "16_plus",
    };
    format!("{prefix}{bucket}")
}

fn is_transition_pair(reference: char, alternate: char) -> bool {
    matches!(
        (
            reference.to_ascii_uppercase(),
            alternate.to_ascii_uppercase()
        ),
        ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C')
    )
}

fn ensure_format_fields(record: &mut RawVcfRecord, additions: &[&str]) {
    let mut keys = record
        .format
        .as_deref()
        .filter(|format| !format.is_empty() && *format != ".")
        .map(|format| format.split(':').map(str::to_string).collect::<Vec<_>>())
        .unwrap_or_default();
    for sample in &mut record.samples {
        let mut values = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        values.resize(keys.len(), ".".to_string());
        *sample = values.join(":");
    }
    for addition in additions {
        if keys.iter().any(|key| key == addition) {
            continue;
        }
        keys.push((*addition).to_string());
        for sample in &mut record.samples {
            sample.push(':');
            sample.push('.');
        }
    }
    record.format = (!keys.is_empty()).then(|| keys.join(":"));
}

fn ensure_ga4gh_headers(headers: &mut Vec<String>) {
    const REQUIRED: [(&str, &str); 8] = [
        (
            "INFO=<ID=BS,",
            "##INFO=<ID=BS,Number=.,Type=Integer,Description=\"Benchmarking superlocus ID for these variants.\">",
        ),
        (
            "FORMAT=<ID=GT,",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        ),
        (
            "FORMAT=<ID=BD,",
            "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">",
        ),
        (
            "FORMAT=<ID=BK,",
            "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">",
        ),
        (
            "FORMAT=<ID=BI,",
            "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">",
        ),
        (
            "FORMAT=<ID=QQ,",
            "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">",
        ),
        (
            "FORMAT=<ID=BVT,",
            "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"High-level variant type (SNP|INDEL).\">",
        ),
        (
            "FORMAT=<ID=BLT,",
            "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"High-level location type (het|homref|hetalt|homalt|nocall).\">",
        ),
    ];
    let insertion_index = headers
        .iter()
        .position(|header| header.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    let mut offset = 0usize;
    for (needle, declaration) in REQUIRED {
        if headers.iter().any(|header| header.contains(needle)) {
            continue;
        }
        headers.insert(insertion_index + offset, declaration.to_string());
        offset += 1;
    }
}

fn validate_ga4gh_qq_fields(headers: &[String], records: &[RawVcfRecord]) -> Result<()> {
    let qq_is_string = headers.iter().any(|header| {
        let Some(body) = header.strip_prefix("##FORMAT=<") else {
            return false;
        };
        let fields = body.split(',').collect::<Vec<_>>();
        fields.contains(&"ID=QQ") && fields.contains(&"Type=String")
    });
    if !qq_is_string {
        return Ok(());
    }

    // Legacy's numeric FORMAT reader sizes its result from the encoded BCF
    // string width. A one-character String QQ is accepted as a missing numeric
    // value, while a wider cell yields multiple values and aborts before any
    // report is published (`BCFHelpers.cpp::getFormatFloat`).
    for record in records {
        let Some(qq_index) = record.format_keys().iter().position(|field| *field == "QQ") else {
            continue;
        };
        let encoded_width = record
            .samples
            .iter()
            .filter_map(|sample| sample.split(':').nth(qq_index))
            .map(str::len)
            .max()
            .unwrap_or(0);
        if encoded_width > 1 {
            bail!(
                "too many QQ fields at {}:{}",
                record.chrom,
                record.pos.saturating_sub(1)
            );
        }
    }
    Ok(())
}

fn canonicalize_ga4gh_header_order(headers: &mut Vec<String>) {
    let Some(pass_index) = headers
        .iter()
        .position(|header| header.starts_with("##FILTER=<ID=PASS,"))
    else {
        return;
    };
    let pass = headers.remove(pass_index);
    let insertion = headers
        .iter()
        .position(|header| !header.starts_with("##fileformat="))
        .unwrap_or(headers.len());
    headers.insert(insertion, pass);
}

fn ensure_info_header(headers: &mut Vec<String>, id: &str, declaration: &str) {
    if headers
        .iter()
        .any(|header| header.contains(&format!("INFO=<ID={id},")))
    {
        return;
    }
    let index = headers
        .iter()
        .position(|header| header.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    headers.insert(index, declaration.to_string());
}

fn propagate_ga4gh_superlocus_for_samples(
    records: &mut [RawVcfRecord],
    samples: BenchmarkSamples,
    preserve_missing_query_qq: bool,
    inherit_same_position_tp_qq: bool,
) {
    let mut same_position_tp_qq = BTreeMap::<usize, String>::new();
    if let Some(query_index) = samples.query {
        for record in records.iter() {
            let query = record.sample_map(query_index);
            let Some(score) = (query.get("BD").map(String::as_str) == Some("TP"))
                .then(|| query.get("QQ").cloned())
                .flatten()
                .filter(|score| {
                    score != "." && score.parse::<f64>().is_ok_and(|value| value.is_finite())
                })
            else {
                continue;
            };
            same_position_tp_qq
                .entry(record.pos)
                .and_modify(|current| {
                    if score.parse::<f64>().unwrap() < current.parse::<f64>().unwrap() {
                        *current = score.clone();
                    }
                })
                .or_insert(score);
        }
    }
    let minimum_tp_qq = records
        .iter()
        .filter_map(|record| {
            let query = record.sample_map(samples.query?);
            (query.get("BD").map(String::as_str) == Some("TP"))
                .then(|| query.get("QQ").cloned())
                .flatten()
        })
        .filter(|score| score != "." && score.parse::<f64>().is_ok_and(|value| value.is_finite()))
        .min_by(|left, right| {
            left.parse::<f64>()
                .unwrap()
                .total_cmp(&right.parse::<f64>().unwrap())
        });

    let block_filters = records
        .iter()
        .flat_map(|record| record.filter.split(';'))
        .filter(|filter| !filter.is_empty() && *filter != "." && *filter != "PASS")
        .map(str::to_string)
        .collect::<BTreeSet<_>>();

    for record in records {
        let truth_tp = samples.truth.is_some_and(|truth| {
            record.sample_map(truth).get("BD").map(String::as_str) == Some("TP")
        });
        if let Some(query) = samples.query {
            let query_sample = record.sample_map(query);
            let query_qq = query_sample.get("QQ").cloned();
            if !preserve_missing_query_qq && query_qq.as_deref().is_none_or(|value| value == ".") {
                let inherited = (truth_tp && inherit_same_position_tp_qq)
                    .then(|| same_position_tp_qq.get(&record.pos))
                    .flatten()
                    .map(String::as_str)
                    .unwrap_or("0");
                set_format_value(record, query, "QQ", inherited);
            }
        }
        let query = samples
            .query
            .map(|query| record.sample_map(query))
            .unwrap_or_default();
        let direct_query_qq = (query.get("BD").map(String::as_str) == Some("TP"))
            .then(|| query.get("QQ").cloned())
            .flatten()
            .filter(|score| {
                score != "." && score.parse::<f64>().is_ok_and(|value| value.is_finite())
            });
        let truth_qq = if truth_tp {
            direct_query_qq.as_ref().or(minimum_tp_qq.as_ref())
        } else {
            None
        };
        if let Some(truth) = samples.truth {
            set_format_value(
                record,
                truth,
                "QQ",
                truth_qq.map(String::as_str).unwrap_or("."),
            );
        }

        if truth_tp
            && query.get("BVT").map(String::as_str) == Some("NOCALL")
            && !block_filters.is_empty()
        {
            let mut filters = record
                .filter
                .split(';')
                .filter(|filter| !filter.is_empty() && *filter != "." && *filter != "PASS")
                .map(str::to_string)
                .collect::<BTreeSet<_>>();
            filters.extend(block_filters.iter().cloned());
            record.filter = filters.into_iter().collect::<Vec<_>>().join(";");
        }
    }
}

fn normalize_integer_like_format_values(record: &mut RawVcfRecord, key: &str) {
    let Some(index) = record.format_keys().iter().position(|field| *field == key) else {
        return;
    };
    for sample in &mut record.samples {
        let mut fields = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        let Some(value) = fields.get_mut(index) else {
            continue;
        };
        if let Ok(number) = value.parse::<f64>()
            && number.is_finite()
            && number.fract() == 0.0
        {
            *value = number.to_string();
            *sample = fields.join(":");
        }
    }
}

/// Legacy `fastainfo` reports the contig length after trimming only terminal
/// N-runs. Internal N-runs remain part of `Subset.Size`.
fn n_trimmed_length(sequence: &str) -> usize {
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
    bytes.len().saturating_sub(leading + trailing)
}

fn legacy_regions_extent(record: &RawVcfRecord) -> String {
    effective_reference_range(record)
        .map(|(start, end, _)| format!("{start}-{end}"))
        .unwrap_or_else(|| format!("{}-{}", record.pos, record.end_pos()))
}

/// Port `QuantifyRegions::annotate`'s per-allele reference span. Coordinates
/// are 1-based inclusive; pure insertions bracket both adjacent reference
/// bases, while the pinned legacy qfy assigns the insertion using its left
/// VCF anchor.
fn effective_reference_range(record: &RawVcfRecord) -> Option<(usize, usize, bool)> {
    let ref_bytes = record.ref_allele.as_bytes();
    let pos_0b = record.pos.saturating_sub(1) as i64;
    let mut updated_start = i64::MAX;
    let mut updated_end = i64::MIN;
    let mut pure_insertion = false;
    let mut has_nucleotide_alt = false;

    for alt in record.alt_allele.split(',') {
        // classifyAlleleString maps a missing ALT to an empty allele and
        // processes it as a reference-consuming deletion. Any other
        // non-nucleotide allele terminates the loop, so later nucleotide
        // ALTs must not expand the range.
        let normalized_alt = if alt.is_empty() || alt == "." {
            String::new()
        } else if alt.bytes().all(is_legacy_nucleotide_base) {
            alt.to_ascii_uppercase()
        } else {
            break;
        };
        let alt_bytes = normalized_alt.as_bytes();
        let mut ref_len = ref_bytes.len();
        let mut alt_len = alt_bytes.len();
        while ref_len > 0 && alt_len > 0 && ref_bytes[ref_len - 1] == alt_bytes[alt_len - 1] {
            ref_len -= 1;
            alt_len -= 1;
        }
        let mut prefix = 0usize;
        while prefix < ref_len && prefix < alt_len && ref_bytes[prefix] == alt_bytes[prefix] {
            prefix += 1;
        }
        let start = pos_0b + prefix as i64;
        let end = pos_0b + ref_len as i64 - 1;
        if !has_nucleotide_alt {
            pure_insertion = true;
        }
        has_nucleotide_alt = true;
        if end >= start {
            updated_start = updated_start.min(start);
            updated_end = updated_end.max(end);
            pure_insertion = false;
        } else {
            updated_start = updated_start.min(start - 1);
            updated_end = updated_end.max(start);
        }
    }

    if !has_nucleotide_alt {
        return None;
    }
    Some((
        (updated_start + 1) as usize,
        (updated_end + 1) as usize,
        pure_insertion,
    ))
}

fn is_legacy_nucleotide_base(base: u8) -> bool {
    matches!(
        base.to_ascii_uppercase(),
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'R'
            | b'Y'
            | b'K'
            | b'M'
            | b'S'
            | b'W'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
            | b'X'
    )
}

fn remove_region_tag(info: &mut String, unwanted: &str) {
    let mut fields = Vec::new();
    for field in info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
    {
        if let Some(regions) = field.strip_prefix("Regions=") {
            let retained = regions
                .split(',')
                .filter(|region| !region.is_empty() && *region != unwanted)
                .collect::<Vec<_>>();
            if !retained.is_empty() {
                fields.push(format!("Regions={}", retained.join(",")));
            }
        } else {
            fields.push(field.to_string());
        }
    }
    *info = if fields.is_empty() {
        ".".to_string()
    } else {
        fields.join(";")
    };
}

fn move_region_to_front(info: &mut String, wanted: &str) {
    let mut fields = Vec::new();
    for field in info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
    {
        if let Some(regions) = field.strip_prefix("Regions=") {
            let mut reordered = regions
                .split(',')
                .filter(|region| !region.is_empty() && *region != wanted)
                .map(str::to_string)
                .collect::<Vec<_>>();
            reordered.insert(0, wanted.to_string());
            fields.push(format!("Regions={}", reordered.join(",")));
        } else {
            fields.push(field.to_string());
        }
    }
    *info = if fields.is_empty() {
        ".".to_string()
    } else {
        fields.join(";")
    };
}

/// Reproduce the old XCMP quantifier rather than consuming GA4GH `BD` fields.
/// XCMP decisions live in record-level `INFO/type`, `kind`, and `ctype`; a
/// finalized hap.py VCF intentionally lacks them, which is why the legacy qfy
/// lane ignores calls inside CONF and labels calls outside CONF as UNK.
#[cfg(test)]
fn reannotate_xcmp_record(
    record: &mut RawVcfRecord,
    has_confidence_regions: bool,
    roc_field: &str,
) {
    reannotate_xcmp_record_for_samples(
        record,
        has_confidence_regions,
        roc_field,
        BenchmarkSamples::POSITIONAL,
    );
}

fn reannotate_xcmp_record_for_samples(
    record: &mut RawVcfRecord,
    has_confidence_regions: bool,
    roc_field: &str,
    samples: BenchmarkSamples,
) {
    let mut decision = info_value(&record.info, "type").unwrap_or_default();
    let mismatch_kind = info_value(&record.info, "kind").unwrap_or_default();
    let comparison_type = info_value(&record.info, "ctype").unwrap_or_default();
    let hap_match = info_flag(&record.info, "HapMatch");
    let import_fail = info_flag(&record.info, "IMPORT_FAIL");
    let query_filtered = info_flag(&record.info, "Q_FILTERED");

    if hap_match && decision != "TP" && !query_filtered {
        decision = "TP".to_string();
    }
    if has_confidence_regions && !has_region(&record.info, "CONF") {
        decision = "UNK".to_string();
    }
    if import_fail {
        decision = "N".to_string();
    }

    let match_kind = if decision == "TP" {
        "gm"
    } else if mismatch_kind == "gtmismatch" {
        "am"
    } else if mismatch_kind == "almismatch" || comparison_type == "hap:mismatch" {
        "lm"
    } else {
        "."
    };

    let record_qual = record.qual.clone();
    for sample_index in 0..record.samples.len() {
        let fields = record.sample_map(sample_index);
        let selected_roc = if roc_field == "QUAL" {
            Some(record_qual.clone())
        } else if let Some(value) = info_value(&record.info, roc_field) {
            Some(value.split(',').next().unwrap_or(".").to_string())
        } else {
            fields.get(roc_field).cloned()
        }
        .unwrap_or_else(|| ".".to_string());
        let gt = fields.get("GT").map(String::as_str).unwrap_or("./.");
        let no_call = gt
            .split(['/', '|'])
            .all(|allele| allele.is_empty() || allele == ".");
        let suppressed =
            import_fail || no_call || (samples.query == Some(sample_index) && query_filtered);
        let sample_decision = if import_fail {
            "N"
        } else if suppressed || decision.is_empty() {
            "."
        } else if samples.truth == Some(sample_index) && decision == "FP" {
            "FN"
        } else {
            decision.as_str()
        };
        set_format_value(record, sample_index, "BD", sample_decision);
        set_format_value(
            record,
            sample_index,
            "BK",
            if suppressed { "." } else { match_kind },
        );
        set_format_value(
            record,
            sample_index,
            "QQ",
            if suppressed {
                "0"
            } else {
                selected_roc.as_str()
            },
        );
    }
}

#[cfg(test)]
fn decorate_quantified_record(
    record: &mut RawVcfRecord,
    annotation_type: &str,
    preserve_info: bool,
    output_vtc: bool,
    has_confidence_regions: bool,
) {
    decorate_quantified_record_for_samples(
        record,
        annotation_type,
        preserve_info,
        output_vtc,
        has_confidence_regions,
        BenchmarkSamples::POSITIONAL,
    );
}

fn decorate_quantified_record_for_samples(
    record: &mut RawVcfRecord,
    annotation_type: &str,
    preserve_info: bool,
    output_vtc: bool,
    has_confidence_regions: bool,
    samples: BenchmarkSamples,
) {
    if output_vtc {
        if annotation_type == "xcmp" {
            let mut decision = info_value(&record.info, "type").unwrap_or_default();
            let mut kind = info_value(&record.info, "kind").unwrap_or_default();
            let ctype = info_value(&record.info, "ctype").unwrap_or_default();
            let query_filtered = info_flag(&record.info, "Q_FILTERED");
            if info_flag(&record.info, "HapMatch") && decision != "TP" && !query_filtered {
                kind = format!("hapmatch__{decision}__{kind}");
                decision = "TP".to_string();
            }
            if has_confidence_regions && !has_region(&record.info, "CONF") {
                decision = "UNK".to_string();
            }
            if info_flag(&record.info, "IMPORT_FAIL") {
                decision = "N".to_string();
                kind = "error".to_string();
            }
            let gtt1 = info_value(&record.info, "gtt1").unwrap_or_else(|| ".".to_string());
            let gtt2 = info_value(&record.info, "gtt2").unwrap_or_else(|| ".".to_string());
            set_info_value(
                &mut record.info,
                "XCMP",
                &format!("{decision}:{kind}:{gtt1}:{gtt2}:{ctype}"),
            );
        }
        let truth = samples
            .truth
            .map(|sample_index| record.sample_map(sample_index))
            .unwrap_or_default();
        let query = samples
            .query
            .map(|sample_index| record.sample_map(sample_index))
            .unwrap_or_default();
        let vtc = legacy_vtc(record, &truth, &query);
        if !vtc.is_empty() {
            set_info_value(&mut record.info, "VTC", &vtc);
        }
        if info_value(&record.info, "Regions").as_deref() == Some("TS_boundary") {
            move_info_field_to_end(&mut record.info, "Regions");
        }
    }

    if annotation_type == "xcmp" && !preserve_info {
        retain_info_fields(
            &mut record.info,
            &["END", "VTC", "Regions", "BS", "XCMP", "IMPORT_FAIL"],
        );
    }
}

fn set_info_value(info: &mut String, key: &str, value: &str) {
    let mut fields = info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .filter(|field| field.split_once('=').map_or(*field, |(name, _)| name) != key)
        .map(str::to_string)
        .collect::<Vec<_>>();
    fields.push(format!("{key}={value}"));
    *info = fields.join(";");
}

fn retain_info_fields(info: &mut String, retained: &[&str]) {
    let fields = info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .filter(|field| {
            let key = field.split_once('=').map_or(*field, |(key, _)| key);
            retained.contains(&key)
        })
        .collect::<Vec<_>>();
    *info = if fields.is_empty() {
        ".".to_string()
    } else {
        fields.join(";")
    };
}

fn move_info_field_to_end(info: &mut String, wanted: &str) {
    let mut fields = info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .map(str::to_string)
        .collect::<Vec<_>>();
    let Some(index) = fields
        .iter()
        .position(|field| field.split_once('=').map_or(field.as_str(), |(key, _)| key) == wanted)
    else {
        return;
    };
    let field = fields.remove(index);
    fields.push(field);
    *info = fields.join(";");
}

fn legacy_vtc(
    record: &RawVcfRecord,
    truth: &BTreeMap<String, String>,
    query: &BTreeMap<String, String>,
) -> String {
    let mut types = BTreeMap::<u8, String>::new();
    for sample in [truth, query] {
        if sample
            .get("BVT")
            .is_none_or(|value| matches!(value.as_str(), "" | "." | "NOCALL"))
        {
            types.insert(0x80, "nocall__nc".to_string());
            continue;
        }
        let alleles = sample
            .get("GT")
            .map(String::as_str)
            .unwrap_or(".")
            .split(['/', '|'])
            .filter_map(|allele| allele.parse::<usize>().ok())
            .collect::<Vec<_>>();
        let mut combined = 0u8;
        for allele in alleles.iter().copied().filter(|allele| *allele > 0) {
            let Some(alternate) = record.alt_allele.split(',').nth(allele - 1) else {
                continue;
            };
            let bits = allele_edit_bits(&record.ref_allele, alternate);
            combined |= bits;
            for bit in [1u8, 2, 4] {
                if bits & bit != 0 {
                    types.insert(bit, format!("nuc__{}", legacy_type_bits(bit)));
                }
            }
            if bits != 0 {
                types.insert(0x10 | bits, format!("al__{}", legacy_type_bits(bits)));
            }
        }
        if combined == 0 {
            continue;
        }
        let location = match sample.get("BLT").map(String::as_str).unwrap_or("") {
            "het" => 0x30,
            "hetalt" => 0x40,
            "hemi" => 0x50,
            "homalt" => 0x90,
            _ => 0xa0,
        };
        let ref_bit = u8::from(alleles.contains(&0)) * 8;
        types.insert(
            location | ref_bit | combined,
            format!(
                "{}__{}",
                sample.get("BLT").map(String::as_str).unwrap_or("unknown"),
                legacy_type_bits(ref_bit | combined)
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

fn info_value(info: &str, key: &str) -> Option<String> {
    info.split(';').find_map(|field| {
        field
            .split_once('=')
            .filter(|(field_key, _)| *field_key == key)
            .map(|(_, value)| value.to_string())
    })
}

fn info_flag(info: &str, key: &str) -> bool {
    info.split(';').any(|field| field == key)
}

fn has_region(info: &str, wanted: &str) -> bool {
    info.split(';')
        .find_map(|field| field.strip_prefix("Regions="))
        .is_some_and(|regions| regions.split(',').any(|region| region == wanted))
}

fn set_format_value(record: &mut RawVcfRecord, sample_index: usize, key: &str, value: &str) {
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

fn merge_region_tags(info: &mut String, additions: &[String]) {
    let mut tags = Vec::<String>::new();
    let mut fields = Vec::<String>::new();
    for field in info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
    {
        if let Some(regions) = field.strip_prefix("Regions=") {
            tags.extend(
                regions
                    .split(',')
                    .filter(|tag| !tag.is_empty())
                    .map(str::to_string),
            );
        } else {
            fields.push(field.to_string());
        }
    }
    for addition in additions {
        if !tags.contains(addition) {
            tags.push(addition.clone());
        }
    }
    let regions = format!("Regions={}", tags.join(","));
    let insertion = fields
        .iter()
        .position(|field| field.starts_with("RegionsExtent="))
        .unwrap_or(fields.len());
    fields.insert(insertion, regions);
    *info = fields.join(";");
}

fn region_size(intervals: &[vcf::BedInterval]) -> usize {
    intervals
        .iter()
        .map(|interval| interval.end.saturating_sub(interval.start))
        .sum()
}

fn merged_regions(intervals: &[vcf::BedInterval]) -> Vec<vcf::BedInterval> {
    let mut sorted = intervals.to_vec();
    sorted.sort_by(|left, right| {
        left.chrom
            .cmp(&right.chrom)
            .then(left.start.cmp(&right.start))
            .then(left.end.cmp(&right.end))
    });
    let mut merged: Vec<vcf::BedInterval> = Vec::with_capacity(sorted.len());
    for interval in sorted {
        if let Some(last) = merged.last_mut()
            && last.chrom == interval.chrom
            && interval.start <= last.end
        {
            last.end = last.end.max(interval.end);
        } else {
            merged.push(interval);
        }
    }
    merged
}

fn region_intersection_size(left: &[vcf::BedInterval], right: &[vcf::BedInterval]) -> usize {
    let left = merged_regions(left);
    let right = merged_regions(right);
    let (mut left_index, mut right_index, mut size) = (0, 0, 0usize);
    while left_index < left.len() && right_index < right.len() {
        let left_interval = &left[left_index];
        let right_interval = &right[right_index];
        match left_interval.chrom.cmp(&right_interval.chrom) {
            std::cmp::Ordering::Less => {
                left_index += 1;
                continue;
            }
            std::cmp::Ordering::Greater => {
                right_index += 1;
                continue;
            }
            std::cmp::Ordering::Equal => {}
        }
        size += left_interval
            .end
            .min(right_interval.end)
            .saturating_sub(left_interval.start.max(right_interval.start));
        if left_interval.end <= right_interval.end {
            left_index += 1;
        } else {
            right_index += 1;
        }
    }
    size
}

fn compact_no_roc_outputs(prefix: &Path) -> Result<()> {
    let all_path = suffixed_report_path(prefix, "roc.all.csv.gz");
    let text = vcf::read_text(&all_path)?;
    let rows = text
        .lines()
        .enumerate()
        .filter(|(index, line)| *index == 0 || line.split(',').nth(6).is_some_and(|qq| qq == "*"))
        .map(|(_, line)| line)
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

fn apply_stratification_levels(
    prefix: &Path,
    levels: &RegionLevels,
    preserve_raw_table: bool,
) -> Result<()> {
    if !levels.values().any(|level| *level > 0) {
        return Ok(());
    }

    let csv_path = suffixed_report_path(prefix, "roc.all.csv.gz");
    let csv = vcf::read_text(&csv_path)?;
    let csv = rewrite_subset_levels(&csv, ',', levels);
    let csv_file = fs::File::create(&csv_path)
        .with_context(|| format!("failed to create {}", csv_path.display()))?;
    let mut encoder = GzEncoder::new(csv_file, Compression::default());
    encoder.write_all(csv.as_bytes())?;
    encoder.finish()?;

    if preserve_raw_table {
        let raw_path = suffixed_report_path(prefix, "roc.tsv");
        if raw_path.exists() {
            let raw = fs::read_to_string(&raw_path)
                .with_context(|| format!("failed to read {}", raw_path.display()))?;
            fs::write(&raw_path, rewrite_subset_levels(&raw, '\t', levels))
                .with_context(|| format!("failed to write {}", raw_path.display()))?;
        }
    }
    Ok(())
}

fn rewrite_subset_levels(text: &str, delimiter: char, levels: &RegionLevels) -> String {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return String::new();
    };
    let header_fields = header.split(delimiter).collect::<Vec<_>>();
    let Some(subset_index) = header_fields.iter().position(|field| *field == "Subset") else {
        return text.to_string();
    };
    let Some(level_index) = header_fields
        .iter()
        .position(|field| *field == "Subset.Level")
    else {
        return text.to_string();
    };

    let separator = delimiter.to_string();
    let mut output = vec![header.to_string()];
    for line in lines {
        let mut fields = line
            .split(delimiter)
            .map(str::to_string)
            .collect::<Vec<_>>();
        if let Some(level) = fields
            .get(subset_index)
            .and_then(|subset| levels.get(subset))
            && let Some(field) = fields.get_mut(level_index)
        {
            *field = format!("{:.6}", *level as f64);
        }
        output.push(fields.join(&separator));
    }
    output.join("\n") + "\n"
}

fn contigs_in_input(records: &[RawVcfRecord]) -> BTreeSet<String> {
    records.iter().map(|record| record.chrom.clone()).collect()
}

fn classify_side(record: &RawVcfRecord, sample_index: usize) -> Option<ClassifiedVariant> {
    let fields = record.sample_map(sample_index);
    let bd = fields.get("BD")?.to_string();
    if bd == "." || bd == "N" {
        return None;
    }
    // XCMP/GA4GH quantification classifies the active genotype alleles, not
    // the record-wide REF/ALT tuple.  This matters for mixed multi-allelic
    // records such as REF=T ALT=C,TATC where the selected allele can be a SNP
    // even though another unselected allele is an insertion.  The quantifier
    // stores that genotype-aware result in BVT/BI/BLT before accumulating.
    let variant_type = fields.get("BVT")?.to_string();
    if !matches!(variant_type.as_str(), "SNP" | "INDEL") {
        return None;
    }
    let info_tokens = fields
        .get("BI")
        .map(String::as_str)
        .unwrap_or(".")
        .split(',')
        .collect::<BTreeSet<_>>();
    let subtypes = if variant_type == "INDEL" {
        info_tokens
            .iter()
            .filter_map(|token| {
                let subtype = token.to_ascii_uppercase();
                INDEL_SUBTYPES
                    .contains(&subtype.as_str())
                    .then_some(subtype)
            })
            .collect()
    } else {
        Vec::new()
    };
    let location_type = fields.get("BLT").map(String::as_str).unwrap_or(".");
    let subsets = parse_subsets(&record.info);
    let fp_class = fp_class(&bd, fields.get("BK").map(String::as_str));
    Some(ClassifiedVariant {
        variant_type,
        subtypes,
        ti: usize::from(info_tokens.contains("ti")),
        tv: usize::from(info_tokens.contains("tv")),
        het: location_type == "het",
        homalt: location_type == "homalt",
        status: bd,
        passes_filter: record.filter == "PASS" || record.filter == ".",
        subsets,
        fp_class,
    })
}

fn fp_class(decision: &str, match_kind: Option<&str>) -> Option<&'static str> {
    if decision != "FP" {
        return None;
    }
    match match_kind {
        Some("am") => Some("gt"),
        Some("lm") => Some("al"),
        _ => None,
    }
}

#[cfg(test)]
fn query_fp_class(record: &RawVcfRecord) -> Option<&'static str> {
    query_fp_class_for_sample(record, 1)
}

fn query_fp_class_for_sample(record: &RawVcfRecord, sample_index: usize) -> Option<&'static str> {
    let fields = record.sample_map(sample_index);
    fp_class(
        fields.get("BD").map(String::as_str).unwrap_or("."),
        fields.get("BK").map(String::as_str),
    )
}

fn parse_subsets(info: &str) -> Vec<String> {
    info.split(';')
        .find_map(|entry| entry.strip_prefix("Regions="))
        .map(|entry| {
            entry
                .split(',')
                .filter(|subset| !subset.is_empty() && *subset != "CONF")
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default()
}

fn register_subsets(
    subsets_present: &mut BTreeMap<String, BTreeSet<String>>,
    classified: &ClassifiedVariant,
) {
    for subset in &classified.subsets {
        subsets_present
            .entry(classified.variant_type.clone())
            .or_default()
            .insert(subset.clone());
    }
}

fn record_truth(counts: &mut QuantifyCountMaps, classified: &ClassifiedVariant) {
    if !matches!(classified.status.as_str(), "TP" | "FN") {
        return;
    }
    add_variant_stats(
        &mut counts
            .by_type
            .entry(classified.variant_type.clone())
            .or_default()
            .truth_total,
        classified,
    );
    for subtype in &classified.subtypes {
        add_variant_stats(
            &mut counts
                .by_subtype
                .entry(classified.variant_type.clone())
                .or_default()
                .entry(subtype.clone())
                .or_default()
                .truth_total,
            classified,
        );
    }
    for subset in &classified.subsets {
        add_variant_stats(
            &mut counts
                .by_subset_type
                .entry(subset.clone())
                .or_default()
                .entry(classified.variant_type.clone())
                .or_default()
                .truth_total,
            classified,
        );
        for subtype in &classified.subtypes {
            add_variant_stats(
                &mut counts
                    .by_subset_subtype
                    .entry(subset.clone())
                    .or_default()
                    .entry(classified.variant_type.clone())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default()
                    .truth_total,
                classified,
            );
        }
    }

    add_variant_stats(
        truth_bucket(
            counts
                .by_type
                .entry(classified.variant_type.clone())
                .or_default(),
            &classified.status,
        ),
        classified,
    );
    for subtype in &classified.subtypes {
        add_variant_stats(
            truth_bucket(
                counts
                    .by_subtype
                    .entry(classified.variant_type.clone())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default(),
                &classified.status,
            ),
            classified,
        );
    }
    for subset in &classified.subsets {
        add_variant_stats(
            truth_bucket(
                counts
                    .by_subset_type
                    .entry(subset.clone())
                    .or_default()
                    .entry(classified.variant_type.clone())
                    .or_default(),
                &classified.status,
            ),
            classified,
        );
        for subtype in &classified.subtypes {
            add_variant_stats(
                truth_bucket(
                    counts
                        .by_subset_subtype
                        .entry(subset.clone())
                        .or_default()
                        .entry(classified.variant_type.clone())
                        .or_default()
                        .entry(subtype.clone())
                        .or_default(),
                    &classified.status,
                ),
                classified,
            );
        }
    }
}

fn record_truth_total_only(counts: &mut QuantifyCountMaps, classified: &ClassifiedVariant) {
    if !matches!(classified.status.as_str(), "TP" | "FN") {
        return;
    }
    add_variant_stats(
        &mut counts
            .by_type
            .entry(classified.variant_type.clone())
            .or_default()
            .truth_total,
        classified,
    );
    for subtype in &classified.subtypes {
        add_variant_stats(
            &mut counts
                .by_subtype
                .entry(classified.variant_type.clone())
                .or_default()
                .entry(subtype.clone())
                .or_default()
                .truth_total,
            classified,
        );
    }
    for subset in &classified.subsets {
        add_variant_stats(
            &mut counts
                .by_subset_type
                .entry(subset.clone())
                .or_default()
                .entry(classified.variant_type.clone())
                .or_default()
                .truth_total,
            classified,
        );
        for subtype in &classified.subtypes {
            add_variant_stats(
                &mut counts
                    .by_subset_subtype
                    .entry(subset.clone())
                    .or_default()
                    .entry(classified.variant_type.clone())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default()
                    .truth_total,
                classified,
            );
        }
    }
}

fn record_truth_filtered(counts: &mut QuantifyCountMaps, classified: &ClassifiedVariant) {
    record_truth_total_only(counts, classified);
    if !matches!(classified.status.as_str(), "TP" | "FN") {
        return;
    }
    add_variant_stats(
        truth_bucket(
            counts
                .by_type
                .entry(classified.variant_type.clone())
                .or_default(),
            &classified.status,
        ),
        classified,
    );
    for subtype in &classified.subtypes {
        add_variant_stats(
            truth_bucket(
                counts
                    .by_subtype
                    .entry(classified.variant_type.clone())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default(),
                &classified.status,
            ),
            classified,
        );
    }
    for subset in &classified.subsets {
        add_variant_stats(
            truth_bucket(
                counts
                    .by_subset_type
                    .entry(subset.clone())
                    .or_default()
                    .entry(classified.variant_type.clone())
                    .or_default(),
                &classified.status,
            ),
            classified,
        );
        for subtype in &classified.subtypes {
            add_variant_stats(
                truth_bucket(
                    counts
                        .by_subset_subtype
                        .entry(subset.clone())
                        .or_default()
                        .entry(classified.variant_type.clone())
                        .or_default()
                        .entry(subtype.clone())
                        .or_default(),
                    &classified.status,
                ),
                classified,
            );
        }
    }
}

fn record_query(counts: &mut QuantifyCountMaps, classified: &ClassifiedVariant) {
    if !matches!(classified.status.as_str(), "TP" | "FP" | "UNK" | "AMBI") {
        return;
    }
    record_query_stats(
        counts
            .by_type
            .entry(classified.variant_type.clone())
            .or_default(),
        classified,
    );
    for subtype in &classified.subtypes {
        record_query_stats(
            counts
                .by_subtype
                .entry(classified.variant_type.clone())
                .or_default()
                .entry(subtype.clone())
                .or_default(),
            classified,
        );
    }
    for subset in &classified.subsets {
        record_query_stats(
            counts
                .by_subset_type
                .entry(subset.clone())
                .or_default()
                .entry(classified.variant_type.clone())
                .or_default(),
            classified,
        );
        for subtype in &classified.subtypes {
            record_query_stats(
                counts
                    .by_subset_subtype
                    .entry(subset.clone())
                    .or_default()
                    .entry(classified.variant_type.clone())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default(),
                classified,
            );
        }
    }
}

fn record_query_stats(stats: &mut QuantifyTypeCounts, classified: &ClassifiedVariant) {
    add_variant_stats(&mut stats.query_total, classified);
    add_variant_stats(query_bucket(stats, &classified.status), classified);
    match classified.fp_class {
        Some("gt") => stats.fp_gt += 1,
        Some("al") => stats.fp_al += 1,
        _ => {}
    }
}

fn truth_bucket<'a>(stats: &'a mut TypeCounts, status: &str) -> &'a mut CountsBucket {
    match status {
        "TP" => &mut stats.truth_tp,
        "FN" => &mut stats.truth_fn,
        _ => &mut stats.truth_total,
    }
}

fn query_bucket<'a>(stats: &'a mut TypeCounts, status: &str) -> &'a mut CountsBucket {
    match status {
        "TP" => &mut stats.query_tp,
        "FP" => &mut stats.query_fp,
        "UNK" | "AMBI" => &mut stats.query_unk,
        _ => &mut stats.query_total,
    }
}

fn add_variant_stats(bucket: &mut CountsBucket, classified: &ClassifiedVariant) {
    bucket.total += 1;
    bucket.ti += classified.ti;
    bucket.tv += classified.tv;
    if classified.het {
        bucket.het += 1;
    }
    if classified.homalt {
        bucket.homalt += 1;
    }
}

fn derive_pass_truth_false_negatives(counts: &mut QuantifyCountMaps) {
    let derive = |stats: &mut TypeCounts| {
        stats.truth_fn = subtract_bucket(&stats.truth_total, &stats.truth_tp);
    };
    for stats in counts.by_type.values_mut() {
        derive(stats);
    }
    for by_subtype in counts.by_subtype.values_mut() {
        for stats in by_subtype.values_mut() {
            derive(stats);
        }
    }
    for by_type in counts.by_subset_type.values_mut() {
        for stats in by_type.values_mut() {
            derive(stats);
        }
    }
    for by_type in counts.by_subset_subtype.values_mut() {
        for by_subtype in by_type.values_mut() {
            for stats in by_subtype.values_mut() {
                derive(stats);
            }
        }
    }
}

fn subtract_bucket(total: &CountsBucket, matched: &CountsBucket) -> CountsBucket {
    CountsBucket {
        total: total.total.saturating_sub(matched.total),
        ti: total.ti.saturating_sub(matched.ti),
        tv: total.tv.saturating_sub(matched.tv),
        het: total.het.saturating_sub(matched.het),
        homalt: total.homalt.saturating_sub(matched.homalt),
    }
}

fn write_quantify_summary(
    path: &Path,
    all_counts: &BTreeMap<String, QuantifyTypeCounts>,
    pass_counts: &BTreeMap<String, QuantifyTypeCounts>,
) -> Result<()> {
    let mut writer = BufWriter::new(
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(
        writer,
        "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio"
    )?;
    if all_counts.is_empty() && pass_counts.is_empty() {
        for variant_type in ["INDEL", "SNP"] {
            writeln!(writer, "{variant_type},ALL,0,0,0,0,0,0,0,0,,,,,,,,")?;
        }
        return writer
            .flush()
            .with_context(|| format!("failed to flush {}", path.display()));
    }
    for variant_type in ["INDEL", "SNP"].into_iter().filter(|variant_type| {
        all_counts.contains_key(*variant_type) || pass_counts.contains_key(*variant_type)
    }) {
        let all = all_counts.get(variant_type).cloned().unwrap_or_default();
        let pass = pass_counts.get(variant_type).cloned().unwrap_or_default();
        write_summary_row(&mut writer, variant_type, "ALL", &all)?;
        write_summary_row(&mut writer, variant_type, "PASS", &pass)?;
    }
    writer
        .flush()
        .with_context(|| format!("failed to flush {}", path.display()))
}

fn write_summary_row<W: Write>(
    writer: &mut W,
    variant_type: &str,
    filter: &str,
    stats: &QuantifyTypeCounts,
) -> Result<()> {
    let truth_titv = if variant_type == "SNP" {
        ti_tv_ratio(stats.truth_total.ti, stats.truth_total.tv)
    } else {
        String::new()
    };
    let query_titv = if variant_type == "SNP" {
        ti_tv_ratio(stats.query_total.ti, stats.query_total.tv)
    } else {
        String::new()
    };
    writeln!(
        writer,
        "{variant_type},{filter},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        stats.truth_total.total,
        stats.truth_tp.total,
        stats.truth_fn.total,
        stats.query_total.total,
        stats.query_fp.total,
        stats.query_unk.total,
        stats.fp_gt,
        stats.fp_al,
        metric_ratio(stats.truth_tp.total, stats.truth_total.total),
        precision_metric(stats),
        metric_ratio(stats.query_unk.total, stats.query_total.total),
        f1_score(
            stats.truth_tp.total,
            stats.truth_total.total,
            stats.query_tp.total,
            stats.query_tp.total + stats.query_fp.total
        ),
        truth_titv,
        query_titv,
        het_hom_ratio(stats.truth_total.het, stats.truth_total.homalt),
        het_hom_ratio(stats.query_total.het, stats.query_total.homalt),
    )?;
    Ok(())
}

fn write_quantify_extended(
    path: &Path,
    all_counts: &QuantifyCountMaps,
    pass_counts: &QuantifyCountMaps,
    subsets_present: &BTreeMap<String, BTreeSet<String>>,
    options: &ExtendedTableOptions<'_>,
) -> Result<()> {
    let mut writer = BufWriter::new(
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    let mut header = report::EXTENDED_HEADER
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>();
    if options.ci_alpha > 0.0 {
        header.extend(
            [
                "METRIC.Recall.Lower",
                "METRIC.Recall.Upper",
                "METRIC.Precision.Lower",
                "METRIC.Precision.Upper",
                "METRIC.Frac_NA.Lower",
                "METRIC.Frac_NA.Upper",
            ]
            .into_iter()
            .map(str::to_string),
        );
    }
    writeln!(writer, "{}", header.join(","))?;
    if all_counts.by_type.is_empty() && pass_counts.by_type.is_empty() {
        for line in report::empty_comparison_extended_lines(options.subset_size) {
            let mut row = line.split(',').map(str::to_string).collect::<Vec<_>>();
            row.resize(header.len(), String::new());
            writeln!(writer, "{}", row.join(","))?;
        }
        return writer
            .flush()
            .with_context(|| format!("failed to flush {}", path.display()));
    }
    let confidence_size_value = options.confidence_size;
    let confidence_size = confidence_size_value
        .map(|size| format!("{:.6}", size as f64))
        .unwrap_or_default();

    for variant_type in ["INDEL", "SNP"].into_iter().filter(|variant_type| {
        all_counts.by_type.contains_key(*variant_type)
            || pass_counts.by_type.contains_key(*variant_type)
    }) {
        let subtype_labels: Vec<String> = if variant_type == "INDEL" {
            INDEL_SUBTYPES
                .iter()
                .map(|value| value.to_string())
                .collect()
        } else {
            Vec::new()
        };
        let subsets = subsets_present
            .get(variant_type)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect::<Vec<_>>();

        write_extended_pair(
            &mut writer,
            variant_type,
            "*",
            "*",
            all_counts
                .by_type
                .get(variant_type)
                .cloned()
                .unwrap_or_default(),
            pass_counts
                .by_type
                .get(variant_type)
                .cloned()
                .unwrap_or_default(),
            options.subset_size,
            0,
            &confidence_size,
            options.qq_field,
            options.ci_alpha,
        )?;

        for subset in &subsets {
            let (subset_size_for_row, subset_confidence_size) = named_subset_report_sizes(
                subset,
                options.subset_size,
                options.whole_reference_size,
                confidence_size_value,
                options.stratification_sizes,
                options.stratification_confidence_sizes,
            );
            let subset_confidence_size = subset_confidence_size
                .map(|size| format!("{:.6}", size as f64))
                .unwrap_or_default();
            let all = all_counts
                .by_subset_type
                .get(subset)
                .and_then(|types| types.get(variant_type))
                .cloned()
                .unwrap_or_default();
            let pass = pass_counts
                .by_subset_type
                .get(subset)
                .and_then(|types| types.get(variant_type))
                .cloned()
                .unwrap_or_default();
            write_extended_pair(
                &mut writer,
                variant_type,
                "*",
                subset,
                all,
                pass,
                subset_size_for_row,
                options
                    .stratification_levels
                    .get(subset)
                    .copied()
                    .unwrap_or(0),
                &subset_confidence_size,
                options.qq_field,
                options.ci_alpha,
            )?;
        }

        if variant_type == "INDEL" {
            for subtype in &subtype_labels {
                let all = all_counts
                    .by_subtype
                    .get(variant_type)
                    .and_then(|subtypes| subtypes.get(subtype))
                    .cloned()
                    .unwrap_or_default();
                let pass = pass_counts
                    .by_subtype
                    .get(variant_type)
                    .and_then(|subtypes| subtypes.get(subtype))
                    .cloned()
                    .unwrap_or_default();
                write_extended_pair(
                    &mut writer,
                    variant_type,
                    subtype,
                    "*",
                    all,
                    pass,
                    options.subset_size,
                    0,
                    &confidence_size,
                    options.qq_field,
                    options.ci_alpha,
                )?;
                for subset in &subsets {
                    let (subset_size_for_row, subset_confidence_size) = named_subset_report_sizes(
                        subset,
                        options.subset_size,
                        options.whole_reference_size,
                        confidence_size_value,
                        options.stratification_sizes,
                        options.stratification_confidence_sizes,
                    );
                    let subset_confidence_size = subset_confidence_size
                        .map(|size| format!("{:.6}", size as f64))
                        .unwrap_or_default();
                    let all = all_counts
                        .by_subset_subtype
                        .get(subset)
                        .and_then(|types| types.get(variant_type))
                        .and_then(|subtypes| subtypes.get(subtype))
                        .cloned()
                        .unwrap_or_default();
                    let pass = pass_counts
                        .by_subset_subtype
                        .get(subset)
                        .and_then(|types| types.get(variant_type))
                        .and_then(|subtypes| subtypes.get(subtype))
                        .cloned()
                        .unwrap_or_default();
                    write_extended_pair(
                        &mut writer,
                        variant_type,
                        subtype,
                        subset,
                        all,
                        pass,
                        subset_size_for_row,
                        options
                            .stratification_levels
                            .get(subset)
                            .copied()
                            .unwrap_or(0),
                        &subset_confidence_size,
                        options.qq_field,
                        options.ci_alpha,
                    )?;
                }
            }
        }
    }

    writer
        .flush()
        .with_context(|| format!("failed to flush {}", path.display()))
}

fn named_subset_report_sizes(
    subset: &str,
    _active_reference_size: usize,
    whole_reference_size: usize,
    confidence_size: Option<usize>,
    stratification_sizes: &BTreeMap<String, usize>,
    stratification_confidence_sizes: &BTreeMap<String, usize>,
) -> (usize, Option<usize>) {
    let subset_size = match subset {
        "TS_boundary" => whole_reference_size,
        "TS_contained" => confidence_size.unwrap_or(0),
        _ => stratification_sizes.get(subset).copied().unwrap_or(0),
    };
    let subset_confidence_size = match subset {
        "TS_boundary" | "TS_contained" => confidence_size,
        _ => confidence_size.map(|_| {
            stratification_confidence_sizes
                .get(subset)
                .copied()
                .unwrap_or(0)
        }),
    };
    (subset_size, subset_confidence_size)
}

#[allow(clippy::too_many_arguments)]
fn write_extended_pair<W: Write>(
    writer: &mut W,
    variant_type: &str,
    subtype: &str,
    subset: &str,
    all: QuantifyTypeCounts,
    pass: QuantifyTypeCounts,
    subset_size: usize,
    subset_level: usize,
    conf_size: &str,
    qq_field: &str,
    ci_alpha: f64,
) -> Result<()> {
    write_extended_row(
        writer,
        variant_type,
        subtype,
        subset,
        "ALL",
        &all,
        subset_size,
        subset_level,
        conf_size,
        qq_field,
        ci_alpha,
    )?;
    write_extended_row(
        writer,
        variant_type,
        subtype,
        subset,
        "PASS",
        &pass,
        subset_size,
        subset_level,
        conf_size,
        qq_field,
        ci_alpha,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_extended_row<W: Write>(
    writer: &mut W,
    variant_type: &str,
    subtype: &str,
    subset: &str,
    filter: &str,
    stats: &QuantifyTypeCounts,
    subset_size: usize,
    subset_level: usize,
    conf_size: &str,
    qq_field: &str,
    ci_alpha: f64,
) -> Result<()> {
    let supports_titv = variant_type == "SNP";
    let mut row = vec![
        variant_type.to_string(),
        subtype.to_string(),
        subset.to_string(),
        filter.to_string(),
        "*".to_string(),
        qq_field.to_string(),
        "*".to_string(),
        metric_ratio(stats.truth_tp.total, stats.truth_total.total),
        precision_metric(stats),
        metric_ratio(stats.query_unk.total, stats.query_total.total),
        f1_score(
            stats.truth_tp.total,
            stats.truth_total.total,
            stats.query_tp.total,
            stats.query_tp.total + stats.query_fp.total,
        ),
        stats.fp_gt.to_string(),
        stats.fp_al.to_string(),
        if subset == "*" {
            subset_size.to_string()
        } else {
            format!("{:.6}", subset_size as f64)
        },
        conf_size.to_string(),
        format!("{:.6}", subset_level as f64),
    ];
    append_stats(&mut row, &stats.truth_total, supports_titv);
    append_stats(&mut row, &stats.truth_tp, supports_titv);
    append_stats(&mut row, &stats.truth_fn, supports_titv);
    append_stats(&mut row, &stats.query_total, supports_titv);
    append_stats(&mut row, &stats.query_tp, supports_titv);
    append_stats(&mut row, &stats.query_fp, supports_titv);
    append_stats(&mut row, &stats.query_unk, supports_titv);
    if ci_alpha > 0.0 {
        roc::append_ci_cells(
            &mut row,
            [
                (stats.truth_tp.total, stats.truth_total.total),
                (
                    stats.query_tp.total,
                    stats.query_tp.total + stats.query_fp.total,
                ),
                (stats.query_unk.total, stats.query_total.total),
            ],
            ci_alpha,
        );
    }
    writeln!(writer, "{}", row.join(","))?;
    Ok(())
}

fn append_stats(row: &mut Vec<String>, stats: &CountsBucket, supports_titv: bool) {
    row.push(stats.total.to_string());
    if supports_titv {
        row.push(format_count(stats.ti));
        row.push(format_count(stats.tv));
    } else {
        row.push(".".to_string());
        row.push(".".to_string());
    }
    row.push(format_count(stats.het));
    row.push(format_count(stats.homalt));
    if supports_titv {
        row.push(ti_tv_ratio(stats.ti, stats.tv));
    } else {
        row.push(String::new());
    }
    row.push(het_hom_ratio(stats.het, stats.homalt));
}

fn format_count(value: usize) -> String {
    format!("{:.6}", value as f64)
}

fn metric_ratio(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        return "0.0".to_string();
    }
    format_metric(numerator as f64 / denominator as f64)
}

fn precision_metric(stats: &TypeCounts) -> String {
    let denominator = stats.query_tp.total + stats.query_fp.total;
    if denominator == 0 && stats.query_total.total > 0 {
        String::new()
    } else {
        metric_ratio(stats.query_tp.total, denominator)
    }
}

fn f1_score(truth_tp: usize, truth_total: usize, query_tp: usize, query_total: usize) -> String {
    let recall = if truth_total == 0 {
        0.0
    } else {
        truth_tp as f64 / truth_total as f64
    };
    let precision = if query_total == 0 {
        0.0
    } else {
        query_tp as f64 / query_total as f64
    };
    if (recall + precision).abs() < f64::EPSILON {
        return String::new();
    }
    format_metric((2.0 * recall * precision) / (recall + precision))
}

fn ti_tv_ratio(ti: usize, tv: usize) -> String {
    if tv == 0 {
        return String::new();
    }
    format_ratio(ti as f64 / tv as f64)
}

fn het_hom_ratio(het: usize, homalt: usize) -> String {
    if homalt == 0 {
        return String::new();
    }
    format_ratio(het as f64 / homalt as f64)
}

fn format_metric(value: f64) -> String {
    report::format_metric(value)
}

fn format_ratio(value: f64) -> String {
    report::python_repr_float(value)
}

fn write_metrics_json(
    prefix: &Path,
    write_counts: bool,
    roc_indices: &roc::MetricIndices,
) -> Result<()> {
    let mut tables = vec![(
        "summary.metrics",
        "summary.metrics",
        suffixed_report_path(prefix, "summary.csv"),
    )];
    if write_counts {
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
    let table_refs = tables
        .iter()
        .map(|(id, label, path)| (*id, *label, path.as_path()))
        .collect::<Vec<_>>();
    let commandline = std::env::args().collect::<Vec<_>>().join(" ");
    metrics_json::write_metrics_gz_for_module_with_indices(
        &suffixed_report_path(prefix, "metrics.json.gz"),
        "qfy.py.comparison",
        "qfy.py",
        &commandline,
        &table_refs,
        Some(&roc_indices.tables),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::GzDecoder;
    use std::io::Read;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEST_ID: AtomicU64 = AtomicU64::new(0);

    fn fixture(file: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/quantify-simple")
            .join(file)
    }

    fn test_root(label: &str) -> PathBuf {
        let id = TEST_ID.fetch_add(1, Ordering::Relaxed);
        let root =
            std::env::temp_dir().join(format!("hap-quantify-{label}-{}-{id}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        root
    }

    fn write_indexed_vcf_text(path: &Path, text: &str) {
        let source = path.with_extension("source.vcf");
        fs::write(&source, text).unwrap();
        let (headers, records) = vcf::load_raw_vcf(&source).unwrap();
        vcf::write_raw_vcf(path, &headers, &records).unwrap();
        fs::remove_file(source).unwrap();
    }

    fn indexed_fixture(root: &Path) -> PathBuf {
        let input = root.join("input.vcf.gz");
        if !input.exists() {
            // The repository fixture intentionally preserves the historical
            // invalid String-typed, multi-character QQ case used by the parity
            // matrix. Unit tests exercising successful quantify behavior need
            // the valid Float declaration the legacy numeric reader accepts.
            let (mut headers, records) = vcf::load_raw_vcf(&fixture("annotated.vcf")).unwrap();
            for header in &mut headers {
                if header.starts_with("##FORMAT=<ID=QQ,") {
                    *header = header.replace("Type=String", "Type=Float");
                }
            }
            vcf::write_raw_vcf(&input, &headers, &records).unwrap();
        }
        input
    }

    fn args(root: &Path) -> QuantifyArgs {
        QuantifyArgs {
            input_vcf: indexed_fixture(root).display().to_string(),
            report_prefix: root.join("result").display().to_string(),
            reference: fixture("ref.fa").display().to_string(),
            annotation_type: Some("ga4gh".to_string()),
            fp_bedfile: None,
            strat_tsv: None,
            strat_regions: Vec::new(),
            strat_fixchr: false,
            write_vcf: false,
            write_counts: true,
            output_vtc: false,
            preserve_info: false,
            adjust_conf_regions: None,
            threads: None,
            bcf: false,
            logfile: None,
            verbose: false,
            quiet: false,
            force_interactive: false,
            roc: "QUAL".to_string(),
            do_roc: false,
            roc_regions: vec!["*".to_string()],
            roc_filter: None,
            roc_delta: 0.5,
            ci_alpha: 0.0,
            no_json: true,
        }
    }

    fn raw_record(pos: usize, bs: usize, regions: &str, query: &str) -> RawVcfRecord {
        RawVcfRecord::from_line(
            &format!(
                "chr1\t{pos}\t.\tA\tG\t50\tPASS\tBS={bs};Regions={regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:.:.:ti:SNP:het:50\t{query}"
            ),
            Path::new("test.vcf"),
        )
        .unwrap()
    }

    #[test]
    fn superlocus_region_flags_follow_recomputed_confidence_across_the_whole_block() {
        let mut records = vec![
            raw_record(2, 7, "CONF", "0/1:UNK:.:ti:SNP:het:50"),
            raw_record(3, 7, "TS_contained", "0/1:UNK:.:ti:SNP:het:50"),
            raw_record(4, 8, "CONF", "0/1:.:.:ti:SNP:het:50"),
        ];

        propagate_superlocus_annotations(&mut records, "xcmp");

        assert!(has_region(&records[0].info, "TS_boundary"));
        assert!(has_region(&records[1].info, "TS_boundary"));
        assert!(has_region(&records[1].info, "TS_contained"));
        assert!(has_region(&records[2].info, "TS_contained"));
        assert_eq!(
            records[0].sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
    }

    #[test]
    fn boundary_propagation_preserves_existing_contained_membership() {
        let mut records = vec![
            raw_record(2, 7, "CONF,TS_contained", "0/1:TP:gm:ti:SNP:het:50"),
            raw_record(3, 7, "TS_boundary", "0/1:UNK:.:ti:SNP:het:50"),
        ];

        propagate_superlocus_annotations(&mut records, "ga4gh");

        assert!(has_region(&records[0].info, "TS_contained"));
        assert!(has_region(&records[0].info, "TS_boundary"));
    }

    #[test]
    fn ga4gh_reannotation_adds_legacy_fields_from_selected_genotype_alleles() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG,AT\t50\tPASS\tBS=3;Regions=CONF\tGT:BD:BK:QQ\t0/2:TP:gm:40\t1/2:TP:am:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();

        reannotate_ga4gh_record(&mut record);

        assert_eq!(record.format.as_deref(), Some("GT:BD:BK:QQ:BI:BVT:BLT"));
        assert_eq!(
            record.samples[0], "0/2:TP:gm:40:i1_5:INDEL:het",
            "the unused SNP allele must not affect truth BVT/BI"
        );
        assert_eq!(record.samples[1], "1/2:TP:am:35:i1_5,ti:INDEL:hetalt");
        assert_eq!(
            record.sample_map(1).get("BK").map(String::as_str),
            Some("am"),
            "GA4GH BK is supplied by RTG and must be preserved"
        );
    }

    #[test]
    fn ga4gh_reannotation_classifies_homref_nocall_and_halfcall() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3\tGT:BD:BK:QQ\t0/0:TP:gm:40\t./.:N:.:.",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        reannotate_ga4gh_record(&mut record);
        assert_eq!(
            record.sample_map(0).get("BVT").map(String::as_str),
            Some("HOMREF")
        );
        assert_eq!(
            record.sample_map(0).get("BLT").map(String::as_str),
            Some("homref")
        );
        assert_eq!(
            record.sample_map(1).get("BVT").map(String::as_str),
            Some("NOCALL")
        );
        assert_eq!(
            record.sample_map(1).get("BLT").map(String::as_str),
            Some("nocall")
        );

        let halfcall = ga4gh_annotation(&record, "./1");
        assert_eq!(halfcall.bvt, "SNP");
        assert_eq!(halfcall.blt, "halfcall");
        assert_eq!(halfcall.bi, "ti");
        assert_eq!(ga4gh_annotation(&record, "1").blt, "halfcall");
    }

    #[test]
    fn ga4gh_reannotation_marks_overwide_genotype_locations_ambiguous() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t.\tPASS\t.\tGT:BD\t0/1/1:N\t.",
            Path::new("rtg.vcf"),
        )
        .unwrap();

        reannotate_ga4gh_record(&mut record);

        assert_eq!(record.format.as_deref(), Some("GT:BD:BI:BVT:BLT:QQ"));
        assert_eq!(record.samples[0], "0/1/1:N:.:UNK:ambi:.");
        assert_eq!(record.samples[1], ".:.:.:UNK:ambi:0");
    }

    #[test]
    fn ga4gh_count_unk_rewrites_both_samples_outside_confidence() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3\tGT:BD:BK:QQ\t0/1:N:.:40\t0/1:FP:.:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        annotate_regions(&mut record, Some(&[]), &RegionMap::new(), true);
        assert_eq!(
            record.sample_map(0).get("BD").map(String::as_str),
            Some("UNK")
        );
        assert_eq!(
            record.sample_map(1).get("BD").map(String::as_str),
            Some("UNK")
        );

        let mut query_only = RawVcfRecord::from_line(
            "chr1\t4\t.\tA\tG\t50\tPASS\tBS=4\tGT:BD:BK:QQ\t./.:.:.:.\t0/1:FP:.:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        annotate_regions(&mut query_only, Some(&[]), &RegionMap::new(), true);
        assert_eq!(
            query_only.sample_map(0).get("BD").map(String::as_str),
            Some("UNK")
        );
        assert_eq!(
            query_only.sample_map(1).get("BD").map(String::as_str),
            Some("UNK")
        );

        let mut filtered_truth_handoff = RawVcfRecord::from_line(
            "chr1\t5\t.\tA\tG\t50\tPASS\tBS=5\tGT:BD:BK:BI:BVT:BLT:QQ\t./.:.:.:.:NOCALL:nocall:.\t0/1:FP:.:ti:SNP:het:35",
            Path::new("xcmp.vcf"),
        )
        .unwrap();
        annotate_regions_for_samples(
            &mut filtered_truth_handoff,
            Some(&[]),
            &RegionMap::new(),
            true,
            BenchmarkSamples::POSITIONAL,
            true,
        );
        assert_eq!(
            filtered_truth_handoff
                .sample_map(0)
                .get("BD")
                .map(String::as_str),
            Some("."),
            "filtered-truth xcmp handoff preserves a missing truth decision"
        );
        assert_eq!(
            filtered_truth_handoff
                .sample_map(1)
                .get("BD")
                .map(String::as_str),
            Some("UNK")
        );

        let mut filtered_truth_only_handoff = RawVcfRecord::from_line(
            "chr1\t6\t.\tA\tG\t50\tPASS\tBS=6\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t./.:.:.:.:NOCALL:nocall:.",
            Path::new("xcmp.vcf"),
        )
        .unwrap();
        annotate_regions_for_samples(
            &mut filtered_truth_only_handoff,
            Some(&[]),
            &RegionMap::new(),
            true,
            BenchmarkSamples::POSITIONAL,
            true,
        );
        assert_eq!(
            filtered_truth_only_handoff
                .sample_map(0)
                .get("BD")
                .map(String::as_str),
            Some("UNK")
        );
        assert_eq!(
            filtered_truth_only_handoff
                .sample_map(1)
                .get("BD")
                .map(String::as_str),
            Some("."),
            "filtered-truth xcmp handoff preserves either NOCALL side"
        );
    }

    #[test]
    fn ga4gh_fp_match_kinds_feed_count_and_roc_classes() {
        let record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t0/1:FP:am:ti:SNP:het:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        let classified = classify_side(&record, 1).unwrap();
        let mut counts = QuantifyCountMaps::default();
        record_query(&mut counts, &classified);

        assert_eq!(counts.by_type["SNP"].fp_gt, 1);
        assert_eq!(counts.by_type["SNP"].fp_al, 0);
        assert_eq!(query_fp_class(&record), Some("gt"));

        let allele_mismatch = RawVcfRecord::from_line(
            "chr1\t4\t.\tA\tG\t50\tPASS\tBS=4\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t0/1:FP:lm:ti:SNP:het:30",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        assert_eq!(query_fp_class(&allele_mismatch), Some("al"));
    }

    #[test]
    fn xcmp_vtc_and_preserve_info_follow_legacy_cleaning_controls() {
        let source = "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3;Regions=CONF;type=TP;kind=match;ctype=simple:match;gtt1=gt_het;gtt2=gt_het;CUSTOM=kept\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:50\t0/1:TP:gm:ti:SNP:het:50";
        let mut cleaned = RawVcfRecord::from_line(source, Path::new("xcmp.vcf")).unwrap();
        decorate_quantified_record(&mut cleaned, "xcmp", false, true, true);
        assert_eq!(
            cleaned.info,
            "BS=3;Regions=CONF;XCMP=TP:match:gt_het:gt_het:simple:match;VTC=nuc__s,al__s,het__rs"
        );

        let mut preserved = RawVcfRecord::from_line(source, Path::new("xcmp.vcf")).unwrap();
        decorate_quantified_record(&mut preserved, "xcmp", true, false, true);
        assert!(preserved.info.contains("CUSTOM=kept"));
        assert!(preserved.info.contains("type=TP"));
    }

    #[test]
    fn adjusted_confidence_padding_uses_truth_records_inside_raw_confidence() {
        let truth = vec![
            RawVcfRecord::from_line(
                "chr1\t3\t.\tA\tAT\t50\tPASS\t.\tGT\t0/1",
                Path::new("truth.vcf"),
            )
            .unwrap(),
        ];
        let confidence = [vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 2,
            end: 3,
        }];

        let padding = truth_confidence_padding(&truth, &confidence);

        assert_eq!(padding.len(), 1);
        assert_eq!(padding[0].chrom, "chr1");
        assert_eq!((padding[0].start, padding[0].end), (2, 4));
    }

    #[test]
    fn confidence_region_precedes_existing_containment_tag() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3;Regions=TS_contained\tGT:BD\t0/1:TP\t0/1:TP",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        let confidence = [vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 2,
            end: 4,
        }];

        annotate_regions(&mut record, Some(&confidence), &RegionMap::new(), true);

        assert!(record.info.contains("Regions=CONF,TS_contained"));
    }

    #[test]
    fn region_updates_stay_before_preserved_extent_metadata() {
        let mut info = "BS=5;RegionsExtent=5-5;Regions=EXTRA".to_string();

        merge_region_tags(&mut info, &["CONF".to_string()]);

        assert_eq!(info, "BS=5;Regions=EXTRA,CONF;RegionsExtent=5-5");
    }

    #[test]
    fn partial_custom_stratification_overlap_marks_the_superlocus_boundary() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t8\t.\tTACG\tT\t50\tPASS\tBS=5;Regions=CONF,TS_contained\tGT:BD\t0/1:TP\t0/1:TP",
            Path::new("xcmp.vcf"),
        )
        .unwrap();
        let mut stratifications = RegionMap::new();
        stratifications.insert(
            "EXTRA".to_string(),
            vec![
                vcf::BedInterval {
                    chrom: "chr1".to_string(),
                    start: 0,
                    end: 7,
                },
                vcf::BedInterval {
                    chrom: "chr1".to_string(),
                    start: 9,
                    end: 100,
                },
            ],
        );
        let confidence = stratifications["EXTRA"].clone();

        annotate_regions(&mut record, Some(&confidence), &stratifications, false);

        assert!(has_region(&record.info, "EXTRA"));
        assert!(has_region(&record.info, "TS_boundary"));
        propagate_superlocus_annotations(std::slice::from_mut(&mut record), "ga4gh");
        assert!(has_region(&record.info, "CONF"));
        assert!(has_region(&record.info, "EXTRA"));
        assert!(has_region(&record.info, "TS_boundary"));
        assert!(has_region(&record.info, "TS_contained"));
    }

    #[test]
    fn insertion_requires_both_reference_anchors_in_confidence() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tATC\t50\tPASS\t.\tGT:BD\t0/1:TP\t0/1:TP",
            Path::new("scmp.vcf"),
        )
        .unwrap();
        let confidence = [vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 2,
            end: 4,
        }];

        annotate_regions(&mut record, Some(&confidence), &RegionMap::new(), true);

        assert!(has_region(&record.info, "CONF"));
        assert_eq!(
            record.sample_map(0).get("BD").map(String::as_str),
            Some("TP")
        );

        for confidence in [
            vcf::BedInterval {
                chrom: "chr1".to_string(),
                start: 2,
                end: 3,
            },
            vcf::BedInterval {
                chrom: "chr1".to_string(),
                start: 3,
                end: 4,
            },
        ] {
            let mut one_anchor = RawVcfRecord::from_line(
                "chr1\t3\t.\tA\tATC\t50\tPASS\t.\tGT:BD\t0/1:TP\t0/1:TP",
                Path::new("scmp.vcf"),
            )
            .unwrap();
            annotate_regions(
                &mut one_anchor,
                Some(std::slice::from_ref(&confidence)),
                &RegionMap::new(),
                true,
            );
            assert!(!has_region(&one_anchor.info, "CONF"));
        }
    }

    #[test]
    fn confidence_prefixed_stratification_is_folded_into_confidence() {
        let root = test_root("conf-vars");
        let base = root.join("base.bed");
        let vars = root.join("vars.bed");
        fs::write(&base, "chr1\t0\t2\n").unwrap();
        fs::write(&vars, "chr1\t2\t4\n").unwrap();
        let mut options = args(&root);
        options.fp_bedfile = Some(base.display().to_string());
        options.strat_regions = vec![format!("CONF_VARS:{}", vars.display())];
        let contigs = ["chr1".to_string()].into_iter().collect();

        let (confidence, regions, _) = load_regions(&options, &contigs).unwrap();

        let confidence = confidence.expect("combined CONF lane");
        assert_eq!(confidence.len(), 2);
        assert_eq!(region_size(&confidence), 4);
        assert!(!regions.contains_key("CONF_VARS"));
    }

    #[test]
    fn pass_truth_false_negatives_include_filtered_query_matches() {
        let mut counts = QuantifyCountMaps::default();
        counts.by_type.insert(
            "SNP".to_string(),
            QuantifyTypeCounts {
                counts: TypeCounts {
                    truth_total: CountsBucket {
                        total: 10,
                        ti: 7,
                        tv: 3,
                        ..CountsBucket::default()
                    },
                    truth_tp: CountsBucket {
                        total: 6,
                        ti: 4,
                        tv: 2,
                        ..CountsBucket::default()
                    },
                    // Explicit FN classification alone misses truth TPs whose
                    // paired query failed FILTER.
                    truth_fn: CountsBucket {
                        total: 2,
                        ti: 1,
                        tv: 1,
                        ..CountsBucket::default()
                    },
                    ..TypeCounts::default()
                },
                ..QuantifyTypeCounts::default()
            },
        );

        derive_pass_truth_false_negatives(&mut counts);

        let stats = &counts.by_type["SNP"].truth_fn;
        assert_eq!(stats.total, 4);
        assert_eq!(stats.ti, 3);
        assert_eq!(stats.tv, 1);
    }

    #[test]
    fn ga4gh_superlocus_propagates_truth_quality_and_query_nocall_filters() {
        let mut records = vec![
            RawVcfRecord::from_line(
                "chr1\t2\t.\tA\tG\t50\tPASS\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:.\t0/1:TP:gm:ti:SNP:het:30",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
            RawVcfRecord::from_line(
                "chr1\t3\t.\tC\tT\t40\tPASS\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:.\t./.:N:.:.:NOCALL:nocall:.",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
            RawVcfRecord::from_line(
                "chr1\t4\t.\tG\tA\t20\tLowQual\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t0/1:FP:.:ti:SNP:het:20",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
        ];

        propagate_superlocus_annotations(&mut records, "ga4gh");

        assert_eq!(
            records[0].sample_map(0).get("QQ").map(String::as_str),
            Some("30")
        );
        assert_eq!(
            records[1].sample_map(0).get("QQ").map(String::as_str),
            Some("30"),
            "truth TP without a paired query TP receives the minimum TP score in BS"
        );
        assert_eq!(
            records[2].sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(records[1].filter, "LowQual");
    }

    #[test]
    fn ga4gh_superlocus_does_not_propagate_nan_query_quality_to_truth() {
        let mut records = vec![
            RawVcfRecord::from_line(
                "chr1\t2\t.\tA\tG\t50\tPASS\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:.\t0/1:TP:gm:ti:SNP:het:nan",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
        ];

        propagate_superlocus_annotations(&mut records, "ga4gh");

        assert_eq!(
            records[0].sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(
            records[0].sample_map(1).get("QQ").map(String::as_str),
            Some("nan")
        );
    }

    #[test]
    fn ga4gh_output_restores_missing_qualities_and_integer_rendering() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t30\t.\tC\t<DEL>\t.\tPASS\t.\tGT:BD\t0/1:N\t./.:N",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        reannotate_ga4gh_record(&mut record);
        propagate_ga4gh_superlocus_for_samples(
            std::slice::from_mut(&mut record),
            BenchmarkSamples::POSITIONAL,
            false,
            false,
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("0")
        );
        set_format_value(&mut record, 1, "QQ", "60.0");
        normalize_integer_like_format_values(&mut record, "QQ");

        assert_eq!(record.format.as_deref(), Some("GT:BD:BI:BVT:BLT:QQ"));
        assert_eq!(
            record.sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("60")
        );

        let mut somatic_record = RawVcfRecord::from_line(
            "chr1\t30\t.\tC\t<DEL>\t.\tPASS\t.\tGT:BD\t0/1:N\t./.:N",
            Path::new("scmp-somatic.vcf"),
        )
        .unwrap();
        reannotate_ga4gh_record(&mut somatic_record);
        propagate_ga4gh_superlocus_for_samples(
            std::slice::from_mut(&mut somatic_record),
            BenchmarkSamples::POSITIONAL,
            true,
            false,
        );
        assert_eq!(
            somatic_record.sample_map(1).get("QQ").map(String::as_str),
            Some("."),
            "scmp-somatic preserves a missing query score"
        );
    }

    #[test]
    fn ga4gh_output_headers_add_only_missing_legacy_declarations() {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string(),
        ];
        ensure_ga4gh_headers(&mut headers);
        ensure_ga4gh_headers(&mut headers);
        canonicalize_ga4gh_header_order(&mut headers);
        for id in ["GT", "BD", "BK", "BI", "QQ", "BVT", "BLT"] {
            assert_eq!(
                headers
                    .iter()
                    .filter(|header| header.contains(&format!("FORMAT=<ID={id},")))
                    .count(),
                1,
                "header {id} must be present exactly once"
            );
        }
        assert!(headers.iter().any(|header| header.contains("INFO=<ID=BS,")));
        assert!(headers[1].starts_with("##FILTER=<ID=PASS,"));
        assert!(headers.last().unwrap().starts_with("#CHROM"));
    }

    #[test]
    fn ga4gh_accepts_a_string_qq_declaration() {
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FORMAT=<ID=QQ,Number=1,Type=String,Description=\"Score\">".to_string(),
        ];
        let mut record = RawVcfRecord::from_line(
            "chr1\t2\t.\tA\tG\t60\tPASS\t.\tGT:QQ\t0/1:7\t0/1:.",
            Path::new("fixture.vcf"),
        )
        .unwrap();
        assert!(validate_ga4gh_qq_fields(&headers, std::slice::from_ref(&record)).is_ok());

        record.samples[0] = "0/1:60".to_string();
        let error = validate_ga4gh_qq_fields(&headers, &[record]).unwrap_err();
        assert_eq!(error.to_string(), "too many QQ fields at chr1:1");
    }

    #[test]
    fn region_header_uses_the_legacy_bcf_compatible_declaration() {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        ensure_info_header(
            &mut headers,
            "Regions",
            "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">",
        );
        assert_eq!(
            headers[1],
            "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">"
        );
    }

    fn read_gzip(path: &Path) -> String {
        let mut text = String::new();
        GzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    }

    fn write_named_sample_vcf(path: &Path, sample_names: &[&str], samples: &[Vec<&str>]) {
        let mut text = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=10>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision\">\n",
            "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Kind\">\n",
            "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Info\">\n",
            "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"Type\">\n",
            "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"Location\">\n",
            "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Quality\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        )
        .to_string();
        for name in sample_names {
            text.push('\t');
            text.push_str(name);
        }
        text.push('\n');
        for (record_index, record_samples) in samples.iter().enumerate() {
            let pos = record_index + 2;
            let qual = 30 + record_index * 10;
            text.push_str(&format!(
                "chr1\t{pos}\t.\tA\tG\t{qual}\tPASS\tBS={pos}\tGT:BD:BK:BI:BVT:BLT:QQ"
            ));
            for sample in record_samples {
                text.push('\t');
                text.push_str(sample);
            }
            text.push('\n');
        }
        write_indexed_vcf_text(path, &text);
    }

    fn summary_snp_all(path: &Path) -> Vec<String> {
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .find(|line| line.starts_with("SNP,ALL,"))
            .unwrap()
            .split(',')
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn compressed_quantifier_inputs_require_their_index() {
        let root = test_root("input-index");
        let (mut headers, records) = vcf::load_raw_vcf(&fixture("annotated.vcf")).unwrap();
        let chrom_header = headers
            .iter()
            .position(|header| header.starts_with("#CHROM"))
            .unwrap();
        headers.insert(
            chrom_header,
            "##INFO=<ID=BS,Number=1,Type=Integer,Description=\"Benchmark superlocus\">".to_string(),
        );
        for (name, index_suffix) in [("input.vcf.gz", ".tbi"), ("input.bcf", ".csi")] {
            let input = root.join(name);
            vcf::write_raw_vcf(&input, &headers, &records).unwrap();
            let index = PathBuf::from(format!("{}{}", input.display(), index_suffix));
            assert!(index.is_file());
            fs::remove_file(index).unwrap();

            let case_root = root.join(format!("case-{name}"));
            let mut options = args(&case_root);
            options.input_vcf = input.display().to_string();
            let error = run(options).unwrap_err();
            assert!(
                error.to_string().contains("index"),
                "unexpected missing-index error: {error:#}"
            );
        }

        let plain_root = root.join("plain");
        let plain_input = plain_root.join("input.vcf");
        fs::create_dir_all(&plain_root).unwrap();
        fs::copy(fixture("annotated.vcf"), &plain_input).unwrap();
        let mut options = args(&plain_root);
        options.input_vcf = plain_input.display().to_string();
        let error = run(options).unwrap_err();
        assert!(error.to_string().contains("compressed and indexed"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn named_truth_and_query_samples_drive_counts_and_rocs_in_any_column() {
        let root = test_root("named-samples");
        let cases = [
            (
                "reversed",
                vec!["QUERY", "TRUTH"],
                vec![
                    vec!["0/1:FP:lm:ti:SNP:het:30", "0/1:FN:lm:ti:SNP:het:30"],
                    vec!["0/1:TP:gm:ti:SNP:het:40", "0/1:TP:gm:ti:SNP:het:40"],
                ],
            ),
            (
                "extra",
                vec!["NOISE", "QUERY", "AUX", "TRUTH"],
                vec![
                    vec![
                        "0/1:TP:gm:ti:SNP:het:99",
                        "0/1:FP:lm:ti:SNP:het:30",
                        "0/1:TP:gm:ti:SNP:het:98",
                        "0/1:FN:lm:ti:SNP:het:30",
                    ],
                    vec![
                        "0/1:TP:gm:ti:SNP:het:99",
                        "0/1:TP:gm:ti:SNP:het:40",
                        "0/1:TP:gm:ti:SNP:het:98",
                        "0/1:TP:gm:ti:SNP:het:40",
                    ],
                ],
            ),
        ];

        for (label, names, samples) in cases {
            let case_root = root.join(label);
            fs::create_dir_all(&case_root).unwrap();
            let input = case_root.join("input.vcf.gz");
            write_named_sample_vcf(&input, &names, &samples);
            let mut options = args(&case_root);
            options.input_vcf = input.display().to_string();
            options.do_roc = true;
            run(options).unwrap();

            let fields = summary_snp_all(&case_root.join("result.summary.csv"));
            assert_eq!(&fields[2..10], ["2", "1", "1", "2", "1", "0", "0", "1"]);
            let roc = read_gzip(&case_root.join("result.roc.all.csv.gz"));
            assert!(roc.lines().any(|line| {
                let fields = line.split(',').collect::<Vec<_>>();
                fields.first() == Some(&"SNP")
                    && fields.get(6) == Some(&"*")
                    && fields.get(51) == Some(&"1")
            }));
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn missing_truth_or_query_names_disable_rocs_without_positional_counts() {
        let root = test_root("missing-named-samples");
        let input = root.join("named.vcf.gz");
        write_named_sample_vcf(
            &input,
            &["FIRST", "SECOND"],
            &[vec!["0/1:TP:gm:ti:SNP:het:30", "0/1:FP:lm:ti:SNP:het:30"]],
        );
        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.do_roc = true;
        run(options).unwrap();

        let summary = fs::read_to_string(root.join("result.summary.csv")).unwrap();
        assert_eq!(summary.lines().count(), 3);
        assert!(
            summary
                .lines()
                .skip(1)
                .all(|line| line.contains(",ALL,0,0,0,0,0,0,0,0,"))
        );
        assert!(!root.join("result.roc.Locations.SNP.csv.gz").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn explicit_zero_roc_delta_uses_the_low_level_default() {
        let root = test_root("zero-roc-delta");
        let input = root.join("named.vcf.gz");
        write_named_sample_vcf(
            &input,
            &["TRUTH", "QUERY"],
            &[
                vec!["0/1:TP:gm:ti:SNP:het:10", "0/1:TP:gm:ti:SNP:het:10"],
                vec!["0/1:TP:gm:ti:SNP:het:10.05", "0/1:TP:gm:ti:SNP:het:10.05"],
            ],
        );
        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.do_roc = true;
        options.roc_delta = 0.0;
        run(options).unwrap();

        let roc = read_gzip(&root.join("result.roc.Locations.SNP.csv.gz"));
        let levels = roc
            .lines()
            .skip(1)
            .filter_map(|line| line.split(',').nth(6))
            .filter(|level| *level != "*")
            .collect::<Vec<_>>();
        assert_eq!(levels, ["10.000000"]);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn preserve_info_output_declares_and_populates_regions_extent() {
        let root = test_root("regions-extent");
        let mut options = args(&root);
        options.annotation_type = Some("xcmp".to_string());
        options.write_vcf = true;
        options.preserve_info = true;
        run(options).unwrap();

        let (headers, records) =
            vcf::load_raw_vcf(&root.join("result.vcf.gz")).expect("quantified VCF");
        assert_eq!(
            headers
                .iter()
                .filter(|line| line.contains("INFO=<ID=RegionsExtent,"))
                .count(),
            1
        );
        assert_eq!(
            records
                .iter()
                .map(|record| info_value(&record.info, "RegionsExtent").unwrap())
                .collect::<Vec<_>>(),
            ["2-2", "8-9"]
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn confidence_and_named_stratifications_drive_counts_and_sizes() {
        let root = test_root("regions");
        let conf = root.join("confidence.bed");
        let early = root.join("early.bed");
        let late = root.join("late.bed");
        let strata = root.join("regions.tsv");
        fs::write(&conf, "1\t0\t5\n").unwrap();
        fs::write(&early, "1\t0\t5\n").unwrap();
        fs::write(&late, "chr1\t5\t10\n").unwrap();
        fs::write(&strata, "EARLY\tearly.bed\n").unwrap();

        let mut options = args(&root);
        options.fp_bedfile = Some(conf.display().to_string());
        options.strat_tsv = Some(strata.display().to_string());
        options.strat_regions = vec![format!("LATE:{}", late.display())];
        options.strat_fixchr = true;
        run(options).unwrap();

        let summary = fs::read_to_string(root.join("result.summary.csv")).unwrap();
        assert!(
            summary.contains("INDEL,ALL,0,0,0,1,0,1,0,0,0.0,,1.0,"),
            "the call outside confidence must be QUERY.UNK: {summary}"
        );
        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        let early_row = extended
            .lines()
            .find(|line| line.starts_with("SNP,*,EARLY,ALL,"))
            .expect("named TSV stratification row");
        let early_fields = early_row.split(',').collect::<Vec<_>>();
        assert_eq!(early_fields[13], "5.000000");
        assert_eq!(early_fields[14], "5.000000");
        assert!(
            extended
                .lines()
                .any(|line| line.starts_with("INDEL,*,LATE,ALL,")),
            "direct NAME:BED stratification must be applied"
        );
        assert!(root.join("result.roc.all.csv.gz").exists());
        assert!(!root.join("result.roc.Locations.SNP.csv.gz").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn four_column_stratifications_add_dynamic_child_lanes_and_parent_membership() {
        let root = test_root("four-column-regions");
        let regions = root.join("hierarchy.bed");
        fs::write(
            &regions,
            "chr1\t1\t2\tcoding_1\nchr1\t4\t5\tunused\nchr1\t7\t9\tcoding_2\n",
        )
        .unwrap();

        let dynamic_root = root.join("dynamic");
        let mut dynamic = args(&dynamic_root);
        dynamic.strat_regions = vec![format!("EXTRA:{}", regions.display())];
        dynamic.write_vcf = true;
        run(dynamic).unwrap();

        let (_, records) = vcf::load_raw_vcf(&dynamic_root.join("result.vcf.gz")).unwrap();
        assert!(has_region(&records[0].info, "EXTRA"));
        assert!(has_region(&records[0].info, "EXTRA_coding_1"));
        assert!(has_region(&records[1].info, "EXTRA"));
        assert!(has_region(&records[1].info, "EXTRA_coding_2"));
        assert!(
            records
                .iter()
                .all(|record| !has_region(&record.info, "EXTRA_unused"))
        );
        assert!(
            !records
                .iter()
                .any(|record| has_region(&record.info, "EXTRA_coding"))
        );

        let extended = fs::read_to_string(dynamic_root.join("result.extended.csv")).unwrap();
        for (variant_type, subset, size, level) in [
            ("SNP", "EXTRA", "4.000000", "0.000000"),
            ("SNP", "EXTRA_coding_1", "1.000000", "1.000000"),
            ("SNP", "EXTRA_coding_2", "2.000000", "1.000000"),
            ("SNP", "EXTRA_unused", "1.000000", "1.000000"),
            ("INDEL", "EXTRA", "4.000000", "0.000000"),
            ("INDEL", "EXTRA_coding_1", "1.000000", "1.000000"),
            ("INDEL", "EXTRA_coding_2", "2.000000", "1.000000"),
            ("INDEL", "EXTRA_unused", "1.000000", "1.000000"),
        ] {
            let row = extended
                .lines()
                .find(|line| line.starts_with(&format!("{variant_type},*,{subset},ALL,")))
                .unwrap_or_else(|| panic!("missing dynamic subset row {variant_type}/{subset}"));
            let fields = row.split(',').collect::<Vec<_>>();
            assert_eq!(fields[13], size, "wrong size for {variant_type}/{subset}");
            assert_eq!(fields[15], level, "wrong level for {variant_type}/{subset}");
        }
        assert!(
            !extended
                .lines()
                .any(|line| { line.split(',').nth(2) == Some("EXTRA_coding") })
        );
        let roc = read_gzip(&dynamic_root.join("result.roc.all.csv.gz"));
        let dynamic_roc = roc
            .lines()
            .find(|line| {
                let fields = line.split(',').collect::<Vec<_>>();
                fields.first() == Some(&"SNP")
                    && fields.get(2) == Some(&"EXTRA_coding_1")
                    && fields.get(6) == Some(&"*")
            })
            .expect("dynamic child ROC row");
        assert_eq!(dynamic_roc.split(',').nth(15), Some("1.000000"));
        for variant_type in ["SNP", "INDEL"] {
            for filter in ["ALL", "PASS"] {
                let empty_roc = roc
                    .lines()
                    .find(|line| {
                        let fields = line.split(',').collect::<Vec<_>>();
                        fields.first() == Some(&variant_type)
                            && fields.get(1) == Some(&"*")
                            && fields.get(2) == Some(&"EXTRA_unused")
                            && fields.get(3) == Some(&filter)
                            && fields.get(6) == Some(&"*")
                    })
                    .unwrap_or_else(|| {
                        panic!("missing empty dynamic ROC row {variant_type}/{filter}")
                    });
                let fields = empty_roc.split(',').collect::<Vec<_>>();
                assert_eq!(fields[13], "1.000000");
                assert_eq!(fields[15], "1.000000");
                assert_eq!(fields[16], "0");
            }
        }

        let fixed_root = root.join("fixed");
        let mut fixed = args(&fixed_root);
        fixed.strat_regions = vec![format!("=EXTRA:{}", regions.display())];
        fixed.write_vcf = true;
        run(fixed).unwrap();
        let (_, fixed_records) = vcf::load_raw_vcf(&fixed_root.join("result.vcf.gz")).unwrap();
        assert!(
            fixed_records
                .iter()
                .all(|record| has_region(&record.info, "EXTRA"))
        );
        assert!(fixed_records.iter().all(|record| {
            !parse_subsets(&record.info)
                .iter()
                .any(|subset| subset.starts_with("EXTRA_"))
        }));

        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn default_qual_roc_preserves_count_tables_and_writes_roc_files() {
        let root = test_root("roc");
        let baseline_root = root.join("baseline");
        let roc_root = root.join("roc");
        let baseline = args(&baseline_root);
        run(baseline).unwrap();
        let mut options = args(&roc_root);
        options.do_roc = true;
        run(options).unwrap();

        assert_eq!(
            fs::read(baseline_root.join("result.summary.csv")).unwrap(),
            fs::read(roc_root.join("result.summary.csv")).unwrap()
        );
        assert_eq!(
            fs::read(baseline_root.join("result.extended.csv")).unwrap(),
            fs::read(roc_root.join("result.extended.csv")).unwrap()
        );
        for suffix in [
            "roc.all.csv.gz",
            "roc.Locations.SNP.csv.gz",
            "roc.Locations.SNP.PASS.csv.gz",
            "roc.Locations.INDEL.csv.gz",
            "roc.Locations.INDEL.PASS.csv.gz",
        ] {
            let path = roc_root.join(format!("result.{suffix}"));
            assert!(
                path.metadata().unwrap().len() > 0,
                "{} is empty",
                path.display()
            );
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn inherited_roc_controls_flow_through_quantify_outputs() {
        let root = test_root("roc-controls");
        let input = root.join("annotated.vcf.gz");
        let source = fs::read_to_string(fixture("annotated.vcf")).unwrap();
        let source = source
            .replace(
                "##FORMAT=<ID=QQ,Number=1,Type=String",
                "##FORMAT=<ID=QQ,Number=1,Type=Float",
            )
            .replace(
                "##FORMAT=<ID=GT",
                "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"ROC score\">\n##FORMAT=<ID=GT",
            )
            .replace("PASS\tBS=2", "LowQual\tBS=2;SCORE=10.0")
            .replace("PASS\tBS=8", "PASS\tBS=8;SCORE=10.4");
        write_indexed_vcf_text(&input, &source);

        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.roc = "SCORE".to_string();
        options.roc_filter = Some("LowQual".to_string());
        options.roc_delta = 0.0;
        options.ci_alpha = 0.05;
        options.do_roc = true;
        options.no_json = false;
        run(options).unwrap();

        let all = read_gzip(&root.join("result.roc.all.csv.gz"));
        let mut lines = all.lines();
        let header = lines.next().unwrap();
        assert!(header.ends_with("METRIC.Frac_NA.Lower,METRIC.Frac_NA.Upper"));
        assert!(lines.any(|line| {
            let fields = line.split(',').collect::<Vec<_>>();
            fields[5] == "SCORE" && fields[6] == "10.000000"
        }));
        assert!(root.join("result.roc.Locations.SNP.SEL.csv.gz").exists());
        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        let mut extended_lines = extended.lines();
        assert!(
            extended_lines
                .next()
                .unwrap()
                .ends_with("METRIC.Frac_NA.Lower,METRIC.Frac_NA.Upper")
        );
        assert!(extended_lines.all(|line| line.split(',').nth(5) == Some("SCORE")));
        let metrics = read_gzip(&root.join("result.metrics.json.gz"));
        assert!(metrics.contains("\"id\":\"roc.Locations.SNP.SEL\""));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn metrics_json_uses_legacy_table_schema_and_column_types() {
        let root = test_root("metrics-json");
        let mut options = args(&root);
        options.do_roc = true;
        options.no_json = false;
        run(options).unwrap();

        let json = read_gzip(&root.join("result.metrics.json.gz"));
        assert!(json.contains("\"name\":\"qfy.py.comparison\""));
        assert!(json.contains("\"module\":\"qfy.py\""));
        for id in [
            "summary.metrics",
            "all.metrics",
            "roc.all",
            "roc.Locations.INDEL.PASS",
            "roc.Locations.SNP.PASS",
            "roc.Locations.SNP",
            "roc.Locations.INDEL",
        ] {
            assert!(
                json.contains(&format!("\"id\":\"{id}\"")),
                "missing legacy metrics table {id}"
            );
        }
        assert!(
            json.contains("\"type\":\"int64\",\"id\":\"TRUTH.TOTAL\",\"label\":\"TRUTH.TOTAL\"")
        );
        assert!(
            json.contains(
                "\"type\":\"double\",\"id\":\"METRIC.Recall\",\"label\":\"METRIC.Recall\""
            )
        );
        assert!(json.contains("\"type\":\"string\",\"id\":\"Type\",\"label\":\"Type\""));
        assert!(!json.contains("\"tool\":\"hap quantify\""));
        assert!(!json.contains("\"truth_total\""));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn metrics_json_omits_all_table_when_extended_counts_are_disabled() {
        let root = test_root("metrics-json-no-counts");
        let mut options = args(&root);
        options.write_counts = false;
        options.no_json = false;
        run(options).unwrap();

        let json = read_gzip(&root.join("result.metrics.json.gz"));
        assert!(json.contains("\"id\":\"summary.metrics\""));
        assert!(json.contains("\"id\":\"roc.all\""));
        assert!(!json.contains("\"id\":\"all.metrics\""));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn no_write_counts_keeps_summary_and_suppresses_only_extended_table() {
        let root = test_root("no-counts");
        let mut options = args(&root);
        options.write_counts = false;
        run(options).unwrap();

        assert!(root.join("result.summary.csv").is_file());
        assert!(!root.join("result.extended.csv").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn dotted_report_prefix_is_preserved_for_every_output_family() {
        let root = test_root("dotted-prefix");
        let mut options = args(&root);
        options.report_prefix = root.join("sample.v1").display().to_string();
        options.write_vcf = true;
        run(options).unwrap();

        for suffix in [
            "summary.csv",
            "extended.csv",
            "vcf.gz",
            "vcf.gz.tbi",
            "roc.all.csv.gz",
        ] {
            assert!(
                root.join(format!("sample.v1.{suffix}")).is_file(),
                "missing dotted-prefix output {suffix}"
            );
        }
        assert!(!root.join("sample.summary.csv").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bcf_logfile_and_verbose_controls_publish_legacy_standalone_artifacts() {
        let root = test_root("bcf-logfile");
        let mut options = args(&root);
        options.write_vcf = true;
        options.bcf = true;
        options.logfile = Some(root.join("qfy.log").display().to_string());
        options.verbose = true;
        options.do_roc = false;
        run(options).unwrap();

        assert!(root.join("result.bcf").is_file());
        assert!(root.join("result.bcf.csi").is_file());
        assert!(!root.join("result.vcf.gz").exists());
        assert!(root.join("qfy.log").is_file());
        let raw_roc = fs::read_to_string(root.join("result.roc.tsv")).unwrap();
        let mut lines = raw_roc.lines();
        let header = lines.next().unwrap();
        assert!(header.starts_with("FP.al\tFP.gt\tFilter\tGenotype\tMETRIC.F1_Score"));
        let columns = header.split('\t').collect::<Vec<_>>();
        let qq = columns.iter().position(|column| *column == "QQ").unwrap();
        assert!(lines.all(|line| line.split('\t').nth(qq) == Some("*")));
        assert!(header.contains("QUERY.TOTAL.hetalt"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn write_vcf_refuses_to_overwrite_its_input_before_writing_reports() {
        let root = test_root("input-overwrite");
        let input = root.join("result.vcf.gz");
        fs::copy(fixture("annotated.vcf.gz"), &input).unwrap();
        let before = fs::read(&input).unwrap();
        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.write_vcf = true;

        let error = run(options).unwrap_err();
        assert!(error.to_string().contains("cannot overwrite input VCF"));
        assert_eq!(fs::read(&input).unwrap(), before);
        assert!(!root.join("result.summary.csv").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn xcmp_mode_rederives_decisions_instead_of_consuming_final_bd_fields() {
        let root = test_root("xcmp-reannotation");
        let input = root.join("annotated.vcf.gz");
        let reference = root.join("ref.fa");
        let confidence = root.join("confidence.bed");
        fs::write(&reference, ">chr1\nNNACGTNN\n").unwrap();
        fs::write(&confidence, "chr1\t0\t3\n").unwrap();
        write_indexed_vcf_text(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision\">\n",
                "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Kind\">\n",
                "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Info\">\n",
                "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"Type\">\n",
                "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"Location\">\n",
                "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Quality\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n",
                // Final BD fields inside CONF are ignored because old XCMP
                // record-level INFO/type is absent.
                "chr1\t2\t.\tA\tG\t60\tPASS\t.\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:60\t0/1:TP:gm:ti:SNP:het:60\n",
                // Outside CONF, XCMP's count_unk rule replaces the absent
                // decision with UNK for both samples; only QUERY contributes.
                "chr1\t5\t.\tG\tA\t50\tPASS\t.\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:50\t0/1:TP:gm:ti:SNP:het:50\n",
            ),
        );

        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.reference = reference.display().to_string();
        options.fp_bedfile = Some(confidence.display().to_string());
        options.annotation_type = Some("xcmp".to_string());
        run(options).unwrap();

        let summary = fs::read_to_string(root.join("result.summary.csv")).unwrap();
        assert!(summary.contains("SNP,ALL,0,0,0,1,0,1,0,0,0.0,,1.0,"));
        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        let snp_all = extended
            .lines()
            .find(|line| line.starts_with("SNP,*,*,ALL,"))
            .unwrap()
            .split(',')
            .collect::<Vec<_>>();
        assert_eq!(
            snp_all[13], "4",
            "terminal Ns are excluded from Subset.Size"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn xcmp_mode_maps_record_level_fp_to_truth_fn_and_query_fp() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t40\tPASS\ttype=FP;kind=gtmismatch;Regions=CONF\tGT:BD:BK:QQ\t0/1:TP:gm:1\t0/1:TP:gm:2",
            Path::new("input.vcf"),
        )
        .unwrap();
        reannotate_xcmp_record(&mut record, true, "QUAL");
        assert_eq!(
            record.sample_map(0).get("BD").map(String::as_str),
            Some("FN")
        );
        assert_eq!(
            record.sample_map(1).get("BD").map(String::as_str),
            Some("FP")
        );
        assert_eq!(
            record.sample_map(0).get("BK").map(String::as_str),
            Some("am")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("40")
        );
    }

    #[test]
    fn xcmp_custom_format_roc_field_is_copied_into_qq() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t40\tPASS\ttype=TP;kind=match\tGT:BD:BK:QQ:GQX\t0/1:TP:gm:1:17.5\t0/1:TP:gm:2:23.5",
            Path::new("input.vcf"),
        )
        .unwrap();
        reannotate_xcmp_record(&mut record, false, "GQX");
        assert_eq!(
            record.sample_map(0).get("QQ").map(String::as_str),
            Some("17.5")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("23.5")
        );
    }

    #[test]
    fn invalid_quantify_controls_fail_before_reading_inputs() {
        let root = test_root("invalid-controls");
        let baseline = args(&root);

        let mut unknown_type = baseline.clone();
        unknown_type.annotation_type = Some("unknown".to_string());
        assert!(
            run(unknown_type)
                .unwrap_err()
                .to_string()
                .contains("unsupported")
        );

        let mut roc_field = baseline.clone();
        roc_field.roc.clear();
        assert!(run(roc_field).unwrap_err().to_string().contains("empty"));

        let mut roc_regions = baseline.clone();
        roc_regions.roc_regions.push(String::new());
        assert!(run(roc_regions).unwrap_err().to_string().contains("empty"));

        let mut roc_delta = baseline.clone();
        roc_delta.roc_delta = -1.0;
        assert!(
            run(roc_delta)
                .unwrap_err()
                .to_string()
                .contains("nonnegative")
        );

        let mut ci = baseline;
        ci.ci_alpha = 1.5;
        assert!(run(ci).unwrap_err().to_string().contains("ci-alpha"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn reported_bed_size_uses_half_open_bed_span() {
        let intervals = [vcf::BedInterval {
            chrom: "chr1".to_string(),
            start: 0,
            end: 16,
        }];
        assert_eq!(region_size(&intervals), 16);
    }

    #[test]
    fn symbolic_alt_has_no_nucleotide_reference_range() {
        let record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\t<DEL>\t40\tPASS\tEND=7\tGT\t0/1",
            Path::new("symbolic.vcf"),
        )
        .unwrap();
        assert_eq!(effective_reference_range(&record), None);
    }

    #[test]
    fn reference_range_stops_at_first_symbolic_alt() {
        let parse = |alt: &str| {
            RawVcfRecord::from_line(
                &format!("chr1\t3\t.\tA\t{alt}\t40\tPASS\t.\tGT\t0/1"),
                Path::new("mixed.vcf"),
            )
            .unwrap()
        };

        // A leading symbolic allele prevents the later insertion from being
        // considered at all.
        assert_eq!(effective_reference_range(&parse("<DEL>,AT")), None);
        // A preceding nucleotide allele is retained, while nucleotide alleles
        // after the symbolic one are ignored.
        assert_eq!(
            effective_reference_range(&parse("AT,<DEL>,AG")),
            effective_reference_range(&parse("AT"))
        );
        assert_eq!(effective_reference_range(&parse("AT")), Some((3, 4, true)));
    }

    #[test]
    fn missing_alt_is_processed_as_an_empty_allele() {
        let parse = |alt: &str| {
            RawVcfRecord::from_line(
                &format!("chr1\t3\t.\tAT\t{alt}\t40\tPASS\t.\tGT\t0/1"),
                Path::new("missing.vcf"),
            )
            .unwrap()
        };

        let reference_span = Some((3, 4, false));
        assert_eq!(effective_reference_range(&parse(".")), reference_span);
        assert_eq!(effective_reference_range(&parse("")), reference_span);
        // Missing is retained before the symbolic ALT terminates processing;
        // the later nucleotide ALT cannot narrow that reference span.
        assert_eq!(
            effective_reference_range(&parse(".,<DEL>,A")),
            reference_span
        );
    }

    #[test]
    fn named_subset_report_sizes_distinguish_boundary_and_confidence_intersection() {
        let stratification_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);
        let stratification_confidence_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);

        assert_eq!(
            named_subset_report_sizes(
                "TS_boundary",
                100,
                140,
                Some(141),
                &stratification_sizes,
                &stratification_confidence_sizes,
            ),
            (140, Some(141))
        );
        assert_eq!(
            named_subset_report_sizes(
                "EXTRA",
                100,
                140,
                Some(141),
                &stratification_sizes,
                &stratification_confidence_sizes,
            ),
            (138, Some(138))
        );
    }

    #[test]
    fn confidence_intersection_uses_unique_interval_union() {
        let interval = |chrom: &str, start, end| vcf::BedInterval {
            chrom: chrom.to_string(),
            start,
            end,
        };
        let subset = vec![interval("chr1", 0, 98), interval("chrX", 0, 40)];
        let confidence = vec![
            interval("chr1", 0, 98),
            interval("chrX", 0, 40),
            interval("chr1", 1, 4),
        ];

        assert_eq!(region_size(&confidence), 141);
        assert_eq!(region_intersection_size(&subset, &confidence), 138);
    }
}
