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
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const INDEL_SUBTYPES: [&str; 9] = [
    "C16_PLUS", "C1_5", "C6_15", "D16_PLUS", "D1_5", "D6_15", "I16_PLUS", "I1_5", "I6_15",
];

type RegionMap = BTreeMap<String, Vec<vcf::BedInterval>>;
type LoadedRegions = (Option<Vec<vcf::BedInterval>>, RegionMap);

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
}

#[derive(Clone, Debug, Default)]
struct QuantifyCountMaps {
    by_type: BTreeMap<String, TypeCounts>,
    by_subset_type: BTreeMap<String, BTreeMap<String, TypeCounts>>,
    by_subtype: BTreeMap<String, BTreeMap<String, TypeCounts>>,
    by_subset_subtype: BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>>,
}

pub fn run(args: QuantifyArgs) -> Result<()> {
    validate_options(&args)?;
    let annotation_type = args.annotation_type.as_deref().unwrap_or("xcmp");
    let prefix = Path::new(&args.report_prefix);
    let output_vcf = suffixed_report_path(prefix, "vcf.gz");
    if args.write_vcf && paths_refer_to_same_file(Path::new(&args.input_vcf), &output_vcf) {
        bail!(
            "cannot overwrite input VCF: {} would be overwritten by output {}",
            args.input_vcf,
            output_vcf.display()
        );
    }
    let (mut headers, mut records) = vcf::load_raw_vcf(Path::new(&args.input_vcf))?;
    let reference = crate::fasta::read_sequences(Path::new(&args.reference))?;
    let reference_contigs = reference.keys().cloned().collect::<BTreeSet<_>>();
    let (confidence, stratifications) = load_regions(&args, &reference_contigs)?;
    for record in &mut records {
        annotate_regions(
            record,
            confidence.as_deref(),
            &stratifications,
            annotation_type == "ga4gh",
        );
        if annotation_type == "xcmp" {
            reannotate_xcmp_record(record, confidence.is_some(), &args.roc);
        } else {
            reannotate_ga4gh_record(record);
        }
    }
    propagate_superlocus_annotations(&mut records, annotation_type);
    let subset_size = contigs_in_input(&records)
        .into_iter()
        .filter_map(|contig| {
            reference
                .get(&contig)
                .map(|sequence| n_trimmed_length(sequence))
        })
        .sum::<usize>();

    let mut all_counts = QuantifyCountMaps::default();
    let mut pass_counts = QuantifyCountMaps::default();
    let mut subsets_present = BTreeMap::<String, BTreeSet<String>>::new();
    let stratification_sizes = stratifications
        .iter()
        .map(|(name, intervals)| (name.clone(), region_size(intervals)))
        .collect::<BTreeMap<_, _>>();
    let confidence_size = confidence.as_deref().map(region_size);

    for record in &records {
        let truth = classify_side(record, 0);
        let query = classify_side(record, 1);

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
            subset_size,
            &stratification_sizes,
            confidence_size,
        )?;
    }
    if args.write_vcf {
        if annotation_type == "ga4gh" {
            ensure_ga4gh_headers(&mut headers);
        }
        let lines = records
            .iter()
            .map(RawVcfRecord::to_line)
            .collect::<Vec<_>>();
        vcf::write_indexed_vcf(&output_vcf, &headers, lines.iter().map(String::as_str))?;
    }
    let mut rows = records
        .iter()
        .enumerate()
        .map(|(index, record)| AnnotatedRow {
            sort_key: (record.chrom.clone(), record.pos, index, 0),
            line: record.to_line(),
            query_pass: record.is_pass(),
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        })
        .collect::<Vec<_>>();
    rows.sort_by(|left, right| left.sort_key.cmp(&right.sort_key));
    let roc_options = roc::RocOptions {
        qq_field: args.roc.clone(),
        ignored_filters: args
            .roc_filter
            .as_deref()
            .into_iter()
            .flat_map(|filters| filters.split(|ch: char| ch.is_whitespace() || ";,".contains(ch)))
            .filter(|filter| !filter.is_empty())
            .map(str::to_string)
            .collect(),
        roc_regions: args.roc_regions.iter().cloned().collect(),
        delta: args.roc_delta,
        ci_alpha: args.ci_alpha,
    };
    let roc_indices = roc::write_roc_files_with_options(
        prefix,
        &rows,
        subset_size,
        confidence_size.unwrap_or(0),
        &roc_options,
    )?;
    if !args.do_roc {
        compact_no_roc_outputs(prefix)?;
    }
    if !args.no_json {
        write_metrics_json(prefix, args.write_counts, args.do_roc, &roc_indices)?;
    }
    Ok(())
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

    let mut regions = BTreeMap::new();
    for (name, path) in paths {
        let intervals = load_region_bed(&path, reference_contigs, args.strat_fixchr)
            .with_context(|| format!("failed to load stratification region {name}"))?;
        // QuantifyRegions::load collapses every lane whose label begins
        // with CONF into the reserved confidence lane. hap.py relies on
        // that behavior for its generated CONF_VARS truth padding.
        if name.starts_with("CONF") {
            confidence.get_or_insert_default().extend(intervals);
        } else {
            regions.insert(name, intervals);
        }
    }
    Ok((confidence, regions))
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
    if name.is_empty() {
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
        intervals.push(vcf::BedInterval { chrom, start, end });
    }
    Ok(intervals)
}

fn annotate_regions(
    record: &mut RawVcfRecord,
    confidence: Option<&[vcf::BedInterval]>,
    stratifications: &RegionMap,
    rewrite_ga4gh_decisions: bool,
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
    let mut additions = stratifications
        .iter()
        .filter(|(_, intervals)| overlaps(intervals))
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();

    if let Some(confidence) = confidence {
        remove_region_tag(&mut record.info, "CONF");
        remove_region_tag(&mut record.info, "TS_boundary");
        if overlaps(confidence) {
            additions.insert(0, "CONF".to_string());
        } else if rewrite_ga4gh_decisions {
            // GA4GHQuantify's `count_unk` rule rewrites both samples outside
            // CONF. Truth-side UNK is omitted from truth counts; query-side
            // UNK supplies QUERY.UNK and ROC unknown counts.
            replace_existing_decision(record, 0, "UNK");
            replace_existing_decision(record, 1, "UNK");
        }
    }
    if !additions.is_empty() {
        merge_region_tags(&mut record.info, &additions);
    }
}

/// `BlockQuantify::count` finalizes region flags per benchmarking superlocus.
/// A block crossing the recomputed CONF boundary marks every member as
/// `TS_boundary`; a wholly confident block marks every member `TS_contained`.
/// Records without a non-negative BS value form singleton blocks.
fn propagate_superlocus_annotations(records: &mut [RawVcfRecord], annotation_type: &str) {
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
        let has_non_confident = records[start..end]
            .iter()
            .any(|record| !has_region(&record.info, "CONF"));
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
            for record in &mut records[start..end] {
                if record.sample_map(0).get("BD").map(String::as_str) != Some("TP") {
                    set_format_value(record, 0, "QQ", ".");
                }
            }
        } else {
            propagate_ga4gh_superlocus(&mut records[start..end]);
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
    ensure_format_fields(record, &["BI", "BVT", "BLT"]);
    for sample_index in 0..record.samples.len() {
        let gt = record
            .sample_map(sample_index)
            .get("GT")
            .cloned()
            .unwrap_or_else(|| "./.".to_string());
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
        "hemi"
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

fn propagate_ga4gh_superlocus(records: &mut [RawVcfRecord]) {
    let minimum_tp_qq = records
        .iter()
        .filter_map(|record| {
            let query = record.sample_map(1);
            (query.get("BD").map(String::as_str) == Some("TP"))
                .then(|| query.get("QQ").cloned())
                .flatten()
        })
        .filter(|score| score != "." && score.parse::<f64>().is_ok())
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
        let truth_tp = record.sample_map(0).get("BD").map(String::as_str) == Some("TP");
        let query = record.sample_map(1);
        let direct_query_qq = (query.get("BD").map(String::as_str) == Some("TP"))
            .then(|| query.get("QQ").cloned())
            .flatten()
            .filter(|score| score != "." && score.parse::<f64>().is_ok());
        let truth_qq = if truth_tp {
            direct_query_qq.as_ref().or(minimum_tp_qq.as_ref())
        } else {
            None
        };
        set_format_value(record, 0, "QQ", truth_qq.map(String::as_str).unwrap_or("."));

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
        if alt.is_empty() || alt == "." || alt.starts_with('<') {
            continue;
        }
        let alt_bytes = alt.as_bytes();
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

/// Reproduce the old XCMP quantifier rather than consuming GA4GH `BD` fields.
/// XCMP decisions live in record-level `INFO/type`, `kind`, and `ctype`; a
/// finalized hap.py VCF intentionally lacks them, which is why the legacy qfy
/// lane ignores calls inside CONF and labels calls outside CONF as UNK.
fn reannotate_xcmp_record(
    record: &mut RawVcfRecord,
    has_confidence_regions: bool,
    roc_field: &str,
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
        let suppressed = import_fail || no_call || (sample_index == 1 && query_filtered);
        let sample_decision = if import_fail {
            "N"
        } else if suppressed || decision.is_empty() {
            "."
        } else if sample_index == 0 && decision == "FP" {
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
    fields.push(format!("Regions={}", tags.join(",")));
    *info = fields.join(";");
}

fn region_size(intervals: &[vcf::BedInterval]) -> usize {
    intervals
        .iter()
        .map(|interval| interval.end.saturating_sub(interval.start))
        .sum()
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
    })
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
    add_variant_stats(
        &mut counts
            .by_type
            .entry(classified.variant_type.clone())
            .or_default()
            .query_total,
        classified,
    );
    add_variant_stats(
        query_bucket(
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
            &mut counts
                .by_subtype
                .entry(classified.variant_type.clone())
                .or_default()
                .entry(subtype.clone())
                .or_default()
                .query_total,
            classified,
        );
        add_variant_stats(
            query_bucket(
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
            &mut counts
                .by_subset_type
                .entry(subset.clone())
                .or_default()
                .entry(classified.variant_type.clone())
                .or_default()
                .query_total,
            classified,
        );
        add_variant_stats(
            query_bucket(
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
                &mut counts
                    .by_subset_subtype
                    .entry(subset.clone())
                    .or_default()
                    .entry(classified.variant_type.clone())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default()
                    .query_total,
                classified,
            );
            add_variant_stats(
                query_bucket(
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
    all_counts: &BTreeMap<String, TypeCounts>,
    pass_counts: &BTreeMap<String, TypeCounts>,
) -> Result<()> {
    let mut writer = BufWriter::new(fs::File::create(path)?);
    writeln!(
        writer,
        "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio"
    )?;
    for variant_type in ["INDEL", "SNP"].into_iter().filter(|variant_type| {
        all_counts.contains_key(*variant_type) || pass_counts.contains_key(*variant_type)
    }) {
        let all = all_counts.get(variant_type).cloned().unwrap_or_default();
        let pass = pass_counts.get(variant_type).cloned().unwrap_or_default();
        write_summary_row(&mut writer, variant_type, "ALL", &all)?;
        write_summary_row(&mut writer, variant_type, "PASS", &pass)?;
    }
    Ok(())
}

fn write_summary_row<W: Write>(
    writer: &mut W,
    variant_type: &str,
    filter: &str,
    stats: &TypeCounts,
) -> Result<()> {
    writeln!(
        writer,
        "{variant_type},{filter},{},{},{},{},{},{},0,0,{},{},{},{},{},{},{},{}",
        stats.truth_total.total,
        stats.truth_tp.total,
        stats.truth_fn.total,
        stats.query_total.total,
        stats.query_fp.total,
        stats.query_unk.total,
        metric_ratio(stats.truth_tp.total, stats.truth_total.total),
        precision_metric(stats),
        metric_ratio(stats.query_unk.total, stats.query_total.total),
        f1_score(
            stats.truth_tp.total,
            stats.truth_total.total,
            stats.query_tp.total,
            stats.query_tp.total + stats.query_fp.total
        ),
        ti_tv_ratio(stats.truth_total.ti, stats.truth_total.tv),
        ti_tv_ratio(stats.query_total.ti, stats.query_total.tv),
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
    subset_size: usize,
    stratification_sizes: &BTreeMap<String, usize>,
    confidence_size: Option<usize>,
) -> Result<()> {
    let mut writer = BufWriter::new(fs::File::create(path)?);
    writeln!(writer, "{}", report::EXTENDED_HEADER.join(","))?;
    let confidence_size_value = confidence_size;
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
            subset_size,
            &confidence_size,
        )?;

        for subset in &subsets {
            let subset_size_for_row = named_subset_size(
                subset,
                subset_size,
                confidence_size_value,
                stratification_sizes,
            );
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
                &confidence_size,
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
                    subset_size,
                    &confidence_size,
                )?;
                for subset in &subsets {
                    let subset_size_for_row = named_subset_size(
                        subset,
                        subset_size,
                        confidence_size_value,
                        stratification_sizes,
                    );
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
                        &confidence_size,
                    )?;
                }
            }
        }
    }

    Ok(())
}

fn named_subset_size(
    subset: &str,
    whole_reference_size: usize,
    confidence_size: Option<usize>,
    stratification_sizes: &BTreeMap<String, usize>,
) -> usize {
    match subset {
        "TS_boundary" => whole_reference_size,
        "TS_contained" => confidence_size.unwrap_or(0),
        _ => stratification_sizes.get(subset).copied().unwrap_or(0),
    }
}

#[allow(clippy::too_many_arguments)]
fn write_extended_pair<W: Write>(
    writer: &mut W,
    variant_type: &str,
    subtype: &str,
    subset: &str,
    all: TypeCounts,
    pass: TypeCounts,
    subset_size: usize,
    conf_size: &str,
) -> Result<()> {
    write_extended_row(
        writer,
        variant_type,
        subtype,
        subset,
        "ALL",
        &all,
        subset_size,
        conf_size,
    )?;
    write_extended_row(
        writer,
        variant_type,
        subtype,
        subset,
        "PASS",
        &pass,
        subset_size,
        conf_size,
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
    stats: &TypeCounts,
    subset_size: usize,
    conf_size: &str,
) -> Result<()> {
    let supports_titv = variant_type == "SNP";
    let mut row = vec![
        variant_type.to_string(),
        subtype.to_string(),
        subset.to_string(),
        filter.to_string(),
        "*".to_string(),
        "QUAL".to_string(),
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
        "0".to_string(),
        "0".to_string(),
        if subset == "*" {
            subset_size.to_string()
        } else {
            format!("{:.6}", subset_size as f64)
        },
        conf_size.to_string(),
        "0.000000".to_string(),
    ];
    append_stats(&mut row, &stats.truth_total, supports_titv);
    append_stats(&mut row, &stats.truth_tp, supports_titv);
    append_stats(&mut row, &stats.truth_fn, supports_titv);
    append_stats(&mut row, &stats.query_total, supports_titv);
    append_stats(&mut row, &stats.query_tp, supports_titv);
    append_stats(&mut row, &stats.query_fp, supports_titv);
    append_stats(&mut row, &stats.query_unk, supports_titv);
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
    if (value.fract()).abs() < f64::EPSILON {
        format!("{value:.1}")
    } else {
        value.to_string()
    }
}

fn write_metrics_json(
    prefix: &Path,
    write_counts: bool,
    do_roc: bool,
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
    tables.push((
        "roc.all",
        "roc.all",
        suffixed_report_path(prefix, "roc.all.csv.gz"),
    ));
    if do_roc {
        for id in [
            "roc.Locations.INDEL",
            "roc.Locations.SNP.PASS",
            "roc.Locations.SNP",
            "roc.Locations.INDEL.PASS",
        ] {
            let path = suffixed_report_path(prefix, &format!("{id}.csv.gz"));
            if path.is_file() {
                tables.push((id, id, path));
            }
        }
        for id in ["roc.Locations.INDEL.SEL", "roc.Locations.SNP.SEL"] {
            let path = suffixed_report_path(prefix, &format!("{id}.csv.gz"));
            if path.exists() {
                tables.push((id, id, path));
            }
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

    fn args(root: &Path) -> QuantifyArgs {
        QuantifyArgs {
            input_vcf: fixture("annotated.vcf").display().to_string(),
            report_prefix: root.join("result").display().to_string(),
            reference: fixture("ref.fa").display().to_string(),
            annotation_type: Some("ga4gh".to_string()),
            fp_bedfile: None,
            strat_tsv: None,
            strat_regions: Vec::new(),
            strat_fixchr: false,
            write_vcf: false,
            write_counts: true,
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
        assert!(has_region(&records[2].info, "TS_contained"));
        assert_eq!(
            records[0].sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
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

        let (confidence, regions) = load_regions(&options, &contigs).unwrap();

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
            TypeCounts {
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
    fn ga4gh_output_headers_add_only_missing_legacy_declarations() {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string(),
        ];
        ensure_ga4gh_headers(&mut headers);
        ensure_ga4gh_headers(&mut headers);
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
        assert!(headers.last().unwrap().starts_with("#CHROM"));
    }

    fn read_gzip(path: &Path) -> String {
        let mut text = String::new();
        GzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
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
        let input = root.join("annotated.vcf");
        let source = fs::read_to_string(fixture("annotated.vcf")).unwrap();
        let source = source
            .replace(
                "##FORMAT=<ID=GT",
                "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"ROC score\">\n##FORMAT=<ID=GT",
            )
            .replace("PASS\tBS=2", "LowQual\tBS=2;SCORE=10.0")
            .replace("PASS\tBS=8", "PASS\tBS=8;SCORE=10.4");
        fs::write(&input, source).unwrap();

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
        let input = root.join("annotated.vcf");
        let reference = root.join("ref.fa");
        let confidence = root.join("confidence.bed");
        fs::write(&reference, ">chr1\nNNACGTNN\n").unwrap();
        fs::write(&confidence, "chr1\t0\t3\n").unwrap();
        fs::write(
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
        )
        .unwrap();

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
}
