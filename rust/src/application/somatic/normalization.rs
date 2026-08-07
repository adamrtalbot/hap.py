//! Cohesive responsibility extracted from the command façade.

use super::allele_frequency::parse_af_bins;
use super::features::csv_join;
use super::features::write_simple_table;
use super::metrics::python2_counter_indices;
use super::{AmbiguousInterval, FilteredRawRecord, QueryClass};
use crate::adapters::report;
use crate::adapters::vcf;
use crate::application::SomaticArgs;
use crate::domain::{Interval, QueryProvenance, RawVcfRecord};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

pub(super) fn validate_args(args: &SomaticArgs) -> Result<()> {
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

pub(super) fn resolve_toggle(enabled: bool, explicitly_disabled: bool) -> bool {
    enabled && !explicitly_disabled
}

pub(super) fn selected_normalizations(args: &SomaticArgs) -> (bool, bool) {
    (
        args.normalize_truth || args.normalize_all,
        args.normalize_query || args.normalize_all,
    )
}

pub(super) fn validate_af_bins(raw: &str) -> Result<()> {
    parse_af_bins(raw).map(|_| ())
}

pub(super) fn classification_bed_chrom(
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

pub(super) fn load_classification_bed(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    fixchr_truth: bool,
) -> Result<Vec<Interval>> {
    let mut intervals = vcf::load_bed(path, &BTreeSet::new())?;
    for interval in &mut intervals {
        interval.chrom = classification_bed_chrom(&interval.chrom, reference_contigs, fixchr_truth);
    }
    Ok(intervals)
}

pub(super) fn load_fp_explanation_bed(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    fixchr_truth: bool,
) -> Result<Vec<AmbiguousInterval>> {
    let text = vcf::read_text(path)
        .with_context(|| format!("failed to read false-positive BED {}", path.display()))?;
    let mut intervals = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 3 {
            bail!(
                "false-positive BED line {} has fewer than 3 columns in {}",
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
            interval: Interval {
                chrom: classification_bed_chrom(fields[0], reference_contigs, fixchr_truth),
                start,
                end,
            },
            label: "FP".to_string(),
            details: fields[3..]
                .iter()
                .map(|value| (*value).to_string())
                .collect(),
        });
    }
    Ok(intervals)
}

pub(super) fn load_ambiguous_beds(
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
                interval: Interval {
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

pub(super) fn record_ambiguous_explanation(
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

pub(super) fn write_legacy_count_table(
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

pub(super) fn classify_query(
    chrom: &str,
    start: usize,
    end: usize,
    fp_regions: &[Interval],
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

pub(super) fn normalize_somatic_records(
    records: Vec<RawVcfRecord>,
    reference_sequences: &BTreeMap<String, String>,
) -> Vec<RawVcfRecord> {
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

pub(super) const MAX_SOMATIC_CONTIG_RECORDS: usize = 10_000_000;

pub(super) fn prepare_somatic_cache(
    source: &Path,
    cache: &Path,
    headers: &[String],
    normalize: bool,
    reference_sequences: Option<&BTreeMap<String, String>>,
) -> Result<()> {
    if !normalize {
        let records = vcf::open_validated_vcf(source)?
            .map(|record| record.map(|record| record.raw().clone()));
        return vcf::write_raw_vcf_iter(cache, headers, records);
    }
    let reference_sequences = reference_sequences.context("normalization requires a reference")?;
    let mut spool = tempfile::NamedTempFile::new().context("failed to create somatic spool")?;
    {
        let mut writer = BufWriter::new(spool.as_file_mut());
        for header in headers {
            writeln!(writer, "{header}")?;
        }
        let mut current_chrom: Option<String> = None;
        let mut contig_records = Vec::new();
        for record in vcf::open_validated_vcf(source)? {
            let record = record?.raw().clone();
            if current_chrom
                .as_deref()
                .is_some_and(|chrom| chrom != record.chrom)
            {
                write_normalized_somatic_contig(
                    &mut writer,
                    &mut contig_records,
                    reference_sequences,
                )?;
            }
            current_chrom = Some(record.chrom.clone());
            contig_records.push(record);
            if contig_records.len() > MAX_SOMATIC_CONTIG_RECORDS {
                bail!(
                    "somatic contig in {} exceeds the {} record active-window limit",
                    source.display(),
                    MAX_SOMATIC_CONTIG_RECORDS
                );
            }
        }
        write_normalized_somatic_contig(&mut writer, &mut contig_records, reference_sequences)?;
        writer.flush()?;
    }
    let records = vcf::open_validated_vcf(spool.path())?
        .map(|record| record.map(|record| record.raw().clone()));
    vcf::write_raw_vcf_iter(cache, headers, records)
}

pub(super) fn write_normalized_somatic_contig(
    writer: &mut dyn Write,
    records: &mut Vec<RawVcfRecord>,
    reference_sequences: &BTreeMap<String, String>,
) -> Result<()> {
    for record in normalize_somatic_records(std::mem::take(records), reference_sequences) {
        let checked = vcf::ValidatedVcfRecord::try_from_raw(record, QueryProvenance::Unavailable)?;
        writeln!(writer, "{}", checked.raw().to_line())?;
    }
    Ok(())
}

pub(super) struct ContigSpool {
    pub(super) chrom: String,
    pub(super) path: tempfile::TempPath,
    pub(super) count: usize,
}

impl ContigSpool {
    pub(super) fn load(&self) -> Result<Vec<FilteredRawRecord>> {
        vcf::open_validated_vcf(&self.path)?
            .map(|record| {
                let record = record?.raw().clone();
                Ok(FilteredRawRecord {
                    key: vcf::VariantKey {
                        chrom: record.chrom.clone(),
                        pos: record.pos,
                        ref_allele: record.ref_allele.clone(),
                        alt_allele: record.alt_allele.clone(),
                    },
                    record,
                })
            })
            .collect()
    }
}

pub(super) fn spool_filtered_contigs(
    path: &Path,
    source_path: &Path,
    options: &RawFilterOptions<'_>,
) -> Result<Vec<ContigSpool>> {
    let mut spools: Vec<ContigSpool> = Vec::new();
    let mut closed = BTreeSet::new();
    let mut active: Option<(String, tempfile::NamedTempFile, usize)> = None;
    for record in vcf::open_validated_vcf(path)? {
        let Some(record) = filter_raw_record(record?.raw().clone(), source_path, options)? else {
            continue;
        };
        if active
            .as_ref()
            .is_none_or(|(chrom, _, _)| chrom != &record.key.chrom)
        {
            if let Some((chrom, mut file, count)) = active.take() {
                file.as_file_mut().flush()?;
                closed.insert(chrom.clone());
                spools.push(ContigSpool {
                    chrom,
                    path: file.into_temp_path(),
                    count,
                });
            }
            if closed.contains(&record.key.chrom) {
                bail!(
                    "somatic records for chromosome {} are not contiguous in {}",
                    record.key.chrom,
                    path.display()
                );
            }
            active = Some((
                record.key.chrom.clone(),
                tempfile::NamedTempFile::new().context("failed to create somatic contig spool")?,
                0,
            ));
        }
        let (chrom, file, count) = active.as_mut().expect("somatic spool was just created");
        writeln!(file.as_file_mut(), "{}", record.record.to_line())?;
        *count += 1;
        if *count > MAX_SOMATIC_CONTIG_RECORDS {
            bail!(
                "somatic contig {} exceeds the {} record active-window limit",
                chrom,
                MAX_SOMATIC_CONTIG_RECORDS
            );
        }
    }
    if let Some((chrom, mut file, count)) = active {
        file.as_file_mut().flush()?;
        spools.push(ContigSpool {
            chrom,
            path: file.into_temp_path(),
            count,
        });
    }
    Ok(spools)
}

pub(super) fn normalize_somatic_alleles(record: &mut RawVcfRecord, reference: &[u8]) {
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

pub(super) struct RawFilterOptions<'a> {
    pub(super) reference_contigs: &'a BTreeSet<String>,
    pub(super) fixchr: bool,
    pub(super) pass_only: bool,
    pub(super) regions: Option<&'a [Interval]>,
    pub(super) targets: Option<&'a [Interval]>,
    pub(super) locations: Option<&'a [vcf::LocationFilter]>,
}

#[cfg(test)]
pub(super) fn filter_raw_records(
    records: Vec<RawVcfRecord>,
    path: &Path,
    options: &RawFilterOptions<'_>,
) -> Result<Vec<FilteredRawRecord>> {
    records
        .into_iter()
        .map(|record| filter_raw_record(record, path, options))
        .filter_map(|result| result.transpose())
        .collect()
}

pub(super) fn filter_raw_record(
    record: RawVcfRecord,
    path: &Path,
    options: &RawFilterOptions<'_>,
) -> Result<Option<FilteredRawRecord>> {
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
    let normalized =
        vcf::ValidatedVcfRecord::try_from_raw(normalized, QueryProvenance::Unavailable)?
            .raw()
            .clone();
    Ok(Some(FilteredRawRecord {
        key,
        record: normalized,
    }))
}

pub(super) fn somatic_chrom(
    chrom: &str,
    _reference_contigs: &BTreeSet<String>,
    fixchr: bool,
) -> String {
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

pub(super) fn calls_terminal_non_ref(record: &RawVcfRecord) -> bool {
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

pub(super) fn somatic_commandline(args: &SomaticArgs) -> String {
    let process_args = std::env::args().collect::<Vec<_>>();
    if process_args.iter().any(|value| value == "somatic") {
        return process_args.join(" ");
    }

    let mut parts = vec!["hap".to_string(), "somatic".to_string()];

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
