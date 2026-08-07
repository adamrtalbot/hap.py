use crate::adapters::{fasta, vcf};
use crate::application::{ValidateArgs as RawValidateArgs, ValidatedValidateArgs as ValidateArgs};
use crate::domain::RawVcfRecord;
use crate::output::OutputTransaction;
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

const WARNING_REFPADDING: usize = 0;
const WARNING_OVERLAP: usize = 1;
const WARNING_SYMALT: usize = 2;
const WARNING_UNCERTAIN_LENGTH: usize = 3;
const WARNING_BCFERROR: usize = 4;
const WARNING_COUNT: usize = 5;

#[derive(Default)]
struct ValidationCounts {
    records: usize,
    ref_records: usize,
    nonref_records: usize,
    haploid: usize,
    diploid: usize,
    polyploid: usize,
    haploid_x: bool,
    diploid_x: bool,
    warnings: [usize; WARNING_COUNT],
}

#[derive(Default)]
struct PreviousRecord {
    chrom: String,
    end: Option<usize>,
    allele_counts: Vec<usize>,
}

#[derive(Default)]
struct RecordGenotypes {
    any_ref: bool,
    any_nonref: bool,
    any_haploid: bool,
    any_diploid: bool,
    any_polyploid: bool,
    any_symbolic: bool,
    any_uncertain: bool,
    allele_counts: Vec<usize>,
}

struct VcfHeader {
    sample_count: usize,
    contigs: BTreeSet<String>,
    info: BTreeSet<String>,
    format: BTreeMap<String, HeaderValueType>,
    filters: BTreeSet<String>,
}

#[derive(Clone, Copy)]
enum HeaderValueType {
    Integer,
    Float,
    Other,
}

struct FormatParseError {
    key: String,
    value_type: HeaderValueType,
    invalid_character: char,
    extreme_value: bool,
}

pub(crate) fn run(args: ValidateArgs) -> Result<()> {
    let stderr = std::io::stderr();
    run_with_diagnostics(args.into_inner(), &mut stderr.lock())
}

fn run_with_diagnostics<W: Write + ?Sized>(
    mut args: RawValidateArgs,
    diagnostics: &mut W,
) -> Result<()> {
    let outputs = [args.output_json.as_deref(), args.errors_bed.as_deref()]
        .into_iter()
        .flatten()
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    if outputs.is_empty() {
        return run_with_diagnostics_inner(args, diagnostics);
    }
    let inputs = [
        Some(args.input.as_str()),
        args.reference.as_deref(),
        args.regions_bedfile.as_deref(),
        args.targets_bedfile.as_deref(),
    ]
    .into_iter()
    .flatten()
    .map(PathBuf::from)
    .collect::<Vec<_>>();
    let transaction = OutputTransaction::files(&inputs, &outputs)?;
    if let Some(output) = args.output_json.as_mut() {
        *output = transaction
            .staged_file(Path::new(output))?
            .to_string_lossy()
            .into_owned();
    }
    if let Some(output) = args.errors_bed.as_mut() {
        *output = transaction
            .staged_file(Path::new(output))?
            .to_string_lossy()
            .into_owned();
    }
    run_with_diagnostics_inner(args, diagnostics).map_err(|error| {
        anyhow::anyhow!(
            "failed to produce validation outputs {}: {error:#}",
            display_paths(&outputs)
        )
    })?;
    transaction.commit()
}

fn run_with_diagnostics_inner<W: Write + ?Sized>(
    args: RawValidateArgs,
    diagnostics: &mut W,
) -> Result<()> {
    let reference_contigs = if let Some(reference) = &args.reference {
        fasta::contig_lengths(Path::new(reference))?
            .into_keys()
            .collect::<BTreeSet<_>>()
    } else {
        BTreeSet::new()
    };
    let regions = args
        .regions_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &reference_contigs))
        .transpose()?;
    let targets = args
        .targets_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &reference_contigs))
        .transpose()?;
    let locations = args
        .locations
        .as_deref()
        .map(|text| vcf::parse_locations(text, &reference_contigs))
        .transpose()?;

    let mut records = vcf::open_validated_vcf(Path::new(&args.input))?;
    let mut header = VcfHeader::from_lines(records.headers())?;
    let location_is_lowercase_x = args
        .locations
        .as_deref()
        .is_some_and(location_selects_lowercase_x);
    let reference_sequences = args
        .reference
        .as_ref()
        .map(|path| fasta::read_sequences(Path::new(path)))
        .transpose()?;

    let mut counts = ValidationCounts::default();
    let mut previous = PreviousRecord::default();
    let mut errors = args
        .errors_bed
        .as_deref()
        .map(|path| {
            let parent = Path::new(path)
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            tempfile::NamedTempFile::new_in(parent)
                .map(|file| (path.to_string(), BufWriter::new(file)))
                .with_context(|| format!("failed to stage {path}"))
        })
        .transpose()?;
    let mut parsed_any_record = false;
    let mut previous_record_failed_to_parse = false;
    let mut reported_extreme_format_value = false;

    for record in &mut records {
        let mut record = record.map_err(restore_legacy_genotype_error)?;
        if record.samples.len() < header.sample_count {
            writeln!(
                diagnostics,
                "[E::vcf_parse_format] Number of columns at {}:{} does not match the number of samples ({} vs {})",
                record.chrom,
                record.pos,
                record.samples.len(),
                header.sample_count
            )?;
            if !parsed_any_record || previous_record_failed_to_parse {
                break;
            }
            previous_record_failed_to_parse = true;
            continue;
        }
        record.try_update(|raw| {
            raw.samples.truncate(header.sample_count);
            Ok(())
        })?;

        if let Some(parse_error) = header.format_parse_error(&record) {
            if parse_error.extreme_value && !reported_extreme_format_value {
                let integer_suffix = match parse_error.value_type {
                    HeaderValueType::Integer => " and set to missing",
                    HeaderValueType::Float | HeaderValueType::Other => "",
                };
                writeln!(
                    diagnostics,
                    "[W::vcf_parse_format] Extreme FORMAT/{} value encountered{} at {}:{}",
                    parse_error.key, integer_suffix, record.chrom, record.pos
                )?;
                reported_extreme_format_value = true;
            }
            writeln!(
                diagnostics,
                "[E::vcf_parse_format] Invalid character '{}' in '{}' FORMAT field at {}:{}",
                parse_error.invalid_character, parse_error.key, record.chrom, record.pos
            )?;
            if !parsed_any_record || previous_record_failed_to_parse {
                break;
            }
            previous_record_failed_to_parse = true;
            continue;
        }
        parsed_any_record = true;
        previous_record_failed_to_parse = false;

        let translation_error = header.translation_error(&record, diagnostics)?;
        if translation_error != 0 {
            if args.check_bcf_errors {
                bail!(
                    "Record at {}:{} will not translate into BCF. Check if the header is incomplete (error code {}). The header must have all contigs present as #contig entries (contrary to the htslib error message, tabix indexing is not sufficient), and all the INFO and FORMAT types must match the values in all records.",
                    record.chrom,
                    record.pos,
                    translation_error
                );
            }
            if counts.warnings[WARNING_BCFERROR] == 0 {
                writeln!(
                    diagnostics,
                    "[W] Record at {}:{} will not translate into BCF. Check if the header is incomplete  (error code {}) -- all records like this are skipped.",
                    record.chrom, record.pos, translation_error
                )?;
            }
            counts.warnings[WARNING_BCFERROR] += 1;
            continue;
        }

        let chrom = if reference_contigs.is_empty() {
            record.chrom.clone()
        } else {
            vcf::normalize_chrom(&record.chrom, &reference_contigs)
        };
        record.try_update(|raw| {
            raw.chrom = chrom;
            Ok(())
        })?;

        if let Some(filters) = locations.as_deref()
            && !filters
                .iter()
                .any(|filter| filter.matches(&record.chrom, record.pos))
        {
            continue;
        }
        if let Some(bed) = regions.as_deref()
            && !bed
                .iter()
                .any(|interval| interval.matches(&record.chrom, record.pos))
        {
            continue;
        }
        if let Some(bed) = targets.as_deref()
            && !bed
                .iter()
                .any(|interval| interval.matches(&record.chrom, record.pos))
        {
            continue;
        }
        if args.apply_filters && !record.is_pass() {
            continue;
        }
        if args.limit_records.is_some_and(|limit| {
            limit != -1 && (limit < 0 || counts.records as u64 >= limit as u64)
        }) {
            break;
        }

        let genotype = inspect_genotypes(&record)?;
        update_warning_counts(
            &record,
            &genotype,
            &mut previous,
            &mut counts,
            args.strict_homref,
            args.all_warnings,
            diagnostics,
        )?;
        update_record_counts(&record, &genotype, &mut counts, location_is_lowercase_x);

        if let Some(reference_sequences) = &reference_sequences
            && let Some(reason) = validate_record(&record, reference_sequences)
            && let Some((_, errors)) = errors.as_mut()
        {
            writeln!(
                errors,
                "{}\t{}\t{}\t{}",
                record.chrom,
                record.pos.saturating_sub(1),
                record.end_pos(),
                reason
            )?;
        }

        if args
            .message_every
            .is_some_and(|every| every > 0 && counts.records.is_multiple_of(every as usize))
        {
            writeln!(diagnostics, "[PROGRESS] {}:{}", record.chrom, record.pos)?;
        }
        counts.records += 1;
    }

    write_warning_summaries(&counts, diagnostics)?;

    if let Some(path) = &args.output_json {
        let json = format!(
            "{{\"OVERLAP\":{},\"REFPADDING\":{},\"SYMALT\":{},\"UNCERTAINLENGTH\":{},\"diploid\":{},\"haploid\":{},\"male\":{},\"nonref\":{},\"polyploid\":{},\"records\":{},\"ref\":{}}}\n",
            counts.warnings[WARNING_OVERLAP],
            counts.warnings[WARNING_REFPADDING],
            counts.warnings[WARNING_SYMALT],
            counts.warnings[WARNING_UNCERTAIN_LENGTH],
            counts.diploid,
            counts.haploid,
            counts.haploid_x && !counts.diploid_x,
            counts.nonref_records,
            counts.polyploid,
            counts.records,
            counts.ref_records,
        );
        fs::write(path, json).with_context(|| format!("failed to write {path}"))?;
    }

    if let Some((path, mut errors)) = errors {
        errors
            .flush()
            .context("failed to flush validation errors BED")?;
        let staged = errors.into_inner().map_err(|error| error.into_error())?;
        staged
            .persist(&path)
            .map_err(|error| error.error)
            .with_context(|| format!("failed to publish {path}"))?;
    }

    Ok(())
}

fn display_paths(paths: &[PathBuf]) -> String {
    paths
        .iter()
        .map(|path| path.display().to_string())
        .collect::<Vec<_>>()
        .join(", ")
}

impl VcfHeader {
    fn from_lines(headers: &[String]) -> Result<Self> {
        let sample_count = headers
            .iter()
            .find(|line| line.starts_with("#CHROM\t"))
            .map(|line| line.split('\t').count().saturating_sub(9))
            .unwrap_or(0);
        if sample_count == 0 {
            bail!("input VCF has no samples; legacy vcfcheck requires at least one sample");
        }
        Ok(Self {
            sample_count,
            contigs: header_ids(headers, "##contig="),
            info: header_ids(headers, "##INFO="),
            format: header_format_types(headers),
            filters: header_ids(headers, "##FILTER="),
        })
    }

    fn translation_error<W: Write + ?Sized>(
        &mut self,
        record: &RawVcfRecord,
        diagnostics: &mut W,
    ) -> Result<usize> {
        let mut error = 0;
        if self.contigs.insert(record.chrom.clone()) {
            writeln!(
                diagnostics,
                "[W::vcf_parse] Contig '{}' is not defined in the header. (Quick workaround: index the file with tabix.)",
                record.chrom
            )?;
            error |= 1;
        }
        for filter in record
            .filter
            .split(';')
            .filter(|filter| *filter != "." && *filter != "PASS")
        {
            if self.filters.insert(filter.to_string()) {
                writeln!(
                    diagnostics,
                    "[W::vcf_parse_filter] FILTER '{filter}' is not defined in the header"
                )?;
                error |= 2;
            }
        }
        for field in record
            .info
            .split(';')
            .filter(|field| *field != "." && !field.is_empty())
        {
            let key = field.split_once('=').map_or(field, |(key, _)| key);
            if self.info.insert(key.to_string()) {
                writeln!(
                    diagnostics,
                    "[W::vcf_parse_info] INFO '{key}' is not defined in the header, assuming Type=String"
                )?;
                error |= 2;
            }
        }
        for key in record.format_keys() {
            if !self.format.contains_key(key) {
                writeln!(
                    diagnostics,
                    "[W::vcf_parse_format] FORMAT '{key}' at {}:{} is not defined in the header, assuming Type=String",
                    record.chrom, record.pos
                )?;
                self.format.insert(key.to_string(), HeaderValueType::Other);
                error |= 2;
            }
        }
        Ok(error)
    }

    fn format_parse_error(&self, record: &RawVcfRecord) -> Option<FormatParseError> {
        let format_keys = record.format_keys();
        record.samples.iter().find_map(|sample| {
            format_keys
                .iter()
                .zip(sample.split(':'))
                .find_map(|(key, value)| {
                    let value_type = *self.format.get(*key)?;
                    let (invalid_character, extreme_value) = value_type.parse_error(value)?;
                    Some(FormatParseError {
                        key: (*key).to_string(),
                        value_type,
                        invalid_character,
                        extreme_value,
                    })
                })
        })
    }
}

impl HeaderValueType {
    fn parse_error(self, value: &str) -> Option<(char, bool)> {
        if value.is_empty() || matches!(self, Self::Other) {
            return None;
        }
        value
            .split(',')
            .filter(|part| *part != ".")
            .find_map(|part| match self {
                Self::Integer => integer_parse_error(part),
                Self::Float => float_parse_error(part),
                Self::Other => None,
            })
    }
}

fn integer_parse_error(value: &str) -> Option<(char, bool)> {
    let unsigned = value
        .strip_prefix('+')
        .or_else(|| value.strip_prefix('-'))
        .unwrap_or(value);
    let invalid = unsigned
        .chars()
        .find(|character| !character.is_ascii_digit());
    invalid.map(|character| {
        (
            character,
            !unsigned.as_bytes().first().is_some_and(u8::is_ascii_digit),
        )
    })
}

fn float_parse_error(value: &str) -> Option<(char, bool)> {
    let unsigned = value
        .strip_prefix('+')
        .or_else(|| value.strip_prefix('-'))
        .unwrap_or(value);
    if value.parse::<f64>().is_ok()
        || ["nan", "inf", "infinity"]
            .iter()
            .any(|special| unsigned.eq_ignore_ascii_case(special))
    {
        return None;
    }

    let bytes = unsigned.as_bytes();
    let mut index = 0;
    while bytes.get(index).is_some_and(u8::is_ascii_digit) {
        index += 1;
    }
    let mut digit_count = index;
    if bytes.get(index) == Some(&b'.') {
        index += 1;
        let decimal_start = index;
        while bytes.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
        digit_count += index - decimal_start;
    }
    if bytes
        .get(index)
        .is_some_and(|byte| matches!(*byte, b'e' | b'E'))
        && digit_count > 0
    {
        index += 1;
        if bytes
            .get(index)
            .is_some_and(|byte| matches!(*byte, b'+' | b'-'))
        {
            index += 1;
        }
        while bytes.get(index).is_some_and(u8::is_ascii_digit) {
            index += 1;
        }
    }
    let invalid = unsigned[index..]
        .chars()
        .next()
        .or_else(|| unsigned.chars().next())?;
    Some((invalid, digit_count == 0))
}

fn header_ids(headers: &[String], prefix: &str) -> BTreeSet<String> {
    headers
        .iter()
        .filter_map(|line| line.strip_prefix(prefix)?.strip_prefix("<ID="))
        .filter_map(|fields| fields.split([',', '>']).next())
        .map(str::to_string)
        .collect()
}

fn header_format_types(headers: &[String]) -> BTreeMap<String, HeaderValueType> {
    headers
        .iter()
        .filter_map(|line| line.strip_prefix("##FORMAT=")?.strip_prefix('<'))
        .filter_map(|fields| {
            let id = header_attribute(fields, "ID")?;
            let value_type = match header_attribute(fields, "Type")? {
                "Integer" => HeaderValueType::Integer,
                "Float" => HeaderValueType::Float,
                _ => HeaderValueType::Other,
            };
            Some((id.to_string(), value_type))
        })
        .collect()
}

fn header_attribute<'a>(fields: &'a str, key: &str) -> Option<&'a str> {
    fields.split(',').find_map(|field| {
        let (field_key, value) = field.split_once('=')?;
        (field_key == key).then_some(value.trim_end_matches('>'))
    })
}

fn location_selects_lowercase_x(location: &str) -> bool {
    location
        .split(',')
        .next()
        .and_then(|location| location.split(':').next())
        == Some("x")
}

/// Keeps the checked adapter boundary while preserving vcfcheck's fatal
/// diagnostic for a genotype that names an absent alternate allele.
fn restore_legacy_genotype_error(error: anyhow::Error) -> anyhow::Error {
    let message = error.to_string();
    let Some((_, coordinate_and_error)) = message.split_once(" at ") else {
        return error;
    };
    let Some((coordinate, _)) = coordinate_and_error.split_once(": allele index ") else {
        return error;
    };
    if !coordinate_and_error.contains(" is out of bounds for ") {
        return error;
    }
    error.context(format!(
        "Call with invalid genotype (non-existent allele) at {coordinate}"
    ))
}

fn inspect_genotypes(record: &RawVcfRecord) -> Result<RecordGenotypes> {
    let mut result = RecordGenotypes::default();
    let alts = record.alt_allele.split(',').collect::<Vec<_>>();
    let Some(gt_index) = record.format_keys().iter().position(|key| *key == "GT") else {
        result.allele_counts.resize(record.samples.len(), 0);
        return Ok(result);
    };

    for sample in &record.samples {
        let gt = sample.split(':').nth(gt_index).unwrap_or(".");
        let alleles = gt.split(['/', '|']).collect::<Vec<_>>();
        if alleles.len() == 1 {
            result.any_haploid = true;
        } else if alleles.len() == 2 {
            if alleles[0] != alleles[1] {
                result.any_diploid = true;
            }
        } else if alleles.len() > 2 {
            result.any_polyploid = true;
        }

        let mut allele_count = 0usize;
        for allele in alleles {
            if allele == "." || allele.is_empty() {
                continue;
            }
            let index = allele.parse::<usize>().with_context(|| {
                format!(
                    "invalid GT allele '{allele}' at {}:{}",
                    record.chrom, record.pos
                )
            })?;
            if index == 0 {
                result.any_ref = true;
                continue;
            }
            result.any_nonref = true;
            let Some(alt) = alts.get(index - 1) else {
                bail!(
                    "Call with invalid genotype (non-existent allele) at {}:{}",
                    record.chrom,
                    record.pos + usize::from(reference_padding(record) > 0)
                );
            };
            if alt.starts_with('<') {
                result.any_symbolic = true;
            }
            if alt.len() > 1 && (alt.contains('*') || alt.contains('.')) {
                result.any_uncertain = true;
            }
            if alt.starts_with('<') || *alt == "." || *alt == "*" || alt.is_empty() {
                continue;
            }
            allele_count += 1;
        }
        result.allele_counts.push(allele_count);
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
fn update_warning_counts<W: Write + ?Sized>(
    record: &RawVcfRecord,
    genotype: &RecordGenotypes,
    previous: &mut PreviousRecord,
    counts: &mut ValidationCounts,
    strict_homref: bool,
    all_warnings: bool,
    diagnostics: &mut W,
) -> Result<()> {
    if previous.chrom != record.chrom {
        previous.end = None;
        previous.allele_counts.clear();
    }

    let ref_padding = reference_padding(record);
    let adjusted_start = record.pos.saturating_sub(1) + usize::from(ref_padding > 0);
    let overlap = previous.end.is_some_and(|end| adjusted_start < end);
    if overlap {
        for (sample, current_count) in genotype.allele_counts.iter().enumerate() {
            let current_count = *current_count + strict_ref_count(record, sample, strict_homref);
            let previous_count = previous.allele_counts.get(sample).copied().unwrap_or(0);
            if current_count + previous_count > 2 {
                emit_warning(
                    counts,
                    WARNING_OVERLAP,
                    all_warnings,
                    diagnostics,
                    format!(
                        "[W] overlapping records at {}:{} for sample {sample}",
                        record.chrom, adjusted_start
                    ),
                )?;
            }
        }
    }

    if ref_padding > 1 {
        emit_warning(
            counts,
            WARNING_REFPADDING,
            all_warnings,
            diagnostics,
            format!(
                "[W] variant at {}:{} has more than one base of reference padding",
                record.chrom,
                adjusted_start + 1
            ),
        )?;
    }
    if genotype.any_symbolic {
        emit_warning(
            counts,
            WARNING_SYMALT,
            all_warnings,
            diagnostics,
            format!(
                "[W] Symbolic / SV ALT alleles at {}:{}",
                record.chrom, adjusted_start
            ),
        )?;
    }
    if genotype.any_uncertain {
        emit_warning(
            counts,
            WARNING_UNCERTAIN_LENGTH,
            all_warnings,
            diagnostics,
            format!(
                "[W] Alleles with uncertain length at {}:{}",
                record.chrom, adjusted_start
            ),
        )?;
    }

    previous.chrom.clone_from(&record.chrom);
    previous.end = Some(record_end(record));
    previous.allele_counts = genotype
        .allele_counts
        .iter()
        .enumerate()
        .map(|(sample, count)| *count + strict_ref_count(record, sample, strict_homref))
        .collect();
    Ok(())
}

fn strict_ref_count(record: &RawVcfRecord, sample: usize, strict: bool) -> usize {
    if !strict {
        return 0;
    }
    record
        .sample_map(sample)
        .get("GT")
        .map(|gt| gt.split(['/', '|']).filter(|allele| *allele == "0").count())
        .unwrap_or(0)
}

fn reference_padding(record: &RawVcfRecord) -> usize {
    let mut max_match = record.ref_allele.len();
    for alt in record.alt_allele.split(',') {
        if alt == "." || alt.starts_with('<') {
            return 0;
        }
        let prefix = record
            .ref_allele
            .bytes()
            .zip(alt.bytes())
            .take_while(|(reference, alternate)| reference == alternate)
            .count();
        max_match = max_match.min(prefix);
    }
    max_match
}

fn record_end(record: &RawVcfRecord) -> usize {
    record
        .info
        .split(';')
        .find_map(|field| field.strip_prefix("END=")?.parse::<usize>().ok())
        .map(|end| end.saturating_sub(1))
        .unwrap_or_else(|| record.end_pos().saturating_sub(1))
}

fn emit_warning<W: Write + ?Sized>(
    counts: &mut ValidationCounts,
    warning: usize,
    all_warnings: bool,
    diagnostics: &mut W,
    message: String,
) -> Result<()> {
    if all_warnings || counts.warnings[warning] == 0 {
        writeln!(diagnostics, "{message}")?;
    }
    counts.warnings[warning] += 1;
    Ok(())
}

fn update_record_counts(
    record: &RawVcfRecord,
    genotype: &RecordGenotypes,
    counts: &mut ValidationCounts,
    location_is_lowercase_x: bool,
) {
    counts.ref_records += usize::from(genotype.any_ref);
    counts.nonref_records += usize::from(genotype.any_nonref);
    counts.haploid += usize::from(genotype.any_haploid);
    counts.diploid += usize::from(genotype.any_diploid);
    counts.polyploid += usize::from(genotype.any_polyploid);
    if record.chrom == "X"
        || record.chrom == "chrX"
        || record.chrom == "chrx"
        || location_is_lowercase_x
    {
        counts.haploid_x |= genotype.any_haploid;
        counts.diploid_x |= genotype.any_diploid || genotype.any_polyploid;
    }
}

fn write_warning_summaries<W: Write + ?Sized>(
    counts: &ValidationCounts,
    diagnostics: &mut W,
) -> Result<()> {
    if counts.warnings[WARNING_BCFERROR] > 0 {
        writeln!(
            diagnostics,
            "[W] Variants that will cause trouble when writing BCF: {}",
            counts.warnings[WARNING_BCFERROR]
        )?;
    }
    for (warning, label) in [
        (
            WARNING_REFPADDING,
            "Variants that have >1 base of reference padding",
        ),
        (
            WARNING_OVERLAP,
            "Variants that overlap on the reference allele",
        ),
        (WARNING_SYMALT, "Variants that have symbolic ALT alleles"),
        (
            WARNING_UNCERTAIN_LENGTH,
            "Variants that have alleles with uncertain length",
        ),
    ] {
        if counts.warnings[warning] > 0 {
            writeln!(diagnostics, "[W] {label}: {}", counts.warnings[warning])?;
        }
    }
    writeln!(
        diagnostics,
        "[I] Total VCF records:         {}",
        counts.records
    )?;
    writeln!(
        diagnostics,
        "[I] Non-reference VCF records: {}",
        counts.nonref_records
    )?;
    Ok(())
}

fn validate_record(
    record: &RawVcfRecord,
    reference_sequences: &BTreeMap<String, String>,
) -> Option<String> {
    let reference = reference_sequences.get(&record.chrom)?;
    // DNA is ASCII; byte slicing avoids allocating a chr-scale Vec<char> per record.
    let bases = reference.as_bytes();
    if record.pos == 0 || record.end_pos() > bases.len() {
        return Some("OUT_OF_RANGE".to_string());
    }
    let observed_bytes = &bases[record.pos - 1..record.end_pos()];
    let ref_bytes = record.ref_allele.as_bytes();
    let case_insensitive_match = observed_bytes.len() == ref_bytes.len()
        && observed_bytes.iter().zip(ref_bytes.iter()).all(|(a, b)| {
            let au = a.to_ascii_uppercase();
            let bu = b.to_ascii_uppercase();
            au == bu || au == b'N' || bu == b'N'
        });
    if !case_insensitive_match {
        return Some(format!(
            "REF_MISMATCH:{}",
            String::from_utf8_lossy(observed_bytes)
        ));
    }
    let gt = record.sample_map(0).get("GT").cloned().unwrap_or_default();
    let ploidy = gt
        .split(['/', '|'])
        .filter(|token| !token.is_empty())
        .count();
    if ploidy > 2 {
        return Some("POLYPLOID_GT".to_string());
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn write_vcf(path: &Path, records: &[&str]) -> Result<()> {
        let mut text = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FILTER=<ID=LowQual,Description=\"Synthetic filter\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n".to_string();
        for record in records {
            text.push_str(record);
            text.push('\n');
        }
        fs::write(path, text)?;
        Ok(())
    }

    fn args(input: &Path, output: &Path) -> crate::application::ValidateArgs {
        crate::application::ValidateArgs {
            input: input.display().to_string(),
            reference: None,
            output_json: Some(output.display().to_string()),
            errors_bed: None,
            locations: None,
            regions_bedfile: None,
            targets_bedfile: None,
            apply_filters: false,
            limit_records: None,
            message_every: None,
            strict_homref: false,
            check_bcf_errors: false,
            all_warnings: false,
        }
    }

    fn run_with_diagnostics(
        args: crate::application::ValidateArgs,
        diagnostics: &mut dyn Write,
    ) -> Result<()> {
        super::run_with_diagnostics(args, diagnostics)
    }

    #[test]
    fn filter_limit_and_progress_match_processed_record_order() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        write_vcf(
            &input,
            &[
                "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
                "chr1\t2\t.\tA\tG\t.\tLowQual\t.\tGT\t0/1",
                "chr1\t3\t.\tA\tT\t.\tPASS\t.\tGT\t1/1",
                "chr1\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0/1",
            ],
        )?;
        let mut args = args(&input, &output);
        args.apply_filters = true;
        args.limit_records = Some(2);
        args.message_every = Some(1);
        let mut diagnostics = Vec::new();
        run_with_diagnostics(args, &mut diagnostics)?;

        let json = fs::read_to_string(output)?;
        assert!(
            json.ends_with('\n'),
            "legacy JsonCpp FastWriter adds a newline"
        );
        assert!(json.contains("\"records\":2"));
        assert!(json.contains("\"nonref\":2"));
        assert!(!json.contains("PROGRESS"));
        let diagnostics = String::from_utf8(diagnostics)?;
        assert!(diagnostics.contains("[PROGRESS] chr1:1"));
        assert!(diagnostics.contains("[PROGRESS] chr1:3"));
        assert!(!diagnostics.contains("chr1:2"));
        assert!(!diagnostics.contains("chr1:4"));
        Ok(())
    }

    #[test]
    fn validation_rejects_output_path_collisions_before_writing() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("result.txt");
        write_vcf(&input, &["chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1"])?;

        let mut colliding_outputs = args(&input, &output);
        colliding_outputs.errors_bed = Some(output.display().to_string());
        assert!(
            run_with_diagnostics(colliding_outputs, &mut Vec::new())
                .unwrap_err()
                .to_string()
                .contains("distinct paths")
        );

        let overwrite_input = args(&input, &input);
        assert!(
            run_with_diagnostics(overwrite_input, &mut Vec::new())
                .unwrap_err()
                .to_string()
                .contains("would overwrite input")
        );
        assert!(fs::read_to_string(&input)?.starts_with("##fileformat"));
        Ok(())
    }

    #[test]
    fn strict_homref_changes_overlap_accounting() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        write_vcf(
            &input,
            &[
                "chr1\t2\t.\tAAA\tA\t.\tPASS\t.\tGT\t0/1",
                "chr1\t2\t.\tAAA\tA\t.\tPASS\t.\tGT\t0/1",
            ],
        )?;

        for (strict, expected) in [(false, 0), (true, 1)] {
            let output = directory.path().join(format!("strict-{strict}.json"));
            let mut args = args(&input, &output);
            args.strict_homref = strict;
            run_with_diagnostics(args, &mut Vec::new())?;
            assert!(fs::read_to_string(output)?.contains(&format!("\"OVERLAP\":{expected}")));
        }
        Ok(())
    }

    #[test]
    fn all_warnings_controls_repeated_diagnostics_not_counts() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        write_vcf(
            &input,
            &[
                "chr1\t2\t.\tA\t<DEL>\t.\tPASS\t.\tGT\t0/1",
                "chr1\t4\t.\tA\t<DUP>\t.\tPASS\t.\tGT\t0/1",
            ],
        )?;

        for (all, expected_occurrences) in [(false, 1), (true, 2)] {
            let output = directory.path().join(format!("warnings-{all}.json"));
            let mut args = args(&input, &output);
            args.all_warnings = all;
            let mut diagnostics = Vec::new();
            run_with_diagnostics(args, &mut diagnostics)?;
            let diagnostics = String::from_utf8(diagnostics)?;
            assert_eq!(
                diagnostics.matches("Symbolic / SV ALT alleles at").count(),
                expected_occurrences
            );
            assert!(fs::read_to_string(output)?.contains("\"SYMALT\":2"));
        }
        Ok(())
    }

    #[test]
    fn unchecked_bcf_translation_errors_are_warned_and_skipped() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        write_vcf(&input, &["chr2\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1"])?;
        let mut diagnostics = Vec::new();
        run_with_diagnostics(args(&input, &output), &mut diagnostics)?;

        assert_eq!(
            fs::read_to_string(output)?,
            "{\"OVERLAP\":0,\"REFPADDING\":0,\"SYMALT\":0,\"UNCERTAINLENGTH\":0,\"diploid\":0,\"haploid\":0,\"male\":false,\"nonref\":0,\"polyploid\":0,\"records\":0,\"ref\":0}\n"
        );
        let diagnostics = String::from_utf8(diagnostics)?;
        assert_eq!(
            diagnostics,
            "[W::vcf_parse] Contig 'chr2' is not defined in the header. (Quick workaround: index the file with tabix.)\n[W] Record at chr2:1 will not translate into BCF. Check if the header is incomplete  (error code 1) -- all records like this are skipped.\n[W] Variants that will cause trouble when writing BCF: 1\n[I] Total VCF records:         0\n[I] Non-reference VCF records: 0\n"
        );
        Ok(())
    }

    #[test]
    fn checked_bcf_translation_errors_are_fatal() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        write_vcf(&input, &["chr2\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1"])?;
        let mut args = args(&input, &output);
        args.check_bcf_errors = true;
        let error = run_with_diagnostics(args, &mut Vec::new()).unwrap_err();
        assert!(error.to_string().contains(
            "Record at chr2:1 will not translate into BCF. Check if the header is incomplete (error code 1)."
        ));
        assert!(!output.exists());
        Ok(())
    }

    #[test]
    fn checked_bcf_translation_accepts_valid_records() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        write_vcf(&input, &["chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1"])?;
        let mut case_args = args(&input, &output);
        case_args.check_bcf_errors = true;
        run_with_diagnostics(case_args, &mut Vec::new())?;

        assert!(fs::read_to_string(output)?.contains("\"records\":1"));
        Ok(())
    }

    #[test]
    fn undefined_field_translation_errors_follow_bcf_mode() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let unchecked_output = directory.path().join("unchecked.json");
        write_vcf(
            &input,
            &[
                "chr1\t1\t.\tA\tC\t.\tLowX\tUNDECLARED=1\tGT:DP\t0/1:10",
                "chr1\t2\t.\tA\tG\t.\tLowX\tUNDECLARED=2\tGT:DP\t0/1:20",
            ],
        )?;

        let mut diagnostics = Vec::new();
        run_with_diagnostics(args(&input, &unchecked_output), &mut diagnostics)?;
        assert!(fs::read_to_string(unchecked_output)?.contains("\"records\":1"));
        assert_eq!(
            String::from_utf8(diagnostics)?,
            "[W::vcf_parse_filter] FILTER 'LowX' is not defined in the header\n[W::vcf_parse_info] INFO 'UNDECLARED' is not defined in the header, assuming Type=String\n[W::vcf_parse_format] FORMAT 'DP' at chr1:1 is not defined in the header, assuming Type=String\n[W] Record at chr1:1 will not translate into BCF. Check if the header is incomplete  (error code 2) -- all records like this are skipped.\n[W] Variants that will cause trouble when writing BCF: 1\n[I] Total VCF records:         1\n[I] Non-reference VCF records: 1\n"
        );

        let checked_output = directory.path().join("checked.json");
        let mut checked_args = args(&input, &checked_output);
        checked_args.check_bcf_errors = true;
        let error = run_with_diagnostics(checked_args, &mut Vec::new()).unwrap_err();
        assert!(error.to_string().contains("(error code 2)"));
        assert!(!checked_output.exists());
        Ok(())
    }

    #[test]
    fn invalid_gt_allele_index_is_fatal_for_both_bcf_modes() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        write_vcf(&input, &["chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t2/2"])?;

        for check_bcf_errors in [false, true] {
            let output = directory
                .path()
                .join(format!("check-{check_bcf_errors}.json"));
            let mut case_args = args(&input, &output);
            case_args.check_bcf_errors = check_bcf_errors;
            let error = run_with_diagnostics(case_args, &mut Vec::new()).unwrap_err();
            let message = error.to_string();
            assert!(message.contains(&output.display().to_string()));
            assert!(message.contains("Call with invalid genotype (non-existent allele) at chr1:1"));
            assert!(!output.exists());
        }
        Ok(())
    }

    #[test]
    fn declared_info_type_and_field_cardinality_are_not_translation_errors() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allelic depth\">\n##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihoods\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t1\t.\tA\tC\t.\tPASS\tDP=bad\tGT\t0/1\nchr1\t2\t.\tA\tG\t.\tPASS\tDP=1,2\tGT:AD:PL\t0/1:3,4:0,10\n",
        )?;
        run_with_diagnostics(args(&input, &output), &mut Vec::new())?;

        assert!(fs::read_to_string(output)?.contains("\"records\":2"));
        Ok(())
    }

    #[test]
    fn malformed_typed_format_value_is_silently_skipped() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:DP\t0/1:bad\nchr1\t2\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:10\n",
        )?;

        for check_bcf_errors in [false, true] {
            let output = directory
                .path()
                .join(format!("check-{check_bcf_errors}.json"));
            let mut case_args = args(&input, &output);
            case_args.check_bcf_errors = check_bcf_errors;
            let mut diagnostics = Vec::new();
            run_with_diagnostics(case_args, &mut diagnostics)?;

            assert!(fs::read_to_string(output)?.contains("\"records\":0"));
            let diagnostics = String::from_utf8(diagnostics)?;
            assert_eq!(
                diagnostics,
                "[W::vcf_parse_format] Extreme FORMAT/DP value encountered and set to missing at chr1:1\n[E::vcf_parse_format] Invalid character 'b' in 'DP' FORMAT field at chr1:1\n[I] Total VCF records:         0\n[I] Non-reference VCF records: 0\n"
            );
        }
        Ok(())
    }

    #[test]
    fn malformed_float_format_diagnostic_matches_htslib() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:AF\t0/1:bad\n",
        )?;
        let mut diagnostics = Vec::new();
        run_with_diagnostics(args(&input, &output), &mut diagnostics)?;

        assert_eq!(
            String::from_utf8(diagnostics)?,
            "[W::vcf_parse_format] Extreme FORMAT/AF value encountered at chr1:1\n[E::vcf_parse_format] Invalid character 'b' in 'AF' FORMAT field at chr1:1\n[I] Total VCF records:         0\n[I] Non-reference VCF records: 0\n"
        );
        Ok(())
    }

    #[test]
    fn malformed_numeric_suffix_emits_only_the_htslib_error() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:DP\t0/1:10x\n",
        )?;
        let mut diagnostics = Vec::new();
        run_with_diagnostics(args(&input, &output), &mut diagnostics)?;

        assert_eq!(
            String::from_utf8(diagnostics)?,
            "[E::vcf_parse_format] Invalid character 'x' in 'DP' FORMAT field at chr1:1\n[I] Total VCF records:         0\n[I] Non-reference VCF records: 0\n"
        );
        Ok(())
    }

    #[test]
    fn malformed_format_recovery_matches_the_legacy_reader() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\nchr1\t1\t.\tA\tC\t.\tPASS\t.\tGT:DP\t0/1:10\nchr1\t2\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:bad\nchr1\t3\t.\tA\tT\t.\tPASS\t.\tGT:DP\t0/1:20\nchr1\t4\t.\tA\tC\t.\tPASS\t.\tGT:DP\t0/1:bad\nchr1\t5\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:also_bad\nchr1\t6\t.\tA\tT\t.\tPASS\t.\tGT:DP\t0/1:30\n",
        )?;
        let mut diagnostics = Vec::new();
        run_with_diagnostics(args(&input, &output), &mut diagnostics)?;

        assert!(fs::read_to_string(output)?.contains("\"records\":2"));
        assert_eq!(
            String::from_utf8(diagnostics)?,
            "[W::vcf_parse_format] Extreme FORMAT/DP value encountered and set to missing at chr1:2\n[E::vcf_parse_format] Invalid character 'b' in 'DP' FORMAT field at chr1:2\n[E::vcf_parse_format] Invalid character 'b' in 'DP' FORMAT field at chr1:4\n[E::vcf_parse_format] Invalid character 'a' in 'DP' FORMAT field at chr1:5\n[I] Total VCF records:         2\n[I] Non-reference VCF records: 2\n"
        );
        Ok(())
    }

    #[test]
    fn extra_sample_columns_are_ignored_like_htslib() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        write_vcf(&input, &["chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1"])?;
        run_with_diagnostics(args(&input, &output), &mut Vec::new())?;

        assert_eq!(
            fs::read_to_string(output)?,
            "{\"OVERLAP\":0,\"REFPADDING\":0,\"SYMALT\":0,\"UNCERTAINLENGTH\":0,\"diploid\":0,\"haploid\":0,\"male\":false,\"nonref\":0,\"polyploid\":0,\"records\":1,\"ref\":1}\n"
        );
        Ok(())
    }

    #[test]
    fn missing_sample_columns_skip_the_record() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        fs::write(
            &input,
            "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\nchr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\n",
        )?;
        let mut diagnostics = Vec::new();
        run_with_diagnostics(args(&input, &output), &mut diagnostics)?;

        assert!(fs::read_to_string(output)?.contains("\"records\":0"));
        assert_eq!(
            String::from_utf8(diagnostics)?,
            "[E::vcf_parse_format] Number of columns at chr1:1 does not match the number of samples (1 vs 2)\n[I] Total VCF records:         0\n[I] Non-reference VCF records: 0\n"
        );
        Ok(())
    }

    #[test]
    fn male_inference_preserves_legacy_chromosome_spelling() -> Result<()> {
        let directory = tempdir()?;
        for (chrom, location, expected_male) in [
            ("X", None, true),
            ("chrX", None, true),
            ("chrx", None, true),
            ("x", None, false),
            ("CHRX", None, false),
            ("x", Some("x"), true),
        ] {
            let input = directory.path().join(format!("{chrom}.vcf"));
            let output = directory.path().join(format!("{chrom}-{location:?}.json"));
            fs::write(
                &input,
                format!(
                    "##fileformat=VCFv4.2\n##contig=<ID={chrom},length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n{chrom}\t1\t.\tA\tC\t.\tPASS\t.\tGT\t1\n"
                ),
            )?;
            let mut case_args = args(&input, &output);
            case_args.locations = location.map(str::to_string);
            run_with_diagnostics(case_args, &mut Vec::new())?;
            assert!(
                fs::read_to_string(output)?.contains(&format!("\"male\":{expected_male}")),
                "unexpected male inference for {chrom} with location {location:?}"
            );
        }
        Ok(())
    }
}
