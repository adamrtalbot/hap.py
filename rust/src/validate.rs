use crate::cli::ValidateArgs;
use crate::{fasta, vcf};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::Path;

const WARNING_REFPADDING: usize = 0;
const WARNING_OVERLAP: usize = 1;
const WARNING_SYMALT: usize = 2;
const WARNING_UNCERTAIN_LENGTH: usize = 3;
const WARNING_COUNT: usize = 4;

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

pub fn run(args: ValidateArgs) -> Result<()> {
    let stderr = std::io::stderr();
    run_with_diagnostics(args, &mut stderr.lock())
}

fn run_with_diagnostics<W: Write>(args: ValidateArgs, diagnostics: &mut W) -> Result<()> {
    if args.check_bcf_errors {
        bail!(
            "--check-bcf-errors true is unsupported: the built-in VCF reader does not perform htslib BCF translation checks"
        );
    }

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

    let (headers, records) = vcf::load_raw_vcf(Path::new(&args.input))?;
    ensure_sample_header(&headers)?;
    let reference_sequences = args
        .reference
        .as_ref()
        .map(|path| fasta::read_sequences(Path::new(path)))
        .transpose()?;

    let mut counts = ValidationCounts::default();
    let mut previous = PreviousRecord::default();
    let mut errors = Vec::new();

    for record in records {
        let chrom = if reference_contigs.is_empty() {
            record.chrom.clone()
        } else {
            vcf::normalize_chrom(&record.chrom, &reference_contigs)
        };
        let record = vcf::RawVcfRecord { chrom, ..record };

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
        update_record_counts(&record, &genotype, &mut counts);

        if let Some(reference_sequences) = &reference_sequences
            && let Some(reason) = validate_record(&record, reference_sequences)
        {
            errors.push(format!(
                "{}\t{}\t{}\t{}",
                record.chrom,
                record.pos.saturating_sub(1),
                record.end_pos(),
                reason
            ));
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

    if let Some(path) = &args.errors_bed {
        fs::write(
            path,
            if errors.is_empty() {
                String::new()
            } else {
                format!("{}\n", errors.join("\n"))
            },
        )
        .with_context(|| format!("failed to write {path}"))?;
    }

    Ok(())
}

fn ensure_sample_header(headers: &[String]) -> Result<()> {
    let sample_count = headers
        .iter()
        .find(|line| line.starts_with("#CHROM\t"))
        .map(|line| line.split('\t').count().saturating_sub(9))
        .unwrap_or(0);
    if sample_count == 0 {
        bail!("input VCF has no samples; legacy vcfcheck requires at least one sample");
    }
    Ok(())
}

fn inspect_genotypes(record: &vcf::RawVcfRecord) -> Result<RecordGenotypes> {
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
            let alt = alts.get(index - 1).ok_or_else(|| {
                anyhow::anyhow!(
                    "call with invalid genotype allele {index} at {}:{}",
                    record.chrom,
                    record.pos
                )
            })?;
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
fn update_warning_counts<W: Write>(
    record: &vcf::RawVcfRecord,
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

fn strict_ref_count(record: &vcf::RawVcfRecord, sample: usize, strict: bool) -> usize {
    if !strict {
        return 0;
    }
    record
        .sample_map(sample)
        .get("GT")
        .map(|gt| gt.split(['/', '|']).filter(|allele| *allele == "0").count())
        .unwrap_or(0)
}

fn reference_padding(record: &vcf::RawVcfRecord) -> usize {
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

fn record_end(record: &vcf::RawVcfRecord) -> usize {
    record
        .info
        .split(';')
        .find_map(|field| field.strip_prefix("END=")?.parse::<usize>().ok())
        .map(|end| end.saturating_sub(1))
        .unwrap_or_else(|| record.end_pos().saturating_sub(1))
}

fn emit_warning<W: Write>(
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
    record: &vcf::RawVcfRecord,
    genotype: &RecordGenotypes,
    counts: &mut ValidationCounts,
) {
    counts.ref_records += usize::from(genotype.any_ref);
    counts.nonref_records += usize::from(genotype.any_nonref);
    counts.haploid += usize::from(genotype.any_haploid);
    counts.diploid += usize::from(genotype.any_diploid);
    counts.polyploid += usize::from(genotype.any_polyploid);
    if record.chrom.eq_ignore_ascii_case("x") || record.chrom.eq_ignore_ascii_case("chrx") {
        counts.haploid_x |= genotype.any_haploid;
        counts.diploid_x |= genotype.any_diploid || genotype.any_polyploid;
    }
}

fn write_warning_summaries<W: Write>(counts: &ValidationCounts, diagnostics: &mut W) -> Result<()> {
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
    record: &vcf::RawVcfRecord,
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
        let mut text = "##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n".to_string();
        for record in records {
            text.push_str(record);
            text.push('\n');
        }
        fs::write(path, text)?;
        Ok(())
    }

    fn args(input: &Path, output: &Path) -> ValidateArgs {
        ValidateArgs {
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
    fn bcf_translation_check_is_rejected_instead_of_ignored() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("check.json");
        write_vcf(&input, &["chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0/1"])?;
        let mut args = args(&input, &output);
        args.check_bcf_errors = true;
        let error = run_with_diagnostics(args, &mut Vec::new()).unwrap_err();
        assert!(error.to_string().contains("unsupported"));
        assert!(!output.exists());
        Ok(())
    }
}
