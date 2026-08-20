//! `hap ftx` — feature-table extractor.
//!
//! Runs the lightweight filtering and optional `bcftools norm` equivalent
//! performed by legacy `ftx.py` before dispatching to a feature-table
//! emitter. This is deliberately separate from `hap pre`: the Python
//! `preprocessVCF` helper did not run hap.py's primitive decomposition or
//! genotype rewriting stages.

use crate::application::{FtxArgs, ValidatedFtxArgs};
use crate::domain::RawVcfRecord;
use crate::engines::partial_credit::RefVar;
use crate::{
    adapters::{fasta, vcf},
    engines::partial_credit,
    output::{FailureOperation, OutputTransaction, fail_operation},
};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, UNIX_EPOCH};

mod bam_depth;
mod common;
mod generic;
mod legacy_callers;
mod strelka_indel;
mod strelka_snv;

static SCRATCH_RUN_ID: AtomicU64 = AtomicU64::new(0);

/// The reference is an argument, and `--normalize` is the only ftx feature
/// that opens one.
fn normalize_reference(args: &FtxArgs) -> Result<Option<&str>> {
    if !args.normalize {
        return Ok(None);
    }
    args.reference
        .as_deref()
        .map(Some)
        .context("no reference file found for --normalize; pass --reference")
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
            let path = parent.join(format!("run-{}-{timestamp}-{id}", std::process::id()));
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
            "failed to allocate a unique ftx scratch directory under {}",
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
            format!("failed to remove scratch directory {}", self.path.display())
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

/// Where a scratch run directory is allocated. Honours `--scratch-prefix` when
/// present (a declared extension mirroring `germline`/`somatic`); otherwise
/// places scratch under the output directory, never `$TMPDIR`.
fn scratch_parent(args: &FtxArgs, output: &Path) -> PathBuf {
    if let Some(prefix) = args.scratch_prefix.as_deref() {
        return PathBuf::from(prefix);
    }
    output_parent(output).join(".hap_scratch")
}

pub(crate) fn run(args: ValidatedFtxArgs) -> Result<()> {
    let mut args = args.into_inner();
    normalize_reference(&args)?;
    let output = ftx_output_path(&args.output);
    // Scratch now lives inside the output directory, so refuse a missing
    // output parent up front rather than materialising it when scratch is
    // created. Mirrors the germline guard in application/compare.rs.
    validate_output_parent(&output)?;
    let inputs = std::iter::once(args.input.as_str())
        .chain(args.reference.as_deref())
        .chain(args.regions_bedfile.as_deref())
        .chain(args.targets_bedfile.as_deref())
        .chain(args.bams.iter().map(String::as_str))
        .map(PathBuf::from)
        .collect::<Vec<_>>();
    let transaction = OutputTransaction::files(&inputs, [&output])?;
    args.output = transaction
        .staged_file(&output)?
        .to_string_lossy()
        .into_owned();
    run_inner(args).map_err(|error| {
        anyhow::anyhow!(
            "failed to produce feature table {}: {error:#}",
            output.display()
        )
    })?;
    transaction.commit()
}

fn output_parent(output: &Path) -> &Path {
    output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."))
}

fn validate_output_parent(output: &Path) -> Result<()> {
    if !output_parent(output).exists() {
        bail!(
            "The output path does not exist. Please specify a valid output path and prefix using -o"
        );
    }
    Ok(())
}

fn ftx_output_path(output: &str) -> PathBuf {
    if output.ends_with(".csv") {
        PathBuf::from(output)
    } else {
        PathBuf::from(format!("{output}.csv"))
    }
}

fn run_inner(args: FtxArgs) -> Result<()> {
    let label = legacy_feature_label(&args.input, args.label.as_deref());

    // Legacy passes the reference only to `bcftools norm`; ordinary feature
    // extraction never opens or validates it. Preserve that lazy behavior so
    // an explicitly missing path is harmless unless normalization is enabled.
    let reference_sequences = match normalize_reference(&args)? {
        Some(path) => fasta::read_sequences(Path::new(path))?,
        None => BTreeMap::new(),
    };
    let reference_contigs: BTreeSet<String> = reference_sequences.keys().cloned().collect();

    // Keep staging scoped to a unique RAII scratch directory: VCF writers can
    // create index sidecars, and both files must disappear on success/error
    // unless --keep-scratch is set. Scratch honours --scratch-prefix and
    // otherwise lives inside the output directory, never $TMPDIR, so a caller's
    // disk accounting and cleanup reach it.
    let output_path = ftx_output_path(&args.output);
    let scratch = ScratchRun::create(&scratch_parent(&args, &output_path), args.keep_scratch)?;
    let temp_path = scratch.path().join("input.vcf.gz");
    let (headers, staged_records) =
        prepare_records(&args, &reference_sequences, &reference_contigs)?;
    vcf::write_raw_vcf(&temp_path, &headers, &staged_records)?;
    let staged_input = temp_path;

    let (headers, records) = vcf::load_raw_vcf(&staged_input)?;

    let bam_depths = bam_normalization_depths(&args.bams)?;
    let lines = emit_feature_table_with_depths(
        &args.features,
        &records,
        &headers,
        &label,
        (!bam_depths.is_empty()).then_some(&bam_depths),
    )?;

    let output = ftx_output_path(&args.output);
    fail_operation(FailureOperation::Writer, &output)?;
    fs::write(&output, format!("{}\n", lines.join("\n")))
        .with_context(|| format!("failed to write {}", output.display()))?;

    scratch.cleanup()?;
    Ok(())
}

fn legacy_feature_label(input: &str, label: Option<&str>) -> String {
    label
        .filter(|label| !label.is_empty())
        .map(str::to_owned)
        .unwrap_or_else(|| {
            Path::new(input)
                .file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .to_string()
        })
}

/// Emit a feature table from records already filtered by the caller.
/// `hap somatic` uses this narrow API after classification so standalone and
/// merged feature extraction share the byte-parity implementation.
pub(crate) fn emit_feature_table_with_depths(
    feature: &str,
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depths: Option<&BTreeMap<String, f64>>,
) -> Result<Vec<String>> {
    emit_feature_table_internal(feature, records, headers, label, depths, false)
}

pub(crate) fn emit_feature_table_for_somatic(
    feature: &str,
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depths: Option<&BTreeMap<String, f64>>,
) -> Result<Vec<String>> {
    emit_feature_table_internal(feature, records, headers, label, depths, true)
}

fn emit_feature_table_internal(
    feature: &str,
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depths: Option<&BTreeMap<String, f64>>,
    somatic_precision: bool,
) -> Result<Vec<String>> {
    let lines = match feature {
        "generic" => generic::emit(records, label),
        "admix.strelka.snv" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, ADMIX_STRELKA_TRUTH_FIELDS)
        }
        "admix.strelka.snv" => strelka_snv::emit_with_depths_precision(
            records,
            headers,
            label,
            depths,
            somatic_precision,
        ),
        "admix.strelka.indel" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, ADMIX_STRELKA_TRUTH_FIELDS)
        }
        "admix.strelka.indel" => strelka_indel::emit_with_depths(records, headers, label, depths),
        "hcc.strelka.snv" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, HCC_STRELKA_SNV_TRUTH_FIELDS)
        }
        "hcc.strelka.snv" => strelka_snv::emit_with_depths_precision(
            records,
            headers,
            label,
            depths,
            somatic_precision,
        ),
        "hcc.strelka.indel" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, HCC_STRELKA_INDEL_TRUTH_FIELDS)
        }
        "hcc.strelka.indel" => strelka_indel::emit_with_depths(records, headers, label, depths),
        "hcc.mutect.snv" | "hcc.mutect.indel" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, HCC_OTHER_SNV_TRUTH_FIELDS)
        }
        "hcc.mutect.snv" | "hcc.mutect.indel" => {
            legacy_callers::emit_mutect_with_depths(records, headers, label, depths)?
        }
        "hcc.varscan2.snv" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, HCC_OTHER_SNV_TRUTH_FIELDS)
        }
        "hcc.varscan2.indel" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, HCC_STRELKA_INDEL_TRUTH_FIELDS)
        }
        "hcc.varscan2.snv" => {
            legacy_callers::emit_varscan_with_depths(records, label, false, depths)
        }
        "hcc.varscan2.indel" => {
            legacy_callers::emit_varscan_with_depths(records, label, true, depths)
        }
        "hcc.pisces.snv" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, HCC_OTHER_SNV_TRUTH_FIELDS)
        }
        "hcc.pisces.indel" if matches!(label, "TP" | "FN") => {
            generic::emit_fields(records, label, HCC_STRELKA_INDEL_TRUTH_FIELDS)
        }
        "hcc.pisces.snv" | "hcc.pisces.indel" => {
            legacy_callers::emit_pisces_with_depths(records, headers, label, depths)?
        }
        other => bail!("unsupported --feature-table: {other}"),
    };
    Ok(lines)
}

/// Return legacy FTX per-contig normalization depths for one or more BAMs.
///
/// This crate-internal entry point keeps BAM parsing private to FTX while
/// allowing `hap somatic` to use the exact same averaging and tripling rules.
pub(crate) fn bam_normalization_depths(paths: &[String]) -> Result<BTreeMap<String, f64>> {
    bam_depth::normalization_depths(paths)
}

const HCC_STRELKA_SNV_TRUTH_FIELDS: &[&str] = &[
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUAL",
    "I.T_ALT_RATE",
    "I.DP_normal",
    "I.DP_tumor",
    "I.tag",
    "I.count",
];

const ADMIX_STRELKA_TRUTH_FIELDS: &[&str] =
    &["CHROM", "POS", "REF", "ALT", "I.editDistance", "S.2.GT"];

const HCC_STRELKA_INDEL_TRUTH_FIELDS: &[&str] = &[
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUAL",
    "S.1.VT",
    "I.T_ALT_RATE",
    "I.DP_normal",
    "I.DP_tumor",
    "I.tag",
    "I.count",
];

const HCC_OTHER_SNV_TRUTH_FIELDS: &[&str] = &[
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUAL",
    "I.MapQrange",
    "I.somatic",
    "I.filtered",
    "S.1.VT",
    "I.T_ALT_RATE",
    "I.DP_normal",
    "I.DP_tumor",
    "I.tag",
    "I.count",
];

/// Mirror the operations in `Tools.bcftools.preprocessVCF` used by legacy
/// `ftx.py`: PASS filtering, optional chr-prefix rewriting, location/region
/// selection, and optional `bcftools norm -f REF -c x -D` semantics.
fn prepare_records(
    args: &FtxArgs,
    reference_sequences: &BTreeMap<String, String>,
    _reference_contigs: &BTreeSet<String>,
) -> Result<(Vec<String>, Vec<RawVcfRecord>)> {
    let (headers, records) = vcf::load_raw_vcf(Path::new(&args.input))?;
    // `--fix-chr` rewrites VCF records before selection, but legacy passes
    // location and BED contigs to bcftools literally.
    let literal_contigs = BTreeSet::new();
    let locations = args
        .location
        .as_deref()
        .map(|value| vcf::parse_locations(value, &literal_contigs))
        .transpose()?;
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

    let mut output = Vec::with_capacity(records.len());
    let mut seen = HashSet::new();
    for mut record in records {
        if crate::application::preprocess::calls_non_ref_allele(&record) {
            continue;
        }
        if args.fixchr {
            record.chrom = legacy_fix_chrom(&record.chrom);
        }
        if !args.include_nonpass && !record.is_pass() {
            continue;
        }
        let effective_end = record.effective_end_pos(Path::new(&args.input))?;
        if !vcf::matches_interval_filters(
            &record.chrom,
            record.pos,
            effective_end,
            regions.as_deref(),
            targets.as_deref(),
            locations.as_deref(),
        ) {
            continue;
        }

        if args.normalize {
            let Some(reference) = reference_sequences.get(&record.chrom) else {
                bail!(
                    "cannot normalize {}:{}: contig is absent from the reference",
                    record.chrom,
                    record.pos
                );
            };
            // `bcftools norm -c x` excludes reference-mismatch records.
            if !record_reference_matches(&record, reference.as_bytes()) {
                continue;
            }
            normalize_record(&mut record, reference.as_bytes());

            // Legacy uses deprecated `-D` (now `-d exact`). Bcftools defines
            // identity by normalized locus + alleles, not by auxiliary VCF
            // columns such as ID, QUAL, or INFO; the first record wins.
            let duplicate_key = (
                record.chrom.clone(),
                record.pos,
                record.ref_allele.clone(),
                record.alt_allele.clone(),
            );
            if !seen.insert(duplicate_key) {
                continue;
            }
        }
        output.push(record);
    }

    Ok((headers, output))
}

/// Legacy's two Perl substitutions prefix a leading numeric/X/Y/M contig,
/// then collapse `chrMT` to `chrM`. Unlike `vcf::normalize_chrom`, this does
/// not depend on which names happen to be present in the reference.
fn legacy_fix_chrom(chrom: &str) -> String {
    let prefixed = match chrom.as_bytes().first() {
        Some(first) if first.is_ascii_digit() || matches!(*first, b'X' | b'Y' | b'M') => {
            format!("chr{chrom}")
        }
        _ => chrom.to_string(),
    };
    prefixed.replacen("chrMT", "chrM", 1)
}

fn record_reference_matches(record: &RawVcfRecord, reference: &[u8]) -> bool {
    let start = record.pos.saturating_sub(1);
    let end = start.saturating_add(record.ref_allele.len());
    reference
        .get(start..end)
        .is_some_and(|observed| observed.eq_ignore_ascii_case(record.ref_allele.as_bytes()))
}

/// Left-align and minimize all concrete alleles while preserving a single
/// multi-allelic record. Normalizing each ALT independently and then padding
/// to their common reference span matches bcftools' representation without
/// splitting the record (`norm -f`, without `-m`).
fn normalize_record(record: &mut RawVcfRecord, reference: &[u8]) {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    if alts
        .iter()
        .any(|alt| alt.is_empty() || *alt == "." || alt.starts_with('<') || *alt == "*")
    {
        return;
    }

    let end = record.end_pos();
    let mut normalized = Vec::with_capacity(alts.len());
    for alt in alts {
        let mut variant = RefVar {
            start: record.pos,
            end,
            alt: alt.to_string(),
        };
        // Normalize without VCF padding first. The ref-padding mode used by
        // hap.py's C++ preprocessor reserves an extra left position and stops
        // one base too early for bcftools normalization in homopolymers.
        partial_credit::left_shift(reference, &mut variant, 1, false);
        normalized.push(pad_normalized_allele(variant, reference));
    }

    let common_start = normalized
        .iter()
        .map(|variant| variant.start)
        .min()
        .unwrap_or(record.pos);
    let common_end = normalized
        .iter()
        .map(|variant| variant.end)
        .max()
        .unwrap_or(end);
    if common_start == 0 || common_end < common_start || common_end > reference.len() {
        return;
    }

    let common_ref = &reference[common_start - 1..common_end];
    let normalized_alts: Vec<String> = normalized
        .into_iter()
        .map(|variant| {
            let mut allele = Vec::new();
            allele.extend_from_slice(&reference[common_start - 1..variant.start - 1]);
            allele.extend_from_slice(variant.alt.as_bytes());
            allele.extend_from_slice(&reference[variant.end..common_end]);
            String::from_utf8_lossy(&allele).to_ascii_uppercase()
        })
        .collect();

    record.pos = common_start;
    record.ref_allele = String::from_utf8_lossy(common_ref).to_ascii_uppercase();
    record.alt_allele = normalized_alts.join(",");
}

fn pad_normalized_allele(mut variant: RefVar, reference: &[u8]) -> RefVar {
    let ref_len = variant.end as i64 - variant.start as i64 + 1;
    if ref_len <= 0 && !variant.alt.is_empty() {
        // Pure insertion: prepend the base immediately left of the insertion
        // boundary, yielding a valid VCF allele.
        if variant.start > 1 {
            let anchor_pos = variant.start - 1;
            let anchor = reference[anchor_pos - 1].to_ascii_uppercase() as char;
            variant.start = anchor_pos;
            variant.end = anchor_pos;
            variant.alt.insert(0, anchor);
        }
    } else if ref_len > 0 && variant.alt.is_empty() {
        // Pure deletion: prefer a left anchor; only use a right anchor when
        // the deletion begins at the first reference base.
        if variant.start > 1 {
            let anchor_pos = variant.start - 1;
            let anchor = reference[anchor_pos - 1].to_ascii_uppercase() as char;
            variant.start = anchor_pos;
            variant.alt.push(anchor);
        } else if variant.end < reference.len() {
            let anchor_pos = variant.end + 1;
            let anchor = reference[anchor_pos - 1].to_ascii_uppercase() as char;
            variant.end = anchor_pos;
            variant.alt.push(anchor);
        }
    }
    variant
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::{FailureOperation, set_failure_operation};
    use std::thread;

    fn args(input: &Path, reference: &Path) -> FtxArgs {
        crate::application::FtxArgs {
            input: input.display().to_string(),
            output: "unused".to_string(),
            location: None,
            regions_bedfile: None,
            targets_bedfile: None,
            include_nonpass: false,
            features: "generic".to_string(),
            label: None,
            bams: Vec::new(),
            reference: Some(reference.display().to_string()),
            normalize: false,
            fixchr: false,
            scratch_prefix: None,
            keep_scratch: false,
        }
    }

    fn run(args: FtxArgs) -> Result<()> {
        super::run(args.validated()?)
    }

    fn fixture(contents: &str, reference: &str) -> (tempfile::TempDir, PathBuf, PathBuf) {
        let scratch = tempfile::tempdir().unwrap();
        let input = scratch.path().join("input.vcf");
        let fasta = scratch.path().join("ref.fa");
        fs::write(&input, contents).unwrap();
        fs::write(&fasta, reference).unwrap();
        (scratch, input, fasta)
    }

    #[test]
    fn scratch_run_is_created_under_the_given_parent() {
        let parent = tempfile::tempdir().unwrap();
        let scratch = ScratchRun::create(parent.path(), false).unwrap();
        assert!(scratch.path().starts_with(parent.path()));
        let run_path = scratch.path().to_path_buf();
        assert!(run_path.is_dir());
        scratch.cleanup().unwrap();
        assert!(!run_path.exists());
    }

    #[test]
    fn keep_scratch_retains_the_run_directory() {
        let parent = tempfile::tempdir().unwrap();
        let scratch = ScratchRun::create(parent.path(), true).unwrap();
        let run_path = scratch.path().to_path_buf();
        scratch.cleanup().unwrap();
        assert!(run_path.is_dir(), "--keep-scratch must retain the run dir");
    }

    #[test]
    fn scratch_prefix_is_honoured_over_the_output_directory() {
        let root = tempfile::tempdir().unwrap();
        let prefix = root.path().join("explicit");
        let mut arguments = args(Path::new("in.vcf"), Path::new("ref.fa"));
        arguments.output = root.path().join("features").display().to_string();
        arguments.scratch_prefix = Some(prefix.display().to_string());
        let parent = scratch_parent(&arguments, &ftx_output_path(&arguments.output));
        assert_eq!(parent, prefix);
    }

    #[test]
    fn concurrent_scratch_runs_are_unique_and_cleanup_sidecars() {
        let parent = tempfile::tempdir().unwrap();
        let first_parent = parent.path().to_path_buf();
        let second_parent = parent.path().to_path_buf();
        let first = thread::spawn(move || ScratchRun::create(&first_parent, false));
        let second = thread::spawn(move || ScratchRun::create(&second_parent, false));
        let first = first.join().expect("first thread panicked").unwrap();
        let second = second.join().expect("second thread panicked").unwrap();
        assert_ne!(first.path(), second.path());

        let first_path = first.path().to_path_buf();
        let second_path = second.path().to_path_buf();
        fs::write(first.path().join("input.vcf.gz"), b"vcf").unwrap();
        fs::write(first.path().join("input.vcf.gz.tbi"), b"index").unwrap();
        fs::write(second.path().join("input.vcf.gz"), b"vcf").unwrap();
        drop(first);
        drop(second);

        assert!(!first_path.exists());
        assert!(!second_path.exists());
    }

    #[test]
    fn injected_csv_writer_preserves_generation_and_cleanup() -> Result<()> {
        let (scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t60\tPASS\t.\n",
            ">chr1\nAAAA\n",
        );
        let output = scratch.path().join("features.csv");
        fs::write(&output, "old-csv")?;
        let mut arguments = args(&input, &reference);
        arguments.output = output.to_string_lossy().into_owned();

        set_failure_operation(Some(FailureOperation::Writer));
        let error = run(arguments).expect_err("injected CSV writer operation must fail");
        set_failure_operation(None);

        assert!(error.to_string().contains(&output.display().to_string()));
        assert_eq!(fs::read_to_string(&output)?, "old-csv");
        assert!(fs::read_dir(scratch.path())?.all(|entry| {
            !entry
                .expect("scratch entry must be readable")
                .file_name()
                .to_string_lossy()
                .contains("hap-rs")
        }));
        Ok(())
    }

    #[test]
    fn empty_feature_label_falls_back_to_input_basename() {
        assert_eq!(
            legacy_feature_label("fixtures/options.vcf", None),
            "options.vcf"
        );
        assert_eq!(
            legacy_feature_label("fixtures/options.vcf", Some("")),
            "options.vcf"
        );
        assert_eq!(
            legacy_feature_label("fixtures/options.vcf", Some(" ")),
            " ",
            "non-empty labels remain truthy in Python"
        );
    }

    #[test]
    fn fix_chr_is_opt_in_and_matches_legacy_perl_rewrite() {
        assert_eq!(legacy_fix_chrom("1"), "chr1");
        assert_eq!(legacy_fix_chrom("X"), "chrX");
        assert_eq!(legacy_fix_chrom("MT"), "chrM");
        assert_eq!(legacy_fix_chrom("chr1"), "chr1");

        let (_scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\t.\tA\tC\t60\tPASS\t.\n",
            ">chr1\nAAAA\n",
        );
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();

        let (_, records) =
            prepare_records(&args(&input, &reference), &sequences, &contigs).unwrap();
        assert_eq!(records[0].chrom, "1");

        let mut fixed = args(&input, &reference);
        fixed.fixchr = true;
        let (_, records) = prepare_records(&fixed, &sequences, &contigs).unwrap();
        assert_eq!(records[0].chrom, "chr1");
    }

    #[test]
    fn legacy_only_called_final_non_ref_alleles_are_removed_by_ftx() {
        let (_scratch, input, reference) = fixture(
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n",
                "chr1\t1\tcalled\tA\tC,<NON_REF>\t60\tPASS\t.\tGT\t0/1\t0|2\n",
            ),
            ">chr1\nAAAA\n",
        );
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();

        let (_, records) =
            prepare_records(&args(&input, &reference), &sequences, &contigs).unwrap();

        assert!(records.is_empty());
    }

    #[test]
    fn normative_uncalled_or_nonfinal_non_ref_alleles_are_retained_by_ftx() {
        let (_scratch, input, reference) = fixture(
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR\n",
                "chr1\t2\tuncalled\tA\tC,<NON_REF>\t60\tPASS\t.\tGT\t0/1\t0/0\n",
                "chr1\t3\tnot-final\tA\t<NON_REF>,C\t60\tPASS\t.\tGT\t0/1\t0/0\n",
            ),
            ">chr1\nAAAA\n",
        );
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();

        let (_, records) =
            prepare_records(&args(&input, &reference), &sequences, &contigs).unwrap();

        assert_eq!(
            records
                .iter()
                .map(|record| record.id.as_str())
                .collect::<Vec<_>>(),
            ["uncalled", "not-final"]
        );
    }

    #[test]
    fn region_uses_symbolic_end_while_target_uses_start_position() {
        let (scratch, input, reference) = fixture(
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "chr1\t3\t.\tC\t<DEL>\t.\tPASS\tEND=5\n",
            ),
            ">chr1\nAACCGGTTAA\n",
        );
        let boundary = scratch.path().join("boundary.bed");
        fs::write(&boundary, "chr1\t4\t5\n").unwrap();
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();

        let mut region_args = args(&input, &reference);
        region_args.regions_bedfile = Some(boundary.display().to_string());
        let (_, region_records) = prepare_records(&region_args, &sequences, &contigs).unwrap();

        let mut target_args = args(&input, &reference);
        target_args.targets_bedfile = Some(boundary.display().to_string());
        let (_, target_records) = prepare_records(&target_args, &sequences, &contigs).unwrap();

        assert_eq!(region_records.len(), 1, "-R must use INFO/END");
        assert!(target_records.is_empty(), "-T must use POS");
    }

    #[test]
    fn normalize_flag_left_aligns_homopolymer_indel() {
        let (_scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrH\t2\t.\tA\tAA\t60\tPASS\t.\n",
            ">chrH\nCAAAAA\n",
        );
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();

        let (_, untouched) =
            prepare_records(&args(&input, &reference), &sequences, &contigs).unwrap();
        assert_eq!(
            (
                untouched[0].pos,
                untouched[0].ref_allele.as_str(),
                untouched[0].alt_allele.as_str()
            ),
            (2, "A", "AA")
        );

        let mut normalized_args = args(&input, &reference);
        normalized_args.normalize = true;
        let (_, normalized) = prepare_records(&normalized_args, &sequences, &contigs).unwrap();
        assert_eq!(
            (
                normalized[0].pos,
                normalized[0].ref_allele.as_str(),
                normalized[0].alt_allele.as_str()
            ),
            (1, "C", "CA")
        );
    }

    #[test]
    fn normalize_matches_bcftools_multi_allelic_padding_and_deduplicates() {
        let (_scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrH\t2\tfirst\tA\tAA,AAA\t60\tPASS\t.\nchrH\t2\tsecond\tA\tAA,AAA\t50\tPASS\tX=1\nchrH\t3\t.\tC\tT\t60\tPASS\t.\n",
            ">chrH\nCAAAAA\n",
        );
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();
        let mut normalized_args = args(&input, &reference);
        normalized_args.normalize = true;

        let (_, records) = prepare_records(&normalized_args, &sequences, &contigs).unwrap();
        assert_eq!(
            records.len(),
            1,
            "duplicate removed and REF mismatch excluded: {:?}",
            records
                .iter()
                .map(RawVcfRecord::to_line)
                .collect::<Vec<_>>()
        );
        assert_eq!(records[0].pos, 1);
        assert_eq!(records[0].ref_allele, "C");
        assert_eq!(records[0].alt_allele, "CA,CAA");
    }

    #[test]
    fn target_is_not_reapplied_after_normalization_moves_a_record() {
        let (scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrH\t2\t.\tA\tAA\t60\tPASS\t.\n",
            ">chrH\nCAAAAA\n",
        );
        let target = scratch.path().join("target.bed");
        fs::write(&target, "chrH\t1\t2\n").unwrap();
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();
        let mut normalized_args = args(&input, &reference);
        normalized_args.normalize = true;
        normalized_args.targets_bedfile = Some(target.display().to_string());

        let (_, records) = prepare_records(&normalized_args, &sequences, &contigs).unwrap();

        assert_eq!(records.len(), 1, "-T selected the pre-normalized POS");
        assert_eq!(records[0].pos, 1, "normalization then moved the record");
        assert_eq!(generic::emit(&records, "tag").len(), 2);
    }

    #[test]
    fn fix_chr_does_not_rewrite_region_or_target_bed_contigs() {
        let (scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t2\t.\tA\tC\t60\tPASS\t.\n",
            ">chr1\nAAAA\n",
        );
        let selector = scratch.path().join("selector.bed");
        fs::write(&selector, "1\t0\t4\n").unwrap();
        let sequences = fasta::read_sequences(&reference).unwrap();
        let contigs = sequences.keys().cloned().collect();

        for use_regions in [true, false] {
            let mut selected_args = args(&input, &reference);
            selected_args.fixchr = true;
            selected_args.normalize = true;
            if use_regions {
                selected_args.regions_bedfile = Some(selector.display().to_string());
            } else {
                selected_args.targets_bedfile = Some(selector.display().to_string());
            }

            let (_, records) = prepare_records(&selected_args, &sequences, &contigs).unwrap();

            assert!(
                records.is_empty(),
                "rewritten record contig must not match the literal BED contig"
            );
        }
    }

    #[test]
    fn bam_depths_override_strelka_header_normalization() {
        let record = RawVcfRecord {
            chrom: "chr1".to_string(),
            pos: 7,
            id: ".".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "C".to_string(),
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            info: "NT=ref;QSS_NT=10".to_string(),
            format: Some("DP:FDP:SDP:AU:CU:GU:TU".to_string()),
            samples: vec![
                "10:0:0:10,10:0,0:0,0:0,0".to_string(),
                "20:0:0:0,0:20,20:0,0:0,0".to_string(),
            ],
            mixed_edit_primitive: false,
            primitive_identity: None,
        };
        let headers = vec!["##MaxDepth_chr1=100".to_string()];
        let depths = BTreeMap::from([("chr1".to_string(), 50.0)]);

        let lines = emit_feature_table_with_depths(
            "hcc.strelka.snv",
            &[record],
            &headers,
            "FP",
            Some(&depths),
        )
        .unwrap();
        let cells: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(cells[18], "0.2");
        assert_eq!(cells[19], "0.4");
    }

    #[test]
    fn missing_explicit_reference_is_ignored_without_normalization() {
        let (scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tA\tC\t60\tPASS\t.\n",
            ">chr1\nAAAA\n",
        );
        let output = scratch.path().join("features");
        let mut run_args = args(&input, &reference);
        run_args.output = output.display().to_string();
        run_args.reference = Some(scratch.path().join("missing.fa").display().to_string());

        run(run_args).unwrap();

        assert!(output.with_extension("csv").is_file());
    }

    #[test]
    fn missing_output_parent_is_not_created() {
        let (scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tA\tC\t60\tPASS\t.\n",
            ">chr1\nAAAA\n",
        );
        let missing_parent = scratch.path().join("missing");
        let output = missing_parent.join("features");
        let mut run_args = args(&input, &reference);
        run_args.output = output.display().to_string();

        assert!(run(run_args).is_err());
        assert!(!missing_parent.exists());
        assert!(!output.with_extension("csv").exists());
    }

    #[test]
    fn missing_explicit_reference_errors_when_normalizing() {
        let (scratch, input, reference) = fixture(
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t2\t.\tA\tC\t60\tPASS\t.\n",
            ">chr1\nAAAA\n",
        );
        let mut run_args = args(&input, &reference);
        run_args.output = scratch.path().join("features").display().to_string();
        run_args.reference = Some(scratch.path().join("missing.fa").display().to_string());
        run_args.normalize = true;

        let error = run(run_args).unwrap_err().to_string();

        assert!(error.contains("failed to read FASTA"), "{error}");
    }
}
