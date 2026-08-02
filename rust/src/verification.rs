use crate::cli::{CompareArgs, PreprocessArgs, QuantifyArgs, SomaticArgs, ValidateArgs};
use crate::fixtures::{
    COMPARE_FIXTURE_CASES, CompareFixtureCase, PREPROCESS_FIXTURE_CASES, PreprocessFixtureCase,
    QUANTIFY_FIXTURE_CASES, QuantifyFixtureCase, SOMATIC_FIXTURE_CASES, SomaticFixtureCase,
    VALIDATE_FIXTURE_CASES, ValidateFixtureCase, ValidateFixtureKind, repo_root,
};
use crate::{compare, preprocess, quantify, somatic, validate};
use anyhow::{Context, Result, bail};
use clap::ValueEnum;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

pub const LEGACY_IMAGE: &str = "community.wave.seqera.io/library/hap.py_rtg-tools@sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6";

pub fn bless_legacy_outputs(image: &str) -> Result<()> {
    for case in COMPARE_FIXTURE_CASES {
        let output_dir = run_legacy_compare_case(case, image, "legacy-bless")?;
        save_expected_compare_outputs(case, &output_dir)?;
    }
    for case in PREPROCESS_FIXTURE_CASES {
        let output_dir = run_legacy_preprocess_case(case, image, "legacy-bless")?;
        save_expected_preprocess_outputs(case, &output_dir)?;
    }
    for case in SOMATIC_FIXTURE_CASES {
        let output_dir = run_legacy_somatic_case(case, image, "legacy-bless")?;
        save_expected_somatic_outputs(case, &output_dir)?;
    }
    for case in VALIDATE_FIXTURE_CASES {
        if let ValidateFixtureKind::SummaryJson = case.kind {
            let output_dir = run_legacy_validate_case(case, image, "legacy-bless")?;
            save_expected_validate_outputs(case, &output_dir)?;
        }
    }
    Ok(())
}

pub fn verify_legacy_outputs(image: &str) -> Result<()> {
    for case in COMPARE_FIXTURE_CASES {
        let output_dir = run_legacy_compare_case(case, image, "legacy-check")?;
        compare_compare_outputs(case, &output_dir)
            .with_context(|| format!("legacy compare verification failed for {}", case.id))?;
    }
    for case in PREPROCESS_FIXTURE_CASES {
        let output_dir = run_legacy_preprocess_case(case, image, "legacy-check")?;
        compare_preprocess_outputs(case, &output_dir)
            .with_context(|| format!("legacy preprocess verification failed for {}", case.id))?;
    }
    for case in SOMATIC_FIXTURE_CASES {
        let output_dir = run_legacy_somatic_case(case, image, "legacy-check")?;
        compare_somatic_outputs(case, &output_dir)
            .with_context(|| format!("legacy somatic verification failed for {}", case.id))?;
    }
    for case in VALIDATE_FIXTURE_CASES {
        if let ValidateFixtureKind::SummaryJson = case.kind {
            let output_dir = run_legacy_validate_case(case, image, "legacy-check")?;
            compare_validate_outputs(case, &output_dir)
                .with_context(|| format!("legacy validate verification failed for {}", case.id))?;
        }
    }
    Ok(())
}

pub fn verify_rust_outputs() -> Result<()> {
    for case in COMPARE_FIXTURE_CASES {
        let output_dir = run_rust_compare_case(case)?;
        compare_compare_outputs(case, &output_dir)
            .with_context(|| format!("rust compare verification failed for {}", case.id))?;
    }
    for case in PREPROCESS_FIXTURE_CASES {
        let output_dir = run_rust_preprocess_case(case)?;
        compare_preprocess_outputs(case, &output_dir)
            .with_context(|| format!("rust preprocess verification failed for {}", case.id))?;
    }
    for case in QUANTIFY_FIXTURE_CASES {
        let output_dir = run_rust_quantify_case(case)?;
        compare_quantify_outputs(case, &output_dir)
            .with_context(|| format!("rust quantify verification failed for {}", case.id))?;
    }
    for case in SOMATIC_FIXTURE_CASES {
        let output_dir = run_rust_somatic_case(case)?;
        compare_somatic_outputs(case, &output_dir)
            .with_context(|| format!("rust somatic verification failed for {}", case.id))?;
    }
    for case in VALIDATE_FIXTURE_CASES {
        let output_dir = run_rust_validate_case(case)?;
        compare_validate_outputs(case, &output_dir)
            .with_context(|| format!("rust validate verification failed for {}", case.id))?;
    }
    Ok(())
}

pub fn verify_realworld_outputs(image: &str) -> Result<()> {
    for case in realworld_cases() {
        let workdir = materialize_realworld_inputs(&case)?;
        run_legacy_realworld_case(&case, image, &workdir)?;
        run_rust_realworld_case(&case, &workdir)?;
        compare_output_trees(&workdir.join("legacy_work"), &workdir.join("rust_work"))
            .with_context(|| format!("real-world parity failed for {}", case.id))?;
    }
    Ok(())
}

fn run_rust_compare_case(case: &CompareFixtureCase) -> Result<PathBuf> {
    let output_dir = prepare_output_dir("rust-check", case.id)?;
    let mut args = CompareArgs::with_paths(
        case.truth_path().display().to_string(),
        case.query_path().display().to_string(),
        case.reference_path().display().to_string(),
        output_dir.join("oracle").display().to_string(),
    );
    args.pass_only = case.pass_only;
    args.regions_bedfile = case
        .restrict_bed_path()
        .map(|path| path.display().to_string());
    args.fp_bedfile = case.fp_bed_path().map(|path| path.display().to_string());
    compare::run(args)?;
    Ok(output_dir)
}

fn run_rust_preprocess_case(case: &PreprocessFixtureCase) -> Result<PathBuf> {
    let output_dir = prepare_output_dir("rust-check", case.id)?;
    preprocess::run(PreprocessArgs {
        input: case.input_path().display().to_string(),
        output: output_dir.join(case.output).display().to_string(),
        version: false,
        reference: Some(case.reference_path().display().to_string()),
        locations: case.locations.map(str::to_string),
        pass_only: case.pass_only,
        filters_only: None,
        regions_bedfile: case
            .regions_bed_path()
            .map(|path| path.display().to_string()),
        targets_bedfile: None,
        fixchr: case.fixchr,
        no_fixchr: matches!(case.fixchr, Some(false)),
        somatic: case.somatic,
        set_gt: case.set_gt,
        filter_nonref: case.filter_nonref,
        convert_gvcf_to_vcf: case.convert_gvcf_to_vcf,
        bcf: false,
        bcftools_norm: false,
        leftshift: true,
        no_leftshift: false,
        decompose: true,
        no_decompose: false,
        gender: crate::cli::PreprocessGender::Auto,
        window_size: 10_000,
        threads: None,
        logfile: None,
        verbose: false,
        quiet: false,
        force_interactive: false,
    })?;
    Ok(output_dir)
}

fn run_rust_quantify_case(case: &QuantifyFixtureCase) -> Result<PathBuf> {
    let output_dir = prepare_output_dir("rust-check", case.id)?;
    quantify::run(QuantifyArgs {
        input_vcf: case.input_vcf_path().display().to_string(),
        report_prefix: output_dir.join("oracle").display().to_string(),
        reference: case.reference_path().display().to_string(),
        // This fixture already carries finalized per-sample BD/BK fields.
        // Select the GA4GH quantifier explicitly; legacy XCMP mode instead
        // re-derives decisions from record-level INFO/type.
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
        do_roc: true,
        roc_regions: vec!["*".to_string()],
        roc_filter: None,
        roc_delta: 0.5,
        ci_alpha: 0.0,
        no_json: true,
    })?;
    Ok(output_dir)
}

fn run_rust_somatic_case(case: &SomaticFixtureCase) -> Result<PathBuf> {
    let output_dir = prepare_output_dir("rust-check", case.id)?;
    somatic::run(SomaticArgs {
        truth: case.truth_path().display().to_string(),
        query: case.query_path().display().to_string(),
        output: output_dir.join("oracle").display().to_string(),
        reference: case.reference_path().display().to_string(),
        location: None,
        regions_bedfile: None,
        targets_bedfile: None,
        fp_bedfile: case.fp_bed_path().map(|path| path.display().to_string()),
        ambiguous_beds: Vec::new(),
        ambi_fp: false,
        no_ambi_fp: false,
        count_unk: case.count_unk,
        no_count_unk: false,
        explain_ambiguous: false,
        include_nonpass: case.include_nonpass,
        fp_region_size: None,
        feature_table: None,
        happy_stats: false,
        bams: Vec::new(),
        normalize_truth: false,
        normalize_query: false,
        normalize_all: false,
        fixchr_truth: None,
        fixchr_query: None,
        fix_chr_truth: None,
        fix_chr_query: None,
        no_fixchr_truth: false,
        no_fixchr_query: false,
        no_order_check: false,
        roc: None,
        af_strat: false,
        af_strat_binsize: "0.2".to_string(),
        af_strat_truth: "I.T_ALT_RATE".to_string(),
        af_strat_query: "T_AF".to_string(),
        count_filtered_fn: false,
        ci_level: 0.95,
        scratch_prefix: None,
        keep_scratch: false,
        cont: false,
        logfile: None,
        verbose: false,
        quiet: false,
    })?;
    Ok(output_dir)
}

fn run_rust_validate_case(case: &ValidateFixtureCase) -> Result<PathBuf> {
    let output_dir = prepare_output_dir("rust-check", case.id)?;
    validate::run(ValidateArgs {
        input: case.input_path().display().to_string(),
        reference: case.reference_path().map(|path| path.display().to_string()),
        output_json: matches!(case.kind, ValidateFixtureKind::SummaryJson)
            .then(|| output_dir.join("check.json").display().to_string()),
        errors_bed: matches!(case.kind, ValidateFixtureKind::ErrorsBed)
            .then(|| output_dir.join("errors.bed").display().to_string()),
        locations: None,
        regions_bedfile: None,
        targets_bedfile: None,
        apply_filters: false,
        limit_records: None,
        message_every: None,
        strict_homref: false,
        check_bcf_errors: false,
        all_warnings: false,
    })?;
    Ok(output_dir)
}

fn run_legacy_compare_case(case: &CompareFixtureCase, image: &str, scope: &str) -> Result<PathBuf> {
    let output_dir = prepare_output_dir(scope, case.id)?;
    let extra_flags = legacy_compare_flags(case);
    run_legacy_shell(
        image,
        &output_dir,
        &format!(
            "hap.py {truth} {query} -r {reference} -o {prefix} --force-interactive -V -X{extra} >{stdout} 2>{stderr}",
            truth = sh_quote(case.truth_path()),
            query = sh_quote(case.query_path()),
            reference = sh_quote(case.reference_path()),
            prefix = sh_quote(output_dir.join("oracle")),
            extra = extra_flags,
            stdout = sh_quote(output_dir.join("stdout.log")),
            stderr = sh_quote(output_dir.join("stderr.log")),
        ),
    )?;
    Ok(output_dir)
}

fn run_legacy_preprocess_case(
    case: &PreprocessFixtureCase,
    image: &str,
    scope: &str,
) -> Result<PathBuf> {
    let output_dir = prepare_output_dir(scope, case.id)?;
    let mut flags = String::new();
    if case.pass_only {
        flags.push_str(" --pass-only");
    }
    if let Some(path) = case.regions_bed_path() {
        flags.push_str(&format!(" -R {}", sh_quote(path)));
    }
    if let Some(locations) = case.locations {
        flags.push_str(&format!(" -l {}", sh_quote(locations)));
    }
    match case.fixchr {
        Some(true) => flags.push_str(" --fixchr"),
        Some(false) => flags.push_str(" --no-fixchr"),
        None => {}
    }
    if case.somatic {
        flags.push_str(" --somatic");
    }
    if let Some(mode) = case.set_gt {
        flags.push_str(&format!(
            " --set-gt {}",
            mode.to_possible_value().unwrap().get_name()
        ));
    }
    if case.filter_nonref {
        flags.push_str(" --filter-nonref");
    }
    if case.convert_gvcf_to_vcf {
        flags.push_str(" --convert-gvcf-to-vcf");
    }
    run_legacy_shell(
        image,
        &output_dir,
        &format!(
            "pre.py {input} {output} -r {reference}{flags} >{stdout} 2>{stderr}",
            input = sh_quote(case.input_path()),
            output = sh_quote(output_dir.join(case.output)),
            reference = sh_quote(case.reference_path()),
            flags = flags,
            stdout = sh_quote(output_dir.join("stdout.log")),
            stderr = sh_quote(output_dir.join("stderr.log")),
        ),
    )?;
    Ok(output_dir)
}

fn run_legacy_somatic_case(case: &SomaticFixtureCase, image: &str, scope: &str) -> Result<PathBuf> {
    let output_dir = prepare_output_dir(scope, case.id)?;
    let mut flags = String::new();
    if case.include_nonpass {
        flags.push_str(" -P");
    }
    if let Some(path) = case.fp_bed_path() {
        flags.push_str(&format!(" -f {}", sh_quote(path)));
    }
    if case.count_unk {
        flags.push_str(" --count-unk");
    }
    run_legacy_shell(
        image,
        &output_dir,
        &format!(
            "som.py {truth} {query} -o {prefix} -r {reference}{flags} >{stdout} 2>{stderr}",
            truth = sh_quote(case.truth_path()),
            query = sh_quote(case.query_path()),
            prefix = sh_quote(output_dir.join("oracle")),
            reference = sh_quote(case.reference_path()),
            flags = flags,
            stdout = sh_quote(output_dir.join("stdout.log")),
            stderr = sh_quote(output_dir.join("stderr.log")),
        ),
    )?;
    Ok(output_dir)
}

fn run_legacy_validate_case(
    case: &ValidateFixtureCase,
    image: &str,
    scope: &str,
) -> Result<PathBuf> {
    let output_dir = prepare_output_dir(scope, case.id)?;
    let shell = match case.kind {
        ValidateFixtureKind::SummaryJson => format!(
            "vcfcheck {input} -o {output} >{stdout} 2>{stderr}",
            input = sh_quote(case.input_path()),
            output = sh_quote(output_dir.join("check.json")),
            stdout = sh_quote(output_dir.join("stdout.log")),
            stderr = sh_quote(output_dir.join("stderr.log")),
        ),
        ValidateFixtureKind::ErrorsBed => format!(
            "true >{stdout} 2>{stderr}",
            stdout = sh_quote(output_dir.join("stdout.log")),
            stderr = sh_quote(output_dir.join("stderr.log")),
        ),
    };
    run_legacy_shell(image, &output_dir, &shell)?;
    Ok(output_dir)
}

fn run_legacy_shell(image: &str, output_dir: &Path, shell_cmd: &str) -> Result<()> {
    let repo = repo_root();
    let repo_str = repo.display().to_string();
    let output_rel = output_dir
        .strip_prefix(&repo)
        .context("failed to compute repo-relative output path")?;
    let status = Command::new("docker")
        .args([
            "run",
            "--rm",
            "-v",
            &format!("{repo_str}:{repo_str}"),
            "-w",
            &repo_str,
            image,
            "sh",
            "-lc",
            &format!(
                "rm -rf {out} && mkdir -p {out} && {cmd}",
                out = sh_quote(output_rel),
                cmd = shell_cmd,
            ),
        ])
        .status()
        .context("failed to invoke docker")?;
    if !status.success() {
        bail!("legacy docker run failed");
    }
    Ok(())
}

fn legacy_compare_flags(case: &CompareFixtureCase) -> String {
    let mut extra_flags = String::new();
    if let Some(path) = case.fp_bed_path() {
        extra_flags.push_str(&format!(" -f {}", sh_quote(path)));
    }
    if let Some(path) = case.restrict_bed_path() {
        extra_flags.push_str(&format!(" -R {}", sh_quote(path)));
    }
    if case.pass_only {
        extra_flags.push_str(" --pass-only");
    }
    extra_flags
}

fn compare_compare_outputs(case: &CompareFixtureCase, output_dir: &Path) -> Result<()> {
    let expected_dir = case.expected_dir();
    let expected_summary = fs::read_to_string(expected_dir.join("summary.csv"))?;
    compare_text_file(
        &expected_dir.join("summary.csv"),
        &output_dir.join("oracle.summary.csv"),
    )?;
    compare_text_file(
        &expected_dir.join("extended.csv"),
        &output_dir.join("oracle.extended.csv"),
    )?;
    compare_normalized_vcf(
        &expected_dir.join("vcf.normalized.vcf"),
        &output_dir.join("oracle.vcf.gz"),
    )?;
    ensure_text_contains(
        &output_dir.join("oracle.runinfo.json"),
        &["\"hap.py\"", "\"final_args\"", "\"runInfo\""],
    )?;
    ensure_gzip_text_contains(
        &output_dir.join("oracle.metrics.json.gz"),
        &[
            "\"hap.py\"",
            "\"metrics\"",
            "\"summary.metrics\"",
            "\"all.metrics\"",
            "\"roc.all\"",
        ],
    )?;
    let mut required_outputs = vec![
        output_dir.join("oracle.vcf.gz.tbi"),
        output_dir.join("oracle.roc.all.csv.gz"),
    ];
    for variant_type in ["SNP", "INDEL"] {
        if expected_summary
            .lines()
            .any(|line| line.starts_with(&format!("{variant_type},")))
        {
            required_outputs
                .push(output_dir.join(format!("oracle.roc.Locations.{variant_type}.csv.gz")));
            required_outputs
                .push(output_dir.join(format!("oracle.roc.Locations.{variant_type}.PASS.csv.gz")));
        }
    }
    for path in required_outputs {
        if !path.exists() {
            bail!("missing expected compare output {}", path.display());
        }
    }
    Ok(())
}

fn compare_preprocess_outputs(case: &PreprocessFixtureCase, output_dir: &Path) -> Result<()> {
    let expected = canonical_preprocess_records(&fs::read_to_string(
        case.expected_dir().join("processed.normalized.vcf"),
    )?);
    let actual =
        canonical_preprocess_records(&crate::vcf::read_text(&output_dir.join(case.output))?);
    if expected != actual {
        bail!("normalized preprocess VCF mismatch");
    }
    Ok(())
}

fn compare_quantify_outputs(case: &QuantifyFixtureCase, output_dir: &Path) -> Result<()> {
    compare_summary_metrics(
        &case.expected_dir().join("summary.csv"),
        &output_dir.join("oracle.summary.csv"),
    )?;
    compare_extended_metrics(
        &case.expected_dir().join("extended.csv"),
        &output_dir.join("oracle.extended.csv"),
    )
}

fn compare_somatic_outputs(case: &SomaticFixtureCase, output_dir: &Path) -> Result<()> {
    compare_somatic_metrics(
        &case.expected_dir().join("stats.csv"),
        &output_dir.join("oracle.stats.csv"),
    )?;
    ensure_text_contains(
        &output_dir.join("oracle.metrics.json"),
        &["\"som.py\"", "\"metrics\"", "\"result\""],
    )?;
    Ok(())
}

fn compare_validate_outputs(case: &ValidateFixtureCase, output_dir: &Path) -> Result<()> {
    match case.kind {
        ValidateFixtureKind::SummaryJson => compare_validate_json(
            &case.expected_dir().join("check.json"),
            &output_dir.join("check.json"),
        ),
        ValidateFixtureKind::ErrorsBed => compare_sorted_text(
            &case.expected_dir().join("errors.bed"),
            &output_dir.join("errors.bed"),
        ),
    }
}

fn ensure_text_contains(path: &Path, needles: &[&str]) -> Result<()> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    for needle in needles {
        if !text.contains(needle) {
            bail!("{} missing expected content {}", path.display(), needle);
        }
    }
    Ok(())
}

fn ensure_gzip_text_contains(path: &Path, needles: &[&str]) -> Result<()> {
    let text = crate::vcf::read_text(path)?;
    for needle in needles {
        if !text.contains(needle) {
            bail!("{} missing expected content {}", path.display(), needle);
        }
    }
    Ok(())
}

fn save_expected_compare_outputs(case: &CompareFixtureCase, output_dir: &Path) -> Result<()> {
    let expected_dir = case.expected_dir();
    fs::create_dir_all(&expected_dir)?;
    fs::copy(
        output_dir.join("oracle.summary.csv"),
        expected_dir.join("summary.csv"),
    )?;
    fs::copy(
        output_dir.join("oracle.extended.csv"),
        expected_dir.join("extended.csv"),
    )?;
    fs::write(
        expected_dir.join("vcf.normalized.vcf"),
        normalized_vcf_records(&output_dir.join("oracle.vcf.gz"))?,
    )?;
    Ok(())
}

fn save_expected_preprocess_outputs(case: &PreprocessFixtureCase, output_dir: &Path) -> Result<()> {
    let expected_dir = case.expected_dir();
    fs::create_dir_all(&expected_dir)?;
    fs::write(
        expected_dir.join("processed.normalized.vcf"),
        canonical_preprocess_records(&crate::vcf::read_text(&output_dir.join(case.output))?),
    )?;
    Ok(())
}

fn save_expected_somatic_outputs(case: &SomaticFixtureCase, output_dir: &Path) -> Result<()> {
    let expected_dir = case.expected_dir();
    fs::create_dir_all(&expected_dir)?;
    fs::copy(
        output_dir.join("oracle.stats.csv"),
        expected_dir.join("stats.csv"),
    )?;
    Ok(())
}

fn save_expected_validate_outputs(case: &ValidateFixtureCase, output_dir: &Path) -> Result<()> {
    let expected_dir = case.expected_dir();
    fs::create_dir_all(&expected_dir)?;
    match case.kind {
        ValidateFixtureKind::SummaryJson => {
            fs::copy(
                output_dir.join("check.json"),
                expected_dir.join("check.json"),
            )?;
        }
        ValidateFixtureKind::ErrorsBed => {
            fs::copy(
                output_dir.join("errors.bed"),
                expected_dir.join("errors.bed"),
            )?;
        }
    }
    Ok(())
}

fn compare_text_file(expected: &Path, actual: &Path) -> Result<()> {
    let expected_text = fs::read_to_string(expected)
        .with_context(|| format!("failed to read {}", expected.display()))?;
    let actual_text = fs::read_to_string(actual)
        .with_context(|| format!("failed to read {}", actual.display()))?;
    if expected_text != actual_text {
        bail!(
            "file mismatch: expected {} and actual {}",
            expected.display(),
            actual.display()
        );
    }
    Ok(())
}

fn compare_normalized_vcf(expected: &Path, actual_gz: &Path) -> Result<()> {
    let expected_text = fs::read_to_string(expected)?;
    let actual_text = normalized_vcf_records(actual_gz)?;
    if expected_text != actual_text {
        bail!("normalized VCF mismatch");
    }
    Ok(())
}

fn compare_sorted_text(expected: &Path, actual: &Path) -> Result<()> {
    let expected_text = sorted_non_empty_text(expected)?;
    let actual_text = sorted_non_empty_text(actual)?;
    if expected_text != actual_text {
        bail!("sorted text mismatch");
    }
    Ok(())
}

fn compare_summary_metrics(expected: &Path, actual: &Path) -> Result<()> {
    let expected = parse_csv_rows(expected)?;
    let actual = parse_csv_rows(actual)?;
    for metric in [
        "TRUTH.TOTAL",
        "QUERY.TOTAL",
        "TRUTH.TP",
        "TRUTH.FN",
        "QUERY.FP",
        "QUERY.UNK",
        "METRIC.Recall",
        "METRIC.Precision",
        "METRIC.Frac_NA",
        "METRIC.F1_Score",
    ] {
        compare_row_metric(&expected, &actual, "Type", metric, 0.001)?;
    }
    Ok(())
}

fn compare_extended_metrics(expected: &Path, actual: &Path) -> Result<()> {
    let expected = parse_csv_rows(expected)?;
    let actual = parse_csv_rows(actual)?;
    for metric in [
        "TRUTH.TOTAL",
        "TRUTH.TP",
        "TRUTH.FN",
        "QUERY.TOTAL",
        "QUERY.FP",
        "QUERY.UNK",
        "METRIC.Recall",
        "METRIC.Precision",
        "METRIC.Frac_NA",
        "METRIC.F1_Score",
    ] {
        compare_row_metric_multi(
            &expected,
            &actual,
            &[
                "Type", "Subtype", "Subset", "Filter", "Genotype", "QQ.Field", "QQ",
            ],
            metric,
            0.001,
        )?;
    }
    Ok(())
}

fn compare_somatic_metrics(expected: &Path, actual: &Path) -> Result<()> {
    let expected = parse_csv_rows(expected)?;
    let actual = parse_csv_rows(actual)?;
    for metric in [
        "total.truth",
        "total.query",
        "tp",
        "fp",
        "fn",
        "unk",
        "recall",
        "precision",
        "fp.region.size",
        "fp.rate",
    ] {
        compare_row_metric(&expected, &actual, "type", metric, 0.001)?;
    }
    Ok(())
}

fn compare_validate_json(expected: &Path, actual: &Path) -> Result<()> {
    let expected_map = parse_simple_json(expected)?;
    let actual_map = parse_simple_json(actual)?;
    for key in [
        "records",
        "nonref",
        "OVERLAP",
        "REFPADDING",
        "SYMALT",
        "UNCERTAINLENGTH",
    ] {
        let left = expected_map
            .get(key)
            .ok_or_else(|| anyhow::anyhow!("missing expected key {key}"))?;
        let right = actual_map
            .get(key)
            .ok_or_else(|| anyhow::anyhow!("missing actual key {key}"))?;
        if left != right {
            bail!("validate json mismatch for {key}: {left} != {right}");
        }
    }
    Ok(())
}

fn compare_row_metric(
    expected: &[BTreeMap<String, String>],
    actual: &[BTreeMap<String, String>],
    label_key: &str,
    metric: &str,
    tolerance: f64,
) -> Result<()> {
    let expected_map = expected
        .iter()
        .map(|row| {
            (
                row.get(label_key).cloned().unwrap_or_default(),
                row.get(metric).cloned().unwrap_or_default(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let actual_map = actual
        .iter()
        .map(|row| {
            (
                row.get(label_key).cloned().unwrap_or_default(),
                row.get(metric).cloned().unwrap_or_default(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    if expected_map.keys().collect::<Vec<_>>() != actual_map.keys().collect::<Vec<_>>() {
        bail!("row labels differ for metric {metric}");
    }
    for label in expected_map.keys() {
        compare_value(
            metric,
            label,
            &expected_map[label],
            &actual_map[label],
            tolerance,
        )?;
    }
    Ok(())
}

fn compare_row_metric_multi(
    expected: &[BTreeMap<String, String>],
    actual: &[BTreeMap<String, String>],
    label_keys: &[&str],
    metric: &str,
    tolerance: f64,
) -> Result<()> {
    let label = |row: &BTreeMap<String, String>| {
        label_keys
            .iter()
            .map(|key| row.get(*key).cloned().unwrap_or_default())
            .collect::<Vec<_>>()
            .join("_")
    };
    let expected_map = expected
        .iter()
        .map(|row| (label(row), row.get(metric).cloned().unwrap_or_default()))
        .collect::<BTreeMap<_, _>>();
    let actual_map = actual
        .iter()
        .map(|row| (label(row), row.get(metric).cloned().unwrap_or_default()))
        .collect::<BTreeMap<_, _>>();
    if expected_map.keys().collect::<Vec<_>>() != actual_map.keys().collect::<Vec<_>>() {
        bail!("row labels differ for extended metric {metric}");
    }
    for row_label in expected_map.keys() {
        compare_value(
            metric,
            row_label,
            &expected_map[row_label],
            &actual_map[row_label],
            tolerance,
        )?;
    }
    Ok(())
}

fn compare_value(
    metric: &str,
    label: &str,
    expected: &str,
    actual: &str,
    tolerance: f64,
) -> Result<()> {
    let left = expected.parse::<f64>();
    let right = actual.parse::<f64>();
    match (left, right) {
        (Ok(left), Ok(right)) => {
            if (left - right).abs() > tolerance {
                bail!("{metric} / {label}: {left} != {right}");
            }
        }
        _ => {
            if expected != actual {
                bail!("{metric} / {label}: {expected} != {actual}");
            }
        }
    }
    Ok(())
}

fn parse_csv_rows(path: &Path) -> Result<Vec<BTreeMap<String, String>>> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("missing CSV header in {}", path.display()))?
        .split(',')
        .map(str::to_string)
        .collect::<Vec<_>>();
    Ok(lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            header
                .iter()
                .cloned()
                .zip(line.split(',').map(str::to_string))
                .collect::<BTreeMap<_, _>>()
        })
        .collect())
}

fn parse_simple_json(path: &Path) -> Result<BTreeMap<String, String>> {
    let text = fs::read_to_string(path)?;
    Ok(text
        .trim()
        .trim_start_matches('{')
        .trim_end_matches('}')
        .split(',')
        .filter_map(|entry| entry.split_once(':'))
        .map(|(key, value)| {
            (
                key.trim().trim_matches('"').to_string(),
                value.trim().trim_matches('"').to_string(),
            )
        })
        .collect())
}

fn normalized_vcf_records(path: &Path) -> Result<String> {
    Ok(normalized_lines(&crate::vcf::read_text(path)?))
}

fn normalized_lines(text: &str) -> String {
    let mut lines: Vec<&str> = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter(|line| !line.is_empty())
        .collect();
    lines.sort_unstable();
    if lines.is_empty() {
        String::new()
    } else {
        format!("{}\n", lines.join("\n"))
    }
}

fn canonical_preprocess_records(text: &str) -> String {
    let mut lines = text
        .lines()
        .filter(|line| !line.starts_with('#'))
        .filter(|line| !line.trim().is_empty())
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 10 {
                return None;
            }
            let gt = fields[9].split(':').next().unwrap_or("./.");
            Some(format!(
                "{}\t{}\t{}\t{}\t{}",
                fields[0], fields[1], fields[3], fields[4], gt
            ))
        })
        .collect::<Vec<_>>();
    lines.sort();
    if lines.is_empty() {
        String::new()
    } else {
        format!("{}\n", lines.join("\n"))
    }
}

fn sorted_non_empty_text(path: &Path) -> Result<String> {
    let mut lines: Vec<_> = fs::read_to_string(path)?
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(str::to_string)
        .collect();
    lines.sort();
    Ok(if lines.is_empty() {
        String::new()
    } else {
        format!("{}\n", lines.join("\n"))
    })
}

fn verification_root() -> PathBuf {
    repo_root().join("target/verification")
}

fn prepare_output_dir(scope: &str, case_id: &str) -> Result<PathBuf> {
    let output_dir = verification_root().join(scope).join(case_id);
    if output_dir.exists() {
        fs::remove_dir_all(&output_dir)?;
    }
    fs::create_dir_all(&output_dir)?;
    Ok(output_dir)
}

fn sh_quote(path: impl AsRef<Path>) -> String {
    let text = path.as_ref().display().to_string();
    format!("'{}'", text.replace('\'', "'\"'\"'"))
}

#[derive(Clone)]
struct RealWorldCase {
    id: &'static str,
    command: &'static str,
    legacy_tool: &'static str,
    inputs: &'static [&'static str],
    args: &'static [&'static str],
}

fn realworld_cases() -> Vec<RealWorldCase> {
    vec![
        RealWorldCase {
            id: "happy-module-germline",
            command: "germline",
            legacy_tool: "hap.py",
            inputs: &[
                "example/happy/PG_NA12878_chr21.vcf.gz",
                "example/happy/PG_NA12878_chr21.vcf.gz.tbi",
                "example/happy/NA12878_chr21.vcf.gz",
                "example/happy/NA12878_chr21.vcf.gz.tbi",
                "example/happy/PG_Conf_chr21.bed.gz",
                "example/happy/PG_Conf_chr21.bed.gz.tbi",
                "example/chr21.fa",
                "example/chr21.fa.fai",
            ],
            args: &[
                "PG_NA12878_chr21.vcf.gz",
                "NA12878_chr21.vcf.gz",
                "--reference",
                "chr21.fa",
                "--threads",
                "1",
                "--false-positives",
                "PG_Conf_chr21.bed.gz",
                "-o",
                "result",
            ],
        },
        RealWorldCase {
            id: "sompy-module-somatic",
            command: "somatic",
            legacy_tool: "som.py",
            inputs: &[
                "example/sompy/PG_admix_truth_snvs.vcf.gz",
                "example/sompy/strelka_admix_snvs.vcf.gz",
                "example/sompy/FP_admix.bed.gz",
                "example/chr21.fa",
                "example/chr21.fa.fai",
            ],
            args: &[
                "PG_admix_truth_snvs.vcf.gz",
                "strelka_admix_snvs.vcf.gz",
                "-o",
                "result",
                "--reference",
                "chr21.fa",
                "--false-positives",
                "FP_admix.bed.gz",
                "-P",
                "--count-unk",
                "--feature-table",
                "hcc.strelka.indel",
            ],
        },
        RealWorldCase {
            id: "prepy-module-pre",
            command: "pre",
            legacy_tool: "pre.py",
            inputs: &[
                "example/happy/NA12878_chr21.vcf.gz",
                "example/happy/NA12878_chr21.vcf.gz.tbi",
                "example/happy/PG_Conf_chr21.bed.gz",
                "example/happy/PG_Conf_chr21.bed.gz.tbi",
                "example/chr21.fa",
                "example/chr21.fa.fai",
            ],
            args: &[
                "--reference",
                "chr21.fa",
                "--threads",
                "1",
                "-R",
                "PG_Conf_chr21.bed.gz",
                "NA12878_chr21.vcf.gz",
                "result.vcf.gz",
            ],
        },
        RealWorldCase {
            id: "ftxpy-module-ftx",
            command: "ftx",
            legacy_tool: "ftx.py",
            inputs: &[
                "example/sompy/strelka_admix_snvs.vcf.gz",
                "example/chr21.fa",
                "example/chr21.fa.fai",
            ],
            args: &[
                "-o",
                "result",
                "--reference",
                "chr21.fa",
                "--feature-table",
                "generic",
                "strelka_admix_snvs.vcf.gz",
            ],
        },
    ]
}

fn materialize_realworld_inputs(case: &RealWorldCase) -> Result<PathBuf> {
    let root = verification_root().join("realworld").join(case.id);
    if root.exists() {
        fs::remove_dir_all(&root)?;
    }
    let legacy_work = root.join("legacy_work");
    let rust_work = root.join("rust_work");
    fs::create_dir_all(&legacy_work)?;
    fs::create_dir_all(&rust_work)?;
    for rel in case.inputs {
        let name = Path::new(rel)
            .file_name()
            .ok_or_else(|| anyhow::anyhow!("invalid input path {rel}"))?;
        let output = Command::new("git")
            .args(["show", &format!("HEAD:{rel}")])
            .output()
            .with_context(|| format!("failed to restore {rel} from git"))?;
        if !output.status.success() {
            bail!("git show failed for {}", rel);
        }
        fs::write(legacy_work.join(name), &output.stdout)?;
        fs::write(rust_work.join(name), &output.stdout)?;
    }
    Ok(root)
}

fn run_legacy_realworld_case(case: &RealWorldCase, image: &str, workdir: &Path) -> Result<()> {
    let input_dir = workdir.join("legacy_work");
    let mut cmd = Command::new("docker");
    cmd.arg("run")
        .arg("--rm")
        .arg("-v")
        .arg(format!("{}:{}", input_dir.display(), input_dir.display()))
        .arg("-w")
        .arg(&input_dir)
        .arg(image)
        .arg(case.legacy_tool);
    for arg in case.args {
        cmd.arg(*arg);
    }
    let status = cmd
        .status()
        .context("failed to run legacy real-world case")?;
    if !status.success() {
        bail!("legacy real-world case {} failed", case.id);
    }
    Ok(())
}

fn run_rust_realworld_case(case: &RealWorldCase, workdir: &Path) -> Result<()> {
    let hap_bin = ensure_hap_binary()?;
    let input_dir = workdir.join("rust_work");
    let mut cmd = Command::new(hap_bin);
    cmd.current_dir(&input_dir).arg(case.command);
    if let Some(timestamp) = legacy_result_timestamp(&workdir.join("legacy_work"))? {
        cmd.env("HAP_FIXED_TIMESTAMP", timestamp);
    }
    for arg in case.args {
        cmd.arg(*arg);
    }
    let status = cmd.status().context("failed to run rust real-world case")?;
    if !status.success() {
        bail!("rust real-world case {} failed", case.id);
    }
    Ok(())
}

fn ensure_hap_binary() -> Result<PathBuf> {
    let path = repo_root().join("target/debug/hap");
    if path.exists() {
        return Ok(path);
    }
    let status = Command::new("cargo")
        .args(["build", "--bin", "hap"])
        .current_dir(repo_root())
        .status()
        .context("failed to build hap binary")?;
    if !status.success() {
        bail!("cargo build --bin hap failed");
    }
    Ok(path)
}

fn legacy_result_timestamp(legacy_dir: &Path) -> Result<Option<String>> {
    let json_path = legacy_dir.join("result.metrics.json");
    if json_path.exists() {
        return Ok(extract_json_string(
            &fs::read_to_string(json_path)?,
            "timestamp",
        ));
    }
    let gzip_path = legacy_dir.join("result.metrics.json.gz");
    if gzip_path.exists() {
        return Ok(extract_json_string(
            &crate::vcf::read_text(&gzip_path)?,
            "timestamp",
        ));
    }
    Ok(None)
}

fn extract_json_string(text: &str, key: &str) -> Option<String> {
    let marker = format!("\"{key}\":");
    let start = text.find(&marker)? + marker.len();
    let tail = text[start..].trim_start();
    if !tail.starts_with('"') {
        return None;
    }
    let tail = &tail[1..];
    let end = tail.find('"')?;
    Some(tail[..end].to_string())
}

fn compare_output_trees(expected_dir: &Path, actual_dir: &Path) -> Result<()> {
    let mut expected = collect_tree(expected_dir)?;
    let mut actual = collect_tree(actual_dir)?;
    expected.sort();
    actual.sort();
    if expected != actual {
        bail!("output file sets differ: {:?} vs {:?}", expected, actual);
    }
    for rel in expected {
        let left = fs::read(expected_dir.join(&rel))?;
        let right = fs::read(actual_dir.join(&rel))?;
        if left != right {
            bail!("file differs: {}", rel.display());
        }
    }
    Ok(())
}

fn collect_tree(root: &Path) -> Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    for entry in fs::read_dir(root)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            let nested = collect_tree(&path)?;
            for rel in nested {
                let sub = path
                    .file_name()
                    .map(PathBuf::from)
                    .unwrap_or_default()
                    .join(rel);
                files.push(sub);
            }
        } else if path
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| name.starts_with("result"))
            .unwrap_or(false)
        {
            files.push(
                path.strip_prefix(root)
                    .map(PathBuf::from)
                    .unwrap_or(path.clone()),
            );
        }
    }
    Ok(files)
}
