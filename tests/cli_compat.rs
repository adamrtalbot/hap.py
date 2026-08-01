use std::process::Command;

fn hap() -> Command {
    Command::new(env!("CARGO_BIN_EXE_hap"))
}

#[test]
fn public_and_legacy_help_surfaces_exit_successfully() {
    let output = hap().arg("--help").output().expect("run hap --help");
    assert_eq!(output.status.code(), Some(0));
    let stdout = String::from_utf8(output.stdout).expect("help is UTF-8");
    assert!(stdout.contains("quantify"));
    assert!(stdout.contains("qfy"));
    assert!(stdout.contains("validate"));
    assert!(stdout.contains("vcfcheck"));

    for alias in ["qfy", "vcfcheck"] {
        let output = hap()
            .args([alias, "--help"])
            .output()
            .expect("run legacy alias help");
        assert_eq!(output.status.code(), Some(0), "{alias}");
    }
}

#[test]
fn germline_version_exits_before_required_input_validation() {
    for version_flag in ["-v", "--version"] {
        let output = hap()
            .args(["germline", version_flag])
            .output()
            .expect("run germline version");
        assert_eq!(output.status.code(), Some(0), "{version_flag}");
        assert_eq!(
            String::from_utf8(output.stdout)
                .expect("version is UTF-8")
                .trim(),
            env!("CARGO_PKG_VERSION")
        );
    }
}

#[test]
fn external_engine_options_are_public_but_engine_values_remain_validated() {
    let output = hap()
        .args([
            "germline",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-r",
            "ref.fa",
            "-o",
            "report",
            "--engine",
            "unsupported",
        ])
        .output()
        .expect("run invalid engine value");
    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).expect("error is UTF-8");
    assert!(stderr.contains("invalid value"), "{stderr}");

    for engine_args in [
        vec![
            "--engine",
            "vcfeval",
            "--engine-vcfeval-path",
            "rtg-custom",
            "--engine-vcfeval-template",
            "reference.sdf",
        ],
        vec!["--engine", "scmp-distance", "--scmp-distance", "42"],
        vec!["--engine", "scmp-distance", "--lose-match-distance", "7"],
    ] {
        let mut command = hap();
        command.args([
            "germline",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-r",
            "ref.fa",
            "-o",
            "report",
        ]);
        let output = command
            .args(engine_args)
            .output()
            .expect("run supported engine options");
        assert_eq!(
            output.status.code(),
            Some(1),
            "accepted legacy spellings should reach runtime validation"
        );
        let stderr = String::from_utf8(output.stderr).expect("error is UTF-8");
        assert!(!stderr.contains("unexpected argument"), "{stderr}");
        assert!(
            !stderr.is_empty(),
            "runtime validation should explain the failure"
        );
    }

    let output = hap()
        .args([
            "germline",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-r",
            "ref.fa",
            "-o",
            "report",
            "--engine",
            "scmp-distance",
            "--scmp-distance",
            "not-a-number",
        ])
        .output()
        .expect("run invalid distance value");
    assert_eq!(output.status.code(), Some(2));
    let stderr = String::from_utf8(output.stderr).expect("error is UTF-8");
    assert!(stderr.contains("invalid value"), "{stderr}");
}

#[test]
fn legacy_aliases_keep_parser_and_runtime_exit_codes_distinct() {
    let output = hap()
        .args(["qfy", "annotated.vcf.gz", "-r", "ref.fa"])
        .output()
        .expect("run qfy without required output prefix");
    assert_eq!(output.status.code(), Some(2));

    let output = hap()
        .args([
            "vcfcheck",
            "does-not-exist.vcf.gz",
            "--output-file",
            "counts.json",
            "--location",
            "chr1:10-20",
        ])
        .output()
        .expect("run vcfcheck legacy spellings");
    assert_eq!(
        output.status.code(),
        Some(1),
        "accepted legacy spellings should reach runtime validation"
    );
    let stderr = String::from_utf8(output.stderr).expect("error is UTF-8");
    assert!(stderr.contains("does-not-exist.vcf.gz"), "{stderr}");
}
