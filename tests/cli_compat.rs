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
    for command in ["germline", "compare"] {
        for version_flag in ["-v", "--version"] {
            let output = hap()
                .args([command, version_flag])
                .output()
                .expect("run germline version");
            assert_eq!(output.status.code(), Some(0), "{command} {version_flag}");
            assert_eq!(output.stdout, b"Hap.py \n", "{command} {version_flag}");
            assert!(output.stderr.is_empty(), "{command} {version_flag}");
        }
    }
}

#[test]
fn pre_and_qfy_versions_require_normal_arguments() {
    for command in ["pre", "preprocess", "prepy", "quantify", "qfy"] {
        for version_flag in ["-v", "--version"] {
            let output = hap()
                .args([command, version_flag])
                .output()
                .expect("run standalone script version");
            assert_eq!(output.status.code(), Some(2), "{command} {version_flag}");
            assert!(output.stdout.is_empty(), "{command} {version_flag}");
        }
    }
}

#[test]
fn pre_and_qfy_versions_match_pinned_script_output_after_validation() {
    for command in ["pre", "preprocess", "prepy"] {
        for version_flag in ["-v", "--version"] {
            let output = hap()
                .args([command, "input.vcf", "output.vcf", version_flag])
                .output()
                .expect("run pre version with required arguments");
            assert_eq!(output.status.code(), Some(0), "{command} {version_flag}");
            assert_eq!(output.stdout, b"pre.py \n", "{command} {version_flag}");
            assert!(output.stderr.is_empty(), "{command} {version_flag}");
        }
    }

    for command in ["quantify", "qfy"] {
        for version_flag in ["-v", "--version"] {
            let output = hap()
                .args([
                    command,
                    "input.vcf",
                    "-o",
                    "report",
                    "-r",
                    "ref.fa",
                    version_flag,
                ])
                .output()
                .expect("run qfy version with required arguments");
            assert_eq!(output.status.code(), Some(0), "{command} {version_flag}");
            assert_eq!(output.stdout, b"qfy.py \n", "{command} {version_flag}");
            assert!(output.stderr.is_empty(), "{command} {version_flag}");
        }
    }
}

#[test]
fn version_shim_does_not_change_other_subcommand_validation() {
    for command in ["ftx", "somatic"] {
        for version_flag in ["-v", "--version"] {
            let output = hap()
                .args([command, version_flag])
                .output()
                .expect("run unsupported subcommand version");
            assert_eq!(output.status.code(), Some(2), "{command} {version_flag}");
            assert!(output.stdout.is_empty(), "{command} {version_flag}");
        }
    }

    let output = hap()
        .args([
            "germline",
            "-r",
            "ref.fa",
            "-o",
            "report",
            "--",
            "--version",
            "query.vcf",
        ])
        .output()
        .expect("run germline with a version-like positional input");
    assert_eq!(
        output.status.code(),
        Some(1),
        "a positional --version should reach runtime validation"
    );
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).expect("error is UTF-8");
    assert!(
        !stderr.is_empty(),
        "runtime validation should explain the failure"
    );
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

#[test]
fn legacy_unknown_arguments_keep_wrapper_specific_exit_codes() {
    let germline = hap()
        .args([
            "germline",
            "truth.vcf",
            "query.vcf",
            "-o",
            "result",
            "--definitely-unknown",
        ])
        .output()
        .expect("run germline unknown option");
    assert_eq!(germline.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&germline.stdout).contains("Usage:"));

    let version = hap()
        .args(["germline", "--version", "--definitely-unknown"])
        .output()
        .expect("run germline version with unknown option");
    assert_eq!(version.status.code(), Some(1));
    assert!(!version.stdout.starts_with(b"Hap.py "));

    for command in [
        vec!["pre", "input.vcf", "output.vcf", "--definitely-unknown"],
        vec![
            "qfy",
            "input.vcf",
            "-o",
            "result",
            "-r",
            "ref.fa",
            "--definitely-unknown",
        ],
    ] {
        let output = hap()
            .args(command)
            .output()
            .expect("run legacy wrapper unknown option");
        assert_eq!(output.status.code(), Some(0));
        assert!(String::from_utf8_lossy(&output.stdout).contains("Usage:"));
    }

    let ftx = hap()
        .args(["ftx", "input.vcf", "-o", "result", "--definitely-unknown"])
        .output()
        .expect("run ftx unknown option");
    assert_eq!(ftx.status.code(), Some(2));
}

#[test]
fn germline_missing_required_arguments_use_legacy_runtime_exit_code() {
    for arguments in [
        vec!["germline", "truth.vcf", "query.vcf"],
        vec!["compare", "truth.vcf", "-o", "result"],
    ] {
        let output = hap()
            .args(arguments)
            .output()
            .expect("run incomplete germline command");
        assert_eq!(output.status.code(), Some(1));
    }
}

#[test]
fn legacy_store_true_switches_reject_following_boolean_tokens() {
    for (arguments, expected_code) in [
        (
            vec![
                "germline",
                "truth.vcf",
                "query.vcf",
                "-o",
                "result",
                "--fixchr",
                "false",
            ],
            1,
        ),
        (
            vec!["pre", "input.vcf", "output.vcf", "--fixchr", "false"],
            0,
        ),
        (
            vec![
                "somatic",
                "truth.vcf",
                "query.vcf",
                "-o",
                "result",
                "--fixchr-truth",
                "false",
            ],
            2,
        ),
    ] {
        let output = hap()
            .args(arguments)
            .output()
            .expect("run store_true switch followed by a Boolean token");
        assert_eq!(output.status.code(), Some(expected_code));
    }
}
