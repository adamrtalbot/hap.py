use std::process::{Command, Output};

fn hap(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_hap"))
        .args(arguments)
        .output()
        .expect("hap process starts")
}

#[test]
fn legacy_aliases_keep_command_specific_unknown_argument_exit_codes() {
    for (command, arguments, expected_code) in [
        (
            "compare",
            vec!["compare", "truth.vcf", "query.vcf", "-o", "report", "--wat"],
            1,
        ),
        (
            "prepy",
            vec!["prepy", "input.vcf", "output.vcf.gz", "--wat"],
            0,
        ),
        (
            "qfy",
            vec![
                "qfy",
                "input.vcf.gz",
                "-o",
                "report",
                "-r",
                "reference.fa",
                "--wat",
            ],
            0,
        ),
    ] {
        let output = hap(&arguments);
        assert_eq!(
            output.status.code(),
            Some(expected_code),
            "legacy exit status for {command}"
        );
        assert!(
            String::from_utf8_lossy(&output.stderr).contains("unexpected argument '--wat'"),
            "legacy diagnostic for {command}"
        );
        assert!(
            String::from_utf8_lossy(&output.stdout).contains("Usage:"),
            "legacy help for {command}"
        );
    }
}

#[test]
fn ordinary_adapter_errors_keep_clap_exit_status() {
    let output = hap(&["validate", "input.vcf", "--wat"]);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("unexpected argument '--wat'"));
}

#[test]
fn legacy_version_spelling_short_circuits_before_file_loading() {
    let germline = hap(&[
        "compare",
        "missing-truth.vcf",
        "missing-query.vcf",
        "-o",
        "report",
        "--version",
    ]);
    assert!(germline.status.success());
    assert_eq!(String::from_utf8_lossy(&germline.stdout), "Hap.py \n");

    let quantify = hap(&[
        "qfy",
        "missing.vcf.gz",
        "-o",
        "report",
        "-r",
        "missing.fa",
        "--version",
    ]);
    assert!(quantify.status.success());
    assert_eq!(String::from_utf8_lossy(&quantify.stdout), "qfy.py \n");
}

#[test]
fn semantic_validation_precedes_input_loading_and_output_creation() {
    let directory = tempfile::tempdir().expect("temporary output directory");
    let output = directory.path().join("must-not-exist.vcf");
    let process = Command::new(env!("CARGO_BIN_EXE_hap"))
        .args([
            "pre",
            "missing-input.vcf",
            output.to_str().unwrap(),
            "--reference",
            "missing-reference.fa",
        ])
        .output()
        .expect("hap process starts");

    assert_eq!(process.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&process.stderr).contains("plain VCF output cannot be indexed")
    );
    assert!(!output.exists());
}
