use std::process::{Command, Output};

fn hap(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_hap"))
        .args(arguments)
        .output()
        .expect("hap process should start")
}

const SUBCOMMANDS: [&str; 6] = ["germline", "somatic", "pre", "quantify", "ftx", "validate"];
const ALIASES: [&str; 6] = ["compare", "preprocess", "prepy", "qfy", "ftxpy", "vcfcheck"];

fn assert_usage_error(invocation: &[&str]) {
    let label = invocation.join(" ");
    let output = hap(invocation);
    assert_eq!(output.status.code(), Some(2), "{label}");
    assert!(output.stdout.is_empty(), "{label}");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.starts_with("error: "), "{label}: {stderr}");
    assert!(stderr.contains("Usage:"), "{label}: {stderr}");
    assert!(
        stderr.contains("For more information, try '--help'"),
        "{label}: {stderr}"
    );
}

#[test]
fn normative_unknown_option_uses_clap_usage_exit_for_every_subcommand() {
    for command in SUBCOMMANDS {
        assert_usage_error(&[command, "--definitely-invalid"]);
    }
}

#[test]
fn normative_missing_required_argument_uses_clap_usage_exit_for_every_subcommand() {
    for command in SUBCOMMANDS {
        assert_usage_error(&[command]);
    }
}

#[test]
fn normative_legacy_aliases_use_clap_usage_exit() {
    for alias in ALIASES {
        assert_usage_error(&[alias, "--definitely-invalid"]);
        assert_usage_error(&[alias]);
    }
}

const VCFEVAL_OPTIONS_IGNORED_WARNING: &str = "warning: --engine-vcfeval-path and --engine-vcfeval-template are ignored; hap-rs reads the FASTA supplied with --reference instead";

/// The warning is emitted after parsing and before dispatch, so an invocation
/// that fails afterwards still carries it. Absent paths keep the test cheap.
fn germline_stderr(extra: &[&str]) -> String {
    let mut invocation = vec![
        "germline",
        "absent-truth.vcf",
        "absent-query.vcf",
        "-o",
        "absent-report",
        "-r",
        "absent-ref.fa",
    ];
    invocation.extend_from_slice(extra);
    String::from_utf8_lossy(&hap(&invocation).stderr).into_owned()
}

#[test]
fn germline_warns_on_stderr_for_each_ignored_vcfeval_option() {
    for option in ["--engine-vcfeval-path", "--engine-vcfeval-template"] {
        let stderr = germline_stderr(&[option, "ignored-value"]);
        assert!(
            stderr.contains(VCFEVAL_OPTIONS_IGNORED_WARNING),
            "{option}: {stderr}"
        );
    }
}

#[test]
fn germline_without_the_vcfeval_options_emits_no_warning() {
    let stderr = germline_stderr(&[]);
    assert!(!stderr.contains("--engine-vcfeval-path"), "{stderr}");
}

/// The reference comes from the command line only. `HGREF` and `HG19` name a
/// file that exists, so a run that still fails proves neither variable is read.
fn hap_with_reference_environment(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_hap"))
        .args(arguments)
        .env("HGREF", env!("CARGO_BIN_EXE_hap"))
        .env("HG19", env!("CARGO_BIN_EXE_hap"))
        .output()
        .expect("hap process should start")
}

#[test]
fn normative_reference_environment_variables_are_never_consulted() {
    let invocations: [&[&str]; 4] = [
        &[
            "germline",
            "absent-truth.vcf",
            "absent-query.vcf",
            "-o",
            "absent-report",
        ],
        &["pre", "absent-input.vcf", "absent-output.vcf.gz"],
        &[
            "ftx",
            "absent-input.vcf",
            "-o",
            "absent-output.csv",
            "--normalize",
        ],
        &[
            "somatic",
            "absent-truth.vcf",
            "absent-query.vcf",
            "-o",
            "absent-report",
            "--normalize-all",
        ],
    ];
    for invocation in invocations {
        let label = invocation.join(" ");
        let output = hap_with_reference_environment(invocation);
        assert!(!output.status.success(), "{label}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(stderr.contains("--reference"), "{label}: {stderr}");
        assert!(!stderr.contains("HGREF"), "{label}: {stderr}");
        assert!(!stderr.contains("HG19"), "{label}: {stderr}");
    }
}

/// Legacy `som.py` compares alleles with no reference at all, so `--reference`
/// is demanded only by the two somatic features that open one. The parser must
/// therefore let the argument stay absent.
#[test]
fn normative_somatic_accepts_an_absent_reference() {
    let output = hap_with_reference_environment(&[
        "somatic",
        "absent-truth.vcf",
        "absent-query.vcf",
        "-o",
        "absent-report",
    ]);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(!stderr.contains("--reference"), "{stderr}");
    assert!(!stderr.contains("Usage:"), "{stderr}");
}
