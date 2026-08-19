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
