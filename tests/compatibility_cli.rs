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
