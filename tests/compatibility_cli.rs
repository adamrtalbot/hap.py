use std::process::{Command, Output};

const LEGACY_SUCCESS_WARNING: &str =
    "success exit status for unknown pre/quantify arguments is deprecated";

fn hap(arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_hap"))
        .args(arguments)
        .output()
        .expect("hap process should start")
}

#[test]
fn legacy_only_pre_and_quantify_unknown_options_exit_success_with_transition_warning() {
    for command in ["pre", "quantify"] {
        let output = hap(&[command, "--definitely-invalid"]);
        assert_eq!(output.status.code(), Some(0), "{command}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert!(
            stderr.contains("unexpected argument"),
            "{command}: {stderr}"
        );
        assert!(
            stderr.contains(LEGACY_SUCCESS_WARNING),
            "{command}: {stderr}"
        );
        assert!(stderr.contains("hap-rs 1.0.0"), "{command}: {stderr}");
    }
}

#[test]
fn normative_pre_and_quantify_missing_arguments_use_standard_nonzero_exit() {
    for command in ["pre", "quantify"] {
        let output = hap(&[command]);
        assert_eq!(output.status.code(), Some(2), "{command}");
        assert!(!String::from_utf8_lossy(&output.stderr).contains(LEGACY_SUCCESS_WARNING));
    }
}

#[test]
fn normative_validate_invalid_option_uses_standard_nonzero_exit() {
    let output = hap(&["validate", "--definitely-invalid"]);
    assert_eq!(output.status.code(), Some(2));
    assert!(!String::from_utf8_lossy(&output.stderr).contains(LEGACY_SUCCESS_WARNING));
}

#[test]
fn legacy_only_germline_invalid_option_retains_failure_exit_one() {
    let output = hap(&["germline", "--definitely-invalid"]);
    assert_eq!(output.status.code(), Some(1));
    assert!(!String::from_utf8_lossy(&output.stderr).contains(LEGACY_SUCCESS_WARNING));
}
