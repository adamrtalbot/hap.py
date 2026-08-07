//! Explicit policies for legacy command-line behavior.

use super::cli::Command;

pub(crate) const LEGACY_SUCCESS_EXIT_DEPRECATION: &str = "warning: success exit status for unknown pre/quantify arguments is deprecated; unknown arguments will exit non-zero in hap-rs 1.0.0";
pub(crate) const VCFEVAL_RUNTIME_DEPRECATION: &str = "warning: --engine-vcfeval-path and --engine-vcfeval-template are deprecated and ignored; use --engine vcfeval --reference <FASTA>; these options will be removed in hap-rs 1.0.0";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum UsageErrorPolicy {
    LegacyFailure,
    DeprecatedSuccess,
}

impl UsageErrorPolicy {
    pub(crate) fn exit_code(self) -> i32 {
        match self {
            Self::LegacyFailure => 1,
            Self::DeprecatedSuccess => 0,
        }
    }

    pub(crate) fn warning(self) -> Option<&'static str> {
        matches!(self, Self::DeprecatedSuccess).then_some(LEGACY_SUCCESS_EXIT_DEPRECATION)
    }
}

pub(crate) fn usage_error_policy(arguments: &[std::ffi::OsString]) -> Option<UsageErrorPolicy> {
    match arguments.get(1).and_then(|value| value.to_str())? {
        "germline" | "compare" => Some(UsageErrorPolicy::LegacyFailure),
        "pre" | "preprocess" | "prepy" | "quantify" | "qfy" => {
            Some(UsageErrorPolicy::DeprecatedSuccess)
        }
        _ => None,
    }
}

pub(crate) fn emit_deprecation_warnings(command: &Command) {
    if let Command::Germline(args) = command
        && (args.engine_vcfeval.is_some() || args.engine_vcfeval_template.is_some())
    {
        eprintln!("{VCFEVAL_RUNTIME_DEPRECATION}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(values: &[&str]) -> Vec<std::ffi::OsString> {
        values.iter().map(std::ffi::OsString::from).collect()
    }

    #[test]
    fn legacy_only_pre_and_quantify_unknown_options_retain_deprecated_success() {
        for command in ["pre", "preprocess", "prepy", "quantify", "qfy"] {
            let policy = usage_error_policy(&args(&["hap", command])).unwrap();
            assert_eq!(policy, UsageErrorPolicy::DeprecatedSuccess);
            assert_eq!(policy.exit_code(), 0);
            assert_eq!(policy.warning(), Some(LEGACY_SUCCESS_EXIT_DEPRECATION));
        }
    }

    #[test]
    fn normative_commands_do_not_override_clap_usage_errors() {
        for command in ["somatic", "ftx", "validate", "vcfcheck"] {
            assert_eq!(usage_error_policy(&args(&["hap", command])), None);
        }
    }

    #[test]
    fn vcfeval_deprecation_warning_has_exact_replacement_and_deadline() {
        assert_eq!(
            VCFEVAL_RUNTIME_DEPRECATION,
            "warning: --engine-vcfeval-path and --engine-vcfeval-template are deprecated and ignored; use --engine vcfeval --reference <FASTA>; these options will be removed in hap-rs 1.0.0"
        );
    }
}
