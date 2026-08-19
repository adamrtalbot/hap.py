//! Runtime warnings for options hap-rs accepts and ignores.

use super::cli::Command;

pub(crate) const VCFEVAL_RUNTIME_OPTIONS_IGNORED: &str = "warning: --engine-vcfeval-path and --engine-vcfeval-template are ignored; hap-rs reads the FASTA supplied with --reference instead";

/// Returns the warning an invocation earns, so the decision is testable without
/// capturing standard error.
pub(crate) fn vcfeval_options_ignored_warning(command: &Command) -> Option<&'static str> {
    match command {
        Command::Germline(args)
            if args.engine_vcfeval.is_some() || args.engine_vcfeval_template.is_some() =>
        {
            Some(VCFEVAL_RUNTIME_OPTIONS_IGNORED)
        }
        _ => None,
    }
}

pub(crate) fn emit_compatibility_warnings(command: &Command) {
    if let Some(warning) = vcfeval_options_ignored_warning(command) {
        eprintln!("{warning}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli_compat::cli::Cli;
    use clap::Parser;

    fn germline_with(path: Option<&str>, template: Option<&str>) -> Command {
        let mut argv = vec![
            "hap",
            "germline",
            "truth.vcf",
            "query.vcf",
            "-o",
            "report",
            "-r",
            "ref.fa",
        ];
        if let Some(path) = path {
            argv.extend(["--engine-vcfeval-path", path]);
        }
        if let Some(template) = template {
            argv.extend(["--engine-vcfeval-template", template]);
        }
        Cli::try_parse_from(argv)
            .expect("germline invocation should parse")
            .command
    }

    #[test]
    fn vcfeval_ignored_options_warning_names_the_reference_replacement() {
        assert_eq!(
            VCFEVAL_RUNTIME_OPTIONS_IGNORED,
            "warning: --engine-vcfeval-path and --engine-vcfeval-template are ignored; hap-rs reads the FASTA supplied with --reference instead"
        );
    }

    #[test]
    fn germline_warns_once_either_vcfeval_option_is_supplied() {
        for command in [
            germline_with(Some("rtg"), None),
            germline_with(None, Some("t.sdf")),
        ] {
            assert_eq!(
                vcfeval_options_ignored_warning(&command),
                Some(VCFEVAL_RUNTIME_OPTIONS_IGNORED)
            );
        }
    }

    #[test]
    fn germline_stays_silent_when_neither_vcfeval_option_is_supplied() {
        assert_eq!(
            vcfeval_options_ignored_warning(&germline_with(None, None)),
            None
        );
    }
}
