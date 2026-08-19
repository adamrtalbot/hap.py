//! Runtime deprecation warnings for options hap-rs accepts and ignores.

use super::cli::Command;

pub(crate) const VCFEVAL_RUNTIME_DEPRECATION: &str = "warning: --engine-vcfeval-path and --engine-vcfeval-template are deprecated and ignored; use --engine vcfeval --reference <FASTA>; these options will be removed in hap-rs 1.0.0";

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

    #[test]
    fn vcfeval_deprecation_warning_has_exact_replacement_and_deadline() {
        assert_eq!(
            VCFEVAL_RUNTIME_DEPRECATION,
            "warning: --engine-vcfeval-path and --engine-vcfeval-template are deprecated and ignored; use --engine vcfeval --reference <FASTA>; these options will be removed in hap-rs 1.0.0"
        );
    }
}
