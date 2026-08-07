//! Governed compatibility behavior at command and input adapter boundaries.
//!
//! Core comparison and reporting code should not choose legacy behavior. This
//! module makes those choices explicit so every emulation has a policy, a
//! warning where applicable, and regression evidence.

use crate::{cli::Command, vcf::LocationFilter};

pub(crate) const LEGACY_SUCCESS_EXIT_DEPRECATION: &str = "warning: success exit status for unknown pre/quantify arguments is deprecated; unknown arguments will exit non-zero in hap-rs 1.0.0";

pub(crate) const VCFEVAL_RUNTIME_DEPRECATION: &str = "warning: --engine-vcfeval-path and --engine-vcfeval-template are deprecated and ignored; use --engine vcfeval --reference <FASTA>; these options will be removed in hap-rs 1.0.0";

/// Whether comma-separated locations are independent streams or a set union.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum LocationStreamPolicy {
    /// pre.py runs block splitting once for every requested location.
    IndependentLegacyStreams,
    /// Normative region selection emits each selected record at most once.
    #[allow(
        dead_code,
        reason = "normative alternative is retained as policy and regression evidence"
    )]
    SetUnion,
}

/// Span construction used when adapting VCF alleles to SCMP `RefVar`s.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ScmpRefVarSpanPolicy {
    /// Pinned SCMP incorrectly derives the reference span from ALT length.
    LegacyAltLength,
    /// A VCF reference span is derived from REF length.
    #[allow(
        dead_code,
        reason = "normative alternative is retained as policy and regression evidence"
    )]
    ReferenceLength,
}

/// Lifecycle of bcftools-compatible per-allele count storage.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum AlleleCountArrayPolicy {
    /// bcftools 1.17 reuses slots without clearing the count array.
    LegacyReuse,
    /// Normative counting starts an empty table for each output record.
    ResetPerTable,
}

/// Governed non-default handling for a command-line usage error.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum UsageErrorPolicy {
    /// Historical hap.py commands reported invalid arguments as failure.
    LegacyFailure,
    /// Historical pre.py and qfy.py reported invalid arguments as success.
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

/// Select the legacy usage-error contract at the CLI adapter boundary.
pub(crate) fn usage_error_policy(arguments: &[std::ffi::OsString]) -> Option<UsageErrorPolicy> {
    match arguments.get(1).and_then(|value| value.to_str())? {
        "germline" | "compare" => Some(UsageErrorPolicy::LegacyFailure),
        "pre" | "preprocess" | "prepy" | "quantify" | "qfy" => {
            Some(UsageErrorPolicy::DeprecatedSuccess)
        }
        _ => None,
    }
}

/// Emit warnings for parsed options that exist only for migration.
pub(crate) fn emit_deprecation_warnings(command: &Command) {
    if let Command::Germline(args) = command
        && (args.engine_vcfeval.is_some() || args.engine_vcfeval_template.is_some())
    {
        eprintln!("{VCFEVAL_RUNTIME_DEPRECATION}");
    }
}

/// Map a record to governed location streams before block splitting.
///
/// The returned stream IDs are the complete compatibility decision consumed by
/// preprocessing: downstream partition expansion is generic and does not
/// inspect location overlap or choose whether records are duplicated.
pub(crate) fn location_stream_groups(
    policy: LocationStreamPolicy,
    filters: &[LocationFilter],
    chrom: &str,
    pos: usize,
) -> Vec<usize> {
    let matching = filters
        .iter()
        .enumerate()
        .filter_map(|(index, filter)| filter.matches(chrom, pos).then_some(index))
        .collect::<Vec<_>>();
    match policy {
        LocationStreamPolicy::IndependentLegacyStreams => matching,
        LocationStreamPolicy::SetUnion if matching.is_empty() => Vec::new(),
        LocationStreamPolicy::SetUnion => vec![0],
    }
}

/// Construct an SCMP end coordinate under the selected compatibility policy.
pub(crate) fn scmp_refvar_end(
    policy: ScmpRefVarSpanPolicy,
    start: i64,
    reference_len: usize,
    alt_len: usize,
) -> anyhow::Result<i64> {
    let span = match policy {
        ScmpRefVarSpanPolicy::LegacyAltLength => alt_len,
        ScmpRefVarSpanPolicy::ReferenceLength => reference_len,
    };
    Ok(start + i64::try_from(span)? - 1)
}

/// Begin a bcftools-compatible allele count table.
pub(crate) fn begin_allele_count_table(policy: AlleleCountArrayPolicy, counts: &mut Vec<usize>) {
    if policy == AlleleCountArrayPolicy::ResetPerTable {
        counts.clear();
    }
}

/// Record one allele occurrence in the governed count-array representation.
pub(crate) fn record_allele_count(counts: &mut Vec<usize>, slot: usize, first_record: bool) {
    if counts.len() <= slot {
        counts.resize(slot + 1, 0);
    }
    if first_record {
        counts[slot] = 1;
    } else {
        counts[slot] += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vcf;

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
    fn legacy_only_overlapping_locations_expand_to_independent_streams() {
        let filters = [
            vcf::LocationFilter::Range {
                chrom: "chr1".into(),
                start: 1,
                end: 100,
            },
            vcf::LocationFilter::Range {
                chrom: "chr1".into(),
                start: 51,
                end: 120,
            },
        ];
        assert_eq!(
            location_stream_groups(
                LocationStreamPolicy::IndependentLegacyStreams,
                &filters,
                "chr1",
                75,
            ),
            [0, 1]
        );
    }

    #[test]
    fn normative_set_union_collapses_overlapping_locations() {
        let filters = [
            vcf::LocationFilter::Contig("chr1".into()),
            vcf::LocationFilter::Contig("chr1".into()),
        ];
        assert_eq!(
            location_stream_groups(LocationStreamPolicy::SetUnion, &filters, "chr1", 75),
            [0]
        );
        assert!(
            location_stream_groups(LocationStreamPolicy::SetUnion, &filters, "chr2", 75).is_empty()
        );
    }

    #[test]
    fn legacy_only_scmp_refvar_span_uses_alt_length() {
        assert_eq!(
            scmp_refvar_end(ScmpRefVarSpanPolicy::LegacyAltLength, 9, 1, 3).unwrap(),
            11
        );
    }

    #[test]
    fn normative_scmp_refvar_span_uses_reference_length() {
        assert_eq!(
            scmp_refvar_end(ScmpRefVarSpanPolicy::ReferenceLength, 9, 1, 3).unwrap(),
            9
        );
    }

    #[test]
    fn legacy_only_bcftools_count_table_reuses_stale_slots() {
        let mut counts = vec![2, 7];
        begin_allele_count_table(AlleleCountArrayPolicy::LegacyReuse, &mut counts);
        record_allele_count(&mut counts, 0, true);
        assert_eq!(counts, [1, 7]);
    }

    #[test]
    fn normative_allele_count_table_clears_all_slots() {
        let mut counts = vec![2, 7];
        begin_allele_count_table(AlleleCountArrayPolicy::ResetPerTable, &mut counts);
        record_allele_count(&mut counts, 0, true);
        assert_eq!(counts, [1]);
    }

    #[test]
    fn vcfeval_deprecation_warning_has_exact_replacement_and_deadline() {
        assert_eq!(
            VCFEVAL_RUNTIME_DEPRECATION,
            "warning: --engine-vcfeval-path and --engine-vcfeval-template are deprecated and ignored; use --engine vcfeval --reference <FASTA>; these options will be removed in hap-rs 1.0.0"
        );
    }
}
