use super::validated::{QueryProvenance, ValidatedVcfRecord};
use super::variant::RawVcfRecord;

/// Counts accumulated for one side and classification of a benchmark.
#[derive(Clone, Debug, Default)]
pub(crate) struct CountsBucket {
    pub(crate) total: usize,
    pub(crate) ti: usize,
    pub(crate) tv: usize,
    pub(crate) het: usize,
    pub(crate) homalt: usize,
}

/// Truth and query counts for one variant grouping.
#[derive(Clone, Debug, Default)]
pub(crate) struct TypeCounts {
    pub(crate) truth_total: CountsBucket,
    pub(crate) truth_tp: CountsBucket,
    pub(crate) truth_fn: CountsBucket,
    pub(crate) query_total: CountsBucket,
    pub(crate) query_tp: CountsBucket,
    pub(crate) query_fp: CountsBucket,
    pub(crate) query_unk: CountsBucket,
}

/// False-positive subclass carried on an FP row.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum FpClass {
    /// Genotype mismatch (`am` block kind).
    Gt,
    /// Allele mismatch (`lm` block kind).
    Al,
}

impl FpClass {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Gt => "gt",
            Self::Al => "al",
        }
    }
}

/// Legacy xcmp block-level comparison classification.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum XcmpCtype {
    SimpleMatch,
    SimpleMismatch,
    HapMatch,
    HapMismatch,
    HapfailMismatch,
}

impl XcmpCtype {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::SimpleMatch => "simple:match",
            Self::SimpleMismatch => "simple:mismatch",
            Self::HapMatch => "hap:match",
            Self::HapMismatch => "hap:mismatch",
            Self::HapfailMismatch => "hapfail:mismatch",
        }
    }
}

/// Lexicographic ordering key for a comparison row. Field order is the sort
/// precedence and must not be reordered.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct SortKey {
    pub(crate) chrom: String,
    pub(crate) pos: usize,
    pub(crate) side_rank: usize,
    pub(crate) type_rank: usize,
}

impl SortKey {
    pub(crate) fn new(chrom: String, pos: usize, side_rank: usize, type_rank: usize) -> Self {
        Self {
            chrom,
            pos,
            side_rank,
            type_rank,
        }
    }
}

/// A comparison record plus the domain facts needed by report engines.
#[derive(Clone, Debug)]
pub(crate) struct AnnotatedRow {
    pub(crate) sort_key: SortKey,
    pub(crate) record: ComparisonRecord,
    /// Whether the originating query variant was PASS-filtered.
    pub(crate) query_pass: bool,
    /// False-positive subclass (`gt` or `al`) when this is an FP row.
    pub(crate) fp_class: Option<FpClass>,
    /// Legacy xcmp block-level comparison classification.
    pub(crate) xcmp_ctype: Option<XcmpCtype>,
    pub(crate) xcmp_hap_match: bool,
}

/// A lossless VCF record whose coordinate, alleles, genotypes, and provenance
/// were checked before entering comparison and reporting code.
#[derive(Clone, Debug)]
pub(crate) struct ComparisonRecord {
    record: ValidatedVcfRecord,
}

impl ComparisonRecord {
    pub(crate) fn checked(raw: RawVcfRecord) -> Self {
        let record = ValidatedVcfRecord::try_from_raw(raw, QueryProvenance::Unavailable)
            .expect("comparison rows preserve checked VCF invariants");
        Self { record }
    }

    pub(crate) fn raw(&self) -> &RawVcfRecord {
        self.record.raw()
    }

    pub(crate) fn into_validated(self) -> ValidatedVcfRecord {
        self.record
    }

    pub(crate) fn try_update<R>(
        &mut self,
        edit: impl FnOnce(&mut RawVcfRecord) -> anyhow::Result<R>,
    ) -> anyhow::Result<R> {
        self.record.try_update(edit)
    }

    pub(crate) fn samples_contain(&self, value: &str) -> bool {
        self.raw()
            .samples
            .iter()
            .any(|sample| sample.contains(value))
    }

    pub(crate) fn replace_sample_fragment(&mut self, from: &str, to: &str) -> anyhow::Result<()> {
        self.try_update(|record| {
            for sample in &mut record.samples {
                *sample = sample.replace(from, to);
            }
            Ok(())
        })
    }

    #[cfg(test)]
    pub(crate) fn fixture(line: &str) -> Self {
        Self::checked(
            RawVcfRecord::from_line(line, std::path::Path::new("comparison-test-fixture"))
                .expect("comparison fixture must be a valid VCF record"),
        )
    }
}

impl From<RawVcfRecord> for ComparisonRecord {
    fn from(record: RawVcfRecord) -> Self {
        Self::checked(record)
    }
}

impl std::ops::Deref for ComparisonRecord {
    type Target = RawVcfRecord;

    fn deref(&self) -> &Self::Target {
        self.raw()
    }
}

impl PartialEq for ComparisonRecord {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other).is_eq()
    }
}

impl Eq for ComparisonRecord {}

impl PartialOrd for ComparisonRecord {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ComparisonRecord {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let left = self.raw();
        let right = other.raw();
        left.chrom
            .cmp(&right.chrom)
            .then(left.pos.cmp(&right.pos))
            .then(left.id.cmp(&right.id))
            .then(left.ref_allele.cmp(&right.ref_allele))
            .then(left.alt_allele.cmp(&right.alt_allele))
            .then(left.qual.cmp(&right.qual))
            .then(left.filter.cmp(&right.filter))
            .then(left.info.cmp(&right.info))
            .then(left.format.cmp(&right.format))
            .then(left.samples.cmp(&right.samples))
    }
}
