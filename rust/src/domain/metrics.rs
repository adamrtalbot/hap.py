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

/// A comparison record plus the domain facts needed by report engines.
#[derive(Clone, Debug)]
pub(crate) struct AnnotatedRow {
    pub(crate) sort_key: (String, usize, usize, usize),
    pub(crate) record: ComparisonRecord,
    /// Whether the originating query variant was PASS-filtered.
    pub(crate) query_pass: bool,
    /// False-positive subclass (`gt` or `al`) when this is an FP row.
    pub(crate) fp_class: Option<&'static str>,
    /// Legacy xcmp block-level comparison classification.
    pub(crate) xcmp_ctype: Option<&'static str>,
    pub(crate) xcmp_hap_match: bool,
}

/// A lossless VCF record whose coordinate, alleles, genotypes, and provenance
/// were checked before entering comparison and reporting code.
#[derive(Clone, Debug)]
pub(crate) struct ComparisonRecord {
    record: ValidatedVcfRecord,
}

impl ComparisonRecord {
    pub(crate) fn checked(raw: super::RawVcfRecord) -> Self {
        let record = ValidatedVcfRecord::try_from_raw(raw, super::QueryProvenance::Unavailable)
            .expect("comparison rows preserve checked VCF invariants");
        Self { record }
    }

    pub(crate) fn raw(&self) -> &super::RawVcfRecord {
        self.record.raw()
    }

    pub(crate) fn validated(&self) -> &ValidatedVcfRecord {
        &self.record
    }

    pub(crate) fn into_validated(self) -> ValidatedVcfRecord {
        self.record
    }

    pub(crate) fn try_update<R>(
        &mut self,
        edit: impl FnOnce(&mut super::RawVcfRecord) -> anyhow::Result<R>,
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
            super::RawVcfRecord::from_line(line, Path::new("comparison-test-fixture"))
                .expect("comparison fixture must be a valid VCF record"),
        )
    }
}

impl From<super::RawVcfRecord> for ComparisonRecord {
    fn from(record: super::RawVcfRecord) -> Self {
        Self::checked(record)
    }
}

impl std::ops::Deref for ComparisonRecord {
    type Target = super::RawVcfRecord;

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
use crate::adapters::vcf::ValidatedVcfRecord;
use std::path::Path;
