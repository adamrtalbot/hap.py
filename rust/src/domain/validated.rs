//! Validated value objects shared by application use cases.
//!
//! Textual CLI values and VCF fields should be converted to these types at
//! their adapters. The application core can then rely on the invariants
//! documented by each type without repeating format-specific validation.

use super::variant::RawVcfRecord;
use anyhow::{Context, Result};
use std::error::Error;
use std::fmt;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// An error raised while constructing a domain value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum DomainError {
    /// A genomic position was zero even though VCF coordinates are one-based.
    ZeroPosition,
    /// A contig name was empty.
    EmptyContig,
    /// A contig name contained whitespace.
    InvalidContig(String),
    /// A coordinate did not use the `contig:position` form.
    InvalidCoordinate(String),
    /// An allele was empty or represented a missing VCF value.
    EmptyAllele,
    /// An allele contained a field or allele-list delimiter.
    InvalidAllele(String),
    /// An allele index did not identify the reference or an alternate allele.
    AlleleIndexOutOfBounds {
        /// The invalid VCF allele index.
        index: usize,
        /// The number of alternate alleles available in the record.
        alternate_allele_count: usize,
    },
    /// Ploidy was zero.
    ZeroPloidy,
    /// A genotype contained no allele slots.
    EmptyGenotype,
    /// A genotype was malformed or used inconsistent phasing separators.
    InvalidGenotype(String),
    /// A phased genotype contained only one allele slot.
    PhasedHaploidGenotype,
    /// Query provenance was neither missing nor a non-negative record index.
    #[allow(dead_code, reason = "constructed by the textual provenance adapter")]
    InvalidQueryProvenance(String),
    /// Query provenance referred beyond the query record collection.
    QuerySourceOutOfBounds {
        /// The invalid zero-based query record index.
        index: usize,
        /// The number of query records available.
        query_record_count: usize,
    },
    /// An output report prefix was empty.
    EmptyOutputPrefix,
}

impl fmt::Display for DomainError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroPosition => write!(formatter, "genomic positions must be at least 1"),
            Self::EmptyContig => write!(formatter, "contig name cannot be empty"),
            Self::InvalidContig(contig) => {
                write!(
                    formatter,
                    "contig name {contig:?} cannot contain whitespace"
                )
            }
            Self::InvalidCoordinate(value) => write!(
                formatter,
                "invalid genomic coordinate {value:?}; expected contig:position with a 1-based position"
            ),
            Self::EmptyAllele => write!(formatter, "allele cannot be empty or '.'"),
            Self::InvalidAllele(allele) => write!(
                formatter,
                "invalid allele {allele:?}; a single allele cannot contain whitespace, a comma, or a tab"
            ),
            Self::AlleleIndexOutOfBounds {
                index,
                alternate_allele_count,
            } => write!(
                formatter,
                "allele index {index} is out of bounds for {alternate_allele_count} alternate allele(s)"
            ),
            Self::ZeroPloidy => write!(formatter, "genotype ploidy must be at least 1"),
            Self::EmptyGenotype => write!(formatter, "genotype must contain at least one allele"),
            Self::InvalidGenotype(value) => write!(
                formatter,
                "invalid genotype {value:?}; use allele indices or '.', separated consistently by '/' or '|'"
            ),
            Self::PhasedHaploidGenotype => {
                write!(formatter, "a haploid genotype cannot be marked as phased")
            }
            Self::InvalidQueryProvenance(value) => write!(
                formatter,
                "invalid query provenance {value:?}; expected '.' or a zero-based record index"
            ),
            Self::QuerySourceOutOfBounds {
                index,
                query_record_count,
            } => write!(
                formatter,
                "query provenance index {index} is out of bounds for {query_record_count} query record(s)"
            ),
            Self::EmptyOutputPrefix => write!(formatter, "output report prefix cannot be empty"),
        }
    }
}

impl Error for DomainError {}

/// A one-based genomic position.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct GenomicPosition(NonZeroUsize);

impl GenomicPosition {
    /// Constructs a position, rejecting zero.
    pub(crate) fn new(value: usize) -> Result<Self, DomainError> {
        NonZeroUsize::new(value)
            .map(Self)
            .ok_or(DomainError::ZeroPosition)
    }

    /// Returns the one-based numeric position.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain accessor"))]
    pub(crate) const fn get(self) -> usize {
        self.0.get()
    }
}

impl fmt::Display for GenomicPosition {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

impl FromStr for GenomicPosition {
    type Err = DomainError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        value
            .parse::<usize>()
            .map_err(|_| DomainError::InvalidCoordinate(value.to_owned()))
            .and_then(Self::new)
    }
}

/// A contig and one-based position identifying a genomic locus.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct GenomicCoordinate {
    contig: String,
    position: GenomicPosition,
}

impl GenomicCoordinate {
    /// Constructs a coordinate from a contig name and checked position.
    pub(crate) fn new(
        contig: impl Into<String>,
        position: GenomicPosition,
    ) -> Result<Self, DomainError> {
        let contig = contig.into();
        if contig.is_empty() {
            return Err(DomainError::EmptyContig);
        }
        if contig.chars().any(char::is_whitespace) {
            return Err(DomainError::InvalidContig(contig));
        }
        Ok(Self { contig, position })
    }

    /// Returns the contig name.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain accessor"))]
    pub(crate) fn contig(&self) -> &str {
        &self.contig
    }

    /// Returns the one-based genomic position.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain accessor"))]
    pub(crate) const fn position(&self) -> GenomicPosition {
        self.position
    }
}

impl fmt::Display for GenomicCoordinate {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}:{}", self.contig, self.position)
    }
}

impl FromStr for GenomicCoordinate {
    type Err = DomainError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let (contig, position) = value
            .rsplit_once(':')
            .ok_or_else(|| DomainError::InvalidCoordinate(value.to_owned()))?;
        let position = position
            .parse::<usize>()
            .map_err(|_| DomainError::InvalidCoordinate(value.to_owned()))?;
        Self::new(contig, GenomicPosition::new(position)?)
    }
}

/// A single called allele, rather than a comma-delimited VCF ALT field.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Allele(String);

impl Allele {
    /// Constructs an allele suitable for use as one REF or ALT value.
    ///
    /// Symbolic, breakend, and spanning-deletion alleles are retained. A VCF
    /// missing value (`.`), delimiters, and whitespace are rejected.
    pub(crate) fn new(value: impl Into<String>) -> Result<Self, DomainError> {
        let value = value.into();
        if value.is_empty() || value == "." {
            return Err(DomainError::EmptyAllele);
        }
        if value.contains(',') || value.chars().any(char::is_whitespace) {
            return Err(DomainError::InvalidAllele(value));
        }
        Ok(Self(value))
    }

    /// Returns the allele text.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain accessor"))]
    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for Allele {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

impl FromStr for Allele {
    type Err = DomainError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::new(value)
    }
}

/// A checked VCF genotype allele index (`0` is the reference allele).
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct AlleleIndex(usize);

impl AlleleIndex {
    /// Constructs an index for a record with the given number of ALT alleles.
    pub(crate) fn new(value: usize, alternate_allele_count: usize) -> Result<Self, DomainError> {
        if value > alternate_allele_count {
            return Err(DomainError::AlleleIndexOutOfBounds {
                index: value,
                alternate_allele_count,
            });
        }
        Ok(Self(value))
    }

    /// Returns the reference allele index, which is always zero.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain constructor"))]
    pub(crate) const fn reference() -> Self {
        Self(0)
    }

    /// Returns the numeric VCF allele index.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain accessor"))]
    pub(crate) const fn get(self) -> usize {
        self.0
    }
}

impl fmt::Display for AlleleIndex {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// The number of chromosome copies represented by a genotype.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct Ploidy(NonZeroUsize);

impl Ploidy {
    /// Constructs a ploidy, rejecting zero.
    pub(crate) fn new(value: usize) -> Result<Self, DomainError> {
        NonZeroUsize::new(value)
            .map(Self)
            .ok_or(DomainError::ZeroPloidy)
    }

    /// Returns the number of allele slots.
    pub(crate) const fn get(self) -> usize {
        self.0.get()
    }
}

impl fmt::Display for Ploidy {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Whether the allele ordering of a polyploid genotype is phased.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum GenotypePhase {
    /// Allele copies have no asserted haplotype ordering (`/`).
    Unphased,
    /// Allele copies have asserted haplotype ordering (`|`).
    Phased,
}

/// A genotype with checked allele indices and explicit phase and ploidy.
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub(crate) struct Genotype {
    alleles: Box<[Option<AlleleIndex>]>,
    phase: GenotypePhase,
    ploidy: Ploidy,
}

impl Genotype {
    /// Constructs a genotype from already-checked allele indices.
    pub(crate) fn new(
        alleles: Vec<Option<AlleleIndex>>,
        phase: GenotypePhase,
    ) -> Result<Self, DomainError> {
        let ploidy = Ploidy::new(alleles.len()).map_err(|_| DomainError::EmptyGenotype)?;
        if ploidy.get() == 1 && phase == GenotypePhase::Phased {
            return Err(DomainError::PhasedHaploidGenotype);
        }
        Ok(Self {
            alleles: alleles.into_boxed_slice(),
            phase,
            ploidy,
        })
    }

    /// Parses a VCF GT value against the number of ALT alleles in its record.
    pub(crate) fn parse(value: &str, alternate_allele_count: usize) -> Result<Self, DomainError> {
        if value.is_empty() {
            return Err(DomainError::EmptyGenotype);
        }
        let has_unphased = value.contains('/');
        let has_phased = value.contains('|');
        if has_unphased && has_phased {
            return Err(DomainError::InvalidGenotype(value.to_owned()));
        }
        let (phase, parts): (GenotypePhase, Vec<&str>) = if has_phased {
            (GenotypePhase::Phased, value.split('|').collect())
        } else if has_unphased {
            (GenotypePhase::Unphased, value.split('/').collect())
        } else {
            (GenotypePhase::Unphased, vec![value])
        };
        if parts.iter().any(|part| part.is_empty()) {
            return Err(DomainError::InvalidGenotype(value.to_owned()));
        }
        let alleles = parts
            .into_iter()
            .map(|part| {
                if part == "." {
                    Ok(None)
                } else {
                    let index = part
                        .parse::<usize>()
                        .map_err(|_| DomainError::InvalidGenotype(value.to_owned()))?;
                    AlleleIndex::new(index, alternate_allele_count).map(Some)
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        Self::new(alleles, phase)
    }

    /// Returns each called or missing allele slot.
    /// Returns the genotype ploidy.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain accessor"))]
    pub(crate) const fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    /// Returns the genotype phasing state.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain accessor"))]
    pub(crate) const fn phase(&self) -> GenotypePhase {
        self.phase
    }

    /// Returns whether every allele slot is missing.
    #[cfg_attr(not(test), allow(dead_code, reason = "domain predicate"))]
    pub(crate) fn is_fully_missing(&self) -> bool {
        self.alleles.iter().all(Option::is_none)
    }
}

impl fmt::Display for Genotype {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let separator = match self.phase {
            GenotypePhase::Unphased => '/',
            GenotypePhase::Phased => '|',
        };
        for (index, allele) in self.alleles.iter().enumerate() {
            if index > 0 {
                separator.fmt(formatter)?;
            }
            match allele {
                Some(allele) => allele.fmt(formatter)?,
                None => formatter.write_str(".")?,
            }
        }
        Ok(())
    }
}

/// A checked zero-based index into the original query record collection.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(crate) struct QueryRecordIndex(usize);

impl QueryRecordIndex {
    /// Constructs an index checked against the number of query records.
    pub(crate) fn new(index: usize, query_record_count: usize) -> Result<Self, DomainError> {
        if index >= query_record_count {
            return Err(DomainError::QuerySourceOutOfBounds {
                index,
                query_record_count,
            });
        }
        Ok(Self(index))
    }

    /// Returns the zero-based record index.
    pub(crate) const fn get(self) -> usize {
        self.0
    }
}

impl fmt::Display for QueryRecordIndex {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.0.fmt(formatter)
    }
}

/// Provenance linking an evaluated record to its original query record.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum QueryProvenance {
    /// The evaluated record has no query source, for example a truth-only call.
    Unavailable,
    /// The zero-based index of the original query record.
    SourceRecord(QueryRecordIndex),
}

impl QueryProvenance {
    /// Constructs checked provenance for a query record collection.
    pub(crate) fn source(index: usize, query_record_count: usize) -> Result<Self, DomainError> {
        QueryRecordIndex::new(index, query_record_count).map(Self::SourceRecord)
    }

    /// Parses `.` or a zero-based source index and checks it against the query.
    #[cfg_attr(not(test), allow(dead_code, reason = "file adapter parser"))]
    pub(crate) fn parse(value: &str, query_record_count: usize) -> Result<Self, DomainError> {
        if value == "." {
            return Ok(Self::Unavailable);
        }
        let index = value
            .parse::<usize>()
            .map_err(|_| DomainError::InvalidQueryProvenance(value.to_owned()))?;
        Self::source(index, query_record_count)
    }

    /// Returns the original query record index, when present.
    pub(crate) const fn source_index(self) -> Option<usize> {
        match self {
            Self::Unavailable => None,
            Self::SourceRecord(index) => Some(index.get()),
        }
    }
}

impl fmt::Display for QueryProvenance {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unavailable => formatter.write_str("."),
            Self::SourceRecord(index) => index.fmt(formatter),
        }
    }
}

/// The encoding selected for an annotated variant output.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub(crate) enum VariantOutputFormat {
    /// Block-gzipped VCF (`.vcf.gz`).
    Vcf,
    /// Binary variant call format (`.bcf`).
    Bcf,
}

impl VariantOutputFormat {
    /// Returns the report suffix for the selected encoding.
    pub(crate) const fn suffix(self) -> &'static str {
        match self {
            Self::Vcf => "vcf.gz",
            Self::Bcf => "bcf",
        }
    }
}

/// A validated plan describing which report artifacts an invocation emits.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OutputPlan {
    report_prefix: PathBuf,
    write_counts: bool,
    write_metrics: bool,
    variant_output: Option<VariantOutputFormat>,
}

impl OutputPlan {
    /// Constructs an output plan. The summary report is always included.
    pub(crate) fn new(
        report_prefix: impl Into<PathBuf>,
        write_counts: bool,
        write_metrics: bool,
        variant_output: Option<VariantOutputFormat>,
    ) -> Result<Self, DomainError> {
        let report_prefix = report_prefix.into();
        if report_prefix.as_os_str().is_empty() {
            return Err(DomainError::EmptyOutputPrefix);
        }
        Ok(Self {
            report_prefix,
            write_counts,
            write_metrics,
            variant_output,
        })
    }

    /// Returns the always-present compact summary path.
    pub(crate) fn summary_path(&self) -> PathBuf {
        self.suffixed_path("summary.csv")
    }

    /// Returns the extended counts path when it is enabled.
    pub(crate) fn counts_path(&self) -> Option<PathBuf> {
        self.write_counts
            .then(|| self.suffixed_path("extended.csv"))
    }

    /// Returns the compressed metrics path when it is enabled.
    pub(crate) fn metrics_path(&self) -> Option<PathBuf> {
        self.write_metrics
            .then(|| self.suffixed_path("metrics.json.gz"))
    }

    /// Returns the annotated variant path when it is enabled.
    pub(crate) fn variant_path(&self) -> Option<PathBuf> {
        self.variant_output
            .map(|format| self.suffixed_path(format.suffix()))
    }

    /// Returns whether any output would overwrite the supplied input path.
    pub(crate) fn conflicts_with_input(&self, input: &Path) -> bool {
        self.summary_path() == input
            || self.counts_path().as_deref() == Some(input)
            || self.metrics_path().as_deref() == Some(input)
            || self.variant_path().as_deref() == Some(input)
    }

    fn suffixed_path(&self, suffix: &str) -> PathBuf {
        let mut path = self.report_prefix.as_os_str().to_os_string();
        path.push(".");
        path.push(suffix);
        PathBuf::from(path)
    }
}

/// A lossless VCF record whose coordinate, alleles, genotypes, and provenance
/// have been checked before entering application and engine code.
#[derive(Clone, Debug)]
pub(crate) struct ValidatedVcfRecord {
    raw: RawVcfRecord,
    coordinate: GenomicCoordinate,
    reference: Allele,
    alternates: Vec<Allele>,
    genotypes: Vec<Option<Genotype>>,
    provenance: QueryProvenance,
}

impl ValidatedVcfRecord {
    /// Converts a raw adapter record into checked domain values.
    pub(crate) fn try_from_raw(raw: RawVcfRecord, provenance: QueryProvenance) -> Result<Self> {
        let position = GenomicPosition::new(raw.pos).map_err(|error| {
            anyhow::anyhow!(
                "invalid VCF position {} on contig {}: {error}",
                raw.pos,
                raw.chrom
            )
        })?;
        let coordinate = GenomicCoordinate::new(raw.chrom.clone(), position)
            .with_context(|| format!("invalid VCF coordinate {}:{}", raw.chrom, raw.pos))?;
        let reference = Allele::new(raw.ref_allele.clone())
            .with_context(|| format!("invalid REF allele at {coordinate}"))?;
        let alternates = if raw.alt_allele == "." {
            Vec::new()
        } else {
            raw.alt_allele
                .split(',')
                .map(|allele| {
                    Allele::new(allele)
                        .with_context(|| format!("invalid ALT allele {allele:?} at {coordinate}"))
                })
                .collect::<Result<Vec<_>>>()?
        };
        let format_keys = raw.format_keys();
        let gt_index = format_keys.iter().position(|key| *key == "GT");
        let genotypes = raw
            .samples
            .iter()
            .map(|sample| {
                let Some(gt_index) = gt_index else {
                    return Ok(None);
                };
                let Some(gt) = sample.split(':').nth(gt_index) else {
                    return Ok(None);
                };
                Genotype::parse(gt, alternates.len())
                    .map(Some)
                    .map_err(|error| anyhow::anyhow!("invalid GT {gt:?} at {coordinate}: {error}"))
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            raw,
            coordinate,
            reference,
            alternates,
            genotypes,
            provenance,
        })
    }

    #[cfg_attr(not(test), allow(dead_code, reason = "checked domain accessor"))]
    pub(crate) fn coordinate(&self) -> &GenomicCoordinate {
        &self.coordinate
    }

    #[cfg_attr(not(test), allow(dead_code, reason = "checked domain accessor"))]
    pub(crate) fn reference(&self) -> &Allele {
        &self.reference
    }

    #[cfg_attr(not(test), allow(dead_code, reason = "checked domain accessor"))]
    pub(crate) fn alternates(&self) -> &[Allele] {
        &self.alternates
    }

    #[cfg_attr(not(test), allow(dead_code, reason = "checked domain accessor"))]
    pub(crate) fn genotypes(&self) -> &[Option<Genotype>] {
        &self.genotypes
    }

    pub(crate) const fn provenance(&self) -> QueryProvenance {
        self.provenance
    }

    pub(crate) fn raw(&self) -> &RawVcfRecord {
        &self.raw
    }

    pub(crate) fn into_raw(self) -> RawVcfRecord {
        self.raw
    }

    /// Applies an edit transactionally and refreshes all checked facts.
    pub(crate) fn try_update<R>(
        &mut self,
        edit: impl FnOnce(&mut RawVcfRecord) -> Result<R>,
    ) -> Result<R> {
        let mut raw = self.raw.clone();
        let result = edit(&mut raw)?;
        *self = Self::try_from_raw(raw, self.provenance)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positions_and_coordinates_reject_invalid_states_and_round_trip() {
        assert_eq!(GenomicPosition::new(0), Err(DomainError::ZeroPosition));
        assert_eq!(
            "chr1:0".parse::<GenomicCoordinate>(),
            Err(DomainError::ZeroPosition)
        );
        assert!(matches!(
            GenomicCoordinate::new("chr 1", GenomicPosition::new(1).unwrap()),
            Err(DomainError::InvalidContig(_))
        ));
        assert!(matches!(
            "chr1".parse::<GenomicCoordinate>(),
            Err(DomainError::InvalidCoordinate(_))
        ));

        let coordinate = "HLA-A:29910247".parse::<GenomicCoordinate>().unwrap();
        assert_eq!(coordinate.contig(), "HLA-A");
        assert_eq!(coordinate.position().get(), 29_910_247);
        assert_eq!(coordinate.to_string(), "HLA-A:29910247");
    }

    #[test]
    fn alleles_accept_vcf_forms_but_not_missing_or_delimited_values() {
        for value in ["A", "ACGT", "*", "<DEL>", "A]chr2:42]"] {
            let allele = Allele::new(value).unwrap();
            assert_eq!(allele.to_string(), value);
        }
        assert_eq!(Allele::new(""), Err(DomainError::EmptyAllele));
        assert_eq!(Allele::new("."), Err(DomainError::EmptyAllele));
        assert!(matches!(
            Allele::new("A,C"),
            Err(DomainError::InvalidAllele(_))
        ));
        assert!(matches!(
            Allele::new("A C"),
            Err(DomainError::InvalidAllele(_))
        ));
    }

    #[test]
    fn allele_indices_are_checked_against_the_record() {
        assert_eq!(AlleleIndex::new(0, 0).unwrap(), AlleleIndex::reference());
        assert_eq!(AlleleIndex::new(2, 2).unwrap().get(), 2);
        assert_eq!(
            AlleleIndex::new(3, 2),
            Err(DomainError::AlleleIndexOutOfBounds {
                index: 3,
                alternate_allele_count: 2,
            })
        );
    }

    #[test]
    fn genotypes_preserve_phase_ploidy_missing_calls_and_display() {
        let cases = [
            ("0", 1, GenotypePhase::Unphased, 1),
            (".", 1, GenotypePhase::Unphased, 1),
            ("0/1", 1, GenotypePhase::Unphased, 2),
            ("1|0", 1, GenotypePhase::Phased, 2),
            ("0/1/2", 2, GenotypePhase::Unphased, 3),
            (".|1", 1, GenotypePhase::Phased, 2),
        ];
        for (text, alt_count, phase, ploidy) in cases {
            let genotype = Genotype::parse(text, alt_count).unwrap();
            assert_eq!(genotype.phase(), phase);
            assert_eq!(genotype.ploidy().get(), ploidy);
            assert_eq!(genotype.to_string(), text);
        }
        assert!(Genotype::parse("./.", 1).unwrap().is_fully_missing());
    }

    #[test]
    fn genotypes_reject_malformed_and_out_of_range_calls() {
        assert_eq!(Genotype::parse("", 1), Err(DomainError::EmptyGenotype));
        assert!(matches!(
            Genotype::parse("0/1|1", 1),
            Err(DomainError::InvalidGenotype(_))
        ));
        assert!(matches!(
            Genotype::parse("0/", 1),
            Err(DomainError::InvalidGenotype(_))
        ));
        assert!(matches!(
            Genotype::parse("-1/0", 1),
            Err(DomainError::InvalidGenotype(_))
        ));
        assert_eq!(
            Genotype::parse("0/2", 1),
            Err(DomainError::AlleleIndexOutOfBounds {
                index: 2,
                alternate_allele_count: 1,
            })
        );
        assert_eq!(
            Genotype::new(vec![Some(AlleleIndex::reference())], GenotypePhase::Phased),
            Err(DomainError::PhasedHaploidGenotype)
        );
    }

    #[test]
    fn query_provenance_checks_bounds_and_round_trips() {
        let source = QueryProvenance::parse("2", 3).unwrap();
        assert_eq!(source.source_index(), Some(2));
        assert_eq!(source.to_string(), "2");
        assert_eq!(QueryProvenance::parse(".", 0).unwrap().source_index(), None);
        assert_eq!(
            QueryProvenance::source(3, 3),
            Err(DomainError::QuerySourceOutOfBounds {
                index: 3,
                query_record_count: 3,
            })
        );
        assert!(matches!(
            QueryProvenance::parse("source-2", 3),
            Err(DomainError::InvalidQueryProvenance(_))
        ));
    }

    #[test]
    fn output_plan_derives_paths_without_creating_outputs() {
        assert_eq!(
            OutputPlan::new("", true, true, Some(VariantOutputFormat::Vcf)),
            Err(DomainError::EmptyOutputPrefix)
        );
        let plan = OutputPlan::new(
            "reports/sample.v1",
            true,
            true,
            Some(VariantOutputFormat::Bcf),
        )
        .unwrap();
        assert_eq!(
            plan.summary_path(),
            Path::new("reports/sample.v1.summary.csv")
        );
        assert_eq!(
            plan.counts_path().as_deref(),
            Some(Path::new("reports/sample.v1.extended.csv"))
        );
        assert_eq!(
            plan.metrics_path().as_deref(),
            Some(Path::new("reports/sample.v1.metrics.json.gz"))
        );
        assert_eq!(
            plan.variant_path().as_deref(),
            Some(Path::new("reports/sample.v1.bcf"))
        );
        assert!(plan.conflicts_with_input(Path::new("reports/sample.v1.bcf")));
        assert!(!plan.conflicts_with_input(Path::new("input.bcf")));
    }
}
