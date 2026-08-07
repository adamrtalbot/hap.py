//! Filesystem-free data shared by comparison and quantification engines.

pub(crate) mod cephes;
mod interval;
mod metrics;
mod statistics;
mod validated;
mod variant;

pub(crate) use interval::Interval;
pub(crate) use metrics::{AnnotatedRow, ComparisonRecord, CountsBucket, TypeCounts};
pub(crate) use statistics::jeffreys_interval;
pub(crate) use validated::{
    Allele, AlleleIndex, GenomicCoordinate, GenomicPosition, Genotype, GenotypePhase, OutputPlan,
    Ploidy, QueryProvenance, QueryRecordIndex, SampleDecision, VariantOutputFormat,
};
pub(crate) use variant::RawVcfRecord;
