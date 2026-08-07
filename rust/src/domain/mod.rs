//! Filesystem-free data shared by comparison and quantification engines.

pub(crate) mod cephes;
mod interval;
mod metrics;
mod statistics;
mod variant;

pub(crate) use interval::Interval;
pub(crate) use metrics::{AnnotatedRow, CountsBucket, TypeCounts};
pub(crate) use statistics::jeffreys_interval;
pub(crate) use variant::RawVcfRecord;
