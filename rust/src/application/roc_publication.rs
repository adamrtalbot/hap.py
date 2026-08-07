//! Publication adapter for filesystem-neutral ROC engine artifacts.

use crate::adapters::report::suffixed_report_path;
use crate::domain::AnnotatedRow;
use crate::engines::roc::{self, MetricIndices, RocOptions};
use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::io::Write;
use std::path::Path;

pub(crate) fn write_roc_files(
    prefix: &Path,
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
) -> Result<MetricIndices> {
    publish(prefix, roc::calculate(rows, subset_size, conf_size)?)
}

pub(crate) fn write_roc_files_with_options(
    prefix: &Path,
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<MetricIndices> {
    publish(
        prefix,
        roc::calculate_with_options(rows, subset_size, conf_size, options)?,
    )
}

fn publish(prefix: &Path, artifacts: roc::Artifacts) -> Result<MetricIndices> {
    for artifact in artifacts.csv {
        let path = suffixed_report_path(prefix, &artifact.suffix);
        if artifact.optional && artifact.rows.is_empty() {
            if path.exists() {
                std::fs::remove_file(&path)
                    .with_context(|| format!("failed to remove stale {}", path.display()))?;
            }
            continue;
        }
        write_gzip_csv(&path, &artifact.header, &artifact.rows)?;
    }
    if let Some(table) = artifacts.raw_table {
        let path = suffixed_report_path(prefix, "roc.tsv");
        std::fs::write(&path, table)
            .with_context(|| format!("failed to write {}", path.display()))?;
    }
    Ok(artifacts.indices)
}

fn write_gzip_csv(path: &Path, header: &str, rows: &[String]) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = GzEncoder::new(file, Compression::default());
    writeln!(writer, "{header}")?;
    for row in rows {
        writeln!(writer, "{row}")?;
    }
    writer.finish()?;
    Ok(())
}
