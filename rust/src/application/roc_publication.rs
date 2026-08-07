//! Publication adapter for filesystem-neutral ROC engine artifacts.

use crate::adapters::report::suffixed_report_path;
use crate::domain::AnnotatedRow;
use crate::engines::roc::{self, MetricIndices, RocOptions};
use crate::output::{FailureOperation, fail_operation};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;

pub(crate) fn write_roc_files_with_options_iter<I, R>(
    prefix: &Path,
    rows: I,
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<MetricIndices>
where
    I: IntoIterator<Item = Result<R>>,
    R: std::borrow::Borrow<AnnotatedRow>,
{
    publish(
        prefix,
        roc::calculate_with_options_iter(rows, subset_size, conf_size, options)?,
    )
}

fn publish(prefix: &Path, artifacts: roc::Artifacts) -> Result<MetricIndices> {
    for artifact in artifacts.csv {
        let path = suffixed_report_path(prefix, &artifact.suffix);
        let Some(source) = artifact.source else {
            debug_assert!(artifact.optional);
            if path.exists() {
                std::fs::remove_file(&path)
                    .with_context(|| format!("failed to remove stale {}", path.display()))?;
            }
            continue;
        };
        publish_file(&source, &path)?;
    }
    if let Some(source) = artifacts.raw_table {
        let path = suffixed_report_path(prefix, "roc.tsv");
        publish_file(&source, &path)?;
    }
    Ok(artifacts.indices)
}

fn publish_file(source: &Path, path: &Path) -> Result<()> {
    fail_operation(FailureOperation::Writer, path)?;
    let input = File::open(source)
        .with_context(|| format!("failed to open ROC artifact {}", source.display()))?;
    let output =
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut reader = BufReader::new(input);
    let mut writer = BufWriter::new(output);
    std::io::copy(&mut reader, &mut writer)
        .with_context(|| format!("failed to publish ROC artifact {}", path.display()))?;
    fail_operation(FailureOperation::Encoder, path)?;
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::{FailureOperation, OutputTransaction, set_failure_operation};
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn injected_roc_operations_preserve_generation_and_cleanup() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let output = directory.path().join("report.roc.all.csv.gz");
        // The engine-owned source spool is outside the destination generation;
        // only the destination and its adjacent transaction staging file belong
        // in this directory.
        let mut source = tempfile::NamedTempFile::new()?;
        source.write_all(b"new-roc")?;
        source.flush()?;

        for operation in [FailureOperation::Writer, FailureOperation::Encoder] {
            fs::write(&output, "old-roc")?;
            let transaction = OutputTransaction::files(Vec::<PathBuf>::new(), [&output])?;
            let staged = transaction.staged_file(&output)?.to_path_buf();
            set_failure_operation(Some(operation));
            let result = publish_file(source.path(), &staged)
                .with_context(|| format!("failed to write ROC table {}", output.display()));
            set_failure_operation(None);
            let error = result.expect_err("injected ROC operation must fail");
            drop(transaction);

            assert!(error.to_string().contains(&output.display().to_string()));
            assert_eq!(fs::read_to_string(&output)?, "old-roc");
            assert_eq!(fs::read_dir(directory.path())?.count(), 1);
        }
        Ok(())
    }
}
