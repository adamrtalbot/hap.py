//! Bounded-memory ROC calculation for germline reports.
//!
//! The engine writes only to an engine-owned temporary directory. Publication
//! of those artifacts to the caller's staged output generation remains an
//! application-adapter responsibility.

use crate::domain::AnnotatedRow;
use anyhow::{Context, Result};
use std::borrow::Borrow;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

mod bounded;

pub(crate) use bounded::{MetricIndices, RocOptions};

/// A calculated CSV artifact held in an engine-owned temporary file.
#[derive(Debug)]
pub(crate) struct CsvArtifact {
    pub(crate) suffix: String,
    pub(crate) source: Option<tempfile::TempPath>,
    pub(crate) optional: bool,
}

/// Filesystem-neutral destination plan returned by the ROC engine.
#[derive(Debug)]
pub(crate) struct Artifacts {
    pub(crate) indices: MetricIndices,
    pub(crate) csv: Vec<CsvArtifact>,
    pub(crate) raw_table: Option<tempfile::TempPath>,
}

/// Calculate ROC reports from a fallible record stream without retaining the
/// complete input or output tables in memory.
pub(crate) fn calculate_with_options_iter<I, R>(
    rows: I,
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<Artifacts>
where
    I: IntoIterator<Item = Result<R>>,
    R: Borrow<AnnotatedRow>,
{
    let directory = tempfile::tempdir().context("failed to create ROC engine workspace")?;
    let prefix = directory.path().join("report");
    let indices =
        bounded::write_roc_files_with_options_iter(&prefix, rows, subset_size, conf_size, options)?;

    let mut candidates = vec![
        ("roc.all.csv.gz".to_string(), false),
        ("roc.Locations.SNP.csv.gz".to_string(), true),
        ("roc.Locations.SNP.PASS.csv.gz".to_string(), true),
        ("roc.Locations.INDEL.csv.gz".to_string(), true),
        ("roc.Locations.INDEL.PASS.csv.gz".to_string(), true),
    ];
    let mut seen = candidates
        .iter()
        .map(|(suffix, _)| suffix.clone())
        .collect::<BTreeSet<_>>();
    for ty in ["SNP", "INDEL"] {
        let suffix = format!("roc.Locations.{ty}.SEL.csv.gz");
        if seen.insert(suffix.clone()) {
            candidates.push((suffix, true));
        }
    }

    let mut csv = Vec::with_capacity(candidates.len());
    for (suffix, optional) in candidates {
        let path = suffixed_path(&prefix, &suffix);
        let source = path
            .exists()
            .then(|| copy_to_temp(&path, "ROC CSV"))
            .transpose()?;
        if !optional && source.is_none() {
            anyhow::bail!("ROC engine did not produce required artifact {suffix}");
        }
        csv.push(CsvArtifact {
            suffix,
            source,
            optional,
        });
    }

    let raw_path = suffixed_path(&prefix, "roc.tsv");
    let raw_table = raw_path
        .exists()
        .then(|| copy_to_temp(&raw_path, "legacy ROC table"))
        .transpose()?;

    Ok(Artifacts {
        indices,
        csv,
        raw_table,
    })
}

fn suffixed_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(".");
    path.push(suffix);
    PathBuf::from(path)
}

fn copy_to_temp(source: &Path, description: &str) -> Result<tempfile::TempPath> {
    let input = File::open(source)
        .with_context(|| format!("failed to open {description} {}", source.display()))?;
    let mut output =
        tempfile::NamedTempFile::new().with_context(|| format!("failed to spool {description}"))?;
    {
        let mut reader = BufReader::new(input);
        let mut writer = BufWriter::new(output.as_file_mut());
        std::io::copy(&mut reader, &mut writer)
            .with_context(|| format!("failed to spool {description}"))?;
        writer.flush()?;
    }
    Ok(output.into_temp_path())
}
