use anyhow::{Context, Result, bail};
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

pub fn read_index(path: &Path) -> Result<BTreeMap<String, usize>> {
    let index_path = index_path(path);
    let text = match fs::read_to_string(&index_path) {
        Ok(text) => text,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            bail!("Fasta file {} is not indexed", path.display())
        }
        Err(error) => {
            return Err(error)
                .with_context(|| format!("failed to read FASTA index {}", index_path.display()));
        }
    };
    let mut contigs = BTreeMap::new();
    for (line_number, line) in text.lines().enumerate() {
        let mut fields = line.split('\t');
        let name = fields.next().unwrap_or_default();
        let Some(raw_length) = fields.next() else {
            bail!(
                "invalid FASTA index line {} in {}",
                line_number + 1,
                index_path.display()
            );
        };
        if name.is_empty() {
            bail!(
                "invalid FASTA index line {} in {}",
                line_number + 1,
                index_path.display()
            );
        }
        let length = raw_length.parse::<usize>().with_context(|| {
            format!(
                "invalid FASTA index length '{}' on line {} in {}",
                raw_length,
                line_number + 1,
                index_path.display()
            )
        })?;
        contigs.insert(name.to_string(), length);
    }
    Ok(contigs)
}

fn index_path(path: &Path) -> PathBuf {
    let mut index = path.as_os_str().to_os_string();
    index.push(".fai");
    PathBuf::from(index)
}

pub fn read_sequences(path: &Path) -> Result<BTreeMap<String, String>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read FASTA {}", path.display()))?;
    let mut contigs = BTreeMap::new();
    let mut current_name: Option<String> = None;
    let mut current_seq = String::new();

    for line in text.lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = current_name.take() {
                insert_contig(&mut contigs, name, std::mem::take(&mut current_seq), path)?;
            }
            let name = rest
                .split_whitespace()
                .next()
                .ok_or_else(|| anyhow::anyhow!("invalid FASTA header in {}", path.display()))?;
            current_name = Some(name.to_string());
        } else {
            if current_name.is_none() && !line.trim().is_empty() {
                bail!(
                    "sequence data precedes the first FASTA header in {}",
                    path.display()
                );
            }
            current_seq.push_str(line.trim());
        }
    }

    if let Some(name) = current_name {
        insert_contig(&mut contigs, name, current_seq, path)?;
    }

    if contigs.is_empty() {
        bail!("no contigs found in FASTA {}", path.display());
    }

    Ok(contigs)
}

fn insert_contig(
    contigs: &mut BTreeMap<String, String>,
    name: String,
    sequence: String,
    path: &Path,
) -> Result<()> {
    if contigs.insert(name.clone(), sequence).is_some() {
        bail!("duplicate FASTA contig '{name}' in {}", path.display());
    }
    Ok(())
}

pub fn contig_lengths(path: &Path) -> Result<BTreeMap<String, usize>> {
    Ok(read_sequences(path)?
        .into_iter()
        .map(|(name, sequence)| (name, sequence.len()))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn read_index_rejects_malformed_names_and_lengths() -> Result<()> {
        let directory = tempdir()?;
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        let index = index_path(&reference);

        fs::write(&index, "chr1\n")?;
        assert!(
            read_index(&reference)
                .unwrap_err()
                .to_string()
                .contains("invalid FASTA index line")
        );

        fs::write(&index, "chr1\tnot-a-length\t6\t5\t6\n")?;
        assert!(
            read_index(&reference)
                .unwrap_err()
                .to_string()
                .contains("invalid FASTA index length")
        );
        Ok(())
    }

    #[test]
    fn read_index_returns_contig_lengths() -> Result<()> {
        let directory = tempdir()?;
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(index_path(&reference), "chr1\t5\t6\t5\t6\n")?;

        assert_eq!(
            read_index(&reference)?,
            BTreeMap::from([("chr1".to_string(), 5)])
        );
        Ok(())
    }

    #[test]
    fn read_index_matches_legacy_permissive_fai_fields() -> Result<()> {
        let directory = tempdir()?;
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        let index = index_path(&reference);

        fs::write(&index, "")?;
        assert_eq!(read_index(&reference)?, BTreeMap::new());

        fs::write(&index, "chr1\t5\n")?;
        assert_eq!(
            read_index(&reference)?,
            BTreeMap::from([("chr1".to_string(), 5)])
        );

        fs::write(&index, "chr1\t1\tinvalid\tgeometry\tis-ignored\nchr1\t5\n")?;
        assert_eq!(
            read_index(&reference)?,
            BTreeMap::from([("chr1".to_string(), 5)])
        );
        Ok(())
    }

    #[test]
    fn read_sequences_rejects_ambiguous_fasta_structure() -> Result<()> {
        let directory = tempdir()?;
        let reference = directory.path().join("ref.fa");

        fs::write(&reference, "ACGT\n>chr1\nACGT\n")?;
        assert!(
            read_sequences(&reference)
                .unwrap_err()
                .to_string()
                .contains("precedes the first FASTA header")
        );

        fs::write(&reference, ">chr1\nAC\n>chr1\nGT\n")?;
        assert!(
            read_sequences(&reference)
                .unwrap_err()
                .to_string()
                .contains("duplicate FASTA contig")
        );
        Ok(())
    }
}
