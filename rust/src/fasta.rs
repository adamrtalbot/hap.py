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
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 2 || fields[0].is_empty() {
            bail!(
                "invalid FASTA index line {} in {}",
                line_number + 1,
                index_path.display()
            );
        }
        let length = fields[1].parse::<usize>().with_context(|| {
            format!(
                "invalid FASTA index length '{}' on line {} in {}",
                fields[1],
                line_number + 1,
                index_path.display()
            )
        })?;
        contigs.insert(fields[0].to_string(), length);
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
                contigs.insert(name, current_seq.clone());
            }
            let name = rest
                .split_whitespace()
                .next()
                .ok_or_else(|| anyhow::anyhow!("invalid FASTA header in {}", path.display()))?;
            current_name = Some(name.to_string());
            current_seq.clear();
        } else {
            current_seq.push_str(line.trim());
        }
    }

    if let Some(name) = current_name {
        contigs.insert(name, current_seq);
    }

    if contigs.is_empty() {
        bail!("no contigs found in FASTA {}", path.display());
    }

    Ok(contigs)
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
}
