use anyhow::{Context, Result, bail};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

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
