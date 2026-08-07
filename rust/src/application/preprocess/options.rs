//! Cohesive preprocessing responsibility.

use crate::adapters::vcf;
use crate::application::{PreprocessArgs, PreprocessGender};
use crate::domain::RawVcfRecord;
use anyhow::{Context, Result, bail};
use std::collections::BTreeSet;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

#[cfg(test)]
pub(super) fn sort_normalized_records(records: &mut [RawVcfRecord]) {
    let mut contig_ranks = std::collections::HashMap::new();
    let mut next_rank = 0usize;
    for record in records.iter() {
        contig_ranks.entry(record.chrom.clone()).or_insert_with(|| {
            let rank = next_rank;
            next_rank += 1;
            rank
        });
    }
    records.sort_by(|left, right| {
        contig_ranks[&left.chrom]
            .cmp(&contig_ranks[&right.chrom])
            .then(left.pos.cmp(&right.pos))
    });
}

pub(super) fn effective_thread_count(threads: Option<usize>) -> usize {
    let available_threads = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1);
    effective_thread_count_with_available(threads, available_threads)
}

pub(super) fn effective_thread_count_with_available(
    threads: Option<usize>,
    available_threads: usize,
) -> usize {
    threads.unwrap_or_else(|| available_threads.max(1))
}

pub(super) fn has_non_reference_genotype(record: &RawVcfRecord) -> bool {
    let Some(format) = record.format.as_deref() else {
        return false;
    };
    let Some(gt_index) = format.split(':').position(|field| field == "GT") else {
        return false;
    };
    record.samples.iter().any(|sample| {
        sample
            .split(':')
            .nth(gt_index)
            .unwrap_or(".")
            .split(['/', '|'])
            .any(|allele| allele.parse::<usize>().is_ok_and(|allele| allele > 0))
    })
}

pub(super) struct PreprocessLogger {
    verbose: bool,
    quiet: bool,
    file: Option<File>,
}

impl PreprocessLogger {
    pub(super) fn new(args: &PreprocessArgs) -> Result<Self> {
        let file = args
            .logfile
            .as_deref()
            .map(|path| {
                if let Some(parent) = Path::new(path)
                    .parent()
                    .filter(|parent| !parent.as_os_str().is_empty())
                {
                    fs::create_dir_all(parent)?;
                }
                File::create(path).with_context(|| format!("failed to create logfile {path}"))
            })
            .transpose()?;
        Ok(Self {
            verbose: args.verbose,
            quiet: args.quiet,
            file,
        })
    }

    pub(super) fn info(&mut self, message: &str) -> Result<()> {
        if !self.verbose || self.quiet {
            return Ok(());
        }
        if let Some(file) = &mut self.file {
            writeln!(file, "INFO {message}")?;
            file.flush()?;
        } else {
            eprintln!("[I] {message}");
        }
        Ok(())
    }
}

/// The Perl prefixing stage in legacy `pre.py` changes record CHROM values
/// without rewriting the input declarations. The following `bcftools view`
/// pass notices each newly used sequence and appends a length-less contig
/// declaration. Mirror that repair while leaving existing declarations and
/// their order untouched.
pub(super) fn ensure_emitted_contig_headers(headers: &mut Vec<String>, emitted_contigs: &[String]) {
    let mut declared: BTreeSet<String> = headers
        .iter()
        .filter_map(|line| {
            line.strip_prefix("##contig=<ID=")
                .and_then(|body| body.split([',', '>']).next())
                .map(str::to_string)
        })
        .collect();
    let mut additions = Vec::new();
    for chrom in emitted_contigs {
        if declared.insert(chrom.clone()) {
            additions.push(format!("##contig=<ID={chrom}>"));
        }
    }
    let insert_at = headers
        .iter()
        .position(|line| line.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    headers.splice(insert_at..insert_at, additions);
}

pub(super) fn ensure_pass_filter_header(headers: &mut Vec<String>) {
    if headers
        .iter()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"))
    {
        return;
    }
    let index = headers
        .iter()
        .position(|line| line.starts_with("##fileformat="))
        .map_or(0, |index| index + 1);
    headers.insert(
        index,
        "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
    );
}

pub(super) fn resolve_reference(explicit: Option<&str>) -> Result<PathBuf> {
    let hg19 = std::env::var_os("HG19").map(PathBuf::from);
    let hgref = std::env::var_os("HGREF").map(PathBuf::from);
    resolve_reference_candidates(
        explicit.map(Path::new),
        hg19.as_deref(),
        hgref.as_deref(),
        Path::new("/opt/hap.py-data/hg19.fa"),
    )
}

pub(super) fn require_output_parent(output: &Path) -> Result<()> {
    if let Some(parent) = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        && !parent.is_dir()
    {
        bail!("output parent does not exist: {}", parent.display());
    }
    Ok(())
}

pub(super) fn require_vcf_sample(headers: &[String]) -> Result<()> {
    let has_sample = headers
        .iter()
        .find(|line| line.starts_with("#CHROM"))
        .is_some_and(|line| {
            line.split('\t')
                .nth(9)
                .is_some_and(|sample| !sample.is_empty())
        });
    if !has_sample {
        bail!("input VCF has no samples");
    }
    Ok(())
}

pub(super) fn resolve_reference_candidates(
    explicit: Option<&Path>,
    hg19: Option<&Path>,
    hgref: Option<&Path>,
    fallback: &Path,
) -> Result<PathBuf> {
    if let Some(path) = explicit {
        return Ok(path.to_path_buf());
    }
    [hg19, hgref, Some(fallback)]
        .into_iter()
        .flatten()
        .find(|path| path.is_file())
        .map(Path::to_path_buf)
        .context("no reference file found; pass --reference or set HG19/HGREF")
}

pub(super) fn has_chr_prefix(contigs: &BTreeSet<String>) -> Option<bool> {
    let plain = (0..23)
        .map(|value| value.to_string())
        .chain(["X".into(), "Y".into(), "MT".into()])
        .filter(|name| contigs.contains(name))
        .count();
    let prefixed = (0..23)
        .map(|value| format!("chr{value}"))
        .chain(["chrX".into(), "chrY".into(), "chrM".into()])
        .filter(|name| contigs.contains(name))
        .count();
    match prefixed.cmp(&plain) {
        std::cmp::Ordering::Greater => Some(true),
        std::cmp::Ordering::Less => Some(false),
        std::cmp::Ordering::Equal => None,
    }
}

pub(super) fn resolve_fixchr(
    requested: Option<bool>,
    reference_contigs: &BTreeSet<String>,
    input_contigs: &BTreeSet<String>,
) -> bool {
    requested.unwrap_or_else(|| {
        has_chr_prefix(reference_contigs) == Some(true)
            && has_chr_prefix(input_contigs) == Some(false)
    })
}

pub(super) fn add_legacy_chr_prefix(chrom: &str) -> String {
    if chrom == "chrMT" {
        return "chrM".to_string();
    }
    if chrom.starts_with("chr") {
        return chrom.to_string();
    }
    let first = chrom.as_bytes().first().copied();
    if !matches!(first, Some(b'0'..=b'9' | b'X' | b'Y' | b'M')) {
        return chrom.to_string();
    }
    if chrom == "MT" || chrom == "M" {
        "chrM".to_string()
    } else {
        format!("chr{chrom}")
    }
}

pub(super) fn passes_filters_only(filter: &str, filters_only: Option<&str>) -> bool {
    let Some(filters_only) = filters_only.filter(|value| !value.is_empty()) else {
        return true;
    };
    if filter.is_empty() || matches!(filter, "." | "PASS") {
        return true;
    }
    let excluded: BTreeSet<&str> = filters_only.split(',').collect();
    filter.split(';').any(|name| !excluded.contains(name))
}

#[cfg(test)]
pub(super) fn resolve_gender(
    requested: PreprocessGender,
    records: &[RawVcfRecord],
) -> PreprocessGender {
    if requested != PreprocessGender::Auto {
        return requested;
    }
    let mut haploid_x = false;
    let mut diploid_x = false;
    for record in records {
        observe_gender(record, &mut haploid_x, &mut diploid_x);
    }
    if haploid_x && !diploid_x {
        PreprocessGender::Male
    } else {
        PreprocessGender::Female
    }
}

pub(super) fn observe_gender(record: &RawVcfRecord, haploid_x: &mut bool, diploid_x: &mut bool) {
    if !matches!(record.chrom.as_str(), "X" | "chrX" | "chrx") {
        return;
    }
    let Some(format) = record.format.as_deref() else {
        return;
    };
    let Some(gt_index) = format.split(':').position(|field| field == "GT") else {
        return;
    };
    for sample in &record.samples {
        let gt = sample.split(':').nth(gt_index).unwrap_or(".");
        // vcfcheck classifies ploidy from the encoded GT vector length,
        // including missing slots. Thus `./1` has ngt == 2 and unequal
        // alleles (-1 and 1), so it is diploid rather than haploid.
        let alleles: Vec<&str> = gt.split(['/', '|']).collect();
        if alleles.len() == 1 {
            *haploid_x = true;
        } else if alleles.len() > 2 || (alleles.len() == 2 && alleles[0] != alleles[1]) {
            *diploid_x = true;
        }
    }
}

pub(crate) fn infer_gender(path: &Path) -> Result<PreprocessGender> {
    let records = vcf::open_validated_vcf(path)?;
    let mut haploid_x = false;
    let mut diploid_x = false;
    for record in records {
        let record = record?;
        observe_gender(record.raw(), &mut haploid_x, &mut diploid_x);
    }
    Ok(if haploid_x && !diploid_x {
        PreprocessGender::Male
    } else {
        PreprocessGender::Female
    })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum SymbolicDeletionMaterialization {
    LeadingAnchor,
    ContigStart,
}
