//! Filesystem boundary for the in-memory comparison engines.

use crate::adapters::{fasta, vcf};
use crate::engines::{scmp, vcfeval};
use anyhow::Result;
use std::path::Path;

pub(crate) fn run_scmp(
    truth: &Path,
    query: &Path,
    reference: &Path,
    mode: scmp::ScmpMode,
    qq_field: &str,
) -> Result<vcf::ValidatedVcf> {
    let (truth_headers, truth_records) = vcf::load_raw_vcf(truth)?;
    let (query_headers, query_records) = vcf::load_raw_vcf(query)?;
    let mut merged = scmp::merge_two_sample_records(
        &truth_headers,
        &truth_records,
        &query_headers,
        &query_records,
    )?;
    let references = fasta::read_sequences(reference)?;
    scmp::annotate_merged_records(
        &mut merged.records,
        &merged.headers,
        &references,
        mode,
        qq_field,
    )?;
    vcf::ValidatedVcf::try_from_raw(merged.headers, merged.records)
}

pub(crate) fn run_vcfeval(
    truth: &Path,
    query: &Path,
    reference: &Path,
    options: vcfeval::Options<'_>,
) -> Result<vcf::ValidatedVcf> {
    let (truth_headers, truth_records) = vcf::load_raw_vcf(truth)?;
    let (query_headers, query_records) = vcf::load_raw_vcf(query)?;
    let references = fasta::read_sequences(reference)?;
    vcfeval::compare_records(
        &truth_headers,
        &truth_records,
        &query_headers,
        &query_records,
        &references,
        options,
    )
}
