//! Bounded-memory transformation of quantified records.

use super::annotations::{
    decorate_quantified_record_for_samples, ga4gh_qq_is_string, legacy_regions_extent,
    normalize_integer_like_format_values, reannotate_ga4gh_record,
    reannotate_xcmp_record_for_samples, set_info_value, validate_ga4gh_qq_record,
};
use super::stratification::{
    annotate_regions_for_samples, benchmark_superlocus,
    propagate_checked_superlocus_annotations_for_samples,
};
use super::{BenchmarkSamples, CompareQuantifyMode, RegionMap};
use crate::adapters::vcf::{ValidatedVcf, ValidatedVcfRecord};
use crate::application::QuantifyArgs;
use anyhow::{Context, Result, bail};
use std::collections::BTreeSet;
use std::io::{BufWriter, Write};

const MAX_QUANTIFY_SUPERLOCUS_RECORDS: usize = 1_000_000;

/// Annotation settings shared by the quantify streaming pass: which annotation
/// dialect to write, the mode toggles that gate legacy quirks, and which
/// samples carry benchmark decisions.
#[derive(Clone, Copy)]
pub(super) struct QuantifyAnnotation<'a> {
    pub(super) annotation_type: &'a str,
    pub(super) mode: CompareQuantifyMode,
    pub(super) benchmark_samples: BenchmarkSamples,
}

pub(super) fn spool_quantified_records<I>(
    records: I,
    headers: &[String],
    annotation: QuantifyAnnotation<'_>,
    args: &QuantifyArgs,
    confidence: Option<&[crate::domain::Interval]>,
    stratifications: &RegionMap,
) -> Result<(tempfile::NamedTempFile, BTreeSet<String>)>
where
    I: IntoIterator<Item = Result<ValidatedVcfRecord>>,
{
    let QuantifyAnnotation {
        annotation_type,
        mode,
        benchmark_samples,
    } = annotation;
    let mut spool = tempfile::NamedTempFile::new().context("failed to create quantify spool")?;
    let qq_is_string = annotation_type == "ga4gh" && ga4gh_qq_is_string(headers);
    let mut input_contigs = BTreeSet::new();
    {
        let mut writer = BufWriter::new(spool.as_file_mut());
        for header in headers {
            writeln!(writer, "{header}")?;
        }
        let mut group = Vec::new();
        let mut group_key: Option<(String, i64)> = None;
        for record in records {
            let mut record = record?;
            if qq_is_string {
                validate_ga4gh_qq_record(&record)?;
            }
            input_contigs.insert(record.chrom.clone());
            record.try_update(|raw| {
                if args.preserve_info {
                    let extent = legacy_regions_extent(raw);
                    set_info_value(&mut raw.info, "RegionsExtent", &extent);
                }
                annotate_regions_for_samples(
                    raw,
                    confidence,
                    stratifications,
                    annotation_type == "ga4gh",
                    benchmark_samples,
                    mode.preserve_missing_nocall_bd,
                );
                if annotation_type == "xcmp" {
                    reannotate_xcmp_record_for_samples(
                        raw,
                        confidence.is_some(),
                        &args.roc,
                        benchmark_samples,
                    );
                } else {
                    reannotate_ga4gh_record(raw);
                }
                Ok(())
            })?;

            let key = benchmark_superlocus(&record.info)
                .map(|superlocus| (record.chrom.clone(), superlocus));
            if !group.is_empty() && (key.is_none() || key != group_key) {
                flush_quantify_group(
                    &mut writer,
                    &mut group,
                    annotation,
                    args,
                    confidence.is_some(),
                )?;
            }
            group_key = key;
            group.push(record);
            if group.len() > MAX_QUANTIFY_SUPERLOCUS_RECORDS {
                bail!(
                    "quantify superlocus in {} exceeds the {} record active-window limit",
                    args.input_vcf,
                    MAX_QUANTIFY_SUPERLOCUS_RECORDS
                );
            }
            if group_key.is_none() {
                flush_quantify_group(
                    &mut writer,
                    &mut group,
                    annotation,
                    args,
                    confidence.is_some(),
                )?;
            }
        }
        flush_quantify_group(
            &mut writer,
            &mut group,
            annotation,
            args,
            confidence.is_some(),
        )?;
        writer.flush()?;
    }
    spool.as_file().sync_all()?;
    Ok((spool, input_contigs))
}

fn flush_quantify_group(
    writer: &mut dyn Write,
    group: &mut Vec<ValidatedVcfRecord>,
    annotation: QuantifyAnnotation<'_>,
    args: &QuantifyArgs,
    has_confidence: bool,
) -> Result<()> {
    let QuantifyAnnotation {
        annotation_type,
        mode,
        benchmark_samples,
    } = annotation;
    if group.is_empty() {
        return Ok(());
    }
    let mut checked = ValidatedVcf::from_parts(Vec::new(), std::mem::take(group));
    propagate_checked_superlocus_annotations_for_samples(
        &mut checked,
        annotation_type,
        benchmark_samples,
        mode.preserve_missing_query_qq,
        mode.inherit_same_position_tp_qq,
    )?;
    let (_, records) = checked.into_parts();
    for mut record in records {
        record.try_update(|raw| {
            decorate_quantified_record_for_samples(
                raw,
                annotation_type,
                args.preserve_info,
                args.output_vtc,
                has_confidence,
                benchmark_samples,
            );
            if annotation_type == "ga4gh" {
                normalize_integer_like_format_values(raw, "QQ");
            }
            Ok(())
        })?;
        writeln!(writer, "{}", record.raw().to_line())?;
    }
    Ok(())
}
