//! Cohesive quantify stratification responsibility.

use super::annotations::{
    effective_reference_range, has_region, info_value, merge_region_tags, move_region_to_front,
    propagate_ga4gh_superlocus_for_samples, remove_region_tag, set_format_value,
};
use super::{BenchmarkSamples, LoadedRegions, RegionLevels, RegionMap};
use crate::adapters::vcf::{self, ValidatedVcf};
use crate::application::QuantifyArgs;
use crate::domain::{Interval, RawVcfRecord};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

pub(super) fn load_regions(
    args: &QuantifyArgs,
    reference_contigs: &BTreeSet<String>,
) -> Result<LoadedRegions> {
    let mut confidence = args
        .fp_bedfile
        .as_deref()
        .map(|path| load_region_bed(Path::new(path), reference_contigs, args.strat_fixchr))
        .transpose()
        .with_context(|| "failed to load qfy confidence regions")?;
    if let Some(truth_vcf) = args.adjust_conf_regions.as_deref() {
        let raw_confidence = confidence.as_deref().ok_or_else(|| {
            anyhow::anyhow!("qfy --adjust-conf-regions requires --false-positives")
        })?;
        let (_, truth_records) = vcf::load_raw_vcf(Path::new(truth_vcf)).with_context(|| {
            format!("failed to load --adjust-conf-regions truth VCF {truth_vcf}")
        })?;
        let padding = truth_confidence_padding(&truth_records, raw_confidence);
        confidence.get_or_insert_default().extend(padding);
    }
    let mut paths = BTreeMap::<String, PathBuf>::new();

    if let Some(tsv_path) = args.strat_tsv.as_deref() {
        let tsv_path = Path::new(tsv_path);
        let text = fs::read_to_string(tsv_path)
            .with_context(|| format!("failed to read stratification TSV {}", tsv_path.display()))?;
        for (line_index, line) in text.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let (name, raw_path) = line.split_once('\t').ok_or_else(|| {
                anyhow::anyhow!(
                    "stratification TSV line {} has no region file",
                    line_index + 1
                )
            })?;
            let path = resolve_stratification_path(raw_path.trim(), tsv_path);
            insert_region_path(&mut paths, name.trim(), path)?;
        }
    }

    for spec in &args.strat_regions {
        let (name, path) = spec.split_once(':').ok_or_else(|| {
            anyhow::anyhow!("invalid --stratification-region '{spec}'; expected NAME:BED")
        })?;
        insert_region_path(&mut paths, name.trim(), PathBuf::from(path.trim()))?;
    }

    let mut regions = BTreeMap::<String, Vec<Interval>>::new();
    let mut levels = BTreeMap::new();
    for (raw_name, path) in paths {
        let (name, fixed_label) = dynamic_region_label(&raw_name);
        // QuantifyRegions::load collapses every lane whose label begins
        // with CONF into the reserved confidence lane. hap.py relies on
        // that behavior for its generated CONF_VARS truth padding.
        if name.starts_with("CONF") {
            let intervals = load_region_bed(&path, reference_contigs, args.strat_fixchr)
                .with_context(|| format!("failed to load stratification region {name}"))?;
            confidence.get_or_insert_default().extend(intervals);
        } else {
            let (loaded, loaded_levels) = load_stratification_bed(
                &path,
                &name,
                fixed_label,
                reference_contigs,
                args.strat_fixchr,
            )
            .with_context(|| format!("failed to load stratification region {name}"))?;
            for (label, intervals) in loaded {
                regions.entry(label).or_default().extend(intervals);
            }
            levels.extend(loaded_levels);
        }
    }
    Ok((confidence, regions, levels))
}

pub(super) fn dynamic_region_label(raw_name: &str) -> (String, bool) {
    if let Some(name) = raw_name.strip_prefix('=') {
        (name.to_string(), true)
    } else {
        (raw_name.to_string(), raw_name.starts_with("CONF"))
    }
}

pub(super) fn truth_confidence_padding(
    records: &[RawVcfRecord],
    targets: &[Interval],
) -> Vec<Interval> {
    let mut records = records.iter().collect::<Vec<_>>();
    records.sort_by(|left, right| left.chrom.cmp(&right.chrom).then(left.pos.cmp(&right.pos)));
    let mut output = Vec::new();
    let mut active: Option<Interval> = None;
    for record in records {
        if !targets
            .iter()
            .any(|target| target.matches(&record.chrom, record.pos))
        {
            continue;
        }
        let (start, end) = effective_reference_range(record)
            .map(|(start, end, _)| (start.saturating_sub(1), end))
            .unwrap_or((
                record.pos.saturating_sub(1),
                record.pos.saturating_sub(1) + record.ref_allele.len(),
            ));
        let overlaps_active = active.as_ref().is_some_and(|interval| {
            interval.chrom == record.chrom && start < interval.end && end > interval.start
        });
        if overlaps_active {
            if let Some(interval) = active.as_mut() {
                interval.start = interval.start.min(start);
                // gvcf2bed intentionally extends with the next record's start,
                // not its end; qfy relies on this historical under-extension.
                interval.end = interval.end.max(start + 1);
            }
        } else {
            if let Some(interval) = active.take() {
                output.push(interval);
            }
            active = Some(Interval {
                chrom: record.chrom.clone(),
                start,
                end,
            });
        }
    }
    if let Some(interval) = active {
        output.push(interval);
    }
    output
}

pub(super) fn resolve_stratification_path(raw_path: &str, tsv_path: &Path) -> PathBuf {
    let path = PathBuf::from(raw_path);
    if path.exists() || path.is_absolute() {
        path
    } else {
        tsv_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(path)
    }
}

pub(super) fn insert_region_path(
    paths: &mut BTreeMap<String, PathBuf>,
    name: &str,
    path: PathBuf,
) -> Result<()> {
    if name.is_empty() || name == "=" {
        bail!("stratification region name cannot be empty");
    }
    if name == "CONF" {
        bail!("stratification region name CONF is reserved for --false-positives");
    }
    if path.as_os_str().is_empty() {
        bail!("no file for stratification region {name}");
    }
    if paths.insert(name.to_string(), path).is_some() {
        bail!("duplicate stratification region ID: {name}");
    }
    Ok(())
}

pub(super) fn load_region_bed(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    fixchr: bool,
) -> Result<Vec<Interval>> {
    Ok(load_region_bed_rows(path, reference_contigs, fixchr)?
        .into_iter()
        .map(|(interval, _)| interval)
        .collect())
}

pub(super) fn load_stratification_bed(
    path: &Path,
    parent_label: &str,
    fixed_label: bool,
    reference_contigs: &BTreeSet<String>,
    fixchr: bool,
) -> Result<(RegionMap, RegionLevels)> {
    let mut regions = RegionMap::new();
    let mut levels = RegionLevels::from([(parent_label.to_string(), 0)]);
    for (interval, dynamic_label) in load_region_bed_rows(path, reference_contigs, fixchr)? {
        regions
            .entry(parent_label.to_string())
            .or_default()
            .push(interval.clone());
        if !fixed_label && let Some(dynamic_label) = dynamic_label {
            // Pinned QuantifyRegions adds the complete fourth-column value
            // as one level-1 child. Its apparent trailing-number hierarchy
            // branch starts at string[size] and is therefore unreachable:
            // `coding_1` is EXTRA_coding_1, never EXTRA_coding + child.
            let label = format!("{parent_label}_{dynamic_label}");
            regions.entry(label.clone()).or_default().push(interval);
            levels.insert(label, 1);
        }
    }
    // Empty BEDs still register their parent lane in legacy QuantifyRegions.
    regions.entry(parent_label.to_string()).or_default();
    Ok((regions, levels))
}

pub(super) fn load_region_bed_rows(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    fixchr: bool,
) -> Result<Vec<(Interval, Option<String>)>> {
    let text =
        vcf::read_text(path).with_context(|| format!("failed to read BED {}", path.display()))?;
    let mut intervals = Vec::new();
    for (line_index, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() < 3 {
            bail!(
                "BED line {} has fewer than 3 columns in {}",
                line_index + 1,
                path.display()
            );
        }
        let chrom = if fixchr {
            vcf::normalize_chrom(fields[0], reference_contigs)
        } else {
            fields[0].to_string()
        };
        let start = fields[1]
            .parse::<usize>()
            .with_context(|| format!("invalid BED start '{}' in {}", fields[1], path.display()))?;
        let end = fields[2]
            .parse::<usize>()
            .with_context(|| format!("invalid BED end '{}' in {}", fields[2], path.display()))?;
        if end < start {
            bail!(
                "BED end precedes start on line {} in {}",
                line_index + 1,
                path.display()
            );
        }
        intervals.push((
            Interval { chrom, start, end },
            fields.get(3).map(|label| (*label).to_string()),
        ));
    }
    Ok(intervals)
}

#[cfg(test)]
pub(super) fn annotate_regions(
    record: &mut RawVcfRecord,
    confidence: Option<&[Interval]>,
    stratifications: &RegionMap,
    rewrite_ga4gh_decisions: bool,
) {
    annotate_regions_for_samples(
        record,
        confidence,
        stratifications,
        rewrite_ga4gh_decisions,
        BenchmarkSamples::POSITIONAL,
        false,
    );
}

pub(super) fn annotate_regions_for_samples(
    record: &mut RawVcfRecord,
    confidence: Option<&[Interval]>,
    stratifications: &RegionMap,
    rewrite_ga4gh_decisions: bool,
    samples: BenchmarkSamples,
    preserve_missing_nocall_bd: bool,
) {
    let effective_range = effective_reference_range(record);
    let record_chrom = record.chrom.clone();
    let record_pos = record.pos;
    let record_end = record.end_pos();
    let overlaps = |intervals: &[Interval]| match effective_range {
        Some((start, end, pure_insertion)) => {
            if pure_insertion {
                // Insertions occupy the gap between two reference anchors.
                // Legacy only marks them confident when both anchors are
                // covered (possibly by different CONF input lanes).
                (start..=end).all(|pos| {
                    intervals
                        .iter()
                        .any(|interval| interval.matches(&record_chrom, pos))
                })
            } else {
                intervals
                    .iter()
                    .any(|interval| interval.overlaps(&record_chrom, start, end))
            }
        }
        None => intervals
            .iter()
            .any(|interval| interval.overlaps(&record_chrom, record_pos, record_end)),
    };
    let fully_covered = |intervals: &[Interval]| {
        let (start, end) = effective_range
            .map(|(start, end, _)| (start, end))
            .unwrap_or((record_pos, record_end));
        (start..=end).all(|pos| {
            intervals
                .iter()
                .any(|interval| interval.matches(&record_chrom, pos))
        })
    };
    let mut additions = Vec::new();
    let mut crosses_stratification_boundary = false;
    for (name, intervals) in stratifications {
        if overlaps(intervals) {
            additions.push(name.clone());
            crosses_stratification_boundary |= !fully_covered(intervals);
        }
    }
    if crosses_stratification_boundary {
        additions.push("TS_boundary".to_string());
    }

    if let Some(confidence) = confidence {
        remove_region_tag(&mut record.info, "CONF");
        remove_region_tag(&mut record.info, "TS_boundary");
        if overlaps(confidence) {
            additions.insert(0, "CONF".to_string());
        } else if rewrite_ga4gh_decisions {
            // GA4GHQuantify's `count_unk` rule rewrites both samples outside
            // CONF. Truth-side UNK is omitted from truth counts; query-side
            // UNK supplies QUERY.UNK and ROC unknown counts.
            if let Some(truth) = samples.truth
                && !preserve_missing_nocall_decision(record, truth, preserve_missing_nocall_bd)
            {
                replace_existing_decision(record, truth, "UNK");
            }
            if let Some(query) = samples.query
                && !preserve_missing_nocall_decision(record, query, preserve_missing_nocall_bd)
            {
                replace_existing_decision(record, query, "UNK");
            }
        }
    }
    if !additions.is_empty() {
        merge_region_tags(&mut record.info, &additions);
    }
    if has_region(&record.info, "CONF") {
        move_region_to_front(&mut record.info, "CONF");
    }
}

/// `BlockQuantify::count` finalizes region flags per benchmarking superlocus.
/// A block crossing the recomputed CONF boundary marks every member as
/// `TS_boundary`; a wholly confident block marks every member `TS_contained`.
/// Records without a non-negative BS value form singleton blocks.
#[cfg(test)]
pub(super) fn propagate_superlocus_annotations(
    records: &mut [RawVcfRecord],
    annotation_type: &str,
) {
    propagate_superlocus_annotations_for_samples(
        records,
        annotation_type,
        BenchmarkSamples::POSITIONAL,
        false,
        false,
    );
}

pub(super) fn propagate_superlocus_annotations_for_samples(
    records: &mut [RawVcfRecord],
    annotation_type: &str,
    samples: BenchmarkSamples,
    preserve_missing_query_qq: bool,
    inherit_same_position_tp_qq: bool,
) {
    let mut start = 0usize;
    while start < records.len() {
        let chrom = records[start].chrom.clone();
        let bs = benchmark_superlocus(&records[start].info);
        let mut end = start + 1;
        if bs.is_some() {
            while end < records.len()
                && records[end].chrom == chrom
                && benchmark_superlocus(&records[end].info) == bs
            {
                end += 1;
            }
        }

        let has_confident = records[start..end]
            .iter()
            .any(|record| has_region(&record.info, "CONF"));
        let has_non_confident = records[start..end].iter().any(|record| {
            !has_region(&record.info, "CONF") || has_region(&record.info, "TS_boundary")
        });
        let region = if has_confident && has_non_confident {
            Some("TS_boundary")
        } else if has_confident {
            Some("TS_contained")
        } else {
            None
        };
        if let Some(region) = region {
            for record in &mut records[start..end] {
                merge_region_tags(&mut record.info, &[region.to_string()]);
            }
        }

        // XCMP has no truth TP in the compatibility lane, so the legacy
        // post-pass clears truth QQ rather than retaining record QUAL.
        if annotation_type == "xcmp" {
            if let Some(truth) = samples.truth {
                for record in &mut records[start..end] {
                    if record.sample_map(truth).get("BD").map(String::as_str) != Some("TP") {
                        set_format_value(record, truth, "QQ", ".");
                    }
                }
            }
        } else {
            propagate_ga4gh_superlocus_for_samples(
                &mut records[start..end],
                samples,
                preserve_missing_query_qq,
                inherit_same_position_tp_qq,
            );
        }
        start = end;
    }
}

pub(super) fn propagate_checked_superlocus_annotations_for_samples(
    input: &mut ValidatedVcf,
    annotation_type: &str,
    samples: BenchmarkSamples,
    preserve_missing_query_qq: bool,
    inherit_same_position_tp_qq: bool,
) -> Result<()> {
    let mut edited = input.records().to_vec();
    propagate_superlocus_annotations_for_samples(
        &mut edited,
        annotation_type,
        samples,
        preserve_missing_query_qq,
        inherit_same_position_tp_qq,
    );
    for (index, replacement) in edited.into_iter().enumerate() {
        input.try_edit_record(index, |record| {
            *record = replacement;
            Ok(())
        })?;
    }
    Ok(())
}

pub(super) fn benchmark_superlocus(info: &str) -> Option<i64> {
    info_value(info, "BS")?
        .split(',')
        .next()?
        .parse::<i64>()
        .ok()
        .filter(|value| *value >= 0)
}

pub(super) fn replace_existing_decision(
    record: &mut RawVcfRecord,
    sample_index: usize,
    decision: &str,
) {
    if record.sample_map(sample_index).contains_key("BD") {
        set_format_value(record, sample_index, "BD", decision);
    }
}

pub(super) fn preserve_missing_nocall_decision(
    record: &RawVcfRecord,
    sample_index: usize,
    enabled: bool,
) -> bool {
    if !enabled {
        return false;
    }
    let sample = record.sample_map(sample_index);
    sample.get("BVT").map(String::as_str) == Some("NOCALL")
        && sample
            .get("BD")
            .is_none_or(|decision| decision.is_empty() || decision == ".")
}
