//! Cohesive quantify reporting responsibility.

use super::{ExtendedTableOptions, INDEL_SUBTYPES, QuantifyCountMaps, QuantifyTypeCounts};
use crate::adapters::metrics_json;
use crate::adapters::report::{self, suffixed_report_path};
use crate::domain::{CountsBucket, TypeCounts};
use crate::engines::roc;
use anyhow::{Context, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

pub(super) fn write_quantify_summary(
    path: &Path,
    all_counts: &BTreeMap<String, QuantifyTypeCounts>,
    pass_counts: &BTreeMap<String, QuantifyTypeCounts>,
) -> Result<()> {
    let mut writer = BufWriter::new(
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(
        writer,
        "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio"
    )?;
    if all_counts.is_empty() && pass_counts.is_empty() {
        for variant_type in ["INDEL", "SNP"] {
            writeln!(writer, "{variant_type},ALL,0,0,0,0,0,0,0,0,,,,,,,,")?;
        }
        return writer
            .flush()
            .with_context(|| format!("failed to flush {}", path.display()));
    }
    for variant_type in ["INDEL", "SNP"].into_iter().filter(|variant_type| {
        all_counts.contains_key(*variant_type) || pass_counts.contains_key(*variant_type)
    }) {
        let all = all_counts.get(variant_type).cloned().unwrap_or_default();
        let pass = pass_counts.get(variant_type).cloned().unwrap_or_default();
        write_summary_row(&mut writer, variant_type, "ALL", &all)?;
        write_summary_row(&mut writer, variant_type, "PASS", &pass)?;
    }
    writer
        .flush()
        .with_context(|| format!("failed to flush {}", path.display()))
}

pub(super) fn write_summary_row<W: Write>(
    writer: &mut W,
    variant_type: &str,
    filter: &str,
    stats: &QuantifyTypeCounts,
) -> Result<()> {
    let truth_titv = if variant_type == "SNP" {
        ti_tv_ratio(stats.truth_total.ti, stats.truth_total.tv)
    } else {
        String::new()
    };
    let query_titv = if variant_type == "SNP" {
        ti_tv_ratio(stats.query_total.ti, stats.query_total.tv)
    } else {
        String::new()
    };
    writeln!(
        writer,
        "{variant_type},{filter},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
        stats.truth_total.total,
        stats.truth_tp.total,
        stats.truth_fn.total,
        stats.query_total.total,
        stats.query_fp.total,
        stats.query_unk.total,
        stats.fp_gt,
        stats.fp_al,
        metric_ratio(stats.truth_tp.total, stats.truth_total.total),
        precision_metric(stats),
        metric_ratio(stats.query_unk.total, stats.query_total.total),
        f1_score(
            stats.truth_tp.total,
            stats.truth_total.total,
            stats.query_tp.total,
            stats.query_tp.total + stats.query_fp.total
        ),
        truth_titv,
        query_titv,
        het_hom_ratio(stats.truth_total.het, stats.truth_total.homalt),
        het_hom_ratio(stats.query_total.het, stats.query_total.homalt),
    )?;
    Ok(())
}

pub(super) fn write_quantify_extended(
    path: &Path,
    all_counts: &QuantifyCountMaps,
    pass_counts: &QuantifyCountMaps,
    subsets_present: &BTreeMap<String, BTreeSet<String>>,
    options: &ExtendedTableOptions<'_>,
) -> Result<()> {
    let mut writer = BufWriter::new(
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    let mut header = report::EXTENDED_HEADER
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>();
    if options.ci_alpha > 0.0 {
        header.extend(
            [
                "METRIC.Recall.Lower",
                "METRIC.Recall.Upper",
                "METRIC.Precision.Lower",
                "METRIC.Precision.Upper",
                "METRIC.Frac_NA.Lower",
                "METRIC.Frac_NA.Upper",
            ]
            .into_iter()
            .map(str::to_string),
        );
    }
    writeln!(writer, "{}", header.join(","))?;
    if all_counts.by_type.is_empty() && pass_counts.by_type.is_empty() {
        for line in report::empty_comparison_extended_lines(options.subset_size) {
            let mut row = line.split(',').map(str::to_string).collect::<Vec<_>>();
            row.resize(header.len(), String::new());
            writeln!(writer, "{}", row.join(","))?;
        }
        return writer
            .flush()
            .with_context(|| format!("failed to flush {}", path.display()));
    }
    let confidence_size_value = options.confidence_size;
    let confidence_size = confidence_size_value
        .map(|size| format!("{:.6}", size as f64))
        .unwrap_or_default();

    for variant_type in ["INDEL", "SNP"].into_iter().filter(|variant_type| {
        all_counts.by_type.contains_key(*variant_type)
            || pass_counts.by_type.contains_key(*variant_type)
    }) {
        let subtype_labels: Vec<String> = if variant_type == "INDEL" {
            INDEL_SUBTYPES
                .iter()
                .map(|value| value.to_string())
                .collect()
        } else {
            Vec::new()
        };
        let subsets = subsets_present
            .get(variant_type)
            .cloned()
            .unwrap_or_default()
            .into_iter()
            .collect::<Vec<_>>();

        write_extended_pair(
            &mut writer,
            variant_type,
            "*",
            "*",
            all_counts
                .by_type
                .get(variant_type)
                .cloned()
                .unwrap_or_default(),
            pass_counts
                .by_type
                .get(variant_type)
                .cloned()
                .unwrap_or_default(),
            options.subset_size,
            0,
            &confidence_size,
            options.qq_field,
            options.ci_alpha,
        )?;

        for subset in &subsets {
            let (subset_size_for_row, subset_confidence_size) = named_subset_report_sizes(
                subset,
                options.subset_size,
                options.whole_reference_size,
                confidence_size_value,
                options.stratification_sizes,
                options.stratification_confidence_sizes,
            );
            let subset_confidence_size = subset_confidence_size
                .map(|size| format!("{:.6}", size as f64))
                .unwrap_or_default();
            let all = all_counts
                .by_subset_type
                .get(subset)
                .and_then(|types| types.get(variant_type))
                .cloned()
                .unwrap_or_default();
            let pass = pass_counts
                .by_subset_type
                .get(subset)
                .and_then(|types| types.get(variant_type))
                .cloned()
                .unwrap_or_default();
            write_extended_pair(
                &mut writer,
                variant_type,
                "*",
                subset,
                all,
                pass,
                subset_size_for_row,
                options
                    .stratification_levels
                    .get(subset)
                    .copied()
                    .unwrap_or(0),
                &subset_confidence_size,
                options.qq_field,
                options.ci_alpha,
            )?;
        }

        if variant_type == "INDEL" {
            for subtype in &subtype_labels {
                let all = all_counts
                    .by_subtype
                    .get(variant_type)
                    .and_then(|subtypes| subtypes.get(subtype))
                    .cloned()
                    .unwrap_or_default();
                let pass = pass_counts
                    .by_subtype
                    .get(variant_type)
                    .and_then(|subtypes| subtypes.get(subtype))
                    .cloned()
                    .unwrap_or_default();
                write_extended_pair(
                    &mut writer,
                    variant_type,
                    subtype,
                    "*",
                    all,
                    pass,
                    options.subset_size,
                    0,
                    &confidence_size,
                    options.qq_field,
                    options.ci_alpha,
                )?;
                for subset in &subsets {
                    let (subset_size_for_row, subset_confidence_size) = named_subset_report_sizes(
                        subset,
                        options.subset_size,
                        options.whole_reference_size,
                        confidence_size_value,
                        options.stratification_sizes,
                        options.stratification_confidence_sizes,
                    );
                    let subset_confidence_size = subset_confidence_size
                        .map(|size| format!("{:.6}", size as f64))
                        .unwrap_or_default();
                    let all = all_counts
                        .by_subset_subtype
                        .get(subset)
                        .and_then(|types| types.get(variant_type))
                        .and_then(|subtypes| subtypes.get(subtype))
                        .cloned()
                        .unwrap_or_default();
                    let pass = pass_counts
                        .by_subset_subtype
                        .get(subset)
                        .and_then(|types| types.get(variant_type))
                        .and_then(|subtypes| subtypes.get(subtype))
                        .cloned()
                        .unwrap_or_default();
                    write_extended_pair(
                        &mut writer,
                        variant_type,
                        subtype,
                        subset,
                        all,
                        pass,
                        subset_size_for_row,
                        options
                            .stratification_levels
                            .get(subset)
                            .copied()
                            .unwrap_or(0),
                        &subset_confidence_size,
                        options.qq_field,
                        options.ci_alpha,
                    )?;
                }
            }
        }
    }

    writer
        .flush()
        .with_context(|| format!("failed to flush {}", path.display()))
}

pub(super) fn named_subset_report_sizes(
    subset: &str,
    _active_reference_size: usize,
    whole_reference_size: usize,
    confidence_size: Option<usize>,
    stratification_sizes: &BTreeMap<String, usize>,
    stratification_confidence_sizes: &BTreeMap<String, usize>,
) -> (usize, Option<usize>) {
    let subset_size = match subset {
        "TS_boundary" => whole_reference_size,
        "TS_contained" => confidence_size.unwrap_or(0),
        _ => stratification_sizes.get(subset).copied().unwrap_or(0),
    };
    let subset_confidence_size = match subset {
        "TS_boundary" | "TS_contained" => confidence_size,
        _ => confidence_size.map(|_| {
            stratification_confidence_sizes
                .get(subset)
                .copied()
                .unwrap_or(0)
        }),
    };
    (subset_size, subset_confidence_size)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write_extended_pair<W: Write>(
    writer: &mut W,
    variant_type: &str,
    subtype: &str,
    subset: &str,
    all: QuantifyTypeCounts,
    pass: QuantifyTypeCounts,
    subset_size: usize,
    subset_level: usize,
    conf_size: &str,
    qq_field: &str,
    ci_alpha: f64,
) -> Result<()> {
    write_extended_row(
        writer,
        variant_type,
        subtype,
        subset,
        "ALL",
        &all,
        subset_size,
        subset_level,
        conf_size,
        qq_field,
        ci_alpha,
    )?;
    write_extended_row(
        writer,
        variant_type,
        subtype,
        subset,
        "PASS",
        &pass,
        subset_size,
        subset_level,
        conf_size,
        qq_field,
        ci_alpha,
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn write_extended_row<W: Write>(
    writer: &mut W,
    variant_type: &str,
    subtype: &str,
    subset: &str,
    filter: &str,
    stats: &QuantifyTypeCounts,
    subset_size: usize,
    subset_level: usize,
    conf_size: &str,
    qq_field: &str,
    ci_alpha: f64,
) -> Result<()> {
    let supports_titv = variant_type == "SNP";
    let mut row = vec![
        variant_type.to_string(),
        subtype.to_string(),
        subset.to_string(),
        filter.to_string(),
        "*".to_string(),
        qq_field.to_string(),
        "*".to_string(),
        metric_ratio(stats.truth_tp.total, stats.truth_total.total),
        precision_metric(stats),
        metric_ratio(stats.query_unk.total, stats.query_total.total),
        f1_score(
            stats.truth_tp.total,
            stats.truth_total.total,
            stats.query_tp.total,
            stats.query_tp.total + stats.query_fp.total,
        ),
        stats.fp_gt.to_string(),
        stats.fp_al.to_string(),
        if subset == "*" {
            subset_size.to_string()
        } else {
            format!("{:.6}", subset_size as f64)
        },
        conf_size.to_string(),
        format!("{:.6}", subset_level as f64),
    ];
    append_stats(&mut row, &stats.truth_total, supports_titv);
    append_stats(&mut row, &stats.truth_tp, supports_titv);
    append_stats(&mut row, &stats.truth_fn, supports_titv);
    append_stats(&mut row, &stats.query_total, supports_titv);
    append_stats(&mut row, &stats.query_tp, supports_titv);
    append_stats(&mut row, &stats.query_fp, supports_titv);
    append_stats(&mut row, &stats.query_unk, supports_titv);
    if ci_alpha > 0.0 {
        report::append_ci_cells(
            &mut row,
            [
                (stats.truth_tp.total, stats.truth_total.total),
                (
                    stats.query_tp.total,
                    stats.query_tp.total + stats.query_fp.total,
                ),
                (stats.query_unk.total, stats.query_total.total),
            ],
            ci_alpha,
        );
    }
    writeln!(writer, "{}", row.join(","))?;
    Ok(())
}

pub(super) fn append_stats(row: &mut Vec<String>, stats: &CountsBucket, supports_titv: bool) {
    row.push(stats.total.to_string());
    if supports_titv {
        row.push(format_count(stats.ti));
        row.push(format_count(stats.tv));
    } else {
        row.push(".".to_string());
        row.push(".".to_string());
    }
    row.push(format_count(stats.het));
    row.push(format_count(stats.homalt));
    if supports_titv {
        row.push(ti_tv_ratio(stats.ti, stats.tv));
    } else {
        row.push(String::new());
    }
    row.push(het_hom_ratio(stats.het, stats.homalt));
}

pub(super) fn format_count(value: usize) -> String {
    format!("{:.6}", value as f64)
}

pub(super) fn metric_ratio(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        return "0.0".to_string();
    }
    format_metric(numerator as f64 / denominator as f64)
}

pub(super) fn precision_metric(stats: &TypeCounts) -> String {
    let denominator = stats.query_tp.total + stats.query_fp.total;
    if denominator == 0 && stats.query_total.total > 0 {
        String::new()
    } else {
        metric_ratio(stats.query_tp.total, denominator)
    }
}

pub(super) fn f1_score(
    truth_tp: usize,
    truth_total: usize,
    query_tp: usize,
    query_total: usize,
) -> String {
    let recall = if truth_total == 0 {
        0.0
    } else {
        truth_tp as f64 / truth_total as f64
    };
    let precision = if query_total == 0 {
        0.0
    } else {
        query_tp as f64 / query_total as f64
    };
    if (recall + precision).abs() < f64::EPSILON {
        return String::new();
    }
    format_metric((2.0 * recall * precision) / (recall + precision))
}

pub(super) fn ti_tv_ratio(ti: usize, tv: usize) -> String {
    if tv == 0 {
        return String::new();
    }
    format_ratio(ti as f64 / tv as f64)
}

pub(super) fn het_hom_ratio(het: usize, homalt: usize) -> String {
    if homalt == 0 {
        return String::new();
    }
    format_ratio(het as f64 / homalt as f64)
}

pub(super) fn format_metric(value: f64) -> String {
    report::format_metric(value)
}

pub(super) fn format_ratio(value: f64) -> String {
    report::python_repr_float(value)
}

pub(super) fn write_metrics_json(
    prefix: &Path,
    write_counts: bool,
    roc_indices: &roc::MetricIndices,
) -> Result<()> {
    let mut tables = vec![(
        "summary.metrics",
        "summary.metrics",
        suffixed_report_path(prefix, "summary.csv"),
    )];
    if write_counts {
        tables.push((
            "all.metrics",
            "all.metrics",
            suffixed_report_path(prefix, "extended.csv"),
        ));
    }
    for id in &roc_indices.table_order {
        let path = suffixed_report_path(prefix, &format!("{id}.csv.gz"));
        if path.is_file() {
            tables.push((id.as_str(), id.as_str(), path));
        }
    }
    let table_refs = tables
        .iter()
        .map(|(id, label, path)| (*id, *label, path.as_path()))
        .collect::<Vec<_>>();
    let commandline = std::env::args().collect::<Vec<_>>().join(" ");
    metrics_json::write_metrics_gz_for_module_with_indices(
        &suffixed_report_path(prefix, "metrics.json.gz"),
        "qfy.py.comparison",
        "qfy.py",
        &commandline,
        &table_refs,
        Some(&roc_indices.tables),
    )
}
