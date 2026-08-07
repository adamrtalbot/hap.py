//! Cohesive ROC rendering responsibility.

use super::EXTENDED_HEADER;
use super::accumulation::{
    EmittedRow, GroupAccum, ObsRecord, introsort_libstdcpp, is_aggregate_filter,
};
use super::model::RowKey;
use crate::adapters::report::{
    append_ci_cells, append_stats_with_missing, f1_score, format_count, het_hom_ratio,
    metric_ratio, ti_tv_ratio,
};
use crate::domain::CountsBucket;
use std::collections::{BTreeMap, HashSet};

#[derive(Clone, Copy, Debug)]
pub(super) enum RowFilter<'a> {
    /// Emit every row from every group — used for `result.roc.all.csv.gz`.
    All,
    /// Legacy `happyroc.py` Locations filter: keep only rows in the base
    /// group `(Type=ty, Subtype=*, Subset=*, Genotype=*, Filter=filter)`
    /// and drop the `QQ="*"` baseline row.
    Locations { ty: &'a str, filter: &'a str },
}

#[derive(Clone, Copy, Debug)]
pub(super) struct RenderConfig<'a> {
    pub(super) subset_size: usize,
    pub(super) whole_reference_size: usize,
    pub(super) conf_size: usize,
    pub(super) subset_sizes: &'a BTreeMap<String, usize>,
    pub(super) subset_confidence_sizes: &'a BTreeMap<String, usize>,
    pub(super) delta: f64,
    pub(super) ci_alpha: f64,
    pub(super) filter_counts_only: bool,
}

/// Build a per-subtype sorted-obs map keyed by `(Type, Subtype, Subset, Filter)`.
///
/// Legacy `RocOutput::write` iterates `(subtype, genotype)` pairs in a
/// fixed order and calls `Roc::getLevels(flag_mask)` per pair, where
/// `getLevels` runs `std::sort` on the shared obs vector. Since
/// `std::sort` is not idempotent, each call leaves obs in a slightly
/// different order at tied levels. The (Subtype, *) row written to the
/// roc.all CSV reflects the obs state AFTER the corresponding
/// `getLevels` call. Genotype iteration is `[het, hetalt, homalt, *]`,
/// so each subtype contributes 4 sorts.
///
/// For each `(Type, Subset, Filter)` wildcard group, sort the obs vector
/// progressively up to the maximum iteration count needed for that Type
/// (3×4=12 for SNP, 10×4=40 for INDEL), snapshotting at each subtype's
/// sort boundary. Each snapshot reproduces the obs state that legacy's
/// `Roc::getLevels` walks for that subtype.
pub(super) fn build_star_sorted(
    groups: &BTreeMap<RowKey, GroupAccum>,
) -> BTreeMap<(String, String, String, String), Vec<ObsRecord>> {
    let mut star_sorted: BTreeMap<(String, String, String, String), Vec<ObsRecord>> =
        BTreeMap::new();
    for (key, accum) in groups {
        if key.subtype == "*" && key.genotype == "*" {
            let max_count = match key.ty.as_str() {
                "SNP" => 12,
                "INDEL" => 40,
                _ => 4,
            };
            let mut obs = accum.obs.clone();
            let mut current_count: usize = 0;
            // Insert per-subtype snapshots at every 4-sort boundary.
            let snapshots: &[(usize, &str)] = match key.ty.as_str() {
                "SNP" => &[(4, "*"), (8, "ti"), (12, "tv")],
                "INDEL" => &[
                    (4, "*"),
                    (8, "I1_5"),
                    (12, "I6_15"),
                    (16, "I16_PLUS"),
                    (20, "D1_5"),
                    (24, "D6_15"),
                    (28, "D16_PLUS"),
                    (32, "C1_5"),
                    (36, "C6_15"),
                    (40, "C16_PLUS"),
                ],
                _ => &[(4, "*")],
            };
            for &(target, subtype_name) in snapshots {
                while current_count < target {
                    introsort_libstdcpp(&mut obs);
                    current_count += 1;
                }
                star_sorted.insert(
                    (
                        key.ty.clone(),
                        subtype_name.to_string(),
                        key.subset.clone(),
                        key.filter.clone(),
                    ),
                    obs.clone(),
                );
            }
            let _ = max_count; // silence unused warning when not needed
        }
    }
    star_sorted
}

pub(super) fn render_rows(
    groups: &BTreeMap<RowKey, GroupAccum>,
    star_sorted: &BTreeMap<(String, String, String, String), Vec<ObsRecord>>,
    row_filter: RowFilter<'_>,
    config: RenderConfig<'_>,
) -> Vec<String> {
    let mut out = Vec::new();
    // `accumulate` pre-seeds empty subtype buckets for types that are present,
    // because legacy reports zero-valued subtype baselines for an observed
    // type. It does not, however, emit a second family of rows for a wholly
    // absent type. Determine presence from actual observations before walking
    // the pre-seeded map.
    let active_types: HashSet<&str> = groups
        .iter()
        .filter(|(_, accum)| !accum.obs.is_empty())
        .map(|(key, _)| key.ty.as_str())
        .collect();
    for (key, accum) in groups {
        if !active_types.contains(key.ty.as_str()) {
            continue;
        }
        match row_filter {
            RowFilter::All => {}
            RowFilter::Locations { ty, filter } => {
                if key.ty != ty
                    || key.subtype != "*"
                    || key.subset != "*"
                    || key.genotype != "*"
                    || key.filter != filter
                {
                    continue;
                }
            }
        }
        let is_filter_tier = !is_aggregate_filter(&key.filter);
        let emitted_rows: Vec<EmittedRow> = if key.subtype == "*" {
            accum.emit_with_delta(config.delta)
        } else if let Some(shared) = star_sorted.get(&(
            key.ty.clone(),
            key.subtype.clone(),
            key.subset.clone(),
            key.filter.clone(),
        )) {
            accum.emit_with_shared_sort_and_delta(shared, &key.subtype, config.delta)
        } else {
            accum.emit_with_delta(config.delta)
        };
        for emitted in emitted_rows {
            if matches!(row_filter, RowFilter::Locations { .. }) && emitted.qq_str == "*" {
                // Legacy's Locations file drops the baseline row.
                continue;
            }
            if is_filter_tier && config.filter_counts_only && emitted.qq_str != "*" {
                // Per-Filter rows in legacy are baseline-only (QQ='*').
                // No numeric thresholds — skip the synthetic 0.0 and any
                // real numeric buckets that may have landed here.
                continue;
            }
            out.push(render_row(
                key,
                &emitted,
                config.subset_size,
                config.whole_reference_size,
                config.conf_size,
                config.subset_sizes,
                config.subset_confidence_sizes,
                is_filter_tier && config.filter_counts_only,
                config.ci_alpha,
            ));
        }
    }
    out
}

#[allow(clippy::too_many_arguments)]
pub(super) fn render_row(
    key: &RowKey,
    emitted: &EmittedRow,
    subset_size: usize,
    whole_reference_size: usize,
    conf_size: usize,
    subset_sizes: &BTreeMap<String, usize>,
    subset_confidence_sizes: &BTreeMap<String, usize>,
    counts_only: bool,
    ci_alpha: f64,
) -> String {
    let counts = &emitted.cum;
    let truth_total = counts.truth_total();
    let query_total = counts.query_total();

    // Per-Filter rows (Filter ∉ {ALL, PASS}): legacy zeroes TRUTH.TOTAL,
    // TRUTH.FN, QUERY.TOTAL and renders METRIC.* as the pandas NaN
    // sentinel (empty cell). Only TRUTH.TP / QUERY.TP / QUERY.FP /
    // QUERY.UNK carry real counts (and their het/homalt substats); the
    // TOTAL / FN blocks emit `0` for the count and `.` for substats.
    let is_filter_tier = counts_only;
    let supports_titv = key.ty == "SNP";
    let missing_titv = if conf_size > 0 { "." } else { "" };
    let mut row = Vec::with_capacity(EXTENDED_HEADER.len());

    row.push(key.ty.clone());
    row.push(key.subtype.clone());
    row.push(key.subset.clone());
    row.push(key.filter.clone());
    row.push(key.genotype.clone());
    row.push(key.qq_field.clone());
    row.push(emitted.qq_str.clone());

    // Metrics — reuse the extended.csv helpers so formatting matches byte
    // for byte (pandas repr for finite floats, `"0.0"` for zero-denominator
    // ratios, empty string for zero-denominator F1).
    if is_filter_tier {
        row.push(String::new()); // METRIC.Recall
        row.push(String::new()); // METRIC.Precision
        row.push(String::new()); // METRIC.Frac_NA
        row.push(String::new()); // METRIC.F1_Score
    } else {
        row.push(metric_ratio(counts.truth_tp.total, truth_total.total));
        // Use precision_ratio: returns empty when query records exist but
        // none are TP/FP (all UNKs), 0.0 when query_total == 0.
        row.push(crate::adapters::report::precision_ratio(
            counts.query_tp.total,
            counts.query_fp.total,
            query_total.total,
        ));
        row.push(metric_ratio(counts.query_unk.total, query_total.total));
        // F1 score: empty when precision is undefined (query records exist
        // but precision denominator is 0 — same condition as precision_ratio
        // returning empty).
        let prec_denom = counts.query_tp.total + counts.query_fp.total;
        if prec_denom == 0 && query_total.total > 0 {
            row.push(String::new());
        } else {
            row.push(f1_score(
                counts.truth_tp.total,
                truth_total.total,
                counts.query_tp.total,
                prec_denom,
            ));
        }
    }

    // FP.gt / FP.al: raw integers (`"11"`, not `"11.000000"`). Legacy
    // emits them at every INDEL stratification row (matches extended.csv's
    // per-subtype FP emission) and at SNP non-base rows. Forward the
    // accumulated counters as-is for every row.
    row.push(counts.fp_gt.to_string());
    row.push(counts.fp_al.to_string());

    // Subset.Size / Subset.IS_CONF.Size / Subset.Level — match the per-row
    // convention of write_extended exactly.
    let (sz, conf) = subset_size_cells(
        &key.subset,
        &key.subtype,
        subset_size,
        whole_reference_size,
        conf_size,
        subset_sizes,
        subset_confidence_sizes,
    );
    row.push(sz);
    row.push(conf);
    row.push("0.000000".to_string());

    // Substat blocks (7 cells each). The baseline `QQ="*"` row gets full
    // ti/tv/het/homalt/ratio cells via `append_stats`, mirroring the
    // extended.csv layout. Per-threshold numeric rows mirror legacy's
    // independent-sweep pattern: ti/tv/het/homalt cells show cumulative
    // values only at levels kept by that substat's own roc-delta sweep;
    // otherwise `.` (count) / `""` (ratio). Ratios are computed in-row
    // from whichever pair is present — if either side is missing, the
    // ratio is blank (legacy's NaN/NaN → NaN rendering via pandas).
    let emit_block = |row: &mut Vec<String>, bucket: &CountsBucket| match &emitted.substats {
        None if is_filter_tier || emitted.qq_str == "*" => {
            append_stats_with_missing(row, bucket, supports_titv, missing_titv)
        }
        None => append_roc_stats(row, bucket, supports_titv, missing_titv),
        Some(avail) => {
            row.push(bucket.total.to_string());
            emit_substat_cell(
                row,
                bucket.ti,
                supports_titv && avail.ti,
                supports_titv,
                missing_titv,
            );
            emit_substat_cell(
                row,
                bucket.tv,
                supports_titv && avail.tv,
                supports_titv,
                missing_titv,
            );
            emit_substat_cell(row, bucket.het, avail.het, true, missing_titv);
            emit_substat_cell(row, bucket.homalt, avail.homalt, true, missing_titv);
            if supports_titv && avail.ti && avail.tv {
                row.push(ti_tv_ratio(bucket.ti, bucket.tv));
            } else {
                row.push(String::new());
            }
            if avail.het && avail.homalt {
                row.push(het_hom_ratio(bucket.het, bucket.homalt));
            } else {
                row.push(String::new());
            }
        }
    };
    // Per-Filter rows render TRUTH.TOTAL / TRUTH.FN / QUERY.TOTAL as a
    // hollow "0 + . + …" block: the count cell is "0", the substat
    // cells are "." (or empty for ratios). Mirror that with
    // `emit_zero_block` instead of the normal `emit_block` for those
    // three groups; TRUTH.TP, QUERY.TP, QUERY.FP, QUERY.UNK still
    // carry their real counts.
    let emit_zero_block = |row: &mut Vec<String>| {
        row.push("0".to_string());
        for _ in 0..4 {
            row.push(".".to_string());
        }
        row.push(String::new());
        row.push(String::new());
    };
    if is_filter_tier {
        emit_zero_block(&mut row);
        emit_block(&mut row, &counts.truth_tp);
        emit_zero_block(&mut row);
        emit_zero_block(&mut row);
        emit_block(&mut row, &counts.query_tp);
        emit_block(&mut row, &counts.query_fp);
        emit_block(&mut row, &counts.query_unk);
    } else {
        emit_block(&mut row, &truth_total);
        emit_block(&mut row, &counts.truth_tp);
        emit_block(&mut row, &counts.truth_fn);
        emit_block(&mut row, &query_total);
        emit_block(&mut row, &counts.query_tp);
        emit_block(&mut row, &counts.query_fp);
        emit_block(&mut row, &counts.query_unk);
    }

    if ci_alpha > 0.0 {
        let observations = [
            (counts.truth_tp.total, truth_total.total),
            (
                counts.query_tp.total,
                counts.query_tp.total + counts.query_fp.total,
            ),
            if is_filter_tier {
                // Filter-tier projection zeroes QUERY.TOTAL before the
                // unknown-fraction CI is calculated. Recall and precision
                // still use the retained TP/FN/FP counts above.
                (0, 0)
            } else {
                (counts.query_unk.total, query_total.total)
            },
        ];
        append_ci_cells(&mut row, observations, ci_alpha);
    }

    row.join(",")
}

pub(super) fn emit_substat_cell(
    row: &mut Vec<String>,
    value: usize,
    kept: bool,
    supported: bool,
    missing: &str,
) {
    if !supported {
        row.push(missing.to_string());
    } else if kept {
        row.push(format_count(value));
    } else {
        row.push(".".to_string());
    }
}

pub(super) fn append_roc_stats(
    row: &mut Vec<String>,
    bucket: &CountsBucket,
    supports_titv: bool,
    missing_titv: &str,
) {
    row.push(bucket.total.to_string());
    if supports_titv {
        row.push(format_count(bucket.ti));
        row.push(format_count(bucket.tv));
    } else {
        row.push(missing_titv.to_string());
        row.push(missing_titv.to_string());
    }
    row.push(format_count(bucket.het));
    row.push(format_count(bucket.homalt));
    row.push(if supports_titv {
        ti_tv_ratio(bucket.ti, bucket.tv)
    } else {
        String::new()
    });
    row.push(het_hom_ratio(bucket.het, bucket.homalt));
}

pub(super) fn subset_size_cells(
    subset: &str,
    _subtype: &str,
    subset_size: usize,
    whole_reference_size: usize,
    conf_size: usize,
    subset_sizes: &BTreeMap<String, usize>,
    subset_confidence_sizes: &BTreeMap<String, usize>,
) -> (String, String) {
    // Derived from the four branches in report::write_extended. Kept in
    // lockstep — if write_extended changes its Subset.Size/IS_CONF.Size
    // convention, this function must follow.
    let is_base_subset = subset == "*";
    let size_cell = if is_base_subset {
        // Subset="*": always the raw subset_size integer, regardless of
        // subtype.
        subset_size.to_string()
    } else if subset == "TS_contained" {
        format_count(conf_size)
    } else if subset == "TS_boundary" {
        format_count(whole_reference_size)
    } else {
        // User-named stratification.
        format_count(subset_sizes.get(subset).copied().unwrap_or(0))
    };
    // The built-in confidence subsets carry the global confidence size.
    // User-named subsets instead carry their interval-union intersection
    // with the confidence regions. Empty cell only when confidence is absent.
    let conf_cell = if conf_size > 0 {
        let size = if matches!(subset, "*" | "TS_boundary" | "TS_contained") {
            conf_size
        } else {
            subset_confidence_sizes.get(subset).copied().unwrap_or(0)
        };
        format_count(size)
    } else {
        String::new()
    };
    (size_cell, conf_cell)
}
