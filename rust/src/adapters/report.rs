use crate::domain::{CountsBucket, TypeCounts};
use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

pub(crate) const EXTENDED_HEADER: [&str; 65] = [
    "Type",
    "Subtype",
    "Subset",
    "Filter",
    "Genotype",
    "QQ.Field",
    "QQ",
    "METRIC.Recall",
    "METRIC.Precision",
    "METRIC.Frac_NA",
    "METRIC.F1_Score",
    "FP.gt",
    "FP.al",
    "Subset.Size",
    "Subset.IS_CONF.Size",
    "Subset.Level",
    "TRUTH.TOTAL",
    "TRUTH.TOTAL.ti",
    "TRUTH.TOTAL.tv",
    "TRUTH.TOTAL.het",
    "TRUTH.TOTAL.homalt",
    "TRUTH.TOTAL.TiTv_ratio",
    "TRUTH.TOTAL.het_hom_ratio",
    "TRUTH.TP",
    "TRUTH.TP.ti",
    "TRUTH.TP.tv",
    "TRUTH.TP.het",
    "TRUTH.TP.homalt",
    "TRUTH.TP.TiTv_ratio",
    "TRUTH.TP.het_hom_ratio",
    "TRUTH.FN",
    "TRUTH.FN.ti",
    "TRUTH.FN.tv",
    "TRUTH.FN.het",
    "TRUTH.FN.homalt",
    "TRUTH.FN.TiTv_ratio",
    "TRUTH.FN.het_hom_ratio",
    "QUERY.TOTAL",
    "QUERY.TOTAL.ti",
    "QUERY.TOTAL.tv",
    "QUERY.TOTAL.het",
    "QUERY.TOTAL.homalt",
    "QUERY.TOTAL.TiTv_ratio",
    "QUERY.TOTAL.het_hom_ratio",
    "QUERY.TP",
    "QUERY.TP.ti",
    "QUERY.TP.tv",
    "QUERY.TP.het",
    "QUERY.TP.homalt",
    "QUERY.TP.TiTv_ratio",
    "QUERY.TP.het_hom_ratio",
    "QUERY.FP",
    "QUERY.FP.ti",
    "QUERY.FP.tv",
    "QUERY.FP.het",
    "QUERY.FP.homalt",
    "QUERY.FP.TiTv_ratio",
    "QUERY.FP.het_hom_ratio",
    "QUERY.UNK",
    "QUERY.UNK.ti",
    "QUERY.UNK.tv",
    "QUERY.UNK.het",
    "QUERY.UNK.homalt",
    "QUERY.UNK.TiTv_ratio",
    "QUERY.UNK.het_hom_ratio",
];

/// Append a legacy report suffix without treating a dotted prefix as an extension.
pub(crate) fn suffixed_report_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(".");
    path.push(suffix);
    PathBuf::from(path)
}

type SubsetSubtypeFpCounts = BTreeMap<String, BTreeMap<String, BTreeMap<String, (usize, usize)>>>;

pub(crate) fn write_summary(
    path: &Path,
    all_counts: &BTreeMap<String, TypeCounts>,
    pass_counts: &BTreeMap<String, TypeCounts>,
    all_fp: &BTreeMap<String, (usize, usize)>,
    pass_fp: &BTreeMap<String, (usize, usize)>,
) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
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
    let mut variant_types: Vec<&String> = all_counts
        .keys()
        .chain(pass_counts.keys())
        .filter(|variant_type| {
            summary_type_has_observations(variant_type, all_counts, pass_counts, all_fp, pass_fp)
        })
        .collect();
    variant_types.sort();
    variant_types.dedup();
    let empty = TypeCounts::default();
    for variant_type in variant_types {
        for filter in ["ALL", "PASS"] {
            let stats = if filter == "ALL" {
                all_counts.get(variant_type).unwrap_or(&empty)
            } else {
                pass_counts.get(variant_type).unwrap_or(&empty)
            };
            let fp_tuple = if filter == "ALL" {
                all_fp.get(variant_type).copied().unwrap_or((0, 0))
            } else {
                pass_fp.get(variant_type).copied().unwrap_or((0, 0))
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
                fp_tuple.0,
                fp_tuple.1,
                metric_ratio(stats.truth_tp.total, stats.truth_total.total),
                precision_ratio(
                    stats.query_tp.total,
                    stats.query_fp.total,
                    stats.query_total.total,
                ),
                metric_ratio(stats.query_unk.total, stats.query_total.total),
                f1_score(
                    stats.truth_tp.total,
                    stats.truth_total.total,
                    stats.query_tp.total,
                    stats.query_tp.total + stats.query_fp.total
                ),
                ti_tv_ratio(stats.truth_total.ti, stats.truth_total.tv),
                ti_tv_ratio(stats.query_total.ti, stats.query_total.tv),
                het_hom_ratio(stats.truth_total.het, stats.truth_total.homalt),
                het_hom_ratio(stats.query_total.het, stats.query_total.homalt),
            )?;
        }
    }
    writer
        .flush()
        .with_context(|| format!("failed to flush {}", path.display()))
}

#[allow(clippy::too_many_arguments)] // Each map is a distinct legacy report axis/tier.
pub(crate) fn write_extended(
    path: &Path,
    all_counts: &BTreeMap<String, TypeCounts>,
    pass_counts: &BTreeMap<String, TypeCounts>,
    all_subtype: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
    pass_subtype: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
    subset_size: usize,
    whole_reference_size: usize,
    conf_size: usize,
    has_conf_regions: bool,
    all_subset: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
    pass_subset: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
    all_subset_subtype: &BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>>,
    pass_subset_subtype: &BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>>,
    all_fp: &BTreeMap<String, (usize, usize)>,
    pass_fp: &BTreeMap<String, (usize, usize)>,
    all_subset_fp: &BTreeMap<String, BTreeMap<String, (usize, usize)>>,
    pass_subset_fp: &BTreeMap<String, BTreeMap<String, (usize, usize)>>,
    all_subtype_fp: &BTreeMap<String, BTreeMap<String, (usize, usize)>>,
    pass_subtype_fp: &BTreeMap<String, BTreeMap<String, (usize, usize)>>,
    all_subset_subtype_fp: &SubsetSubtypeFpCounts,
    pass_subset_subtype_fp: &SubsetSubtypeFpCounts,
) -> Result<()> {
    let mut writer = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(writer, "{}", EXTENDED_HEADER.join(","))?;

    if all_counts.is_empty() && pass_counts.is_empty() {
        for row in empty_comparison_extended_lines(subset_size) {
            writeln!(writer, "{row}")?;
        }
        return writer
            .flush()
            .with_context(|| format!("failed to flush {}", path.display()));
    }

    let mut variant_types: Vec<&String> = all_counts.keys().chain(pass_counts.keys()).collect();
    variant_types.sort();
    variant_types.dedup();
    let pick = |filter: &str,
                counts: &BTreeMap<String, TypeCounts>,
                pass: &BTreeMap<String, TypeCounts>,
                ty: &str|
     -> TypeCounts {
        let src = if filter == "ALL" { counts } else { pass };
        src.get(ty).cloned().unwrap_or_default()
    };
    let pick_nested = |filter: &str,
                       all: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
                       pass: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
                       outer: &str,
                       inner: &str|
     -> TypeCounts {
        let src = if filter == "ALL" { all } else { pass };
        src.get(outer)
            .and_then(|inner_map| inner_map.get(inner))
            .cloned()
            .unwrap_or_default()
    };
    // Legacy `qfy.py` emits extended.csv rows in the lex-ASC tuple order
    // `(Type, Subtype, Subset, Filter)`. Each Subtype's full Subset/Filter
    // rectangle is emitted before moving to the next Subtype — i.e. the
    // outer key is Subtype, the inner keys are Subset then Filter. Match
    // that ordering exactly so byte-MD5 nf-test comparisons pass.
    for variant_type in &variant_types {
        let is_indel = *variant_type == "INDEL";
        let supports_titv = *variant_type == "SNP";
        let subtypes: Vec<&str> = if is_indel {
            vec![
                "*", "C16_PLUS", "C1_5", "C6_15", "D16_PLUS", "D1_5", "D6_15", "I16_PLUS", "I1_5",
                "I6_15",
            ]
        } else {
            vec!["*"]
        };
        for subtype in &subtypes {
            let is_base_subtype = *subtype == "*";
            let subsets = if has_conf_regions {
                &["*", "TS_boundary", "TS_contained"][..]
            } else {
                &["*"][..]
            };
            for subset in subsets {
                let subset = *subset;
                let is_base_subset = subset == "*";
                // Legacy emission gates:
                //   • base subset: always (INDEL subtype rows are seeded for
                //     every INDEL_SUBTYPE).
                //   • named subset: only if that subset carries this Type.
                //     Once the Type axis exists, every INDEL subtype/filter
                //     row is seeded, including zero-count subtype rows.
                if !is_base_subset {
                    let has_any = all_subset
                        .get(subset)
                        .and_then(|m| m.get(variant_type.as_str()))
                        .is_some()
                        || pass_subset
                            .get(subset)
                            .and_then(|m| m.get(variant_type.as_str()))
                            .is_some();
                    if !has_any {
                        continue;
                    }
                }

                for filter in ["ALL", "PASS"] {
                    // Resolve the right TypeCounts for this 4-tuple.
                    let stats: TypeCounts = match (is_base_subtype, is_base_subset) {
                        (true, true) => pick(filter, all_counts, pass_counts, variant_type),
                        (true, false) => {
                            pick_nested(filter, all_subset, pass_subset, subset, variant_type)
                        }
                        (false, true) => {
                            pick_nested(filter, all_subtype, pass_subtype, variant_type, subtype)
                        }
                        (false, false) => {
                            let src = if filter == "ALL" {
                                all_subset_subtype
                            } else {
                                pass_subset_subtype
                            };
                            src.get(subset)
                                .and_then(|m| m.get(variant_type.as_str()))
                                .and_then(|m| m.get(*subtype))
                                .cloned()
                                .unwrap_or_default()
                        }
                    };

                    // Same indexing for FP.gt / FP.al.
                    let fp_tuple: (usize, usize) = match (is_base_subtype, is_base_subset) {
                        (true, true) => {
                            if filter == "ALL" {
                                all_fp.get(variant_type.as_str()).copied().unwrap_or((0, 0))
                            } else {
                                pass_fp
                                    .get(variant_type.as_str())
                                    .copied()
                                    .unwrap_or((0, 0))
                            }
                        }
                        (true, false) => {
                            let src = if filter == "ALL" {
                                all_subset_fp
                            } else {
                                pass_subset_fp
                            };
                            src.get(subset)
                                .and_then(|m| m.get(variant_type.as_str()))
                                .copied()
                                .unwrap_or((0, 0))
                        }
                        (false, true) => {
                            let src = if filter == "ALL" {
                                all_subtype_fp
                            } else {
                                pass_subtype_fp
                            };
                            src.get(variant_type.as_str())
                                .and_then(|m| m.get(*subtype))
                                .copied()
                                .unwrap_or((0, 0))
                        }
                        (false, false) => {
                            let src = if filter == "ALL" {
                                all_subset_subtype_fp
                            } else {
                                pass_subset_subtype_fp
                            };
                            src.get(subset)
                                .and_then(|m| m.get(variant_type.as_str()))
                                .and_then(|m| m.get(*subtype))
                                .copied()
                                .unwrap_or((0, 0))
                        }
                    };

                    // Subset.Size cell varies by (subtype × subset). Match
                    // legacy:
                    //   subset="*"             → raw `subset_size`
                    //                              (formatted as integer
                    //                              when subtype="*",
                    //                              else still raw integer)
                    //   subset="TS_boundary"   → format_count(whole_reference_size)
                    //   subset="TS_contained"  → format_count(conf_size)
                    let size_cell = match subset {
                        "*" => subset_size.to_string(),
                        "TS_contained" => format_count(conf_size),
                        "TS_boundary" => format_count(whole_reference_size),
                        _ => format_count(subset_size),
                    };
                    // Subset.IS_CONF.Size cell:
                    //   subtype="*", subset="*"  → format_count(conf_size)
                    //                              when conf_size > 0
                    //                              else empty
                    //   subtype="*", subset!="*" → format_count(conf_size)
                    //                              when conf_size > 0
                    //                              else empty
                    //   subtype!="*", subset="*" → empty cell at base subset
                    //                              when conf_size==0; else
                    //                              format_count(conf_size)
                    //   subtype!="*", subset!="*"→ format_count(conf_size)
                    //                              when conf_size > 0
                    let conf_cell = if conf_size > 0 {
                        format_count(conf_size)
                    } else {
                        String::new()
                    };

                    let mut row = vec![
                        variant_type.to_string(),
                        subtype.to_string(),
                        subset.to_string(),
                        filter.to_string(),
                        "*".to_string(),
                        "QUAL".to_string(),
                        "*".to_string(),
                        metric_ratio(stats.truth_tp.total, stats.truth_total.total),
                        precision_ratio(
                            stats.query_tp.total,
                            stats.query_fp.total,
                            stats.query_total.total,
                        ),
                        metric_ratio(stats.query_unk.total, stats.query_total.total),
                        f1_score(
                            stats.truth_tp.total,
                            stats.truth_total.total,
                            stats.query_tp.total,
                            stats.query_tp.total + stats.query_fp.total,
                        ),
                        fp_tuple.0.to_string(),
                        fp_tuple.1.to_string(),
                        size_cell,
                        conf_cell,
                        "0.000000".to_string(),
                    ];

                    // Substat blocks: SNP rows expose ti/tv; INDEL rows
                    // (including INDEL subtypes) leave them empty.
                    let titv = supports_titv && is_base_subtype;
                    let missing_titv = if has_conf_regions { "." } else { "" };
                    append_extended_stats(&mut row, &stats.truth_total, titv, missing_titv);
                    append_extended_stats(&mut row, &stats.truth_tp, titv, missing_titv);
                    append_extended_stats(&mut row, &stats.truth_fn, titv, missing_titv);
                    append_extended_stats(&mut row, &stats.query_total, titv, missing_titv);
                    append_extended_stats(&mut row, &stats.query_tp, titv, missing_titv);
                    append_extended_stats(&mut row, &stats.query_fp, titv, missing_titv);
                    append_extended_stats(&mut row, &stats.query_unk, titv, missing_titv);

                    writeln!(writer, "{}", row.join(","))?;
                }
            }
        }
    }
    writer
        .flush()
        .with_context(|| format!("failed to flush {}", path.display()))
}

fn summary_type_has_observations(
    variant_type: &str,
    all_counts: &BTreeMap<String, TypeCounts>,
    pass_counts: &BTreeMap<String, TypeCounts>,
    all_fp: &BTreeMap<String, (usize, usize)>,
    pass_fp: &BTreeMap<String, (usize, usize)>,
) -> bool {
    let counted = |counts: &BTreeMap<String, TypeCounts>| {
        counts.get(variant_type).is_some_and(|stats| {
            stats.truth_total.total > 0
                || stats.query_total.total > 0
                || stats.truth_tp.total > 0
                || stats.truth_fn.total > 0
                || stats.query_tp.total > 0
                || stats.query_fp.total > 0
                || stats.query_unk.total > 0
        })
    };
    let classified_fp = |counts: &BTreeMap<String, (usize, usize)>| {
        counts
            .get(variant_type)
            .is_some_and(|(gt, allele)| *gt > 0 || *allele > 0)
    };
    counted(all_counts) || counted(pass_counts) || classified_fp(all_fp) || classified_fp(pass_fp)
}

pub(crate) fn empty_comparison_extended_lines(subset_size: usize) -> Vec<String> {
    ["INDEL", "SNP"]
        .into_iter()
        .map(|variant_type| {
            let mut row = vec![String::new(); EXTENDED_HEADER.len()];
            row[0] = variant_type.to_string();
            row[1] = "*".to_string();
            row[2] = "*".to_string();
            row[3] = "ALL".to_string();
            row[4] = "*".to_string();
            row[5] = "nan".to_string();
            row[6] = "*".to_string();
            for index in [11, 12, 16, 23, 30, 37, 44, 51, 58] {
                row[index] = "0".to_string();
            }
            row[13] = format!("{subset_size}.0");
            row.join(",")
        })
        .collect()
}

fn append_extended_stats(
    row: &mut Vec<String>,
    stats: &CountsBucket,
    supports_titv: bool,
    missing_titv: &str,
) {
    append_stats_with_missing(row, stats, supports_titv, missing_titv);
}

pub(crate) fn append_stats_with_missing(
    row: &mut Vec<String>,
    stats: &CountsBucket,
    supports_titv: bool,
    missing_titv: &str,
) {
    row.push(stats.total.to_string());
    if supports_titv {
        row.push(format_count(stats.ti));
        row.push(format_count(stats.tv));
    } else {
        row.push(missing_titv.to_string());
        row.push(missing_titv.to_string());
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

pub(crate) fn format_count(value: usize) -> String {
    format!("{:.6}", value as f64)
}

pub(crate) fn metric_ratio(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        return "0.0".to_string();
    }
    format_metric(numerator as f64 / denominator as f64)
}

/// Legacy quantify renders METRIC.Precision as the pandas NaN sentinel
/// (empty cell) when there are query records on the row but none of them
/// classify as TP or FP — i.e. all of them are UNKs. The 0/0 divide
/// surfaces as NaN through pandas because there's data on the row, just
/// nothing in the precision denominator. When the row has no query
/// records at all (`query_total == 0`) legacy keeps 0.0, matching the
/// flat metric_ratio for unobserved cells.
pub(crate) fn precision_ratio(query_tp: usize, query_fp: usize, query_total: usize) -> String {
    let denom = query_tp + query_fp;
    if denom == 0 {
        if query_total == 0 {
            return "0.0".to_string();
        }
        return String::new();
    }
    format_metric(query_tp as f64 / denom as f64)
}

pub(crate) fn f1_score(
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

pub(crate) fn ti_tv_ratio(ti: usize, tv: usize) -> String {
    if tv == 0 {
        return String::new();
    }
    format_ratio(ti as f64 / tv as f64)
}

pub(crate) fn het_hom_ratio(het: usize, homalt: usize) -> String {
    if homalt == 0 {
        return String::new();
    }
    format_ratio(het as f64 / homalt as f64)
}

pub(crate) fn format_metric(value: f64) -> String {
    // Legacy data flow for METRIC.* columns: C++ writes the raw double via
    // `std::to_string` (printf "%.6f" — six fractional decimals), Python
    // pandas re-parses that string with its lossy `xstrtod` (a power-of-2
    // multiply/divide walk that introduces specific 1-ULP differences vs.
    // a correctly-rounded strtod), then pandas writes the resulting f64
    // via Python 2's 12-significant-digit `str(float)`. Mirroring those
    // three steps in Rust reproduces the exact decimal text legacy emits
    // without exposing the parser's adjacent-double representation.
    if !value.is_finite() {
        return python_repr_float(value);
    }
    let formatted = format!("{value:.6}");
    let lossy = pandas_xstrtod(&formatted);
    python_repr_float(lossy)
}

/// Rust port of pandas 0.24 `xstrtod` (the lossy parser used by
/// `pandas.to_numeric`). Pandas applies the scale factor by walking the
/// binary representation of the (positive) exponent magnitude and
/// repeatedly multiplying or dividing by `p10 = 10, 100, 10000, ...`.
/// The intermediate rounding produces results that differ by 1 ULP from
/// `f64::from_str` on roughly half of inputs — for byte parity we have
/// to reproduce that walk.
pub(crate) fn pandas_xstrtod(s: &str) -> f64 {
    let bytes = s.as_bytes();
    let mut idx = 0;
    let mut sign = 1.0_f64;
    if idx < bytes.len() && bytes[idx] == b'-' {
        sign = -1.0;
        idx += 1;
    } else if idx < bytes.len() && bytes[idx] == b'+' {
        idx += 1;
    }
    let mut number: f64 = 0.0;
    while idx < bytes.len() && bytes[idx].is_ascii_digit() {
        number = number * 10.0 + (bytes[idx] - b'0') as f64;
        idx += 1;
    }
    let mut num_decimals: i32 = 0;
    if idx < bytes.len() && bytes[idx] == b'.' {
        idx += 1;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            number = number * 10.0 + (bytes[idx] - b'0') as f64;
            num_decimals += 1;
            idx += 1;
        }
    }
    let mut exponent: i32 = -num_decimals;
    if idx < bytes.len() && (bytes[idx] == b'e' || bytes[idx] == b'E') {
        idx += 1;
        let mut neg = false;
        if idx < bytes.len() && bytes[idx] == b'-' {
            neg = true;
            idx += 1;
        } else if idx < bytes.len() && bytes[idx] == b'+' {
            idx += 1;
        }
        let mut exp_digits = 0i32;
        while idx < bytes.len() && bytes[idx].is_ascii_digit() {
            exp_digits = exp_digits * 10 + (bytes[idx] - b'0') as i32;
            idx += 1;
        }
        if neg {
            exponent -= exp_digits;
        } else {
            exponent += exp_digits;
        }
    }
    number *= sign;
    let mut p10: f64 = 10.0;
    let mut n = exponent.unsigned_abs() as i32;
    while n != 0 {
        if n & 1 == 1 {
            if exponent < 0 {
                number /= p10;
            } else {
                number *= p10;
            }
        }
        n >>= 1;
        p10 *= p10;
    }
    number
}

pub(crate) fn format_ratio(value: f64) -> String {
    python_repr_float(value)
}

pub(crate) fn full_repr_float(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            "-inf".to_string()
        } else {
            "inf".to_string()
        };
    }
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0.0".to_string()
        } else {
            "0.0".to_string()
        };
    }
    // Python's shortest representation uses fixed notation exactly within
    // this numeric interval. Direct comparison preserves adjacent f64 values;
    // `log10` can round them across either boundary.
    let magnitude = value.abs();
    if !(1e-4..1e16).contains(&magnitude) {
        let scientific = format!("{value:e}");
        let (mantissa, exponent_text) = scientific.split_once('e').unwrap();
        let exponent_value = exponent_text.parse::<i32>().unwrap();
        return format!("{mantissa}e{exponent_value:+03}");
    }
    let mut rendered = value.to_string();
    if !rendered.contains('.') {
        rendered.push_str(".0");
    }
    rendered
}

// Python 2's str(float), used by the other legacy adapters, keeps 12
// significant digits. It uses fixed notation for rounded exponents in
// [-4, 10] and scientific notation otherwise. Integer-valued fixed numbers
// retain `.0`.
pub(crate) fn python_repr_float(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            "-inf".to_string()
        } else {
            "inf".to_string()
        };
    }
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0.0".to_string()
        } else {
            "0.0".to_string()
        };
    }
    let scientific = format!("{value:.11e}");
    let (mantissa, exponent_text) = scientific.split_once('e').unwrap();
    let exponent = exponent_text.parse::<i32>().unwrap();
    if !(-4..11).contains(&exponent) {
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        return format!("{mantissa}e{exponent:+03}");
    }
    let decimal_places = usize::try_from(11 - exponent).unwrap_or(0);
    let mut rendered = format!("{value:.decimal_places$}");
    if rendered.contains('.') {
        while rendered.ends_with('0') {
            rendered.pop();
        }
        if rendered.ends_with('.') {
            rendered.push('0');
        }
    } else {
        rendered.push_str(".0");
    }
    rendered
}

pub(crate) fn append_ci_cells(
    row: &mut Vec<String>,
    observations: [(usize, usize); 3],
    alpha: f64,
) {
    for (successes, trials) in observations {
        let (lower, upper) = crate::domain::jeffreys_interval(successes, trials, alpha);
        row.push(python_repr_float(lower));
        row.push(python_repr_float(upper));
    }
}

#[cfg(test)]
mod format_tests {
    use super::{
        SubsetSubtypeFpCounts, format_ratio, full_repr_float, python_repr_float, write_extended,
        write_summary,
    };
    use crate::domain::TypeCounts;
    use std::collections::BTreeMap;

    #[test]
    fn integer_valued_float_keeps_trailing_zero() {
        assert_eq!(python_repr_float(1.0), "1.0");
        assert_eq!(python_repr_float(0.0), "0.0");
        assert_eq!(python_repr_float(-1.0), "-1.0");
    }

    #[test]
    fn short_decimal_stays_short() {
        assert_eq!(python_repr_float(0.1), "0.1");
        assert_eq!(python_repr_float(0.861224), "0.861224");
    }

    #[test]
    fn germline_ratio_matches_python_two_twelve_significant_digits() {
        let v: f64 = 0.9651790000000001;
        assert_eq!(format_ratio(v), "0.965179");
        let w: f64 = 1.58980044345898;
        assert_eq!(format_ratio(w), "1.58980044346");
        assert_eq!(format_ratio(0.00009499999999999999), "9.5e-05");
        assert_eq!(format_ratio(1_234_567_890_123.0), "1.23456789012e+12");
        assert_eq!(format_ratio(1.713487071977638), "1.71348707198");
        assert_eq!(format_ratio(2.1964285714285716), "2.19642857143");
    }

    #[test]
    fn python_two_ratio_notation_uses_the_rounded_exponent() {
        assert_eq!(format_ratio(99_999_999_999.0), "99999999999.0");
        assert_eq!(format_ratio(100_000_000_000.0), "1e+11");
        assert_eq!(format_ratio(9.99999999999e-5), "9.99999999999e-05");
        assert_eq!(format_ratio(9.999999999999e-5), "0.0001");
    }

    #[test]
    fn full_repr_notation_uses_adjacent_numeric_thresholds() {
        let below_upper = f64::from_bits(1e16_f64.to_bits() - 1);
        let below_lower = f64::from_bits(1e-4_f64.to_bits() - 1);

        assert_eq!(full_repr_float(below_upper), "9999999999999998.0");
        assert_eq!(full_repr_float(1e16), "1e+16");
        assert_eq!(full_repr_float(below_lower), "9.999999999999999e-05");
        assert_eq!(full_repr_float(1e-4), "0.0001");
    }

    #[test]
    fn pandas_xstrtod_matches_legacy_reference() {
        use super::pandas_xstrtod;
        // (input, expected pandas-parsed f64). These were captured by
        // running `pandas.to_numeric` inside the legacy Wave container
        // (pandas 0.19.2 / Python 2.7) — the exact pipeline that drives
        // legacy METRIC.* CSV emission.
        let cases: &[(&str, f64)] = &[
            ("0.877140", f64::from_bits(0x3fec1187e7c06e19)), // predecessor of 0.87714
            ("0.958635", f64::from_bits(0x3feead234eb9a177)), // exact
            ("0.916079", f64::from_bits(0x3fed5084e831ad22)), // successor of "0.916079"
            ("0.844803", f64::from_bits(0x3feb08a04d12018b)), // successor
            ("0.020633", f64::from_bits(0x3f9520d130df9bdd)), // successor
            ("0.196971", f64::from_bits(0x3fc9365881a15550)), // exact
            ("0.298002", f64::from_bits(0x3fd31276fb09203a)), // exact
            ("0.992971", f64::from_bits(0x3fefc66b1e5c0b99)), // predecessor
            ("0.90076", f64::from_bits(0x3fecd306a2b17050)),  // exact
        ];
        for (s, expected) in cases {
            let got = pandas_xstrtod(s);
            assert_eq!(
                got.to_bits(),
                expected.to_bits(),
                "pandas_xstrtod({s}) = {got} (hex {:x}), expected {expected} (hex {:x})",
                got.to_bits(),
                expected.to_bits()
            );
        }
    }

    #[test]
    fn format_metric_matches_legacy_byte_for_byte() {
        use super::format_metric;
        assert_eq!(format_metric(7839.0 / 8937.0), "0.87714");
        // Precision = 7949/8292 → "0.958635" → pandas-xstrtod → "0.958635"
        assert_eq!(format_metric(7949.0 / 8292.0), "0.958635");
        // F1 ≈ 0.916078519... → C++ "0.916079" → Python 2 "0.916079".
        let r = 7839.0 / 8937.0;
        let p = 7949.0 / 8292.0;
        let f1 = 2.0 * r * p / (r + p);
        assert_eq!(format_metric(f1), "0.916079");
        // Frac_NA = 3520/11812 → "0.298002" → exact f64 → "0.298002"
        assert_eq!(format_metric(3520.0 / 11812.0), "0.298002");
        assert_eq!(format_metric(0.976271), "0.976271");
    }

    #[test]
    fn summary_omits_zero_observation_unknown_type() {
        let output = tempfile::NamedTempFile::new().expect("temporary summary output");
        let mut indel = TypeCounts::default();
        indel.truth_total.total = 1;
        let counts = BTreeMap::from([
            ("INDEL".to_string(), indel),
            ("UNK".to_string(), TypeCounts::default()),
        ]);
        write_summary(
            output.path(),
            &counts,
            &counts,
            &BTreeMap::new(),
            &BTreeMap::new(),
        )
        .expect("write summary");
        let summary = std::fs::read_to_string(output.path()).expect("read summary");
        assert!(summary.lines().any(|line| line.starts_with("INDEL,ALL,")));
        assert!(!summary.lines().any(|line| line.starts_with("UNK,")));
    }

    #[test]
    fn extended_omits_indel_subtype_rows_for_subsets_without_indels() {
        let output = tempfile::NamedTempFile::new().expect("temporary extended output");
        let counts = BTreeMap::from([
            ("INDEL".to_string(), TypeCounts::default()),
            ("SNP".to_string(), TypeCounts::default()),
        ]);
        let subset_counts = BTreeMap::from([
            (
                "TS_boundary".to_string(),
                BTreeMap::from([("SNP".to_string(), TypeCounts::default())]),
            ),
            (
                "TS_contained".to_string(),
                BTreeMap::from([("INDEL".to_string(), TypeCounts::default())]),
            ),
        ]);
        let nested_counts = BTreeMap::<String, BTreeMap<String, TypeCounts>>::new();
        let subset_subtype_counts =
            BTreeMap::<String, BTreeMap<String, BTreeMap<String, TypeCounts>>>::new();
        let fp_counts = BTreeMap::<String, (usize, usize)>::new();
        let nested_fp_counts = BTreeMap::<String, BTreeMap<String, (usize, usize)>>::new();
        let subset_subtype_fp_counts = SubsetSubtypeFpCounts::new();

        write_extended(
            output.path(),
            &counts,
            &counts,
            &nested_counts,
            &nested_counts,
            100,
            100,
            80,
            true,
            &subset_counts,
            &subset_counts,
            &subset_subtype_counts,
            &subset_subtype_counts,
            &fp_counts,
            &fp_counts,
            &nested_fp_counts,
            &nested_fp_counts,
            &nested_fp_counts,
            &nested_fp_counts,
            &subset_subtype_fp_counts,
            &subset_subtype_fp_counts,
        )
        .expect("write extended report");

        let output = std::fs::read_to_string(output.path()).expect("read extended report");
        let indel_subtype_rows = output
            .lines()
            .skip(1)
            .filter(|line| line.starts_with("INDEL,") && !line.starts_with("INDEL,*,"))
            .collect::<Vec<_>>();
        assert_eq!(
            indel_subtype_rows
                .iter()
                .filter(|line| line.split(',').nth(2) == Some("TS_contained"))
                .count(),
            18,
            "present INDEL subsets retain their seeded subtype/filter rows"
        );
        assert!(
            indel_subtype_rows
                .iter()
                .all(|line| line.split(',').nth(2) != Some("TS_boundary")),
            "a SNP-only subset must not seed empty INDEL subtype rows"
        );
    }
}
