use crate::compare::{AnnotatedRow, TypeCounts};
use crate::vcf;
use anyhow::Result;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

pub const EXTENDED_HEADER: [&str; 65] = [
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

#[derive(Clone, Debug, Default)]
pub struct CountsBucket {
    pub total: usize,
    pub ti: usize,
    pub tv: usize,
    pub het: usize,
    pub homalt: usize,
}

type SubsetSubtypeFpCounts = BTreeMap<String, BTreeMap<String, BTreeMap<String, (usize, usize)>>>;

pub fn write_summary(
    path: &Path,
    all_counts: &BTreeMap<String, TypeCounts>,
    pass_counts: &BTreeMap<String, TypeCounts>,
    all_fp: &BTreeMap<String, (usize, usize)>,
    pass_fp: &BTreeMap<String, (usize, usize)>,
) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,FP.al,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio"
    )?;
    let mut variant_types: Vec<&String> = all_counts.keys().chain(pass_counts.keys()).collect();
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
    Ok(())
}

#[allow(clippy::too_many_arguments)] // Each map is a distinct legacy report axis/tier.
pub fn write_extended(
    path: &Path,
    all_counts: &BTreeMap<String, TypeCounts>,
    pass_counts: &BTreeMap<String, TypeCounts>,
    all_subtype: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
    pass_subtype: &BTreeMap<String, BTreeMap<String, TypeCounts>>,
    subset_size: usize,
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
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(writer, "{}", EXTENDED_HEADER.join(","))?;

    let mut variant_types: Vec<&String> = all_counts.keys().chain(pass_counts.keys()).collect();
    variant_types.sort();
    variant_types.dedup();
    let empty = TypeCounts::default();

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
    let _ = &empty; // silence unused var warning in case of dead path

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
                //   • base subtype + base subset: always.
                //   • base subtype + non-base subset: only if the subset
                //     actually carries data for this Type.
                //   • non-base subtype + base subset: always (subtype rows
                //     are seeded for every INDEL_SUBTYPE).
                //   • non-base subtype + non-base subset: always
                //     (cross-product rows are seeded zero-counts).
                if is_base_subtype && !is_base_subset {
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
                    //   subset="TS_boundary"   → format_count(subset_size)
                    //   subset="TS_contained"  → format_count(conf_size)
                    let size_cell = match subset {
                        "*" => subset_size.to_string(),
                        "TS_contained" => format_count(conf_size),
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
    Ok(())
}

pub fn write_vcf(path: &Path, headers: &[String], rows: &[AnnotatedRow]) -> Result<()> {
    vcf::write_indexed_vcf(path, headers, rows.iter().map(|row| row.line.as_str()))
}

pub(crate) fn append_stats(row: &mut Vec<String>, stats: &CountsBucket, supports_titv: bool) {
    append_stats_with_missing(row, stats, supports_titv, ".");
}

fn append_extended_stats(
    row: &mut Vec<String>,
    stats: &CountsBucket,
    supports_titv: bool,
    missing_titv: &str,
) {
    append_stats_with_missing(row, stats, supports_titv, missing_titv);
}

fn append_stats_with_missing(
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
    // via Python's shortest_repr. Mirroring those three steps in Rust
    // reproduces the exact decimal text legacy emits — which is generally
    // 1 ULP off from a plain `(value * 1e6).round() / 1e6` round-trip.
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
fn pandas_xstrtod(s: &str) -> f64 {
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

fn format_ratio(value: f64) -> String {
    python_repr_float(value)
}

// Match Python's repr() for floats: shortest-roundtrip representation that
// always includes a decimal point. Rust's Debug formatter for f64 uses the
// same shortest-roundtrip algorithm as Python's repr (Grisu/Ryu) and
// preserves ".0" for integer-valued floats — exactly what pandas emits
// when DataFrame.to_csv renders a float column with the default format.
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
    let s = format!("{value:?}");
    // Python pandas zero-pads scientific notation exponents to ≥2 digits:
    // Rust:  "9.499999999999999e-5" / "1.0e10"
    // Python: "9.499999999999999e-05" / "1.0e+10"
    // Pandas always emits sign (e+xx or e-xx) with 2-digit minimum.
    if let Some(e_pos) = s.find('e') {
        let (mantissa, exponent) = s.split_at(e_pos);
        let exp_str = &exponent[1..]; // skip 'e'
        let (sign, digits) = if let Some(stripped) = exp_str.strip_prefix('-') {
            ("-", stripped)
        } else if let Some(stripped) = exp_str.strip_prefix('+') {
            ("+", stripped)
        } else {
            ("+", exp_str)
        };
        if digits.len() < 2 {
            return format!("{mantissa}e{sign}0{digits}");
        } else {
            return format!("{mantissa}e{sign}{digits}");
        }
    }
    s
}

#[cfg(test)]
mod format_tests {
    use super::python_repr_float;

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
    fn long_repr_preserves_all_digits() {
        // These values are chosen so their f64 shortest-repr has trailing 9s/0s
        // the way Python's repr does.
        let v: f64 = 0.9651790000000001;
        assert_eq!(python_repr_float(v), "0.9651790000000001");
        let w: f64 = 1.58980044345898;
        assert_eq!(python_repr_float(w), "1.58980044345898");
    }

    #[test]
    fn pandas_xstrtod_matches_legacy_oracle() {
        use super::pandas_xstrtod;
        // (input, expected pandas-parsed f64). These were captured by
        // running `pandas.to_numeric` inside the legacy Wave container
        // (pandas 0.24.2 / Python 2.7) — the exact pipeline that drives
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
        // Recall = 7839/8937 → C++ "0.877140" → pandas-xstrtod → "0.8771399999999999"
        assert_eq!(format_metric(7839.0 / 8937.0), "0.8771399999999999");
        // Precision = 7949/8292 → "0.958635" → pandas-xstrtod → "0.958635"
        assert_eq!(format_metric(7949.0 / 8292.0), "0.958635");
        // F1 ≈ 0.916078519... → "0.916079" → "0.9160790000000001"
        let r = 7839.0 / 8937.0;
        let p = 7949.0 / 8292.0;
        let f1 = 2.0 * r * p / (r + p);
        assert_eq!(format_metric(f1), "0.9160790000000001");
        // Frac_NA = 3520/11812 → "0.298002" → exact f64 → "0.298002"
        assert_eq!(format_metric(3520.0 / 11812.0), "0.298002");
    }
}
