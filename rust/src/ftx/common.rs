//! Shared helpers for the `hap ftx` feature-table extractors.
//!
//! Every extractor emits a CSV whose float cells must byte-match what
//! pandas' `DataFrame.to_csv(...)` produces for the legacy Python
//! extractor. `format_python_float` is the single entry point for that —
//! other feature tables should never reach for `format!("{v}")` directly.

use crate::report;

/// Escapes a CSV cell only when the content demands it. Mirrors the CSV
/// writer used by pandas: delimiters, quotes, and embedded line endings all
/// require quoting, with quotes doubled inside the cell.
pub(super) fn csv_escape(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

/// FILTER column rendering matches legacy `vcfextract.vcfExtract` semantics:
/// `PASS` and `.` both resolve to an empty list which pandas then serialises as
/// an empty cell; anything else is passed through comma-joined (multi-filter
/// records already use that separator in the VCF field).
pub(super) fn render_filter(value: &str) -> String {
    if value == "PASS" || value == "." || value.is_empty() {
        String::new()
    } else {
        value.to_string()
    }
}

/// Formats an `f64` the way pandas `DataFrame.to_csv` does for a float
/// column: integer-valued floats keep a trailing `.0`, non-integer floats
/// use the shortest round-trippable representation, and `NaN` renders as
/// an empty cell. Infinities survive as `inf` / `-inf` because that's
/// Python `repr()`'s choice for the same values.
pub(super) fn format_python_float(value: f64) -> String {
    if value.is_nan() {
        return String::new();
    }
    report::python_repr_float(value)
}

/// Numeric INFO values pass through htslib during legacy preprocessing,
/// which canonicalises negative zero to positive zero before pandas sees it.
pub(super) fn format_info_float(value: f64) -> String {
    format_python_float(if value == 0.0 { 0.0 } else { value })
}

/// Returns the ordered list of scoring-feature names from a Strelka VCF
/// header. The line looks like `##snv_scoring_features=GQX,EVS_LOG,...`
/// (Strelka's VQSR output); absent → empty Vec.
///
/// Whitespace around names is stripped to match the Python extractor's
/// tolerance of hand-edited inputs. The same helper serves indel scoring
/// features (`##indel_scoring_features=`) via the `prefix` parameter.
pub(super) fn parse_scoring_features(headers: &[String], prefix: &str) -> Vec<String> {
    let needle = format!("##{prefix}=");
    for header in headers {
        if let Some(rest) = header.strip_prefix(&needle) {
            return rest
                .split(',')
                .map(|n| n.trim().to_string())
                .filter(|n| !n.is_empty())
                .collect();
        }
    }
    Vec::new()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_python_float_matches_pandas_repr() {
        assert_eq!(format_python_float(26.0), "26.0");
        assert_eq!(format_python_float(0.0), "0.0");
        assert_eq!(format_python_float(-1.0), "-1.0");
        assert_eq!(format_python_float(0.5), "0.5");
        assert_eq!(format_python_float(0.1), "0.1");
        assert_eq!(format_python_float(26.0 / 135.5), "0.1918819188191882");
        assert_eq!(format_python_float(77.0 / 135.5), "0.5682656826568265");
    }

    #[test]
    fn format_python_float_nan_is_empty_cell() {
        assert_eq!(format_python_float(f64::NAN), "");
    }

    #[test]
    fn format_python_float_infinite_renders_python_repr() {
        assert_eq!(format_python_float(f64::INFINITY), "inf");
        assert_eq!(format_python_float(f64::NEG_INFINITY), "-inf");
    }

    #[test]
    fn info_float_collapses_negative_zero_like_htslib() {
        assert_eq!(format_python_float(-0.0), "-0.0");
        assert_eq!(format_info_float(-0.0), "0.0");
    }

    #[test]
    fn parse_scoring_features_absent_returns_empty() {
        let headers = vec!["##fileformat=VCFv4.1".to_string()];
        assert!(parse_scoring_features(&headers, "snv_scoring_features").is_empty());
    }

    #[test]
    fn parse_scoring_features_splits_on_comma() {
        let headers = vec![
            "##fileformat=VCFv4.1".to_string(),
            "##snv_scoring_features=GQX,EVS_LOG,AD_RATE".to_string(),
        ];
        assert_eq!(
            parse_scoring_features(&headers, "snv_scoring_features"),
            vec![
                "GQX".to_string(),
                "EVS_LOG".to_string(),
                "AD_RATE".to_string()
            ]
        );
    }

    #[test]
    fn parse_scoring_features_trims_whitespace() {
        let headers = vec!["##snv_scoring_features= GQX , EVS_LOG".to_string()];
        assert_eq!(
            parse_scoring_features(&headers, "snv_scoring_features"),
            vec!["GQX".to_string(), "EVS_LOG".to_string()]
        );
    }

    #[test]
    fn csv_escape_quotes_on_comma_and_double_quote() {
        assert_eq!(csv_escape("abc"), "abc");
        assert_eq!(csv_escape("a,b"), "\"a,b\"");
        assert_eq!(csv_escape("a\"b"), "\"a\"\"b\"");
    }

    #[test]
    fn csv_escape_quotes_embedded_line_endings() {
        assert_eq!(csv_escape("first\nsecond"), "\"first\nsecond\"");
        assert_eq!(csv_escape("first\rsecond"), "\"first\rsecond\"");
    }

    #[test]
    fn render_filter_collapses_pass_and_missing() {
        assert_eq!(render_filter("PASS"), "");
        assert_eq!(render_filter("."), "");
        assert_eq!(render_filter(""), "");
        assert_eq!(render_filter("LowQscore"), "LowQscore");
        assert_eq!(render_filter("FAIL1,FAIL2"), "FAIL1,FAIL2");
    }
}
