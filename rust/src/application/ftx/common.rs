//! Shared helpers for the `hap ftx` feature-table extractors.
//!
//! Every extractor emits a CSV whose float cells must byte-match what
//! pandas' `DataFrame.to_csv(...)` produces for the legacy Python
//! extractor. `format_python_float` is the single entry point for that —
//! other feature tables should never reach for `format!("{v}")` directly.

use std::collections::BTreeMap;

use crate::adapters::report;

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
    report::full_repr_float(value)
}

/// Numeric INFO values pass through htslib during legacy preprocessing,
/// which canonicalises negative zero to positive zero before pandas sees it.
pub(super) fn format_info_float(value: f64) -> String {
    format_python_float(if value == 0.0 { 0.0 } else { value })
}

/// Legacy keeps every scoring-feature column from every matching header, but
/// its index-to-name lookup is overwritten by each later header. Those two
/// views intentionally differ when malformed VCFs repeat the metadata line.
#[derive(Debug, Default, Eq, PartialEq)]
pub(super) struct ScoringFeatures {
    pub(super) columns: Vec<String>,
    pub(super) names_by_index: Vec<String>,
}

impl ScoringFeatures {
    pub(super) fn render_evsf(&self, raw: Option<&str>) -> Vec<String> {
        let mut values = self
            .names_by_index
            .iter()
            .map(|name| (name.as_str(), 0.0))
            .collect::<BTreeMap<_, _>>();

        // `vcfExtract.field` produces a list only when the source contains a
        // comma. Enumerating a scalar raises in legacy and leaves defaults.
        if let Some(raw) = raw.filter(|value| value.contains(',')) {
            for (index, value) in raw.split(',').enumerate() {
                if let Some(name) = self.names_by_index.get(index)
                    && let Ok(value) = value.parse::<f64>()
                {
                    values.insert(name, value);
                }
            }
        }

        self.columns
            .iter()
            .map(|name| {
                values
                    .get(name.as_str())
                    .copied()
                    .map(format_info_float)
                    .unwrap_or_default()
            })
            .collect()
    }
}

/// Parse literal Strelka/Pisces scoring metadata using the legacy loop.
pub(super) fn parse_scoring_features(headers: &[String], prefix: &str) -> ScoringFeatures {
    let needle = format!("##{prefix}");
    let mut parsed = ScoringFeatures::default();
    for header in headers {
        if !header.contains(&needle) {
            continue;
        }
        let Some((_, value)) = header.split_once('=') else {
            continue;
        };
        for (index, name) in value.split(',').enumerate() {
            let name = name.to_string();
            parsed.columns.push(name.clone());
            if index < parsed.names_by_index.len() {
                parsed.names_by_index[index] = name;
            } else {
                parsed.names_by_index.push(name);
            }
        }
    }
    parsed
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
        assert_eq!(format_python_float(1.0 / 3.0), "0.3333333333333333");
    }

    #[test]
    fn format_python_float_matches_pandas_notation_threshold_bytes() {
        let below_upper = f64::from_bits(1e16_f64.to_bits() - 1);
        let below_lower = f64::from_bits(1e-4_f64.to_bits() - 1);

        assert_eq!(format_python_float(1e12), "1000000000000.0");
        assert_eq!(format_python_float(1e15), "1000000000000000.0");
        assert_eq!(format_python_float(below_upper), "9999999999999998.0");
        assert_eq!(format_python_float(1e16), "1e+16");
        assert_eq!(format_python_float(below_lower), "9.999999999999999e-05");
        assert_eq!(format_python_float(1e-4), "0.0001");
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
        assert!(
            parse_scoring_features(&headers, "snv_scoring_features")
                .columns
                .is_empty()
        );
    }

    #[test]
    fn parse_scoring_features_splits_on_comma() {
        let headers = vec![
            "##fileformat=VCFv4.1".to_string(),
            "##snv_scoring_features=GQX,EVS_LOG,AD_RATE".to_string(),
        ];
        assert_eq!(
            parse_scoring_features(&headers, "snv_scoring_features").columns,
            vec![
                "GQX".to_string(),
                "EVS_LOG".to_string(),
                "AD_RATE".to_string()
            ]
        );
    }

    #[test]
    fn parse_scoring_features_preserves_whitespace() {
        let headers = vec!["##snv_scoring_features= GQX , EVS_LOG".to_string()];
        assert_eq!(
            parse_scoring_features(&headers, "snv_scoring_features").columns,
            vec![" GQX ".to_string(), " EVS_LOG".to_string()]
        );
    }

    #[test]
    fn parse_scoring_features_appends_every_matching_header() {
        let headers = vec![
            "##snv_scoring_features=first,second".to_string(),
            "##snv_scoring_features=later".to_string(),
        ];
        let parsed = parse_scoring_features(&headers, "snv_scoring_features");
        assert_eq!(
            parsed.columns,
            vec![
                "first".to_string(),
                "second".to_string(),
                "later".to_string()
            ]
        );
        assert_eq!(
            parsed.names_by_index,
            vec!["later".to_string(), "second".to_string()]
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
