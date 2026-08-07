//! Strelka-specific parsing helpers shared across subcommands.
//!
//! Both `hap somatic` (comparison / feature-table emission) and
//! `hap ftx` (standalone feature extraction) need to read Strelka's
//! per-chromosome depth headers, pull INFO/FORMAT values, and compute
//! tier-1 allele frequencies. Keeping the helpers in one module prevents
//! subtle drift between the two call sites.

use std::collections::BTreeMap;

/// Returns the INFO value for `key`, or `None` when absent.
///
/// Bare flag entries (`SOMATIC` with no `=value`) resolve to
/// `Some("True")` — matching the legacy `vcfextract.vcfExtract` behaviour
/// that stores Number=0 info entries as Python `True`.
pub(crate) fn info_value(info: &str, key: &str) -> Option<String> {
    for part in info.split(';') {
        if let Some((k, v)) = part.split_once('=') {
            if k == key {
                return Some(v.to_string());
            }
        } else if part == key {
            return Some("True".to_string());
        }
    }
    None
}

/// Returns the INFO value parsed as an `f64`, or `None` when absent or
/// unparseable.
pub(crate) fn info_float(info: &str, key: &str) -> Option<f64> {
    info_value(info, key).and_then(|v| v.parse::<f64>().ok())
}

/// Parses the first comma-separated field of a FORMAT value as an `f64`.
///
/// Strelka encodes tier-1/tier-2 counts as `"<t1>,<t2>"`; the extractors
/// usually only want the tier-1 number, so this helper short-circuits the
/// `.split(',').next()` + parse pattern.
pub(crate) fn parse_first_number(value: Option<&String>) -> Option<f64> {
    value
        .and_then(|v| v.split(',').next())
        .and_then(|v| v.parse::<f64>().ok())
}

/// Reads Strelka's per-chromosome depth headers.
///
/// Strelka writes one of `##maxdepth_<chrom>=N`, `##meandepth_<chrom>=N`,
/// or `##depth_<chrom>=N` for every contig it called on. The legacy
/// Python lowercases the header line before matching, so uppercase contig
/// names in headers still hit. The returned map is keyed by the contig
/// name as written in the header (case preserved for the actual chrom).
pub(crate) fn parse_depths(headers: &[String]) -> BTreeMap<String, f64> {
    let mut out = BTreeMap::new();
    for header in headers {
        let lower = header.to_lowercase();
        if !(lower.starts_with("##maxdepth_")
            || lower.starts_with("##meandepth_")
            || lower.starts_with("##depth_"))
        {
            continue;
        }
        if let Some((left, right)) = header.trim_start_matches('#').split_once('=')
            && let Some((_, chrom)) = left.split_once('_')
            && let Ok(v) = right.parse::<f64>()
        {
            out.insert(chrom.to_string(), v);
        }
    }
    out
}

/// Computes tier-1 allele frequency from Strelka TIR/TAR-style counts.
///
/// Returns 0.0 when both counts are zero (matches legacy `if denom == 0`
/// short-circuit). Used by both the Strelka indel row (`TIR`/`TAR`) and
/// any future extractor that maps to a similar tier-1 fraction.
pub(crate) fn af_from_tir_tar(sample: &BTreeMap<String, String>, tir: &str, tar: &str) -> f64 {
    let alt = parse_first_number(sample.get(tir)).unwrap_or(0.0);
    let refc = parse_first_number(sample.get(tar)).unwrap_or(0.0);
    if alt + refc == 0.0 {
        0.0
    } else {
        alt / (alt + refc)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn info_value_finds_keyed_entry() {
        let info = "NT=ref;QSS_NT=47;SOMATIC";
        assert_eq!(info_value(info, "NT"), Some("ref".to_string()));
        assert_eq!(info_value(info, "QSS_NT"), Some("47".to_string()));
    }

    #[test]
    fn info_value_flags_become_true() {
        let info = "SOMATIC;NT=ref";
        assert_eq!(info_value(info, "SOMATIC"), Some("True".to_string()));
    }

    #[test]
    fn info_value_missing_returns_none() {
        assert_eq!(info_value("NT=ref", "QSS_NT"), None);
        assert_eq!(info_value("", "NT"), None);
    }

    #[test]
    fn info_float_parses_numeric() {
        assert_eq!(info_float("MQ=57.82;DP=112", "MQ"), Some(57.82));
        assert_eq!(info_float("MQ=57.82;DP=112", "MQ0"), None);
        assert_eq!(info_float("NT=ref", "NT"), None);
    }

    #[test]
    fn parse_first_number_handles_tier_tuples() {
        let v1 = Some("26,29".to_string());
        let v2 = Some("3".to_string());
        let v3 = Some("bad".to_string());
        assert_eq!(parse_first_number(v1.as_ref()), Some(26.0));
        assert_eq!(parse_first_number(v2.as_ref()), Some(3.0));
        assert_eq!(parse_first_number(v3.as_ref()), None);
        assert_eq!(parse_first_number(None), None);
    }

    #[test]
    fn parse_depths_accepts_all_three_prefixes() {
        let headers = vec![
            "##maxdepth_chr1=42.5".to_string(),
            "##meandepth_chr2=17.1".to_string(),
            "##depth_chrX=9".to_string(),
            "##unrelated=nope".to_string(),
        ];
        let depths = parse_depths(&headers);
        assert_eq!(depths.get("chr1"), Some(&42.5));
        assert_eq!(depths.get("chr2"), Some(&17.1));
        assert_eq!(depths.get("chrX"), Some(&9.0));
        assert!(!depths.contains_key("unrelated"));
    }

    #[test]
    fn af_from_tir_tar_zero_denominator_returns_zero() {
        let mut sample = BTreeMap::new();
        sample.insert("TIR".to_string(), "0,0".to_string());
        sample.insert("TAR".to_string(), "0,0".to_string());
        assert_eq!(af_from_tir_tar(&sample, "TIR", "TAR"), 0.0);
    }

    #[test]
    fn af_from_tir_tar_normal_ratio() {
        let mut sample = BTreeMap::new();
        sample.insert("TIR".to_string(), "3,4".to_string());
        sample.insert("TAR".to_string(), "7,9".to_string());
        assert_eq!(af_from_tir_tar(&sample, "TIR", "TAR"), 0.3);
    }
}
