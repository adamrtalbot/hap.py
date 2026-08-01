//! Legacy MuTect, Pisces, and VarScan2 feature-table extractors.

use std::collections::BTreeMap;

use crate::strelka;
use crate::vcf::RawVcfRecord;

use super::common::{
    csv_escape, format_info_float, format_python_float, parse_scoring_features, render_filter,
};

pub(super) fn emit_pisces_with_depths(
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depth_override: Option<&BTreeMap<String, f64>>,
) -> Vec<String> {
    let scoring = parse_scoring_features(headers, "snv_scoring_features");
    let header_depths = strelka::parse_depths(headers);
    let depths = depth_override.unwrap_or(&header_depths);
    let mut header = String::from(",CHROM,POS,REF,ALT,FILTER,GQX,EVS,T_DP,T_DP_RATE,T_AF,tag");
    for feature in &scoring {
        header.push_str(",E.");
        header.push_str(feature);
    }
    let mut lines = vec![header];
    for (index, record) in records.iter().enumerate() {
        let tumor = record.sample_map(0);
        let t_dp = number(tumor.get("DP"), 0.0);
        let cells = [
            index.to_string(),
            csv_escape(&record.chrom),
            record.pos.to_string(),
            csv_escape(&record.ref_allele),
            csv_escape(&record.alt_allele),
            csv_escape(&render_filter(&record.filter)),
            format_python_float(number(tumor.get("GQX"), 0.0)),
            format_info_float(strelka::info_float(&record.info, "EVS").unwrap_or(-1.0)),
            format_python_float(t_dp),
            format_python_float(
                depths
                    .get(&record.chrom)
                    .filter(|depth| **depth != 0.0)
                    .map(|depth| t_dp / depth)
                    .unwrap_or(0.0),
            ),
            format_python_float(number(tumor.get("VF"), 0.0)),
            csv_escape(label),
        ];
        let mut line = cells.join(",");
        // Legacy allocates scoring-feature columns but never populates them.
        line.push_str(&",".repeat(scoring.len()));
        lines.push(line);
    }
    lines
}

pub(super) fn emit_varscan_with_depths(
    records: &[RawVcfRecord],
    label: &str,
    indel: bool,
    depths: Option<&BTreeMap<String, f64>>,
) -> Vec<String> {
    const HEADER: &str = ",CHROM,POS,REF,ALT,FILTER,SSC,GPV,SPV,N_DP,T_DP,N_DP_RATE,T_DP_RATE,N_GT,T_GT,N_GQ,T_GQ,N_AD,T_AD,N_FA,T_FA,N_ALT_RATE,T_ALT_RATE,tag";
    let depth_rates_are_float = records.iter().any(|record| {
        depths.is_some_and(|depths| depths.get(&record.chrom).is_some_and(|depth| *depth != 0.0))
    });
    let mut lines = vec![HEADER.to_string()];
    for (index, record) in records.iter().enumerate() {
        let normal = record.sample_map(0);
        let tumor = record.sample_map(1);
        let n_dp = int_sample(&normal, "DP", 0) as f64;
        let t_dp = int_sample(&tumor, "DP", 0) as f64;
        let n_rd = int_sample(&normal, "RD", 0) as f64;
        let t_rd = int_sample(&tumor, "RD", 0) as f64;
        let n_ad = int_sample(&normal, "AD", 0) as f64;
        let t_ad = int_sample(&tumor, "AD", 0) as f64;
        let n_rate = ratio(n_ad, n_rd);
        let t_rate = ratio(t_ad, t_rd);
        let n_fa = varscan_frequency(normal.get("FREQ"), indel);
        let t_fa = varscan_frequency(tumor.get("FREQ"), indel);
        let depth = depths.and_then(|depths| depths.get(&record.chrom)).copied();
        lines.push(
            [
                index.to_string(),
                csv_escape(&record.chrom),
                record.pos.to_string(),
                csv_escape(&record.ref_allele),
                csv_escape(&record.alt_allele),
                csv_escape(&render_filter(&record.filter)),
                info_scalar(&record.info, "SSC"),
                info_scalar(&record.info, "GPV"),
                info_scalar(&record.info, "SPV"),
                format_python_float(n_dp),
                format_python_float(t_dp),
                format_depth_rate(n_dp, depth, depth_rates_are_float),
                format_depth_rate(t_dp, depth, depth_rates_are_float),
                int_sample(&normal, "GT", 0).to_string(),
                int_sample(&tumor, "GT", 0).to_string(),
                int_sample(&normal, "GQ", 0).to_string(),
                int_sample(&tumor, "GQ", 0).to_string(),
                int_sample(&normal, "AD", 0).to_string(),
                int_sample(&tumor, "AD", 0).to_string(),
                n_fa,
                t_fa,
                format_python_float(n_rate),
                format_python_float(t_rate),
                csv_escape(label),
            ]
            .join(","),
        );
    }
    lines
}

pub(super) fn emit_mutect_with_depths(
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depths: Option<&BTreeMap<String, f64>>,
) -> Vec<String> {
    const HEADER: &str = ",CHROM,POS,REF,ALT,FILTER,TLOD,NLOD,DBSNP,N_DP,T_DP,N_DP_RATE,T_DP_RATE,N_GT,T_GT,N_AD,T_AD,N_QSS,T_QSS,N_AF,T_AF,tag";
    let (tumor_index, normal_index) = mutect_sample_indices(headers);
    let depth_rates_are_float = records.iter().any(|record| {
        depths.is_some_and(|depths| depths.get(&record.chrom).is_some_and(|depth| *depth != 0.0))
    });
    let mut lines = vec![HEADER.to_string()];
    for (index, record) in records.iter().enumerate() {
        let normal = record.sample_map(normal_index);
        let tumor = record.sample_map(tumor_index);
        let n_dp = coerced_int(normal.get("DP"), 0) as f64;
        let t_dp = coerced_int(tumor.get("DP"), 0) as f64;
        let n_ad = numeric_list(normal.get("AD"), record.alt_allele.split(',').count() + 1);
        let t_ad = numeric_list(tumor.get("AD"), record.alt_allele.split(',').count() + 1);
        let n_qss = numeric_list(normal.get("QSS"), record.alt_allele.split(',').count() + 1);
        let t_qss = numeric_list(tumor.get("QSS"), record.alt_allele.split(',').count() + 1);
        let depth = depths.and_then(|depths| depths.get(&record.chrom)).copied();
        lines.push(
            [
                index.to_string(),
                csv_escape(&record.chrom),
                record.pos.to_string(),
                csv_escape(&record.ref_allele),
                csv_escape(&record.alt_allele),
                csv_escape(&render_filter(&record.filter)),
                format_python_float(coerced_info_int(&record.info, "TLOD", 0) as f64),
                format_python_float(coerced_info_int(&record.info, "NLOD", 0) as f64),
                coerced_info_int(&record.info, "DB", 0).to_string(),
                format_python_float(n_dp),
                format_python_float(t_dp),
                format_depth_rate(n_dp, depth, depth_rates_are_float),
                format_depth_rate(t_dp, depth, depth_rates_are_float),
                coerced_int(normal.get("GT"), -1).to_string(),
                coerced_int(tumor.get("GT"), -1).to_string(),
                render_python_list(&n_ad),
                render_python_list(&t_ad),
                render_python_list(&n_qss),
                render_python_list(&t_qss),
                format_python_float(allele_fraction(&n_ad)),
                format_python_float(allele_fraction(&t_ad)),
                csv_escape(label),
            ]
            .join(","),
        );
    }
    lines
}

fn mutect_sample_indices(headers: &[String]) -> (usize, usize) {
    let samples: Vec<&str> = headers
        .iter()
        .find(|line| line.starts_with("#CHROM\t"))
        .map(|line| line.split('\t').skip(9).collect())
        .unwrap_or_default();
    let command = headers.iter().find(|line| {
        let lower = line.to_ascii_lowercase();
        lower.contains("gatkcommandline") && lower.contains("id=mutect")
    });
    let find_option = |name: &str| -> Option<usize> {
        let command = command?;
        let start = command.find(name)? + name.len();
        let value = command[start..]
            .split(|ch: char| ch.is_whitespace() || matches!(ch, '"' | '>' | ','))
            .next()?;
        samples.iter().position(|sample| *sample == value)
    };
    (
        find_option("tumor_sample_name=").unwrap_or(0),
        find_option("normal_sample_name=").unwrap_or(1),
    )
}

fn info_scalar(info: &str, key: &str) -> String {
    strelka::info_value(info, key)
        .map(|value| {
            if let Ok(value) = value.parse::<i64>() {
                value.to_string()
            } else if let Ok(value) = value.parse::<f64>() {
                format_python_float(value)
            } else {
                csv_escape(&value)
            }
        })
        .unwrap_or_default()
}

fn int_sample(sample: &std::collections::BTreeMap<String, String>, key: &str, default: i64) -> i64 {
    coerced_int(sample.get(key), default)
}

fn coerced_int(value: Option<&String>, default: i64) -> i64 {
    match value {
        None => default,
        Some(value) if value == "True" => 1,
        Some(value) if value == "False" => 0,
        Some(value) => value
            .parse::<i64>()
            .or_else(|_| value.parse::<f64>().map(|number| number as i64))
            .unwrap_or(-1),
    }
}

fn coerced_info_int(info: &str, key: &str, default: i64) -> i64 {
    let value = strelka::info_value(info, key);
    coerced_int(value.as_ref(), default)
}

fn number(value: Option<&String>, default: f64) -> f64 {
    value
        .and_then(|value| value.parse().ok())
        .unwrap_or(default)
}

fn ratio(alt: f64, reference: f64) -> f64 {
    if alt + reference == 0.0 {
        0.0
    } else {
        alt / (alt + reference)
    }
}

fn depth_rate(depth: f64, normalization: Option<f64>) -> f64 {
    normalization
        .filter(|normalization| *normalization != 0.0)
        .map(|normalization| depth / normalization)
        .unwrap_or(0.0)
}

fn format_depth_rate(depth: f64, normalization: Option<f64>, float_column: bool) -> String {
    if float_column {
        format_python_float(depth_rate(depth, normalization))
    } else {
        // Legacy initializes these fields with integer zero. Pandas retains
        // that integer dtype when no row was normalized, so `to_csv` writes
        // `0` rather than `0.0`.
        "0".to_string()
    }
}

fn varscan_frequency(value: Option<&String>, indel: bool) -> String {
    let Some(value) = value else {
        return "0".to_string();
    };
    if indel {
        coerced_int(Some(value), 0).to_string()
    } else {
        value
            .parse::<f64>()
            .map(format_python_float)
            .unwrap_or_default()
    }
}

fn numeric_list(value: Option<&String>, expected: usize) -> Vec<f64> {
    let Some(value) = value.filter(|value| value.contains(',')) else {
        return vec![0.0; expected];
    };
    value
        .split(',')
        .map(|part| part.parse::<f64>().unwrap_or(0.0))
        .collect()
}

fn render_python_list(values: &[f64]) -> String {
    let rendered = values
        .iter()
        .map(|value| {
            if value.fract() == 0.0 {
                (*value as i64).to_string()
            } else {
                format_python_float(*value)
            }
        })
        .collect::<Vec<_>>()
        .join(", ");
    csv_escape(&format!("[{rendered}]"))
}

fn allele_fraction(values: &[f64]) -> f64 {
    let reference = values.first().copied().unwrap_or(0.0);
    let alt: f64 = values.iter().skip(1).sum();
    ratio(alt, reference)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(info: &str, format: &str, samples: &[&str]) -> RawVcfRecord {
        RawVcfRecord {
            chrom: "chr1".to_string(),
            pos: 9,
            id: ".".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "C".to_string(),
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            info: info.to_string(),
            format: Some(format.to_string()),
            samples: samples.iter().map(|v| v.to_string()).collect(),
        }
    }

    #[test]
    fn pisces_scoring_columns_are_intentionally_blank() {
        let lines = emit_pisces_with_depths(
            &[record("EVS=2", "DP:VF:GQX", &["20:0.25:30"])],
            &["##snv_scoring_features=x,y".to_string()],
            "FP",
            None,
        );
        assert!(lines[1].ends_with(",FP,,"));
    }

    #[test]
    fn varscan_matches_legacy_gt_and_frequency_coercion() {
        let lines = emit_varscan_with_depths(
            &[record(
                "SSC=12",
                "GT:GQ:DP:RD:AD:FREQ",
                &["0/0:30:10:9:1:10%", "0/1:40:20:5:15:75%"],
            )],
            "FP",
            false,
            None,
        );
        let cells: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(&cells[13..15], &["-1", "-1"]);
        assert_eq!(&cells[19..21], &["", ""]);
        assert_eq!(cells[21], "0.1");
        assert_eq!(cells[22], "0.75");
        assert_eq!(&cells[11..13], &["0", "0"]);
    }

    #[test]
    fn mutect_defaults_to_tumor_then_normal_and_renders_lists() {
        let lines = emit_mutect_with_depths(
            &[record(
                "DB;TLOD=12.7;NLOD=3.2",
                "GT:DP:AD:QSS",
                &["0/1:20:5,15:1,2", "0/0:10:9,1:3,4"],
            )],
            &[],
            "FP",
            None,
        );
        assert!(lines[1].contains(",12.0,3.0,1,10.0,20.0,"));
        assert!(lines[1].contains("\"[9, 1]\",\"[5, 15]\""));
        let cells: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(&cells[11..13], &["0", "0"]);
    }
}
