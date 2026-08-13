//! Legacy MuTect, Pisces, and VarScan2 feature-table extractors.

use std::collections::{BTreeMap, BTreeSet};

use anyhow::{Result, bail};

use crate::domain::RawVcfRecord;
use crate::engines::strelka;

use super::common::{
    csv_escape, format_info_float, format_python_float, parse_scoring_features, render_filter,
};

pub(super) fn emit_pisces_with_depths(
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depth_override: Option<&BTreeMap<String, f64>>,
) -> Result<Vec<String>> {
    let scoring = parse_scoring_features(headers, "snv_scoring_features");
    let header_depths = strelka::parse_depths(headers);
    let depths = depth_override.unwrap_or(&header_depths);
    let mut header = String::from(",CHROM,POS,REF,ALT,FILTER,GQX,EVS,T_DP,T_DP_RATE,T_AF,tag");
    for feature in &scoring.columns {
        header.push_str(",E.");
        header.push_str(feature);
    }
    let mut lines = vec![header];
    for (index, record) in records.iter().enumerate() {
        let tumor = record.sample_map(0);
        let t_dp = required_number(tumor.get("DP"), "Pisces DP")?;
        let gqx = required_number(tumor.get("GQX"), "Pisces GQX")?;
        let vf = required_number(tumor.get("VF"), "Pisces VF")?;
        let cells = [
            index.to_string(),
            csv_escape(&record.chrom),
            record.pos.to_string(),
            csv_escape(&record.ref_allele),
            csv_escape(&record.alt_allele),
            csv_escape(&render_filter(&record.filter)),
            format_python_float(gqx),
            format_info_float(strelka::info_float(&record.info, "EVS").unwrap_or(-1.0)),
            format_python_float(t_dp),
            format_python_float(
                depths
                    .get(&record.chrom)
                    .filter(|depth| **depth != 0.0)
                    .map(|depth| t_dp / depth)
                    .unwrap_or(0.0),
            ),
            format_python_float(vf),
            csv_escape(label),
        ];
        let mut line = cells.join(",");
        // Legacy allocates scoring-feature columns but never populates them.
        line.push_str(&",".repeat(scoring.columns.len()));
        lines.push(line);
    }
    Ok(lines)
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
) -> Result<Vec<String>> {
    const HEADER: &str = ",CHROM,POS,REF,ALT,FILTER,TLOD,NLOD,DBSNP,N_DP,T_DP,N_DP_RATE,T_DP_RATE,N_GT,T_GT,N_AD,T_AD,N_QSS,T_QSS,N_AF,T_AF,tag";
    let (tumor_index, normal_index) = mutect_sample_indices(headers);
    let depth_rates_are_float = records.iter().any(|record| {
        depths.is_some_and(|depths| depths.get(&record.chrom).is_some_and(|depth| *depth != 0.0))
    });
    let mut lines = vec![HEADER.to_string()];
    let mut scalar_list_warnings = BTreeSet::new();
    for (index, record) in records.iter().enumerate() {
        let normal = record.sample_map(normal_index);
        let tumor = record.sample_map(tumor_index);
        let n_dp = coerced_int(normal.get("DP"), 0) as f64;
        let t_dp = coerced_int(tumor.get("DP"), 0) as f64;
        let alt_count = record.alt_allele.split(',').count();
        let expected = alt_count + 1;
        let mut record_lists = BTreeMap::new();
        let n_qss = mutect_list_value(
            &mut record_lists,
            &mut scalar_list_warnings,
            normal_index,
            MutectListField::Qss,
            normal.get("QSS"),
            expected,
        )?;
        let t_qss = mutect_list_value(
            &mut record_lists,
            &mut scalar_list_warnings,
            tumor_index,
            MutectListField::Qss,
            tumor.get("QSS"),
            expected,
        )?;
        let n_ad = mutect_list_value(
            &mut record_lists,
            &mut scalar_list_warnings,
            normal_index,
            MutectListField::Ad,
            normal.get("AD"),
            expected,
        )?
        .into_list();
        let t_ad = mutect_list_value(
            &mut record_lists,
            &mut scalar_list_warnings,
            tumor_index,
            MutectListField::Ad,
            tumor.get("AD"),
            expected,
        )?
        .into_list();
        let required_ad_len = if record.alt_allele == "." {
            1
        } else {
            expected
        };
        if n_ad.len() < required_ad_len || t_ad.len() < required_ad_len {
            bail!("MuTect AD has fewer values than REF plus ALT alleles");
        }
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
                n_qss.render(),
                t_qss.render(),
                format_python_float(allele_fraction(&n_ad, &record.alt_allele)),
                format_python_float(allele_fraction(&t_ad, &record.alt_allele)),
                csv_escape(label),
            ]
            .join(","),
        );
    }
    Ok(lines)
}

fn mutect_sample_indices(headers: &[String]) -> (usize, usize) {
    let samples: Vec<&str> = headers
        .iter()
        .find(|line| line.starts_with("#CHROM\t"))
        .map(|line| line.split('\t').skip(9).collect())
        .unwrap_or_default();
    let find_option = |command: &str, name: &str| -> Option<usize> {
        let start = command.find(name)? + name.len();
        // The Python regex captures every non-whitespace character. In
        // particular, it retains the closing quote when an option is last in
        // CommandLineOptions, causing that sample lookup to fall back to its
        // legacy default. Do not trim VCF-header delimiters here.
        let value = command[start..].split_whitespace().next()?;
        samples.iter().position(|sample| *sample == value)
    };
    let mut tumor_index = 0;
    let mut normal_index = 1;
    for command in headers.iter().filter(|line| {
        let lower = line.to_ascii_lowercase();
        lower.contains("gatkcommandline") && lower.contains("id=mutect")
    }) {
        if let Some(index) = find_option(command, "tumor_sample_name=") {
            tumor_index = index;
        }
        if let Some(index) = find_option(command, "normal_sample_name=") {
            normal_index = index;
        }
    }
    (tumor_index, normal_index)
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

fn required_number(value: Option<&String>, field: &str) -> Result<f64> {
    let Some(value) = value else {
        bail!("missing required {field}");
    };
    value
        .parse()
        .map_err(|_| anyhow::anyhow!("invalid required {field}: {value}"))
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

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum MutectListField {
    Ad,
    Qss,
}

impl MutectListField {
    fn name(self) -> &'static str {
        match self {
            Self::Ad => "AD",
            Self::Qss => "QSS",
        }
    }
}

#[derive(Clone, Debug)]
enum MutectListValue {
    ScalarZero,
    List(Vec<f64>),
}

impl MutectListValue {
    fn into_list(self) -> Vec<f64> {
        match self {
            Self::List(values) => values,
            Self::ScalarZero => unreachable!("AD is rejected before scalar output"),
        }
    }

    fn render(&self) -> String {
        match self {
            Self::ScalarZero => "0".to_string(),
            Self::List(values) => render_python_list(values),
        }
    }
}

fn mutect_list_value(
    record_cache: &mut BTreeMap<(usize, MutectListField), MutectListValue>,
    scalar_warnings: &mut BTreeSet<(usize, MutectListField)>,
    sample_index: usize,
    field: MutectListField,
    value: Option<&String>,
    expected: usize,
) -> Result<MutectListValue> {
    let key = (sample_index, field);
    if let Some(value) = record_cache.get(&key) {
        return Ok(value.clone());
    }

    let parsed = match value {
        None if field == MutectListField::Ad => {
            bail!("missing required MuTect {}", field.name())
        }
        None => MutectListValue::ScalarZero,
        Some(value) if value.contains(',') && field == MutectListField::Ad => {
            let values = value
                .split(',')
                .map(|part| {
                    part.parse::<f64>().map_err(|_| {
                        anyhow::anyhow!("invalid MuTect {} value: {part}", field.name())
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            MutectListValue::List(values)
        }
        Some(value) if value.contains(',') => {
            MutectListValue::List(numeric_list(Some(value), expected))
        }
        Some(_) if scalar_warnings.insert(key) => MutectListValue::List(vec![0.0; expected]),
        Some(_) => bail!(
            "repeated scalar MuTect {} for sample {}",
            field.name(),
            sample_index + 1
        ),
    };
    record_cache.insert(key, parsed.clone());
    Ok(parsed)
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

fn allele_fraction(values: &[f64], alt_allele: &str) -> f64 {
    let reference = values.first().copied().unwrap_or(0.0);
    let alt = if alt_allele == "." {
        0.0
    } else {
        values
            .iter()
            .skip(1)
            .take(alt_allele.split(',').count())
            .sum()
    };
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
            mixed_edit_primitive: false,
        }
    }

    #[test]
    fn pisces_scoring_columns_are_intentionally_blank() {
        let lines = emit_pisces_with_depths(
            &[record("EVS=2", "DP:VF:GQX", &["20:0.25:30"])],
            &["##snv_scoring_features=x,y".to_string()],
            "FP",
            None,
        )
        .unwrap();
        assert!(lines[1].ends_with(",FP,,"));
    }

    #[test]
    fn pisces_keeps_all_literal_scoring_columns_blank() {
        let lines = emit_pisces_with_depths(
            &[record("EVS=2", "DP:VF:GQX", &["20:0.25:30"])],
            &[
                "##snv_scoring_features= old ".to_string(),
                "##snv_scoring_features=new,".to_string(),
            ],
            "FP",
            None,
        )
        .unwrap();
        assert!(lines[0].ends_with(",E. old ,E.new,E."), "{}", lines[0]);
        assert!(lines[1].ends_with(",FP,,,"), "{}", lines[1]);
    }

    #[test]
    fn pisces_missing_or_malformed_required_numbers_are_rejected() {
        for sample in ["0.25:30", "bad:0.25:30"] {
            let format = if sample.starts_with("bad") {
                "DP:VF:GQX"
            } else {
                "VF:GQX"
            };
            let result =
                emit_pisces_with_depths(&[record("EVS=2", format, &[sample])], &[], "FP", None);
            assert!(result.is_err(), "legacy extraction must fail before a row");
        }
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
        )
        .unwrap();
        assert!(lines[1].contains(",12.0,3.0,1,10.0,20.0,"));
        assert!(lines[1].contains("\"[9, 1]\",\"[5, 15]\""));
        let cells: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(&cells[11..13], &["0", "0"]);
    }

    #[test]
    fn mutect_missing_ad_is_rejected() {
        let result = emit_mutect_with_depths(
            &[record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:QSS",
                &["0/1:20:1,2", "0/0:10:3,4"],
            )],
            &[],
            "FP",
            None,
        );
        assert!(result.is_err(), "legacy indexes the missing scalar AD");
    }

    #[test]
    fn mutect_af_ignores_ad_values_beyond_alt_cardinality() {
        let lines = emit_mutect_with_depths(
            &[record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:AD:QSS",
                &["0/1:20:5,15,80:1,2", "0/0:10:9,1,90:3,4"],
            )],
            &[],
            "FP",
            None,
        )
        .unwrap();
        assert!(lines[1].contains(",0.1,0.75,FP"), "{}", lines[1]);
    }

    #[test]
    fn mutect_missing_qss_stays_scalar_zero() {
        let lines = emit_mutect_with_depths(
            &[record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:AD",
                &["0/1:20:5,15", "0/0:10:9,1"],
            )],
            &[],
            "FP",
            None,
        )
        .unwrap();

        assert!(
            lines[1].contains("\"[9, 1]\",\"[5, 15]\",0,0,0.1,0.75,FP"),
            "{}",
            lines[1]
        );
    }

    #[test]
    fn mutect_second_scalar_ad_for_a_sample_is_rejected() {
        let records = [
            record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:AD:QSS",
                &["0/1:20:7:1,2", "0/0:10:9:3,4"],
            ),
            record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:AD:QSS",
                &["0/1:20:8:1,2", "0/0:10:10:3,4"],
            ),
        ];

        assert!(emit_mutect_with_depths(&records, &[], "FP", None).is_err());
    }

    #[test]
    fn mutect_second_scalar_qss_for_a_sample_is_rejected() {
        let records = [
            record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:AD:QSS",
                &["0/1:20:5,15:7", "0/0:10:9,1:9"],
            ),
            record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:AD:QSS",
                &["0/1:20:5,15:8", "0/0:10:9,1:10"],
            ),
        ];

        assert!(emit_mutect_with_depths(&records, &[], "FP", None).is_err());
    }

    #[test]
    fn mutect_duplicate_sample_selection_repairs_a_scalar_only_once() {
        let samples = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR";
        let trailing_normal = concat!(
            "##GATKCommandLine=<ID=MuTect,Version=1.1.7,",
            "CommandLineOptions=\"tumor_sample_name=TUMOR normal_sample_name=NORMAL\">"
        );
        let lines = emit_mutect_with_depths(
            &[record(
                "DB;TLOD=12;NLOD=3",
                "GT:DP:AD:QSS",
                &["0/0:10:9,1:3,4", "0/1:20:7:8"],
            )],
            &[trailing_normal.to_string(), samples.to_string()],
            "FP",
            None,
        )
        .unwrap();

        assert!(
            lines[1].contains("\"[0, 0]\",\"[0, 0]\",\"[0, 0]\",\"[0, 0]\",0.0,0.0,FP"),
            "{}",
            lines[1]
        );
    }

    #[test]
    fn later_mutect_command_line_headers_take_precedence() {
        let samples = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR";
        let headers = vec![
            concat!(
                "##GATKCommandLine=<ID=MuTect,CommandLineOptions=\"",
                "tumor_sample_name=NORMAL normal_sample_name=TUMOR emit_mode=VCF\">"
            )
            .to_string(),
            concat!(
                "##GATKCommandLine=<ID=MuTect,CommandLineOptions=\"",
                "tumor_sample_name=TUMOR normal_sample_name=NORMAL emit_mode=VCF\">"
            )
            .to_string(),
            samples.to_string(),
        ];

        assert_eq!(mutect_sample_indices(&headers), (1, 0));
    }

    #[test]
    fn mutect_header_sample_order_preserves_legacy_token_boundaries() {
        let samples = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR";
        let trailing_normal = concat!(
            "##GATKCommandLine=<ID=MuTect,Version=1.1.7,",
            "CommandLineOptions=\"tumor_sample_name=TUMOR normal_sample_name=NORMAL\">"
        );
        let trailing_headers = [trailing_normal.to_string(), samples.to_string()];
        assert_eq!(
            mutect_sample_indices(&trailing_headers),
            (1, 1),
            "legacy retains the closing quote on the final normal-sample token"
        );
        let lines = emit_mutect_with_depths(
            &[record(
                "DB;TLOD=17.4;NLOD=2.1",
                "GT:DP:AD:QSS",
                &["0/0:12:11,1:7,1", "0/1:24:6,18:2,8"],
            )],
            &trailing_headers,
            "reference",
            None,
        )
        .unwrap();
        assert!(
            lines[1].contains(
                ",24.0,24.0,0,0,-1,-1,\"[6, 18]\",\"[6, 18]\",\"[2, 8]\",\"[2, 8]\",0.75,0.75,reference"
            ),
            "both legacy prefixes resolve to the second sample: {}",
            lines[1]
        );

        let followed_normal = concat!(
            "##GATKCommandLine=<ID=MuTect,Version=1.1.7,",
            "CommandLineOptions=\"tumor_sample_name=TUMOR normal_sample_name=NORMAL emit_mode=VCF\">"
        );
        assert_eq!(
            mutect_sample_indices(&[followed_normal.to_string(), samples.to_string()]),
            (1, 0),
            "a following option terminates the normal-sample token with whitespace"
        );
    }
}
