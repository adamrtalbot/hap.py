//! Strelka indel feature tables used by both admix and HCC modes.

use std::collections::BTreeMap;

use crate::domain::RawVcfRecord;
use crate::engines::strelka;

use super::common::{
    ScoringFeatures, csv_escape, format_info_float, format_python_float, parse_scoring_features,
    render_filter,
};

const FIXED_COLUMNS: &[&str] = &[
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "LENGTH",
    "INDELTYPE",
    "FILTER",
    "NT",
    "NT_REF",
    "EVS",
    "QSI_NT",
    "N_DP",
    "T_DP",
    "N_DP_RATE",
    "T_DP_RATE",
    "N_BCN",
    "T_BCN",
    "N_FDP",
    "T_FDP",
    "N_AF",
    "T_AF",
    "SGT",
    "RC",
    "RU",
    "RU_LEN",
    "IC",
    "IHP",
    "MQ",
    "MQ0",
    "tag",
];

pub(super) fn emit_with_depths(
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depth_override: Option<&BTreeMap<String, f64>>,
) -> Vec<String> {
    let scoring_features = parse_scoring_features(headers, "indel_scoring_features");
    let header_depths = strelka::parse_depths(headers);
    let depths = depth_override.unwrap_or(&header_depths);
    let mut header = format!(",{}", FIXED_COLUMNS.join(","));
    for feature in &scoring_features.columns {
        header.push_str(",E.");
        header.push_str(feature);
    }
    let mut lines = vec![header];
    lines.extend(
        records
            .iter()
            .enumerate()
            .map(|(index, record)| render_row(record, depths, &scoring_features, index, label)),
    );
    lines
}

fn render_row(
    record: &RawVcfRecord,
    depths: &BTreeMap<String, f64>,
    scoring_features: &ScoringFeatures,
    index: usize,
    label: &str,
) -> String {
    let normal = record.sample_map(0);
    let tumor = record.sample_map(1);
    let info = &record.info;
    let ref_len = record.ref_allele.len();
    let mut min_len = ref_len;
    let mut max_len = ref_len;
    let mut indel_type = 0;
    for alt in record.alt_allele.split(',') {
        if alt.len() > ref_len {
            indel_type |= 1;
        } else {
            indel_type |= 2;
        }
        min_len = min_len.min(alt.len());
        max_len = max_len.max(alt.len());
    }

    let nt = strelka::info_value(info, "NT").unwrap_or_default();
    let n_dp = first(&normal, "DP");
    let t_dp = first(&tumor, "DP");
    let depth = depths
        .get(&record.chrom)
        .copied()
        .filter(|depth| *depth != 0.0);
    let ru = strelka::info_value(info, "RU").unwrap_or_default();

    let mut cells = vec![
        index.to_string(),
        csv_escape(&record.chrom),
        record.pos.to_string(),
        csv_escape(&record.ref_allele),
        csv_escape(&record.alt_allele),
        (max_len - min_len).to_string(),
        indel_type.to_string(),
        csv_escape(&render_filter(&record.filter)),
        csv_escape(&nt),
        usize::from(nt == "ref").to_string(),
        format_info_float(strelka::info_float(info, "SomaticEVS").unwrap_or(-1.0)),
        int_info(info, "QSI_NT").to_string(),
        format_python_float(n_dp),
        format_python_float(t_dp),
        format_python_float(depth.map(|value| n_dp / value).unwrap_or(0.0)),
        format_python_float(depth.map(|value| t_dp / value).unwrap_or(0.0)),
        format_python_float(first(&normal, "BCN50")),
        format_python_float(first(&tumor, "BCN50")),
        format_python_float(first(&normal, "FDP50")),
        format_python_float(first(&tumor, "FDP50")),
        format_python_float(strelka::af_from_tir_tar(&normal, "TIR", "TAR")),
        format_python_float(strelka::af_from_tir_tar(&tumor, "TIR", "TAR")),
        csv_escape(&strelka::info_value(info, "SGT").unwrap_or_default()),
        int_info(info, "RC").to_string(),
        csv_escape(&ru),
        ru.len().to_string(),
        int_info(info, "IC").to_string(),
        int_info(info, "IHP").to_string(),
        format_info_float(strelka::info_float(info, "MQ").unwrap_or(0.0)),
        format_info_float(strelka::info_float(info, "MQ0").unwrap_or(0.0)),
        csv_escape(label),
    ];

    cells.extend(scoring_features.render_evsf(strelka::info_value(info, "EVSF").as_deref()));
    cells.join(",")
}

fn first(sample: &BTreeMap<String, String>, field: &str) -> f64 {
    strelka::parse_first_number(sample.get(field)).unwrap_or(0.0)
}

fn int_info(info: &str, field: &str) -> i64 {
    strelka::info_value(info, field)
        .and_then(|value| {
            value
                .parse::<i64>()
                .or_else(|_| value.parse::<f64>().map(|number| number as i64))
                .ok()
        })
        .unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record() -> RawVcfRecord {
        RawVcfRecord {
            chrom: "chr1".to_string(),
            pos: 10,
            id: ".".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "AT".to_string(),
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            info: "NT=ref;QSI_NT=42;SomaticEVS=17.5;SGT=ref->het;RC=3;RU=T;IC=2;IHP=4;MQ=55;EVSF=1.5,bad".to_string(),
            format: Some("DP:TAR:TIR:BCN50:FDP50".to_string()),
            samples: vec!["10:8,9:2,3:0.5:1".to_string(), "20:5,6:15,16:0.25:2".to_string()],
            mixed_edit_primitive: false,
            primitive_identity: None,
        }
    }

    #[test]
    fn renders_legacy_strelka_indel_shape_and_defaults() {
        let mut depths = BTreeMap::new();
        depths.insert("chr1".to_string(), 40.0);
        let scoring = ScoringFeatures {
            columns: vec!["one".to_string(), "two".to_string()],
            names_by_index: vec!["one".to_string(), "two".to_string()],
        };
        let row = render_row(&record(), &depths, &scoring, 0, "FP");
        let cells: Vec<&str> = row.split(',').collect();
        assert_eq!(
            &cells[1..11],
            &["chr1", "10", "A", "AT", "1", "1", "", "ref", "1", "17.5"]
        );
        assert_eq!(cells[13], "20.0");
        assert_eq!(cells[14], "0.25");
        assert_eq!(cells[15], "0.5");
        assert_eq!(cells[20], "0.2");
        assert_eq!(cells[21], "0.75");
        assert_eq!(cells[31], "1.5");
        assert_eq!(cells[32], "0.0");
    }

    #[test]
    fn later_scoring_headers_remap_values_but_keep_earlier_columns() {
        let headers = vec![
            "##indel_scoring_features=old_first,shared_second".to_string(),
            "##indel_scoring_features=new_first".to_string(),
        ];

        let lines = emit_with_depths(&[record()], &headers, "FP", None);

        assert!(
            lines[0].ends_with(",E.old_first,E.shared_second,E.new_first"),
            "{}",
            lines[0]
        );
        assert!(lines[1].ends_with(",,0.0,1.5"), "{}", lines[1]);
    }
}
