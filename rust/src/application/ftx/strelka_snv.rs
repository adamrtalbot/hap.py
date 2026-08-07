//! Strelka SNV feature-table extractor — `hap ftx --feature-table
//! hcc.strelka.snv`. Port of `src/python/Somatic/Strelka.py::
//! extractStrelkaSNVFeatures` with pandas `to_csv` output semantics.
//!
//! Sample identification follows legacy: column 0 = NORMAL, column 1 =
//! TUMOR (not name-based). Per-chromosome depth normalisation reads
//! Strelka's `##MaxDepth_<chrom>=N` / `##MeanDepth_...` / `##Depth_...`
//! header lines; if absent, `*_DP_RATE` stays at 0.0.

use std::collections::BTreeMap;

use crate::domain::RawVcfRecord;
use crate::engines::strelka;

use super::common::{ScoringFeatures, csv_escape, format_info_float, parse_scoring_features};

/// Fixed column order, excluding the leading pandas index cell and the
/// dynamic `E.<scoring-feature>` tail. The `SomaticEVS` column is kept
/// because pandas emits it, but legacy never populates it (the
/// `extractStrelkaSNVFeatures` qrec dict has no `"SomaticEVS"` key — the
/// column ends up as NaN → empty cell).
const FIXED_COLUMNS: &[&str] = &[
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "NT",
    "NT_REF",
    "QSS_NT",
    "FILTER",
    "SomaticEVS",
    "EVS",
    "VQSR",
    "N_FDP_RATE",
    "T_FDP_RATE",
    "N_SDP_RATE",
    "T_SDP_RATE",
    "N_DP",
    "T_DP",
    "N_DP_RATE",
    "T_DP_RATE",
    "N_AF",
    "T_AF",
    "MQ",
    "MQ0",
    "SNVSB",
    "ReadPosRankSum",
    "tag",
];

#[cfg(test)]
pub(super) fn emit_with_depths(
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depth_override: Option<&BTreeMap<String, f64>>,
) -> Vec<String> {
    emit_with_depths_precision(records, headers, label, depth_override, false)
}

pub(super) fn emit_with_depths_precision(
    records: &[RawVcfRecord],
    headers: &[String],
    label: &str,
    depth_override: Option<&BTreeMap<String, f64>>,
    full_precision: bool,
) -> Vec<String> {
    let scoring_features = parse_scoring_features(headers, "snv_scoring_features");
    let header_depths = strelka::parse_depths(headers);
    let avg_depth = depth_override.unwrap_or(&header_depths);
    let float_columns = FloatColumns::from_records(records, avg_depth);

    let mut header_line = String::from(",");
    header_line.push_str(&FIXED_COLUMNS.join(","));
    for name in &scoring_features.columns {
        header_line.push_str(",E.");
        header_line.push_str(name);
    }
    let mut lines = vec![header_line];

    for (index, record) in records.iter().enumerate() {
        lines.push(render_row(
            record,
            avg_depth,
            &scoring_features,
            &float_columns,
            full_precision,
            index,
            label,
        ));
    }
    lines
}

fn render_row(
    record: &RawVcfRecord,
    avg_depth: &BTreeMap<String, f64>,
    scoring_features: &ScoringFeatures,
    float_columns: &FloatColumns,
    full_precision: bool,
    index: usize,
    label: &str,
) -> String {
    let n = record.sample_map(0);
    let t = record.sample_map(1);
    let info = record.info.as_str();

    let nt = strelka::info_value(info, "NT").unwrap_or_default();
    let nt_ref = if nt == "ref" { 1 } else { 0 };
    let qss_nt = strelka::info_value(info, "QSS_NT")
        .and_then(|value| {
            value
                .parse::<i64>()
                .or_else(|_| value.parse::<f64>().map(|number| number as i64))
                .ok()
        })
        .unwrap_or(0);
    let vqsr = strelka::info_float(info, "VQSR").unwrap_or(-1.0);
    // `vcfExtract` always creates the requested SomaticEVS key, even when
    // absent (with a None value). The legacy key-presence branch therefore
    // returns -1.0 rather than falling back to the older EVS field.
    let evs = strelka::info_float(info, "SomaticEVS").unwrap_or(-1.0);

    let n_dp = strelka::parse_first_number(n.get("DP")).unwrap_or(0.0);
    let t_dp = strelka::parse_first_number(t.get("DP")).unwrap_or(0.0);
    let n_fdp = strelka::parse_first_number(n.get("FDP")).unwrap_or(0.0);
    let t_fdp = strelka::parse_first_number(t.get("FDP")).unwrap_or(0.0);
    let n_sdp = strelka::parse_first_number(n.get("SDP")).unwrap_or(0.0);
    let t_sdp = strelka::parse_first_number(t.get("SDP")).unwrap_or(0.0);

    let n_fdp_rate = if n_dp != 0.0 { n_fdp / n_dp } else { 0.0 };
    let t_fdp_rate = if t_dp != 0.0 { t_fdp / t_dp } else { 0.0 };
    let n_sdp_rate = if n_dp + n_sdp != 0.0 {
        n_sdp / (n_dp + n_sdp)
    } else {
        0.0
    };
    let t_sdp_rate = if t_dp + t_sdp != 0.0 {
        t_sdp / (t_dp + t_sdp)
    } else {
        0.0
    };

    let norm = avg_depth
        .get(&record.chrom)
        .copied()
        .filter(|depth| *depth != 0.0);
    let n_dp_rate = norm.map(|d| n_dp / d).unwrap_or(0.0);
    let t_dp_rate = norm.map(|d| t_dp / d).unwrap_or(0.0);

    let (n_ref_t1, n_alt_t1) = tier1_ref_alt(&n, &record.ref_allele, &record.alt_allele);
    let (t_ref_t1, t_alt_t1) = tier1_ref_alt(&t, &record.ref_allele, &record.alt_allele);
    let n_af = if n_ref_t1 + n_alt_t1 != 0.0 {
        n_alt_t1 / (n_ref_t1 + n_alt_t1)
    } else {
        0.0
    };
    let t_af = if t_ref_t1 + t_alt_t1 != 0.0 {
        t_alt_t1 / (t_ref_t1 + t_alt_t1)
    } else {
        0.0
    };

    let mq_cell = info_float_or_zero(info, "MQ");
    let mq0_cell = info_float_or_zero(info, "MQ0");
    // These two fields are copied directly into the legacy row after its
    // missing-feature repair. A present VCF float therefore renders as a
    // pandas float, while an absent value is the literal integer `0` that
    // the repair inserted (unlike MQ/MQ0, which are always float-cast).
    let snvsb_cell = info_float_or_integer_zero(info, "SNVSB");
    let rprs_cell = info_float_or_integer_zero(info, "ReadPosRankSum");
    let filter_cell = csv_escape(&render_strelka_filter(&record.filter));

    let mut cells: Vec<String> =
        Vec::with_capacity(FIXED_COLUMNS.len() + scoring_features.columns.len() + 1);
    cells.push(index.to_string());
    cells.push(record.chrom.clone());
    cells.push(record.pos.to_string());
    cells.push(csv_escape(&record.ref_allele));
    cells.push(csv_escape(&record.alt_allele));
    cells.push(csv_escape(&nt));
    cells.push(nt_ref.to_string());
    cells.push(qss_nt.to_string());
    cells.push(filter_cell);
    cells.push(String::new()); // SomaticEVS: qrec never writes it
    cells.push(format_staged_float(evs, full_precision));
    cells.push(format_staged_float(vqsr, full_precision));
    cells.push(format_inferred_rate(
        n_fdp_rate,
        float_columns.n_fdp,
        full_precision,
    ));
    cells.push(format_inferred_rate(
        t_fdp_rate,
        float_columns.t_fdp,
        full_precision,
    ));
    cells.push(format_inferred_rate(
        n_sdp_rate,
        float_columns.n_sdp,
        full_precision,
    ));
    cells.push(format_inferred_rate(
        t_sdp_rate,
        float_columns.t_sdp,
        full_precision,
    ));
    cells.push(format_staged_float(n_dp, full_precision));
    cells.push(format_staged_float(t_dp, full_precision));
    cells.push(format_inferred_rate(
        n_dp_rate,
        float_columns.dp_rate,
        full_precision,
    ));
    cells.push(format_inferred_rate(
        t_dp_rate,
        float_columns.dp_rate,
        full_precision,
    ));
    cells.push(format_inferred_rate(
        n_af,
        float_columns.n_af,
        full_precision,
    ));
    cells.push(format_inferred_rate(
        t_af,
        float_columns.t_af,
        full_precision,
    ));
    cells.push(if full_precision {
        info_float_or_zero_precision(info, "MQ", true)
    } else {
        mq_cell
    });
    cells.push(if full_precision {
        info_float_or_zero_precision(info, "MQ0", true)
    } else {
        mq0_cell
    });
    cells.push(snvsb_cell);
    cells.push(rprs_cell);
    cells.push(csv_escape(label));

    // E.<scoring-feature> columns from INFO/EVSF, one per name in the
    // `##snv_scoring_features` header. Missing / unparseable → 0.0;
    // pandas float-coerces the column so even all-missing renders as
    // "0.0" when any other row had a real float.
    cells.extend(scoring_features.render_evsf(strelka::info_value(info, "EVSF").as_deref()));

    cells.join(",")
}

#[derive(Default)]
struct FloatColumns {
    n_fdp: bool,
    t_fdp: bool,
    n_sdp: bool,
    t_sdp: bool,
    dp_rate: bool,
    n_af: bool,
    t_af: bool,
}

impl FloatColumns {
    fn from_records(records: &[RawVcfRecord], avg_depth: &BTreeMap<String, f64>) -> Self {
        let mut columns = Self::default();
        for record in records {
            let normal = record.sample_map(0);
            let tumor = record.sample_map(1);
            let n_dp = strelka::parse_first_number(normal.get("DP")).unwrap_or(0.0);
            let t_dp = strelka::parse_first_number(tumor.get("DP")).unwrap_or(0.0);
            let n_sdp = strelka::parse_first_number(normal.get("SDP")).unwrap_or(0.0);
            let t_sdp = strelka::parse_first_number(tumor.get("SDP")).unwrap_or(0.0);
            columns.n_fdp |= n_dp != 0.0;
            columns.t_fdp |= t_dp != 0.0;
            columns.n_sdp |= n_dp + n_sdp != 0.0;
            columns.t_sdp |= t_dp + t_sdp != 0.0;
            columns.dp_rate |= avg_depth
                .get(&record.chrom)
                .is_some_and(|depth| *depth != 0.0);

            let (n_ref, n_alt) = tier1_ref_alt(&normal, &record.ref_allele, &record.alt_allele);
            let (t_ref, t_alt) = tier1_ref_alt(&tumor, &record.ref_allele, &record.alt_allele);
            columns.n_af |= n_ref + n_alt != 0.0;
            columns.t_af |= t_ref + t_alt != 0.0;
        }
        columns
    }
}

fn format_inferred_rate(value: f64, float_column: bool, full_precision: bool) -> String {
    if float_column {
        format_staged_float(value, full_precision)
    } else {
        "0".to_string()
    }
}

fn format_staged_float(value: f64, full_precision: bool) -> String {
    if full_precision {
        (if value == 0.0 { 0.0 } else { value }).to_string()
    } else {
        format_info_float(value)
    }
}

fn info_float_or_zero(info: &str, key: &str) -> String {
    info_float_or_zero_precision(info, key, false)
}

fn info_float_or_zero_precision(info: &str, key: &str, full_precision: bool) -> String {
    format_staged_float(
        strelka::info_float(info, key).unwrap_or(0.0),
        full_precision,
    )
}

fn info_float_or_integer_zero(info: &str, key: &str) -> String {
    strelka::info_float(info, key)
        .map(format_info_float)
        .unwrap_or_else(|| "0".to_string())
}

/// Legacy joins `rec["FILTER"]` (a list split on commas by vcfExtract)
/// with commas. Standard semicolon-separated FILTER tokens are therefore
/// preserved verbatim; PASS and missing collapse to the empty list.
fn render_strelka_filter(filter: &str) -> String {
    if filter == "PASS" || filter == "." || filter.is_empty() {
        String::new()
    } else {
        filter.to_string()
    }
}

fn tier1_ref_alt(
    sample: &BTreeMap<String, String>,
    ref_allele: &str,
    alt_allele: &str,
) -> (f64, f64) {
    (
        tier1_for_base(sample, ref_allele),
        tier1_for_alt_list(sample, alt_allele),
    )
}

fn tier1_for_base(sample: &BTreeMap<String, String>, base: &str) -> f64 {
    if base.len() != 1 {
        return 0.0;
    }
    let key = format!("{}U", base.to_ascii_uppercase());
    strelka::parse_first_number(sample.get(&key)).unwrap_or(0.0)
}

/// Mirrors the legacy `try: sum += rec['S.x.' + a + 'U'][0] except: [0,0]`
/// short-circuit — any non-SNV alt or missing FORMAT field aborts the
/// sum, matching Python's broad `except` catch.
fn tier1_for_alt_list(sample: &BTreeMap<String, String>, alts: &str) -> f64 {
    let mut sum = 0.0;
    for a in alts.split(',') {
        if a.len() != 1 {
            return 0.0;
        }
        let key = format!("{}U", a.to_ascii_uppercase());
        match strelka::parse_first_number(sample.get(&key)) {
            Some(v) => sum += v,
            None => return 0.0,
        }
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_record(info: &str, fmt: &str, n: &str, t: &str) -> RawVcfRecord {
        RawVcfRecord {
            chrom: "chr21".to_string(),
            pos: 9412105,
            id: ".".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            info: info.to_string(),
            format: Some(fmt.to_string()),
            samples: vec![n.to_string(), t.to_string()],
        }
    }

    fn render_test_row(
        record: &RawVcfRecord,
        avg_depth: &BTreeMap<String, f64>,
        scoring_features: &ScoringFeatures,
        index: usize,
        label: &str,
    ) -> String {
        render_row(
            record,
            avg_depth,
            scoring_features,
            &FloatColumns::from_records(std::slice::from_ref(record), avg_depth),
            false,
            index,
            label,
        )
    }

    #[test]
    fn strelka_first_admix_record_matches_legacy_byte_exact() {
        // Derived from example/sompy/strelka_admix_snvs.vcf.gz row 0:
        //   chr21:9412105 A>T PASS NT=ref QSS_NT=61 VQSR=15.09 MQ=57.82 MQ0=8 SNVSB=0 ReadPosRankSum=-0.12
        //   Normal DP=26 FDP=0 SDP=0 AU=26,29 CU=0,0 GU=0,0 TU=0,0
        //   Tumor  DP=77 FDP=0 SDP=0 AU=55,56 CU=0,0 GU=0,0 TU=22,27
        //   With ##MaxDepth_chr21=135.50 → N_DP_RATE=26/135.5, T_DP_RATE=77/135.5
        let rec = mk_record(
            "SOMATIC;QSS=61;NT=ref;QSS_NT=61;SGT=AA->AT;DP=112;MQ=57.82;MQ0=8;\
             ReadPosRankSum=-0.12;SNVSB=0;VQSR=15.09",
            "DP:FDP:SDP:SUBDP:AU:CU:GU:TU",
            "26:0:0:0:26,29:0,0:0,0:0,0",
            "77:0:0:0:55,56:0,0:0,0:22,27",
        );
        let mut avg = BTreeMap::new();
        avg.insert("chr21".to_string(), 135.5);
        let line = render_test_row(
            &rec,
            &avg,
            &ScoringFeatures::default(),
            0,
            "strelka_admix_snvs.vcf.gz",
        );
        let expected = "0,chr21,9412105,A,T,ref,1,61,,,-1.0,15.09,0.0,0.0,0.0,0.0,\
                        26.0,77.0,0.191881918819,0.568265682657,0.0,\
                        0.285714285714,57.82,8.0,0.0,-0.12,strelka_admix_snvs.vcf.gz";
        assert_eq!(line, expected);
    }

    #[test]
    fn strelka_tier1_ref_alt_sums_alt_bases() {
        let sample: BTreeMap<String, String> = [
            ("AU", "26,29"),
            ("CU", "0,0"),
            ("GU", "0,0"),
            ("TU", "22,27"),
        ]
        .iter()
        .map(|(k, v)| (k.to_string(), v.to_string()))
        .collect();
        assert_eq!(tier1_for_base(&sample, "A"), 26.0);
        assert_eq!(tier1_for_alt_list(&sample, "T"), 22.0);
        assert_eq!(tier1_for_alt_list(&sample, "T,C"), 22.0); // 22+0
    }

    #[test]
    fn strelka_tier1_alt_non_single_base_returns_zero() {
        let sample: BTreeMap<String, String> = [("AU", "26,29"), ("TU", "22,27")]
            .iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        assert_eq!(tier1_for_alt_list(&sample, "CG"), 0.0);
    }

    #[test]
    fn strelka_filter_pass_and_missing_render_empty() {
        assert_eq!(render_strelka_filter("PASS"), "");
        assert_eq!(render_strelka_filter("."), "");
        assert_eq!(render_strelka_filter(""), "");
        assert_eq!(
            render_strelka_filter("LowQscore;HighDepth"),
            "LowQscore;HighDepth"
        );
    }

    #[test]
    fn strelka_missing_numeric_features_render_zero() {
        let rec = mk_record(
            "NT=ref;QSS_NT=10;SNVSB=0;VQSR=1.0", // no MQ
            "DP:FDP:SDP:AU:CU:GU:TU",
            "10:0:0:10,10:0,0:0,0:0,0",
            "10:0:0:0,0:0,0:0,0:10,10",
        );
        let avg: BTreeMap<String, f64> = BTreeMap::new();
        let line = render_test_row(&rec, &avg, &ScoringFeatures::default(), 0, "t");
        // MQ column is the 23rd field (0-based index 22); find it.
        let fields: Vec<&str> = line.split(',').collect();
        assert_eq!(fields[22], "0.0");
        assert_eq!(fields[23], "0.0");
        assert_eq!(fields[24], "0.0");
        assert_eq!(fields[25], "0");
    }

    #[test]
    fn later_scoring_headers_remap_values_but_keep_earlier_columns() {
        let rec = mk_record(
            "NT=ref;QSS_NT=10;EVSF=1.5,2.5",
            "DP:FDP:SDP:AU:CU:GU:TU",
            "10:0:0:10,10:0,0:0,0:0,0",
            "10:0:0:0,0:0,0:0,0:10,10",
        );
        let headers = vec![
            "##snv_scoring_features=old_first,shared_second".to_string(),
            "##snv_scoring_features=new_first".to_string(),
        ];

        let lines = emit_with_depths(&[rec], &headers, "tag", None);

        assert!(
            lines[0].ends_with(",E.old_first,E.shared_second,E.new_first"),
            "{}",
            lines[0]
        );
        assert!(lines[1].ends_with(",,2.5,1.5"), "{}", lines[1]);
    }

    #[test]
    fn scalar_evsf_is_not_treated_as_a_feature_list() {
        let rec = mk_record(
            "NT=ref;QSS_NT=10;EVSF=1.5",
            "DP:FDP:SDP:AU:CU:GU:TU",
            "10:0:0:10,10:0,0:0,0:0,0",
            "10:0:0:0,0:0,0:0,0:10,10",
        );

        let lines = emit_with_depths(
            &[rec],
            &["##snv_scoring_features=only".to_string()],
            "tag",
            None,
        );

        assert!(lines[1].ends_with(",0.0"), "{}", lines[1]);
    }

    #[test]
    fn missing_bias_features_keep_legacy_integer_zero_shape() {
        let rec = mk_record(
            "NT=ref;QSS_NT=10",
            "DP:FDP:SDP:AU:CU:GU:TU",
            "10:0:0:10,10:0,0:0,0:0,0",
            "10:0:0:0,0:0,0:0,0:10,10",
        );
        let line = render_test_row(&rec, &BTreeMap::new(), &ScoringFeatures::default(), 0, "t");
        let fields: Vec<&str> = line.split(',').collect();
        assert_eq!(fields[22], "0.0");
        assert_eq!(fields[23], "0.0");
        assert_eq!(fields[24], "0");
        assert_eq!(fields[25], "0");
    }

    #[test]
    fn older_evs_does_not_bypass_missing_somatic_evs_branch() {
        let rec = mk_record(
            "NT=ref;QSS_NT=10;EVS=17.5",
            "DP:FDP:SDP:AU:CU:GU:TU",
            "10:0:0:10,10:0,0:0,0:0,0",
            "10:0:0:0,0:0,0:0,0:10,10",
        );
        let line = render_test_row(&rec, &BTreeMap::new(), &ScoringFeatures::default(), 0, "t");
        let fields: Vec<&str> = line.split(',').collect();
        assert_eq!(fields[10], "-1.0");
    }

    #[test]
    fn comma_filter_is_csv_quoted_without_rewriting_semicolons() {
        let mut rec = mk_record(
            "NT=ref;QSS_NT=10",
            "DP:FDP:SDP:AU:CU:GU:TU",
            "10:0:0:10,10:0,0:0,0:0,0",
            "10:0:0:0,0:0,0:0,0:10,10",
        );
        rec.filter = "First,Second".to_string();
        let line = render_test_row(&rec, &BTreeMap::new(), &ScoringFeatures::default(), 0, "t");
        assert!(line.contains(",\"First,Second\","));
    }
}
