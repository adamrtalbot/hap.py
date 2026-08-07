//! Cohesive responsibility extracted from the command façade.

use crate::application::ftx;
use crate::domain::RawVcfRecord;
use crate::engines::strelka;
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

pub(super) fn renumber_feature_rows(groups: &[Vec<String>]) -> Vec<String> {
    let mut out = Vec::new();
    for group in groups {
        for (idx, row) in group.iter().enumerate() {
            let tail = row.split_once(',').map(|(_, tail)| tail).unwrap_or(row);
            out.push(format!("{idx},{tail}"));
        }
    }
    out
}

pub(super) struct CallerFeatureTable {
    pub(super) header: String,
    pub(super) tp: Vec<String>,
    pub(super) fp: Vec<String>,
    pub(super) fn_rows: Vec<String>,
    pub(super) ambi: Vec<String>,
    pub(super) unk: Vec<String>,
}

pub(super) struct CallerRecordGroups<'a> {
    pub(super) tp_truth: &'a [RawVcfRecord],
    pub(super) tp_query: &'a [RawVcfRecord],
    pub(super) fn_truth: &'a [RawVcfRecord],
    pub(super) fp_query: &'a [RawVcfRecord],
    pub(super) ambi_query: &'a [RawVcfRecord],
    pub(super) unk_query: &'a [RawVcfRecord],
}

pub(super) struct ParsedFeatureTable {
    columns: BTreeSet<String>,
    rows: Vec<BTreeMap<String, String>>,
}

impl ParsedFeatureTable {
    fn from_lines(lines: Vec<String>) -> Result<Self> {
        let mut lines = lines.into_iter();
        let header = parse_csv_line(&lines.next().context("feature table has no header")?);
        let named_columns = header.iter().skip(1).cloned().collect::<Vec<_>>();
        let columns = named_columns.iter().cloned().collect();
        let rows = lines
            .map(|line| {
                let cells = parse_csv_line(&line);
                named_columns
                    .iter()
                    .cloned()
                    .zip(cells.into_iter().skip(1))
                    .collect()
            })
            .collect();
        Ok(Self { columns, rows })
    }
}

pub(super) fn build_caller_feature_table(
    feature: &str,
    truth_headers: &[String],
    query_headers: &[String],
    depths: Option<&BTreeMap<String, f64>>,
    check_order: bool,
    groups: &CallerRecordGroups<'_>,
) -> Result<CallerFeatureTable> {
    let truth_tp = ParsedFeatureTable::from_lines(ftx::emit_feature_table_for_somatic(
        feature,
        groups.tp_truth,
        truth_headers,
        "TP",
        depths,
    )?)?;
    let query_tp = ParsedFeatureTable::from_lines(ftx::emit_feature_table_for_somatic(
        feature,
        groups.tp_query,
        query_headers,
        "TP_r",
        depths,
    )?)?;
    if truth_tp.rows.len() != query_tp.rows.len() {
        bail!(
            "cannot merge TP features: truth and query lengths differ ({} != {})",
            truth_tp.rows.len(),
            query_tp.rows.len()
        );
    }
    if check_order {
        for (index, (truth, query)) in truth_tp.rows.iter().zip(&query_tp.rows).enumerate() {
            for column in ["CHROM", "POS"] {
                if truth.get(column) != query.get(column) {
                    bail!(
                        "cannot merge TP features: inputs are out of order at row {index} ({column})"
                    );
                }
            }
        }
    }

    let (tp_columns, tp) = merge_caller_tp_tables(&truth_tp, &query_tp);
    let truth_fn = ParsedFeatureTable::from_lines(ftx::emit_feature_table_for_somatic(
        feature,
        groups.fn_truth,
        truth_headers,
        "FN",
        depths,
    )?)?;
    let query_fp = ParsedFeatureTable::from_lines(ftx::emit_feature_table_for_somatic(
        feature,
        groups.fp_query,
        query_headers,
        "FP",
        depths,
    )?)?;
    let query_ambi = ParsedFeatureTable::from_lines(ftx::emit_feature_table_for_somatic(
        feature,
        groups.ambi_query,
        query_headers,
        "AMBI",
        depths,
    )?)?;
    let query_unk = ParsedFeatureTable::from_lines(ftx::emit_feature_table_for_somatic(
        feature,
        groups.unk_query,
        query_headers,
        "UNK",
        depths,
    )?)?;
    let (fn_columns, fn_rows) = suffix_truth_columns(truth_fn, &tp_columns);

    let mut all_columns = tp_columns;
    all_columns.extend(fn_columns);
    all_columns.extend(query_fp.columns.iter().cloned());
    all_columns.extend(query_ambi.columns.iter().cloned());
    all_columns.extend(query_unk.columns.iter().cloned());
    let columns = ordered_somatic_feature_columns(&all_columns);
    let header = format!(",{}", columns.join(","));

    Ok(CallerFeatureTable {
        header,
        tp: render_feature_group(&tp, &columns),
        fp: render_feature_group(&query_fp.rows, &columns),
        fn_rows: render_feature_group(&fn_rows, &columns),
        ambi: render_feature_group(&query_ambi.rows, &columns),
        unk: render_feature_group(&query_unk.rows, &columns),
    })
}

pub(super) fn merge_caller_tp_tables(
    truth: &ParsedFeatureTable,
    query: &ParsedFeatureTable,
) -> (BTreeSet<String>, Vec<BTreeMap<String, String>>) {
    let shared_keys = ["CHROM", "POS", "tag"];
    let union = truth
        .columns
        .union(&query.columns)
        .cloned()
        .collect::<BTreeSet<_>>();
    let mut columns = BTreeSet::new();
    for column in &union {
        if shared_keys.contains(&column.as_str()) {
            columns.insert(column.clone());
        } else if truth.columns.contains(column) && query.columns.contains(column) {
            columns.insert(column.clone());
            columns.insert(format!("{column}.truth"));
        } else {
            columns.insert(column.clone());
        }
    }

    let rows = truth
        .rows
        .iter()
        .zip(&query.rows)
        .map(|(truth_row, query_row)| {
            let mut row = BTreeMap::new();
            for column in &union {
                if shared_keys.contains(&column.as_str()) {
                    if let Some(value) = truth_row.get(column).or_else(|| query_row.get(column)) {
                        row.insert(column.clone(), value.clone());
                    }
                } else if let (Some(truth_value), Some(query_value)) =
                    (truth_row.get(column), query_row.get(column))
                {
                    row.insert(column.clone(), query_value.clone());
                    row.insert(format!("{column}.truth"), truth_value.clone());
                } else if let Some(value) = query_row.get(column).or_else(|| truth_row.get(column))
                {
                    row.insert(column.clone(), value.clone());
                }
            }
            row
        })
        .collect();
    (columns, rows)
}

pub(super) fn suffix_truth_columns(
    table: ParsedFeatureTable,
    tp_columns: &BTreeSet<String>,
) -> (BTreeSet<String>, Vec<BTreeMap<String, String>>) {
    let rename = |column: &str| {
        let truth_name = format!("{column}.truth");
        if tp_columns.contains(&truth_name) {
            truth_name
        } else {
            column.to_string()
        }
    };
    let columns = table.columns.iter().map(|column| rename(column)).collect();
    let rows = table
        .rows
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|(column, value)| (rename(&column), value))
                .collect()
        })
        .collect();
    (columns, rows)
}

pub(super) fn ordered_somatic_feature_columns(columns: &BTreeSet<String>) -> Vec<String> {
    let first = [
        "CHROM",
        "POS",
        "tag",
        "REF",
        "REF.truth",
        "ALT",
        "ALT.truth",
    ];
    first
        .iter()
        .filter(|column| columns.contains(**column))
        .map(|column| (*column).to_string())
        .chain(
            columns
                .iter()
                .filter(|column| !first.contains(&column.as_str()))
                .cloned(),
        )
        .collect()
}

pub(super) fn render_feature_group(
    rows: &[BTreeMap<String, String>],
    columns: &[String],
) -> Vec<String> {
    rows.iter()
        .enumerate()
        .map(|(index, row)| {
            csv_join(
                std::iter::once(index.to_string()).chain(columns.iter().map(|column| {
                    format_somatic_feature_cell(column, row.get(column).map_or("", String::as_str))
                })),
            )
        })
        .collect()
}

pub(super) fn format_somatic_feature_cell(column: &str, value: &str) -> String {
    const TEXT_COLUMNS: &[&str] = &[
        "CHROM",
        "tag",
        "REF",
        "REF.truth",
        "ALT",
        "ALT.truth",
        "FILTER",
        "NT",
        "SGT",
        "S.2.GT",
    ];
    if value.is_empty() || column == "POS" || TEXT_COLUMNS.contains(&column) {
        value.to_string()
    } else if let Ok(number) = value.parse::<f64>() {
        format!("{number:.8}")
    } else {
        value.to_string()
    }
}

pub(super) fn render_generic_tp_row(
    index: usize,
    truth: &RawVcfRecord,
    query: &RawVcfRecord,
) -> String {
    csv_join([
        index.to_string(),
        query.chrom.clone(),
        query.pos.to_string(),
        "TP".to_string(),
        query.ref_allele.clone(),
        truth.ref_allele.clone(),
        query.alt_allele.clone(),
        truth.alt_allele.clone(),
        normalized_filter(&query.filter),
        normalized_filter(&truth.filter),
        format_feature_number_or_text(&normalized_qual(&query.qual)),
        format_feature_number_or_text(&normalized_qual(&truth.qual)),
    ])
}

pub(super) fn render_generic_fn_row(index: usize, truth: &RawVcfRecord) -> String {
    csv_join([
        index.to_string(),
        truth.chrom.clone(),
        truth.pos.to_string(),
        "FN".to_string(),
        String::new(),
        truth.ref_allele.clone(),
        String::new(),
        truth.alt_allele.clone(),
        String::new(),
        normalized_filter(&truth.filter),
        String::new(),
        format_feature_number_or_text(&normalized_qual(&truth.qual)),
    ])
}

pub(super) fn render_generic_query_row(index: usize, query: &RawVcfRecord, tag: &str) -> String {
    csv_join([
        index.to_string(),
        query.chrom.clone(),
        query.pos.to_string(),
        tag.to_string(),
        query.ref_allele.clone(),
        String::new(),
        query.alt_allele.clone(),
        String::new(),
        normalized_filter(&query.filter),
        String::new(),
        format_feature_number_or_text(&normalized_qual(&query.qual)),
        String::new(),
    ])
}

pub(super) const STRELKA_HCC_INDEL_HEADER: &str = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,EVS,FILTER,I.DP_normal,I.DP_tumor,I.T_ALT_RATE,I.count,I.tag,IC,IHP,INDELTYPE,LENGTH,MQ,MQ0,NT,NT_REF,N_AF,N_BCN,N_DP,N_DP_RATE,N_FDP,QSI_NT,QUAL,RC,RU,RU_LEN,S.1.VT,SGT,T_AF,T_BCN,T_DP,T_DP_RATE,T_FDP";

pub(super) fn render_strelka_hcc_indel_tp_row(
    index: usize,
    truth: &RawVcfRecord,
    query: &RawVcfRecord,
    avg_depth: &BTreeMap<String, f64>,
) -> String {
    let query_row = strelka_indel_row(query, avg_depth);
    let truth_row = generic_truth_row(truth);
    merge_feature_rows(index, "TP", &query_row, &truth_row)
}

pub(super) fn render_strelka_hcc_indel_fn_row(index: usize, truth: &RawVcfRecord) -> String {
    let query_row = blank_strelka_query_row(truth);
    let truth_row = generic_truth_row(truth);
    merge_feature_rows(index, "FN", &query_row, &truth_row)
}

pub(super) fn render_strelka_hcc_indel_query_row(
    index: usize,
    query: &RawVcfRecord,
    tag: &str,
    avg_depth: &BTreeMap<String, f64>,
) -> String {
    let query_row = strelka_indel_row(query, avg_depth);
    merge_feature_rows(index, tag, &query_row, &BTreeMap::new())
}

pub(super) fn merge_feature_rows(
    index: usize,
    tag: &str,
    query_row: &BTreeMap<String, String>,
    truth_row: &BTreeMap<String, String>,
) -> String {
    let cols = vec![
        index.to_string(),
        query_row.get("CHROM").cloned().unwrap_or_default(),
        query_row.get("POS").cloned().unwrap_or_default(),
        tag.to_string(),
        query_row.get("REF").cloned().unwrap_or_default(),
        truth_row.get("REF").cloned().unwrap_or_default(),
        query_row.get("ALT").cloned().unwrap_or_default(),
        truth_row.get("ALT").cloned().unwrap_or_default(),
        query_row.get("EVS").cloned().unwrap_or_default(),
        query_row.get("FILTER").cloned().unwrap_or_default(),
        truth_row.get("I.DP_normal").cloned().unwrap_or_default(),
        truth_row.get("I.DP_tumor").cloned().unwrap_or_default(),
        truth_row.get("I.T_ALT_RATE").cloned().unwrap_or_default(),
        truth_row.get("I.count").cloned().unwrap_or_default(),
        truth_row.get("I.tag").cloned().unwrap_or_default(),
        query_row.get("IC").cloned().unwrap_or_default(),
        query_row.get("IHP").cloned().unwrap_or_default(),
        query_row.get("INDELTYPE").cloned().unwrap_or_default(),
        query_row.get("LENGTH").cloned().unwrap_or_default(),
        query_row.get("MQ").cloned().unwrap_or_default(),
        query_row.get("MQ0").cloned().unwrap_or_default(),
        query_row.get("NT").cloned().unwrap_or_default(),
        query_row.get("NT_REF").cloned().unwrap_or_default(),
        query_row.get("N_AF").cloned().unwrap_or_default(),
        query_row.get("N_BCN").cloned().unwrap_or_default(),
        query_row.get("N_DP").cloned().unwrap_or_default(),
        query_row.get("N_DP_RATE").cloned().unwrap_or_default(),
        query_row.get("N_FDP").cloned().unwrap_or_default(),
        query_row.get("QSI_NT").cloned().unwrap_or_default(),
        truth_row.get("QUAL").cloned().unwrap_or_default(),
        query_row.get("RC").cloned().unwrap_or_default(),
        query_row.get("RU").cloned().unwrap_or_default(),
        query_row.get("RU_LEN").cloned().unwrap_or_default(),
        truth_row.get("S.1.VT").cloned().unwrap_or_default(),
        query_row.get("SGT").cloned().unwrap_or_default(),
        query_row.get("T_AF").cloned().unwrap_or_default(),
        query_row.get("T_BCN").cloned().unwrap_or_default(),
        query_row.get("T_DP").cloned().unwrap_or_default(),
        query_row.get("T_DP_RATE").cloned().unwrap_or_default(),
        query_row.get("T_FDP").cloned().unwrap_or_default(),
    ];
    csv_join(cols)
}

pub(super) fn csv_join<I>(cols: I) -> String
where
    I: IntoIterator<Item = String>,
{
    cols.into_iter()
        .map(|value| csv_escape(&value))
        .collect::<Vec<_>>()
        .join(",")
}

pub(super) fn csv_escape(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

pub(super) fn blank_strelka_query_row(truth: &RawVcfRecord) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("CHROM".to_string(), truth.chrom.clone());
    row.insert("POS".to_string(), truth.pos.to_string());
    row.insert("REF".to_string(), String::new());
    row.insert("ALT".to_string(), String::new());
    for key in [
        "EVS",
        "FILTER",
        "IC",
        "IHP",
        "INDELTYPE",
        "LENGTH",
        "MQ",
        "MQ0",
        "NT",
        "NT_REF",
        "N_AF",
        "N_BCN",
        "N_DP",
        "N_DP_RATE",
        "N_FDP",
        "QSI_NT",
        "RC",
        "RU",
        "RU_LEN",
        "SGT",
        "T_AF",
        "T_BCN",
        "T_DP",
        "T_DP_RATE",
        "T_FDP",
    ] {
        row.insert(key.to_string(), String::new());
    }
    row
}

pub(super) fn generic_truth_row(record: &RawVcfRecord) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    row.insert("REF".to_string(), record.ref_allele.clone());
    row.insert("ALT".to_string(), record.alt_allele.clone());
    row.insert(
        "QUAL".to_string(),
        format_feature_number_or_text(&normalized_qual(&record.qual)),
    );
    let sample0 = record.sample_map(0);
    for key in [
        "I.T_ALT_RATE",
        "I.DP_normal",
        "I.DP_tumor",
        "I.tag",
        "I.count",
    ] {
        row.insert(
            key.to_string(),
            format_feature_number_or_text(
                &strelka::info_value(&record.info, key.trim_start_matches("I."))
                    .unwrap_or_default(),
            ),
        );
    }
    row.insert(
        "S.1.VT".to_string(),
        sample0.get("VT").cloned().unwrap_or_default(),
    );
    row
}

pub(super) fn format_feature_number_or_text(value: &str) -> String {
    value
        .parse::<f64>()
        .map(|number| format!("{number:.8}"))
        .unwrap_or_else(|_| value.to_string())
}

pub(super) fn strelka_indel_row(
    record: &RawVcfRecord,
    avg_depth: &BTreeMap<String, f64>,
) -> BTreeMap<String, String> {
    let mut row = BTreeMap::new();
    let n = record.sample_map(0);
    let t = record.sample_map(1);
    let ref_len = record.ref_allele.len();
    let alt_lens: Vec<usize> = record.alt_allele.split(',').map(str::len).collect();
    let max_len = alt_lens
        .iter()
        .copied()
        .max()
        .unwrap_or(ref_len)
        .max(ref_len);
    let min_len = alt_lens
        .iter()
        .copied()
        .min()
        .unwrap_or(ref_len)
        .min(ref_len);
    let mut indel_type = 0;
    for alt in record.alt_allele.split(',') {
        if alt.len() > ref_len {
            indel_type |= 1;
        } else {
            indel_type |= 2;
        }
    }
    let nt = strelka::info_value(&record.info, "NT").unwrap_or_default();
    let n_dp = strelka::parse_first_number(n.get("DP")).unwrap_or(0.0);
    let t_dp = strelka::parse_first_number(t.get("DP")).unwrap_or(0.0);
    let norm = avg_depth.get(&record.chrom).copied().unwrap_or(0.0);

    row.insert("CHROM".to_string(), record.chrom.clone());
    row.insert("POS".to_string(), record.pos.to_string());
    row.insert("REF".to_string(), record.ref_allele.clone());
    row.insert("ALT".to_string(), record.alt_allele.clone());
    row.insert(
        "EVS".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "SomaticEVS")
                .or_else(|| strelka::info_float(&record.info, "EVS"))
                .unwrap_or(-1.0)
        ),
    );
    row.insert("FILTER".to_string(), normalized_filter(&record.filter));
    row.insert(
        "IC".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "IC").unwrap_or(0.0)
        ),
    );
    row.insert(
        "IHP".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "IHP").unwrap_or(0.0)
        ),
    );
    row.insert("INDELTYPE".to_string(), format!("{:.8}", indel_type as f64));
    row.insert(
        "LENGTH".to_string(),
        format!("{:.8}", (max_len - min_len) as f64),
    );
    row.insert(
        "MQ".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "MQ").unwrap_or(0.0)
        ),
    );
    row.insert(
        "MQ0".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "MQ0").unwrap_or(0.0)
        ),
    );
    row.insert("NT".to_string(), nt.clone());
    row.insert(
        "NT_REF".to_string(),
        format!("{:.8}", if nt == "ref" { 1.0 } else { 0.0 }),
    );
    row.insert(
        "N_AF".to_string(),
        format!("{:.8}", strelka::af_from_tir_tar(&n, "TIR", "TAR")),
    );
    row.insert(
        "N_BCN".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(n.get("BCN50")).unwrap_or(0.0)
        ),
    );
    row.insert("N_DP".to_string(), format!("{:.8}", n_dp));
    row.insert(
        "N_DP_RATE".to_string(),
        format!("{:.8}", if norm > 0.0 { n_dp / norm } else { 0.0 }),
    );
    row.insert(
        "N_FDP".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(n.get("FDP50")).unwrap_or(0.0)
        ),
    );
    row.insert(
        "QSI_NT".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "QSI_NT").unwrap_or(0.0)
        ),
    );
    row.insert(
        "RC".to_string(),
        format!(
            "{:.8}",
            strelka::info_float(&record.info, "RC").unwrap_or(0.0)
        ),
    );
    row.insert(
        "RU".to_string(),
        strelka::info_value(&record.info, "RU").unwrap_or_default(),
    );
    row.insert(
        "RU_LEN".to_string(),
        format!(
            "{:.8}",
            strelka::info_value(&record.info, "RU")
                .map(|v| v.len() as f64)
                .unwrap_or(0.0)
        ),
    );
    row.insert(
        "SGT".to_string(),
        strelka::info_value(&record.info, "SGT").unwrap_or_default(),
    );
    row.insert(
        "T_AF".to_string(),
        format!("{:.8}", strelka::af_from_tir_tar(&t, "TIR", "TAR")),
    );
    row.insert(
        "T_BCN".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(t.get("BCN50")).unwrap_or(0.0)
        ),
    );
    row.insert("T_DP".to_string(), format!("{:.8}", t_dp));
    row.insert(
        "T_DP_RATE".to_string(),
        format!("{:.8}", if norm > 0.0 { t_dp / norm } else { 0.0 }),
    );
    row.insert(
        "T_FDP".to_string(),
        format!(
            "{:.8}",
            strelka::parse_first_number(t.get("FDP50")).unwrap_or(0.0)
        ),
    );
    row
}

pub(super) fn normalized_filter(filter: &str) -> String {
    if filter == "PASS" || filter == "." {
        String::new()
    } else {
        filter.to_string()
    }
}

pub(super) fn normalized_qual(qual: &str) -> String {
    if qual == "." {
        String::new()
    } else {
        qual.to_string()
    }
}

pub(super) fn write_simple_table(path: &Path, body: &str) -> Result<()> {
    fs::write(path, body).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}
