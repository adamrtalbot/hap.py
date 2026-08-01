//! Generic feature table — the default `hap ftx --feature-table generic`
//! output shape. Seven columns derived entirely from VCF fixed columns,
//! one row per input record.

use crate::strelka;
use crate::vcf::RawVcfRecord;

use super::common::{csv_escape, format_python_float, render_filter};

/// Emits the generic CSV as an ordered line buffer. The caller is
/// responsible for joining with `\n` and writing to disk.
///
/// `label` is the feature-label tag column. `targets` narrows the output
/// to records overlapping any interval when provided.
pub(super) fn emit(records: &[RawVcfRecord], label: &str) -> Vec<String> {
    emit_fields(
        records,
        label,
        &["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"],
    )
}

/// Legacy uses `GenericFeatures.collectFeatures` for truth-side `TP` and
/// `FN` inputs. Those tables have caller-specific schemas rather than the
/// query caller's feature schema.
pub(super) fn emit_fields(records: &[RawVcfRecord], label: &str, fields: &[&str]) -> Vec<String> {
    let mut lines = vec![format!(",{},tag", fields.join(","))];
    for (index, record) in records.iter().enumerate() {
        let mut cells = Vec::with_capacity(fields.len() + 2);
        cells.push(index.to_string());
        cells.extend(fields.iter().map(|field| render_field(record, field)));
        cells.push(csv_escape(label));
        lines.push(cells.join(","));
    }
    lines
}

fn render_field(record: &RawVcfRecord, field: &str) -> String {
    match field {
        "CHROM" => csv_escape(&record.chrom),
        "POS" => record.pos.to_string(),
        "REF" => csv_escape(&record.ref_allele),
        "ALT" => csv_escape(&record.alt_allele),
        "QUAL" => render_qual(&record.qual),
        "FILTER" => csv_escape(&render_filter(&record.filter)),
        field if field.starts_with("I.") => strelka::info_value(&record.info, &field[2..])
            .map(|value| render_inferred_value(&value))
            .unwrap_or_default(),
        field if field.starts_with("S.") => render_sample_field(record, field),
        literal => csv_escape(literal),
    }
}

fn render_sample_field(record: &RawVcfRecord, field: &str) -> String {
    let mut parts = field.split('.');
    let _ = parts.next();
    let Some(sample_number) = parts.next().and_then(|v| v.parse::<usize>().ok()) else {
        return String::new();
    };
    let Some(key) = parts.next() else {
        return String::new();
    };
    record
        .sample_map(sample_number.saturating_sub(1))
        .get(key)
        .map(|value| render_inferred_value(value))
        .unwrap_or_default()
}

fn render_inferred_value(value: &str) -> String {
    if value.contains(',') {
        return csv_escape(
            &value
                .split(',')
                .map(render_scalar_value)
                .collect::<Vec<_>>()
                .join(","),
        );
    }
    csv_escape(&render_scalar_value(value))
}

fn render_scalar_value(value: &str) -> String {
    if let Ok(value) = value.parse::<i64>() {
        value.to_string()
    } else if let Ok(value) = value.parse::<f64>() {
        format_python_float(value)
    } else {
        value.to_string()
    }
}

fn render_qual(value: &str) -> String {
    if value == "." || value.is_empty() {
        String::new()
    } else {
        value
            .parse::<f64>()
            .map(format_python_float)
            .unwrap_or_else(|_| value.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qual_uses_pandas_float_column_rendering() {
        assert_eq!(render_qual("60"), "60.0");
        assert_eq!(render_qual("47.25"), "47.25");
        assert_eq!(render_qual("."), "");
    }

    fn record(info: &str) -> RawVcfRecord {
        RawVcfRecord {
            chrom: "chr1".to_string(),
            pos: 7,
            id: ".".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "C,G".to_string(),
            qual: "60".to_string(),
            filter: "PASS".to_string(),
            info: info.to_string(),
            format: Some("GT".to_string()),
            samples: vec!["0/1".to_string(), "1/1".to_string()],
        }
    }

    #[test]
    fn caller_specific_truth_schema_extracts_info_and_sample_fields() {
        let lines = emit_fields(
            &[record("T_ALT_RATE=0.25;DP_normal=17;tag=truth")],
            "TP",
            &["CHROM", "ALT", "I.T_ALT_RATE", "I.DP_normal", "S.2.GT"],
        );
        assert_eq!(
            lines,
            vec![
                ",CHROM,ALT,I.T_ALT_RATE,I.DP_normal,S.2.GT,tag",
                "0,chr1,\"C,G\",0.25,17,1/1,TP",
            ]
        );
    }
}
