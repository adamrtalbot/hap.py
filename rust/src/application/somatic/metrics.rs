//! Metrics JSON and confidence calculations.

use anyhow::{Context, Result};
use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

pub(super) fn ratio(numerator: usize, denominator: usize) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        numerator as f64 / denominator as f64
    }
}

pub(super) fn py_float(value: f64) -> String {
    crate::adapters::report::full_repr_float(value)
}

pub(super) fn write_legacy_metrics_json(
    path: &Path,
    commandline: &str,
    stats_csv: &Path,
    ci_alpha: f64,
    ambiguous_classes: Option<&BTreeMap<String, usize>>,
    ambiguous_reasons: Option<&BTreeMap<String, usize>>,
) -> Result<()> {
    let text = fs::read_to_string(stats_csv)
        .with_context(|| format!("failed to read {}", stats_csv.display()))?;
    let mut lines = text.lines();
    let headers: Vec<String> = lines
        .next()
        .unwrap_or_default()
        .split(',')
        .map(str::to_string)
        .collect();
    let rows: Vec<Vec<String>> = lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split(',').map(str::to_string).collect())
        .collect();
    let metric_headers = headers[..headers.len().saturating_sub(2)].to_vec();
    let metric_rows: Vec<Vec<String>> = rows
        .iter()
        .map(|row| row[..row.len().saturating_sub(2)].to_vec())
        .collect();

    let mut body = String::new();
    body.push_str("{\"runInfo\": [{\"value\": ");
    body.push_str(&json_string(commandline));
    body.push_str(", \"key\": \"commandline\"}], \"metrics\": [");
    for (id, column, counts) in [
        ("ambiclasses", "class", ambiguous_classes),
        ("ambireasons", "reason", ambiguous_reasons),
    ] {
        if let Some(counts) = counts.filter(|counts| !counts.is_empty()) {
            body.push_str(&count_metric_json(id, column, counts));
            body.push_str(", ");
        }
    }
    body.push_str("{\"data\": [");
    body.push_str(&column_json(
        "types",
        "types",
        "string",
        &metric_rows.iter().map(|r| r[0].clone()).collect::<Vec<_>>(),
        true,
    ));
    for (i, header) in metric_headers.iter().enumerate().skip(1) {
        let exact = exact_somatic_metric_values(header, &metric_headers, &metric_rows, ci_alpha);
        let exact_metric = exact.is_some();
        let values = exact.unwrap_or_else(|| {
            metric_rows
                .iter()
                .map(|r| r.get(i).cloned().unwrap_or_default())
                .collect::<Vec<_>>()
        });
        let kind = if exact_metric {
            "double"
        } else {
            infer_type(&values)
        };
        body.push_str(", ");
        body.push_str(&column_json(header, header, kind, &values, false));
    }
    body.push_str(
        "], \"properties\": [], \"type\": \"Table\", \"id\": \"result\", \"label\": \"result\"}], ",
    );
    body.push_str("\"version\": \"\", \"sampleInfo\": [], \"name\": \"som.py.comparison\", \"parameters\": [], ");
    body.push_str("\"timestamp\": ");
    body.push_str(&json_string(&iso_timestamp_now()));
    body.push_str(", \"metadata\": {\"required\": {\"version\": \"\", \"id\": \"haplotypes\", \"module\": \"som.py\", \"description\": ");
    let executable = commandline.split_whitespace().next().unwrap_or("hap");
    body.push_str(&json_string(&format!(
        "{executable} generated this JSON file via command line {commandline}",
    )));
    body.push_str("}}}");
    fs::write(path, body).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

pub(super) fn exact_somatic_metric_values(
    column: &str,
    headers: &[String],
    rows: &[Vec<String>],
    ci_alpha: f64,
) -> Option<Vec<String>> {
    let metric = matches!(
        column,
        "recall"
            | "recall_lower"
            | "recall_upper"
            | "recall2"
            | "precision"
            | "precision_lower"
            | "precision_upper"
            | "na"
            | "ambiguous"
            | "fp.rate"
            | "recall.filtered"
            | "precision.filtered"
            | "fp.rate.filtered"
            | "na.filtered"
            | "ambiguous.filtered"
    );
    if !metric {
        return None;
    }

    let index = |name: &str| headers.iter().position(|header| header == name);
    let count = |row: &[String], name: &str| {
        index(name)
            .and_then(|column| row.get(column))
            .filter(|value| !value.is_empty() && value.as_str() != ".")
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| value.is_finite() && *value >= 0.0)
            .map(|value| value as usize)
    };
    let ratio = |numerator: usize, denominator: usize| {
        (denominator != 0).then(|| (numerator as f64 / denominator as f64).to_string())
    };
    let fp_rate = |false_positives: usize, region_size: usize| {
        if region_size == 0 {
            (false_positives != 0).then(|| "inf".to_string())
        } else {
            Some((1_000_000.0 * false_positives as f64 / region_size as f64).to_string())
        }
    };

    Some(
        rows.iter()
            .map(|row| {
                let Some(tp) = count(row, "tp") else {
                    return String::new();
                };
                let fp = count(row, "fp").unwrap_or(0);
                let fn_count = count(row, "fn").unwrap_or(0);
                let truth_total = count(row, "total.truth").unwrap_or(0);
                let query_total = count(row, "total.query").unwrap_or(0);
                let unk = count(row, "unk").unwrap_or(0);
                let ambi = count(row, "ambi").unwrap_or(0);
                let fp_region_size = count(row, "fp.region.size").unwrap_or(0);
                let (recall, recall_lower, recall_upper) = jeffreys_ci(tp, tp + fn_count, ci_alpha);
                let (precision, precision_lower, precision_upper) =
                    jeffreys_ci(tp, tp + fp, ci_alpha);
                let exact = match column {
                    "recall" => Some(recall.to_string()),
                    "recall_lower" => Some(recall_lower.to_string()),
                    "recall_upper" => Some(recall_upper.to_string()),
                    "recall2" => ratio(tp, truth_total),
                    "precision" => Some(precision.to_string()),
                    "precision_lower" => Some(precision_lower.to_string()),
                    "precision_upper" => Some(precision_upper.to_string()),
                    "na" => ratio(unk, query_total),
                    "ambiguous" => ratio(ambi, query_total),
                    "fp.rate" => fp_rate(fp, fp_region_size),
                    filtered => {
                        let Some(filtered_tp) = count(row, "tp.filtered") else {
                            return String::new();
                        };
                        let filtered_fp = count(row, "fp.filtered").unwrap_or(0);
                        let filtered_unk = count(row, "unk.filtered").unwrap_or(0);
                        let filtered_ambi = count(row, "ambi.filtered").unwrap_or(0);
                        let unfiltered_tp = tp.saturating_sub(filtered_tp);
                        let unfiltered_fp = fp.saturating_sub(filtered_fp);
                        match filtered {
                            "recall.filtered" => ratio(unfiltered_tp, tp + fn_count),
                            "precision.filtered" => {
                                ratio(unfiltered_tp, unfiltered_tp + unfiltered_fp)
                            }
                            "fp.rate.filtered" => fp_rate(unfiltered_fp, fp_region_size),
                            "na.filtered" => ratio(unk.saturating_sub(filtered_unk), query_total),
                            "ambiguous.filtered" => {
                                ratio(ambi.saturating_sub(filtered_ambi), query_total)
                            }
                            _ => None,
                        }
                    }
                };
                exact.unwrap_or_default()
            })
            .collect(),
    )
}

pub(super) fn count_metric_json(
    id: &str,
    column: &str,
    counts: &BTreeMap<String, usize>,
) -> String {
    let legacy_indices = python2_counter_indices(counts);
    let indices = counts
        .keys()
        .map(|key| legacy_indices[key].to_string())
        .collect::<Vec<_>>();
    let labels = counts.keys().cloned().collect::<Vec<_>>();
    let values = counts
        .values()
        .map(|count| count.to_string())
        .collect::<Vec<_>>();
    format!(
        "{{\"data\": [{}, {}, {}], \"properties\": [], \"type\": \"Table\", \"id\": {}, \"label\": {}}}",
        column_json("types", "types", "string", &indices, true),
        column_json(column, column, "string", &labels, false),
        column_json("count", "count", "int64", &values, false),
        json_string(id),
        json_string(id),
    )
}

pub(super) fn python2_counter_indices(counts: &BTreeMap<String, usize>) -> BTreeMap<String, usize> {
    fn hash(value: &str) -> u64 {
        let bytes = value.as_bytes();
        let mut hash = bytes.first().copied().unwrap_or_default() as u64 * 128;
        for byte in bytes {
            hash = hash.wrapping_mul(1_000_003) ^ u64::from(*byte);
        }
        hash ^= bytes.len() as u64;
        if hash == u64::MAX { u64::MAX - 1 } else { hash }
    }

    fn insert(table: &mut [Option<String>], key: String) {
        let mask = table.len() - 1;
        let hashed = hash(&key);
        let mut slot = hashed as usize & mask;
        let mut perturb = hashed;
        while table[slot].is_some() {
            slot = slot
                .wrapping_mul(5)
                .wrapping_add(perturb as usize)
                .wrapping_add(1)
                & mask;
            perturb >>= 5;
        }
        table[slot] = Some(key);
    }

    let mut table = vec![None; 8];
    let mut used = 0usize;
    for key in counts.keys() {
        insert(&mut table, key.clone());
        used += 1;
        if used * 3 >= table.len() * 2 {
            let minimum = if used > 50_000 { used * 2 } else { used * 4 };
            let mut size = 8usize;
            while size <= minimum {
                size *= 2;
            }
            let old = std::mem::replace(&mut table, vec![None; size]);
            for key in old.into_iter().flatten() {
                insert(&mut table, key);
            }
        }
    }

    table
        .into_iter()
        .flatten()
        .enumerate()
        .map(|(index, key)| (key, index))
        .collect()
}

pub(super) fn infer_type(values: &[String]) -> &'static str {
    let is_missing = |value: &str| value.is_empty() || value == ".";
    let populated = values
        .iter()
        .filter(|value| !is_missing(value))
        .collect::<Vec<_>>();
    let has_missing = populated.len() != values.len();
    if !populated.is_empty() && populated.iter().all(|value| value.parse::<i64>().is_ok()) {
        // pandas promotes integer columns containing NaN to float64.
        if has_missing { "double" } else { "int64" }
    } else if !populated.is_empty() && populated.iter().all(|value| value.parse::<f64>().is_ok()) {
        "double"
    } else {
        "string"
    }
}

pub(super) fn column_json(
    id: &str,
    label: &str,
    kind: &str,
    values: &[String],
    numeric_strings: bool,
) -> String {
    let rendered = values
        .iter()
        .map(|value| match kind {
            "int64" | "double" if value.is_empty() || value == "." => "null".to_string(),
            "int64" | "double" if value.parse::<f64>().is_ok_and(|number| !number.is_finite()) => {
                "null".to_string()
            }
            // json.dumps serializes the numeric binary64 value, not the CSV
            // cell's quoted text, so parse the cell back to f64 and emit the
            // number here.
            "double" => value
                .parse::<f64>()
                .map(json_float)
                .unwrap_or_else(|_| value.to_string()),
            "int64" => value.to_string(),
            _ if numeric_strings && value.parse::<i64>().is_ok() => value.to_string(),
            _ => json_string(value),
        })
        .collect::<Vec<_>>()
        .join(", ");
    format!(
        "{{\"values\": [{rendered}], \"type\": {}, \"id\": {}, \"label\": {}}}",
        json_string(kind),
        json_string(id),
        json_string(label)
    )
}

pub(super) fn json_float(value: f64) -> String {
    let mut rendered = value.to_string();
    if !rendered.contains(['.', 'e', 'E']) {
        rendered.push_str(".0");
    }
    rendered
}

pub(super) fn json_string(value: &str) -> String {
    let escaped = value
        .replace('\\', "\\\\")
        .replace('"', "\\\"")
        .replace('\n', "\\n");
    format!("\"{escaped}\"")
}

/// Current UTC time formatted the same way Python's `datetime.datetime.now().isoformat()`
/// renders it (microsecond precision, no timezone suffix) — matches legacy `som.py` JSON.
pub(super) fn iso_timestamp_now() -> String {
    let duration = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default();
    let secs = duration.as_secs() as i64;
    let micros = duration.subsec_micros();

    let days = secs.div_euclid(86_400);
    let time_of_day = secs.rem_euclid(86_400);
    let hour = time_of_day / 3_600;
    let minute = (time_of_day / 60) % 60;
    let second = time_of_day % 60;

    // Civil-from-days algorithm (Howard Hinnant, public domain).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1_460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = if m <= 2 { y + 1 } else { y };

    format!("{year:04}-{m:02}-{d:02}T{hour:02}:{minute:02}:{second:02}.{micros:06}")
}

pub(super) fn jeffreys_ci(x: usize, n: usize, alpha: f64) -> (f64, f64, f64) {
    if n == 0 {
        return (0.0, 0.0, 1.0);
    }
    let p = x as f64 / n as f64;
    let lower = if x == n {
        (alpha / 2.0).powf(1.0 / n as f64)
    } else if x <= 1 {
        0.0
    } else {
        beta_ppf(alpha / 2.0, x as f64 + 0.5, (n - x) as f64 + 0.5)
    };
    let upper = if x == 0 {
        1.0 - (alpha / 2.0).powf(1.0 / n as f64)
    } else if x >= n.saturating_sub(1) {
        1.0
    } else {
        beta_isf(alpha / 2.0, x as f64 + 0.5, (n - x) as f64 + 0.5)
    };
    (p, lower, upper)
}

pub(super) fn beta_isf(q: f64, a: f64, b: f64) -> f64 {
    beta_ppf(1.0 - q, a, b)
}

pub(super) fn beta_ppf(p: f64, a: f64, b: f64) -> f64 {
    // Bit-exact port of scipy 1.2.1's beta.ppf via cephes_incbi.
    // Used for `Tools/ci.py::jeffreys` confidence intervals so sompy
    // recall_lower / precision_lower / etc. match legacy hap.py output.
    crate::domain::cephes::incbi(a, b, p)
}
