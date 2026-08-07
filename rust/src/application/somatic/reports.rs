//! Cohesive responsibility extracted from the command façade.

use super::allele_frequency::{parse_af_bins, round_four, rounded_metric};
use super::features::{csv_join, parse_csv_line, write_simple_table};
use super::metrics::{jeffreys_ci, py_float, ratio};
use super::{
    AmbiguousInterval, FilteredCounts, FilteredRawRecord, SOM_VERSION, SOMATIC_ROC_CHUNK,
    SOMATIC_ROC_MERGE_FAN_IN, SomaticCounts, StatsRowContext, somatic_roc_config,
};
use crate::adapters::vcf;
use crate::domain::{Interval, RawVcfRecord};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

pub(super) fn write_happy_style_summary(
    path: &Path,
    feature_header: &str,
    feature_rows: &Path,
    feature_table: &str,
) -> Result<()> {
    let mut lines = vec![
        ",Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio".to_string()
    ];
    let headers = parse_csv_line(feature_header);
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = csv_column_index(&headers, "FILTER")?;
    let ref_index = csv_column_index(&headers, "REF")?;
    let alt_index = csv_column_index(&headers, "ALT")?;
    let truth_ref_index = csv_column_index(&headers, "REF.truth")?;
    let happy_type = match feature_table.rsplit('.').next() {
        Some("snv") => "SNP",
        Some("indel") => "INDEL",
        _ => "NA",
    };
    let mut truth_total = 0usize;
    let mut truth_tp = 0usize;
    let mut truth_fn = 0usize;
    let mut query = [(0usize, 0usize); 2];
    for row in BufReader::new(File::open(feature_rows)?).lines() {
        let row = parse_csv_line(&row?);
        if nonempty_csv_field(&row, truth_ref_index) {
            truth_total += 1;
            match row.get(tag_index).map(String::as_str) {
                Some("TP") => truth_tp += 1,
                Some("FN") => truth_fn += 1,
                _ => {}
            }
        }
        if nonempty_csv_field(&row, ref_index) && nonempty_csv_field(&row, alt_index) {
            let pass = row
                .get(filter_index)
                .is_none_or(|value| value.is_empty() || value == "." || value == "PASS");
            for (index, include) in [pass, true].into_iter().enumerate() {
                if include {
                    match row.get(tag_index).map(String::as_str) {
                        Some("FP") => query[index].0 += 1,
                        Some("UNK" | "AMBI") => query[index].1 += 1,
                        _ => {}
                    }
                }
            }
        }
    }

    for (filter_index, filter) in ["PASS", "ALL"].into_iter().enumerate() {
        let (query_fp, query_unk) = query[filter_index];
        let query_total = truth_tp + query_fp + query_unk;
        let recall = rounded_metric(truth_tp, truth_total);
        let precision = rounded_metric(truth_tp, truth_tp + query_fp);
        let frac_na = rounded_metric(query_unk, query_total);
        let f1 = match (recall, precision) {
            (Some(recall), Some(precision)) if recall + precision > 0.0 => {
                Some(round_four(2.0 * recall * precision / (recall + precision)))
            }
            _ => None,
        };
        lines.push(csv_join([
            "0".to_string(),
            happy_type.to_string(),
            filter.to_string(),
            truth_total.to_string(),
            truth_tp.to_string(),
            truth_fn.to_string(),
            query_total.to_string(),
            query_fp.to_string(),
            query_unk.to_string(),
            "NA".to_string(),
            render_summary_metric(recall),
            render_summary_metric(precision),
            render_summary_metric(frac_na),
            render_summary_metric(f1),
            "NA".to_string(),
            "NA".to_string(),
            "NA".to_string(),
            "NA".to_string(),
        ]));
    }
    write_simple_table(path, &format!("{}\n", lines.join("\n")))
}

pub(super) fn render_summary_metric(value: Option<f64>) -> String {
    value.map_or_else(String::new, py_float)
}

pub(super) fn write_happy_style_extended(
    path: &Path,
    feature_header: &str,
    feature_rows: &Path,
    feature_table: &str,
    bin_sizes: &str,
    truth_af_field: &str,
    query_af_field: &str,
) -> Result<()> {
    let mut lines = vec![
        "Type,Subtype,Subset,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio".to_string()
    ];
    let headers = parse_csv_line(feature_header);
    let truth_af_index = headers
        .iter()
        .position(|header| header == truth_af_field)
        .with_context(|| {
            format!("truth AF feature '{truth_af_field}' is not in the feature table")
        })?;
    let query_af_index = headers
        .iter()
        .position(|header| header == query_af_field)
        .with_context(|| {
            format!("query AF feature '{query_af_field}' is not in the feature table")
        })?;
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = csv_column_index(&headers, "FILTER")?;
    let ref_index = csv_column_index(&headers, "REF")?;
    let alt_index = csv_column_index(&headers, "ALT")?;
    let truth_ref_index = csv_column_index(&headers, "REF.truth")?;
    let happy_type = match feature_table.rsplit('.').next() {
        Some("snv") => "SNP",
        Some("indel") => "INDEL",
        _ => "NA",
    };

    for (start, end) in parse_af_bins(bin_sizes)? {
        let inclusive_last = end >= 1.0;
        let subset = if inclusive_last {
            format!("[{start:.2},1.00]")
        } else {
            let end = if end.is_nan() {
                "nan".to_string()
            } else {
                format!("{end:.2}")
            };
            format!("[{start:.2},{end})")
        };
        let in_bin = |value: &str| {
            value.parse::<f64>().is_ok_and(|value| {
                value >= start && (value < end || (inclusive_last && value <= 1.0))
            })
        };
        let mut truth_total = 0usize;
        let mut truth_tp = 0usize;
        let mut truth_fn = 0usize;
        let mut query = [(0usize, 0usize); 2];
        for row in BufReader::new(File::open(feature_rows)?).lines() {
            let row = parse_csv_line(&row?);
            if nonempty_csv_field(&row, truth_ref_index)
                && row.get(truth_af_index).is_some_and(|value| in_bin(value))
            {
                truth_total += 1;
                match row.get(tag_index).map(String::as_str) {
                    Some("TP") => truth_tp += 1,
                    Some("FN") => truth_fn += 1,
                    _ => {}
                }
            }
            if nonempty_csv_field(&row, ref_index)
                && nonempty_csv_field(&row, alt_index)
                && row.get(query_af_index).is_some_and(|value| in_bin(value))
            {
                let pass = row
                    .get(filter_index)
                    .is_none_or(|value| value.is_empty() || value == "." || value == "PASS");
                for (index, include) in [pass, true].into_iter().enumerate() {
                    if include {
                        match row.get(tag_index).map(String::as_str) {
                            Some("FP") => query[index].0 += 1,
                            Some("UNK" | "AMBI") => query[index].1 += 1,
                            _ => {}
                        }
                    }
                }
            }
        }

        // Preserve the legacy bin-major PASS, ALL append order.
        for (filter_index, filter) in ["PASS", "ALL"].into_iter().enumerate() {
            let (query_fp, query_unk) = query[filter_index];
            let query_total = truth_tp + query_fp + query_unk;
            let recall = rounded_metric(truth_tp, truth_total);
            let precision = rounded_metric(truth_tp, truth_tp + query_fp);
            let frac_na = rounded_metric(query_unk, query_total);
            let f1 = match (recall, precision) {
                (Some(recall), Some(precision)) if recall + precision > 0.0 => {
                    Some(round_four(2.0 * recall * precision / (recall + precision)))
                }
                _ => None,
            };
            lines.push(csv_join([
                happy_type.to_string(),
                "*".to_string(),
                subset.clone(),
                filter.to_string(),
                truth_total.to_string(),
                truth_tp.to_string(),
                truth_fn.to_string(),
                query_total.to_string(),
                query_fp.to_string(),
                query_unk.to_string(),
                "NA".to_string(),
                render_extended_metric(recall),
                render_extended_metric(precision),
                render_extended_metric(frac_na),
                render_extended_metric(f1),
                "NA".to_string(),
                "NA".to_string(),
                "NA".to_string(),
                "NA".to_string(),
            ]));
        }
    }
    write_simple_table(path, &format!("{}\n", lines.join("\n")))
}

pub(super) fn render_extended_metric(value: Option<f64>) -> String {
    value.map_or_else(|| "NA".to_string(), py_float)
}

pub(super) fn csv_column_index(headers: &[String], name: &str) -> Result<usize> {
    headers
        .iter()
        .position(|header| header == name)
        .with_context(|| format!("feature table is missing required column '{name}'"))
}

pub(super) fn nonempty_csv_field(row: &[String], index: usize) -> bool {
    row.get(index)
        .is_some_and(|value| !value.is_empty() && value != ".")
}

pub(super) fn calculate_af_stats(
    feature_header: &str,
    feature_rows: &Path,
    bin_sizes: &str,
    truth_af_field: &str,
    query_af_field: &str,
    type_label: Option<&str>,
) -> Result<Vec<(f64, f64, SomaticCounts, FilteredCounts)>> {
    let headers = parse_csv_line(feature_header);
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = csv_column_index(&headers, "FILTER")?;
    let truth_af_index = csv_column_index(&headers, truth_af_field)?;
    let query_af_index = csv_column_index(&headers, query_af_field)?;
    let mut output = Vec::new();
    for (start, end) in parse_af_bins(bin_sizes)? {
        let in_bin = |row: &[String], index: usize| {
            row.get(index)
                .and_then(|value| value.parse::<f64>().ok())
                .is_some_and(|value| value >= start && value <= 1.0 && value < end)
                || (end >= 1.0
                    && row
                        .get(index)
                        .and_then(|value| value.parse::<f64>().ok())
                        .is_some_and(|value| value == 1.0))
        };

        let mut counts = SomaticCounts::default();
        let mut filtered = FilteredCounts::default();
        for row in BufReader::new(File::open(feature_rows)?).lines() {
            let row = parse_csv_line(&row?);
            if type_label.is_some_and(|expected| feature_row_type(&headers, &row) != Some(expected))
            {
                continue;
            }
            let tag = row.get(tag_index).map(String::as_str).unwrap_or_default();
            let filtered_call = row.get(filter_index).is_some_and(|value| !value.is_empty());
            match tag {
                "TP" if in_bin(&row, truth_af_index) => {
                    counts.tp += 1;
                    if filtered_call {
                        filtered.tp += 1;
                    }
                }
                "FN" if in_bin(&row, truth_af_index) => counts.fn_count += 1,
                "FP" if in_bin(&row, query_af_index) => {
                    counts.fp += 1;
                    if filtered_call {
                        filtered.fp += 1;
                    }
                }
                "UNK" if in_bin(&row, query_af_index) => {
                    counts.unk += 1;
                    if filtered_call {
                        filtered.unk += 1;
                    }
                }
                "AMBI" if in_bin(&row, query_af_index) => {
                    counts.ambi += 1;
                    if filtered_call {
                        filtered.ambi += 1;
                    }
                }
                _ => {}
            }
        }
        counts.truth_total = counts.tp + counts.fn_count;
        counts.query_total = counts.tp + counts.fp + counts.unk + counts.ambi;
        output.push((start, end, counts, filtered));
    }
    Ok(output)
}

pub(super) fn feature_rows_for_type(
    feature_rows: &Path,
    feature_header: &str,
    type_label: Option<&str>,
) -> Result<Option<tempfile::NamedTempFile>> {
    let headers = parse_csv_line(feature_header);
    let mut output =
        tempfile::NamedTempFile::new().context("failed to create feature type spool")?;
    let mut count = 0usize;
    for row in BufReader::new(File::open(feature_rows)?).lines() {
        let row = row?;
        if type_label.is_none() || feature_row_type(&headers, &parse_csv_line(&row)) == type_label {
            writeln!(output.as_file_mut(), "{row}")?;
            count += 1;
        }
    }
    output.as_file_mut().flush()?;
    Ok((count > 0).then_some(output))
}

pub(super) fn feature_row_type(headers: &[String], row: &[String]) -> Option<&'static str> {
    let tag_index = csv_column_index(headers, "tag").ok()?;
    let tag = row.get(tag_index)?.as_str();
    let (reference_field, alternate_field) = if matches!(tag, "TP" | "FN") {
        ("REF.truth", "ALT.truth")
    } else {
        ("REF", "ALT")
    };
    let reference = row.get(csv_column_index(headers, reference_field).ok()?)?;
    let alternate = row.get(csv_column_index(headers, alternate_field).ok()?)?;
    let (reference, alternate) = if reference.is_empty() || alternate.is_empty() {
        (
            row.get(csv_column_index(headers, "REF").ok()?)?,
            row.get(csv_column_index(headers, "ALT").ok()?)?,
        )
    } else {
        (reference, alternate)
    };
    feature_allele_type(reference, alternate)
}

pub(super) fn feature_allele_type(reference: &str, alternate: &str) -> Option<&'static str> {
    if alternate.is_empty() || alternate == "." {
        return None;
    }
    if alternate == "*" || alternate.starts_with('<') || alternate.contains(['[', ']']) {
        return Some("others");
    }
    if reference.len() == 1 && alternate.len() == 1 {
        Some("SNVs")
    } else if reference.len() > 1 && alternate.len() == reference.len() {
        Some("MNPs")
    } else {
        Some("indels")
    }
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(super) enum SomaticRocTag {
    Tp,
    Fp,
    Fn,
}

impl SomaticRocTag {
    fn code(self) -> u8 {
        match self {
            Self::Tp => 0,
            Self::Fp => 1,
            Self::Fn => 2,
        }
    }

    fn from_code(code: u8) -> Result<Self> {
        match code {
            0 => Ok(Self::Tp),
            1 => Ok(Self::Fp),
            2 => Ok(Self::Fn),
            _ => bail!("invalid somatic ROC tag code {code}"),
        }
    }
}

struct SomaticRocSpool {
    buffer: Vec<(f64, SomaticRocTag, u64)>,
    chunks: Vec<tempfile::TempPath>,
    next_serial: u64,
    totals: [usize; 3],
}

impl SomaticRocSpool {
    fn new() -> Self {
        Self {
            buffer: Vec::new(),
            chunks: Vec::new(),
            next_serial: 0,
            totals: [0; 3],
        }
    }

    fn push(&mut self, score: f64, tag: SomaticRocTag) -> Result<()> {
        self.totals[tag.code() as usize] += 1;
        self.buffer.push((score, tag, self.next_serial));
        self.next_serial += 1;
        if self.buffer.len() >= SOMATIC_ROC_CHUNK {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer.sort_by(|left, right| {
            left.0
                .total_cmp(&right.0)
                .then(left.1.cmp(&right.1))
                .then(left.2.cmp(&right.2))
        });
        let mut chunk =
            tempfile::NamedTempFile::new().context("failed to create somatic ROC chunk")?;
        {
            let mut writer = BufWriter::new(chunk.as_file_mut());
            for (score, tag, serial) in self.buffer.drain(..) {
                writeln!(writer, "{}\t{}\t{serial}", score.to_bits(), tag.code())?;
            }
            writer.flush()?;
        }
        self.chunks.push(chunk.into_temp_path());
        Ok(())
    }

    fn finish(mut self) -> Result<(SomaticRocMerge, [usize; 3])> {
        self.flush()?;
        let chunks = collapse_somatic_roc_chunks(self.chunks)?;
        Ok((SomaticRocMerge::open(chunks)?, self.totals))
    }
}

struct SomaticRocMerge {
    _chunks: Vec<tempfile::TempPath>,
    readers: Vec<std::io::Lines<BufReader<File>>>,
    current: Vec<Option<(f64, SomaticRocTag, u64)>>,
    heap: std::collections::BinaryHeap<std::cmp::Reverse<(u64, u8, u64, usize)>>,
}

impl SomaticRocMerge {
    fn open(chunks: Vec<tempfile::TempPath>) -> Result<Self> {
        let mut readers = Vec::with_capacity(chunks.len());
        for chunk in &chunks {
            readers.push(BufReader::new(File::open(chunk)?).lines());
        }
        let mut merge = Self {
            current: (0..readers.len()).map(|_| None).collect(),
            readers,
            heap: std::collections::BinaryHeap::new(),
            _chunks: chunks,
        };
        for index in 0..merge.readers.len() {
            merge.advance(index)?;
        }
        Ok(merge)
    }

    fn advance(&mut self, index: usize) -> Result<()> {
        let Some(line) = self.readers[index].next() else {
            return Ok(());
        };
        let line = line?;
        let mut fields = line.split('\t');
        let score = f64::from_bits(
            fields
                .next()
                .context("somatic ROC score missing")?
                .parse()?,
        );
        let tag =
            SomaticRocTag::from_code(fields.next().context("somatic ROC tag missing")?.parse()?)?;
        let serial = fields
            .next()
            .context("somatic ROC serial missing")?
            .parse()?;
        self.heap.push(std::cmp::Reverse((
            somatic_sortable_f64(score),
            tag.code(),
            serial,
            index,
        )));
        self.current[index] = Some((score, tag, serial));
        Ok(())
    }
}

impl Iterator for SomaticRocMerge {
    type Item = Result<(f64, SomaticRocTag)>;

    fn next(&mut self) -> Option<Self::Item> {
        let std::cmp::Reverse((_, _, _, index)) = self.heap.pop()?;
        let (score, tag, _) = self.current[index]
            .take()
            .expect("somatic ROC merge entry has a row");
        if let Err(error) = self.advance(index) {
            return Some(Err(error));
        }
        Some(Ok((score, tag)))
    }
}

fn somatic_sortable_f64(value: f64) -> u64 {
    let bits = value.to_bits();
    if bits & (1 << 63) == 0 {
        bits ^ (1 << 63)
    } else {
        !bits
    }
}

fn collapse_somatic_roc_chunks(
    mut chunks: Vec<tempfile::TempPath>,
) -> Result<Vec<tempfile::TempPath>> {
    while chunks.len() > SOMATIC_ROC_MERGE_FAN_IN {
        let mut merged = Vec::with_capacity(chunks.len().div_ceil(SOMATIC_ROC_MERGE_FAN_IN));
        let mut remaining = chunks.into_iter();
        loop {
            let batch = remaining
                .by_ref()
                .take(SOMATIC_ROC_MERGE_FAN_IN)
                .collect::<Vec<_>>();
            if batch.is_empty() {
                break;
            }
            let mut output = tempfile::NamedTempFile::new()
                .context("failed to create merged somatic ROC chunk")?;
            {
                let mut writer = BufWriter::new(output.as_file_mut());
                let merge = SomaticRocMerge::open(batch)?;
                for (serial, row) in merge.enumerate() {
                    let (score, tag) = row?;
                    writeln!(writer, "{}\t{}\t{serial}", score.to_bits(), tag.code())?;
                }
                writer.flush()?;
            }
            merged.push(output.into_temp_path());
        }
        chunks = merged;
    }
    Ok(chunks)
}

pub(super) fn write_somatic_roc(
    path: &Path,
    feature_header: &str,
    feature_rows: &Path,
    roc_name: &str,
) -> Result<()> {
    let config = somatic_roc_config(roc_name)
        .with_context(|| format!("unsupported somatic ROC mode '{roc_name}'"))?;
    let headers = parse_csv_line(feature_header);
    let score_index = csv_column_index(&headers, config.score)?;
    let tag_index = csv_column_index(&headers, "tag")?;
    let filter_index = headers.iter().position(|field| field == "FILTER");
    let nt_index = headers.iter().position(|field| field == "NT");
    let mut observations = SomaticRocSpool::new();

    for line in BufReader::new(File::open(feature_rows)?).lines() {
        let line = line?;
        let fields = parse_csv_line(&line);
        let Some(tag) = fields.get(tag_index).and_then(|value| {
            let lower = value.to_ascii_lowercase();
            if lower.starts_with("tp") {
                Some(SomaticRocTag::Tp)
            } else if lower.starts_with("fp") {
                Some(SomaticRocTag::Fp)
            } else if lower.starts_with("fn") {
                Some(SomaticRocTag::Fn)
            } else {
                None
            }
        }) else {
            continue;
        };
        let mut score = fields
            .get(score_index)
            .and_then(|value| value.parse::<f64>().ok())
            .unwrap_or(0.0);
        if config.zero_score_unless_nt_ref
            && nt_index
                .and_then(|index| fields.get(index))
                .is_some_and(|value| value != "ref")
        {
            score = 0.0;
        }
        if let (Some(index), Some(filter_name)) = (filter_index, config.filter_name) {
            let mut filters = fields
                .get(index)
                .map(String::as_str)
                .unwrap_or_default()
                .split([';', ','])
                .filter(|value| !value.is_empty() && *value != "." && *value != "PASS")
                .collect::<Vec<_>>();
            filters.retain(|value| *value != filter_name);
            if !filters.is_empty() {
                score = f64::MIN_POSITIVE;
            }
        }
        observations.push(score, tag)?;
    }
    let (observations, totals) = observations.finish()?;
    let mut tp = totals[SomaticRocTag::Tp.code() as usize];
    let mut fp = totals[SomaticRocTag::Fp.code() as usize];
    let mut fn_count = totals[SomaticRocTag::Fn.code() as usize];
    let mut roc_rows =
        tempfile::NamedTempFile::new().context("failed to create somatic ROC output spool")?;
    let mut integer_scores = true;
    let mut integer_precision = true;
    let mut integer_recall = true;
    let mut previous = None;
    for observation in observations {
        let (score, tag) = observation?;
        if previous != Some(score) {
            let precision = if tp + fp == 0 {
                1.0
            } else {
                tp as f64 / (tp + fp) as f64
            };
            let recall = if tp + fn_count == 0 {
                0.0
            } else {
                tp as f64 / (tp + fn_count) as f64
            };
            previous = Some(score);
            let score = cpp_default_six(score);
            let precision = cpp_default_six(precision);
            let recall = cpp_default_six(recall);
            integer_scores &= score.parse::<i64>().is_ok();
            integer_precision &= precision.parse::<i64>().is_ok();
            integer_recall &= recall.parse::<i64>().is_ok();
            writeln!(
                roc_rows.as_file_mut(),
                "{score}\t{tp}\t{fp}\t{fn_count}\t{precision}\t{recall}"
            )?;
        }
        match tag {
            SomaticRocTag::Tp => {
                tp = tp.saturating_sub(1);
                fn_count += 1;
            }
            SomaticRocTag::Fp => fp = fp.saturating_sub(1),
            SomaticRocTag::Fn => {}
        }
    }
    roc_rows.as_file_mut().flush()?;
    let mut output = BufWriter::new(
        File::create(path).with_context(|| format!("failed to create {}", path.display()))?,
    );
    writeln!(output, ",{},tp,fp,fn,precision,recall", config.score)?;
    for (row_index, row) in BufReader::new(File::open(roc_rows.path())?)
        .lines()
        .enumerate()
    {
        let row = row?;
        let mut fields = row.split('\t');
        let score = fields
            .next()
            .context("somatic ROC score missing")?
            .to_string();
        let tp = fields.next().context("somatic ROC TP missing")?;
        let fp = fields.next().context("somatic ROC FP missing")?;
        let fn_count = fields.next().context("somatic ROC FN missing")?;
        let precision = fields
            .next()
            .context("somatic ROC precision missing")?
            .to_string();
        let recall = fields
            .next()
            .context("somatic ROC recall missing")?
            .to_string();
        let score = if integer_scores {
            score
        } else {
            score
                .parse::<f64>()
                .map(|value| format!("{value:.8}"))
                .unwrap_or(score)
        };
        let precision = if integer_precision {
            precision
        } else {
            precision
                .parse::<f64>()
                .map(|value| format!("{value:.8}"))
                .unwrap_or(precision)
        };
        let recall = if integer_recall {
            recall
        } else {
            recall
                .parse::<f64>()
                .map(|value| format!("{value:.8}"))
                .unwrap_or(recall)
        };
        writeln!(
            output,
            "{row_index},{score},{tp},{fp},{fn_count},{precision},{recall}"
        )?;
    }
    output.flush()?;
    Ok(())
}

pub(super) fn cpp_default_six(value: f64) -> String {
    if value == 0.0 {
        return "0".to_string();
    }
    let exponent = value.abs().log10().floor() as i32;
    if !(-4..6).contains(&exponent) {
        let scientific = format!("{value:.5e}");
        let Some((mantissa, exponent)) = scientific.split_once('e') else {
            return scientific;
        };
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        return format!("{mantissa}e{exponent}");
    }
    let decimals = (5 - exponent).max(0) as usize;
    let fixed = format!("{value:.decimals$}");
    if fixed.contains('.') {
        fixed
            .trim_end_matches('0')
            .trim_end_matches('.')
            .to_string()
    } else {
        fixed
    }
}

pub(super) fn feature_rows_for_af_roc(
    feature_header: &str,
    feature_rows: &Path,
    start: f64,
    end: f64,
    truth_af_field: &str,
    query_af_field: &str,
) -> Result<Option<tempfile::NamedTempFile>> {
    let headers = parse_csv_line(feature_header);
    let tag_index = csv_column_index(&headers, "tag")?;
    let truth_index = csv_column_index(&headers, truth_af_field)?;
    let query_index = csv_column_index(&headers, query_af_field)?;
    let in_bin = |fields: &[String], index: usize| {
        fields
            .get(index)
            .and_then(|value| value.parse::<f64>().ok())
            .is_some_and(|value| value >= start && (value < end || (end >= 1.0 && value == 1.0)))
    };
    let mut output = tempfile::NamedTempFile::new().context("failed to create feature AF spool")?;
    let mut count = 0usize;
    for line in BufReader::new(File::open(feature_rows)?).lines() {
        let line = line?;
        let fields = parse_csv_line(&line);
        let keep = match fields.get(tag_index).map(String::as_str) {
            Some("TP" | "FN") => in_bin(&fields, truth_index),
            Some("FP") => in_bin(&fields, query_index),
            _ => false,
        };
        if keep {
            writeln!(output.as_file_mut(), "{line}")?;
            count += 1;
        }
    }
    output.as_file_mut().flush()?;
    Ok((count > 0).then_some(output))
}

pub(super) fn raw_type_label(record: &RawVcfRecord) -> Option<&'static str> {
    let alternates = record.alt_allele.split(',').collect::<Vec<_>>();
    if alternates.iter().all(|alternate| *alternate == ".") {
        return None;
    }
    if alternates.iter().any(|alternate| {
        alternate == &"*" || alternate.starts_with('<') || alternate.contains(['[', ']'])
    }) {
        return Some("others");
    }

    let reference_len = record.ref_allele.len();
    if reference_len == 1 && alternates.iter().all(|alternate| alternate.len() == 1) {
        Some("SNVs")
    } else if reference_len > 1
        && alternates
            .iter()
            .all(|alternate| alternate.len() == reference_len)
    {
        Some("MNPs")
    } else {
        Some("indels")
    }
}

#[cfg(test)]
pub(super) fn contigs_in_truth(truth: &[FilteredRawRecord]) -> BTreeSet<String> {
    truth
        .iter()
        .map(|record| record.key.chrom.clone())
        .collect()
}

pub(super) fn automatic_fp_intervals<'a>(
    fp_regions: &'a [Interval],
    ambiguous_regions: &'a [AmbiguousInterval],
) -> impl Iterator<Item = &'a Interval> {
    fp_regions.iter().chain(
        ambiguous_regions
            .iter()
            .filter(|entry| entry.label == "FP")
            .map(|entry| &entry.interval),
    )
}

pub(super) fn has_automatic_fp_bases(
    fp_regions: &[Interval],
    ambiguous_regions: &[AmbiguousInterval],
) -> bool {
    automatic_fp_intervals(fp_regions, ambiguous_regions)
        .any(|interval| interval.end > interval.start)
}

pub(super) fn fp_region_size_requires_reference(
    requested: Option<&str>,
    fp_regions: &[Interval],
    ambiguous_regions: &[AmbiguousInterval],
) -> bool {
    requested
        .and_then(|value| value.parse::<usize>().ok())
        .is_none()
        && !has_automatic_fp_bases(fp_regions, ambiguous_regions)
}

pub(super) fn validate_legacy_fp_location_denominator(
    requested_size: Option<&str>,
    location: Option<&str>,
    has_automatic_fp_bases: bool,
) -> Result<()> {
    if !has_automatic_fp_bases
        || requested_size
            .and_then(|size| size.parse::<i64>().ok())
            .is_some()
    {
        return Ok(());
    }
    let Some(location) = location.filter(|location| !location.is_empty()) else {
        return Ok(());
    };
    let Some((_, rest)) = location.split_once(':') else {
        return Ok(());
    };
    if rest.is_empty() {
        return Ok(());
    }

    // This deliberately mirrors som.py's `partition("_")` typo. Bcftools
    // accepts `chr:start-end`, then the denominator code tries to parse the
    // entire `start-end` token as an integer and the run exits after writing
    // any requested feature/ROC artifacts but before stats or metrics.
    let (start, end) = rest.split_once('_').unwrap_or((rest, ""));
    for bound in [start, end].into_iter().filter(|bound| !bound.is_empty()) {
        bound
            .parse::<i64>()
            .with_context(|| format!("invalid literal for int() with base 10: '{bound}'"))?;
    }
    Ok(())
}

#[cfg(test)]
pub(super) fn calculate_fp_region_size(
    requested: Option<&str>,
    fp_regions: &[Interval],
    ambiguous_regions: &[AmbiguousInterval],
    locations: Option<&[vcf::LocationFilter]>,
    reference_sequences: &BTreeMap<String, String>,
    truth: &[FilteredRawRecord],
) -> usize {
    calculate_fp_region_size_for_contigs(
        requested,
        fp_regions,
        ambiguous_regions,
        locations,
        reference_sequences,
        &contigs_in_truth(truth),
    )
}

pub(super) fn calculate_fp_region_size_for_contigs(
    requested: Option<&str>,
    fp_regions: &[Interval],
    ambiguous_regions: &[AmbiguousInterval],
    locations: Option<&[vcf::LocationFilter]>,
    reference_sequences: &BTreeMap<String, String>,
    truth_contigs: &BTreeSet<String>,
) -> usize {
    if let Some(size) = requested.and_then(|value| value.parse::<usize>().ok()) {
        return size;
    }

    if has_automatic_fp_bases(fp_regions, ambiguous_regions) {
        // BedIntervalTree stores each value as a list (`[label, ...]`), but
        // legacy som.py's location branch compares that list directly with
        // the string "FP". The comparison never succeeds, so any `-l` used
        // with labeled FP bases produces the historical zero denominator.
        if locations.is_some() {
            return 0;
        }
        return automatic_fp_intervals(fp_regions, ambiguous_regions)
            .map(|interval| interval.end.saturating_sub(interval.start))
            .sum();
    }

    if let Some(locations) = locations {
        return locations
            .iter()
            .map(|location| match location {
                vcf::LocationFilter::Contig(chrom) => reference_sequences
                    .get(chrom)
                    .map_or(0, |sequence| sequence.len()),
                vcf::LocationFilter::Range { chrom, start, end } => {
                    reference_sequences.get(chrom).map_or(0, |sequence| {
                        end.min(&sequence.len())
                            .saturating_sub(start.saturating_sub(1))
                    })
                }
            })
            .sum();
    }

    truth_contigs
        .iter()
        .filter_map(|contig| {
            reference_sequences
                .get(contig)
                .map(|sequence| sequence.len())
        })
        .sum()
}

pub(super) fn pair_exact_records(
    truth: &[FilteredRawRecord],
    query: &[FilteredRawRecord],
) -> (Vec<Option<usize>>, Vec<Option<usize>>) {
    let mut query_by_key: BTreeMap<vcf::VariantKey, VecDeque<usize>> = BTreeMap::new();
    for (index, record) in query.iter().enumerate() {
        query_by_key
            .entry(record.key.clone())
            .or_default()
            .push_back(index);
    }

    let mut truth_matches = vec![None; truth.len()];
    let mut query_matches = vec![None; query.len()];
    for (truth_index, truth_record) in truth.iter().enumerate() {
        let Some(query_index) = query_by_key
            .get_mut(&truth_record.key)
            .and_then(VecDeque::pop_front)
        else {
            continue;
        };
        truth_matches[truth_index] = Some(query_index);
        query_matches[query_index] = Some(truth_index);
    }
    (truth_matches, query_matches)
}

pub(super) fn render_row(
    index: usize,
    label: &str,
    counts: SomaticCounts,
    context: &StatsRowContext<'_>,
) -> String {
    render_row_standard_order(index, label, counts, context, false)
}

pub(super) fn render_row_af_without_bins(
    index: usize,
    label: &str,
    counts: SomaticCounts,
    context: &StatsRowContext<'_>,
) -> String {
    render_row_standard_order(index, label, counts, context, true)
}

pub(super) fn render_row_standard_order(
    index: usize,
    label: &str,
    counts: SomaticCounts,
    context: &StatsRowContext<'_>,
    blank_undefined_ratios: bool,
) -> String {
    let (recall, recall_lower, recall_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fn_count, context.ci_alpha);
    let (precision, precision_lower, precision_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fp, context.ci_alpha);
    let version = SOM_VERSION;
    let mut columns = vec![
        index.to_string(),
        label.to_string(),
        counts.truth_total.to_string(),
        counts.query_total.to_string(),
        counts.tp.to_string(),
        counts.fp.to_string(),
        counts.fn_count.to_string(),
        counts.unk.to_string(),
        counts.ambi.to_string(),
    ];
    if context.include_filtered_columns {
        if let Some(filtered) = context.filtered {
            columns.extend([
                py_float(filtered.fp as f64),
                py_float(filtered.tp as f64),
                py_float(filtered.unk as f64),
                py_float(filtered.ambi as f64),
            ]);
        } else {
            columns.extend([String::new(), String::new(), String::new(), String::new()]);
        }
    }
    columns.extend([
        py_float(recall),
        py_float(recall_lower),
        py_float(recall_upper),
        formatted_ratio(counts.tp, counts.truth_total, blank_undefined_ratios),
        py_float(precision),
        py_float(precision_lower),
        py_float(precision_upper),
        formatted_ratio(counts.unk, counts.query_total, blank_undefined_ratios),
        formatted_ratio(counts.ambi, counts.query_total, blank_undefined_ratios),
        context.fp_region_size.to_string(),
        fp_rate_or_blank(counts.fp, context.fp_region_size),
    ]);
    if context.include_filtered_columns {
        if let Some(filtered) = context.filtered {
            let unfiltered_tp = counts.tp.saturating_sub(filtered.tp);
            let unfiltered_fp = counts.fp.saturating_sub(filtered.fp);
            columns.extend([
                formatted_ratio(
                    unfiltered_tp,
                    counts.tp + counts.fn_count,
                    blank_undefined_ratios,
                ),
                formatted_ratio(
                    unfiltered_tp,
                    unfiltered_tp + unfiltered_fp,
                    blank_undefined_ratios,
                ),
                fp_rate_or_blank(unfiltered_fp, context.fp_region_size),
                formatted_ratio(
                    counts.unk.saturating_sub(filtered.unk),
                    counts.query_total,
                    blank_undefined_ratios,
                ),
                formatted_ratio(
                    counts.ambi.saturating_sub(filtered.ambi),
                    counts.query_total,
                    blank_undefined_ratios,
                ),
            ]);
        } else {
            columns.extend([
                String::new(),
                String::new(),
                String::new(),
                String::new(),
                String::new(),
            ]);
        }
    }
    columns.push(version.to_string());
    columns.push(context.commandline.to_string());
    columns.join(",")
}

pub(super) fn formatted_ratio(
    numerator: usize,
    denominator: usize,
    blank_undefined: bool,
) -> String {
    if blank_undefined {
        ratio_or_blank(numerator, denominator)
    } else {
        py_float(ratio(numerator, denominator))
    }
}

/// Render the column order produced by the legacy pandas concat used when
/// allele-frequency stratification is enabled.  Adding the bin rows causes
/// pandas 0.x to alphabetize the count columns while leaving subsequently
/// calculated metrics in assignment order; that accidental ordering is part
/// of the byte-for-byte CSV contract.
pub(super) fn render_row_af(
    index: usize,
    label: &str,
    counts: SomaticCounts,
    context: &StatsRowContext<'_>,
) -> String {
    let (recall, recall_lower, recall_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fn_count, context.ci_alpha);
    let (precision, precision_lower, precision_upper) =
        jeffreys_ci(counts.tp, counts.tp + counts.fp, context.ci_alpha);
    let mut columns = vec![index.to_string(), counts.ambi.to_string()];
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.ambi as f64)),
        );
    }
    columns.extend([counts.fn_count.to_string(), counts.fp.to_string()]);
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.fp as f64)),
        );
    }
    columns.extend([
        counts.query_total.to_string(),
        counts.truth_total.to_string(),
        counts.tp.to_string(),
    ]);
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.tp as f64)),
        );
    }
    columns.extend([label.to_string(), counts.unk.to_string()]);
    if context.include_filtered_columns {
        columns.push(
            context
                .filtered
                .map_or_else(String::new, |counts| py_float(counts.unk as f64)),
        );
    }
    columns.extend([
        py_float(recall),
        py_float(recall_lower),
        py_float(recall_upper),
        ratio_or_blank(counts.tp, counts.truth_total),
        py_float(precision),
        py_float(precision_lower),
        py_float(precision_upper),
        ratio_or_blank(counts.unk, counts.query_total),
        ratio_or_blank(counts.ambi, counts.query_total),
        context.fp_region_size.to_string(),
        fp_rate_or_blank(counts.fp, context.fp_region_size),
    ]);
    if context.include_filtered_columns {
        if let Some(filtered) = context.filtered {
            let unfiltered_tp = counts.tp.saturating_sub(filtered.tp);
            let unfiltered_fp = counts.fp.saturating_sub(filtered.fp);
            columns.extend([
                ratio_or_blank(unfiltered_tp, counts.tp + counts.fn_count),
                ratio_or_blank(unfiltered_tp, unfiltered_tp + unfiltered_fp),
                fp_rate_or_blank(unfiltered_fp, context.fp_region_size),
                ratio_or_blank(counts.unk.saturating_sub(filtered.unk), counts.query_total),
                ratio_or_blank(
                    counts.ambi.saturating_sub(filtered.ambi),
                    counts.query_total,
                ),
            ]);
        } else {
            columns.extend(std::iter::repeat_n(String::new(), 5));
        }
    }
    columns.push(SOM_VERSION.to_string());
    columns.push(context.commandline.to_string());
    columns.join(",")
}

pub(super) fn ratio_or_blank(numerator: usize, denominator: usize) -> String {
    if denominator == 0 {
        String::new()
    } else {
        py_float(numerator as f64 / denominator as f64)
    }
}

pub(super) fn fp_rate_or_blank(false_positives: usize, region_size: usize) -> String {
    if region_size == 0 {
        if false_positives == 0 {
            String::new()
        } else {
            "inf".to_string()
        }
    } else {
        py_float(1_000_000.0 * false_positives as f64 / region_size as f64)
    }
}

pub(super) fn stats_header(count_filtered_fn: bool) -> String {
    let mut columns = vec![
        "",
        "type",
        "total.truth",
        "total.query",
        "tp",
        "fp",
        "fn",
        "unk",
        "ambi",
    ];
    if count_filtered_fn {
        columns.extend([
            "fp.filtered",
            "tp.filtered",
            "unk.filtered",
            "ambi.filtered",
        ]);
    }
    columns.extend([
        "recall",
        "recall_lower",
        "recall_upper",
        "recall2",
        "precision",
        "precision_lower",
        "precision_upper",
        "na",
        "ambiguous",
        "fp.region.size",
        "fp.rate",
    ]);
    if count_filtered_fn {
        columns.extend([
            "recall.filtered",
            "precision.filtered",
            "fp.rate.filtered",
            "na.filtered",
            "ambiguous.filtered",
        ]);
    }
    columns.extend(["sompyversion", "sompycmd"]);
    columns.join(",")
}

pub(super) fn stats_header_af(count_filtered_fn: bool) -> String {
    let mut columns = vec![""];
    columns.push("ambi");
    if count_filtered_fn {
        columns.push("ambi.filtered");
    }
    columns.extend(["fn", "fp"]);
    if count_filtered_fn {
        columns.push("fp.filtered");
    }
    columns.extend(["total.query", "total.truth", "tp"]);
    if count_filtered_fn {
        columns.push("tp.filtered");
    }
    columns.extend(["type", "unk"]);
    if count_filtered_fn {
        columns.push("unk.filtered");
    }
    columns.extend([
        "recall",
        "recall_lower",
        "recall_upper",
        "recall2",
        "precision",
        "precision_lower",
        "precision_upper",
        "na",
        "ambiguous",
        "fp.region.size",
        "fp.rate",
    ]);
    if count_filtered_fn {
        columns.extend([
            "recall.filtered",
            "precision.filtered",
            "fp.rate.filtered",
            "na.filtered",
            "ambiguous.filtered",
        ]);
    }
    columns.extend(["sompyversion", "sompycmd"]);
    columns.join(",")
}

pub(super) fn filtered_counts_for_type(
    enabled: bool,
    feature_table: Option<&str>,
    label: &str,
    filtered_by_type: &BTreeMap<&'static str, FilteredCounts>,
) -> Option<FilteredCounts> {
    if !enabled {
        return None;
    }
    if feature_table == Some("generic") {
        return Some(filtered_by_type.get(label).copied().unwrap_or_default());
    }
    let selected_type = match feature_table.and_then(|name| name.rsplit('.').next()) {
        Some("snv") => "SNVs",
        Some("indel") => "indels",
        _ => return None,
    };
    if label == selected_type {
        Some(
            filtered_by_type
                .get(selected_type)
                .copied()
                .unwrap_or_default(),
        )
    } else {
        None
    }
}
