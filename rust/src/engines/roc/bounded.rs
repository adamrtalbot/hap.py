//! ROC enumeration for `hap germline`.
//!
//! Emits `roc.all` plus the non-empty per-type ROC CSVs per prefix:
//!   result.roc.all.csv.gz, result.roc.Locations.SNP.csv.gz,
//!   result.roc.Locations.SNP.PASS.csv.gz, result.roc.Locations.INDEL.csv.gz,
//!   result.roc.Locations.INDEL.PASS.csv.gz.
//!
//! Each file mirrors the 65-column schema of `result.extended.csv` but adds
//! per-QQ cumulative rows: for every grouping tuple `(Type, Subtype, Subset,
//! Filter, Genotype="*", QQ.Field="QUAL")` we emit a `QQ="*"` baseline row
//! whose counts sum every contribution in the group, followed by one row per
//! distinct numeric QQ threshold carrying the cumulative counts at `QQ >=
//! threshold`. Row order within a group is lexicographic ascending over the
//! QQ column string (`"*"` sorts before the six-decimal-formatted numeric
//! strings, exactly matching pandas `sort_values` on an object column).

use crate::adapters::report::{
    EXTENDED_HEADER, append_stats_with_missing, empty_comparison_extended_lines, f1_score,
    format_count, full_repr_float, het_hom_ratio, metric_ratio, suffixed_report_path, ti_tv_ratio,
};
use crate::domain::{AnnotatedRow, CountsBucket, FpClass};
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::borrow::Borrow;
use std::cmp::{Ordering, Reverse};
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;
use std::sync::Arc;

const INDEL_SUBTYPES: [&str; 9] = [
    "C16_PLUS", "C1_5", "C6_15", "D16_PLUS", "D1_5", "D6_15", "I16_PLUS", "I1_5", "I6_15",
];
const ROC_OBSERVATION_CHUNK: usize = 16_384;
const ROC_INDEX_ENTRY_BYTES: u64 = 16;
// The full GIAB profile produces roughly 62 million observations: a 1 GiB
// fixed-width index. Keeping that index in memory avoids billions of tiny
// random writes to the temporary filesystem while leaving the much larger
// variable-width observation payload disk-backed. The production pipeline
// gives hap.py 24 GiB, so a 2 GiB bound remains conservative.
const ROC_IN_MEMORY_INDEX_LIMIT_BYTES: u64 = 2 * 1024 * 1024 * 1024;
const MAX_RENDERED_ROC_THRESHOLDS: usize = 2_500_000;
const MAX_ROC_METRIC_INDEX_KEYS: usize = 500_000;

/// Write `roc.all` and each non-empty Locations ROC file alongside `prefix`.
#[derive(Clone, Debug, Default)]
pub(crate) struct MetricIndices {
    pub(crate) tables: BTreeMap<String, Vec<usize>>,
    /// Table order produced by Python 2.7's insertion-ordered hash table
    /// iteration in `happyroc.roc`. This is data-dependent because location
    /// tables are inserted when their first raw ROC row is encountered.
    pub(crate) table_order: Vec<String>,
}

/// Controls inherited from qfy's ROC command line.  The default deliberately
/// remains identical to the historical hap.py invocation.
#[derive(Clone, Debug)]
pub(crate) struct RocOptions {
    pub(crate) threads: usize,
    pub(crate) qq_field: String,
    /// Optional FORMAT/INFO field used for thresholds while `qq_field`
    /// remains the user-facing label in metrics and ROC tables.
    pub(crate) score_field: Option<String>,
    pub(crate) ignored_filters: HashSet<String>,
    pub(crate) roc_regions: HashSet<String>,
    pub(crate) delta: f64,
    pub(crate) ci_alpha: f64,
    /// Preserve qfy's private C++ quantifier table. Legacy qfy removes this
    /// intermediate unless `--verbose` is active.
    pub(crate) preserve_raw_table: bool,
    /// Include threshold rows in the private table. This follows qfy's
    /// `--roc`/`--no-roc` switch independently of the public compacting pass.
    pub(crate) output_rocs: bool,
    /// Full N-trimmed FASTA size used by the legacy TS_boundary lane.
    /// `None` preserves the historical caller contract where `subset_size`
    /// is also the complete reference size.
    pub(crate) whole_reference_size: Option<usize>,
    /// Per-user-named stratification interval-union sizes.
    pub(crate) subset_sizes: BTreeMap<String, usize>,
    /// Per-named-subset intersection with the confidence regions.
    pub(crate) subset_confidence_sizes: BTreeMap<String, usize>,
}

impl Default for RocOptions {
    fn default() -> Self {
        Self {
            threads: 1,
            qq_field: "QUAL".to_string(),
            score_field: None,
            ignored_filters: HashSet::new(),
            roc_regions: HashSet::from(["*".to_string()]),
            delta: 0.5,
            ci_alpha: 0.0,
            preserve_raw_table: false,
            output_rocs: true,
            whole_reference_size: None,
            subset_sizes: BTreeMap::new(),
            subset_confidence_sizes: BTreeMap::new(),
        }
    }
}

/// Write ROC artifacts and return the original pandas row indices used by
/// legacy qfy's metrics JSON tables.
#[cfg(test)]
fn write_roc_files(
    prefix: &Path,
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
) -> Result<MetricIndices> {
    write_roc_files_with_options(prefix, rows, subset_size, conf_size, &RocOptions::default())
}

/// Write ROC artifacts with the qfy controls supplied by the caller.
#[cfg(test)]
fn write_roc_files_with_options(
    prefix: &Path,
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<MetricIndices> {
    write_roc_files_with_options_iter(
        prefix,
        rows.iter().map(Ok::<_, anyhow::Error>),
        subset_size,
        conf_size,
        options,
    )
}

pub(crate) fn write_roc_files_with_options_iter<I, R>(
    prefix: &Path,
    rows: I,
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<MetricIndices>
where
    I: IntoIterator<Item = Result<R>>,
    R: Borrow<AnnotatedRow>,
{
    if options.qq_field.is_empty() {
        bail!("ROC field cannot be empty");
    }
    if !options.delta.is_finite() || options.delta < 0.0 {
        bail!("ROC delta must be finite and nonnegative");
    }
    if !options.ci_alpha.is_finite()
        || (options.ci_alpha != 0.0 && !(0.0 < options.ci_alpha && options.ci_alpha < 1.0))
    {
        bail!("ROC CI alpha must be 0 (disabled) or strictly between 0 and 1");
    }
    let groups = accumulate_impl(rows, options)?;

    if !groups.values().any(|group| group.records > 0) {
        let all = RenderedRows::from_lines(empty_comparison_extended_lines(subset_size))?;
        write_gzip_csv(
            &suffixed_report_path(prefix, "roc.all.csv.gz"),
            &roc_header(options.ci_alpha),
            &all,
        )?;
        let indices = vec![1, 0];
        return Ok(MetricIndices {
            tables: BTreeMap::from([
                ("summary.metrics".to_string(), indices.clone()),
                ("all.metrics".to_string(), indices.clone()),
                ("roc.all".to_string(), indices),
            ]),
            table_order: vec!["roc.all".to_string()],
        });
    }

    if options.preserve_raw_table {
        write_legacy_roc_table(prefix, &groups, subset_size, conf_size, options)?;
        clear_all_sorted_indices(&groups);
    }

    // Compute per-subtype star_sorted snapshots ONCE (heavy operation: up to
    // 10×4 sorts × full obs vector clone for INDEL groups). Pass by reference
    // to render_rows so the cost isn't paid 5 times.
    let star_sorted = build_star_sorted(&groups);

    let header = roc_header(options.ci_alpha);
    let render_config = RenderConfig {
        subset_size,
        whole_reference_size: options.whole_reference_size.unwrap_or(subset_size),
        conf_size,
        subset_sizes: &options.subset_sizes,
        subset_confidence_sizes: &options.subset_confidence_sizes,
        delta: options.delta,
        ci_alpha: options.ci_alpha,
        filter_counts_only: options.roc_regions.contains("*"),
    };
    // `roc.all`: every row from every group, no filtering.
    let all = render_rows_parallel(&groups, &star_sorted, render_config, options.threads)?;
    write_gzip_csv(
        &suffixed_report_path(prefix, "roc.all.csv.gz"),
        &header,
        &all,
    )?;

    // `roc.Locations.<TYPE>[.PASS]`: matches legacy `happyroc.py` Locations
    // filter — keeps only `(Type, Subtype=*, Subset=*, Genotype=*, Filter,
    // QQ != *)` rows. Each file is the per-QQ cumulative threshold sweep
    // for the corresponding Type under a single Filter.
    // Locations tables are an exact row subset of roc.all. Re-filter the
    // already-rendered ordered stream instead of re-running the same four
    // libstdc++-compatible sorts for each public table.
    // The public Locations tables are all projections of the same ordered
    // `roc.all` stream.  Classify a line once while it is warm in the page
    // cache instead of rescanning the (potentially multi-gigabyte) rendered
    // spool once per table.  This preserves both the source order and the
    // bounded on-disk representation used by metric-index publication.
    let locations = filter_locations_rows(&all, !options.ignored_filters.is_empty())?;
    let snp = locations.snp;
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.SNP.csv.gz"),
        &header,
        &snp,
    )?;
    let snp_pass = locations.snp_pass;
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.SNP.PASS.csv.gz"),
        &header,
        &snp_pass,
    )?;
    let indel = locations.indel;
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.INDEL.csv.gz"),
        &header,
        &indel,
    )?;
    let indel_pass = locations.indel_pass;
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.INDEL.PASS.csv.gz"),
        &header,
        &indel_pass,
    )?;

    let mut selective = Vec::new();
    if let Some(snp_sel) = locations.snp_sel {
        write_optional_gzip_csv(
            &suffixed_report_path(prefix, "roc.Locations.SNP.SEL.csv.gz"),
            &header,
            &snp_sel,
        )?;
        selective.push(("roc.Locations.SNP.SEL".to_string(), snp_sel, "SNP"));
        let indel_sel = locations.indel_sel.expect("SEL location spools are paired");
        write_optional_gzip_csv(
            &suffixed_report_path(prefix, "roc.Locations.INDEL.SEL.csv.gz"),
            &header,
            &indel_sel,
        )?;
        selective.push(("roc.Locations.INDEL.SEL".to_string(), indel_sel, "INDEL"));
    }

    build_metric_indices(
        &groups,
        MetricRows {
            all: &all,
            snp: &snp,
            snp_pass: &snp_pass,
            indel: &indel,
            indel_pass: &indel_pass,
            selective: &selective,
        },
        options.delta,
        options.output_rocs,
        options.threads,
    )
}

struct RenderedRows {
    path: tempfile::TempPath,
    len: usize,
}

struct LocationRows {
    snp: RenderedRows,
    snp_pass: RenderedRows,
    indel: RenderedRows,
    indel_pass: RenderedRows,
    snp_sel: Option<RenderedRows>,
    indel_sel: Option<RenderedRows>,
}

fn filter_locations_rows(all: &RenderedRows, include_selective: bool) -> Result<LocationRows> {
    let mut outputs = (0..if include_selective { 6 } else { 4 })
        .map(|_| {
            tempfile::NamedTempFile::new().context("failed to create Locations ROC report spool")
        })
        .collect::<Result<Vec<_>>>()?;
    let mut lengths = vec![0usize; outputs.len()];
    let mut writers = outputs
        .iter_mut()
        .map(|output| BufWriter::new(output.as_file_mut()))
        .collect::<Vec<_>>();
    for line in all.lines()? {
        let line = line?;
        let mut fields = line.splitn(8, ',');
        let ty = fields.next();
        let subtype = fields.next();
        let subset = fields.next();
        let filter = fields.next();
        let genotype = fields.next();
        let _qq_field = fields.next();
        let qq = fields.next();
        if subtype != Some("*")
            || subset != Some("*")
            || genotype != Some("*")
            || qq == Some("*")
            || qq.is_none()
        {
            continue;
        }
        let index = match (ty, filter) {
            (Some("SNP"), Some("ALL")) => Some(0),
            (Some("SNP"), Some("PASS")) => Some(1),
            (Some("INDEL"), Some("ALL")) => Some(2),
            (Some("INDEL"), Some("PASS")) => Some(3),
            (Some("SNP"), Some("SEL")) if include_selective => Some(4),
            (Some("INDEL"), Some("SEL")) if include_selective => Some(5),
            _ => None,
        };
        if let Some(index) = index {
            writeln!(writers[index], "{line}")?;
            lengths[index] += 1;
        }
    }
    for writer in &mut writers {
        writer.flush()?;
    }
    drop(writers);
    let mut rows = outputs
        .into_iter()
        .zip(lengths)
        .map(|(output, len)| RenderedRows {
            path: output.into_temp_path(),
            len,
        });
    Ok(LocationRows {
        snp: rows.next().expect("SNP spool exists"),
        snp_pass: rows.next().expect("SNP PASS spool exists"),
        indel: rows.next().expect("INDEL spool exists"),
        indel_pass: rows.next().expect("INDEL PASS spool exists"),
        snp_sel: include_selective.then(|| rows.next().expect("SNP SEL spool exists")),
        indel_sel: include_selective.then(|| rows.next().expect("INDEL SEL spool exists")),
    })
}

impl RenderedRows {
    fn from_lines(lines: impl IntoIterator<Item = String>) -> Result<Self> {
        let mut file =
            tempfile::NamedTempFile::new().context("failed to create rendered ROC report spool")?;
        let mut len = 0usize;
        for line in lines {
            writeln!(file.as_file_mut(), "{line}")?;
            len += 1;
        }
        file.as_file_mut().flush()?;
        Ok(Self {
            path: file.into_temp_path(),
            len,
        })
    }

    fn lines(&self) -> Result<std::io::Lines<BufReader<File>>> {
        Ok(BufReader::new(File::open(&self.path)?).lines())
    }
}

struct MetricRows<'a> {
    all: &'a RenderedRows,
    snp: &'a RenderedRows,
    snp_pass: &'a RenderedRows,
    indel: &'a RenderedRows,
    indel_pass: &'a RenderedRows,
    selective: &'a [(String, RenderedRows, &'a str)],
}

type MaskedMetricLevels = BTreeMap<(String, String), Vec<f64>>;

fn metric_subtype_flags(ty: &str) -> &'static [(&'static str, Option<&'static str>)] {
    if ty == "SNP" {
        &[("*", None), ("ti", Some("ti")), ("tv", Some("tv"))]
    } else {
        &[
            ("*", None),
            ("I1_5", Some("I1_5")),
            ("I6_15", Some("I6_15")),
            ("I16_PLUS", Some("I16_PLUS")),
            ("D1_5", Some("D1_5")),
            ("D6_15", Some("D6_15")),
            ("D16_PLUS", Some("D16_PLUS")),
            ("C1_5", Some("C1_5")),
            ("C6_15", Some("C6_15")),
            ("C16_PLUS", Some("C16_PLUS")),
        ]
    }
}

fn build_metric_level_sets(
    groups: &[(&RowKey, &BoundedGroupAccum)],
    delta: f64,
    output_rocs: bool,
    threads: usize,
) -> Result<Vec<MaskedMetricLevels>> {
    if groups.len() <= 1 || threads <= 1 {
        return groups
            .iter()
            .map(|(key, accum)| {
                if !output_rocs || !is_aggregate_filter(&key.filter) {
                    Ok(BTreeMap::new())
                } else {
                    legacy_masked_level_sets(
                        &accum.observations,
                        metric_subtype_flags(&key.ty),
                        delta,
                    )
                }
            })
            .collect();
    }

    // Two simultaneous compact level vectors keep peak RSS close to the
    // existing bounded profile while overlapping the independent ALL/PASS
    // and SNP/INDEL metric walks.
    let worker_count = threads.max(1).min(2).min(groups.len());
    let next_group = std::sync::atomic::AtomicUsize::new(0);
    let (sender, receiver) = std::sync::mpsc::sync_channel(worker_count);
    let mut results = (0..groups.len()).map(|_| None).collect::<Vec<_>>();
    let mut first_error = None;
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let sender = sender.clone();
            let next_group = &next_group;
            handles.push(scope.spawn(move || {
                loop {
                    let index = next_group.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some((key, accum)) = groups.get(index).copied() else {
                        break;
                    };
                    let result = if !output_rocs || !is_aggregate_filter(&key.filter) {
                        Ok(BTreeMap::new())
                    } else {
                        legacy_masked_level_sets(
                            &accum.observations,
                            metric_subtype_flags(&key.ty),
                            delta,
                        )
                    };
                    if sender.send((index, result)).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(sender);
        for (index, result) in receiver {
            match result {
                Ok(levels) if first_error.is_none() => results[index] = Some(levels),
                Ok(_) => {}
                Err(error) if first_error.is_none() => first_error = Some(error),
                Err(_) => {}
            }
        }
        for handle in handles {
            if handle.join().is_err() && first_error.is_none() {
                first_error = Some(anyhow::anyhow!("ROC metric worker panicked"));
            }
        }
    });
    if let Some(error) = first_error {
        return Err(error);
    }
    Ok(results
        .into_iter()
        .map(|result| result.expect("each ROC metric group returned a result"))
        .collect())
}

fn build_metric_indices(
    groups: &BTreeMap<RowKey, BoundedGroupAccum>,
    rows: MetricRows<'_>,
    delta: f64,
    output_rocs: bool,
    threads: usize,
) -> Result<MetricIndices> {
    let mut rocs = BTreeMap::<String, (&RowKey, &BoundedGroupAccum)>::new();
    let active_types = groups
        .iter()
        .filter(|(_, accum)| accum.records > 0)
        .map(|(key, _)| key.ty.as_str())
        .collect::<HashSet<_>>();
    for (key, accum) in groups {
        if key.subtype != "*" || !active_types.contains(key.ty.as_str()) {
            continue;
        }
        let name = if key.subset != "*" {
            format!("s|{}:{}:{}", key.subset, key.ty, key.filter)
        } else if is_aggregate_filter(&key.filter) {
            format!("a:{}:{}", key.ty, key.filter)
        } else {
            format!("f:{}:{}", key.ty, key.filter)
        };
        rocs.insert(name, (key, accum));
    }

    let roc_groups = rocs.into_values().collect::<Vec<_>>();
    let masked_level_sets = build_metric_level_sets(&roc_groups, delta, output_rocs, threads)?;
    let mut table = LegacyUnorderedRows::default();
    for ((key, _accum), masked_levels) in roc_groups.into_iter().zip(masked_level_sets) {
        let subtype_flags = metric_subtype_flags(&key.ty);
        let counts_only = !output_rocs || !is_aggregate_filter(&key.filter);
        for (subtype, _) in subtype_flags {
            for (genotype, _) in [
                ("het", Some("het")),
                ("hetalt", Some("hetalt")),
                ("homalt", Some("homalt")),
                ("*", None),
            ] {
                let baseline = legacy_row_key(&key.ty, subtype, &key.filter, &key.subset, "*");
                if !matches!(*subtype, "ti" | "tv") {
                    table.set(baseline, genotype == "*")?;
                } else if genotype == "*" {
                    table.set(
                        legacy_row_key(&key.ty, "*", &key.filter, &key.subset, "*"),
                        false,
                    )?;
                }
                if counts_only {
                    continue;
                }
                let mask_key = ((*subtype).to_string(), genotype.to_string());
                for level in masked_levels.get(&mask_key).into_iter().flatten().copied() {
                    let qq = format!("{level:.6}");
                    if !matches!(*subtype, "ti" | "tv") && genotype == "*" {
                        table.set(
                            legacy_row_key(&key.ty, subtype, &key.filter, &key.subset, &qq),
                            true,
                        )?;
                    } else if matches!(*subtype, "ti" | "tv") && genotype == "*" {
                        table.set(
                            legacy_row_key(&key.ty, "*", &key.filter, &key.subset, &qq),
                            false,
                        )?;
                    } else if !matches!(*subtype, "ti" | "tv") {
                        // Preserve the extra empty field in the legacy
                        // genotype aggregation prefix. These rows are removed
                        // by dropRowsWithMissing(Type), but their insertion can
                        // trigger an unordered_map rehash and therefore affects
                        // the pandas indices of retained rows.
                        table.set(
                            format!(
                                "{}\t{}\t*\t\t{}\t{}\t{}",
                                key.ty, subtype, key.filter, key.subset, qq
                            ),
                            false,
                        )?;
                    }
                }
            }
        }
    }

    let raw = table.retained_order();
    let raw_positions = raw
        .iter()
        .enumerate()
        .map(|(index, key)| (key.clone(), index))
        .collect::<BTreeMap<_, _>>();
    // qfy's `--no-roc` prevents C++ from placing threshold rows in the raw
    // table. The public CSV is compacted to the same baseline-only set before
    // JSON serialization, so its retained pandas indices must be computed
    // against that compact set rather than falling back to 0..N.
    let all_indices = indices_for_lines(metric_lines(rows.all, !output_rocs)?, &raw_positions)?;

    let mut available_tables = BTreeSet::from(["roc.all"]);
    for (id, lines) in [
        ("roc.Locations.SNP", rows.snp),
        ("roc.Locations.SNP.PASS", rows.snp_pass),
        ("roc.Locations.INDEL", rows.indel),
        ("roc.Locations.INDEL.PASS", rows.indel_pass),
    ] {
        if lines.len != 0 {
            available_tables.insert(id);
        }
    }
    for (id, lines, _) in rows.selective {
        if lines.len != 0 {
            available_tables.insert(id.as_str());
        }
    }
    let table_order = legacy_python_table_order(&raw, &available_tables);

    let mut tables = BTreeMap::new();
    tables.insert("roc.all".to_string(), all_indices.clone());
    let mut all_metrics = Vec::new();
    let mut summary_metrics = Vec::new();
    for (line, index) in metric_lines(rows.all, !output_rocs)?.zip(&all_indices) {
        let line = line?;
        let fields = line.split(',').collect::<Vec<_>>();
        if fields.get(6) == Some(&"*") && matches!(fields.get(3), Some(&"ALL") | Some(&"PASS")) {
            all_metrics.push(*index);
        }
        if fields.get(1) == Some(&"*")
            && fields.get(2) == Some(&"*")
            && matches!(fields.get(3), Some(&"ALL") | Some(&"PASS"))
            && fields.get(4) == Some(&"*")
            && fields.get(6) == Some(&"*")
        {
            summary_metrics.push(*index);
        }
    }
    tables.insert("all.metrics".to_string(), all_metrics);
    tables.insert("summary.metrics".to_string(), summary_metrics);
    for (id, lines, ty, filter) in [
        ("roc.Locations.SNP", rows.snp, "SNP", "ALL"),
        ("roc.Locations.SNP.PASS", rows.snp_pass, "SNP", "PASS"),
        ("roc.Locations.INDEL", rows.indel, "INDEL", "ALL"),
        ("roc.Locations.INDEL.PASS", rows.indel_pass, "INDEL", "PASS"),
    ] {
        let mut local = BTreeMap::new();
        let mut next = 0usize;
        for raw_key in &raw {
            let fields = raw_key.split('\t').collect::<Vec<_>>();
            if fields.len() == 6
                && fields[0] == ty
                && fields[1] == "*"
                && fields[2] == "*"
                && fields[3] == filter
                && fields[4] == "*"
                && fields[5] != "*"
            {
                local.insert(raw_key.clone(), next);
                next += 1;
            }
        }
        tables.insert(
            id.to_string(),
            indices_for_lines(lines.lines()?.map(|line| line.map_err(Into::into)), &local)?,
        );
    }
    for (id, lines, ty) in rows.selective {
        let mut local = BTreeMap::new();
        let mut next = 0usize;
        for raw_key in &raw {
            let fields = raw_key.split('\t').collect::<Vec<_>>();
            if fields.len() == 6
                && fields[0] == *ty
                && fields[1] == "*"
                && fields[2] == "*"
                && fields[3] == "SEL"
                && fields[4] == "*"
                && fields[5] != "*"
            {
                local.insert(raw_key.clone(), next);
                next += 1;
            }
        }
        tables.insert(
            id.clone(),
            indices_for_lines(lines.lines()?.map(|line| line.map_err(Into::into)), &local)?,
        );
    }
    Ok(MetricIndices {
        tables,
        table_order,
    })
}

fn metric_lines(
    rows: &RenderedRows,
    baseline_only: bool,
) -> Result<Box<dyn Iterator<Item = Result<String>>>> {
    Ok(Box::new(rows.lines()?.filter_map(move |line| match line {
        Ok(line) if !baseline_only || line.split(',').nth(6) == Some("*") => Some(Ok(line)),
        Ok(_) => None,
        Err(error) => Some(Err(error.into())),
    })))
}

/// Reproduce the iteration order of the pinned Python 2.7 dictionary used by
/// `qfy.py` to append `res` tables to metrics JSON. The dictionary keys and
/// hashes are fixed by the reference (`hash_randomization=0`); insertion order is
/// determined by the first matching row in the C++ unordered ROC table.
fn legacy_python_table_order(raw: &[String], available_tables: &BTreeSet<&str>) -> Vec<String> {
    let mut insertions = Vec::new();
    for raw_key in raw {
        let fields = raw_key.split('\t').collect::<Vec<_>>();
        if fields.len() == 6
            && matches!(fields[0], "SNP" | "INDEL")
            && fields[1] == "*"
            && fields[2] == "*"
            && fields[4] == "*"
            && fields[5] != "*"
        {
            let suffix = match fields[3] {
                "ALL" => "",
                "PASS" => ".PASS",
                "SEL" => ".SEL",
                _ => "",
            };
            if matches!(fields[3], "ALL" | "PASS" | "SEL") {
                let table = format!("roc.Locations.{}{suffix}", fields[0]);
                if available_tables.contains(table.as_str()) {
                    insertions.push(table);
                }
            }
        }
        insertions.push("roc.all".to_string());
    }

    python27_dict_iteration_order(&insertions)
}

fn python27_dict_iteration_order(insertions: &[String]) -> Vec<String> {
    #[derive(Clone)]
    struct Entry {
        key: String,
        hash: u64,
    }

    fn python_key(table: &str) -> &str {
        table.strip_prefix("roc.").unwrap_or(table)
    }

    fn hash(table: &str) -> u64 {
        match python_key(table) {
            "all" => 1_453_079_729_202_098_178,
            "Locations.SNP" => 14_809_960_256_853_955_918,
            "Locations.INDEL" => 8_577_287_029_120_378_309,
            "Locations.SNP.PASS" => 15_867_700_624_972_327_818,
            "Locations.INDEL.PASS" => 9_569_443_920_230_493_239,
            "Locations.SNP.SEL" => 2_407_900_041_572_956_280,
            "Locations.INDEL.SEL" => 15_566_398_430_702_022_267,
            key => panic!("unsupported legacy metrics table key: {key}"),
        }
    }

    fn insert(slots: &mut [Option<Entry>], entry: Entry) {
        let mask = slots.len() - 1;
        let mut index = entry.hash as usize & mask;
        let mut perturb = entry.hash;
        loop {
            if slots[index].is_none() {
                slots[index] = Some(entry);
                return;
            }
            index = index
                .wrapping_mul(5)
                .wrapping_add(perturb as usize)
                .wrapping_add(1)
                & mask;
            perturb >>= 5;
        }
    }

    let mut slots: Vec<Option<Entry>> = vec![None; 8];
    let mut used = 0usize;
    for key in insertions {
        if slots.iter().flatten().any(|entry| entry.key == *key) {
            continue;
        }
        insert(
            &mut slots,
            Entry {
                key: key.clone(),
                hash: hash(key),
            },
        );
        used += 1;

        // CPython 2.7 resizes once the table reaches two-thirds full and,
        // below 50k keys, asks for four times the number of live entries.
        if used * 3 >= slots.len() * 2 {
            let minimum = used * 4;
            let mut size = 8usize;
            while size <= minimum {
                size *= 2;
            }
            let old = std::mem::replace(&mut slots, vec![None; size]);
            for entry in old.into_iter().flatten() {
                insert(&mut slots, entry);
            }
        }
    }

    slots.into_iter().flatten().map(|entry| entry.key).collect()
}

fn legacy_row_key(ty: &str, subtype: &str, filter: &str, subset: &str, qq: &str) -> String {
    format!("{ty}\t{subtype}\t*\t{filter}\t{subset}\t{qq}")
}

fn csv_row_key(line: &str) -> String {
    let fields = line.split(',').take(7).collect::<Vec<_>>();
    if fields.len() < 7 {
        return String::new();
    }
    legacy_row_key(fields[0], fields[1], fields[3], fields[2], fields[6])
}

fn indices_for_lines(
    lines: impl Iterator<Item = Result<String>>,
    positions: &BTreeMap<String, usize>,
) -> Result<Vec<usize>> {
    lines
        .map(|line| {
            let line = line?;
            Ok(positions.get(&csv_row_key(&line)).copied().unwrap_or(0))
        })
        .collect()
}

/// Build every subtype/genotype level set for one legacy ROC group in a
/// single sequential pass. Metric-index construction only needs the ordered
/// threshold values; it does not consume cumulative counts, so repeatedly
/// introsorting and randomly rereading the full observation payload (up to
/// forty times for INDEL) cannot affect its result.
fn legacy_masked_level_sets(
    observations: &ObservationStore,
    subtype_flags: &[(&str, Option<&str>)],
    delta: f64,
) -> Result<BTreeMap<(String, String), Vec<f64>>> {
    let genotypes = [
        ("het", Some("het")),
        ("hetalt", Some("hetalt")),
        ("homalt", Some("homalt")),
        ("*", None),
    ];
    let mut kept = BTreeMap::<(String, String), Vec<f64>>::new();
    for (subtype, _) in subtype_flags {
        for (genotype, _) in genotypes {
            kept.insert(((*subtype).to_string(), genotype.to_string()), Vec::new());
        }
    }

    // Retain only the fields needed by the metric-index threshold walk and
    // sort them once. The former implementation built and sorted up to forty
    // multi-million-element f64 vectors for one INDEL group.
    let mut observations_by_level = Vec::<(f64, u16, u8)>::new();
    for observation in observations.unsorted()? {
        let observation = observation?;
        let mut subtype_mask = 0u16;
        for (index, (_, subtype_flag)) in subtype_flags.iter().enumerate() {
            let subtype_matches = match subtype_flag {
                None => true,
                Some("ti") => observation.ti_flag,
                Some("tv") => observation.tv_flag,
                Some(value) => observation.has_subtype(value),
            };
            if subtype_matches {
                subtype_mask |= 1u16 << index;
            }
        }
        observations_by_level.push((observation.level, subtype_mask, observation.blt));
    }
    observations_by_level
        .sort_by(|left, right| left.0.partial_cmp(&right.0).unwrap_or(Ordering::Equal));

    for (level, subtype_mask, observation_genotype) in observations_by_level {
        for (subtype_index, (subtype, _)) in subtype_flags.iter().enumerate() {
            if subtype_mask & (1u16 << subtype_index) == 0 {
                continue;
            }
            for (genotype_index, (genotype, _)) in genotypes.iter().enumerate() {
                if genotype_index != 3 && observation_genotype != genotype_index as u8 + 1 {
                    continue;
                }
                let levels = kept
                    .get_mut(&((**subtype).to_string(), (*genotype).to_string()))
                    .expect("legacy metric level mask was preseeded");
                if levels
                    .last()
                    .is_none_or(|previous| (level - previous).abs() > delta)
                {
                    levels.push(level);
                    if levels.len() > MAX_RENDERED_ROC_THRESHOLDS {
                        bail!(
                            "ROC output exceeds the {MAX_RENDERED_ROC_THRESHOLDS} rendered-threshold resource limit"
                        );
                    }
                }
            }
        }
    }
    Ok(kept)
}

#[derive(Default)]
struct LegacyUnorderedRows {
    bucket_count: usize,
    buckets: BTreeMap<usize, VecDeque<String>>,
    bucket_order: VecDeque<usize>,
    seen: HashSet<String>,
    typed: HashSet<String>,
}

impl LegacyUnorderedRows {
    fn set(&mut self, key: String, has_type: bool) -> Result<()> {
        if self.seen.contains(&key) {
            if has_type {
                self.typed.insert(key);
            }
            return Ok(());
        }
        if self.seen.len() >= MAX_ROC_METRIC_INDEX_KEYS {
            bail!("ROC metric index exceeds the {MAX_ROC_METRIC_INDEX_KEYS} key resource limit");
        }
        if self.bucket_count == 0 {
            self.bucket_count = 1;
        }
        if self.seen.len() + 1 > self.bucket_count {
            self.rehash(next_legacy_bucket(self.bucket_count, self.seen.len() + 1));
        }
        self.seen.insert(key.clone());
        if has_type {
            self.typed.insert(key.clone());
        }
        let bucket = (legacy_string_hash(&key) % self.bucket_count as u64) as usize;
        if let Some(entries) = self.buckets.get_mut(&bucket) {
            entries.push_front(key);
        } else {
            self.buckets.insert(bucket, VecDeque::from([key]));
            self.bucket_order.push_front(bucket);
        }
        Ok(())
    }

    fn rehash(&mut self, bucket_count: usize) {
        let old = self.all_rows();
        self.bucket_count = bucket_count;
        self.buckets.clear();
        self.bucket_order.clear();
        for key in old {
            let bucket = (legacy_string_hash(&key) % bucket_count as u64) as usize;
            if let Some(entries) = self.buckets.get_mut(&bucket) {
                entries.push_front(key);
            } else {
                self.buckets.insert(bucket, VecDeque::from([key]));
                self.bucket_order.push_front(bucket);
            }
        }
    }

    fn all_rows(&self) -> Vec<String> {
        self.bucket_order
            .iter()
            .flat_map(|bucket| self.buckets[bucket].iter().cloned())
            .collect()
    }

    fn retained_order(&self) -> Vec<String> {
        self.all_rows()
            .into_iter()
            .filter(|key| self.typed.contains(key))
            .collect()
    }
}

fn next_legacy_bucket(current: usize, required: usize) -> usize {
    const PRIMES: &[usize] = &[
        13, 29, 59, 127, 257, 541, 1109, 2357, 5087, 10273, 20753, 42043, 85229, 172933, 351061,
        712697, 1447153, 2938679, 5967347, 12117689, 24607243, 49969847,
    ];
    let target = required.max(current.saturating_mul(2));
    PRIMES
        .iter()
        .copied()
        .find(|prime| *prime >= target)
        .unwrap_or_else(|| target.next_power_of_two())
}

/// libstdc++ `_Hash_bytes`, used by `std::hash<std::string>` in the pinned
/// legacy image (64-bit MurmurHash2, seed 0xc70f6907).
fn legacy_string_hash(value: &str) -> u64 {
    const MULTIPLIER: u64 = 0xc6a4_a793_5bd1_e995;
    const SHIFT: u32 = 47;
    let bytes = value.as_bytes();
    let mut hash = 0xc70f_6907_u64 ^ (bytes.len() as u64).wrapping_mul(MULTIPLIER);
    let mut chunks = bytes.chunks_exact(8);
    for chunk in &mut chunks {
        let mut key = u64::from_le_bytes(chunk.try_into().expect("eight-byte hash chunk"));
        key = key.wrapping_mul(MULTIPLIER);
        key ^= key >> SHIFT;
        key = key.wrapping_mul(MULTIPLIER);
        hash ^= key;
        hash = hash.wrapping_mul(MULTIPLIER);
    }
    for (index, byte) in chunks.remainder().iter().enumerate() {
        hash ^= u64::from(*byte) << (index * 8);
    }
    if !chunks.remainder().is_empty() {
        hash = hash.wrapping_mul(MULTIPLIER);
    }
    hash ^= hash >> SHIFT;
    hash = hash.wrapping_mul(MULTIPLIER);
    hash ^ (hash >> SHIFT)
}

/// The tab-separated table emitted by the legacy C++ quantifier before
/// happyroc turns it into the public CSV reports. qfy normally unlinks this
/// intermediate; `--verbose` deliberately leaves it behind.
#[derive(Default)]
struct LegacyRawTable {
    order: LegacyUnorderedRows,
    rows: BTreeMap<String, BTreeMap<String, String>>,
}

impl LegacyRawTable {
    fn set(&mut self, row: &str, column: &str, value: String, has_type: bool) -> Result<()> {
        self.order
            .set(row.to_string(), has_type)
            .with_context(|| format!("failed to admit legacy ROC row '{row}'"))?;
        self.rows
            .entry(row.to_string())
            .or_default()
            .insert(column.to_string(), value);
        Ok(())
    }

    fn write(&self, path: &Path) -> Result<()> {
        let retained = self.order.retained_order();
        let columns = retained
            .iter()
            .filter_map(|key| self.rows.get(key))
            .flat_map(|row| row.keys().cloned())
            .collect::<BTreeSet<_>>();
        let mut writer = BufWriter::new(
            std::fs::File::create(path)
                .with_context(|| format!("failed to create {}", path.display()))?,
        );
        writeln!(
            writer,
            "{}",
            columns.iter().cloned().collect::<Vec<_>>().join("\t")
        )?;
        for key in retained {
            let row = &self.rows[&key];
            writeln!(
                writer,
                "{}",
                columns
                    .iter()
                    .map(|column| row.get(column).map(String::as_str).unwrap_or("."))
                    .collect::<Vec<_>>()
                    .join("\t")
            )?;
        }
        writer
            .flush()
            .with_context(|| format!("failed to flush {}", path.display()))
    }
}

fn write_legacy_roc_table(
    prefix: &Path,
    groups: &BTreeMap<RowKey, BoundedGroupAccum>,
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<()> {
    // ROCOutput iterates a std::map keyed by its internal ROC name, not by
    // the final report axes. Recreate those names so insertion/rehash order
    // in LegacyUnorderedRows is byte-identical to libstdc++.
    let mut rocs = BTreeMap::<String, (&RowKey, &BoundedGroupAccum)>::new();
    for (key, accum) in groups {
        if key.subtype != "*" {
            continue;
        }
        let name = if key.subset != "*" {
            format!("s|{}:{}:{}", key.subset, key.ty, key.filter)
        } else if is_aggregate_filter(&key.filter) {
            format!("a:{}:{}", key.ty, key.filter)
        } else {
            format!("f:{}:{}", key.ty, key.filter)
        };
        rocs.insert(name, (key, accum));
    }

    let mut table = LegacyRawTable::default();
    for (_, (key, accum)) in rocs {
        let subtypes: &[&str] = if key.ty == "SNP" {
            &["*", "ti", "tv"]
        } else {
            &[
                "*", "I1_5", "I6_15", "I16_PLUS", "D1_5", "D6_15", "D16_PLUS", "C1_5", "C6_15",
                "C16_PLUS",
            ]
        };
        let counts_only =
            !is_aggregate_filter(&key.filter) && options.roc_regions.contains(key.subset.as_str());
        // getLevels sorts the shared observation vector in-place on every
        // subtype/genotype call. Preserve that progressive tied-level order.
        let mut sorted_obs = accum.observations.small_records();
        let mut sort_passes = 0usize;
        for subtype in subtypes {
            for genotype in ["het", "hetalt", "homalt", "*"] {
                let totals = legacy_totals_store(&accum.observations, subtype, genotype)?;
                add_legacy_level(
                    &mut table,
                    key,
                    RowLabel {
                        subtype,
                        genotype,
                        qq: "*",
                    },
                    &totals,
                    counts_only,
                    &RowSizes {
                        subset_size,
                        whole_reference_size: options.whole_reference_size.unwrap_or(subset_size),
                        conf_size,
                        subset_sizes: &options.subset_sizes,
                        subset_confidence_sizes: &options.subset_confidence_sizes,
                    },
                )?;
                if !counts_only && options.output_rocs {
                    sort_passes += 1;
                    let levels = if accum.observations.is_disk_backed() {
                        legacy_levels_store(
                            &accum.observations,
                            subtype,
                            genotype,
                            options.delta,
                            sort_passes,
                        )?
                    } else {
                        introsort_libstdcpp(&mut sorted_obs);
                        legacy_levels(&sorted_obs, subtype, genotype, options.delta)
                    };
                    for (qq, level) in levels {
                        add_legacy_level(
                            &mut table,
                            key,
                            RowLabel {
                                subtype,
                                genotype,
                                qq: &qq,
                            },
                            &level,
                            false,
                            &RowSizes {
                                subset_size,
                                whole_reference_size: options
                                    .whole_reference_size
                                    .unwrap_or(subset_size),
                                conf_size,
                                subset_sizes: &options.subset_sizes,
                                subset_confidence_sizes: &options.subset_confidence_sizes,
                            },
                        )?;
                    }
                }
            }
        }
    }
    table.write(&suffixed_report_path(prefix, "roc.tsv"))
}

fn legacy_obs_matches(record: &ObsRecord, subtype: &str, genotype: &str) -> bool {
    let subtype_matches = match subtype {
        "*" => true,
        "ti" => record.ti_flag,
        "tv" => record.tv_flag,
        value => record.has_subtype(value),
    };
    subtype_matches && (genotype == "*" || record.blt_is(genotype))
}

fn legacy_totals(obs: &[ObsRecord], subtype: &str, genotype: &str) -> Cumul {
    let mut totals = Cumul::default();
    for record in obs
        .iter()
        .filter(|record| legacy_obs_matches(record, subtype, genotype))
    {
        totals.add(&record.counts);
    }
    totals
}

fn legacy_totals_store(
    observations: &ObservationStore,
    subtype: &str,
    genotype: &str,
) -> Result<Cumul> {
    if !observations.is_disk_backed() {
        return Ok(legacy_totals(
            &observations.small_records(),
            subtype,
            genotype,
        ));
    }
    let mut totals = Cumul::default();
    for observation in observations.sorted()? {
        let observation = observation?;
        if legacy_obs_matches(&observation, subtype, genotype) {
            totals.add(&observation.counts);
        }
    }
    Ok(totals)
}

fn legacy_levels(
    sorted_obs: &[ObsRecord],
    subtype: &str,
    genotype: &str,
    delta: f64,
) -> Vec<(String, Cumul)> {
    let mut running = Cumul::default();
    let mut prefixes = Vec::new();
    for record in sorted_obs
        .iter()
        .filter(|record| legacy_obs_matches(record, subtype, genotype))
    {
        running.add(&record.counts);
        prefixes.push((record.level, running.clone()));
    }
    let Some((_, final_counts)) = prefixes.last().cloned() else {
        return Vec::new();
    };
    let mut candidates = prefixes
        .into_iter()
        .map(|(level, below)| {
            let above = Cumul {
                truth_tp: sub_buckets(&final_counts.truth_tp, &below.truth_tp),
                truth_fn: sum_buckets(&final_counts.truth_fn, &below.truth_tp),
                query_tp: sub_buckets(&final_counts.query_tp, &below.query_tp),
                query_fp: sub_buckets(&final_counts.query_fp, &below.query_fp),
                query_unk: sub_buckets(&final_counts.query_unk, &below.query_unk),
                fp_gt: final_counts.fp_gt.saturating_sub(below.fp_gt),
                fp_al: final_counts.fp_al.saturating_sub(below.fp_al),
            };
            (level, above)
        })
        .collect::<Vec<_>>();
    let first = candidates[0].0;
    let mut previous = first;
    let mut is_first = true;
    candidates.retain(|(level, _)| {
        let keep = if delta < f64::EPSILON {
            format!("{level:.6}") != format!("{previous:.6}")
        } else if is_first {
            is_first = false;
            true
        } else {
            (*level - previous).abs() > delta
        };
        if keep {
            previous = *level;
        }
        keep
    });
    candidates
        .into_iter()
        .map(|(level, counts)| (format!("{level:.6}"), counts))
        .collect()
}

fn legacy_levels_store(
    observations: &ObservationStore,
    subtype: &str,
    genotype: &str,
    delta: f64,
    sort_passes: usize,
) -> Result<Vec<(String, Cumul)>> {
    let totals = legacy_totals_store(observations, subtype, genotype)?;
    let mut running = Cumul::default();
    let mut previous = None;
    let mut levels = Vec::new();
    for observation in observations.sorted_with_passes(sort_passes)? {
        let observation = observation?;
        if !legacy_obs_matches(&observation, subtype, genotype) {
            continue;
        }
        running.add(&observation.counts);
        if previous.is_none_or(|value: f64| (observation.level - value).abs() > delta) {
            previous = Some(observation.level);
            levels.push((
                format!("{:.6}", observation.level as f32 as f64),
                Cumul {
                    truth_tp: sub_buckets(&totals.truth_tp, &running.truth_tp),
                    truth_fn: sum_buckets(&totals.truth_fn, &running.truth_tp),
                    query_tp: sub_buckets(&totals.query_tp, &running.query_tp),
                    query_fp: sub_buckets(&totals.query_fp, &running.query_fp),
                    query_unk: sub_buckets(&totals.query_unk, &running.query_unk),
                    fp_gt: totals.fp_gt.saturating_sub(running.fp_gt),
                    fp_al: totals.fp_al.saturating_sub(running.fp_al),
                },
            ));
            if levels.len() > MAX_RENDERED_ROC_THRESHOLDS {
                bail!(
                    "verbose ROC table exceeds the {MAX_RENDERED_ROC_THRESHOLDS} rendered-threshold resource limit"
                );
            }
        }
    }
    Ok(levels)
}

/// The reference/confidence sizing columns every ROC row renders against.
struct RowSizes<'a> {
    subset_size: usize,
    whole_reference_size: usize,
    conf_size: usize,
    subset_sizes: &'a BTreeMap<String, usize>,
    subset_confidence_sizes: &'a BTreeMap<String, usize>,
}

/// The subtype/genotype/QQ triple identifying one legacy raw-table level.
struct RowLabel<'a> {
    subtype: &'a str,
    genotype: &'a str,
    qq: &'a str,
}

fn add_legacy_level(
    table: &mut LegacyRawTable,
    key: &RowKey,
    label: RowLabel<'_>,
    counts: &Cumul,
    counts_only: bool,
    sizes: &RowSizes<'_>,
) -> Result<()> {
    let RowLabel {
        subtype,
        genotype,
        qq,
    } = label;
    let &RowSizes {
        subset_size,
        whole_reference_size,
        conf_size,
        subset_sizes,
        subset_confidence_sizes,
    } = sizes;
    let is_baseline = qq == "*";
    let row_key = if is_baseline || genotype == "*" {
        legacy_row_key(&key.ty, subtype, &key.filter, &key.subset, qq)
    } else {
        // Retain ROCOutput.cpp's historical extra-tab bug. These rows are
        // later dropped for lacking Type, but can trigger a table rehash.
        format!(
            "{}\t{}\t*\t\t{}\t{}\t{}",
            key.ty, subtype, key.filter, key.subset, qq
        )
    };

    if !matches!(subtype, "ti" | "tv") && genotype == "*" {
        table.set(&row_key, "QQ", qq.to_string(), true)?;
        for (column, value) in [
            ("Type", key.ty.as_str()),
            ("Subtype", subtype),
            ("Genotype", genotype),
            ("Subset", key.subset.as_str()),
            ("Filter", key.filter.as_str()),
            ("QQ.Field", key.qq_field.as_str()),
        ] {
            table.set(&row_key, column, value.to_string(), true)?;
        }
        for (column, value) in legacy_primary_counts(counts, counts_only) {
            table.set(&row_key, column, value, true)?;
        }
        let (size, conf) = subset_size_cells(
            &key.subset,
            subtype,
            subset_size,
            whole_reference_size,
            conf_size,
            subset_sizes,
            subset_confidence_sizes,
        );
        table.set(&row_key, "Subset.Size", legacy_count_string(&size), true)?;
        table.set(
            &row_key,
            "Subset.IS_CONF.Size",
            legacy_count_string(&conf),
            true,
        )?;
        table.set(&row_key, "Subset.Level", "0.000000".to_string(), true)?;
    } else if genotype == "*" && matches!(subtype, "ti" | "tv") {
        let aggregate = legacy_row_key(&key.ty, "*", &key.filter, &key.subset, qq);
        for (metric, count) in legacy_count_buckets(counts, counts_only) {
            table.set(
                &aggregate,
                &format!("{metric}.{subtype}"),
                legacy_usize(count),
                false,
            )?;
        }
    } else if genotype != "*" && !matches!(subtype, "ti" | "tv") {
        for (metric, count) in legacy_count_buckets(counts, counts_only) {
            table.set(
                &row_key,
                &format!("{metric}.{genotype}"),
                legacy_usize(count),
                false,
            )?;
        }
    }
    Ok(())
}

fn legacy_count_buckets(counts: &Cumul, counts_only: bool) -> Vec<(&'static str, usize)> {
    let mut values = vec![
        ("TRUTH.TP", counts.truth_tp.total),
        ("QUERY.TP", counts.query_tp.total),
        ("QUERY.FP", counts.query_fp.total),
        ("QUERY.UNK", counts.query_unk.total),
    ];
    if !counts_only {
        values.extend([
            ("TRUTH.FN", counts.truth_fn.total),
            ("TRUTH.TOTAL", counts.truth_total().total),
            ("QUERY.TOTAL", counts.query_total().total),
        ]);
    }
    values
}

fn legacy_primary_counts(counts: &Cumul, counts_only: bool) -> Vec<(&'static str, String)> {
    let mut values = legacy_count_buckets(counts, counts_only)
        .into_iter()
        .map(|(column, value)| (column, legacy_usize(value)))
        .collect::<Vec<_>>();
    values.extend([
        ("FP.al", legacy_usize(counts.fp_al)),
        ("FP.gt", legacy_usize(counts.fp_gt)),
    ]);
    if !counts_only {
        let truth_total = counts.truth_total().total;
        let query_total = counts.query_total().total;
        let recall = if truth_total == 0 {
            0.0
        } else {
            counts.truth_tp.total as f64 / truth_total as f64
        };
        let precision = if query_total == 0 {
            0.0
        } else {
            counts.query_tp.total as f64 / (counts.query_tp.total + counts.query_fp.total) as f64
        };
        let frac_na = if query_total == 0 {
            0.0
        } else {
            counts.query_unk.total as f64 / query_total as f64
        };
        let f1 = 2.0 * precision * recall / (precision + recall);
        values.extend([
            ("METRIC.Recall", legacy_f64(recall)),
            ("METRIC.Precision", legacy_f64(precision)),
            ("METRIC.F1_Score", legacy_f64(f1)),
            ("METRIC.Frac_NA", legacy_f64(frac_na)),
        ]);
    }
    values
}

fn legacy_usize(value: usize) -> String {
    format!("{:.6}", value as f64)
}

fn legacy_count_string(value: &str) -> String {
    value
        .parse::<f64>()
        .map(legacy_f64)
        .unwrap_or_else(|_| "0.000000".to_string())
}

fn legacy_f64(value: f64) -> String {
    if value.is_nan() {
        "-nan".to_string()
    } else {
        format!("{value:.6}")
    }
}

// ---------------------------------------------------------------------------
// Grouping
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct RowKey {
    // Legacy sort order: (Type, Subtype, Subset, Filter, Genotype, QQ.Field)
    // all ascending. Derive Ord via field order.
    ty: String,
    subtype: String,
    subset: String,
    filter: String,
    genotype: String,
    qq_field: String,
}

impl RowKey {
    #[cfg(test)]
    fn new(ty: &str, subtype: &str, subset: &str, filter: &str) -> Self {
        Self::new_with_qq_field(ty, subtype, subset, filter, "QUAL")
    }

    fn new_with_qq_field(
        ty: &str,
        subtype: &str,
        subset: &str,
        filter: &str,
        qq_field: &str,
    ) -> Self {
        Self {
            ty: ty.to_string(),
            subtype: subtype.to_string(),
            subset: subset.to_string(),
            filter: filter.to_string(),
            genotype: "*".to_string(),
            qq_field: qq_field.to_string(),
        }
    }
}

/// Aggregated counts for one cumulative ROC row (or baseline).
#[derive(Clone, Debug, Default)]
struct Cumul {
    truth_tp: CountsBucket,
    truth_fn: CountsBucket,
    query_tp: CountsBucket,
    query_fp: CountsBucket,
    query_unk: CountsBucket,
    fp_gt: usize,
    fp_al: usize,
}

/// Per-record observation that mirrors a single legacy `roc::Observation`
/// entry. Each call into `emit_contributions` produces exactly one such
/// observation (truth-side TP/FN, or query-side TP2/FP/UNK), which we
/// then sort + walk via the legacy `Roc::getLevels` algorithm to compute
/// the strict-above cumulative for each kept QQ row.
///
/// `level` carries the f32-rounded QQ value so that obs at "the same"
/// QQ from a parity standpoint cluster identically to legacy. `counts`
/// is a single-record contribution (exactly one of the truth/query
/// fields is non-zero). `subtypes` is the set of subtype strings this
/// obs contributes to (always includes `"*"`); legacy's per-Type ROC
/// vector groups all subtypes together and uses flag-mask filtering
/// during `getLevels`, so we mirror that by storing every obs once at
/// the (Type, Subset, Filter) level and filtering by subtype during
/// emission. This is what reproduces legacy's "single sort, walk per
/// subtype" semantics — sorting per-subtype-slice independently
/// produces different cluster orderings across slices and breaks
/// bit-parity on the per-subtype rows.
#[derive(Clone, Debug)]
struct ObsRecord {
    level: f64,
    counts: Cumul,
    subtype_bits: u16,
    /// True if the original BI tag included `ti` — needed to mirror
    /// legacy's `OBS_FLAG_TI` mask which is set per-record from the BI
    /// string, *not* from non-zero count fields. Filter-failed query
    /// phantoms have all-zero counts but still carry their flag bits.
    ti_flag: bool,
    tv_flag: bool,
    blt: u8,
}

impl ObsRecord {
    fn has_subtype(&self, subtype: &str) -> bool {
        subtype_bit(subtype).is_some_and(|bit| self.subtype_bits & bit != 0)
    }

    fn blt_is(&self, value: &str) -> bool {
        self.blt == encode_blt(Some(value))
    }
}

fn subtype_bit(subtype: &str) -> Option<u16> {
    if subtype == "*" {
        Some(1)
    } else {
        INDEL_SUBTYPES
            .iter()
            .position(|candidate| *candidate == subtype)
            .map(|index| 1u16 << (index + 1))
    }
}

fn encode_subtypes(subtypes: &[String]) -> u16 {
    subtypes
        .iter()
        .filter_map(|subtype| subtype_bit(subtype))
        .fold(0, |bits, bit| bits | bit)
}

fn subtype_names(bits: u16) -> impl Iterator<Item = &'static str> {
    std::iter::once("*")
        .chain(INDEL_SUBTYPES)
        .enumerate()
        .filter_map(move |(index, subtype)| (bits & (1u16 << index) != 0).then_some(subtype))
}

fn encode_blt(blt: Option<&str>) -> u8 {
    match blt {
        Some("het") => 1,
        Some("hetalt") => 2,
        Some("homalt") => 3,
        _ => 0,
    }
}

fn decode_blt(blt: u8) -> &'static str {
    match blt {
        1 => "het",
        2 => "hetalt",
        3 => "homalt",
        _ => ".",
    }
}

#[derive(Clone, Debug)]
struct StoredObs {
    serial: u64,
    record: ObsRecord,
}

#[derive(Default)]
struct ObservationStore {
    buffer: Vec<StoredObs>,
    chunks: Vec<tempfile::TempPath>,
    disk: Option<DiskObservationStore>,
    sorted_indices: std::sync::Mutex<BTreeMap<usize, Arc<Vec<DiskIndexEntry>>>>,
    disk_sorted_indices: std::sync::Mutex<BTreeMap<usize, tempfile::TempPath>>,
    next_serial: u64,
    error: Option<String>,
}

struct DiskObservationStore {
    data: tempfile::TempPath,
    index: tempfile::TempPath,
    len: usize,
}

impl ObservationStore {
    fn push(&mut self, record: ObsRecord) {
        self.buffer.push(StoredObs {
            serial: self.next_serial,
            record,
        });
        self.next_serial += 1;
        if self.buffer.len() >= ROC_OBSERVATION_CHUNK
            && let Err(error) = self.flush_chunk()
        {
            self.error.get_or_insert_with(|| error.to_string());
        }
    }

    fn flush_chunk(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        let mut chunk =
            tempfile::NamedTempFile::new().context("failed to create ROC observation chunk")?;
        {
            let mut writer = BufWriter::new(chunk.as_file_mut());
            for observation in self.buffer.drain(..) {
                write_observation(&mut writer, &observation)?;
            }
            writer.flush()?;
        }
        self.chunks.push(chunk.into_temp_path());
        Ok(())
    }

    fn finish(&mut self) -> Result<()> {
        if let Some(error) = self.error.take() {
            bail!("failed to spool ROC observations: {error}");
        }
        if !self.chunks.is_empty() {
            self.flush_chunk()?;
            let mut data = tempfile::NamedTempFile::new()
                .context("failed to create ROC observation data spool")?;
            let mut index = tempfile::NamedTempFile::new()
                .context("failed to create ROC observation index spool")?;
            let mut len = 0usize;
            for chunk in std::mem::take(&mut self.chunks) {
                for line in BufReader::new(File::open(&chunk)?).lines() {
                    let line = line?;
                    let observation = parse_observation(&line)?;
                    let encoded = format!("{line}\n");
                    data.as_file_mut().write_all(encoded.as_bytes())?;
                    write_disk_index_entry(
                        index.as_file_mut(),
                        DiskIndexEntry {
                            observation_bits: encode_compact_observation(&observation.record)?,
                            level_bits: observation.record.level.to_bits(),
                        },
                    )?;
                    len += 1;
                }
            }
            data.as_file_mut().flush()?;
            index.as_file_mut().flush()?;
            self.disk = Some(DiskObservationStore {
                data: data.into_temp_path(),
                index: index.into_temp_path(),
                len,
            });
        }
        Ok(())
    }

    fn is_disk_backed(&self) -> bool {
        self.disk.is_some() || !self.chunks.is_empty()
    }

    fn small_records(&self) -> Vec<ObsRecord> {
        self.buffer
            .iter()
            .map(|observation| observation.record.clone())
            .collect()
    }

    fn unsorted(&self) -> Result<Box<dyn Iterator<Item = Result<ObsRecord>> + '_>> {
        if let Some(disk) = &self.disk {
            let lines = BufReader::new(File::open(&disk.data)?).lines();
            return Ok(Box::new(lines.map(|line| {
                parse_observation(&line?).map(|observation| observation.record)
            })));
        }
        Ok(Box::new(
            self.buffer
                .iter()
                .map(|observation| Ok(observation.record.clone())),
        ))
    }

    fn sorted(&self) -> Result<Box<dyn Iterator<Item = Result<ObsRecord>>>> {
        self.sorted_with_passes(4)
    }

    fn sorted_with_passes(
        &self,
        passes: usize,
    ) -> Result<Box<dyn Iterator<Item = Result<ObsRecord>>>> {
        if let Some(disk) = &self.disk {
            let index_bytes = (disk.len as u64).saturating_mul(ROC_INDEX_ENTRY_BYTES);
            if index_bytes <= ROC_IN_MEMORY_INDEX_LIMIT_BYTES {
                let cached = self
                    .sorted_indices
                    .lock()
                    .unwrap_or_else(|error| error.into_inner())
                    .range(..=passes)
                    .next_back()
                    .map(|(completed, entries)| (*completed, Arc::clone(entries)));
                let (completed_passes, mut entries) = if let Some((completed, entries)) = cached {
                    (completed, entries.as_ref().clone())
                } else {
                    let mut index_file = File::open(&disk.index)?;
                    (0, read_disk_index_entries(&mut index_file, disk.len)?)
                };
                for _ in completed_passes..passes {
                    introsort_libstdcpp_by(&mut entries, |left, right| {
                        f64::from_bits(left.level_bits) < f64::from_bits(right.level_bits)
                    });
                }
                let entries = Arc::new(entries);
                let mut sorted_indices = self
                    .sorted_indices
                    .lock()
                    .unwrap_or_else(|error| error.into_inner());
                sorted_indices
                    .entry(passes)
                    .or_insert_with(|| Arc::clone(&entries));
                // Passes 4/8/12 feed the main, ti, and tv sweeps. Later
                // subtype passes are requested monotonically, so retaining
                // only the newest one bounds each worker to four snapshots.
                sorted_indices.retain(|completed, _| *completed <= 12 || *completed == passes);
                drop(sorted_indices);
                return Ok(Box::new(MemoryObservationIter {
                    entries,
                    position: 0,
                }));
            }

            let (completed_passes, source) = self
                .disk_sorted_indices
                .lock()
                .unwrap_or_else(|error| error.into_inner())
                .range(..=passes)
                .next_back()
                .map(|(completed, path)| (*completed, path.to_path_buf()))
                .unwrap_or_else(|| (0, disk.index.to_path_buf()));
            let index = copy_temp_path(&source)?;
            let mut index_file = File::options().read(true).write(true).open(&index)?;
            for _ in completed_passes..passes {
                disk_introsort_libstdcpp(&mut index_file, disk.len)?;
            }
            let mut sorted_indices = self
                .disk_sorted_indices
                .lock()
                .unwrap_or_else(|error| error.into_inner());
            if !sorted_indices.contains_key(&passes) {
                sorted_indices.insert(passes, copy_temp_path(&index)?);
            }
            return Ok(Box::new(DiskObservationIter {
                _index_path: index,
                sorted_index: index_file,
                position: 0,
                len: disk.len,
            }));
        }
        let mut observations = self.small_records();
        for _ in 0..passes {
            introsort_libstdcpp(&mut observations);
        }
        Ok(Box::new(observations.into_iter().map(Ok)))
    }

    fn clear_sorted_indices(&self) {
        self.sorted_indices
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clear();
        self.disk_sorted_indices
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clear();
    }
}

fn copy_temp_path(path: &Path) -> Result<tempfile::TempPath> {
    let mut copy = tempfile::NamedTempFile::new().context("failed to copy ROC spool")?;
    std::io::copy(&mut File::open(path)?, copy.as_file_mut())?;
    copy.as_file_mut().flush()?;
    Ok(copy.into_temp_path())
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct DiskIndexEntry {
    observation_bits: u64,
    level_bits: u64,
}

fn encode_compact_observation(observation: &ObsRecord) -> Result<u64> {
    let mut encoded = 0u64;
    let mut shift = 0u32;
    let mut push_count = |value: usize| -> Result<()> {
        if value > 1 {
            bail!("single ROC observation count {value} cannot be compacted");
        }
        encoded |= (value as u64) << shift;
        shift += 1;
        Ok(())
    };
    for bucket in [
        &observation.counts.truth_tp,
        &observation.counts.truth_fn,
        &observation.counts.query_tp,
        &observation.counts.query_fp,
        &observation.counts.query_unk,
    ] {
        for value in [
            bucket.total,
            bucket.ti,
            bucket.tv,
            bucket.het,
            bucket.homalt,
        ] {
            push_count(value)?;
        }
    }
    push_count(observation.counts.fp_gt)?;
    push_count(observation.counts.fp_al)?;
    drop(push_count);

    encoded |= u64::from(observation.subtype_bits) << shift;
    shift += 10;
    if observation.ti_flag {
        encoded |= 1u64 << shift;
    }
    shift += 1;
    if observation.tv_flag {
        encoded |= 1u64 << shift;
    }
    shift += 1;
    encoded |= u64::from(observation.blt) << shift;
    Ok(encoded)
}

fn decode_compact_observation(observation_bits: u64, level_bits: u64) -> Result<ObsRecord> {
    let mut shift = 0u32;
    let mut next_count = || {
        let value = ((observation_bits >> shift) & 1) as usize;
        shift += 1;
        value
    };
    let mut next_bucket = || CountsBucket {
        total: next_count(),
        ti: next_count(),
        tv: next_count(),
        het: next_count(),
        homalt: next_count(),
    };
    let truth_tp = next_bucket();
    let truth_fn = next_bucket();
    let query_tp = next_bucket();
    let query_fp = next_bucket();
    let query_unk = next_bucket();
    drop(next_bucket);
    let counts = Cumul {
        truth_tp,
        truth_fn,
        query_tp,
        query_fp,
        query_unk,
        fp_gt: next_count(),
        fp_al: next_count(),
    };
    drop(next_count);

    let subtype_bits = ((observation_bits >> shift) & 0x03ff) as u16;
    shift += 10;
    let ti_flag = observation_bits & (1u64 << shift) != 0;
    shift += 1;
    let tv_flag = observation_bits & (1u64 << shift) != 0;
    shift += 1;
    let blt = ((observation_bits >> shift) & 0b11) as u8;
    Ok(ObsRecord {
        level: f64::from_bits(level_bits),
        counts,
        subtype_bits,
        ti_flag,
        tv_flag,
        blt,
    })
}

fn write_disk_index_entry(file: &mut File, entry: DiskIndexEntry) -> Result<()> {
    file.write_all(&entry.observation_bits.to_le_bytes())?;
    file.write_all(&entry.level_bits.to_le_bytes())?;
    Ok(())
}

fn read_disk_index_entries(file: &mut File, len: usize) -> Result<Vec<DiskIndexEntry>> {
    file.seek(SeekFrom::Start(0))?;
    let byte_len = len
        .checked_mul(ROC_INDEX_ENTRY_BYTES as usize)
        .context("ROC index is too large to address")?;
    let mut bytes = vec![0u8; byte_len];
    file.read_exact(&mut bytes)?;
    Ok(bytes
        .chunks_exact(ROC_INDEX_ENTRY_BYTES as usize)
        .map(|entry| DiskIndexEntry {
            observation_bits: u64::from_le_bytes(
                entry[0..8].try_into().expect("observation bytes"),
            ),
            level_bits: u64::from_le_bytes(entry[8..16].try_into().expect("level bytes")),
        })
        .collect())
}

#[cfg(test)]
fn write_disk_index_entries(file: &mut File, entries: &[DiskIndexEntry]) -> Result<()> {
    file.seek(SeekFrom::Start(0))?;
    for entry in entries {
        write_disk_index_entry(file, *entry)?;
    }
    file.set_len((entries.len() as u64).saturating_mul(ROC_INDEX_ENTRY_BYTES))?;
    file.flush()?;
    Ok(())
}

fn read_disk_index_entry(file: &mut File, position: usize) -> Result<DiskIndexEntry> {
    file.seek(SeekFrom::Start(position as u64 * ROC_INDEX_ENTRY_BYTES))?;
    let mut bytes = [0u8; ROC_INDEX_ENTRY_BYTES as usize];
    file.read_exact(&mut bytes)?;
    Ok(DiskIndexEntry {
        observation_bits: u64::from_le_bytes(bytes[0..8].try_into().expect("observation bytes")),
        level_bits: u64::from_le_bytes(bytes[8..16].try_into().expect("level bytes")),
    })
}

fn write_disk_index_at(file: &mut File, position: usize, entry: DiskIndexEntry) -> Result<()> {
    file.seek(SeekFrom::Start(position as u64 * ROC_INDEX_ENTRY_BYTES))?;
    write_disk_index_entry(file, entry)
}

fn swap_disk_index(file: &mut File, left: usize, right: usize) -> Result<()> {
    if left == right {
        return Ok(());
    }
    let left_entry = read_disk_index_entry(file, left)?;
    let right_entry = read_disk_index_entry(file, right)?;
    write_disk_index_at(file, left, right_entry)?;
    write_disk_index_at(file, right, left_entry)
}

fn disk_level(file: &mut File, position: usize) -> Result<f64> {
    Ok(f64::from_bits(
        read_disk_index_entry(file, position)?.level_bits,
    ))
}

struct DiskObservationIter {
    _index_path: tempfile::TempPath,
    sorted_index: File,
    position: usize,
    len: usize,
}

struct MemoryObservationIter {
    entries: Arc<Vec<DiskIndexEntry>>,
    position: usize,
}

impl Iterator for MemoryObservationIter {
    type Item = Result<ObsRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        let entry = self.entries.get(self.position).copied()?;
        self.position += 1;
        Some(decode_compact_observation(
            entry.observation_bits,
            entry.level_bits,
        ))
    }
}

impl Iterator for DiskObservationIter {
    type Item = Result<ObsRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.len {
            return None;
        }
        let result = (|| {
            let entry = read_disk_index_entry(&mut self.sorted_index, self.position)?;
            self.position += 1;
            decode_compact_observation(entry.observation_bits, entry.level_bits)
        })();
        Some(result)
    }
}

fn sortable_f64(value: f64) -> u64 {
    let bits = value.to_bits();
    if bits & (1 << 63) == 0 {
        bits ^ (1 << 63)
    } else {
        !bits
    }
}

fn write_observation(writer: &mut dyn Write, observation: &StoredObs) -> Result<()> {
    let counts = &observation.record.counts;
    let subtypes = subtype_names(observation.record.subtype_bits)
        .collect::<Vec<_>>()
        .join("|");
    write!(
        writer,
        "{}\t{}\t{}",
        sortable_f64(observation.record.level),
        observation.serial,
        observation.record.level.to_bits()
    )?;
    for bucket in [
        &counts.truth_tp,
        &counts.truth_fn,
        &counts.query_tp,
        &counts.query_fp,
        &counts.query_unk,
    ] {
        write!(
            writer,
            "\t{}\t{}\t{}\t{}\t{}",
            bucket.total, bucket.ti, bucket.tv, bucket.het, bucket.homalt
        )?;
    }
    writeln!(
        writer,
        "\t{}\t{}\t{}\t{}\t{}\t{}",
        counts.fp_gt,
        counts.fp_al,
        subtypes,
        u8::from(observation.record.ti_flag),
        u8::from(observation.record.tv_flag),
        decode_blt(observation.record.blt)
    )?;
    Ok(())
}

fn parse_observation(line: &str) -> Result<StoredObs> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() != 34 {
        bail!(
            "ROC observation spool has {} fields, expected 34",
            fields.len()
        );
    }
    let mut index = 1usize;
    let serial = fields[index].parse()?;
    index += 1;
    let level = f64::from_bits(fields[index].parse()?);
    index += 1;
    let mut next_bucket = || -> Result<CountsBucket> {
        let bucket = CountsBucket {
            total: fields[index].parse()?,
            ti: fields[index + 1].parse()?,
            tv: fields[index + 2].parse()?,
            het: fields[index + 3].parse()?,
            homalt: fields[index + 4].parse()?,
        };
        index += 5;
        Ok(bucket)
    };
    let counts = Cumul {
        truth_tp: next_bucket()?,
        truth_fn: next_bucket()?,
        query_tp: next_bucket()?,
        query_fp: next_bucket()?,
        query_unk: next_bucket()?,
        fp_gt: fields[index].parse()?,
        fp_al: fields[index + 1].parse()?,
    };
    index += 2;
    let subtype_bits = fields[index]
        .split('|')
        .filter_map(subtype_bit)
        .fold(0u16, |bits, bit| bits | bit);
    let ti_flag = fields[index + 1] == "1";
    let tv_flag = fields[index + 2] == "1";
    let blt = encode_blt(Some(fields[index + 3]));
    Ok(StoredObs {
        serial,
        record: ObsRecord {
            level,
            counts,
            subtype_bits,
            ti_flag,
            tv_flag,
            blt,
        },
    })
}

impl Cumul {
    fn add(&mut self, other: &Cumul) {
        add_bucket(&mut self.truth_tp, &other.truth_tp);
        add_bucket(&mut self.truth_fn, &other.truth_fn);
        add_bucket(&mut self.query_tp, &other.query_tp);
        add_bucket(&mut self.query_fp, &other.query_fp);
        add_bucket(&mut self.query_unk, &other.query_unk);
        self.fp_gt += other.fp_gt;
        self.fp_al += other.fp_al;
    }

    fn truth_total(&self) -> CountsBucket {
        sum_buckets(&self.truth_tp, &self.truth_fn)
    }

    fn query_total(&self) -> CountsBucket {
        let a = sum_buckets(&self.query_tp, &self.query_fp);
        sum_buckets(&a, &self.query_unk)
    }
}

fn add_bucket(dst: &mut CountsBucket, src: &CountsBucket) {
    dst.total += src.total;
    dst.ti += src.ti;
    dst.tv += src.tv;
    dst.het += src.het;
    dst.homalt += src.homalt;
}

fn sum_buckets(a: &CountsBucket, b: &CountsBucket) -> CountsBucket {
    CountsBucket {
        total: a.total + b.total,
        ti: a.ti + b.ti,
        tv: a.tv + b.tv,
        het: a.het + b.het,
        homalt: a.homalt + b.homalt,
    }
}

/// Per-group accumulator. Tracks two parallel views:
///   • `baseline`: total counts across the group (used for `QQ="*"` row
///     and for grand-total subtraction in cumulative formulas).
///   • `obs`: one `ObsRecord` per legacy `addObs` call (one per
///     truth/query side per record). Sorted via libstdc++-equivalent
///     `introsort_libstdcpp`, walked to build `target` (cumulative
///     through each obs), then filtered with roc-delta to produce kept
///     ROC rows. This is the only path that can reproduce legacy's
///     std::sort instability on equal-level (truth, query) pairs at
///     bit-exact parity.
/// `threshold_window` is only the bounded, pre-spill fast path. Once the
/// observation chunk reaches `ROC_OBSERVATION_CHUNK`, observations are
/// indexed on disk and threshold/substat state is derived through bounded
/// streaming passes; no per-threshold map grows with a whole-genome input.
#[derive(Default)]
struct BoundedGroupAccum {
    baseline: Cumul,
    observations: ObservationStore,
    threshold_window: BTreeMap<String, NumericBucket>,
    records: usize,
}

struct NumericBucket {
    /// Numeric QQ value used for DESC cumulation ordering. All contributions
    /// that round to the same `{:.6}` string share a bucket and are merged
    /// (matches what the legacy pandas pipeline would do when it groups by
    /// the formatted QQ string).
    qq: f64,
    counts: Cumul,
    /// Substat-flag presence at this bucket. True if any obs (real or phantom)
    /// at this level was tagged with the corresponding bi. Used by
    /// substat_kept_levels to decide whether the substat sweep keeps this
    /// level. Mirrors legacy's BlockQuantify::observe which sets OBS_FLAG_TI/TV
    /// on the FN2 phantom obs from filter-failed query.TP records too.
    has_ti: bool,
    has_tv: bool,
}

/// Substat availability at a given numeric QQ row. Each Option carries the
/// cumulative count for that substat if the independent per-substat
/// roc-delta sweep kept this level; `None` means the substat cell renders
/// as the legacy na_rep (`.` for count cells, `""` for ratio cells).
#[derive(Clone, Debug, Default)]
struct SubstatAvail {
    ti: bool,
    tv: bool,
    het: bool,
    homalt: bool,
}

/// One row emitted by a `BoundedGroupAccum`. For the baseline row (qq_str = "*"),
/// `substats` is `None` and every substat cell is taken from `cum`. For
/// numeric QQ rows, `substats` marks which substat cells the independent
/// sweeps kept at this level; cells where the substat sweep didn't keep
/// render as na_rep.
#[derive(Clone, Debug)]
struct EmittedRow {
    qq_str: String,
    cum: Cumul,
    substats: Option<SubstatAvail>,
}

struct SubstatSnapshotCursor {
    _path: tempfile::TempPath,
    lines: std::io::Lines<BufReader<File>>,
    current: Option<ObsRecord>,
}

impl SubstatSnapshotCursor {
    fn build(
        source: &ObservationStore,
        sort_passes: usize,
        include: impl Fn(&ObsRecord) -> bool,
    ) -> Result<Self> {
        let mut file = tempfile::NamedTempFile::new()
            .context("failed to create ROC substat snapshot spool")?;
        let mut running = Cumul::default();
        let mut previous = None;
        let mut serial = 0u64;
        for observation in source.sorted_with_passes(sort_passes)? {
            let observation = observation?;
            if !include(&observation) {
                continue;
            }
            running.add(&observation.counts);
            if previous != Some(observation.level) {
                write_observation(
                    file.as_file_mut(),
                    &StoredObs {
                        serial,
                        record: ObsRecord {
                            level: observation.level,
                            counts: running.clone(),
                            subtype_bits: subtype_bit("*").expect("wildcard subtype is encoded"),
                            ti_flag: false,
                            tv_flag: false,
                            blt: 0,
                        },
                    },
                )?;
                previous = Some(observation.level);
                serial += 1;
            }
        }
        file.as_file_mut().flush()?;
        let path = file.into_temp_path();
        let mut cursor = Self {
            lines: BufReader::new(File::open(&path)?).lines(),
            current: None,
            _path: path,
        };
        cursor.advance()?;
        Ok(cursor)
    }

    fn advance(&mut self) -> Result<()> {
        self.current = self
            .lines
            .next()
            .transpose()?
            .map(|line| parse_observation(&line).map(|stored| stored.record))
            .transpose()?;
        Ok(())
    }

    fn counts_at(&mut self, level: f64) -> Result<Option<Cumul>> {
        while self
            .current
            .as_ref()
            .is_some_and(|record| record.level < level)
        {
            self.advance()?;
        }
        Ok(self
            .current
            .as_ref()
            .filter(|record| rendered_roc_level(record.level) == rendered_roc_level(level))
            .map(|record| record.counts.clone()))
    }
}

fn rendered_roc_level(level: f64) -> String {
    format!("{:.6}", level as f32 as f64)
}

enum EmittedRows {
    Memory(std::vec::IntoIter<EmittedRow>),
    Disk(Box<DiskEmittedRows>),
}

struct DiskEmittedRows {
    baseline: Option<EmittedRow>,
    merge: EmittedRowMerge,
}

impl Iterator for EmittedRows {
    type Item = Result<EmittedRow>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Memory(rows) => rows.next().map(Ok),
            Self::Disk(rows) => rows.baseline.take().map(Ok).or_else(|| rows.merge.next()),
        }
    }
}

struct EmittedRowSpool {
    buffer: Vec<(String, u64, EmittedRow)>,
    chunks: Vec<tempfile::TempPath>,
    serial: u64,
}

impl EmittedRowSpool {
    fn new() -> Self {
        Self {
            buffer: Vec::new(),
            chunks: Vec::new(),
            serial: 0,
        }
    }

    fn push(&mut self, row: EmittedRow) -> Result<()> {
        if self.serial as usize >= MAX_RENDERED_ROC_THRESHOLDS {
            bail!(
                "ROC output exceeds the {MAX_RENDERED_ROC_THRESHOLDS} rendered-threshold resource limit"
            );
        }
        self.buffer.push((row.qq_str.clone(), self.serial, row));
        self.serial += 1;
        if self.buffer.len() >= ROC_OBSERVATION_CHUNK {
            self.flush()?;
        }
        Ok(())
    }

    fn flush(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }
        self.buffer
            .sort_by(|left, right| left.0.cmp(&right.0).then(left.1.cmp(&right.1)));
        let mut chunk =
            tempfile::NamedTempFile::new().context("failed to create rendered ROC row chunk")?;
        {
            let mut writer = BufWriter::new(chunk.as_file_mut());
            for (_, serial, row) in self.buffer.drain(..) {
                write_emitted_row(&mut writer, serial, &row)?;
            }
            writer.flush()?;
        }
        self.chunks.push(chunk.into_temp_path());
        Ok(())
    }

    fn finish(mut self) -> Result<EmittedRowMerge> {
        self.flush()?;
        while self.chunks.len() > 32 {
            let mut merged = Vec::new();
            let mut remaining = self.chunks.into_iter();
            loop {
                let batch = remaining.by_ref().take(32).collect::<Vec<_>>();
                if batch.is_empty() {
                    break;
                }
                let mut output =
                    tempfile::NamedTempFile::new().context("failed to merge rendered ROC rows")?;
                {
                    let mut writer = BufWriter::new(output.as_file_mut());
                    let mut merge = EmittedRowMerge::open(batch)?;
                    while let Some(row) = merge.next_keyed() {
                        let (serial, row) = row?;
                        write_emitted_row(&mut writer, serial, &row)?;
                    }
                    writer.flush()?;
                }
                merged.push(output.into_temp_path());
            }
            self.chunks = merged;
        }
        EmittedRowMerge::open(self.chunks)
    }
}

struct EmittedRowMerge {
    _chunks: Vec<tempfile::TempPath>,
    readers: Vec<std::io::Lines<BufReader<File>>>,
    current: Vec<Option<(String, u64, EmittedRow)>>,
    heap: BinaryHeap<Reverse<(String, u64, usize)>>,
}

impl EmittedRowMerge {
    fn open(chunks: Vec<tempfile::TempPath>) -> Result<Self> {
        let readers = chunks
            .iter()
            .map(|path| Ok(BufReader::new(File::open(path)?).lines()))
            .collect::<Result<Vec<_>>>()?;
        let mut merge = Self {
            current: (0..readers.len()).map(|_| None).collect(),
            readers,
            heap: BinaryHeap::new(),
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
        let (serial, row) = parse_emitted_row(&line?)?;
        self.heap.push(Reverse((row.qq_str.clone(), serial, index)));
        self.current[index] = Some((row.qq_str.clone(), serial, row));
        Ok(())
    }

    fn next_keyed(&mut self) -> Option<Result<(u64, EmittedRow)>> {
        let Reverse((_, _, index)) = self.heap.pop()?;
        let (_, serial, row) = self.current[index]
            .take()
            .expect("rendered ROC merge entry has a row");
        if let Err(error) = self.advance(index) {
            return Some(Err(error));
        }
        Some(Ok((serial, row)))
    }
}

impl Iterator for EmittedRowMerge {
    type Item = Result<EmittedRow>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_keyed().map(|row| row.map(|(_, row)| row))
    }
}

fn write_emitted_row(writer: &mut dyn Write, serial: u64, row: &EmittedRow) -> Result<()> {
    write!(writer, "{}\t{serial}", row.qq_str)?;
    for bucket in [
        &row.cum.truth_tp,
        &row.cum.truth_fn,
        &row.cum.query_tp,
        &row.cum.query_fp,
        &row.cum.query_unk,
    ] {
        write!(
            writer,
            "\t{}\t{}\t{}\t{}\t{}",
            bucket.total, bucket.ti, bucket.tv, bucket.het, bucket.homalt
        )?;
    }
    let substats = row.substats.clone().unwrap_or_default();
    writeln!(
        writer,
        "\t{}\t{}\t{}\t{}\t{}\t{}",
        row.cum.fp_gt,
        row.cum.fp_al,
        u8::from(substats.ti),
        u8::from(substats.tv),
        u8::from(substats.het),
        u8::from(substats.homalt)
    )?;
    Ok(())
}

fn parse_emitted_row(line: &str) -> Result<(u64, EmittedRow)> {
    let fields = line.split('\t').collect::<Vec<_>>();
    if fields.len() != 33 {
        bail!(
            "rendered ROC spool has {} fields, expected 33",
            fields.len()
        );
    }
    let mut index = 2usize;
    let mut bucket = || -> Result<CountsBucket> {
        let result = CountsBucket {
            total: fields[index].parse()?,
            ti: fields[index + 1].parse()?,
            tv: fields[index + 2].parse()?,
            het: fields[index + 3].parse()?,
            homalt: fields[index + 4].parse()?,
        };
        index += 5;
        Ok(result)
    };
    let cum = Cumul {
        truth_tp: bucket()?,
        truth_fn: bucket()?,
        query_tp: bucket()?,
        query_fp: bucket()?,
        query_unk: bucket()?,
        fp_gt: fields[index].parse()?,
        fp_al: fields[index + 1].parse()?,
    };
    index += 2;
    Ok((
        fields[1].parse()?,
        EmittedRow {
            qq_str: fields[0].to_string(),
            cum,
            substats: Some(SubstatAvail {
                ti: fields[index] == "1",
                tv: fields[index + 1] == "1",
                het: fields[index + 2] == "1",
                homalt: fields[index + 3] == "1",
            }),
        },
    ))
}

impl BoundedGroupAccum {
    fn add(
        &mut self,
        qq: Option<f64>,
        counts: &Cumul,
        subtypes: &[String],
        bi: Option<&str>,
        blt: Option<&str>,
        track_thresholds: bool,
    ) {
        self.baseline.add(counts);
        self.records += 1;
        if !track_thresholds {
            return;
        }
        // Legacy `BlockQuantify::observe` maps NaN-level obs to level=0
        // (`if(std::isnan(qq)) qq = 0`) before any further processing.
        // Mirror that by treating None / non-finite QQ as 0 for sorting
        // purposes (and bucket-key computation).
        let q = qq.filter(|q| q.is_finite()).unwrap_or(0.0);
        // Legacy's ROC table carries QQ values with f32 precision (its
        // upstream `quantify` binary stores thresholds as floats). Round
        // through f32 before formatting so neighbouring f64 inputs that
        // share an f32 representation collapse into the same threshold
        // bucket — matching legacy's threshold cardinality and its
        // characteristic 7-sig-fig formatted strings (`1001.340027`
        // rather than `1001.340000`).
        let q32 = q as f32 as f64;
        let key = rendered_roc_level(q32);
        if !self.observations.is_disk_backed() {
            let entry = self
                .threshold_window
                .entry(key)
                .or_insert_with(|| NumericBucket {
                    qq: q32,
                    counts: Cumul::default(),
                    has_ti: false,
                    has_tv: false,
                });
            entry.counts.add(counts);
            // Track substat-flag presence while the in-memory fast path is
            // active. Disk-backed groups derive the same flags while merging.
            if counts.truth_tp.ti > 0
                || counts.truth_fn.ti > 0
                || counts.query_tp.ti > 0
                || counts.query_fp.ti > 0
                || counts.query_unk.ti > 0
            {
                entry.has_ti = true;
            }
            if counts.truth_tp.tv > 0
                || counts.truth_fn.tv > 0
                || counts.query_tp.tv > 0
                || counts.query_fp.tv > 0
                || counts.query_unk.tv > 0
            {
                entry.has_tv = true;
            }
            if let Some(bi_str) = bi {
                for tok in bi_str.split(',') {
                    match tok {
                        "ti" => entry.has_ti = true,
                        "tv" => entry.has_tv = true,
                        _ => {}
                    }
                }
            }
        }
        // Per-record obs for the libstdc++-style sort + walk that
        // reproduces legacy's strict-above cumulative formula at
        // each kept QQ row. Legacy applies a level==0 remap for
        // observations whose `final_dt` is FN, FN2, or N (all
        // non-positive decision types): the obs is stored at
        // `numeric_limits<double>::min()` (a tiny positive subnormal
        // ≈ 2.225e-308) so it sorts immediately before any
        // positive-level obs — preventing FN/N obs from clustering
        // with real level=0 TP/TP2/FP/UNK obs. We identify FN/N
        // type by the populated count fields:
        //   • truth_fn populated → FN type (or filter-failed TP→FN)
        //   • all-zero counts   → phantom from filter-failed query
        //                          (final_dt is FN2 or N)
        let is_fn_n_type = counts.truth_fn.total > 0
            || (counts.truth_tp.total == 0
                && counts.query_tp.total == 0
                && counts.query_fp.total == 0
                && counts.query_unk.total == 0);
        let obs_level = if q32 == 0.0 && is_fn_n_type {
            f64::MIN_POSITIVE
        } else {
            q32
        };
        let mut ti_flag = false;
        let mut tv_flag = false;
        if let Some(bi_str) = bi {
            for tok in bi_str.split(',') {
                match tok {
                    "ti" => ti_flag = true,
                    "tv" => tv_flag = true,
                    _ => {}
                }
            }
        }
        // Records with non-zero .ti / .tv counts but no `bi` tag (rare —
        // shouldn't happen since sample_bucket only sets ti/tv from bi)
        // also count for substat sweeps. Mirror that conservatively.
        if !ti_flag
            && (counts.truth_tp.ti > 0
                || counts.truth_fn.ti > 0
                || counts.query_tp.ti > 0
                || counts.query_fp.ti > 0
                || counts.query_unk.ti > 0)
        {
            ti_flag = true;
        }
        if !tv_flag
            && (counts.truth_tp.tv > 0
                || counts.truth_fn.tv > 0
                || counts.query_tp.tv > 0
                || counts.query_fp.tv > 0
                || counts.query_unk.tv > 0)
        {
            tv_flag = true;
        }
        self.observations.push(ObsRecord {
            level: obs_level,
            counts: counts.clone(),
            subtype_bits: encode_subtypes(subtypes),
            ti_flag,
            tv_flag,
            blt: encode_blt(blt),
        });
        if self.observations.is_disk_backed() {
            self.threshold_window.clear();
        }
    }

    /// Emit (qq_string, cumulative_counts) pairs in the legacy output order:
    /// baseline `"*"` first, then numeric QQ strings sorted lex-ASC (which
    /// matches what pandas does when sorting an object column that mixes
    /// `"*"` with numeric strings).
    ///
    /// At each numeric threshold `t`, `TRUTH.TOTAL` is constant (equal to
    /// baseline's TP+FN), `TRUTH.TP` is the cumulative count of TPs with
    /// `qq >= t`, and `TRUTH.FN` is derived as `TRUTH.TOTAL − TRUTH.TP(t)`
    /// per substat field — matching the legacy semantics where TPs that
    /// fall below the threshold are reclassified as FN rather than vanishing.
    ///
    /// A legacy-compatible roc-delta filter of 0.5 (the default from
    /// `qfy.py`'s `--roc-delta`) drops rows whose QQ is within 0.5 of the
    /// last kept row when walking ASC from the lowest threshold — matching
    /// the C++ `Roc::getLevels` / `RocOutput` behaviour.
    ///
    /// Substat cells (ti/tv/het/homalt and their ratios) mirror legacy's
    /// independent-sweep pattern: each substat has its own roc-delta
    /// filter applied only to levels where that substat contributed. At a
    /// kept main-sweep level, a substat cell shows the cumulative value
    /// only if the substat's sweep also kept that exact level; otherwise
    /// `.` (count) or `""` (ratio). This mirrors legacy's `_addLevel` +
    /// `getLevels(flag_mask)` × `dropRowsWithMissing("Type")` pipeline
    /// where substat-only levels get filtered out and substat values are
    /// written only into prefixes the main sweep already visited.
    /// Standard emit — sorts the bounded in-memory observation chunk with
    /// libstdc++ introsort; spilled groups use the external merge path.
    /// Used for the `Subtype="*"` group and as a fallback for callers
    /// that don't have a pre-sorted obs vector to share.
    #[cfg(test)]
    fn emit(&self) -> Vec<EmittedRow> {
        self.emit_with_delta(0.5).expect("emit test rows")
    }

    #[cfg(test)]
    fn emit_with_delta(&self, delta: f64) -> Result<Vec<EmittedRow>> {
        self.stream_with_delta(delta)?.collect()
    }

    fn stream_with_delta(&self, delta: f64) -> Result<EmittedRows> {
        if self.observations.is_disk_backed() {
            self.emit_disk_with_delta_from(&self.observations, None, 4, delta)
        } else {
            Ok(EmittedRows::Memory(
                self.emit_internal(None, None, delta).into_iter(),
            ))
        }
    }

    /// Emit using a pre-sorted obs vector (typically from the `(*)`
    /// sibling in the same `(Type, Subset, Filter)` group). When
    /// `my_subtype` is non-empty and not `"*"`, the shared obs are
    /// filtered to entries whose `subtypes` field contains `my_subtype`
    /// — matching legacy's "single per-Type sort, walk per subtype
    /// with flag mask" semantics. This is what reproduces legacy's
    /// per-cluster cluster-ordering at non-(*)-subtype slices.
    fn emit_with_shared_sort_and_delta(
        &self,
        shared_sorted: &[ObsRecord],
        my_subtype: &str,
        delta: f64,
    ) -> Result<EmittedRows> {
        Ok(EmittedRows::Memory(
            self.emit_internal(Some(shared_sorted), Some(my_subtype), delta)
                .into_iter(),
        ))
    }

    fn emit_internal(
        &self,
        shared_sorted: Option<&[ObsRecord]>,
        my_subtype: Option<&str>,
        delta: f64,
    ) -> Vec<EmittedRow> {
        let truth_total_const = sum_buckets(&self.baseline.truth_tp, &self.baseline.truth_fn);

        let mut out = Vec::with_capacity(1 + self.threshold_window.len());
        out.push(EmittedRow {
            qq_str: "*".to_string(),
            cum: self.baseline.clone(),
            substats: None,
        });

        // ----- Main sweep (per-record obs + libstdc++ introsort) -----
        //
        // Legacy `Roc::getLevels` walks the obs vector AFTER `std::sort`
        // and snapshots `last` (running cumulative) into `target` after
        // adding each obs. After sort + roc-delta filter, the kept obs
        // at level L is the FIRST obs at L in sorted order; the row's
        // cumulative cells are computed via:
        //   tp  = total_tp  − target[first_at_L].tp()      (strict above)
        //   fn  = total_fn  + target[first_at_L].tp()      (TPs at-or-below)
        //   tp2 = total_tp2 − target[first_at_L].tp2()
        //   fp  = total_fp  − target[first_at_L].fp()
        //   unk = total_unk − target[first_at_L].unk()
        // Whether row.tp excludes truth obs at L (strict) or includes
        // them (legacy "first=query" outcome) hinges on whether the
        // first obs at the level cluster is a truth.TP or a query.TP2 —
        // an order determined by libstdc++'s `std::sort` partition
        // behavior. Reproducing that order bit-exactly requires the
        // `introsort_libstdcpp` port below.
        // The level==0 → MIN_POSITIVE remap for FN/N obs is now done
        // at insertion time (`BoundedGroupAccum::add`), matching legacy's
        // `BlockQuantify::observe`. No synthetic injection needed —
        // the FN obs whose `final_dt` is FN/N at level==0 already
        // sit at MIN_POSITIVE, so they form the "0.000000"-formatted
        // baseline row naturally after sort + roc-delta.
        //
        // When a `shared_sorted` vector is provided (typically the
        // `(*)`-sibling's already-sorted obs), filter it to obs that
        // belong to `my_subtype` and use the resulting view directly.
        // This mirrors legacy's "sort once per Type, walk with flag
        // mask" pipeline and ensures the cluster-ordering outcome at
        // any per-subtype slice matches the global sort.
        let (mut sorted_obs, presorted): (Vec<ObsRecord>, bool) = match (shared_sorted, my_subtype)
        {
            (Some(sorted), Some(st)) if st != "*" => {
                let filtered: Vec<ObsRecord> = sorted
                    .iter()
                    .filter(|o| o.has_subtype(st))
                    .cloned()
                    .collect();
                (filtered, true)
            }
            _ => (self.observations.small_records(), false),
        };

        // Match legacy's `OBS_FLAG_TI` / `OBS_FLAG_TV` masks: the flag is
        // set per-record from the BI string at insertion time, *not*
        // derived from non-zero count fields. Filter-failed query
        // phantoms carry their flag bits despite all-zero counts and
        // contribute to the obs sweep order legacy walks.
        fn obs_has_ti(o: &ObsRecord) -> bool {
            o.ti_flag
        }
        fn obs_has_tv(o: &ObsRecord) -> bool {
            o.tv_flag
        }
        // Optimization: only do per-substat snapshots for groups that have
        // ti or tv contributions (SNP groups).
        let needs_substat_snapshots =
            !presorted && (sorted_obs.iter().any(|o| obs_has_ti(o) || obs_has_tv(o)));
        // BASE row main counts: 4 sorts. After 4 sorts, sorted_obs is the
        // 4-sort state. We use it directly for the main loop, but we also
        // need to do 8/12 sorts for substats. Solution: clone the 4-sort
        // state once for substat work; main BASE uses sorted_obs directly.
        if !presorted {
            for _ in 0..4 {
                introsort_libstdcpp(&mut sorted_obs);
            }
        }
        // Build per-substat indexes. Use ONE shared clone (substat_obs) that
        // we progressively sort: 4 more sorts for target_8 (total 8 sorts),
        // 4 more for target_12 (total 12 sorts). Saves cloning overhead.
        let (target_8_index, total_8, target_12_index, total_12) = if needs_substat_snapshots {
            let mut substat_obs = sorted_obs.clone();
            // 8-sort state.
            for _ in 0..4 {
                introsort_libstdcpp(&mut substat_obs);
            }
            let mut last8 = Cumul::default();
            let mut idx_8: std::collections::HashMap<u64, Cumul> = std::collections::HashMap::new();
            let mut total_ti = Cumul::default();
            for o in substat_obs.iter().filter(|o| obs_has_ti(o)) {
                last8.add(&o.counts);
                total_ti.add(&o.counts);
                idx_8
                    .entry(o.level.to_bits())
                    .or_insert_with(|| last8.clone());
            }
            // 12-sort state.
            for _ in 0..4 {
                introsort_libstdcpp(&mut substat_obs);
            }
            let mut last12 = Cumul::default();
            let mut idx_12: std::collections::HashMap<u64, Cumul> =
                std::collections::HashMap::new();
            let mut total_tv = Cumul::default();
            for o in substat_obs.iter().filter(|o| obs_has_tv(o)) {
                last12.add(&o.counts);
                total_tv.add(&o.counts);
                idx_12
                    .entry(o.level.to_bits())
                    .or_insert_with(|| last12.clone());
            }
            (idx_8, total_ti, idx_12, total_tv)
        } else {
            (
                std::collections::HashMap::new(),
                Cumul::default(),
                std::collections::HashMap::new(),
                Cumul::default(),
            )
        };

        // Build target: cumulative AFTER each obs (inclusive) on sorted_obs
        // (which is now in 4-sort state for non-presorted, or shared-sort
        // state for presorted).
        let mut last = Cumul::default();
        let mut target_levels: Vec<f64> = Vec::with_capacity(sorted_obs.len());
        let mut target_cums: Vec<Cumul> = Vec::with_capacity(sorted_obs.len());
        for obs in &sorted_obs {
            last.add(&obs.counts);
            target_levels.push(obs.level);
            target_cums.push(last.clone());
        }
        let total = last.clone();

        // Apply roc-delta=0.5 filter: walk ASC, keep first row, then
        // rows where |level − prev_kept_level| > delta.
        let mut kept: Vec<(f64, Cumul)> = Vec::new();
        let mut prev: Option<f64> = None;
        for (i, level) in target_levels.iter().enumerate() {
            let push = match prev {
                None => true,
                Some(p) => (level - p).abs() > delta,
            };
            if push {
                prev = Some(*level);
                kept.push((*level, target_cums[i].clone()));
            }
        }

        // Build numeric ROC rows using legacy's strict-above cumulative
        // formula at the kept obs.
        let total_truth_fn = sub_buckets(&truth_total_const, &total.truth_tp);
        let mut numeric_rows: Vec<EmittedRow> = Vec::with_capacity(kept.len());
        for (level, cum_through) in &kept {
            let mut truth_tp = sub_buckets(&total.truth_tp, &cum_through.truth_tp);
            let mut truth_fn = sum_buckets(&total_truth_fn, &cum_through.truth_tp);
            let mut query_tp = sub_buckets(&total.query_tp, &cum_through.query_tp);
            let mut query_fp = sub_buckets(&total.query_fp, &cum_through.query_fp);
            let mut query_unk = sub_buckets(&total.query_unk, &cum_through.query_unk);

            // Override .ti sub-fields using snapshot_8 cum_through.
            if let Some(c8) = target_8_index.get(&level.to_bits()) {
                truth_tp.ti = total_8.truth_tp.ti.saturating_sub(c8.truth_tp.ti);
                truth_fn.ti = total_truth_fn.ti + c8.truth_tp.ti;
                query_tp.ti = total_8.query_tp.ti.saturating_sub(c8.query_tp.ti);
                query_fp.ti = total_8.query_fp.ti.saturating_sub(c8.query_fp.ti);
                query_unk.ti = total_8.query_unk.ti.saturating_sub(c8.query_unk.ti);
            }

            // Override .tv sub-fields using snapshot_12 cum_through.
            if let Some(c12) = target_12_index.get(&level.to_bits()) {
                truth_tp.tv = total_12.truth_tp.tv.saturating_sub(c12.truth_tp.tv);
                truth_fn.tv = total_truth_fn.tv + c12.truth_tp.tv;
                query_tp.tv = total_12.query_tp.tv.saturating_sub(c12.query_tp.tv);
                query_fp.tv = total_12.query_fp.tv.saturating_sub(c12.query_fp.tv);
                query_unk.tv = total_12.query_unk.tv.saturating_sub(c12.query_unk.tv);
            }

            let fp_gt = total.fp_gt.saturating_sub(cum_through.fp_gt);
            let fp_al = total.fp_al.saturating_sub(cum_through.fp_al);
            let cum = Cumul {
                truth_tp,
                truth_fn,
                query_tp,
                query_fp,
                query_unk,
                fp_gt,
                fp_al,
            };
            // Format level via f32 → f64 → "%.6f" to match legacy's
            // bucket-key cardinality (e.g. `1001.340027` rather than
            // `1001.340000`).
            let level_f32 = *level as f32 as f64;
            let qq_str = rendered_roc_level(level_f32);
            numeric_rows.push(EmittedRow {
                qq_str,
                cum,
                substats: None,
            });
        }

        // ----- Per-substat sweeps (ti/tv only — het/homalt suppressed) -----
        // Legacy's per-substat sweep filters obs by flag mask and applies
        // its own roc-delta. We mirror this with the bucket map (already
        // collapsed to {:.6}-precision keys), since the per-substat sweep
        // doesn't have the truth/query-cluster-order ambiguity that
        // affects the main count columns.
        let zero_key = format!("{:.6}", 0.0_f64);
        let synthetic_zero = NumericBucket {
            qq: 0.0,
            counts: Cumul::default(),
            has_ti: false,
            has_tv: false,
        };
        let zero_already_bucket = self.threshold_window.contains_key(&zero_key);
        let need_synthetic_zero_bucket = truth_total_const.total > 0 && !zero_already_bucket;
        let mut ordered: Vec<(&String, &NumericBucket)> = self.threshold_window.iter().collect();
        if need_synthetic_zero_bucket {
            ordered.push((&zero_key, &synthetic_zero));
        }
        ordered.sort_by(|(_, a), (_, b)| b.qq.partial_cmp(&a.qq).unwrap_or(Ordering::Equal));
        let ti_kept = substat_kept_levels(&ordered, delta, bucket_has_ti);
        let tv_kept = substat_kept_levels(&ordered, delta, bucket_has_tv);

        // Attach substat availability to each numeric row by keying on
        // the row's QQ string against the `_kept` sets.
        for row in &mut numeric_rows {
            row.substats = Some(SubstatAvail {
                ti: ti_kept.contains(&row.qq_str),
                tv: tv_kept.contains(&row.qq_str),
                het: false,
                homalt: false,
            });
        }

        // Re-sort by the lex-ASC QQ string so output ordering matches
        // pandas' object-column sort (legacy roc.all.csv.gz convention).
        numeric_rows.sort_by(|a, b| a.qq_str.cmp(&b.qq_str));
        out.extend(numeric_rows);
        out
    }

    fn emit_disk_with_delta_from(
        &self,
        source: &ObservationStore,
        subtype: Option<&str>,
        sort_passes: usize,
        delta: f64,
    ) -> Result<EmittedRows> {
        let truth_total_const = sum_buckets(&self.baseline.truth_tp, &self.baseline.truth_fn);
        let mut total = Cumul::default();
        let mut total_ti = Cumul::default();
        let mut total_tv = Cumul::default();
        for observation in source.sorted_with_passes(sort_passes)? {
            let observation = observation?;
            if subtype.is_some_and(|subtype| !observation.has_subtype(subtype)) {
                continue;
            }
            total.add(&observation.counts);
            if observation.ti_flag {
                total_ti.add(&observation.counts);
            }
            if observation.tv_flag {
                total_tv.add(&observation.counts);
            }
        }

        let total_truth_fn = sub_buckets(&truth_total_const, &total.truth_tp);
        let mut numeric_rows = EmittedRowSpool::new();
        let mut running = Cumul::default();
        let mut ti_snapshots = SubstatSnapshotCursor::build(source, 8, |row| {
            row.ti_flag && subtype.is_none_or(|subtype| row.has_subtype(subtype))
        })?;
        let mut tv_snapshots = SubstatSnapshotCursor::build(source, 12, |row| {
            row.tv_flag && subtype.is_none_or(|subtype| row.has_subtype(subtype))
        })?;
        let mut previous_main: Option<f64> = None;
        let mut previous_ti: Option<f64> = None;
        let mut previous_tv: Option<f64> = None;
        let mut current: Option<(f64, Cumul)> = None;

        let mut emit_level = |level: f64, cum_through: Cumul| -> Result<()> {
            let ti_through = ti_snapshots.counts_at(level)?;
            let tv_through = tv_snapshots.counts_at(level)?;
            let has_ti = ti_through.is_some();
            let has_tv = tv_through.is_some();
            let keep_main = previous_main.is_none_or(|value| (level - value).abs() > delta);
            let keep_ti = has_ti && previous_ti.is_none_or(|value| (level - value).abs() > delta);
            let keep_tv = has_tv && previous_tv.is_none_or(|value| (level - value).abs() > delta);
            if keep_ti {
                previous_ti = Some(level);
            }
            if keep_tv {
                previous_tv = Some(level);
            }
            if !keep_main {
                return Ok(());
            }
            previous_main = Some(level);
            let mut truth_tp = sub_buckets(&total.truth_tp, &cum_through.truth_tp);
            let mut truth_fn = sum_buckets(&total_truth_fn, &cum_through.truth_tp);
            let mut query_tp = sub_buckets(&total.query_tp, &cum_through.query_tp);
            let mut query_fp = sub_buckets(&total.query_fp, &cum_through.query_fp);
            let mut query_unk = sub_buckets(&total.query_unk, &cum_through.query_unk);
            if let Some(cumulative) = ti_through {
                truth_tp.ti = total_ti.truth_tp.ti.saturating_sub(cumulative.truth_tp.ti);
                truth_fn.ti = total_truth_fn.ti + cumulative.truth_tp.ti;
                query_tp.ti = total_ti.query_tp.ti.saturating_sub(cumulative.query_tp.ti);
                query_fp.ti = total_ti.query_fp.ti.saturating_sub(cumulative.query_fp.ti);
                query_unk.ti = total_ti
                    .query_unk
                    .ti
                    .saturating_sub(cumulative.query_unk.ti);
            }
            if let Some(cumulative) = tv_through {
                truth_tp.tv = total_tv.truth_tp.tv.saturating_sub(cumulative.truth_tp.tv);
                truth_fn.tv = total_truth_fn.tv + cumulative.truth_tp.tv;
                query_tp.tv = total_tv.query_tp.tv.saturating_sub(cumulative.query_tp.tv);
                query_fp.tv = total_tv.query_fp.tv.saturating_sub(cumulative.query_fp.tv);
                query_unk.tv = total_tv
                    .query_unk
                    .tv
                    .saturating_sub(cumulative.query_unk.tv);
            }
            numeric_rows.push(EmittedRow {
                qq_str: rendered_roc_level(level),
                cum: Cumul {
                    truth_tp,
                    truth_fn,
                    query_tp,
                    query_fp,
                    query_unk,
                    fp_gt: total.fp_gt.saturating_sub(cum_through.fp_gt),
                    fp_al: total.fp_al.saturating_sub(cum_through.fp_al),
                },
                substats: Some(SubstatAvail {
                    ti: keep_ti,
                    tv: keep_tv,
                    het: false,
                    homalt: false,
                }),
            })?;
            Ok(())
        };

        for observation in source.sorted_with_passes(sort_passes)? {
            let observation = observation?;
            if subtype.is_some_and(|subtype| !observation.has_subtype(subtype)) {
                continue;
            }
            if current
                .as_ref()
                .is_some_and(|entry| entry.0 != observation.level)
                && let Some(entry) = current.take()
            {
                emit_level(entry.0, entry.1)?;
            }
            running.add(&observation.counts);
            current.get_or_insert_with(|| (observation.level, running.clone()));
        }
        if let Some(entry) = current {
            emit_level(entry.0, entry.1)?;
        }
        Ok(EmittedRows::Disk(Box::new(DiskEmittedRows {
            baseline: Some(EmittedRow {
                qq_str: "*".to_string(),
                cum: self.baseline.clone(),
                substats: None,
            }),
            merge: numeric_rows.finish()?,
        })))
    }
}

fn substat_kept_levels(
    ordered_desc: &[(&String, &NumericBucket)],
    delta: f64,
    has_substat: impl Fn(&NumericBucket) -> bool,
) -> HashSet<String> {
    let mut candidates: Vec<(&String, f64)> = ordered_desc
        .iter()
        .filter(|(_, b)| has_substat(b))
        .map(|(k, b)| (*k, b.qq))
        .collect();
    candidates.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
    let mut kept = HashSet::new();
    let mut prev: Option<f64> = None;
    for (key, q) in candidates {
        let push = match prev {
            None => true,
            Some(p) => (q - p).abs() > delta,
        };
        if push {
            prev = Some(q);
            kept.insert(key.clone());
        }
    }
    kept
}

fn bucket_has_ti(b: &NumericBucket) -> bool {
    b.has_ti
}

fn bucket_has_tv(b: &NumericBucket) -> bool {
    b.has_tv
}

fn sub_buckets(a: &CountsBucket, b: &CountsBucket) -> CountsBucket {
    // Saturating subtraction: clamps to 0 if somehow cum ever exceeds the
    // baseline total (can't happen if accumulation is consistent, but this
    // keeps usize-underflow at bay during incremental development).
    CountsBucket {
        total: a.total.saturating_sub(b.total),
        ti: a.ti.saturating_sub(b.ti),
        tv: a.tv.saturating_sub(b.tv),
        het: a.het.saturating_sub(b.het),
        homalt: a.homalt.saturating_sub(b.homalt),
    }
}

// ---------------------------------------------------------------------------
// libstdc++ std::sort emulation (introsort)
// ---------------------------------------------------------------------------
//
// Legacy `Roc::getLevels` (src/c++/lib/tools/Roc.cpp in the upstream
// hap.py source preserved at commit 9a3ac91) calls `std::sort` on the
// per-ROC observation vector with a `level < level` comparator. To
// reproduce legacy's ROC cumulative cells bit-exactly we mirror
// libstdc++-v3's std::sort algorithm: introsort (median-of-three
// quicksort) with a 16-element threshold switch to insertion sort,
// and a depth-limit heapsort fallback.
//
// Validated against gcc-15 libstdc++ on synthetic equal-level-cluster
// inputs: 8333-element shuffled-level dataset with 33% truth-only
// records and 67% truth+query pairs sorts byte-identically.
//
// Reference: libstdc++-v3/include/bits/stl_algo.h, gcc-15.

const ROC_SORT_THRESHOLD: usize = 16;

#[inline]
fn lg_floor(n: usize) -> usize {
    if n <= 1 {
        0
    } else {
        (usize::BITS - 1 - n.leading_zeros()) as usize
    }
}

fn introsort_libstdcpp(arr: &mut [ObsRecord]) {
    introsort_libstdcpp_by(arr, |left, right| left.level < right.level);
}

fn introsort_libstdcpp_by<T, F>(arr: &mut [T], less: F)
where
    F: Fn(&T, &T) -> bool + Copy,
{
    let n = arr.len();
    if n > 1 {
        let depth_limit = lg_floor(n) * 2;
        introsort_loop(arr, 0, n, depth_limit, less);
        final_insertion_sort(arr, less);
    }
}

fn introsort_loop<T, F>(
    arr: &mut [T],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
    less: F,
) where
    F: Fn(&T, &T) -> bool + Copy,
{
    while last - first > ROC_SORT_THRESHOLD {
        if depth_limit == 0 {
            heapsort_range(arr, first, last, less);
            return;
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(arr, first, last, less);
        introsort_loop(arr, cut, last, depth_limit, less);
        last = cut;
    }
}

fn unguarded_partition_pivot<T, F>(arr: &mut [T], first: usize, last: usize, less: F) -> usize
where
    F: Fn(&T, &T) -> bool + Copy,
{
    let mid = first + (last - first) / 2;
    move_median_to_first(arr, first, first + 1, mid, last - 1, less);
    unguarded_partition(arr, first + 1, last, first, less)
}

fn move_median_to_first<T, F>(arr: &mut [T], result: usize, a: usize, b: usize, c: usize, less: F)
where
    F: Fn(&T, &T) -> bool,
{
    if less(&arr[a], &arr[b]) {
        if less(&arr[b], &arr[c]) {
            arr.swap(result, b);
        } else if less(&arr[a], &arr[c]) {
            arr.swap(result, c);
        } else {
            arr.swap(result, a);
        }
    } else if less(&arr[a], &arr[c]) {
        arr.swap(result, a);
    } else if less(&arr[b], &arr[c]) {
        arr.swap(result, c);
    } else {
        arr.swap(result, b);
    }
}

fn unguarded_partition<T, F>(
    arr: &mut [T],
    mut first: usize,
    mut last: usize,
    pivot: usize,
    less: F,
) -> usize
where
    F: Fn(&T, &T) -> bool,
{
    loop {
        while less(&arr[first], &arr[pivot]) {
            first += 1;
        }
        last -= 1;
        while less(&arr[pivot], &arr[last]) {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        arr.swap(first, last);
        first += 1;
    }
}

fn final_insertion_sort<T, F>(arr: &mut [T], less: F)
where
    F: Fn(&T, &T) -> bool + Copy,
{
    let n = arr.len();
    if n > ROC_SORT_THRESHOLD {
        insertion_sort_range(arr, 0, ROC_SORT_THRESHOLD, less);
        unguarded_insertion_sort_range(arr, ROC_SORT_THRESHOLD, n, less);
    } else {
        insertion_sort_range(arr, 0, n, less);
    }
}

fn insertion_sort_range<T, F>(arr: &mut [T], first: usize, last: usize, less: F)
where
    F: Fn(&T, &T) -> bool,
{
    if first == last {
        return;
    }
    for i in (first + 1)..last {
        if less(&arr[i], &arr[first]) {
            arr[first..=i].rotate_right(1);
        } else {
            let mut j = i;
            while j > first && less(&arr[j], &arr[j - 1]) {
                arr.swap(j, j - 1);
                j -= 1;
            }
        }
    }
}

fn unguarded_insertion_sort_range<T, F>(arr: &mut [T], first: usize, last: usize, less: F)
where
    F: Fn(&T, &T) -> bool,
{
    for i in first..last {
        let mut j = i;
        while j > 0 && less(&arr[j], &arr[j - 1]) {
            arr.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn heapsort_range<T, F>(arr: &mut [T], first: usize, last: usize, less: F)
where
    F: Fn(&T, &T) -> bool + Copy,
{
    let n = last - first;
    if n < 2 {
        return;
    }
    for i in (0..n / 2).rev() {
        sift_down(arr, first + i, first, last, less);
    }
    for i in (1..n).rev() {
        arr.swap(first, first + i);
        sift_down(arr, first, first, first + i, less);
    }
}

fn sift_down<T, F>(arr: &mut [T], start: usize, first: usize, last: usize, less: F)
where
    F: Fn(&T, &T) -> bool,
{
    let mut root = start;
    loop {
        let lc = first + 2 * (root - first) + 1;
        if lc >= last {
            break;
        }
        let rc = lc + 1;
        let mut child = lc;
        if rc < last && less(&arr[lc], &arr[rc]) {
            child = rc;
        }
        if less(&arr[root], &arr[child]) {
            arr.swap(root, child);
            root = child;
        } else {
            break;
        }
    }
}

fn disk_introsort_libstdcpp(file: &mut File, len: usize) -> Result<()> {
    if len > 1 {
        disk_introsort_loop(file, 0, len, lg_floor(len) * 2)?;
        disk_final_insertion_sort(file, len)?;
    }
    Ok(())
}

fn disk_introsort_loop(
    file: &mut File,
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
) -> Result<()> {
    while last - first > ROC_SORT_THRESHOLD {
        if depth_limit == 0 {
            return disk_heapsort_range(file, first, last);
        }
        depth_limit -= 1;
        let cut = disk_unguarded_partition_pivot(file, first, last)?;
        disk_introsort_loop(file, cut, last, depth_limit)?;
        last = cut;
    }
    Ok(())
}

fn disk_unguarded_partition_pivot(file: &mut File, first: usize, last: usize) -> Result<usize> {
    let mid = first + (last - first) / 2;
    disk_move_median_to_first(file, first, first + 1, mid, last - 1)?;
    disk_unguarded_partition(file, first + 1, last, first)
}

fn disk_move_median_to_first(
    file: &mut File,
    result: usize,
    a: usize,
    b: usize,
    c: usize,
) -> Result<()> {
    let av = disk_level(file, a)?;
    let bv = disk_level(file, b)?;
    let cv = disk_level(file, c)?;
    let selected = if av < bv {
        if bv < cv {
            b
        } else if av < cv {
            c
        } else {
            a
        }
    } else if av < cv {
        a
    } else if bv < cv {
        c
    } else {
        b
    };
    swap_disk_index(file, result, selected)
}

fn disk_unguarded_partition(
    file: &mut File,
    mut first: usize,
    mut last: usize,
    pivot: usize,
) -> Result<usize> {
    loop {
        while disk_level(file, first)? < disk_level(file, pivot)? {
            first += 1;
        }
        last -= 1;
        while disk_level(file, pivot)? < disk_level(file, last)? {
            last -= 1;
        }
        if first >= last {
            return Ok(first);
        }
        swap_disk_index(file, first, last)?;
        first += 1;
    }
}

fn disk_final_insertion_sort(file: &mut File, len: usize) -> Result<()> {
    if len > ROC_SORT_THRESHOLD {
        disk_insertion_sort_range(file, 0, ROC_SORT_THRESHOLD)?;
        disk_unguarded_insertion_sort_range(file, ROC_SORT_THRESHOLD, len)
    } else {
        disk_insertion_sort_range(file, 0, len)
    }
}

fn disk_insertion_sort_range(file: &mut File, first: usize, last: usize) -> Result<()> {
    if first == last {
        return Ok(());
    }
    for i in (first + 1)..last {
        if disk_level(file, i)? < disk_level(file, first)? {
            let entry = read_disk_index_entry(file, i)?;
            for position in (first..i).rev() {
                let shifted = read_disk_index_entry(file, position)?;
                write_disk_index_at(file, position + 1, shifted)?;
            }
            write_disk_index_at(file, first, entry)?;
        } else {
            let mut position = i;
            while position > first && disk_level(file, position)? < disk_level(file, position - 1)?
            {
                swap_disk_index(file, position, position - 1)?;
                position -= 1;
            }
        }
    }
    Ok(())
}

fn disk_unguarded_insertion_sort_range(file: &mut File, first: usize, last: usize) -> Result<()> {
    for i in first..last {
        let mut position = i;
        while position > 0 && disk_level(file, position)? < disk_level(file, position - 1)? {
            swap_disk_index(file, position, position - 1)?;
            position -= 1;
        }
    }
    Ok(())
}

fn disk_heapsort_range(file: &mut File, first: usize, last: usize) -> Result<()> {
    let len = last - first;
    if len < 2 {
        return Ok(());
    }
    for index in (0..len / 2).rev() {
        disk_sift_down(file, first + index, first, last)?;
    }
    for index in (1..len).rev() {
        swap_disk_index(file, first, first + index)?;
        disk_sift_down(file, first, first, first + index)?;
    }
    Ok(())
}

fn disk_sift_down(file: &mut File, start: usize, first: usize, last: usize) -> Result<()> {
    let mut root = start;
    loop {
        let left = first + 2 * (root - first) + 1;
        if left >= last {
            break;
        }
        let right = left + 1;
        let child = if right < last && disk_level(file, left)? < disk_level(file, right)? {
            right
        } else {
            left
        };
        if disk_level(file, root)? < disk_level(file, child)? {
            swap_disk_index(file, root, child)?;
            root = child;
        } else {
            break;
        }
    }
    Ok(())
}

#[cfg(test)]
fn accumulate(rows: &[AnnotatedRow]) -> BTreeMap<RowKey, BoundedGroupAccum> {
    accumulate_impl(
        rows.iter().map(Ok::<_, anyhow::Error>),
        &RocOptions::default(),
    )
    .expect("in-memory ROC rows are infallible")
}

fn accumulate_impl<I, R>(
    rows: I,
    options: &RocOptions,
) -> Result<BTreeMap<RowKey, BoundedGroupAccum>>
where
    I: IntoIterator<Item = Result<R>>,
    R: Borrow<AnnotatedRow>,
{
    let mut groups: BTreeMap<RowKey, BoundedGroupAccum> = BTreeMap::new();
    let track_thresholds = options.output_rocs || options.preserve_raw_table;

    // Named stratifications are configured lanes, not merely observed axes.
    // QuantifyRegions registers each one when it loads the BED, so an empty
    // fourth-column child still receives zero-valued baseline ROC rows for
    // every active variant type. Built-in TS_* lanes remain observation-led.
    let mut observed_subsets = options
        .subset_sizes
        .keys()
        .cloned()
        .collect::<BTreeSet<_>>();
    // Track which non-PASS filter tags appeared per variant Type so we
    // pre-seed empty subtype rows only for filters that variant type
    // actually carries (e.g. SB filter is SNP-only in our chr21 data —
    // we shouldn't seed INDEL SB rows that legacy doesn't emit).
    let mut observed_filters_per_ty: std::collections::BTreeMap<
        String,
        std::collections::BTreeSet<String>,
    > = std::collections::BTreeMap::new();
    observed_subsets.insert("*".to_string());
    for row in rows {
        let row = row?;
        let row = row.borrow();
        emit_contributions_with_options(
            row,
            options,
            |key, qq, counts, subtypes: &[String], bi: Option<&str>, blt: Option<&str>| {
                observed_subsets.insert(key.subset.clone());
                if key.filter != "ALL" && key.filter != "PASS" {
                    observed_filters_per_ty
                        .entry(key.ty.clone())
                        .or_default()
                        .insert(key.filter.clone());
                }
                groups
                    .entry(key)
                    .or_default()
                    .add(qq, counts, subtypes, bi, blt, track_thresholds);
            },
        );
    }

    // Pre-seed baseline rows for every (type, subtype, subset, filter)
    // expected by legacy's RocOutput. Legacy's `_addLevel` initialises
    // the output table with a baseline entry per (subtype × genotype)
    // cartesian even when the corresponding ROC has zero observations —
    // that's why legacy's roc.all contains C1_5/C6_15/C16_PLUS baseline
    // rows (all zero) while rust previously only emitted subtypes with
    // real records. Using `or_default` preserves any accumulated data.
    for (ty, subtypes) in [
        ("SNP", &["*"][..]),
        (
            "INDEL",
            &[
                "*", "C1_5", "C6_15", "C16_PLUS", "D1_5", "D6_15", "D16_PLUS", "I1_5", "I6_15",
                "I16_PLUS",
            ][..],
        ),
    ] {
        for subtype in subtypes {
            // Per-Filter pre-seed: legacy emits one row per (Type,
            // Subtype, Subset='*', Filter=<observed-tag>) — including
            // empty-bucket subtypes like C16_PLUS — at QQ=*, but only
            // for filter tags that this Type actually carries (SB on
            // chr21 is SNP-only; INDEL must not carry SB rows).
            if let Some(filters) = observed_filters_per_ty.get(ty) {
                for tag in filters {
                    groups
                        .entry(RowKey::new_with_qq_field(
                            ty,
                            subtype,
                            "*",
                            tag,
                            &options.qq_field,
                        ))
                        .or_default();
                }
            }
            for subset in &observed_subsets {
                let mut filters = vec!["ALL", "PASS"];
                if !options.ignored_filters.is_empty() {
                    filters.push("SEL");
                }
                for filter in filters {
                    groups
                        .entry(RowKey::new_with_qq_field(
                            ty,
                            subtype,
                            subset,
                            filter,
                            &options.qq_field,
                        ))
                        .or_default();
                }
            }
        }
    }

    for group in groups.values_mut() {
        group.observations.finish()?;
    }
    Ok(groups)
}

#[cfg(test)]
fn accumulate_with_options(
    rows: &[AnnotatedRow],
    options: &RocOptions,
) -> BTreeMap<RowKey, BoundedGroupAccum> {
    accumulate_impl(rows.iter().map(Ok::<_, anyhow::Error>), options)
        .expect("in-memory ROC rows are infallible")
}

fn roc_header(ci_alpha: f64) -> String {
    let mut header = EXTENDED_HEADER
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>();
    if ci_alpha > 0.0 {
        header.extend(
            [
                "METRIC.Recall.Lower",
                "METRIC.Recall.Upper",
                "METRIC.Precision.Lower",
                "METRIC.Precision.Upper",
                "METRIC.Frac_NA.Lower",
                "METRIC.Frac_NA.Upper",
            ]
            .into_iter()
            .map(str::to_string),
        );
    }
    header.join(",")
}

fn is_aggregate_filter(filter: &str) -> bool {
    matches!(filter, "ALL" | "PASS" | "SEL")
}

// ---------------------------------------------------------------------------
// Per-row contribution emission
// ---------------------------------------------------------------------------

/// Lightweight VCF sample parser (reads the handful of FORMAT keys roc.rs
/// cares about). Mirrors `compare::SampleView` without introducing a
/// cross-module dep on the private struct.
struct Sample<'a> {
    format_keys: &'a [&'a str],
    parts: &'a [&'a str],
    gt: Option<&'a str>,
    bd: Option<&'a str>,
    bi: Option<&'a str>,
    bvt: Option<&'a str>,
    blt: Option<&'a str>,
    qq: Option<&'a str>,
}

impl<'a> Sample<'a> {
    fn new(format_keys: &'a [&'a str], parts: &'a [&'a str]) -> Self {
        let lookup = |name: &str| -> Option<&'a str> {
            format_keys
                .iter()
                .position(|key| *key == name)
                .and_then(|index| parts.get(index).copied())
        };
        Self {
            format_keys,
            parts,
            gt: lookup("GT"),
            bd: lookup("BD"),
            bi: lookup("BI"),
            bvt: lookup("BVT"),
            blt: lookup("BLT"),
            qq: lookup("QQ"),
        }
    }

    fn variant_type(&self) -> Option<&'a str> {
        self.bvt.filter(|value| matches!(*value, "SNP" | "INDEL"))
    }

    /// Per-side ROC threshold value: legacy hap.py uses each sample's
    /// FORMAT.QQ for its own ROC sweep, **not** the record's QUAL column.
    /// This matters on multi-allelic indels where xcmp emits `QUAL=0`
    /// at the record level but writes the matched query's quality into
    /// each sample's QQ field — meaning truth-side and query-side ROC
    /// thresholds can differ for the same record. Reading QUAL would
    /// drop ~800 chr21 records into the QQ=0 bucket and shift their
    /// cumulative TPs out of the > QUAL>0 thresholds, breaking ROC
    /// parity at the upper end.
    fn roc_qq(&self) -> Option<f64> {
        self.qq.and_then(|raw| raw.parse::<f64>().ok())
    }

    fn roc_value(&self, field: &str, record_qual: &str, info: &str) -> Option<f64> {
        // QUAL has already been propagated into FORMAT/QQ by xcmp, including
        // the truth-side superlocus minimum adjustment. Reading QQ here is
        // therefore required for parity rather than using the record column.
        if field == "QUAL" || field == "QQ" {
            return self.roc_qq();
        }
        if field == "." {
            return None;
        }
        if let Some(value) = info.split(';').find_map(|entry| {
            entry
                .split_once('=')
                .filter(|(key, _)| *key == field)
                .map(|(_, value)| value)
        }) {
            return parse_roc_number(value);
        }
        self.format_keys
            .iter()
            .position(|key| *key == field)
            .and_then(|index| self.parts.get(index).copied())
            .and_then(parse_roc_number)
            .or_else(|| {
                (field == "QUAL")
                    .then(|| parse_roc_number(record_qual))
                    .flatten()
            })
    }
}

fn parse_roc_number(raw: &str) -> Option<f64> {
    raw.split(',')
        .next()
        .filter(|value| !matches!(*value, "" | "."))
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite())
}

fn emit_contributions_with_options<
    F: FnMut(RowKey, Option<f64>, &Cumul, &[String], Option<&str>, Option<&str>),
>(
    row: &AnnotatedRow,
    options: &RocOptions,
    mut emit: F,
) {
    let line = row.record.raw().to_line();
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return;
    }

    let info = fields[7];
    let subsets = extract_subsets(info);

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let truth_parts: Vec<&str> = fields[9].split(':').collect();
    let query_parts: Vec<&str> = fields[10].split(':').collect();
    let truth = Sample::new(&format_keys, &truth_parts);
    let query = Sample::new(&format_keys, &query_parts);

    // Per-side ROC threshold: read FORMAT.QQ from each sample. The
    // record-level QUAL column would lose multi-allelic indels whose
    // xcmp output sets QUAL=0 while still writing the matched per-side
    // quality into FORMAT.QQ.
    let score_field = options.score_field.as_deref().unwrap_or(&options.qq_field);
    let truth_qq = truth.roc_value(score_field, fields[5], info);
    let query_qq = query.roc_value(score_field, fields[5], info);

    let filter_tags = fields[6]
        .split(';')
        .filter(|tag| !tag.is_empty() && *tag != "." && *tag != "PASS")
        .collect::<Vec<_>>();
    let ignores =
        |tag: &str| options.ignored_filters.contains("*") || options.ignored_filters.contains(tag);
    let selectively_filtered = filter_tags.iter().any(|tag| !ignores(tag));

    // Selective filtering adds the legacy SEL tier: filters named by
    // --roc-filter are ignored for SEL, while PASS continues to apply every
    // record filter and ALL applies none.
    let mut tiers = vec![("ALL", false), ("PASS", !row.query_pass)];
    if !options.ignored_filters.is_empty() {
        tiers.push(("SEL", selectively_filtered));
    }
    for (filter, filtered_out) in tiers {
        // Truth side
        if let Some(truth_bvt) = truth.variant_type() {
            let effective = if filtered_out && truth.bd == Some("TP") {
                Some("FN")
            } else {
                truth.bd
            };
            if matches!(effective, Some("TP") | Some("FN")) {
                let bucket = sample_bucket(&truth);
                let mut counts = Cumul::default();
                if effective == Some("TP") {
                    counts.truth_tp = bucket;
                } else {
                    counts.truth_fn = bucket;
                }
                emit_for_axes(
                    AxisContribution {
                        ty: truth_bvt,
                        bi: truth.bi,
                        blt: truth.blt,
                        subsets: &subsets,
                        filter,
                        qq_field: &options.qq_field,
                        qq: truth_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
        }

        // Query side. When the row passes the filter, contribute the
        // BD-typed bucket. When it fails the filter (PASS tier only),
        // emit a zero-counts phantom at every (subtype, subset) the
        // record's BI / region tags imply — legacy's BlockQuantify
        // demotes filter-failed obs but still adds them to the ROC
        // vector with the obs's flag bits, so the per-subtype level
        // bucket exists at the demoted level even though its TTP/QTP
        // counts come entirely from cumulative-above passing records.
        if let Some(query_bvt) = query.variant_type() {
            if !filtered_out && matches!(query.bd, Some("TP") | Some("FP") | Some("UNK")) {
                let bucket = sample_bucket(&query);
                let mut counts = Cumul::default();
                match query.bd {
                    Some("TP") => counts.query_tp = bucket,
                    Some("FP") => {
                        counts.query_fp = bucket;
                        if row.fp_class == Some(FpClass::Gt) {
                            counts.fp_gt = 1;
                        } else if row.fp_class == Some(FpClass::Al) {
                            counts.fp_al = 1;
                        }
                    }
                    Some("UNK") => counts.query_unk = bucket,
                    _ => {}
                }
                emit_for_axes(
                    AxisContribution {
                        ty: query_bvt,
                        bi: query.bi,
                        blt: query.blt,
                        subsets: &subsets,
                        filter,
                        qq_field: &options.qq_field,
                        qq: query_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            } else if filtered_out && matches!(query.bd, Some("TP") | Some("FP") | Some("UNK")) {
                // Filter-failed query phantom (no count contribution).
                // Pass bi through phantom_bi for substat detection.
                let counts = Cumul::default();
                let phantom_bi = if query.bd == Some("TP") {
                    query.bi
                } else {
                    None
                };
                emit_for_axes(
                    AxisContribution {
                        ty: query_bvt,
                        bi: phantom_bi,
                        blt: (query.bd == Some("TP")).then_some(query.blt).flatten(),
                        subsets: &subsets,
                        filter,
                        qq_field: &options.qq_field,
                        qq: query_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
        }
    }

    // Per-Filter rows: legacy emits one (Type, Subtype, Subset='*',
    // Filter=<tag>, QQ='*') row per non-PASS filter tag. Each query
    // record with FILTER=<tag> contributes its BD bucket; truth-side
    // TPs whose paired query carries that filter contribute too. The
    // rendering for these rows zeros TRUTH.TOTAL/FN, QUERY.TOTAL and
    // METRIC.* (handled in render_row); only the *.TP / *.FP / *.UNK
    // buckets carry real counts. The filter column is fields[6]; an
    // empty/`.`/`PASS` filter contributes nothing here (already covered
    // by ALL/PASS tiers above).
    let star_subset = ["*".to_string()];
    if !filter_tags.is_empty() {
        for tag in filter_tags {
            let output_tag = if ignores(tag) {
                format!("SEL_IGN_{tag}")
            } else {
                tag.to_string()
            };
            // Truth side: TP only (FN truths have no query → no filter).
            if let Some(truth_bvt) = truth.variant_type()
                && truth.bd == Some("TP")
            {
                let bucket = sample_bucket(&truth);
                let counts = Cumul {
                    truth_tp: bucket,
                    ..Cumul::default()
                };
                emit_for_axes(
                    AxisContribution {
                        ty: truth_bvt,
                        bi: truth.bi,
                        blt: truth.blt,
                        subsets: &star_subset,
                        filter: &output_tag,
                        qq_field: &options.qq_field,
                        qq: truth_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
            // Query side: all BD types contribute.
            if let Some(query_bvt) = query.variant_type()
                && matches!(query.bd, Some("TP") | Some("FP") | Some("UNK"))
            {
                let bucket = sample_bucket(&query);
                let mut counts = Cumul::default();
                match query.bd {
                    Some("TP") => counts.query_tp = bucket,
                    Some("FP") => {
                        counts.query_fp = bucket;
                        if row.fp_class == Some(FpClass::Gt) {
                            counts.fp_gt = 1;
                        } else if row.fp_class == Some(FpClass::Al) {
                            counts.fp_al = 1;
                        }
                    }
                    Some("UNK") => counts.query_unk = bucket,
                    _ => {}
                }
                emit_for_axes(
                    AxisContribution {
                        ty: query_bvt,
                        bi: query.bi,
                        blt: query.blt,
                        subsets: &star_subset,
                        filter: &output_tag,
                        qq_field: &options.qq_field,
                        qq: query_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
        }
    }
}

struct AxisContribution<'a> {
    ty: &'a str,
    bi: Option<&'a str>,
    blt: Option<&'a str>,
    subsets: &'a [String],
    filter: &'a str,
    qq_field: &'a str,
    qq: Option<f64>,
    counts: &'a Cumul,
}

fn emit_for_axes<F: FnMut(RowKey, Option<f64>, &Cumul, &[String], Option<&str>, Option<&str>)>(
    axes: AxisContribution<'_>,
    emit: &mut F,
) {
    let subtypes = compute_subtypes(axes.ty, axes.bi);
    for subtype in &subtypes {
        for subset in axes.subsets {
            let key =
                RowKey::new_with_qq_field(axes.ty, subtype, subset, axes.filter, axes.qq_field);
            emit(key, axes.qq, axes.counts, &subtypes, axes.bi, axes.blt);
        }
    }
}

fn compute_subtypes(variant_type: &str, bi: Option<&str>) -> Vec<String> {
    let mut out = vec!["*".to_string()];
    if variant_type == "INDEL"
        && let Some(bi) = bi
    {
        // Multi-allelic INDELs emit comma-joined BI (e.g. `i1_5,i6_15`).
        // Each indel-class primitive contributes to its own subtype bucket.
        // The `ti`/`tv` decoration on complex INDELs is SNP-side and does
        // not spawn an extra bucket here — same fanout rule extended.csv
        // uses (compare::SampleView::variant_type_and_subtypes).
        for tok in bi.split(',') {
            if matches!(tok, "ti" | "tv") {
                continue;
            }
            let upper = tok.to_uppercase();
            if INDEL_SUBTYPES.contains(&upper.as_str()) && !out.contains(&upper) {
                out.push(upper);
            }
        }
    }
    out
}

fn sample_bucket(sample: &Sample<'_>) -> CountsBucket {
    let mut bucket = CountsBucket {
        total: 1,
        ..CountsBucket::default()
    };
    if let Some("SNP") = sample.bvt
        && let Some(bi) = sample.bi
    {
        // Legacy `roc::makeObservationFlags` splits BI on `,` and OR's
        // the matching flags. A multi-allelic SNP with bi="ti,tv"
        // therefore counts as both a transition AND a transversion in
        // the per-substat sweeps. Mirror that by setting both bits.
        for tok in bi.split(',') {
            match tok {
                "ti" => bucket.ti = 1,
                "tv" => bucket.tv = 1,
                _ => {}
            }
        }
    }
    // Mirror compare::add_sample_stats's general het/homalt rule so ROC
    // counts agree with extended.csv: het = exactly one allele equals the
    // reference index 0 (covers 0|2, 0/3, …); homalt = both alleles equal
    // and non-zero (1/1, 2|2, 3/3, …). Hetalt (1|2, 2/1, …) lands in
    // neither bucket. Literal "0/1" / "1/1" matching under-counted truth-
    // side multi-allelic GTs that bcftools merge preserves verbatim.
    if let Some(gt) = sample.gt {
        let alleles: Vec<Option<usize>> = gt
            .split(['/', '|'])
            .map(|part| part.parse::<usize>().ok())
            .collect();
        if let [Some(left), Some(right)] = alleles.as_slice() {
            let zero_count = usize::from(*left == 0) + usize::from(*right == 0);
            if zero_count == 1 {
                bucket.het = 1;
            } else if *left != 0 && left == right {
                bucket.homalt = 1;
            }
        }
    }
    bucket
}

fn extract_subsets(info: &str) -> Vec<String> {
    let mut out = vec!["*".to_string()];
    if let Some(tail) = info.split(";Regions=").nth(1) {
        // The Regions= value runs to end-of-INFO (we're already past the
        // FORMAT column, so there's no trailing field; but be defensive
        // against future INFO tags after Regions= by stopping at ';').
        let tag_str = tail.split(';').next().unwrap_or("");
        for tag in tag_str.split(',') {
            if !tag.is_empty() && tag != "CONF" && !out.iter().any(|value| value == tag) {
                out.push(tag.to_string());
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Row rendering
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
enum RowFilter<'a> {
    /// Emit every row from every group — used for `result.roc.all.csv.gz`.
    All,
    Bundle {
        ty: &'a str,
        subset: &'a str,
        filter: &'a str,
        genotype: &'a str,
        qq_field: &'a str,
    },
}

#[derive(Clone, Copy, Debug)]
struct RenderConfig<'a> {
    subset_size: usize,
    whole_reference_size: usize,
    conf_size: usize,
    subset_sizes: &'a BTreeMap<String, usize>,
    subset_confidence_sizes: &'a BTreeMap<String, usize>,
    delta: f64,
    ci_alpha: f64,
    filter_counts_only: bool,
}

/// Build a per-subtype sorted-obs map keyed by `(Type, Subtype, Subset, Filter)`.
///
/// Legacy `RocOutput::write` iterates `(subtype, genotype)` pairs in a
/// fixed order and calls `Roc::getLevels(flag_mask)` per pair, where
/// `getLevels` runs `std::sort` on the shared obs vector. Since
/// `std::sort` is not idempotent, each call leaves obs in a slightly
/// different order at tied levels. The (Subtype, *) row written to the
/// roc.all CSV reflects the obs state AFTER the corresponding
/// `getLevels` call. Genotype iteration is `[het, hetalt, homalt, *]`,
/// so each subtype contributes 4 sorts.
///
/// For each `(Type, Subset, Filter)` wildcard group, sort the obs vector
/// progressively up to the maximum iteration count needed for that Type
/// (3×4=12 for SNP, 10×4=40 for INDEL), snapshotting at each subtype's
/// sort boundary. Each snapshot reproduces the obs state that legacy's
/// `Roc::getLevels` walks for that subtype.
fn build_star_sorted(
    groups: &BTreeMap<RowKey, BoundedGroupAccum>,
) -> BTreeMap<(String, String, String, String), Vec<ObsRecord>> {
    let mut star_sorted: BTreeMap<(String, String, String, String), Vec<ObsRecord>> =
        BTreeMap::new();
    for (key, accum) in groups {
        if key.subtype == "*" && key.genotype == "*" && !accum.observations.is_disk_backed() {
            let mut obs = accum.observations.small_records();
            let mut current_count: usize = 0;
            // Insert per-subtype snapshots at every 4-sort boundary.
            let snapshots: &[(usize, &str)] = match key.ty.as_str() {
                "SNP" => &[(4, "*"), (8, "ti"), (12, "tv")],
                "INDEL" => &[
                    (4, "*"),
                    (8, "I1_5"),
                    (12, "I6_15"),
                    (16, "I16_PLUS"),
                    (20, "D1_5"),
                    (24, "D6_15"),
                    (28, "D16_PLUS"),
                    (32, "C1_5"),
                    (36, "C6_15"),
                    (40, "C16_PLUS"),
                ],
                _ => &[(4, "*")],
            };
            for &(target, subtype_name) in snapshots {
                while current_count < target {
                    introsort_libstdcpp(&mut obs);
                    current_count += 1;
                }
                star_sorted.insert(
                    (
                        key.ty.clone(),
                        subtype_name.to_string(),
                        key.subset.clone(),
                        key.filter.clone(),
                    ),
                    obs.clone(),
                );
            }
        }
    }
    star_sorted
}

fn legacy_subtype_sort_passes(ty: &str, subtype: &str) -> usize {
    let order: &[&str] = if ty == "SNP" {
        &["*", "ti", "tv"]
    } else {
        &[
            "*", "I1_5", "I6_15", "I16_PLUS", "D1_5", "D6_15", "D16_PLUS", "C1_5", "C6_15",
            "C16_PLUS",
        ]
    };
    order
        .iter()
        .position(|candidate| *candidate == subtype)
        .map_or(4, |index| (index + 1) * 4)
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct RenderBundleKey {
    ty: String,
    subset: String,
    filter: String,
    genotype: String,
    qq_field: String,
}

fn render_rows_parallel(
    groups: &BTreeMap<RowKey, BoundedGroupAccum>,
    star_sorted: &BTreeMap<(String, String, String, String), Vec<ObsRecord>>,
    config: RenderConfig<'_>,
    threads: usize,
) -> Result<RenderedRows> {
    let active_types = groups
        .iter()
        .filter(|(_, accum)| accum.records > 0)
        .map(|(key, _)| key.ty.as_str())
        .collect::<HashSet<_>>();
    let bundles = groups
        .keys()
        .filter(|key| active_types.contains(key.ty.as_str()))
        .map(|key| RenderBundleKey {
            ty: key.ty.clone(),
            subset: key.subset.clone(),
            filter: key.filter.clone(),
            genotype: key.genotype.clone(),
            qq_field: key.qq_field.clone(),
        })
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    if bundles.len() <= 1 || threads <= 1 {
        let result = render_rows(groups, star_sorted, RowFilter::All, config);
        clear_all_sorted_indices(groups);
        return result;
    }

    let worker_count = threads.max(1).min(bundles.len());
    let next_bundle = std::sync::atomic::AtomicUsize::new(0);
    let (sender, receiver) = std::sync::mpsc::sync_channel(worker_count);
    let mut rendered = (0..bundles.len()).map(|_| None).collect::<Vec<_>>();
    let mut first_error = None;
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let sender = sender.clone();
            let next_bundle = &next_bundle;
            let bundles = &bundles;
            handles.push(scope.spawn(move || {
                loop {
                    let index = next_bundle.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    let Some(bundle) = bundles.get(index) else {
                        break;
                    };
                    let result = render_rows(
                        groups,
                        star_sorted,
                        RowFilter::Bundle {
                            ty: &bundle.ty,
                            subset: &bundle.subset,
                            filter: &bundle.filter,
                            genotype: &bundle.genotype,
                            qq_field: &bundle.qq_field,
                        },
                        config,
                    );
                    clear_bundle_sorted_indices(groups, bundle);
                    if sender.send((index, result)).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(sender);
        for (index, result) in receiver {
            match result {
                Ok(rows) if first_error.is_none() => rendered[index] = Some(rows),
                Ok(_) => {}
                Err(error) if first_error.is_none() => first_error = Some(error),
                Err(_) => {}
            }
        }
        for handle in handles {
            if handle.join().is_err() && first_error.is_none() {
                first_error = Some(anyhow::anyhow!("ROC render worker panicked"));
            }
        }
    });
    if let Some(error) = first_error {
        return Err(error);
    }
    merge_rendered_rows(
        rendered
            .into_iter()
            .map(|rows| rows.expect("each ROC render bundle returned a result"))
            .collect(),
    )
}

fn clear_bundle_sorted_indices(
    groups: &BTreeMap<RowKey, BoundedGroupAccum>,
    bundle: &RenderBundleKey,
) {
    for (key, accum) in groups {
        if key.ty == bundle.ty
            && key.subtype == "*"
            && key.subset == bundle.subset
            && key.filter == bundle.filter
            && key.genotype == bundle.genotype
            && key.qq_field == bundle.qq_field
        {
            accum.observations.clear_sorted_indices();
        }
    }
}

fn clear_all_sorted_indices(groups: &BTreeMap<RowKey, BoundedGroupAccum>) {
    for accum in groups.values() {
        accum.observations.clear_sorted_indices();
    }
}

type RenderedSortKey = [String; 7];

fn rendered_sort_key(line: &str) -> Result<RenderedSortKey> {
    line.split(',')
        .take(7)
        .map(str::to_string)
        .collect::<Vec<_>>()
        .try_into()
        .map_err(|fields: Vec<String>| {
            anyhow::anyhow!(
                "rendered ROC row has {} key fields, expected 7",
                fields.len()
            )
        })
}

fn merge_rendered_rows(rows: Vec<RenderedRows>) -> Result<RenderedRows> {
    let mut readers = rows
        .iter()
        .map(|rows| rows.lines())
        .collect::<Result<Vec<_>>>()?;
    let mut current = (0..readers.len()).map(|_| None).collect::<Vec<_>>();
    let mut heap = BinaryHeap::<Reverse<(RenderedSortKey, usize)>>::new();
    let advance = |index: usize,
                   readers: &mut Vec<std::io::Lines<BufReader<File>>>,
                   current: &mut Vec<Option<String>>,
                   heap: &mut BinaryHeap<Reverse<(RenderedSortKey, usize)>>|
     -> Result<()> {
        if let Some(line) = readers[index].next() {
            let line = line?;
            heap.push(Reverse((rendered_sort_key(&line)?, index)));
            current[index] = Some(line);
        }
        Ok(())
    };
    for index in 0..readers.len() {
        advance(index, &mut readers, &mut current, &mut heap)?;
    }

    let mut output = tempfile::NamedTempFile::new()
        .context("failed to create merged rendered ROC report spool")?;
    let mut len = 0usize;
    {
        let mut writer = BufWriter::new(output.as_file_mut());
        while let Some(Reverse((_, index))) = heap.pop() {
            let line = current[index]
                .take()
                .expect("rendered ROC merge entry has a current row");
            writeln!(writer, "{line}")?;
            len += 1;
            advance(index, &mut readers, &mut current, &mut heap)?;
        }
        writer.flush()?;
    }
    Ok(RenderedRows {
        path: output.into_temp_path(),
        len,
    })
}

fn render_rows(
    groups: &BTreeMap<RowKey, BoundedGroupAccum>,
    star_sorted: &BTreeMap<(String, String, String, String), Vec<ObsRecord>>,
    row_filter: RowFilter<'_>,
    config: RenderConfig<'_>,
) -> Result<RenderedRows> {
    let mut output =
        tempfile::NamedTempFile::new().context("failed to create rendered ROC report spool")?;
    let mut out = BufWriter::new(output.as_file_mut());
    let mut row_count = 0usize;
    // `accumulate` pre-seeds empty subtype buckets for types that are present,
    // because legacy reports zero-valued subtype baselines for an observed
    // type. It does not, however, emit a second family of rows for a wholly
    // absent type. Determine presence from actual observations before walking
    // the pre-seeded map.
    let active_types: HashSet<&str> = groups
        .iter()
        .filter(|(_, accum)| accum.records > 0)
        .map(|(key, _)| key.ty.as_str())
        .collect();
    for (key, accum) in groups {
        if !active_types.contains(key.ty.as_str()) {
            continue;
        }
        match row_filter {
            RowFilter::All => {}
            RowFilter::Bundle {
                ty,
                subset,
                filter,
                genotype,
                qq_field,
            } => {
                if key.ty != ty
                    || key.subset != subset
                    || key.filter != filter
                    || key.genotype != genotype
                    || key.qq_field != qq_field
                {
                    continue;
                }
            }
        }
        let is_filter_tier = !is_aggregate_filter(&key.filter);
        let wildcard_key = RowKey {
            ty: key.ty.clone(),
            subtype: "*".to_string(),
            subset: key.subset.clone(),
            filter: key.filter.clone(),
            genotype: key.genotype.clone(),
            qq_field: key.qq_field.clone(),
        };
        let emitted_rows = if key.subtype == "*" {
            accum.stream_with_delta(config.delta)?
        } else if let Some(wildcard) = groups
            .get(&wildcard_key)
            .filter(|wildcard| wildcard.observations.is_disk_backed())
        {
            accum.emit_disk_with_delta_from(
                &wildcard.observations,
                Some(&key.subtype),
                legacy_subtype_sort_passes(&key.ty, &key.subtype),
                config.delta,
            )?
        } else if let Some(shared) = star_sorted.get(&(
            key.ty.clone(),
            key.subtype.clone(),
            key.subset.clone(),
            key.filter.clone(),
        )) {
            accum.emit_with_shared_sort_and_delta(shared, &key.subtype, config.delta)?
        } else {
            accum.stream_with_delta(config.delta)?
        };
        for emitted in emitted_rows {
            let emitted = emitted?;
            if is_filter_tier && config.filter_counts_only && emitted.qq_str != "*" {
                // Per-Filter rows in legacy are baseline-only (QQ='*').
                // No numeric thresholds — skip the synthetic 0.0 and any
                // real numeric buckets that may have landed here.
                continue;
            }
            writeln!(
                out,
                "{}",
                render_row(
                    key,
                    &emitted,
                    &RowSizes {
                        subset_size: config.subset_size,
                        whole_reference_size: config.whole_reference_size,
                        conf_size: config.conf_size,
                        subset_sizes: config.subset_sizes,
                        subset_confidence_sizes: config.subset_confidence_sizes,
                    },
                    is_filter_tier && config.filter_counts_only,
                    config.ci_alpha,
                )
            )?;
            row_count += 1;
            if row_count > MAX_RENDERED_ROC_THRESHOLDS {
                bail!(
                    "ROC report exceeds the {MAX_RENDERED_ROC_THRESHOLDS} rendered-row resource limit"
                );
            }
        }
    }
    out.flush()?;
    drop(out);
    Ok(RenderedRows {
        path: output.into_temp_path(),
        len: row_count,
    })
}

fn render_row(
    key: &RowKey,
    emitted: &EmittedRow,
    sizes: &RowSizes<'_>,
    counts_only: bool,
    ci_alpha: f64,
) -> String {
    let &RowSizes {
        subset_size,
        whole_reference_size,
        conf_size,
        subset_sizes,
        subset_confidence_sizes,
    } = sizes;
    let counts = &emitted.cum;
    let truth_total = counts.truth_total();
    let query_total = counts.query_total();

    // Per-Filter rows (Filter ∉ {ALL, PASS}): legacy zeroes TRUTH.TOTAL,
    // TRUTH.FN, QUERY.TOTAL and renders METRIC.* as the pandas NaN
    // sentinel (empty cell). Only TRUTH.TP / QUERY.TP / QUERY.FP /
    // QUERY.UNK carry real counts (and their het/homalt substats); the
    // TOTAL / FN blocks emit `0` for the count and `.` for substats.
    let is_filter_tier = counts_only;
    let supports_titv = key.ty == "SNP";
    let missing_titv = if conf_size > 0 { "." } else { "" };
    let mut row = Vec::with_capacity(EXTENDED_HEADER.len());

    row.push(key.ty.clone());
    row.push(key.subtype.clone());
    row.push(key.subset.clone());
    row.push(key.filter.clone());
    row.push(key.genotype.clone());
    row.push(key.qq_field.clone());
    row.push(emitted.qq_str.clone());

    // Metrics — reuse the extended.csv helpers so formatting matches byte
    // for byte (pandas repr for finite floats, `"0.0"` for zero-denominator
    // ratios, empty string for zero-denominator F1).
    if is_filter_tier {
        row.push(String::new()); // METRIC.Recall
        row.push(String::new()); // METRIC.Precision
        row.push(String::new()); // METRIC.Frac_NA
        row.push(String::new()); // METRIC.F1_Score
    } else {
        row.push(metric_ratio(counts.truth_tp.total, truth_total.total));
        // Use precision_ratio: returns empty when query records exist but
        // none are TP/FP (all UNKs), 0.0 when query_total == 0.
        row.push(crate::adapters::report::precision_ratio(
            counts.query_tp.total,
            counts.query_fp.total,
            query_total.total,
        ));
        row.push(metric_ratio(counts.query_unk.total, query_total.total));
        // F1 score: empty when precision is undefined (query records exist
        // but precision denominator is 0 — same condition as precision_ratio
        // returning empty).
        let prec_denom = counts.query_tp.total + counts.query_fp.total;
        if prec_denom == 0 && query_total.total > 0 {
            row.push(String::new());
        } else {
            row.push(f1_score(
                counts.truth_tp.total,
                truth_total.total,
                counts.query_tp.total,
                prec_denom,
            ));
        }
    }

    // FP.gt / FP.al: raw integers (`"11"`, not `"11.000000"`). Legacy
    // emits them at every INDEL stratification row (matches extended.csv's
    // per-subtype FP emission) and at SNP non-base rows. Forward the
    // accumulated counters as-is for every row.
    row.push(counts.fp_gt.to_string());
    row.push(counts.fp_al.to_string());

    // Subset.Size / Subset.IS_CONF.Size / Subset.Level — match the per-row
    // convention of write_extended exactly.
    let (sz, conf) = subset_size_cells(
        &key.subset,
        &key.subtype,
        subset_size,
        whole_reference_size,
        conf_size,
        subset_sizes,
        subset_confidence_sizes,
    );
    row.push(sz);
    row.push(conf);
    row.push("0.000000".to_string());

    // Substat blocks (7 cells each). The baseline `QQ="*"` row gets full
    // ti/tv/het/homalt/ratio cells via `append_stats`, mirroring the
    // extended.csv layout. Per-threshold numeric rows mirror legacy's
    // independent-sweep pattern: ti/tv/het/homalt cells show cumulative
    // values only at levels kept by that substat's own roc-delta sweep;
    // otherwise `.` (count) / `""` (ratio). Ratios are computed in-row
    // from whichever pair is present — if either side is missing, the
    // ratio is blank (legacy's NaN/NaN → NaN rendering via pandas).
    let emit_block = |row: &mut Vec<String>, bucket: &CountsBucket| match &emitted.substats {
        None if is_filter_tier || emitted.qq_str == "*" => {
            append_stats_with_missing(row, bucket, supports_titv, missing_titv)
        }
        None => append_roc_stats(row, bucket, supports_titv, missing_titv),
        Some(avail) => {
            row.push(bucket.total.to_string());
            emit_substat_cell(
                row,
                bucket.ti,
                supports_titv && avail.ti,
                supports_titv,
                missing_titv,
            );
            emit_substat_cell(
                row,
                bucket.tv,
                supports_titv && avail.tv,
                supports_titv,
                missing_titv,
            );
            emit_substat_cell(row, bucket.het, avail.het, true, missing_titv);
            emit_substat_cell(row, bucket.homalt, avail.homalt, true, missing_titv);
            if supports_titv && avail.ti && avail.tv {
                row.push(ti_tv_ratio(bucket.ti, bucket.tv));
            } else {
                row.push(String::new());
            }
            if avail.het && avail.homalt {
                row.push(het_hom_ratio(bucket.het, bucket.homalt));
            } else {
                row.push(String::new());
            }
        }
    };
    // Per-Filter rows render TRUTH.TOTAL / TRUTH.FN / QUERY.TOTAL as a
    // hollow "0 + . + …" block: the count cell is "0", the substat
    // cells are "." (or empty for ratios). Mirror that with
    // `emit_zero_block` instead of the normal `emit_block` for those
    // three groups; TRUTH.TP, QUERY.TP, QUERY.FP, QUERY.UNK still
    // carry their real counts.
    let emit_zero_block = |row: &mut Vec<String>| {
        row.push("0".to_string());
        for _ in 0..4 {
            row.push(".".to_string());
        }
        row.push(String::new());
        row.push(String::new());
    };
    if is_filter_tier {
        emit_zero_block(&mut row);
        emit_block(&mut row, &counts.truth_tp);
        emit_zero_block(&mut row);
        emit_zero_block(&mut row);
        emit_block(&mut row, &counts.query_tp);
        emit_block(&mut row, &counts.query_fp);
        emit_block(&mut row, &counts.query_unk);
    } else {
        emit_block(&mut row, &truth_total);
        emit_block(&mut row, &counts.truth_tp);
        emit_block(&mut row, &counts.truth_fn);
        emit_block(&mut row, &query_total);
        emit_block(&mut row, &counts.query_tp);
        emit_block(&mut row, &counts.query_fp);
        emit_block(&mut row, &counts.query_unk);
    }

    if ci_alpha > 0.0 {
        let observations = [
            (counts.truth_tp.total, truth_total.total),
            (
                counts.query_tp.total,
                counts.query_tp.total + counts.query_fp.total,
            ),
            if is_filter_tier {
                // Filter-tier projection zeroes QUERY.TOTAL before the
                // unknown-fraction CI is calculated. Recall and precision
                // still use the retained TP/FN/FP counts above.
                (0, 0)
            } else {
                (counts.query_unk.total, query_total.total)
            },
        ];
        append_ci_cells(&mut row, observations, ci_alpha);
    }

    row.join(",")
}

pub(crate) fn append_ci_cells(
    row: &mut Vec<String>,
    observations: [(usize, usize); 3],
    alpha: f64,
) {
    for (successes, trials) in observations {
        let (lower, upper) = jeffreys_interval(successes, trials, alpha);
        row.push(format_ci(lower));
        row.push(format_ci(upper));
    }
}

fn format_ci(value: f64) -> String {
    full_repr_float(value)
}

/// Modified Jeffreys interval used by legacy Tools/ci.py.
pub(crate) fn jeffreys_interval(x: usize, n: usize, alpha: f64) -> (f64, f64) {
    if n == 0 {
        return (0.0, 1.0);
    }
    let lower = if x == n {
        crate::domain::cephes::legacy_pow(alpha / 2.0, 1.0 / n as f64)
    } else if x <= 1 {
        0.0
    } else {
        crate::domain::cephes::incbi(x as f64 + 0.5, (n - x) as f64 + 0.5, alpha / 2.0)
    };
    let upper = if x == 0 {
        1.0 - crate::domain::cephes::legacy_pow(alpha / 2.0, 1.0 / n as f64)
    } else if x >= n - 1 {
        1.0
    } else {
        crate::domain::cephes::incbi(x as f64 + 0.5, (n - x) as f64 + 0.5, 1.0 - alpha / 2.0)
    };
    (lower.max(0.0), upper.min(1.0))
}

fn emit_substat_cell(
    row: &mut Vec<String>,
    value: usize,
    kept: bool,
    supported: bool,
    missing: &str,
) {
    if !supported {
        row.push(missing.to_string());
    } else if kept {
        row.push(format_count(value));
    } else {
        row.push(".".to_string());
    }
}

fn append_roc_stats(
    row: &mut Vec<String>,
    bucket: &CountsBucket,
    supports_titv: bool,
    missing_titv: &str,
) {
    row.push(bucket.total.to_string());
    if supports_titv {
        row.push(format_count(bucket.ti));
        row.push(format_count(bucket.tv));
    } else {
        row.push(missing_titv.to_string());
        row.push(missing_titv.to_string());
    }
    row.push(format_count(bucket.het));
    row.push(format_count(bucket.homalt));
    row.push(if supports_titv {
        ti_tv_ratio(bucket.ti, bucket.tv)
    } else {
        String::new()
    });
    row.push(het_hom_ratio(bucket.het, bucket.homalt));
}

fn subset_size_cells(
    subset: &str,
    _subtype: &str,
    subset_size: usize,
    whole_reference_size: usize,
    conf_size: usize,
    subset_sizes: &BTreeMap<String, usize>,
    subset_confidence_sizes: &BTreeMap<String, usize>,
) -> (String, String) {
    // Derived from the four branches in report::write_extended. Kept in
    // lockstep — if write_extended changes its Subset.Size/IS_CONF.Size
    // convention, this function must follow.
    let is_base_subset = subset == "*";
    let size_cell = if is_base_subset {
        // Subset="*": always the raw subset_size integer, regardless of
        // subtype.
        subset_size.to_string()
    } else if subset == "TS_contained" {
        format_count(conf_size)
    } else if subset == "TS_boundary" {
        format_count(whole_reference_size)
    } else {
        // User-named stratification.
        format_count(subset_sizes.get(subset).copied().unwrap_or(0))
    };
    // The built-in confidence subsets carry the global confidence size.
    // User-named subsets instead carry their interval-union intersection
    // with the confidence regions. Empty cell only when confidence is absent.
    let conf_cell = if conf_size > 0 {
        let size = if matches!(subset, "*" | "TS_boundary" | "TS_contained") {
            conf_size
        } else {
            subset_confidence_sizes.get(subset).copied().unwrap_or(0)
        };
        format_count(size)
    } else {
        String::new()
    };
    (size_cell, conf_cell)
}

// ---------------------------------------------------------------------------
// Gzipped CSV output
// ---------------------------------------------------------------------------

fn write_gzip_csv(path: &Path, header: &str, rows: &RenderedRows) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = GzEncoder::new(file, Compression::default());
    writeln!(writer, "{header}")?;
    for row in rows.lines()? {
        writeln!(writer, "{}", row?)?;
    }
    writer.finish()?;
    Ok(())
}

fn write_optional_gzip_csv(path: &Path, header: &str, rows: &RenderedRows) -> Result<()> {
    if rows.len == 0 {
        // Reusing an output prefix must not retain a stale table from a prior
        // run where this variant type was present.
        if path.exists() {
            std::fs::remove_file(path)
                .with_context(|| format!("failed to remove stale {}", path.display()))?;
        }
        return Ok(());
    }
    write_gzip_csv(path, header, rows)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::SortKey;

    #[test]
    fn substat_cursor_matches_public_roc_bucket_identity() -> Result<()> {
        let mut source = ObservationStore::default();
        source.push(ObsRecord {
            level: f64::MIN_POSITIVE,
            counts: Cumul {
                truth_fn: CountsBucket {
                    total: 1,
                    tv: 1,
                    ..CountsBucket::default()
                },
                ..Cumul::default()
            },
            subtype_bits: subtype_bit("*").expect("wildcard subtype is encoded"),
            ti_flag: false,
            tv_flag: true,
            blt: 0,
        });

        let mut cursor = SubstatSnapshotCursor::build(&source, 8, |row| row.tv_flag)?;
        let counts = cursor
            .counts_at(0.0)?
            .expect("private zero levels share one rendered ROC bucket");

        assert_eq!(rendered_roc_level(0.0), "0.000000");
        assert_eq!(rendered_roc_level(f64::MIN_POSITIVE), "0.000000");
        assert_eq!(counts.truth_fn.tv, 1);
        Ok(())
    }

    #[test]
    fn half_call_is_not_counted_as_heterozygous_in_roc_stats() {
        let format = ["GT", "BD", "BI", "BVT", "BLT", "QQ"];
        let values = ["./1", "TP", "tv", "SNP", "halfcall", "60"];
        let sample = Sample::new(&format, &values);
        let bucket = sample_bucket(&sample);
        assert_eq!(bucket.het, 0);
        assert_eq!(bucket.homalt, 0);
        assert_eq!(bucket.tv, 1);
    }

    #[test]
    fn legacy_string_hash_matches_pinned_libstdcpp() {
        let key = "SNP\t*\t*\tPASS\tTS_contained\t379.290009";
        assert_eq!(legacy_string_hash(key), 0x1dac_92aa_2553_6c1f);
        assert_eq!(legacy_string_hash(key) % 10_273, 5_747);
    }

    #[test]
    fn python27_metrics_dictionary_preserves_pinned_iteration_order() {
        let five_table_insertion = [
            "roc.all",
            "roc.Locations.SNP",
            "roc.Locations.INDEL.PASS",
            "roc.Locations.INDEL",
            "roc.Locations.SNP.PASS",
        ]
        .map(str::to_string);
        assert_eq!(
            python27_dict_iteration_order(&five_table_insertion),
            [
                "roc.Locations.SNP.PASS",
                "roc.all",
                "roc.Locations.INDEL",
                "roc.Locations.SNP",
                "roc.Locations.INDEL.PASS",
            ]
        );

        // The sixth insertion crosses CPython 2.7's two-thirds load factor,
        // so this also pins the resize and rehash path used by SEL reports.
        let seven_table_insertion = [
            "roc.all",
            "roc.Locations.INDEL.SEL",
            "roc.Locations.SNP",
            "roc.Locations.SNP.SEL",
            "roc.Locations.INDEL.PASS",
            "roc.Locations.SNP.PASS",
            "roc.Locations.INDEL",
        ]
        .map(str::to_string);
        assert_eq!(
            python27_dict_iteration_order(&seven_table_insertion),
            [
                "roc.all",
                "roc.Locations.INDEL",
                "roc.Locations.SNP.PASS",
                "roc.Locations.SNP",
                "roc.Locations.INDEL.PASS",
                "roc.Locations.SNP.SEL",
                "roc.Locations.INDEL.SEL",
            ]
        );
    }

    #[test]
    fn introsort_depth_floor_uses_target_pointer_width() {
        assert_eq!(lg_floor(0), 0);
        assert_eq!(lg_floor(1), 0);
        assert_eq!(lg_floor(2), 1);
        assert_eq!(lg_floor(3), 1);
        assert_eq!(lg_floor(16), 4);
        assert_eq!(lg_floor(usize::MAX), usize::BITS as usize - 1);
    }

    #[test]
    fn in_memory_disk_index_sort_matches_random_access_path() -> Result<()> {
        let original = (0..257usize)
            .map(|index| DiskIndexEntry {
                observation_bits: index as u64 * 101,
                level_bits: (((index * 37) % 13) as f64 - 6.0).to_bits(),
            })
            .collect::<Vec<_>>();

        for passes in [1, 4, 12, 40] {
            let mut expected_file = tempfile::NamedTempFile::new()?;
            write_disk_index_entries(expected_file.as_file_mut(), &original)?;
            for _ in 0..passes {
                disk_introsort_libstdcpp(expected_file.as_file_mut(), original.len())?;
            }
            let expected = read_disk_index_entries(expected_file.as_file_mut(), original.len())?;

            let mut actual = original.clone();
            for _ in 0..passes {
                introsort_libstdcpp_by(&mut actual, |left, right| {
                    f64::from_bits(left.level_bits) < f64::from_bits(right.level_bits)
                });
            }
            assert_eq!(actual, expected, "sort order differs after {passes} passes");
        }
        Ok(())
    }

    #[test]
    fn compact_observation_round_trips_every_encoded_axis() -> Result<()> {
        let observation = ObsRecord {
            level: 42.25,
            counts: Cumul {
                truth_tp: CountsBucket {
                    total: 1,
                    ti: 1,
                    het: 1,
                    ..CountsBucket::default()
                },
                query_fp: CountsBucket {
                    total: 1,
                    tv: 1,
                    homalt: 1,
                    ..CountsBucket::default()
                },
                fp_gt: 1,
                ..Cumul::default()
            },
            subtype_bits: encode_subtypes(&[
                "*".to_string(),
                "I1_5".to_string(),
                "D16_PLUS".to_string(),
            ]),
            ti_flag: true,
            tv_flag: true,
            blt: encode_blt(Some("homalt")),
        };
        let encoded = encode_compact_observation(&observation)?;
        let decoded = decode_compact_observation(encoded, observation.level.to_bits())?;

        assert_eq!(decoded.level.to_bits(), observation.level.to_bits());
        assert_eq!(encode_compact_observation(&decoded)?, encoded);
        assert!(decoded.has_subtype("*"));
        assert!(decoded.has_subtype("I1_5"));
        assert!(decoded.has_subtype("D16_PLUS"));
        assert!(decoded.blt_is("homalt"));
        Ok(())
    }

    fn annotated(
        chrom: &str,
        pos: usize,
        qual: &str,
        samples: [&str; 2],
        regions: &str,
        query_pass: bool,
        fp_class: Option<FpClass>,
    ) -> AnnotatedRow {
        let [truth_sample, query_sample] = samples;
        let regions_tag = if regions.is_empty() {
            String::new()
        } else {
            format!(";Regions={regions}")
        };
        let line = format!(
            "{chrom}\t{pos}\t.\tA\tT\t{qual}\t.\tBS=1{regions_tag}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_sample}\t{query_sample}"
        );
        AnnotatedRow {
            sort_key: SortKey::new(chrom.to_string(), pos, 1, 0),
            record: crate::domain::ComparisonRecord::fixture(&line),
            query_pass,
            fp_class,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }
    }

    // Per-side ROC threshold pin: legacy reads each sample's FORMAT.QQ
    // for its own ROC sweep. On chr21 multi-allelic indels the record
    // QUAL is `0` while TRUTH.QQ carries the matched per-side quality
    // (e.g. `829.15`); reading QUAL bins these into the QQ=0 bucket and
    // shifts ~800 truth-TPs out of upper QQ thresholds. Reading
    // FORMAT.QQ keeps each side's records in the right cumulative
    // bucket.
    #[test]
    fn truth_and_query_use_per_side_format_qq() {
        let rows = vec![
            // Multi-allelic-style row: record QUAL=0 but per-side QQ
            // diverge (truth=500, query=0). Truth contributions should
            // land in the 500.000000 bucket, query in 0.000000.
            annotated(
                "chr1",
                100,
                "0",
                ["0/1:TP:gm:tv:SNP:het:500", "0/1:TP:gm:tv:SNP:het:0"],
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");

        // Truth-side bucket at qq=500.0 must hold the TP. (Reading
        // record QUAL=0 would put it in the 0.000000 bucket instead.)
        let bucket_500 = accum
            .threshold_window
            .get("500.000000")
            .expect("missing 500.000000 bucket");
        assert_eq!(bucket_500.counts.truth_tp.total, 1);
        assert_eq!(bucket_500.counts.query_tp.total, 0);

        // Query-side bucket at qq=0.0 must hold the matching query TP.
        let bucket_0 = accum
            .threshold_window
            .get("0.000000")
            .expect("missing 0.000000 bucket");
        assert_eq!(bucket_0.counts.truth_tp.total, 0);
        assert_eq!(bucket_0.counts.query_tp.total, 1);
    }

    #[test]
    fn cumulates_a_single_snp_group_across_qq_thresholds() {
        // Four rows, all SNP TP/TP matches, in subset="*". QUAL values are
        // 10, 20, 30, and ".". The "." row lands only in the baseline.
        let rows = vec![
            annotated(
                "chr1",
                100,
                "10",
                ["0/1:TP:gm:tv:SNP:het:10", "0/1:TP:gm:tv:SNP:het:10"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "20",
                ["1/1:TP:gm:ti:SNP:homalt:20", "1/1:TP:gm:ti:SNP:homalt:20"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                ["0/1:TP:gm:ti:SNP:het:30", "0/1:TP:gm:ti:SNP:het:30"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                ".",
                ["0/1:TP:gm:tv:SNP:het:.", "0/1:TP:gm:tv:SNP:het:."],
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        // Baseline sums all 4 rows; the "." row maps to level=0 (mirroring
        // legacy's `if(std::isnan(qq)) qq = 0` in BlockQuantify::observe),
        // producing a fourth numeric threshold at "0.000000". Lex-ASC
        // string order on integer-valued QQs happens to match numeric order
        // here.
        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        assert_eq!(
            qq_strs,
            vec!["*", "0.000000", "10.000000", "20.000000", "30.000000"]
        );

        // Cumulative semantics (strict-above): row at level L reports
        //   tp = total − cum_through_first_at_L
        // The obs vector contains BOTH a truth-side and a query-side
        // observation per row, so 8 records total. The exact tp values
        // at tied levels depend on libstdc++ tie-break ordering between
        // the truth and query siblings, which is implementation-defined
        // and not stable for std::sort. We only assert structural
        // properties: monotonic decreasing across thresholds, baseline
        // holds total, and the highest threshold reaches the boundary.
        let tp_totals: Vec<usize> = emitted.iter().map(|r| r.cum.truth_tp.total).collect();
        assert_eq!(tp_totals[0], 4, "baseline truth_tp should sum all rows");
        for win in tp_totals.windows(2) {
            assert!(win[0] >= win[1], "tp_totals must monotonically decrease");
        }
        let query_totals: Vec<usize> = emitted.iter().map(|r| r.cum.query_tp.total).collect();
        assert_eq!(query_totals[0], 4);
        for win in query_totals.windows(2) {
            assert!(win[0] >= win[1]);
        }
    }

    #[test]
    fn splits_contributions_across_axes_for_pass_snp_with_region() {
        // One SNP PASS row in built-in and named regions contributes to every
        // non-CONF region axis.
        let rows = vec![annotated(
            "chr1",
            100,
            "42",
            ["0/1:TP:gm:tv:SNP:het:42", "0/1:TP:gm:tv:SNP:het:42"],
            "CONF,TS_contained,EXTRA",
            true,
            None,
        )];
        let groups = accumulate(&rows);
        // accumulate now pre-seeds empty baseline entries for every
        // expected (type, subtype, subset, filter) combo so legacy's
        // empty-subtype baseline rows appear in roc.all. Restrict the
        // assertion to SNP groups with the row's observed axes.
        let snp_combos: std::collections::BTreeSet<(String, String, String)> = groups
            .iter()
            .filter(|(k, accum)| {
                k.ty == "SNP" && k.subtype == "*" && accum.baseline.truth_tp.total > 0
            })
            .map(|(k, _)| (k.subset.clone(), k.filter.clone(), k.subtype.clone()))
            .collect();
        let expected: std::collections::BTreeSet<(String, String, String)> = [
            ("*", "ALL", "*"),
            ("*", "PASS", "*"),
            ("TS_contained", "ALL", "*"),
            ("TS_contained", "PASS", "*"),
            ("EXTRA", "ALL", "*"),
            ("EXTRA", "PASS", "*"),
        ]
        .into_iter()
        .map(|(a, b, c)| (a.to_string(), b.to_string(), c.to_string()))
        .collect();
        let combos = snp_combos;
        assert_eq!(combos, expected);
    }

    #[test]
    fn render_emits_expected_cell_formats() {
        // INDEL baseline row with all zero counts: ti/tv cells should be
        // empty, TiTv_ratio should be empty, denom-zero metrics should render
        // "0.0" (matching `metric_ratio`).
        let key = RowKey::new("INDEL", "*", "*", "ALL");
        let emitted = EmittedRow {
            qq_str: "*".to_string(),
            cum: Cumul::default(),
            substats: None,
        };
        let subset_confidence_sizes = BTreeMap::new();
        let subset_sizes = BTreeMap::new();
        let rendered = render_row(
            &key,
            &emitted,
            &RowSizes {
                subset_size: 100,
                whole_reference_size: 140,
                conf_size: 50,
                subset_sizes: &subset_sizes,
                subset_confidence_sizes: &subset_confidence_sizes,
            },
            false,
            0.0,
        );
        let cells: Vec<&str> = rendered.split(',').collect();
        assert_eq!(cells.len(), 65, "expected 65 columns, got {}", cells.len());
        assert_eq!(cells[0], "INDEL");
        assert_eq!(cells[1], "*");
        assert_eq!(cells[6], "*");
        // Zero-denominator metrics render as "0.0" (legacy metric_ratio).
        assert_eq!(cells[7], "0.0");
        // FP.gt / FP.al are raw integers on the base row.
        assert_eq!(cells[11], "0");
        assert_eq!(cells[12], "0");
        // Subset.Size = raw subset_size integer at Subset="*".
        assert_eq!(cells[13], "100");
        // Confidence regions make unsupported INDEL ti/tv cells use `.`.
        assert_eq!(cells[17], ".");
        assert_eq!(cells[18], ".");
        // TiTv_ratio: empty for INDEL.
        assert_eq!(cells[21], "");

        // Het/hom ratio: both zero → empty.
        let het_hom = het_hom_ratio(0, 0);
        assert_eq!(het_hom, "");

        let without_confidence = render_row(
            &key,
            &emitted,
            &RowSizes {
                subset_size: 100,
                whole_reference_size: 140,
                conf_size: 0,
                subset_sizes: &subset_sizes,
                subset_confidence_sizes: &subset_confidence_sizes,
            },
            false,
            0.0,
        );
        let cells = without_confidence.split(',').collect::<Vec<_>>();
        assert_eq!(cells[17], "");
        assert_eq!(cells[18], "");
    }

    #[test]
    fn subset_size_cells_use_boundary_reference_and_named_confidence_intersection() {
        let subset_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);
        let subset_confidence_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);

        assert_eq!(
            subset_size_cells(
                "TS_boundary",
                "*",
                100,
                140,
                141,
                &subset_sizes,
                &subset_confidence_sizes,
            ),
            ("140.000000".to_string(), "141.000000".to_string())
        );
        assert_eq!(
            subset_size_cells(
                "EXTRA",
                "*",
                100,
                140,
                141,
                &subset_sizes,
                &subset_confidence_sizes,
            ),
            ("138.000000".to_string(), "138.000000".to_string())
        );
        assert_eq!(
            subset_size_cells(
                "*",
                "*",
                100,
                140,
                141,
                &subset_sizes,
                &subset_confidence_sizes,
            ),
            ("100".to_string(), "141.000000".to_string())
        );
    }

    #[test]
    fn absent_variant_type_has_no_roc_rows_or_location_files() {
        let rows = vec![annotated(
            "chr1",
            100,
            "42",
            ["0/1:TP:gm:i1_5:INDEL:het:42", "0/1:TP:gm:i1_5:INDEL:het:42"],
            "",
            true,
            None,
        )];
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("result");
        write_roc_files(&prefix, &rows, 100, 0).unwrap();

        let all = crate::adapters::vcf::read_text(&suffixed_report_path(&prefix, "roc.all.csv.gz"))
            .unwrap();
        assert!(all.lines().skip(1).all(|line| line.starts_with("INDEL,")));
        assert!(suffixed_report_path(&prefix, "roc.Locations.INDEL.csv.gz").exists());
        assert!(suffixed_report_path(&prefix, "roc.Locations.INDEL.PASS.csv.gz").exists());
        assert!(!suffixed_report_path(&prefix, "roc.Locations.SNP.csv.gz").exists());
        assert!(!suffixed_report_path(&prefix, "roc.Locations.SNP.PASS.csv.gz").exists());
    }

    #[test]
    fn configured_unobserved_subset_has_zero_baselines_for_each_active_type() {
        let rows = vec![
            annotated(
                "chr1",
                100,
                "42",
                ["0/1:TP:gm:tv:SNP:het:42", "0/1:TP:gm:tv:SNP:het:42"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "41",
                ["0/1:TP:gm:i1_5:INDEL:het:41", "0/1:TP:gm:i1_5:INDEL:het:41"],
                "",
                true,
                None,
            ),
        ];
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("result");
        let options = RocOptions {
            subset_sizes: BTreeMap::from([("EXTRA_unused".to_string(), 7)]),
            ..RocOptions::default()
        };
        write_roc_files_with_options(&prefix, &rows, 100, 0, &options).unwrap();

        let all = crate::adapters::vcf::read_text(&suffixed_report_path(&prefix, "roc.all.csv.gz"))
            .unwrap();
        let unused = all
            .lines()
            .skip(1)
            .filter(|line| line.split(',').nth(2) == Some("EXTRA_unused"))
            .map(|line| line.split(',').collect::<Vec<_>>())
            .collect::<Vec<_>>();

        let actual_axes = unused
            .iter()
            .map(|fields| (fields[0], fields[1], fields[3], fields[6]))
            .collect::<BTreeSet<_>>();
        let mut expected_axes = BTreeSet::new();
        for filter in ["ALL", "PASS"] {
            expected_axes.insert(("SNP", "*", filter, "*"));
            expected_axes.insert(("INDEL", "*", filter, "*"));
            for subtype in INDEL_SUBTYPES {
                expected_axes.insert(("INDEL", subtype, filter, "*"));
            }
        }
        assert_eq!(actual_axes, expected_axes);
        assert!(unused.iter().all(|fields| fields[6] == "*"));
        assert!(unused.iter().all(|fields| fields[13] == "7.000000"));
        assert!(unused.iter().all(|fields| fields[16] == "0"));
    }

    #[test]
    fn threaded_roc_render_is_byte_deterministic() {
        let rows = vec![
            annotated(
                "chr1",
                100,
                "42",
                ["0/1:TP:gm:tv:SNP:het:42", "0/1:TP:gm:tv:SNP:het:42"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "41",
                ["0/1:TP:gm:i1_5:INDEL:het:41", "0/1:TP:gm:i1_5:INDEL:het:41"],
                "",
                true,
                None,
            ),
        ];
        let dir = tempfile::tempdir().unwrap();
        let serial_prefix = dir.path().join("serial");
        let threaded_prefix = dir.path().join("threaded");
        let serial_options = RocOptions {
            threads: 1,
            ..RocOptions::default()
        };
        let threaded_options = RocOptions {
            threads: 6,
            ..serial_options.clone()
        };
        let serial =
            write_roc_files_with_options(&serial_prefix, &rows, 100, 0, &serial_options).unwrap();
        let threaded =
            write_roc_files_with_options(&threaded_prefix, &rows, 100, 0, &threaded_options)
                .unwrap();

        assert_eq!(serial.tables, threaded.tables);
        assert_eq!(serial.table_order, threaded.table_order);
        for suffix in [
            "roc.all.csv.gz",
            "roc.Locations.SNP.csv.gz",
            "roc.Locations.SNP.PASS.csv.gz",
            "roc.Locations.INDEL.csv.gz",
            "roc.Locations.INDEL.PASS.csv.gz",
        ] {
            assert_eq!(
                std::fs::read(suffixed_report_path(&serial_prefix, suffix)).unwrap(),
                std::fs::read(suffixed_report_path(&threaded_prefix, suffix)).unwrap(),
                "threaded output differs for {suffix}"
            );
        }
    }

    #[test]
    fn pass_tier_demotes_tp_when_query_filtered() {
        // Truth/query TP but query_pass=false. Expected: ALL counts truth
        // as TP, PASS counts truth as FN and query contributes nothing.
        let rows = vec![annotated(
            "chr1",
            100,
            "10",
            ["0/1:TP:gm:tv:SNP:het:10", "0/1:TP:gm:tv:SNP:het:10"],
            "",
            false,
            None,
        )];
        let groups = accumulate(&rows);

        let all_key = RowKey::new("SNP", "*", "*", "ALL");
        let pass_key = RowKey::new("SNP", "*", "*", "PASS");
        let all = groups.get(&all_key).expect("ALL group missing").emit();
        let pass = groups.get(&pass_key).expect("PASS group missing").emit();

        // ALL baseline: 1 TP on truth + 1 TP on query.
        assert_eq!(all[0].cum.truth_tp.total, 1);
        assert_eq!(all[0].cum.truth_fn.total, 0);
        assert_eq!(all[0].cum.query_tp.total, 1);

        // PASS baseline: truth demoted to FN, query contributes nothing.
        assert_eq!(pass[0].cum.truth_tp.total, 0);
        assert_eq!(pass[0].cum.truth_fn.total, 1);
        assert_eq!(pass[0].cum.query_tp.total, 0);
    }

    #[test]
    fn pass_tier_ignores_filtered_queries_without_a_roc_decision() {
        let rows = vec![
            annotated(
                "chr1",
                100,
                "10",
                ["0/0:.:.:.:NOCALL:homref:.", "0/1:UNK:gm:tv:SNP:het:10"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "100",
                ["0/0:.:.:.:NOCALL:homref:.", "0/1:.:gm:tv:SNP:het:100"],
                "",
                false,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let pass_key = RowKey::new("SNP", "*", "*", "PASS");
        let pass = groups.get(&pass_key).expect("PASS group missing");

        // Legacy rocEvaluate calls addROCValue only for TP/FP/UNK query
        // decisions. A filtered record with BD=. must not create an N
        // observation or shift the roc-delta threshold sequence.
        assert_eq!(pass.observations.next_serial, 1);
        let qq: Vec<_> = pass.emit().into_iter().map(|row| row.qq_str).collect();
        assert_eq!(qq, vec!["*", "10.000000"]);
    }

    #[test]
    fn roc_delta_keeps_well_spaced_rows() {
        // Three SNP contributions: TP/TP matches at QQ=50 and QQ=30, plus a
        // truth-only FN at QQ=40. All three numeric thresholds are >0.5
        // apart, so legacy's roc_delta=0.5 filter keeps all of them — the
        // rust emit must match. This pins the absence of a 7-tuple dedup
        // (which legacy doesn't apply) and preserves the derived-FN
        // semantics across all kept rows.
        let rows = vec![
            annotated(
                "chr1",
                100,
                "50",
                ["0/1:TP:gm:tv:SNP:het:50", "0/1:TP:gm:tv:SNP:het:50"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "40",
                ["0/1:FN:gm:tv:SNP:het:40", "0/1:FN:gm:tv:SNP:het:40"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                ["0/1:TP:gm:tv:SNP:het:30", "0/1:TP:gm:tv:SNP:het:30"],
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        // Expected: baseline + three numeric rows (lex-ASC by QQ string).
        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        assert_eq!(qq_strs, vec!["*", "30.000000", "40.000000", "50.000000"]);

        // Baseline cum: 2 truth TP + 1 truth FN + 2 query TP (raw).
        // Numeric rows use legacy's strict-above formula `tp = total −
        // cum_through_first_at_L`. Exact per-row tp values depend on
        // libstdc++ cluster-tie-break ordering between truth and query
        // observations at the same level, which is implementation-
        // defined. Assert structural properties only.
        let tp_totals: Vec<usize> = emitted.iter().map(|r| r.cum.truth_tp.total).collect();
        assert_eq!(tp_totals[0], 2, "baseline truth_tp should sum the two TPs");
        for win in tp_totals.windows(2) {
            assert!(win[0] >= win[1]);
        }
        let fn_totals: Vec<usize> = emitted.iter().map(|r| r.cum.truth_fn.total).collect();
        assert_eq!(fn_totals[0], 1, "baseline truth_fn should sum the one FN");
        for win in fn_totals.windows(2) {
            assert!(win[0] <= win[1]);
        }
        let qtp_totals: Vec<usize> = emitted.iter().map(|r| r.cum.query_tp.total).collect();
        assert_eq!(qtp_totals[0], 2);
        for win in qtp_totals.windows(2) {
            assert!(win[0] >= win[1]);
        }
    }

    #[test]
    fn roc_delta_drops_rows_within_half_unit() {
        // Five SNP query-TP contributions at closely-spaced QQs. Legacy
        // --roc-delta 0.5 keeps only rows where QQ moves >0.5 from the
        // last kept row (ASC walk, first always kept). Expected kept
        // numeric QQs: 30.0, 30.6, 31.2 — 30.3 and 30.9 are within 0.5
        // of the prior kept row and must be dropped. This is the
        // dominant source of rust/legacy roc.all row-count divergence
        // before the fix (~6700 → ~2572 rows per SNP group).
        let rows = vec![
            annotated(
                "chr1",
                100,
                "30.0",
                ["0/1:TP:gm:tv:SNP:het:30", "0/1:TP:gm:tv:SNP:het:30"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "30.3",
                ["0/1:TP:gm:tv:SNP:het:30.3", "0/1:TP:gm:tv:SNP:het:30.3"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30.6",
                ["0/1:TP:gm:tv:SNP:het:30.6", "0/1:TP:gm:tv:SNP:het:30.6"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                "30.9",
                ["0/1:TP:gm:tv:SNP:het:30.9", "0/1:TP:gm:tv:SNP:het:30.9"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                500,
                "31.2",
                ["0/1:TP:gm:tv:SNP:het:31.2", "0/1:TP:gm:tv:SNP:het:31.2"],
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        // 31.2 rounds through f32 to 31.2000007629... which formats to
        // "31.200001" at 6 decimals — legacy exhibits the same artifact.
        assert_eq!(qq_strs, vec!["*", "30.000000", "30.600000", "31.200001"]);
    }

    #[test]
    fn no_dedup_when_every_qq_changes_the_cumulative_tuple() {
        // Four SNP TP/TP matches at distinct QQs. Every threshold moves both
        // cum.truth_tp.total and cum.query_tp.total by 1, so the 7-tuple is
        // different at every row. Dedup must keep all four numeric rows
        // alongside the baseline — guards against over-dedup regressions.
        let rows = vec![
            annotated(
                "chr1",
                100,
                "10",
                ["0/1:TP:gm:tv:SNP:het:10", "0/1:TP:gm:tv:SNP:het:10"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "20",
                ["0/1:TP:gm:ti:SNP:het:20", "0/1:TP:gm:ti:SNP:het:20"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                ["1/1:TP:gm:ti:SNP:homalt:30", "1/1:TP:gm:ti:SNP:homalt:30"],
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                "40",
                ["0/1:TP:gm:tv:SNP:het:40", "0/1:TP:gm:tv:SNP:het:40"],
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        assert_eq!(emitted.len(), 5);
        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        assert_eq!(
            qq_strs,
            vec!["*", "10.000000", "20.000000", "30.000000", "40.000000",]
        );
    }

    #[test]
    fn custom_roc_field_and_delta_control_thresholds() {
        let mut first = annotated(
            "chr1",
            100,
            "90",
            ["0/1:TP:gm:tv:SNP:het:90", "0/1:TP:gm:tv:SNP:het:90"],
            "",
            true,
            None,
        );
        first
            .record
            .try_update(|record| {
                record.info = record.info.replace("BS=1", "BS=1;SCORE=10.0");
                Ok(())
            })
            .unwrap();
        let mut second = annotated(
            "chr1",
            200,
            "80",
            ["0/1:TP:gm:tv:SNP:het:80", "0/1:TP:gm:tv:SNP:het:80"],
            "",
            true,
            None,
        );
        second
            .record
            .try_update(|record| {
                record.info = record.info.replace("BS=1", "BS=1;SCORE=10.4");
                Ok(())
            })
            .unwrap();
        let options = RocOptions {
            qq_field: "SCORE".to_string(),
            delta: 0.0,
            ..RocOptions::default()
        };
        let groups = accumulate_with_options(&[first, second], &options);
        let key = RowKey::new_with_qq_field("SNP", "*", "*", "ALL", "SCORE");
        let rows = groups[&key].emit_with_delta(options.delta).unwrap();
        assert_eq!(
            rows.iter()
                .map(|row| row.qq_str.as_str())
                .collect::<Vec<_>>(),
            vec!["*", "10.000000", "10.400000"]
        );
        assert!(!groups[&key].threshold_window.contains_key("90.000000"));
        let sorted = build_star_sorted(&groups);
        for subtype in ["*", "ti", "tv"] {
            assert!(
                sorted.contains_key(&(
                    "SNP".to_string(),
                    subtype.to_string(),
                    "*".to_string(),
                    "ALL".to_string(),
                )),
                "custom ROC fields need the same subtype sort snapshots as QUAL"
            );
        }
    }

    #[test]
    fn roc_threshold_source_can_differ_from_reported_field() {
        let row = annotated(
            "chr1",
            100,
            "90",
            ["0/1:TP:gm:tv:SNP:het:90", "0/1:TP:gm:tv:SNP:het:90"],
            "",
            true,
            None,
        );
        let options = RocOptions {
            qq_field: "INFO.SCORE".to_string(),
            score_field: Some("QQ".to_string()),
            delta: 0.0,
            ..RocOptions::default()
        };

        let groups = accumulate_with_options(&[row], &options);
        let key = RowKey::new_with_qq_field("SNP", "*", "*", "ALL", "INFO.SCORE");
        let rows = groups[&key].emit_with_delta(options.delta).unwrap();

        assert_eq!(
            rows.iter()
                .map(|row| row.qq_str.as_str())
                .collect::<Vec<_>>(),
            vec!["*", "90.000000"]
        );
    }

    #[test]
    fn ignored_filter_creates_selective_tier_and_renamed_filter_counts() {
        let mut row = annotated(
            "chr1",
            100,
            "20",
            ["0/1:TP:gm:tv:SNP:het:20", "0/1:TP:gm:tv:SNP:het:20"],
            "",
            false,
            None,
        );
        row.record
            .try_update(|record| {
                record.filter = "LowQual".to_string();
                Ok(())
            })
            .unwrap();
        let options = RocOptions {
            ignored_filters: HashSet::from(["LowQual".to_string()]),
            ..RocOptions::default()
        };
        let groups = accumulate_with_options(&[row], &options);
        let pass = &groups[&RowKey::new_with_qq_field("SNP", "*", "*", "PASS", "QUAL")];
        let selective = &groups[&RowKey::new_with_qq_field("SNP", "*", "*", "SEL", "QUAL")];
        assert_eq!(pass.baseline.truth_fn.total, 1);
        assert_eq!(pass.baseline.query_tp.total, 0);
        assert_eq!(selective.baseline.truth_tp.total, 1);
        assert_eq!(selective.baseline.query_tp.total, 1);
        let ignored =
            &groups[&RowKey::new_with_qq_field("SNP", "*", "*", "SEL_IGN_LowQual", "QUAL")];
        assert_eq!(ignored.baseline.query_tp.total, 1);
    }

    #[test]
    fn roc_regions_preserve_legacy_filter_sweep_switch() {
        let mut row = annotated(
            "chr1",
            100,
            "20",
            ["0/1:TP:gm:tv:SNP:het:20", "0/1:TP:gm:tv:SNP:het:20"],
            "",
            false,
            None,
        );
        row.record
            .try_update(|record| {
                record.filter = "LowQual".to_string();
                Ok(())
            })
            .unwrap();
        let options = RocOptions {
            roc_regions: HashSet::from(["TS_contained".to_string()]),
            ..RocOptions::default()
        };
        let groups = accumulate_with_options(&[row], &options);
        let sorted = build_star_sorted(&groups);
        let rendered = render_rows(
            &groups,
            &sorted,
            RowFilter::All,
            RenderConfig {
                subset_size: 100,
                whole_reference_size: 100,
                conf_size: 0,
                subset_sizes: &BTreeMap::new(),
                subset_confidence_sizes: &BTreeMap::new(),
                delta: 0.5,
                ci_alpha: 0.0,
                filter_counts_only: options.roc_regions.contains("*"),
            },
        )
        .unwrap();
        assert!(rendered.lines().unwrap().any(|line| {
            let line = line.unwrap();
            let fields = line.split(',').collect::<Vec<_>>();
            fields[3] == "LowQual" && fields[6] != "*"
        }));
    }

    #[test]
    fn confidence_interval_columns_and_modified_jeffreys_edges() {
        let header = roc_header(0.05);
        assert_eq!(header.split(',').count(), EXTENDED_HEADER.len() + 6);
        let (empty_lower, empty_upper) = jeffreys_interval(0, 0, 0.05);
        assert_eq!((empty_lower, empty_upper), (0.0, 1.0));
        let (lower, upper) = jeffreys_interval(0, 10, 0.05);
        assert_eq!(lower, 0.0);
        assert!((upper - (1.0 - 0.025_f64.powf(0.1))).abs() < 1e-14);
        let (lower, upper) = jeffreys_interval(5, 10, 0.05);
        assert!((lower - 0.223_528_670_252_705_2).abs() < 1e-12);
        assert!((upper - 0.776_471_329_747_294_7).abs() < 1e-12);
        let (lower, upper) = jeffreys_interval(1037, 1156, 0.05);
        assert_eq!(lower.to_bits(), 0x3fec_1d15_8e8d_75af);
        assert_eq!(upper.to_bits(), 0x3fed_3c11_2b31_7194);

        let (lower, _) = jeffreys_interval(1233, 1233, 0.05);
        assert_eq!(lower.to_bits(), 0x3fef_e787_2242_48a7);
        let (_, upper) = jeffreys_interval(0, 136, 0.05);
        assert_eq!(upper.to_bits(), 0x3f9b_66db_9060_1320);
    }

    #[test]
    fn confidence_interval_csv_renders_full_repr() {
        assert_eq!(format_ci(5.280_579_842_943_484e-5), "5.280579842943484e-05");
        assert_eq!(
            format_ci(0.000_814_921_550_822_522_7),
            "0.0008149215508225227"
        );
    }

    #[test]
    fn filter_tier_unknown_fraction_ci_uses_zero_query_total() {
        let key = RowKey::new("SNP", "*", "*", "LowMQ");
        let emitted = EmittedRow {
            qq_str: "*".to_string(),
            cum: Cumul {
                truth_tp: CountsBucket {
                    total: 2,
                    ..CountsBucket::default()
                },
                query_tp: CountsBucket {
                    total: 1,
                    ..CountsBucket::default()
                },
                query_fp: CountsBucket {
                    total: 3,
                    ..CountsBucket::default()
                },
                query_unk: CountsBucket {
                    total: 4,
                    ..CountsBucket::default()
                },
                ..Cumul::default()
            },
            substats: None,
        };

        let subset_confidence_sizes = BTreeMap::new();
        let subset_sizes = BTreeMap::new();
        let rendered = render_row(
            &key,
            &emitted,
            &RowSizes {
                subset_size: 100,
                whole_reference_size: 140,
                conf_size: 50,
                subset_sizes: &subset_sizes,
                subset_confidence_sizes: &subset_confidence_sizes,
            },
            true,
            0.05,
        );
        let cells = rendered.split(',').collect::<Vec<_>>();
        assert_eq!(
            &cells[cells.len() - 6..],
            [
                "0.15811388300841897",
                "1.0",
                "0.0",
                "0.7162483204365873",
                "0.0",
                "1.0",
            ]
        );
    }
}
