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
    EXTENDED_HEADER, append_ci_cells, append_stats_with_missing, empty_comparison_extended_lines,
    f1_score, format_count, het_hom_ratio, metric_ratio, ti_tv_ratio,
};
use crate::domain::{AnnotatedRow, CountsBucket};
use anyhow::{Result, bail};
use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet, HashSet, VecDeque};

mod accumulation;
mod contributions;
mod legacy;
mod model;
mod rendering;

#[cfg(test)]
use crate::domain::jeffreys_interval;
use accumulation::*;
use contributions::*;
use legacy::*;
use model::*;
use rendering::*;

const INDEL_SUBTYPES: [&str; 9] = [
    "C16_PLUS", "C1_5", "C6_15", "D16_PLUS", "D1_5", "D6_15", "I16_PLUS", "I1_5", "I6_15",
];

/// Write `roc.all` and each non-empty Locations ROC file alongside `prefix`.
#[derive(Clone, Debug, Default)]
pub(crate) struct MetricIndices {
    pub tables: BTreeMap<String, Vec<usize>>,
    /// Table order produced by Python 2.7's insertion-ordered hash table
    /// iteration in `happyroc.roc`. This is data-dependent because location
    /// tables are inserted when their first raw ROC row is encountered.
    pub table_order: Vec<String>,
}

/// One CSV artifact calculated by the engine and published by an adapter.
#[derive(Clone, Debug)]
pub(crate) struct CsvArtifact {
    pub suffix: String,
    pub header: String,
    pub rows: Vec<String>,
    pub optional: bool,
}

/// Filesystem-neutral result of ROC calculation.
#[derive(Clone, Debug)]
pub(crate) struct Artifacts {
    pub indices: MetricIndices,
    pub csv: Vec<CsvArtifact>,
    pub raw_table: Option<String>,
}

/// Controls inherited from qfy's ROC command line.  The default deliberately
/// remains identical to the historical hap.py invocation.
#[derive(Clone, Debug)]
pub(crate) struct RocOptions {
    pub qq_field: String,
    /// Optional FORMAT/INFO field used for thresholds while `qq_field`
    /// remains the user-facing label in metrics and ROC tables.
    pub score_field: Option<String>,
    pub ignored_filters: HashSet<String>,
    pub roc_regions: HashSet<String>,
    pub delta: f64,
    pub ci_alpha: f64,
    /// Preserve qfy's private C++ quantifier table. Legacy qfy removes this
    /// intermediate unless `--verbose` is active.
    pub preserve_raw_table: bool,
    /// Include threshold rows in the private table. This follows qfy's
    /// `--roc`/`--no-roc` switch independently of the public compacting pass.
    pub output_rocs: bool,
    /// Full N-trimmed FASTA size used by the legacy TS_boundary lane.
    /// `None` preserves the historical caller contract where `subset_size`
    /// is also the complete reference size.
    pub whole_reference_size: Option<usize>,
    /// Per-user-named stratification interval-union sizes.
    pub subset_sizes: BTreeMap<String, usize>,
    /// Per-named-subset intersection with the confidence regions.
    pub subset_confidence_sizes: BTreeMap<String, usize>,
}

impl Default for RocOptions {
    fn default() -> Self {
        Self {
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

/// Calculate ROC artifacts without accessing the filesystem.
pub(crate) fn calculate(
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
) -> Result<Artifacts> {
    calculate_with_options(rows, subset_size, conf_size, &RocOptions::default())
}

/// Calculate ROC artifacts with qfy controls supplied by the caller.
pub(crate) fn calculate_with_options(
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<Artifacts> {
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
    let groups = accumulate_with_options(rows, options);

    if !groups.values().any(|group| !group.obs.is_empty()) {
        let all = empty_comparison_extended_lines(subset_size);
        let csv = vec![CsvArtifact {
            suffix: "roc.all.csv.gz".to_string(),
            header: roc_header(options.ci_alpha),
            rows: all,
            optional: false,
        }];
        let indices = vec![1, 0];
        return Ok(Artifacts {
            indices: MetricIndices {
                tables: BTreeMap::from([
                    ("summary.metrics".to_string(), indices.clone()),
                    ("all.metrics".to_string(), indices.clone()),
                    ("roc.all".to_string(), indices),
                ]),
                table_order: vec!["roc.all".to_string()],
            },
            csv,
            raw_table: None,
        });
    }

    let raw_table = options
        .preserve_raw_table
        .then(|| build_legacy_roc_table(&groups, subset_size, conf_size, options));

    // Compute per-subtype star_sorted snapshots ONCE (heavy operation: up to
    // 10×4 sorts × full obs vector clone for INDEL groups). Pass by reference
    // to render_rows so the cost isn't paid 5 times.
    let star_sorted = build_star_sorted(&groups);

    let header = roc_header(options.ci_alpha);
    let mut csv = Vec::new();
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
    let all = render_rows(&groups, &star_sorted, RowFilter::All, render_config);
    csv.push(CsvArtifact {
        suffix: "roc.all.csv.gz".to_string(),
        header: header.clone(),
        rows: all.clone(),
        optional: false,
    });

    // `roc.Locations.<TYPE>[.PASS]`: matches legacy `happyroc.py` Locations
    // filter — keeps only `(Type, Subtype=*, Subset=*, Genotype=*, Filter,
    // QQ != *)` rows. Each file is the per-QQ cumulative threshold sweep
    // for the corresponding Type under a single Filter.
    let snp = render_rows(
        &groups,
        &star_sorted,
        RowFilter::Locations {
            ty: "SNP",
            filter: "ALL",
        },
        render_config,
    );
    csv.push(CsvArtifact {
        suffix: "roc.Locations.SNP.csv.gz".to_string(),
        header: header.clone(),
        rows: snp.clone(),
        optional: true,
    });
    let snp_pass = render_rows(
        &groups,
        &star_sorted,
        RowFilter::Locations {
            ty: "SNP",
            filter: "PASS",
        },
        render_config,
    );
    csv.push(CsvArtifact {
        suffix: "roc.Locations.SNP.PASS.csv.gz".to_string(),
        header: header.clone(),
        rows: snp_pass.clone(),
        optional: true,
    });
    let indel = render_rows(
        &groups,
        &star_sorted,
        RowFilter::Locations {
            ty: "INDEL",
            filter: "ALL",
        },
        render_config,
    );
    csv.push(CsvArtifact {
        suffix: "roc.Locations.INDEL.csv.gz".to_string(),
        header: header.clone(),
        rows: indel.clone(),
        optional: true,
    });
    let indel_pass = render_rows(
        &groups,
        &star_sorted,
        RowFilter::Locations {
            ty: "INDEL",
            filter: "PASS",
        },
        render_config,
    );
    csv.push(CsvArtifact {
        suffix: "roc.Locations.INDEL.PASS.csv.gz".to_string(),
        header: header.clone(),
        rows: indel_pass.clone(),
        optional: true,
    });

    let mut selective = Vec::new();
    if !options.ignored_filters.is_empty() {
        for ty in ["SNP", "INDEL"] {
            let lines = render_rows(
                &groups,
                &star_sorted,
                RowFilter::Locations { ty, filter: "SEL" },
                render_config,
            );
            let id = format!("roc.Locations.{ty}.SEL");
            csv.push(CsvArtifact {
                suffix: format!("roc.Locations.{ty}.SEL.csv.gz"),
                header: header.clone(),
                rows: lines.clone(),
                optional: true,
            });
            selective.push((id, lines, ty));
        }
    }

    let indices = build_metric_indices(
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
    );
    Ok(Artifacts {
        indices,
        csv,
        raw_table,
    })
}

#[cfg(test)]
mod test_suite;
