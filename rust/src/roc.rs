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

use crate::compare::{AnnotatedRow, suffixed_report_path};
use crate::report::{
    CountsBucket, EXTENDED_HEADER, append_stats, f1_score, format_count, het_hom_ratio,
    metric_ratio, ti_tv_ratio,
};
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::cmp::Ordering;
use std::collections::{BTreeMap, HashSet, VecDeque};
use std::io::Write;
use std::path::Path;

const INDEL_SUBTYPES: [&str; 9] = [
    "C16_PLUS", "C1_5", "C6_15", "D16_PLUS", "D1_5", "D6_15", "I16_PLUS", "I1_5", "I6_15",
];

/// Write `roc.all` and each non-empty Locations ROC file alongside `prefix`.
#[derive(Clone, Debug, Default)]
pub struct MetricIndices {
    pub tables: BTreeMap<String, Vec<usize>>,
}

/// Controls inherited from qfy's ROC command line.  The default deliberately
/// remains identical to the historical hap.py invocation.
#[derive(Clone, Debug)]
pub struct RocOptions {
    pub qq_field: String,
    pub ignored_filters: HashSet<String>,
    pub roc_regions: HashSet<String>,
    pub delta: f64,
    pub ci_alpha: f64,
}

impl Default for RocOptions {
    fn default() -> Self {
        Self {
            qq_field: "QUAL".to_string(),
            ignored_filters: HashSet::new(),
            roc_regions: HashSet::from(["*".to_string()]),
            delta: 0.5,
            ci_alpha: 0.0,
        }
    }
}

/// Write ROC artifacts and return the original pandas row indices used by
/// legacy qfy's metrics JSON tables.
pub fn write_roc_files(
    prefix: &Path,
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
) -> Result<MetricIndices> {
    write_roc_files_with_options(prefix, rows, subset_size, conf_size, &RocOptions::default())
}

/// Write ROC artifacts with the qfy controls supplied by the caller.
pub fn write_roc_files_with_options(
    prefix: &Path,
    rows: &[AnnotatedRow],
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> Result<MetricIndices> {
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

    // Compute per-subtype star_sorted snapshots ONCE (heavy operation: up to
    // 10×4 sorts × full obs vector clone for INDEL groups). Pass by reference
    // to render_rows so the cost isn't paid 5 times.
    let star_sorted = build_star_sorted(&groups);

    let header = roc_header(options.ci_alpha);
    let render_config = RenderConfig {
        subset_size,
        conf_size,
        delta: options.delta,
        ci_alpha: options.ci_alpha,
        filter_counts_only: options.roc_regions.contains("*"),
    };
    // `roc.all`: every row from every group, no filtering.
    let all = render_rows(&groups, &star_sorted, RowFilter::All, render_config);
    write_gzip_csv(
        &suffixed_report_path(prefix, "roc.all.csv.gz"),
        &header,
        &all,
    )?;

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
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.SNP.csv.gz"),
        &header,
        &snp,
    )?;
    let snp_pass = render_rows(
        &groups,
        &star_sorted,
        RowFilter::Locations {
            ty: "SNP",
            filter: "PASS",
        },
        render_config,
    );
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.SNP.PASS.csv.gz"),
        &header,
        &snp_pass,
    )?;
    let indel = render_rows(
        &groups,
        &star_sorted,
        RowFilter::Locations {
            ty: "INDEL",
            filter: "ALL",
        },
        render_config,
    );
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.INDEL.csv.gz"),
        &header,
        &indel,
    )?;
    let indel_pass = render_rows(
        &groups,
        &star_sorted,
        RowFilter::Locations {
            ty: "INDEL",
            filter: "PASS",
        },
        render_config,
    );
    write_optional_gzip_csv(
        &suffixed_report_path(prefix, "roc.Locations.INDEL.PASS.csv.gz"),
        &header,
        &indel_pass,
    )?;

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
            write_optional_gzip_csv(
                &suffixed_report_path(prefix, &format!("roc.Locations.{ty}.SEL.csv.gz")),
                &header,
                &lines,
            )?;
            selective.push((id, lines, ty));
        }
    }

    Ok(build_metric_indices(
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
    ))
}

struct MetricRows<'a> {
    all: &'a [String],
    snp: &'a [String],
    snp_pass: &'a [String],
    indel: &'a [String],
    indel_pass: &'a [String],
    selective: &'a [(String, Vec<String>, &'a str)],
}

fn build_metric_indices(
    groups: &BTreeMap<RowKey, GroupAccum>,
    rows: MetricRows<'_>,
    delta: f64,
) -> MetricIndices {
    let mut rocs = BTreeMap::<String, (&RowKey, &GroupAccum)>::new();
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

    let mut table = LegacyUnorderedRows::default();
    for (_, (key, accum)) in rocs {
        let subtype_flags: &[(&str, Option<&str>)] = if key.ty == "SNP" {
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
        };
        let counts_only = !is_aggregate_filter(&key.filter);
        for (subtype, subtype_flag) in subtype_flags {
            for (genotype, genotype_flag) in [
                ("het", Some("het")),
                ("hetalt", Some("hetalt")),
                ("homalt", Some("homalt")),
                ("*", None),
            ] {
                let baseline = legacy_row_key(&key.ty, subtype, &key.filter, &key.subset, "*");
                if !matches!(*subtype, "ti" | "tv") {
                    table.set(baseline, genotype == "*");
                } else if genotype == "*" {
                    table.set(
                        legacy_row_key(&key.ty, "*", &key.filter, &key.subset, "*"),
                        false,
                    );
                }
                if counts_only {
                    continue;
                }
                for level in legacy_masked_levels(&accum.obs, *subtype_flag, genotype_flag, delta) {
                    let qq = format!("{level:.6}");
                    if !matches!(*subtype, "ti" | "tv") && genotype == "*" {
                        table.set(
                            legacy_row_key(&key.ty, subtype, &key.filter, &key.subset, &qq),
                            true,
                        );
                    } else if matches!(*subtype, "ti" | "tv") && genotype == "*" {
                        table.set(
                            legacy_row_key(&key.ty, "*", &key.filter, &key.subset, &qq),
                            false,
                        );
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
                        );
                    }
                }
            }
        }
    }

    let raw = table.retained_order();
    let raw_positions = raw
        .iter()
        .enumerate()
        .map(|(index, key)| (key.as_str(), index))
        .collect::<BTreeMap<_, _>>();
    let all_indices = indices_for_lines(rows.all, &raw_positions);

    let mut tables = BTreeMap::new();
    tables.insert("roc.all".to_string(), all_indices.clone());
    tables.insert(
        "all.metrics".to_string(),
        rows.all
            .iter()
            .zip(&all_indices)
            .filter(|(line, _)| {
                let fields = line.split(',').collect::<Vec<_>>();
                fields.get(6) == Some(&"*") && matches!(fields.get(3), Some(&"ALL") | Some(&"PASS"))
            })
            .map(|(_, index)| *index)
            .collect(),
    );
    tables.insert(
        "summary.metrics".to_string(),
        rows.all
            .iter()
            .zip(&all_indices)
            .filter(|(line, _)| {
                let fields = line.split(',').collect::<Vec<_>>();
                fields.get(1) == Some(&"*")
                    && fields.get(2) == Some(&"*")
                    && matches!(fields.get(3), Some(&"ALL") | Some(&"PASS"))
                    && fields.get(4) == Some(&"*")
                    && fields.get(6) == Some(&"*")
            })
            .map(|(_, index)| *index)
            .collect(),
    );
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
                local.insert(raw_key.as_str(), next);
                next += 1;
            }
        }
        tables.insert(id.to_string(), indices_for_lines(lines, &local));
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
                local.insert(raw_key.as_str(), next);
                next += 1;
            }
        }
        tables.insert(id.clone(), indices_for_lines(lines, &local));
    }
    MetricIndices { tables }
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

fn indices_for_lines(lines: &[String], positions: &BTreeMap<&str, usize>) -> Vec<usize> {
    lines
        .iter()
        .map(|line| {
            positions
                .get(csv_row_key(line).as_str())
                .copied()
                .unwrap_or(0)
        })
        .collect()
}

fn legacy_masked_levels(
    obs: &[ObsRecord],
    subtype: Option<&str>,
    genotype: Option<&str>,
    delta: f64,
) -> Vec<f64> {
    let mut levels = obs
        .iter()
        .filter(|record| {
            let subtype_matches = match subtype {
                None => true,
                Some("ti") => record.ti_flag,
                Some("tv") => record.tv_flag,
                Some(value) => record.subtypes.iter().any(|candidate| candidate == value),
            };
            subtype_matches && genotype.is_none_or(|value| record.blt.as_deref() == Some(value))
        })
        .map(|record| record.level)
        .collect::<Vec<_>>();
    levels.sort_by(|left, right| left.partial_cmp(right).unwrap_or(Ordering::Equal));
    let mut kept = Vec::new();
    let mut previous = None;
    for level in levels {
        if previous.is_none_or(|value: f64| (level - value).abs() > delta) {
            kept.push(level);
            previous = Some(level);
        }
    }
    kept
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
    fn set(&mut self, key: String, has_type: bool) {
        if has_type {
            self.typed.insert(key.clone());
        }
        if self.seen.contains(&key) {
            return;
        }
        if self.bucket_count == 0 {
            self.bucket_count = 1;
        }
        if self.seen.len() + 1 > self.bucket_count {
            self.rehash(next_legacy_bucket(self.bucket_count, self.seen.len() + 1));
        }
        self.seen.insert(key.clone());
        let bucket = legacy_string_hash(&key) as usize % self.bucket_count;
        if let Some(entries) = self.buckets.get_mut(&bucket) {
            entries.push_front(key);
        } else {
            self.buckets.insert(bucket, VecDeque::from([key]));
            self.bucket_order.push_front(bucket);
        }
    }

    fn rehash(&mut self, bucket_count: usize) {
        let old = self.all_rows();
        self.bucket_count = bucket_count;
        self.buckets.clear();
        self.bucket_order.clear();
        for key in old {
            let bucket = legacy_string_hash(&key) as usize % bucket_count;
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
    subtypes: Vec<String>,
    /// True if the original BI tag included `ti` — needed to mirror
    /// legacy's `OBS_FLAG_TI` mask which is set per-record from the BI
    /// string, *not* from non-zero count fields. Filter-failed query
    /// phantoms have all-zero counts but still carry their flag bits.
    ti_flag: bool,
    tv_flag: bool,
    blt: Option<String>,
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
///   • `numeric_buckets`: per-`{:.6}` bucket aggregate, retained for
///     the per-substat (ti/tv/het/homalt) roc-delta sweeps and for
///     bookkeeping helpers like `bucket_has_ti`.
#[derive(Default)]
struct GroupAccum {
    baseline: Cumul,
    obs: Vec<ObsRecord>,
    numeric_buckets: BTreeMap<String, NumericBucket>,
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

/// One row emitted by a `GroupAccum`. For the baseline row (qq_str = "*"),
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

impl GroupAccum {
    fn add(
        &mut self,
        qq: Option<f64>,
        counts: &Cumul,
        subtypes: &[String],
        bi: Option<&str>,
        blt: Option<&str>,
    ) {
        self.baseline.add(counts);
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
        let key = format!("{q32:.6}");
        let entry = self
            .numeric_buckets
            .entry(key)
            .or_insert_with(|| NumericBucket {
                qq: q32,
                counts: Cumul::default(),
                has_ti: false,
                has_tv: false,
            });
        entry.counts.add(counts);
        // Track substat-flag presence: any obs at this level with ti/tv
        // contribution flags this bucket for substat sweep. Detected from
        // either counts.{*}.ti/tv (sample_bucket sets these from bi) OR
        // bi parameter (handles filter-failed phantoms whose counts are
        // zero but who carry bi flags per legacy's BlockQuantify::observe).
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
        self.obs.push(ObsRecord {
            level: obs_level,
            counts: counts.clone(),
            subtypes: subtypes.to_vec(),
            ti_flag,
            tv_flag,
            blt: blt.map(str::to_string),
        });
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
    /// Standard emit — sorts self.obs with libstdc++ introsort.
    /// Used for the `Subtype="*"` group and as a fallback for callers
    /// that don't have a pre-sorted obs vector to share.
    #[cfg(test)]
    fn emit(&self) -> Vec<EmittedRow> {
        self.emit_with_delta(0.5)
    }

    fn emit_with_delta(&self, delta: f64) -> Vec<EmittedRow> {
        self.emit_internal(None, None, delta)
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
    ) -> Vec<EmittedRow> {
        self.emit_internal(Some(shared_sorted), Some(my_subtype), delta)
    }

    fn emit_internal(
        &self,
        shared_sorted: Option<&[ObsRecord]>,
        my_subtype: Option<&str>,
        delta: f64,
    ) -> Vec<EmittedRow> {
        let truth_total_const = sum_buckets(&self.baseline.truth_tp, &self.baseline.truth_fn);

        let mut out = Vec::with_capacity(1 + self.numeric_buckets.len());
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
        // at insertion time (`GroupAccum::add`), matching legacy's
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
                    .filter(|o| o.subtypes.iter().any(|x| x == st))
                    .cloned()
                    .collect();
                (filtered, true)
            }
            _ => (self.obs.clone(), false),
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
            let mut truth_fn = add_buckets_total(&total_truth_fn, &cum_through.truth_tp);
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
            let qq_str = format!("{level_f32:.6}");
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
        let zero_already_bucket = self.numeric_buckets.contains_key(&zero_key);
        let need_synthetic_zero_bucket = truth_total_const.total > 0 && !zero_already_bucket;
        let mut ordered: Vec<(&String, &NumericBucket)> = self.numeric_buckets.iter().collect();
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

fn add_buckets_total(a: &CountsBucket, b: &CountsBucket) -> CountsBucket {
    CountsBucket {
        total: a.total + b.total,
        ti: a.ti + b.ti,
        tv: a.tv + b.tv,
        het: a.het + b.het,
        homalt: a.homalt + b.homalt,
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
        (usize::BITS - 1 - (n as u64).leading_zeros()) as usize
    }
}

fn introsort_libstdcpp(arr: &mut [ObsRecord]) {
    let n = arr.len();
    if n > 1 {
        let depth_limit = lg_floor(n) * 2;
        introsort_loop(arr, 0, n, depth_limit);
        final_insertion_sort(arr);
    }
}

fn introsort_loop(arr: &mut [ObsRecord], first: usize, mut last: usize, mut depth_limit: usize) {
    while last - first > ROC_SORT_THRESHOLD {
        if depth_limit == 0 {
            heapsort_range(arr, first, last);
            return;
        }
        depth_limit -= 1;
        let cut = unguarded_partition_pivot(arr, first, last);
        introsort_loop(arr, cut, last, depth_limit);
        last = cut;
    }
}

fn unguarded_partition_pivot(arr: &mut [ObsRecord], first: usize, last: usize) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(arr, first, first + 1, mid, last - 1);
    unguarded_partition(arr, first + 1, last, first)
}

fn move_median_to_first(arr: &mut [ObsRecord], result: usize, a: usize, b: usize, c: usize) {
    if arr[a].level < arr[b].level {
        if arr[b].level < arr[c].level {
            arr.swap(result, b);
        } else if arr[a].level < arr[c].level {
            arr.swap(result, c);
        } else {
            arr.swap(result, a);
        }
    } else if arr[a].level < arr[c].level {
        arr.swap(result, a);
    } else if arr[b].level < arr[c].level {
        arr.swap(result, c);
    } else {
        arr.swap(result, b);
    }
}

fn unguarded_partition(
    arr: &mut [ObsRecord],
    mut first: usize,
    mut last: usize,
    pivot: usize,
) -> usize {
    loop {
        while arr[first].level < arr[pivot].level {
            first += 1;
        }
        last -= 1;
        while arr[pivot].level < arr[last].level {
            last -= 1;
        }
        if first >= last {
            return first;
        }
        arr.swap(first, last);
        first += 1;
    }
}

fn final_insertion_sort(arr: &mut [ObsRecord]) {
    let n = arr.len();
    if n > ROC_SORT_THRESHOLD {
        insertion_sort_range(arr, 0, ROC_SORT_THRESHOLD);
        unguarded_insertion_sort_range(arr, ROC_SORT_THRESHOLD, n);
    } else {
        insertion_sort_range(arr, 0, n);
    }
}

fn insertion_sort_range(arr: &mut [ObsRecord], first: usize, last: usize) {
    if first == last {
        return;
    }
    for i in (first + 1)..last {
        if arr[i].level < arr[first].level {
            arr[first..=i].rotate_right(1);
        } else {
            let mut j = i;
            while j > first && arr[j].level < arr[j - 1].level {
                arr.swap(j, j - 1);
                j -= 1;
            }
        }
    }
}

fn unguarded_insertion_sort_range(arr: &mut [ObsRecord], first: usize, last: usize) {
    for i in first..last {
        let mut j = i;
        while j > 0 && arr[j].level < arr[j - 1].level {
            arr.swap(j, j - 1);
            j -= 1;
        }
    }
}

fn heapsort_range(arr: &mut [ObsRecord], first: usize, last: usize) {
    let n = last - first;
    if n < 2 {
        return;
    }
    for i in (0..n / 2).rev() {
        sift_down(arr, first + i, first, last);
    }
    for i in (1..n).rev() {
        arr.swap(first, first + i);
        sift_down(arr, first, first, first + i);
    }
}

fn sift_down(arr: &mut [ObsRecord], start: usize, first: usize, last: usize) {
    let mut root = start;
    loop {
        let lc = first + 2 * (root - first) + 1;
        if lc >= last {
            break;
        }
        let rc = lc + 1;
        let mut child = lc;
        if rc < last && arr[lc].level < arr[rc].level {
            child = rc;
        }
        if arr[root].level < arr[child].level {
            arr.swap(root, child);
            root = child;
        } else {
            break;
        }
    }
}

#[cfg(test)]
fn accumulate(rows: &[AnnotatedRow]) -> BTreeMap<RowKey, GroupAccum> {
    accumulate_impl(rows, &RocOptions::default())
}

fn accumulate_impl(rows: &[AnnotatedRow], options: &RocOptions) -> BTreeMap<RowKey, GroupAccum> {
    let mut groups: BTreeMap<RowKey, GroupAccum> = BTreeMap::new();

    // Populate first from actual contributions so we learn which subsets
    // (TS_boundary, TS_contained) the dataset has any records in.
    let mut observed_subsets: std::collections::BTreeSet<String> =
        std::collections::BTreeSet::new();
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
                    .add(qq, counts, subtypes, bi, blt);
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

    groups
}

fn accumulate_with_options(
    rows: &[AnnotatedRow],
    options: &RocOptions,
) -> BTreeMap<RowKey, GroupAccum> {
    accumulate_impl(rows, options)
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
        self.bvt.filter(|value| *value != "NOCALL")
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
    let fields: Vec<&str> = row.line.split('\t').collect();
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
    let truth_qq = truth.roc_value(&options.qq_field, fields[5], info);
    let query_qq = query.roc_value(&options.qq_field, fields[5], info);

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
                        if row.fp_class == Some("gt") {
                            counts.fp_gt = 1;
                        } else if row.fp_class == Some("al") {
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
                        if row.fp_class == Some("gt") {
                            counts.fp_gt = 1;
                        } else if row.fp_class == Some("al") {
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
            if tag == "TS_boundary" || tag == "TS_contained" {
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
    /// Legacy `happyroc.py` Locations filter: keep only rows in the base
    /// group `(Type=ty, Subtype=*, Subset=*, Genotype=*, Filter=filter)`
    /// and drop the `QQ="*"` baseline row.
    Locations { ty: &'a str, filter: &'a str },
}

#[derive(Clone, Copy, Debug)]
struct RenderConfig {
    subset_size: usize,
    conf_size: usize,
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
    groups: &BTreeMap<RowKey, GroupAccum>,
) -> BTreeMap<(String, String, String, String), Vec<ObsRecord>> {
    let mut star_sorted: BTreeMap<(String, String, String, String), Vec<ObsRecord>> =
        BTreeMap::new();
    for (key, accum) in groups {
        if key.subtype == "*" && key.genotype == "*" && key.qq_field == "QUAL" {
            let max_count = match key.ty.as_str() {
                "SNP" => 12,
                "INDEL" => 40,
                _ => 4,
            };
            let mut obs = accum.obs.clone();
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
            let _ = max_count; // silence unused warning when not needed
        }
    }
    star_sorted
}

fn render_rows(
    groups: &BTreeMap<RowKey, GroupAccum>,
    star_sorted: &BTreeMap<(String, String, String, String), Vec<ObsRecord>>,
    row_filter: RowFilter<'_>,
    config: RenderConfig,
) -> Vec<String> {
    let mut out = Vec::new();
    // `accumulate` pre-seeds empty subtype buckets for types that are present,
    // because legacy reports zero-valued subtype baselines for an observed
    // type. It does not, however, emit a second family of rows for a wholly
    // absent type. Determine presence from actual observations before walking
    // the pre-seeded map.
    let active_types: HashSet<&str> = groups
        .iter()
        .filter(|(_, accum)| !accum.obs.is_empty())
        .map(|(key, _)| key.ty.as_str())
        .collect();
    for (key, accum) in groups {
        if !active_types.contains(key.ty.as_str()) {
            continue;
        }
        match row_filter {
            RowFilter::All => {}
            RowFilter::Locations { ty, filter } => {
                if key.ty != ty
                    || key.subtype != "*"
                    || key.subset != "*"
                    || key.genotype != "*"
                    || key.filter != filter
                {
                    continue;
                }
            }
        }
        let is_filter_tier = !is_aggregate_filter(&key.filter);
        let emitted_rows: Vec<EmittedRow> = if key.subtype == "*" {
            accum.emit_with_delta(config.delta)
        } else if let Some(shared) = star_sorted.get(&(
            key.ty.clone(),
            key.subtype.clone(),
            key.subset.clone(),
            key.filter.clone(),
        )) {
            accum.emit_with_shared_sort_and_delta(shared, &key.subtype, config.delta)
        } else {
            accum.emit_with_delta(config.delta)
        };
        for emitted in emitted_rows {
            if matches!(row_filter, RowFilter::Locations { .. }) && emitted.qq_str == "*" {
                // Legacy's Locations file drops the baseline row.
                continue;
            }
            if is_filter_tier && config.filter_counts_only && emitted.qq_str != "*" {
                // Per-Filter rows in legacy are baseline-only (QQ='*').
                // No numeric thresholds — skip the synthetic 0.0 and any
                // real numeric buckets that may have landed here.
                continue;
            }
            out.push(render_row(
                key,
                &emitted,
                config.subset_size,
                config.conf_size,
                is_filter_tier && config.filter_counts_only,
                config.ci_alpha,
            ));
        }
    }
    out
}

fn render_row(
    key: &RowKey,
    emitted: &EmittedRow,
    subset_size: usize,
    conf_size: usize,
    counts_only: bool,
    ci_alpha: f64,
) -> String {
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
        row.push(crate::report::precision_ratio(
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
    let (sz, conf) = subset_size_cells(&key.subset, &key.subtype, subset_size, conf_size);
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
        None if is_filter_tier || emitted.qq_str == "*" => append_stats(row, bucket, supports_titv),
        None => append_roc_stats(row, bucket, supports_titv),
        Some(avail) => {
            row.push(bucket.total.to_string());
            emit_substat_cell(row, bucket.ti, supports_titv && avail.ti, supports_titv);
            emit_substat_cell(row, bucket.tv, supports_titv && avail.tv, supports_titv);
            emit_substat_cell(row, bucket.het, avail.het, true);
            emit_substat_cell(row, bucket.homalt, avail.homalt, true);
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
        for (successes, trials) in [
            (counts.truth_tp.total, truth_total.total),
            (
                counts.query_tp.total,
                counts.query_tp.total + counts.query_fp.total,
            ),
            (counts.query_unk.total, query_total.total),
        ] {
            let (lower, upper) = jeffreys_interval(successes, trials, ci_alpha);
            row.push(format_ci(lower));
            row.push(format_ci(upper));
        }
    }

    row.join(",")
}

fn format_ci(value: f64) -> String {
    if value == 0.0 || value == 1.0 {
        format!("{value:.1}")
    } else {
        value.to_string()
    }
}

/// Modified Jeffreys interval used by legacy Tools/ci.py.
fn jeffreys_interval(x: usize, n: usize, alpha: f64) -> (f64, f64) {
    if n == 0 {
        return (0.0, 1.0);
    }
    let lower = if x == n {
        (alpha / 2.0).powf(1.0 / n as f64)
    } else if x <= 1 {
        0.0
    } else {
        inverse_regularized_beta(alpha / 2.0, x as f64 + 0.5, (n - x) as f64 + 0.5)
    };
    let upper = if x == 0 {
        1.0 - (alpha / 2.0).powf(1.0 / n as f64)
    } else if x >= n - 1 {
        1.0
    } else {
        inverse_regularized_beta(1.0 - alpha / 2.0, x as f64 + 0.5, (n - x) as f64 + 0.5)
    };
    (lower.max(0.0), upper.min(1.0))
}

fn inverse_regularized_beta(probability: f64, a: f64, b: f64) -> f64 {
    let mut low = 0.0;
    let mut high = 1.0;
    for _ in 0..80 {
        let middle = (low + high) / 2.0;
        if regularized_beta(middle, a, b) < probability {
            low = middle;
        } else {
            high = middle;
        }
    }
    (low + high) / 2.0
}

fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_continued_fraction(x, a, b) / a
    } else {
        1.0 - front * beta_continued_fraction(1.0 - x, b, a) / b
    }
}

fn beta_continued_fraction(x: f64, a: f64, b: f64) -> f64 {
    const EPSILON: f64 = 3.0e-14;
    const FLOOR: f64 = 1.0e-300;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FLOOR {
        d = FLOOR;
    }
    d = 1.0 / d;
    let mut result = d;
    for m in 1..=200 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let mut aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FLOOR {
            d = FLOOR;
        }
        c = 1.0 + aa / c;
        if c.abs() < FLOOR {
            c = FLOOR;
        }
        d = 1.0 / d;
        result *= d * c;

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FLOOR {
            d = FLOOR;
        }
        c = 1.0 + aa / c;
        if c.abs() < FLOOR {
            c = FLOOR;
        }
        d = 1.0 / d;
        let delta = d * c;
        result *= delta;
        if (delta - 1.0).abs() < EPSILON {
            break;
        }
    }
    result
}

fn log_gamma(value: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if value < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * value).sin().ln()
            - log_gamma(1.0 - value);
    }
    let z = value - 1.0;
    let mut sum = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        sum += coefficient / (z + index as f64);
    }
    let t = z + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + sum.ln()
}

fn emit_substat_cell(row: &mut Vec<String>, value: usize, kept: bool, supported: bool) {
    if !supported {
        row.push(".".to_string());
    } else if kept {
        row.push(format_count(value));
    } else {
        row.push(".".to_string());
    }
}

fn append_roc_stats(row: &mut Vec<String>, bucket: &CountsBucket, supports_titv: bool) {
    row.push(bucket.total.to_string());
    if supports_titv {
        row.push(format_count(bucket.ti));
        row.push(format_count(bucket.tv));
    } else {
        row.push(".".to_string());
        row.push(".".to_string());
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
    subtype: &str,
    subset_size: usize,
    conf_size: usize,
) -> (String, String) {
    // Derived from the four branches in report::write_extended. Kept in
    // lockstep — if write_extended changes its Subset.Size/IS_CONF.Size
    // convention, this function must follow.
    let is_base_subset = subset == "*";
    let is_base_subtype = subtype == "*";
    let size_cell = if is_base_subset {
        // Subset="*": always the raw subset_size integer, regardless of
        // subtype.
        subset_size.to_string()
    } else if subset == "TS_contained" {
        format_count(conf_size)
    } else {
        // TS_boundary
        format_count(subset_size)
    };
    // Legacy emits Subset.IS_CONF.Size as format_count(conf_size)
    // on every row when conf_size > 0 — independent of subtype or
    // subset selection. Empty cell only when conf_size==0.
    let _ = is_base_subset;
    let _ = is_base_subtype;
    let conf_cell = if conf_size > 0 {
        format_count(conf_size)
    } else {
        String::new()
    };
    (size_cell, conf_cell)
}

// ---------------------------------------------------------------------------
// Gzipped CSV output
// ---------------------------------------------------------------------------

fn write_gzip_csv(path: &Path, header: &str, rows: &[String]) -> Result<()> {
    let file = std::fs::File::create(path)
        .with_context(|| format!("failed to create {}", path.display()))?;
    let mut writer = GzEncoder::new(file, Compression::default());
    writeln!(writer, "{header}")?;
    for row in rows {
        writeln!(writer, "{row}")?;
    }
    writer.finish()?;
    Ok(())
}

fn write_optional_gzip_csv(path: &Path, header: &str, rows: &[String]) -> Result<()> {
    if rows.is_empty() {
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

#[cfg(test)]
mod tests {
    use super::*;

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

    #[allow(clippy::too_many_arguments)] // Keeps row fixtures legible at each call site.
    fn annotated(
        chrom: &str,
        pos: usize,
        qual: &str,
        truth_sample: &str,
        query_sample: &str,
        regions: &str,
        query_pass: bool,
        fp_class: Option<&'static str>,
    ) -> AnnotatedRow {
        let regions_tag = if regions.is_empty() {
            String::new()
        } else {
            format!(";Regions={regions}")
        };
        let line = format!(
            "{chrom}\t{pos}\t.\tA\tT\t{qual}\t.\tBS=1{regions_tag}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_sample}\t{query_sample}"
        );
        AnnotatedRow {
            sort_key: (chrom.to_string(), pos, 1, 0),
            line,
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
                "0/1:TP:gm:tv:SNP:het:500",
                "0/1:TP:gm:tv:SNP:het:0",
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
            .numeric_buckets
            .get("500.000000")
            .expect("missing 500.000000 bucket");
        assert_eq!(bucket_500.counts.truth_tp.total, 1);
        assert_eq!(bucket_500.counts.query_tp.total, 0);

        // Query-side bucket at qq=0.0 must hold the matching query TP.
        let bucket_0 = accum
            .numeric_buckets
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
                "0/1:TP:gm:tv:SNP:het:10",
                "0/1:TP:gm:tv:SNP:het:10",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "20",
                "1/1:TP:gm:ti:SNP:homalt:20",
                "1/1:TP:gm:ti:SNP:homalt:20",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                "0/1:TP:gm:ti:SNP:het:30",
                "0/1:TP:gm:ti:SNP:het:30",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                ".",
                "0/1:TP:gm:tv:SNP:het:.",
                "0/1:TP:gm:tv:SNP:het:.",
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
        // One SNP PASS row in region=TS_contained → contributions to (ALL,
        // PASS) × (*, TS_contained) × (Subtype="*" only, since SNP has no
        // non-* subtype) = 4 keys.
        let rows = vec![annotated(
            "chr1",
            100,
            "42",
            "0/1:TP:gm:tv:SNP:het:42",
            "0/1:TP:gm:tv:SNP:het:42",
            "CONF,TS_contained",
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
        let rendered = render_row(&key, &emitted, 100, 50, false, 0.0);
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
        // INDEL ti/tv cells are unsupported and use the legacy `.` marker.
        assert_eq!(cells[17], ".");
        assert_eq!(cells[18], ".");
        // TiTv_ratio: empty for INDEL.
        assert_eq!(cells[21], "");

        // Het/hom ratio: both zero → empty.
        let het_hom = het_hom_ratio(0, 0);
        assert_eq!(het_hom, "");
    }

    #[test]
    fn absent_variant_type_has_no_roc_rows_or_location_files() {
        let rows = vec![annotated(
            "chr1",
            100,
            "42",
            "0/1:TP:gm:i1_5:INDEL:het:42",
            "0/1:TP:gm:i1_5:INDEL:het:42",
            "",
            true,
            None,
        )];
        let dir = tempfile::tempdir().unwrap();
        let prefix = dir.path().join("result");
        write_roc_files(&prefix, &rows, 100, 0).unwrap();

        let all = crate::vcf::read_text(&suffixed_report_path(&prefix, "roc.all.csv.gz")).unwrap();
        assert!(all.lines().skip(1).all(|line| line.starts_with("INDEL,")));
        assert!(suffixed_report_path(&prefix, "roc.Locations.INDEL.csv.gz").exists());
        assert!(suffixed_report_path(&prefix, "roc.Locations.INDEL.PASS.csv.gz").exists());
        assert!(!suffixed_report_path(&prefix, "roc.Locations.SNP.csv.gz").exists());
        assert!(!suffixed_report_path(&prefix, "roc.Locations.SNP.PASS.csv.gz").exists());
    }

    #[test]
    fn pass_tier_demotes_tp_when_query_filtered() {
        // Truth/query TP but query_pass=false. Expected: ALL counts truth
        // as TP, PASS counts truth as FN and query contributes nothing.
        let rows = vec![annotated(
            "chr1",
            100,
            "10",
            "0/1:TP:gm:tv:SNP:het:10",
            "0/1:TP:gm:tv:SNP:het:10",
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
                "0/0:.:.:.:NOCALL:homref:.",
                "0/1:UNK:gm:tv:SNP:het:10",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "100",
                "0/0:.:.:.:NOCALL:homref:.",
                "0/1:.:gm:tv:SNP:het:100",
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
        assert_eq!(pass.obs.len(), 1);
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
                "0/1:TP:gm:tv:SNP:het:50",
                "0/1:TP:gm:tv:SNP:het:50",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "40",
                "0/1:FN:gm:tv:SNP:het:40",
                "0/1:FN:gm:tv:SNP:het:40",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                "0/1:TP:gm:tv:SNP:het:30",
                "0/1:TP:gm:tv:SNP:het:30",
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
                "0/1:TP:gm:tv:SNP:het:30",
                "0/1:TP:gm:tv:SNP:het:30",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "30.3",
                "0/1:TP:gm:tv:SNP:het:30.3",
                "0/1:TP:gm:tv:SNP:het:30.3",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30.6",
                "0/1:TP:gm:tv:SNP:het:30.6",
                "0/1:TP:gm:tv:SNP:het:30.6",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                "30.9",
                "0/1:TP:gm:tv:SNP:het:30.9",
                "0/1:TP:gm:tv:SNP:het:30.9",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                500,
                "31.2",
                "0/1:TP:gm:tv:SNP:het:31.2",
                "0/1:TP:gm:tv:SNP:het:31.2",
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
                "0/1:TP:gm:tv:SNP:het:10",
                "0/1:TP:gm:tv:SNP:het:10",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "20",
                "0/1:TP:gm:ti:SNP:het:20",
                "0/1:TP:gm:ti:SNP:het:20",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                "1/1:TP:gm:ti:SNP:homalt:30",
                "1/1:TP:gm:ti:SNP:homalt:30",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                "40",
                "0/1:TP:gm:tv:SNP:het:40",
                "0/1:TP:gm:tv:SNP:het:40",
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
            "0/1:TP:gm:tv:SNP:het:90",
            "0/1:TP:gm:tv:SNP:het:90",
            "",
            true,
            None,
        );
        first.line = first.line.replace("BS=1", "BS=1;SCORE=10.0");
        let mut second = annotated(
            "chr1",
            200,
            "80",
            "0/1:TP:gm:tv:SNP:het:80",
            "0/1:TP:gm:tv:SNP:het:80",
            "",
            true,
            None,
        );
        second.line = second.line.replace("BS=1", "BS=1;SCORE=10.4");
        let options = RocOptions {
            qq_field: "SCORE".to_string(),
            delta: 0.0,
            ..RocOptions::default()
        };
        let groups = accumulate_with_options(&[first, second], &options);
        let key = RowKey::new_with_qq_field("SNP", "*", "*", "ALL", "SCORE");
        let rows = groups[&key].emit_with_delta(options.delta);
        assert_eq!(
            rows.iter()
                .map(|row| row.qq_str.as_str())
                .collect::<Vec<_>>(),
            vec!["*", "10.000000", "10.400000"]
        );
        assert!(!groups[&key].numeric_buckets.contains_key("90.000000"));
    }

    #[test]
    fn ignored_filter_creates_selective_tier_and_renamed_filter_counts() {
        let mut row = annotated(
            "chr1",
            100,
            "20",
            "0/1:TP:gm:tv:SNP:het:20",
            "0/1:TP:gm:tv:SNP:het:20",
            "",
            false,
            None,
        );
        row.line = row.line.replacen("\t.\tBS=1", "\tLowQual\tBS=1", 1);
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
            "0/1:TP:gm:tv:SNP:het:20",
            "0/1:TP:gm:tv:SNP:het:20",
            "",
            false,
            None,
        );
        row.line = row.line.replacen("\t.\tBS=1", "\tLowQual\tBS=1", 1);
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
                conf_size: 0,
                delta: 0.5,
                ci_alpha: 0.0,
                filter_counts_only: options.roc_regions.contains("*"),
            },
        );
        assert!(rendered.iter().any(|line| {
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
    }
}
