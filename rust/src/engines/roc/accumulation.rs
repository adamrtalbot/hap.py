//! Cohesive ROC accumulation responsibility.

use super::*;

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(super) struct RowKey {
    // Legacy sort order: (Type, Subtype, Subset, Filter, Genotype, QQ.Field)
    // all ascending. Derive Ord via field order.
    pub(super) ty: String,
    pub(super) subtype: String,
    pub(super) subset: String,
    pub(super) filter: String,
    pub(super) genotype: String,
    pub(super) qq_field: String,
}

impl RowKey {
    #[cfg(test)]
    pub(super) fn new(ty: &str, subtype: &str, subset: &str, filter: &str) -> Self {
        Self::new_with_qq_field(ty, subtype, subset, filter, "QUAL")
    }

    pub(super) fn new_with_qq_field(
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
pub(super) struct Cumul {
    pub(super) truth_tp: CountsBucket,
    pub(super) truth_fn: CountsBucket,
    pub(super) query_tp: CountsBucket,
    pub(super) query_fp: CountsBucket,
    pub(super) query_unk: CountsBucket,
    pub(super) fp_gt: usize,
    pub(super) fp_al: usize,
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
pub(super) struct ObsRecord {
    pub(super) level: f64,
    pub(super) counts: Cumul,
    pub(super) subtypes: Vec<String>,
    /// True if the original BI tag included `ti` — needed to mirror
    /// legacy's `OBS_FLAG_TI` mask which is set per-record from the BI
    /// string, *not* from non-zero count fields. Filter-failed query
    /// phantoms have all-zero counts but still carry their flag bits.
    pub(super) ti_flag: bool,
    pub(super) tv_flag: bool,
    pub(super) blt: Option<String>,
}

impl Cumul {
    pub(super) fn add(&mut self, other: &Cumul) {
        add_bucket(&mut self.truth_tp, &other.truth_tp);
        add_bucket(&mut self.truth_fn, &other.truth_fn);
        add_bucket(&mut self.query_tp, &other.query_tp);
        add_bucket(&mut self.query_fp, &other.query_fp);
        add_bucket(&mut self.query_unk, &other.query_unk);
        self.fp_gt += other.fp_gt;
        self.fp_al += other.fp_al;
    }

    pub(super) fn truth_total(&self) -> CountsBucket {
        sum_buckets(&self.truth_tp, &self.truth_fn)
    }

    pub(super) fn query_total(&self) -> CountsBucket {
        let a = sum_buckets(&self.query_tp, &self.query_fp);
        sum_buckets(&a, &self.query_unk)
    }
}

pub(super) fn add_bucket(dst: &mut CountsBucket, src: &CountsBucket) {
    dst.total += src.total;
    dst.ti += src.ti;
    dst.tv += src.tv;
    dst.het += src.het;
    dst.homalt += src.homalt;
}

pub(super) fn sum_buckets(a: &CountsBucket, b: &CountsBucket) -> CountsBucket {
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
pub(super) struct GroupAccum {
    pub(super) baseline: Cumul,
    pub(super) obs: Vec<ObsRecord>,
    pub(super) numeric_buckets: BTreeMap<String, NumericBucket>,
}

pub(super) struct NumericBucket {
    /// Numeric QQ value used for DESC cumulation ordering. All contributions
    /// that round to the same `{:.6}` string share a bucket and are merged
    /// (matches what the legacy pandas pipeline would do when it groups by
    /// the formatted QQ string).
    pub(super) qq: f64,
    pub(super) counts: Cumul,
    /// Substat-flag presence at this bucket. True if any obs (real or phantom)
    /// at this level was tagged with the corresponding bi. Used by
    /// substat_kept_levels to decide whether the substat sweep keeps this
    /// level. Mirrors legacy's BlockQuantify::observe which sets OBS_FLAG_TI/TV
    /// on the FN2 phantom obs from filter-failed query.TP records too.
    pub(super) has_ti: bool,
    pub(super) has_tv: bool,
}

/// Substat availability at a given numeric QQ row. Each Option carries the
/// cumulative count for that substat if the independent per-substat
/// roc-delta sweep kept this level; `None` means the substat cell renders
/// as the legacy na_rep (`.` for count cells, `""` for ratio cells).
#[derive(Clone, Debug, Default)]
pub(super) struct SubstatAvail {
    pub(super) ti: bool,
    pub(super) tv: bool,
    pub(super) het: bool,
    pub(super) homalt: bool,
}

/// One row emitted by a `GroupAccum`. For the baseline row (qq_str = "*"),
/// `substats` is `None` and every substat cell is taken from `cum`. For
/// numeric QQ rows, `substats` marks which substat cells the independent
/// sweeps kept at this level; cells where the substat sweep didn't keep
/// render as na_rep.
#[derive(Clone, Debug)]
pub(super) struct EmittedRow {
    pub(super) qq_str: String,
    pub(super) cum: Cumul,
    pub(super) substats: Option<SubstatAvail>,
}

impl GroupAccum {
    pub(super) fn add(
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
    pub(super) fn emit(&self) -> Vec<EmittedRow> {
        self.emit_with_delta(0.5)
    }

    pub(super) fn emit_with_delta(&self, delta: f64) -> Vec<EmittedRow> {
        self.emit_internal(None, None, delta)
    }

    /// Emit using a pre-sorted obs vector (typically from the `(*)`
    /// sibling in the same `(Type, Subset, Filter)` group). When
    /// `my_subtype` is non-empty and not `"*"`, the shared obs are
    /// filtered to entries whose `subtypes` field contains `my_subtype`
    /// — matching legacy's "single per-Type sort, walk per subtype
    /// with flag mask" semantics. This is what reproduces legacy's
    /// per-cluster cluster-ordering at non-(*)-subtype slices.
    pub(super) fn emit_with_shared_sort_and_delta(
        &self,
        shared_sorted: &[ObsRecord],
        my_subtype: &str,
        delta: f64,
    ) -> Vec<EmittedRow> {
        self.emit_internal(Some(shared_sorted), Some(my_subtype), delta)
    }

    pub(super) fn emit_internal(
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

pub(super) fn substat_kept_levels(
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

pub(super) fn bucket_has_ti(b: &NumericBucket) -> bool {
    b.has_ti
}

pub(super) fn bucket_has_tv(b: &NumericBucket) -> bool {
    b.has_tv
}

pub(super) fn sub_buckets(a: &CountsBucket, b: &CountsBucket) -> CountsBucket {
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

pub(super) fn add_buckets_total(a: &CountsBucket, b: &CountsBucket) -> CountsBucket {
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

pub(super) const ROC_SORT_THRESHOLD: usize = 16;

#[inline]
pub(super) fn lg_floor(n: usize) -> usize {
    if n <= 1 {
        0
    } else {
        (usize::BITS - 1 - n.leading_zeros()) as usize
    }
}

pub(super) fn introsort_libstdcpp(arr: &mut [ObsRecord]) {
    let n = arr.len();
    if n > 1 {
        let depth_limit = lg_floor(n) * 2;
        introsort_loop(arr, 0, n, depth_limit);
        final_insertion_sort(arr);
    }
}

pub(super) fn introsort_loop(
    arr: &mut [ObsRecord],
    first: usize,
    mut last: usize,
    mut depth_limit: usize,
) {
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

pub(super) fn unguarded_partition_pivot(arr: &mut [ObsRecord], first: usize, last: usize) -> usize {
    let mid = first + (last - first) / 2;
    move_median_to_first(arr, first, first + 1, mid, last - 1);
    unguarded_partition(arr, first + 1, last, first)
}

pub(super) fn move_median_to_first(
    arr: &mut [ObsRecord],
    result: usize,
    a: usize,
    b: usize,
    c: usize,
) {
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

pub(super) fn unguarded_partition(
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

pub(super) fn final_insertion_sort(arr: &mut [ObsRecord]) {
    let n = arr.len();
    if n > ROC_SORT_THRESHOLD {
        insertion_sort_range(arr, 0, ROC_SORT_THRESHOLD);
        unguarded_insertion_sort_range(arr, ROC_SORT_THRESHOLD, n);
    } else {
        insertion_sort_range(arr, 0, n);
    }
}

pub(super) fn insertion_sort_range(arr: &mut [ObsRecord], first: usize, last: usize) {
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

pub(super) fn unguarded_insertion_sort_range(arr: &mut [ObsRecord], first: usize, last: usize) {
    for i in first..last {
        let mut j = i;
        while j > 0 && arr[j].level < arr[j - 1].level {
            arr.swap(j, j - 1);
            j -= 1;
        }
    }
}

pub(super) fn heapsort_range(arr: &mut [ObsRecord], first: usize, last: usize) {
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

pub(super) fn sift_down(arr: &mut [ObsRecord], start: usize, first: usize, last: usize) {
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
pub(super) fn accumulate(rows: &[AnnotatedRow]) -> BTreeMap<RowKey, GroupAccum> {
    accumulate_impl(rows, &RocOptions::default())
}

pub(super) fn accumulate_impl(
    rows: &[AnnotatedRow],
    options: &RocOptions,
) -> BTreeMap<RowKey, GroupAccum> {
    let mut groups: BTreeMap<RowKey, GroupAccum> = BTreeMap::new();

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

pub(super) fn accumulate_with_options(
    rows: &[AnnotatedRow],
    options: &RocOptions,
) -> BTreeMap<RowKey, GroupAccum> {
    accumulate_impl(rows, options)
}

pub(super) fn roc_header(ci_alpha: f64) -> String {
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

pub(super) fn is_aggregate_filter(filter: &str) -> bool {
    matches!(filter, "ALL" | "PASS" | "SEL")
}

// ---------------------------------------------------------------------------
// Per-row contribution emission
// ---------------------------------------------------------------------------

/// Lightweight VCF sample parser (reads the handful of FORMAT keys roc.rs
/// cares about). Mirrors `compare::SampleView` without introducing a
/// cross-module dep on the private struct.
pub(super) struct Sample<'a> {
    pub(super) format_keys: &'a [&'a str],
    pub(super) parts: &'a [&'a str],
    pub(super) gt: Option<&'a str>,
    pub(super) bd: Option<&'a str>,
    pub(super) bi: Option<&'a str>,
    pub(super) bvt: Option<&'a str>,
    pub(super) blt: Option<&'a str>,
    pub(super) qq: Option<&'a str>,
}

impl<'a> Sample<'a> {
    pub(super) fn new(format_keys: &'a [&'a str], parts: &'a [&'a str]) -> Self {
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

    pub(super) fn variant_type(&self) -> Option<&'a str> {
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
    pub(super) fn roc_qq(&self) -> Option<f64> {
        self.qq.and_then(|raw| raw.parse::<f64>().ok())
    }

    pub(super) fn roc_value(&self, field: &str, record_qual: &str, info: &str) -> Option<f64> {
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

pub(super) fn parse_roc_number(raw: &str) -> Option<f64> {
    raw.split(',')
        .next()
        .filter(|value| !matches!(*value, "" | "."))
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite())
}
