//! Cohesive ROC legacy responsibility.

use super::*;

pub(super) struct MetricRows<'a> {
    pub(super) all: &'a [String],
    pub(super) snp: &'a [String],
    pub(super) snp_pass: &'a [String],
    pub(super) indel: &'a [String],
    pub(super) indel_pass: &'a [String],
    pub(super) selective: &'a [(String, Vec<String>, &'a str)],
}

pub(super) fn build_metric_indices(
    groups: &BTreeMap<RowKey, GroupAccum>,
    rows: MetricRows<'_>,
    delta: f64,
    output_rocs: bool,
) -> MetricIndices {
    let mut rocs = BTreeMap::<String, (&RowKey, &GroupAccum)>::new();
    let active_types = groups
        .iter()
        .filter(|(_, accum)| !accum.obs.is_empty())
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
        let counts_only = !output_rocs || !is_aggregate_filter(&key.filter);
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
    // qfy's `--no-roc` prevents C++ from placing threshold rows in the raw
    // table. The public CSV is compacted to the same baseline-only set before
    // JSON serialization, so its retained pandas indices must be computed
    // against that compact set rather than falling back to 0..N.
    let compact_all = (!output_rocs).then(|| {
        rows.all
            .iter()
            .filter(|line| line.split(',').nth(6) == Some("*"))
            .cloned()
            .collect::<Vec<_>>()
    });
    let metric_all = compact_all.as_deref().unwrap_or(rows.all);
    let all_indices = indices_for_lines(metric_all, &raw_positions);

    let mut available_tables = BTreeSet::from(["roc.all"]);
    for (id, lines) in [
        ("roc.Locations.SNP", rows.snp),
        ("roc.Locations.SNP.PASS", rows.snp_pass),
        ("roc.Locations.INDEL", rows.indel),
        ("roc.Locations.INDEL.PASS", rows.indel_pass),
    ] {
        if !lines.is_empty() {
            available_tables.insert(id);
        }
    }
    for (id, lines, _) in rows.selective {
        if !lines.is_empty() {
            available_tables.insert(id.as_str());
        }
    }
    let table_order = legacy_python_table_order(&raw, &available_tables);

    let mut tables = BTreeMap::new();
    tables.insert("roc.all".to_string(), all_indices.clone());
    tables.insert(
        "all.metrics".to_string(),
        metric_all
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
        metric_all
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
    MetricIndices {
        tables,
        table_order,
    }
}

/// Reproduce the iteration order of the pinned Python 2.7 dictionary used by
/// `qfy.py` to append `res` tables to metrics JSON. The dictionary keys and
/// hashes are fixed by the reference (`hash_randomization=0`); insertion order is
/// determined by the first matching row in the C++ unordered ROC table.
pub(super) fn legacy_python_table_order(
    raw: &[String],
    available_tables: &BTreeSet<&str>,
) -> Vec<String> {
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

pub(super) fn python27_dict_iteration_order(insertions: &[String]) -> Vec<String> {
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

pub(super) fn legacy_row_key(
    ty: &str,
    subtype: &str,
    filter: &str,
    subset: &str,
    qq: &str,
) -> String {
    format!("{ty}\t{subtype}\t*\t{filter}\t{subset}\t{qq}")
}

pub(super) fn csv_row_key(line: &str) -> String {
    let fields = line.split(',').take(7).collect::<Vec<_>>();
    if fields.len() < 7 {
        return String::new();
    }
    legacy_row_key(fields[0], fields[1], fields[3], fields[2], fields[6])
}

pub(super) fn indices_for_lines(lines: &[String], positions: &BTreeMap<&str, usize>) -> Vec<usize> {
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

pub(super) fn legacy_masked_levels(
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
pub(super) struct LegacyUnorderedRows {
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
        let bucket = (legacy_string_hash(&key) % self.bucket_count as u64) as usize;
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

pub(super) fn next_legacy_bucket(current: usize, required: usize) -> usize {
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
pub(super) fn legacy_string_hash(value: &str) -> u64 {
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
pub(super) struct LegacyRawTable {
    order: LegacyUnorderedRows,
    rows: BTreeMap<String, BTreeMap<String, String>>,
}

impl LegacyRawTable {
    fn set(&mut self, row: &str, column: &str, value: String, has_type: bool) {
        self.order.set(row.to_string(), has_type);
        self.rows
            .entry(row.to_string())
            .or_default()
            .insert(column.to_string(), value);
    }

    fn render(&self) -> String {
        use std::fmt::Write as _;

        let retained = self.order.retained_order();
        let columns = retained
            .iter()
            .filter_map(|key| self.rows.get(key))
            .flat_map(|row| row.keys().cloned())
            .collect::<BTreeSet<_>>();
        let mut output = String::new();
        writeln!(
            output,
            "{}",
            columns.iter().cloned().collect::<Vec<_>>().join("\t")
        )
        .expect("writing to a String cannot fail");
        for key in retained {
            let row = &self.rows[&key];
            writeln!(
                output,
                "{}",
                columns
                    .iter()
                    .map(|column| row.get(column).map(String::as_str).unwrap_or("."))
                    .collect::<Vec<_>>()
                    .join("\t")
            )
            .expect("writing to a String cannot fail");
        }
        output
    }
}

pub(super) fn build_legacy_roc_table(
    groups: &BTreeMap<RowKey, GroupAccum>,
    subset_size: usize,
    conf_size: usize,
    options: &RocOptions,
) -> String {
    // ROCOutput iterates a std::map keyed by its internal ROC name, not by
    // the final report axes. Recreate those names so insertion/rehash order
    // in LegacyUnorderedRows is byte-identical to libstdc++.
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
        let mut sorted_obs = accum.obs.clone();
        for subtype in subtypes {
            for genotype in ["het", "hetalt", "homalt", "*"] {
                let totals = legacy_totals(&accum.obs, subtype, genotype);
                add_legacy_level(
                    &mut table,
                    key,
                    subtype,
                    genotype,
                    "*",
                    &totals,
                    counts_only,
                    subset_size,
                    options.whole_reference_size.unwrap_or(subset_size),
                    conf_size,
                    &options.subset_sizes,
                    &options.subset_confidence_sizes,
                );
                if !counts_only && options.output_rocs {
                    introsort_libstdcpp(&mut sorted_obs);
                    for (qq, level) in legacy_levels(&sorted_obs, subtype, genotype, options.delta)
                    {
                        add_legacy_level(
                            &mut table,
                            key,
                            subtype,
                            genotype,
                            &qq,
                            &level,
                            false,
                            subset_size,
                            options.whole_reference_size.unwrap_or(subset_size),
                            conf_size,
                            &options.subset_sizes,
                            &options.subset_confidence_sizes,
                        );
                    }
                }
            }
        }
    }
    table.render()
}

pub(super) fn legacy_obs_matches(record: &ObsRecord, subtype: &str, genotype: &str) -> bool {
    let subtype_matches = match subtype {
        "*" => true,
        "ti" => record.ti_flag,
        "tv" => record.tv_flag,
        value => record.subtypes.iter().any(|candidate| candidate == value),
    };
    subtype_matches && (genotype == "*" || record.blt.as_deref() == Some(genotype))
}

pub(super) fn legacy_totals(obs: &[ObsRecord], subtype: &str, genotype: &str) -> Cumul {
    let mut totals = Cumul::default();
    for record in obs
        .iter()
        .filter(|record| legacy_obs_matches(record, subtype, genotype))
    {
        totals.add(&record.counts);
    }
    totals
}

pub(super) fn legacy_levels(
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
                truth_fn: add_buckets_total(&final_counts.truth_fn, &below.truth_tp),
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

#[allow(clippy::too_many_arguments)]
pub(super) fn add_legacy_level(
    table: &mut LegacyRawTable,
    key: &RowKey,
    subtype: &str,
    genotype: &str,
    qq: &str,
    counts: &Cumul,
    counts_only: bool,
    subset_size: usize,
    whole_reference_size: usize,
    conf_size: usize,
    subset_sizes: &BTreeMap<String, usize>,
    subset_confidence_sizes: &BTreeMap<String, usize>,
) {
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
        table.set(&row_key, "QQ", qq.to_string(), true);
        for (column, value) in [
            ("Type", key.ty.as_str()),
            ("Subtype", subtype),
            ("Genotype", genotype),
            ("Subset", key.subset.as_str()),
            ("Filter", key.filter.as_str()),
            ("QQ.Field", key.qq_field.as_str()),
        ] {
            table.set(&row_key, column, value.to_string(), true);
        }
        for (column, value) in legacy_primary_counts(counts, counts_only) {
            table.set(&row_key, column, value, true);
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
        table.set(&row_key, "Subset.Size", legacy_count_string(&size), true);
        table.set(
            &row_key,
            "Subset.IS_CONF.Size",
            legacy_count_string(&conf),
            true,
        );
        table.set(&row_key, "Subset.Level", "0.000000".to_string(), true);
    } else if genotype == "*" && matches!(subtype, "ti" | "tv") {
        let aggregate = legacy_row_key(&key.ty, "*", &key.filter, &key.subset, qq);
        for (metric, count) in legacy_count_buckets(counts, counts_only) {
            table.set(
                &aggregate,
                &format!("{metric}.{subtype}"),
                legacy_usize(count),
                false,
            );
        }
    } else if genotype != "*" && !matches!(subtype, "ti" | "tv") {
        for (metric, count) in legacy_count_buckets(counts, counts_only) {
            table.set(
                &row_key,
                &format!("{metric}.{genotype}"),
                legacy_usize(count),
                false,
            );
        }
    }
}

pub(super) fn legacy_count_buckets(
    counts: &Cumul,
    counts_only: bool,
) -> Vec<(&'static str, usize)> {
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

pub(super) fn legacy_primary_counts(
    counts: &Cumul,
    counts_only: bool,
) -> Vec<(&'static str, String)> {
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

pub(super) fn legacy_usize(value: usize) -> String {
    format!("{:.6}", value as f64)
}

pub(super) fn legacy_count_string(value: &str) -> String {
    value
        .parse::<f64>()
        .map(legacy_f64)
        .unwrap_or_else(|_| "0.000000".to_string())
}

pub(super) fn legacy_f64(value: f64) -> String {
    if value.is_nan() {
        "-nan".to_string()
    } else {
        format!("{value:.6}")
    }
}

// ---------------------------------------------------------------------------
// Grouping
// ---------------------------------------------------------------------------
