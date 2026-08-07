//! Native implementation of the RTG `vcfeval -m ga4gh` behavior used by
//! hap.py.
//!
//! The path search follows the RTG 3.12.1 model: each side may include or
//! exclude a variant, unphased diploid orientations are explored, and the
//! reconciled path containing the greatest total number of truth and query
//! variants wins. A second haploid pass identifies allele-only matches and a
//! final proximity pass applies GA4GH loose-match annotations.

use crate::domain::RawVcfRecord;
use crate::engines::scmp;
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::time::{SystemTime, UNIX_EPOCH};

const RTG_MAX_PATHS: usize = 50_000;
const PATH_PADDING: usize = 32;

#[derive(Clone, Copy, Debug)]
pub(crate) struct Options<'a> {
    pub roc_field: &'a str,
    pub loose_match_distance: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Side {
    Truth,
    Query,
}

#[derive(Clone, Debug)]
struct Call {
    record_index: usize,
    side: Side,
    chrom: String,
    pos: usize,
    ref_allele: String,
    alts: Vec<String>,
    gt: Vec<Option<usize>>,
    skipped: bool,
}

impl Call {
    fn end(&self) -> usize {
        self.pos + self.ref_allele.len().saturating_sub(1)
    }

    fn non_ref_alleles(&self) -> impl Iterator<Item = usize> + '_ {
        self.gt
            .iter()
            .flatten()
            .copied()
            .filter(|allele| *allele > 0)
    }
}

#[derive(Clone, Debug)]
struct Cluster {
    chrom: String,
    start: usize,
    end: usize,
    calls: Vec<usize>,
}

#[derive(Clone, Debug)]
struct Edit {
    pos: usize,
    reference: String,
    alternate: String,
}

#[derive(Clone, Debug, Default)]
struct SearchState {
    left: Vec<Edit>,
    right: Vec<Edit>,
    included: BTreeSet<usize>,
}

#[derive(Clone, Debug, Default)]
struct Candidate {
    included: BTreeSet<usize>,
    count: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct Verdict {
    gt_match: bool,
    allele_match: bool,
    loose_match: bool,
    skipped: bool,
    sync_start: usize,
}

pub(crate) fn compare_records(
    truth_headers: &[String],
    truth_records: &[RawVcfRecord],
    query_headers: &[String],
    query_records: &[RawVcfRecord],
    references: &BTreeMap<String, String>,
    options: Options<'_>,
) -> Result<(Vec<String>, Vec<RawVcfRecord>)> {
    let (mut merged, query_sources) =
        merge_records_rtg(truth_headers, truth_records, query_headers, query_records)?;
    let calls = extract_calls(&merged.records)?;
    let clusters = build_clusters(&calls);
    let mut verdicts = vec![Verdict::default(); calls.len()];

    for cluster in &clusters {
        let reference = references
            .get(&cluster.chrom)
            .ok_or_else(|| anyhow::anyhow!("reference has no contig {}", cluster.chrom))?;
        compare_cluster(cluster, &calls, reference, &mut verdicts)?;
    }
    apply_allele_matches(&calls, &references, &mut verdicts)?;
    apply_loose_matches(&calls, options.loose_match_distance, &mut verdicts);

    let qq = query_scores(
        &merged.records,
        &query_records,
        &query_sources,
        options.roc_field,
    )?;
    annotate_records(&mut merged.records, &calls, &verdicts, &qq);
    let headers = output_headers(query_headers, references, options.loose_match_distance);
    Ok((headers, merged.records))
}

/// RTG combines truth and query records only when their start and REF agree.
/// Records merely sharing a coordinate remain separate (query sorts first),
/// which is essential for decomposed-vs-MNP path matches.
fn merge_records_rtg(
    truth_headers: &[String],
    truth: &[RawVcfRecord],
    query_headers: &[String],
    query: &[RawVcfRecord],
) -> Result<(scmp::MergedVcf, Vec<Option<usize>>)> {
    let mut used_query = vec![false; query.len()];
    let mut tagged = Vec::with_capacity(truth.len() + query.len());
    for truth_record in truth {
        let match_index = query.iter().enumerate().position(|(index, query_record)| {
            !used_query[index]
                && query_record.chrom == truth_record.chrom
                && query_record.pos == truth_record.pos
                && query_record.ref_allele == truth_record.ref_allele
        });
        if let Some(index) = match_index {
            used_query[index] = true;
            let merged = scmp::merge_two_sample_records(
                truth_headers,
                std::slice::from_ref(truth_record),
                query_headers,
                std::slice::from_ref(&query[index]),
            )?;
            tagged.extend(
                merged
                    .records
                    .into_iter()
                    .map(|record| (record, 0usize, Some(index))),
            );
        } else {
            tagged.push((simple_two_sample_record(truth_record, Side::Truth), 1, None));
        }
    }
    for (index, query_record) in query.iter().enumerate() {
        if !used_query[index] {
            tagged.push((
                simple_two_sample_record(query_record, Side::Query),
                0,
                Some(index),
            ));
        }
    }
    tagged.sort_by(|(left, left_side, _), (right, right_side, _)| {
        left.chrom
            .cmp(&right.chrom)
            .then(left.pos.cmp(&right.pos))
            .then(left_side.cmp(right_side))
            .then(left.ref_allele.len().cmp(&right.ref_allele.len()))
    });
    let (records, query_sources) = tagged
        .into_iter()
        .map(|(record, _, query_source)| (record, query_source))
        .unzip();
    Ok((
        scmp::MergedVcf {
            headers: Vec::new(),
            records,
        },
        query_sources,
    ))
}

fn simple_two_sample_record(record: &RawVcfRecord, side: Side) -> RawVcfRecord {
    let gt = record
        .sample_map(0)
        .get("GT")
        .cloned()
        .unwrap_or_else(|| ".".to_string());
    RawVcfRecord {
        chrom: record.chrom.clone(),
        pos: record.pos,
        id: record.id.clone(),
        ref_allele: record.ref_allele.clone(),
        alt_allele: record.alt_allele.clone(),
        qual: record.qual.clone(),
        filter: if side == Side::Query {
            record.filter.clone()
        } else {
            ".".to_string()
        },
        info: record.info.clone(),
        format: Some("GT".to_string()),
        samples: match side {
            Side::Truth => vec![gt, ".".to_string()],
            Side::Query => vec![".".to_string(), gt],
        },
    }
}

fn extract_calls(records: &[RawVcfRecord]) -> Result<Vec<Call>> {
    let mut calls = Vec::new();
    for (record_index, record) in records.iter().enumerate() {
        let alts = record
            .alt_allele
            .split(',')
            .map(str::to_string)
            .collect::<Vec<_>>();
        for (sample_index, side) in [(0, Side::Truth), (1, Side::Query)] {
            let sample = record.sample_map(sample_index);
            let Some(gt) = sample.get("GT") else {
                continue;
            };
            let alleles = gt
                .split(['/', '|'])
                .map(|value| value.parse::<usize>().ok())
                .collect::<Vec<_>>();
            if !alleles.iter().flatten().any(|allele| *allele > 0) {
                continue;
            }
            let invalid_index = alleles.iter().flatten().any(|allele| *allele > alts.len());
            if invalid_index {
                bail!(
                    "GT allele index exceeds ALT count at {}:{}",
                    record.chrom,
                    record.pos
                );
            }
            let skipped = alleles.len() > 2
                || alleles.iter().any(Option::is_none)
                || alleles
                    .iter()
                    .flatten()
                    .filter(|allele| **allele > 0)
                    .any(|allele| {
                        let alt = &alts[*allele - 1];
                        alt == "."
                            || alt == "*"
                            || alt.starts_with('<')
                            || !alt.bytes().all(|base| base.is_ascii_alphabetic())
                    });
            calls.push(Call {
                record_index,
                side,
                chrom: record.chrom.clone(),
                pos: record.pos,
                ref_allele: record.ref_allele.clone(),
                alts: alts.clone(),
                gt: alleles,
                skipped,
            });
        }
    }
    Ok(calls)
}

fn build_clusters(calls: &[Call]) -> Vec<Cluster> {
    let mut order = (0..calls.len()).collect::<Vec<_>>();
    order.sort_by(|left, right| {
        calls[*left]
            .chrom
            .cmp(&calls[*right].chrom)
            .then(calls[*left].pos.cmp(&calls[*right].pos))
            .then(calls[*left].end().cmp(&calls[*right].end()))
    });
    let mut clusters = Vec::new();
    let mut current: Option<Cluster> = None;
    for call_id in order {
        let call = &calls[call_id];
        match current.as_mut() {
            Some(cluster)
                if cluster.chrom == call.chrom && call.pos <= cluster.end.saturating_add(1) =>
            {
                cluster.end = cluster.end.max(call.end());
                cluster.calls.push(call_id);
            }
            _ => {
                if let Some(cluster) = current.take() {
                    clusters.push(cluster);
                }
                current = Some(Cluster {
                    chrom: call.chrom.clone(),
                    start: call.pos,
                    end: call.end(),
                    calls: vec![call_id],
                });
            }
        }
    }
    if let Some(cluster) = current {
        clusters.push(cluster);
    }
    clusters
}

fn compare_cluster(
    cluster: &Cluster,
    calls: &[Call],
    reference: &str,
    verdicts: &mut [Verdict],
) -> Result<()> {
    for call_id in &cluster.calls {
        if calls[*call_id].skipped {
            verdicts[*call_id].skipped = true;
        }
    }
    let truth = cluster
        .calls
        .iter()
        .copied()
        .filter(|call_id| calls[*call_id].side == Side::Truth && !calls[*call_id].skipped)
        .collect::<Vec<_>>();
    let query = cluster
        .calls
        .iter()
        .copied()
        .filter(|call_id| calls[*call_id].side == Side::Query && !calls[*call_id].skipped)
        .collect::<Vec<_>>();
    if truth.is_empty() && query.is_empty() {
        return Ok(());
    }

    let region_start = cluster.start.saturating_sub(PATH_PADDING).max(1);
    let region_end = cluster
        .end
        .saturating_add(PATH_PADDING)
        .min(reference.len());
    let truth_paths = enumerate_paths(&truth, calls, reference, region_start, region_end)?;
    let query_paths = enumerate_paths(&query, calls, reference, region_start, region_end)?;

    let mut best_score = 0usize;
    let mut best_delta = usize::MAX;
    let mut best_truth = BTreeSet::new();
    let mut best_query = BTreeSet::new();
    for (signature, truth_candidate) in &truth_paths {
        let Some(query_candidate) = query_paths.get(signature) else {
            continue;
        };
        let score = truth_candidate.count + query_candidate.count;
        let delta = truth_candidate.count.abs_diff(query_candidate.count);
        if score > best_score || (score == best_score && delta < best_delta) {
            best_score = score;
            best_delta = delta;
            best_truth = truth_candidate.included.clone();
            best_query = query_candidate.included.clone();
        }
    }
    let truth_sync = truth.iter().map(|call_id| calls[*call_id].pos).min();
    let query_sync = query.iter().map(|call_id| calls[*call_id].pos).min();
    let sync_start = match (truth_sync, query_sync) {
        (Some(truth), Some(query)) => truth.max(query),
        (Some(start), None) | (None, Some(start)) => start,
        (None, None) => cluster.start,
    };
    for call_id in truth.into_iter().chain(query) {
        verdicts[call_id].gt_match = best_truth.contains(&call_id) || best_query.contains(&call_id);
        verdicts[call_id].sync_start = sync_start;
    }
    Ok(())
}

fn enumerate_paths(
    call_ids: &[usize],
    calls: &[Call],
    reference: &str,
    region_start: usize,
    region_end: usize,
) -> Result<BTreeMap<String, Candidate>> {
    let mut states = vec![SearchState::default()];
    for call_id in call_ids {
        let call = &calls[*call_id];
        let choices = call_choices(call);
        let projected = checked_state_count(states.len(), choices.len())?;
        let mut next = Vec::with_capacity(projected);
        for state in &states {
            next.push(state.clone());
            for (left, right) in &choices {
                let mut included = state.clone();
                if let Some(edit) = allele_edit(call, *left) {
                    included.left.push(edit);
                }
                if let Some(edit) = allele_edit(call, *right) {
                    included.right.push(edit);
                }
                included.included.insert(*call_id);
                next.push(included);
            }
        }
        states = next;
    }

    let mut paths: BTreeMap<String, Candidate> = BTreeMap::new();
    for state in states {
        let Some(signature) = render_signature(
            reference,
            region_start,
            region_end,
            &state.left,
            &state.right,
        )?
        else {
            continue;
        };
        let count = state.included.len();
        match paths.get_mut(&signature) {
            Some(candidate) if candidate.count == count => {}
            Some(candidate) if candidate.count > count => {}
            _ => {
                paths.insert(
                    signature,
                    Candidate {
                        included: state.included,
                        count,
                    },
                );
            }
        }
    }
    Ok(paths)
}

fn checked_state_count(state_count: usize, choice_count: usize) -> Result<usize> {
    let branch_count = choice_count
        .checked_add(1)
        .context("vcfeval path branch count overflow")?;
    let projected = state_count
        .checked_mul(branch_count)
        .context("vcfeval path count overflow")?;
    if projected > RTG_MAX_PATHS {
        bail!("vcfeval cluster requires {projected} paths, exceeding the limit of {RTG_MAX_PATHS}");
    }
    Ok(projected)
}

fn call_choices(call: &Call) -> Vec<(usize, usize)> {
    let called = call.gt.iter().flatten().copied().collect::<Vec<_>>();
    let (left, right) = match called.as_slice() {
        [allele] => (*allele, *allele),
        [left, right] => (*left, *right),
        _ => return Vec::new(),
    };
    let mut choices = vec![(left, right)];
    if left != right {
        choices.push((right, left));
    }
    choices
}

fn allele_edit(call: &Call, allele: usize) -> Option<Edit> {
    (allele > 0).then(|| Edit {
        pos: call.pos,
        reference: call.ref_allele.clone(),
        alternate: call.alts[allele - 1].clone(),
    })
}

fn render_signature(
    reference: &str,
    region_start: usize,
    region_end: usize,
    left: &[Edit],
    right: &[Edit],
) -> Result<Option<String>> {
    let Some(left) = apply_edits(reference, region_start, region_end, left)? else {
        return Ok(None);
    };
    let Some(right) = apply_edits(reference, region_start, region_end, right)? else {
        return Ok(None);
    };
    let (first, second) = if left <= right {
        (left, right)
    } else {
        (right, left)
    };
    Ok(Some(format!("{first}\u{0}{second}")))
}

fn apply_edits(
    reference: &str,
    region_start: usize,
    region_end: usize,
    edits: &[Edit],
) -> Result<Option<String>> {
    if region_start == 0 || region_start > region_end || region_end > reference.len() {
        bail!("invalid vcfeval reference slice {region_start}-{region_end}");
    }
    let mut edits = edits.to_vec();
    edits.sort_by(|left, right| {
        left.pos
            .cmp(&right.pos)
            .then(left.reference.len().cmp(&right.reference.len()))
    });
    let mut cursor = region_start;
    let mut output = String::new();
    for edit in edits {
        let edit_end = edit.pos + edit.reference.len().saturating_sub(1);
        if edit.pos < cursor || edit.pos < region_start || edit_end > region_end {
            return Ok(None);
        }
        let observed = &reference.as_bytes()[edit.pos - 1..edit_end];
        if !observed.eq_ignore_ascii_case(edit.reference.as_bytes()) {
            return Ok(None);
        }
        output.push_str(
            reference
                .get(cursor - 1..edit.pos - 1)
                .context("reference positions are not UTF-8 boundaries")?,
        );
        output.push_str(&edit.alternate.to_ascii_uppercase());
        cursor = edit_end + 1;
    }
    output.push_str(
        reference
            .get(cursor - 1..region_end)
            .context("reference positions are not UTF-8 boundaries")?,
    );
    Ok(Some(output.to_ascii_uppercase()))
}

fn apply_allele_matches(
    calls: &[Call],
    references: &BTreeMap<String, String>,
    verdicts: &mut [Verdict],
) -> Result<()> {
    for (truth_id, truth) in calls.iter().enumerate() {
        if truth.side != Side::Truth || truth.skipped || verdicts[truth_id].gt_match {
            continue;
        }
        for (query_id, query) in calls.iter().enumerate() {
            if query.side != Side::Query
                || query.skipped
                || verdicts[query_id].gt_match
                || truth.chrom != query.chrom
            {
                continue;
            }
            let reference = references
                .get(&truth.chrom)
                .ok_or_else(|| anyhow::anyhow!("reference has no contig {}", truth.chrom))?;
            if calls_share_allele(truth, query, reference)? {
                verdicts[truth_id].allele_match = true;
                verdicts[query_id].allele_match = true;
                let sync = truth.pos.min(query.pos);
                verdicts[truth_id].sync_start = sync;
                verdicts[query_id].sync_start = sync;
            }
        }
    }
    Ok(())
}

fn calls_share_allele(truth: &Call, query: &Call, reference: &str) -> Result<bool> {
    let start = truth.pos.min(query.pos).saturating_sub(PATH_PADDING).max(1);
    let end = truth
        .end()
        .max(query.end())
        .saturating_add(PATH_PADDING)
        .min(reference.len());
    for truth_allele in truth.non_ref_alleles() {
        let truth_edit = allele_edit(truth, truth_allele).expect("non-reference allele has edit");
        let Some(truth_sequence) = apply_edits(reference, start, end, &[truth_edit])? else {
            continue;
        };
        for query_allele in query.non_ref_alleles() {
            let query_edit =
                allele_edit(query, query_allele).expect("non-reference allele has edit");
            if apply_edits(reference, start, end, &[query_edit])?.as_ref() == Some(&truth_sequence)
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn apply_loose_matches(calls: &[Call], distance: usize, verdicts: &mut [Verdict]) {
    let distance = distance as i128;
    for call_id in 0..calls.len() {
        if verdicts[call_id].gt_match || verdicts[call_id].allele_match || verdicts[call_id].skipped
        {
            continue;
        }
        let call = &calls[call_id];
        let opposite = match call.side {
            Side::Truth => Side::Query,
            Side::Query => Side::Truth,
        };
        verdicts[call_id].loose_match = calls.iter().enumerate().any(|(other_id, other)| {
            if other.side != opposite || other.chrom != call.chrom || verdicts[other_id].skipped {
                return false;
            }
            let start = call.pos as i128 - 1;
            let end = call.end() as i128;
            let other_start = other.pos as i128 - 1;
            let other_end = other.end() as i128;
            start < other_end + distance && other_start - distance < end
        });
    }
}

fn query_scores(
    merged: &[RawVcfRecord],
    query: &[RawVcfRecord],
    query_sources: &[Option<usize>],
    field: &str,
) -> Result<BTreeMap<usize, String>> {
    if merged.len() != query_sources.len() {
        bail!("vcfeval query provenance does not match merged records");
    }
    let mut scores = BTreeMap::new();
    for (record_index, record) in merged.iter().enumerate() {
        if !record.sample_map(1).get("GT").is_some_and(|gt| {
            gt.split(['/', '|'])
                .filter_map(|value| value.parse::<usize>().ok())
                .any(|allele| allele > 0)
        }) {
            continue;
        }
        let Some(source_index) = query_sources[record_index] else {
            continue;
        };
        let source = query.get(source_index).with_context(|| {
            format!("vcfeval query provenance index {source_index} is out of bounds")
        })?;
        let value = if field == "QUAL" {
            Some(source.qual.clone())
        } else if let Some(key) = field.strip_prefix("INFO.") {
            info_value(&source.info, key).map(str::to_string)
        } else {
            let key = field.strip_prefix("FORMAT.").unwrap_or(field);
            source.sample_map(0).get(key).cloned()
        };
        if let Some(value) = value
            && value != "."
            && !value.is_empty()
        {
            let score = value
                .split(',')
                .next()
                .unwrap_or(&value)
                .parse::<f64>()
                .with_context(|| format!("invalid vcfeval ROC score '{value}'"))?;
            scores.insert(record_index, format!("{score:.1}"));
        }
    }
    Ok(scores)
}

fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';').find_map(|entry| {
        let (candidate, value) = entry.split_once('=')?;
        (candidate == key).then_some(value)
    })
}

fn annotate_records(
    records: &mut Vec<RawVcfRecord>,
    calls: &[Call],
    verdicts: &[Verdict],
    qq: &BTreeMap<usize, String>,
) {
    let mut by_record: BTreeMap<usize, [Option<usize>; 2]> = BTreeMap::new();
    for (call_id, call) in calls.iter().enumerate() {
        by_record.entry(call.record_index).or_insert([None, None])[match call.side {
            Side::Truth => 0,
            Side::Query => 1,
        }] = Some(call_id);
    }
    let mut output = Vec::with_capacity(records.len());
    for (record_index, mut record) in records.drain(..).enumerate() {
        let Some(call_ids) = by_record.get(&record_index) else {
            continue;
        };
        let mut sync = BTreeSet::new();
        let mut samples = Vec::with_capacity(2);
        let include_qq = qq.contains_key(&record_index) && call_ids[1].is_some();
        let query_only = call_ids[0].is_none() && call_ids[1].is_some();
        let skipped = call_ids
            .iter()
            .flatten()
            .all(|call_id| verdicts[*call_id].skipped);
        for (sample_index, call_id) in call_ids.iter().enumerate() {
            let gt = record
                .sample_map(sample_index)
                .get("GT")
                .cloned()
                .map(|gt| if gt == "./." { ".".to_string() } else { gt })
                .unwrap_or_else(|| ".".to_string());
            let Some(call_id) = call_id else {
                samples.push(".".to_string());
                continue;
            };
            let verdict = verdicts[*call_id];
            if verdict.sync_start > 0 {
                sync.insert(verdict.sync_start);
            }
            let decision = if verdict.skipped {
                "N"
            } else if verdict.gt_match {
                "TP"
            } else if calls[*call_id].side == Side::Truth {
                "FN"
            } else {
                "FP"
            };
            let kind = if verdict.gt_match {
                "gm"
            } else if verdict.allele_match {
                "am"
            } else if verdict.loose_match {
                "lm"
            } else {
                "."
            };
            let score = if sample_index == 1 && verdict.skipped {
                "0"
            } else if sample_index == 1 {
                qq.get(&record_index).map(String::as_str).unwrap_or(".")
            } else {
                "."
            };
            let sample = if skipped {
                format!("{gt}:{decision}")
            } else if query_only {
                format!("{gt}:{score}:{decision}:{kind}")
            } else if include_qq {
                format!("{gt}:{decision}:{kind}:{score}")
            } else {
                format!("{gt}:{decision}:{kind}")
            };
            samples.push(trim_trailing_missing(&sample));
        }
        record.id = ".".to_string();
        record.qual = ".".to_string();
        if call_ids[1].is_none() {
            record.filter = ".".to_string();
        }
        record.info = if sync.is_empty() {
            ".".to_string()
        } else {
            format!(
                "BS={}",
                sync.into_iter()
                    .map(|value| value.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        };
        record.format = Some(if skipped {
            "GT:BD".to_string()
        } else if query_only {
            "GT:QQ:BD:BK".to_string()
        } else if include_qq {
            "GT:BD:BK:QQ".to_string()
        } else {
            "GT:BD:BK".to_string()
        });
        record.samples = samples;
        output.push(record);
    }
    *records = output;
}

fn trim_trailing_missing(sample: &str) -> String {
    let mut fields = sample.split(':').collect::<Vec<_>>();
    while fields.len() > 1 && fields.last() == Some(&".") {
        fields.pop();
    }
    fields.join(":")
}

fn output_headers(
    query_headers: &[String],
    references: &BTreeMap<String, String>,
    loose_match_distance: usize,
) -> Vec<String> {
    let mut headers = vec![
        "##fileformat=VCFv4.2".to_string(),
        format!("##fileDate={}", utc_date()),
        "##source=hap-rs native vcfeval (RTG Tools 3.12.1 compatible)".to_string(),
    ];
    for (name, sequence) in references {
        headers.push(format!("##contig=<ID={name},length={}>", sequence.len()));
    }
    for line in query_headers
        .iter()
        .filter(|line| line.starts_with("##FILTER=<"))
    {
        if !headers.contains(line) {
            headers.push(line.clone());
        }
    }
    headers.extend([
        "##INFO=<ID=BS,Number=.,Type=Integer,Description=\"Benchmarking superlocus ID for these variants\">".to_string(),
        "##INFO=<ID=CALL_WEIGHT,Number=1,Type=Float,Description=\"Call weight (equivalent number of truth variants). When unspecified, assume 1.0\">".to_string(),
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">".to_string(),
        format!("##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type). (Loose match distance is {loose_match_distance})\">"),
        "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">".to_string(),
        "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation\">".to_string(),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string(),
    ]);
    headers
}

fn utc_date() -> String {
    let days = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs()
        / 86_400;
    let (year, month, day) = civil_from_days(days as i64);
    format!("{year:04}{month:02}{day:02}")
}

// Howard Hinnant's public-domain civil calendar conversion.
fn civil_from_days(days_since_epoch: i64) -> (i64, i64, i64) {
    let z = days_since_epoch + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let day_of_era = z - era * 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let mut year = year_of_era + era * 400;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_prime = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * month_prime + 2) / 5 + 1;
    let month = month_prime + if month_prime < 10 { 3 } else { -9 };
    year += i64::from(month <= 2);
    (year, month, day)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn call(side: Side, record_index: usize, pos: usize, alternate: &str) -> Call {
        Call {
            record_index,
            side,
            chrom: "chr1".to_string(),
            pos,
            ref_allele: "A".to_string(),
            alts: vec![alternate.to_string()],
            gt: vec![Some(1), Some(1)],
            skipped: false,
        }
    }

    #[test]
    fn civil_date_matches_unix_epoch() {
        assert_eq!(civil_from_days(0), (1970, 1, 1));
        assert_eq!(civil_from_days(20_671), (2026, 8, 6));
    }

    #[test]
    fn duplicate_equivalent_calls_do_not_create_an_impossible_match_set() -> Result<()> {
        let calls = vec![
            call(Side::Truth, 0, 2, "C"),
            call(Side::Truth, 1, 2, "C"),
            call(Side::Query, 2, 2, "C"),
        ];
        let cluster = Cluster {
            chrom: "chr1".to_string(),
            start: 2,
            end: 2,
            calls: vec![0, 1, 2],
        };
        let mut verdicts = vec![Verdict::default(); calls.len()];

        compare_cluster(&cluster, &calls, "AAAA", &mut verdicts)?;

        assert_eq!(
            verdicts[..2]
                .iter()
                .filter(|verdict| verdict.gt_match)
                .count(),
            1
        );
        assert!(verdicts[2].gt_match);
        Ok(())
    }

    #[test]
    fn path_budget_exhaustion_is_an_explicit_error() {
        assert_eq!(
            checked_state_count(RTG_MAX_PATHS, 0).unwrap(),
            RTG_MAX_PATHS
        );
        let error = checked_state_count(RTG_MAX_PATHS, 1).unwrap_err();
        assert!(error.to_string().contains("exceeding the limit"));
    }

    #[test]
    fn query_scores_follow_source_provenance_at_duplicate_coordinates() -> Result<()> {
        let path = Path::new("scores.vcf");
        let low = RawVcfRecord::from_line("chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t1/1", path)?;
        let high = RawVcfRecord::from_line("chr1\t2\t.\tA\tG\t99\tPASS\t.\tGT\t1/1", path)?;
        let query = vec![low.clone(), high.clone()];
        let merged = vec![
            simple_two_sample_record(&high, Side::Query),
            simple_two_sample_record(&low, Side::Query),
        ];

        let scores = query_scores(&merged, &query, &[Some(1), Some(0)], "QUAL")?;

        assert_eq!(scores.get(&0).map(String::as_str), Some("99.0"));
        assert_eq!(scores.get(&1).map(String::as_str), Some("10.0"));
        Ok(())
    }

    #[test]
    fn equivalent_shifted_insertions_share_an_allele() -> Result<()> {
        let truth = Call {
            record_index: 0,
            side: Side::Truth,
            chrom: "chr1".to_string(),
            pos: 2,
            ref_allele: "A".to_string(),
            alts: vec!["AA".to_string()],
            gt: vec![Some(1), Some(1)],
            skipped: false,
        };
        let query = Call {
            record_index: 1,
            side: Side::Query,
            chrom: "chr1".to_string(),
            pos: 3,
            ref_allele: "A".to_string(),
            alts: vec!["AA".to_string()],
            gt: vec![Some(1), Some(1)],
            skipped: false,
        };
        assert!(calls_share_allele(&truth, &query, "CAAAAG")?);
        Ok(())
    }

    #[test]
    fn loose_match_uses_strict_untrimmed_span_distance() {
        let call = |side, pos| Call {
            record_index: pos,
            side,
            chrom: "chr1".to_string(),
            pos,
            ref_allele: "A".to_string(),
            alts: vec!["C".to_string()],
            gt: vec![Some(1), Some(1)],
            skipped: false,
        };
        let calls = vec![call(Side::Truth, 10), call(Side::Query, 12)];
        let mut verdicts = vec![Verdict::default(); 2];
        apply_loose_matches(&calls, 2, &mut verdicts);
        assert!(verdicts[0].loose_match);
        assert!(verdicts[1].loose_match);
    }

    #[test]
    fn symbolic_calls_are_unscored_and_have_no_match_metadata() -> Result<()> {
        let mut records = vec![RawVcfRecord::from_line(
            "chr1\t7\t.\tA\t<DEL>\t20\tPASS\tSCORE=97\tGT\t1/1\t1/1",
            Path::new("symbolic.vcf"),
        )?];
        let calls = extract_calls(&records)?;
        let verdicts = vec![
            Verdict {
                skipped: true,
                ..Verdict::default()
            };
            calls.len()
        ];
        let scores = BTreeMap::from([(0, "97.0".to_string())]);

        annotate_records(&mut records, &calls, &verdicts, &scores);

        assert_eq!(records[0].info, ".");
        assert_eq!(records[0].format.as_deref(), Some("GT:BD"));
        assert_eq!(records[0].samples, ["1/1:N", "1/1:N"]);
        Ok(())
    }
}
