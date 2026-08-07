//! Pure-Rust implementation of the legacy two-sample `scmp` comparison.

use crate::vcf::{self, RawVcfRecord};
use anyhow::{Result, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ScmpMode {
    Alleles,
    Distance { max_distance: i64 },
}

#[derive(Clone, Debug, Default)]
pub struct MergedVcf {
    pub headers: Vec<String>,
    pub records: Vec<RawVcfRecord>,
}

/// Header declarations added by legacy `scmp`.
pub const SCMP_HEADER_DECLARATIONS: [&str; 9] = [
    "##INFO=<ID=BS,Number=.,Type=Integer,Description=\"Benchmarking superlocus ID for these variants.\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">",
    "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">",
    "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">",
    "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">",
    "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">",
    "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"High-level variant type (SNP|INDEL).\">",
    "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"High-level location type (het|homref|hetalt|homalt|nocall).\">",
];

pub fn compare_files(
    truth: &Path,
    query: &Path,
    reference: &Path,
    mode: ScmpMode,
    qq_field: &str,
    output: &Path,
) -> Result<()> {
    let (truth_headers, truth_records) = vcf::load_raw_vcf(truth)?;
    let (query_headers, query_records) = vcf::load_raw_vcf(query)?;
    let mut merged = merge_two_sample_records(
        &truth_headers,
        &truth_records,
        &query_headers,
        &query_records,
    )?;
    let references = crate::fasta::read_sequences(reference)?;
    annotate_merged_records(
        &mut merged.records,
        &merged.headers,
        &references,
        mode,
        qq_field,
    )?;
    vcf::write_raw_vcf(output, &merged.headers, &merged.records)
}

pub fn merge_two_sample_records(
    truth_headers: &[String],
    truth_records: &[RawVcfRecord],
    query_headers: &[String],
    query_records: &[RawVcfRecord],
) -> Result<MergedVcf> {
    let headers = merged_headers(truth_headers, query_headers);
    let contig_order = contig_order(&headers);
    let truth_numbers = field_numbers(truth_headers);
    let query_numbers = field_numbers(query_headers);
    let mut records = Vec::with_capacity(truth_records.len() + query_records.len());
    let mut truth_index = 0;
    let mut query_index = 0;

    while truth_index < truth_records.len() || query_index < query_records.len() {
        match (
            truth_records.get(truth_index),
            query_records.get(query_index),
        ) {
            (Some(truth), Some(query)) => {
                let order = record_order(truth, query, &contig_order);
                let chrom = if order == std::cmp::Ordering::Greater {
                    query.chrom.clone()
                } else {
                    truth.chrom.clone()
                };
                let pos = if order == std::cmp::Ordering::Greater {
                    query.pos
                } else {
                    truth.pos
                };
                let truth_group =
                    take_coordinate_group(truth_records, &mut truth_index, &chrom, pos);
                let query_group =
                    take_coordinate_group(query_records, &mut query_index, &chrom, pos);
                merge_coordinate_group(
                    &mut records,
                    truth_group,
                    query_group,
                    &truth_numbers,
                    &query_numbers,
                )?;
            }
            (Some(truth), None) => {
                let chrom = truth.chrom.clone();
                let pos = truth.pos;
                let truth_group =
                    take_coordinate_group(truth_records, &mut truth_index, &chrom, pos);
                merge_coordinate_group(
                    &mut records,
                    truth_group,
                    &[],
                    &truth_numbers,
                    &query_numbers,
                )?;
            }
            (None, Some(query)) => {
                let chrom = query.chrom.clone();
                let pos = query.pos;
                let query_group =
                    take_coordinate_group(query_records, &mut query_index, &chrom, pos);
                merge_coordinate_group(
                    &mut records,
                    &[],
                    query_group,
                    &truth_numbers,
                    &query_numbers,
                )?;
            }
            (None, None) => break,
        }
    }
    Ok(MergedVcf { headers, records })
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum FieldNumber {
    A,
    R,
    G,
    Other,
}

fn merged_headers(truth: &[String], query: &[String]) -> Vec<String> {
    let mut headers = Vec::new();
    let mut seen = BTreeSet::new();
    let mut seen_declarations = BTreeSet::new();
    for (source, side_name) in [(truth, "TRUTH"), (query, "QUERY")] {
        for raw_line in source {
            if raw_line.starts_with("#CHROM") {
                continue;
            }
            if raw_line.starts_with("##fileformat=") {
                if side_name == "TRUTH" && seen.insert("##fileformat=VCFv4.2".to_string()) {
                    headers.push("##fileformat=VCFv4.2".to_string());
                }
                continue;
            }
            if raw_line.starts_with("##PEDIGREE=")
                || (side_name == "QUERY"
                    && (raw_line.starts_with("##reference=") || raw_line.starts_with("##source=")))
            {
                // The legacy merge keeps query-side generic provenance but
                // suppresses these three input-identifying metadata families.
                continue;
            }
            let line = rename_sample_info_header(raw_line, side_name);
            if let Some(identity) = header_identity(&line)
                && !seen_declarations.insert(identity)
            {
                continue;
            }
            if seen.insert(line.clone()) {
                headers.push(line);
            }
        }
    }
    for declaration in SCMP_HEADER_DECLARATIONS {
        let id = declaration_id(declaration).unwrap_or_default();
        let kind = if declaration.starts_with("##INFO") {
            "##INFO=<ID="
        } else {
            "##FORMAT=<ID="
        };
        if !headers
            .iter()
            .any(|line| line.starts_with(&format!("{kind}{id},")))
        {
            headers.push(declaration.to_string());
        }
    }
    headers.push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string());
    headers
}

fn rename_sample_info_header(line: &str, side_name: &str) -> String {
    if line.starts_with("##INFO=<ID=SAMPLE_") {
        line.replacen("ID=SAMPLE_", &format!("ID={side_name}_"), 1)
    } else {
        line.to_string()
    }
}

fn header_identity(line: &str) -> Option<String> {
    let kind = if line.starts_with("##INFO=<") {
        "INFO"
    } else if line.starts_with("##FORMAT=<") {
        "FORMAT"
    } else if line.starts_with("##FILTER=<") {
        "FILTER"
    } else if line.starts_with("##contig=<") {
        "contig"
    } else {
        return None;
    };
    Some(format!("{kind}:{}", declaration_id(line)?))
}

fn declaration_id(line: &str) -> Option<&str> {
    let rest = line.split("ID=").nth(1)?;
    rest.split([',', '>']).next()
}

fn field_numbers(headers: &[String]) -> BTreeMap<String, FieldNumber> {
    let mut result = BTreeMap::new();
    for line in headers {
        if !line.starts_with("##FORMAT=<") && !line.starts_with("##INFO=<") {
            continue;
        }
        let Some(id) = declaration_id(line) else {
            continue;
        };
        let number = line
            .split("Number=")
            .nth(1)
            .and_then(|rest| rest.split([',', '>']).next())
            .unwrap_or(".");
        let number = match number {
            "A" => FieldNumber::A,
            "R" => FieldNumber::R,
            "G" => FieldNumber::G,
            _ => FieldNumber::Other,
        };
        result.insert(id.to_string(), number);
    }
    result
}

fn contig_order(headers: &[String]) -> BTreeMap<String, usize> {
    let mut order = BTreeMap::new();
    for line in headers {
        let Some(rest) = line.strip_prefix("##contig=<ID=") else {
            continue;
        };
        let Some(name) = rest.split([',', '>']).next() else {
            continue;
        };
        let next = order.len();
        order.entry(name.to_string()).or_insert(next);
    }
    order
}

fn record_order(
    left: &RawVcfRecord,
    right: &RawVcfRecord,
    contigs: &BTreeMap<String, usize>,
) -> std::cmp::Ordering {
    match (contigs.get(&left.chrom), contigs.get(&right.chrom)) {
        (Some(left_rank), Some(right_rank)) => left_rank.cmp(right_rank),
        _ => left.chrom.cmp(&right.chrom),
    }
    .then_with(|| left.pos.cmp(&right.pos))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum MergeClass {
    Substitution,
    Indel,
    Other,
}

fn merge_class(record: &RawVcfRecord) -> MergeClass {
    let ref_len = record.ref_allele.len();
    let mut saw_indel = false;
    for alt in record.alt_allele.split(',') {
        if alt == "." || alt == "*" || alt.starts_with('<') {
            return MergeClass::Other;
        }
        saw_indel |= alt.len() != ref_len;
    }
    if saw_indel {
        MergeClass::Indel
    } else {
        MergeClass::Substitution
    }
}

fn take_coordinate_group<'a>(
    records: &'a [RawVcfRecord],
    index: &mut usize,
    chrom: &str,
    pos: usize,
) -> &'a [RawVcfRecord] {
    let start = *index;
    while records
        .get(*index)
        .is_some_and(|record| record.chrom == chrom && record.pos == pos)
    {
        *index += 1;
    }
    &records[start..*index]
}

fn merge_coordinate_group(
    output: &mut Vec<RawVcfRecord>,
    truth: &[RawVcfRecord],
    query: &[RawVcfRecord],
    truth_numbers: &BTreeMap<String, FieldNumber>,
    query_numbers: &BTreeMap<String, FieldNumber>,
) -> Result<()> {
    for class in [
        MergeClass::Substitution,
        MergeClass::Indel,
        MergeClass::Other,
    ] {
        let mut truth_ids: Vec<usize> = truth
            .iter()
            .enumerate()
            .filter_map(|(index, record)| (merge_class(record) == class).then_some(index))
            .collect();
        let mut query_ids: Vec<usize> = query
            .iter()
            .enumerate()
            .filter_map(|(index, record)| (merge_class(record) == class).then_some(index))
            .collect();
        let mut allele_counts = Vec::<usize>::new();

        while !truth_ids.is_empty() || !query_ids.is_empty() {
            // bcftools 1.17 rebuilds an allele table for every emitted line,
            // selects the most frequent ALT (first wins ties), and chooses the
            // first record from each input carrying that ALT. Its count array
            // is reused without clearing newly reintroduced slots; this
            // historical quirk affects duplicate indel positions.
            crate::compatibility::begin_allele_count_table(
                crate::compatibility::AlleleCountArrayPolicy::LegacyReuse,
                &mut allele_counts,
            );
            let mut allele_order = Vec::<AlleleKey>::new();
            let mut first_record = true;
            for record in truth_ids
                .iter()
                .map(|id| &truth[*id])
                .chain(query_ids.iter().map(|id| &query[*id]))
            {
                for key in allele_keys(record) {
                    let slot = match allele_order.iter().position(|existing| existing == &key) {
                        Some(slot) => slot,
                        None => {
                            allele_order.push(key);
                            allele_order.len() - 1
                        }
                    };
                    crate::compatibility::record_allele_count(
                        &mut allele_counts,
                        slot,
                        first_record,
                    );
                }
                first_record = false;
            }
            let selected_slot = (1..allele_order.len()).fold(0, |best, slot| {
                if allele_counts[slot] > allele_counts[best] {
                    slot
                } else {
                    best
                }
            });
            let selected = &allele_order[selected_slot];
            let truth_slot = truth_ids
                .iter()
                .position(|id| record_has_allele(&truth[*id], selected))
                .or((!truth_ids.is_empty()).then_some(0));
            let query_slot = query_ids
                .iter()
                .position(|id| record_has_allele(&query[*id], selected))
                .or((!query_ids.is_empty()).then_some(0));

            match (truth_slot, query_slot) {
                (Some(truth_slot), Some(query_slot)) => {
                    let truth_id = truth_ids.remove(truth_slot);
                    let query_id = query_ids.remove(query_slot);
                    output.push(merge_record_pair_ordered(
                        &truth[truth_id],
                        &query[query_id],
                        truth_numbers,
                        query_numbers,
                        &allele_order,
                    )?);
                }
                (Some(truth_slot), None) => {
                    output.push(single_side_record(
                        &truth[truth_ids.remove(truth_slot)],
                        0,
                        truth_numbers,
                        query_numbers,
                    ));
                }
                (None, Some(query_slot)) => {
                    output.push(single_side_record(
                        &query[query_ids.remove(query_slot)],
                        1,
                        truth_numbers,
                        query_numbers,
                    ));
                }
                (None, None) => unreachable!("non-empty coordinate group"),
            }
        }
    }
    Ok(())
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct AlleleKey {
    pos: usize,
    reference: String,
    alt: String,
}

fn allele_keys(record: &RawVcfRecord) -> Vec<AlleleKey> {
    record
        .alt_allele
        .split(',')
        .map(|alt| {
            let mut allele = record.clone();
            allele.alt_allele = alt.to_string();
            minimize_record_alleles(&mut allele);
            AlleleKey {
                pos: allele.pos,
                reference: allele.ref_allele,
                alt: allele.alt_allele,
            }
        })
        .collect()
}

fn record_has_allele(record: &RawVcfRecord, selected: &AlleleKey) -> bool {
    allele_keys(record).iter().any(|key| key == selected)
}

fn single_side_record(
    source: &RawVcfRecord,
    side: usize,
    truth_numbers: &BTreeMap<String, FieldNumber>,
    query_numbers: &BTreeMap<String, FieldNumber>,
) -> RawVcfRecord {
    let mut record = source.clone();
    let mut keys: Vec<String> = record
        .format_keys()
        .into_iter()
        .map(str::to_string)
        .collect();
    if keys.is_empty() {
        keys.push("GT".to_string());
        record.format = Some("GT".to_string());
    }
    let called = source
        .samples
        .first()
        .cloned()
        .unwrap_or_else(|| keys.iter().map(|_| ".").collect::<Vec<_>>().join(":"));
    let missing = keys
        .iter()
        .map(|key| if key == "GT" { "./." } else { "." })
        .collect::<Vec<_>>()
        .join(":");
    record.samples = if side == 0 {
        vec![called, missing]
    } else {
        vec![missing, called]
    };
    record.info =
        rename_sample_info_values(&record.info, if side == 0 { "TRUTH" } else { "QUERY" });
    canonicalize_bcftools_info(&mut record, truth_numbers, query_numbers);
    minimize_record_alleles(&mut record);
    record
}

fn merge_record_pair_ordered(
    truth: &RawVcfRecord,
    query: &RawVcfRecord,
    truth_numbers: &BTreeMap<String, FieldNumber>,
    query_numbers: &BTreeMap<String, FieldNumber>,
    allele_order: &[AlleleKey],
) -> Result<RawVcfRecord> {
    let merged_ref = if truth.ref_allele.len() >= query.ref_allele.len() {
        truth.ref_allele.clone()
    } else {
        query.ref_allele.clone()
    };
    let mut preferred_alts = Vec::new();
    for selected in allele_order {
        append_selected_alt(truth, &merged_ref, selected, &mut preferred_alts)?;
        append_selected_alt(query, &merged_ref, selected, &mut preferred_alts)?;
    }
    let (alts, truth_map) = merged_alts(truth, &merged_ref, preferred_alts)?;
    let (alts, query_map) = merged_alts(query, &merged_ref, alts)?;

    let mut format_keys: Vec<String> = truth
        .format_keys()
        .into_iter()
        .chain(query.format_keys())
        .map(str::to_string)
        .collect();
    if !format_keys.iter().any(|key| key == "GT") {
        format_keys.insert(0, "GT".to_string());
    }
    let mut seen = BTreeSet::new();
    format_keys.retain(|key| seen.insert(key.clone()));
    let truth_sample = remap_sample(truth, &format_keys, &truth_map, alts.len(), truth_numbers)?;
    let query_sample = remap_sample(query, &format_keys, &query_map, alts.len(), query_numbers)?;

    let mut merged = RawVcfRecord {
        chrom: truth.chrom.clone(),
        pos: truth.pos,
        id: merge_scalar(&truth.id, &query.id, ";"),
        ref_allele: merged_ref,
        alt_allele: alts.join(","),
        qual: max_qual(&truth.qual, &query.qual),
        filter: merge_filter(&truth.filter, &query.filter),
        info: merge_info(
            truth,
            query,
            &truth_map,
            &query_map,
            alts.len(),
            truth_numbers,
            query_numbers,
        ),
        format: Some(format_keys.join(":")),
        samples: vec![truth_sample, query_sample],
    };
    canonicalize_bcftools_info(&mut merged, truth_numbers, query_numbers);
    minimize_record_alleles(&mut merged);
    Ok(merged)
}

/// Reproduce the INFO layout that bcftools merge writes after combining
/// samples. Allele-count tags are regenerated from the merged GTs, while
/// merge-sensitive fields are emitted after ordinary annotations.
fn canonicalize_bcftools_info(
    record: &mut RawVcfRecord,
    truth_numbers: &BTreeMap<String, FieldNumber>,
    query_numbers: &BTreeMap<String, FieldNumber>,
) {
    let emit_allele_count = truth_numbers.contains_key("AC") || query_numbers.contains_key("AC");
    let emit_allele_number = truth_numbers.contains_key("AN") || query_numbers.contains_key("AN");
    let mut ordinary = Vec::new();
    let mut depth = Vec::new();
    let mut frequency = Vec::new();
    let mut genotype_numbered = Vec::new();
    for entry in record
        .info
        .split(';')
        .filter(|entry| !entry.is_empty() && *entry != ".")
    {
        let key = entry.split_once('=').map_or(entry, |(key, _)| key);
        if matches!(key, "AC" | "AN") {
            continue;
        }
        let number = merged_field_number(key, truth_numbers, query_numbers);
        let entry = match number {
            FieldNumber::A => truncate_info_values(entry, record.alt_allele.split(',').count()),
            FieldNumber::G => {
                let sample_index = usize::from(key.starts_with("QUERY_"));
                let ploidy = record
                    .sample_map(sample_index)
                    .get("GT")
                    .map_or(2, |gt| gt.split(['/', '|']).count());
                let allele_count = record.alt_allele.split(',').count() + 1;
                let cardinality = match ploidy {
                    1 => allele_count,
                    2 => allele_count * (allele_count + 1) / 2,
                    _ => usize::MAX,
                };
                truncate_info_values(entry, cardinality)
            }
            _ => entry.to_string(),
        };
        let target = if key == "DP" {
            &mut depth
        } else if key == "AF" {
            &mut frequency
        } else if number == FieldNumber::G {
            &mut genotype_numbered
        } else {
            &mut ordinary
        };
        target.push(entry);
    }

    // bcftools emits query SAMPLE_* annotations in merged-header order.
    // The governed preprocess header uses the same lexical order (AD, DP,
    // GQ, GQX, GT, MQ, PL, VF); Number=G fields such as PL stay deferred.
    let mut query_sample_info = Vec::new();
    ordinary.retain(|entry| {
        let key = entry.split_once('=').map_or(entry.as_str(), |(key, _)| key);
        if key.starts_with("QUERY_") {
            query_sample_info.push(entry.clone());
            false
        } else {
            true
        }
    });
    query_sample_info.sort_by(|left, right| {
        let left = left.split_once('=').map_or(left.as_str(), |(key, _)| key);
        let right = right.split_once('=').map_or(right.as_str(), |(key, _)| key);
        left.cmp(right)
    });
    ordinary.extend(query_sample_info);

    ordinary.extend(depth);
    ordinary.extend(frequency);
    ordinary.extend(genotype_numbered);

    let mut counts = vec![0usize; record.alt_allele.split(',').count()];
    let mut allele_number = 0usize;
    for sample_index in 0..record.samples.len() {
        if let Some(gt) = record.sample_map(sample_index).get("GT") {
            for allele in gt.split(['/', '|']) {
                let Ok(index) = allele.parse::<usize>() else {
                    continue;
                };
                allele_number += 1;
                if index > 0
                    && let Some(count) = counts.get_mut(index - 1)
                {
                    *count += 1;
                }
            }
        }
    }
    if emit_allele_number {
        ordinary.push(format!("AN={allele_number}"));
    }
    if emit_allele_count {
        ordinary.push(format!(
            "AC={}",
            counts
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        ));
    }
    record.info = ordinary.join(";");
}

fn truncate_info_values(entry: &str, cardinality: usize) -> String {
    let Some((key, value)) = entry.split_once('=') else {
        return entry.to_string();
    };
    let values: Vec<&str> = value.split(',').collect();
    if values.len() <= cardinality {
        return entry.to_string();
    }
    format!("{key}={}", values[..cardinality].join(","))
}

fn merged_field_number(
    key: &str,
    truth_numbers: &BTreeMap<String, FieldNumber>,
    query_numbers: &BTreeMap<String, FieldNumber>,
) -> FieldNumber {
    if let Some(suffix) = key.strip_prefix("TRUTH_") {
        return truth_numbers
            .get(&format!("SAMPLE_{suffix}"))
            .copied()
            .unwrap_or(FieldNumber::Other);
    }
    if let Some(suffix) = key.strip_prefix("QUERY_") {
        return query_numbers
            .get(&format!("SAMPLE_{suffix}"))
            .copied()
            .unwrap_or(FieldNumber::Other);
    }
    truth_numbers
        .get(key)
        .or_else(|| query_numbers.get(key))
        .copied()
        .unwrap_or(FieldNumber::Other)
}

fn append_selected_alt(
    record: &RawVcfRecord,
    merged_ref: &str,
    selected: &AlleleKey,
    output: &mut Vec<String>,
) -> Result<()> {
    let Some(suffix) = merged_ref.strip_prefix(&record.ref_allele) else {
        bail!(
            "incompatible REF alleles at {}:{}",
            record.chrom,
            record.pos
        );
    };
    for (alt, key) in record.alt_allele.split(',').zip(allele_keys(record)) {
        if &key != selected {
            continue;
        }
        let expanded = if alt.starts_with('<') || alt == "*" || alt == "." {
            alt.to_string()
        } else {
            format!("{alt}{suffix}")
        };
        if !output.contains(&expanded) {
            output.push(expanded);
        }
    }
    Ok(())
}

fn minimize_record_alleles(record: &mut RawVcfRecord) {
    let mut alts: Vec<String> = record.alt_allele.split(',').map(str::to_string).collect();
    if alts
        .iter()
        .any(|alt| alt == "." || alt == "*" || alt.starts_with('<'))
    {
        return;
    }

    while record.ref_allele.len() > 1
        && alts.iter().all(|alt| alt.len() > 1)
        && alts
            .iter()
            .all(|alt| alt.as_bytes().last() == record.ref_allele.as_bytes().last())
    {
        record.ref_allele.pop();
        for alt in &mut alts {
            alt.pop();
        }
    }
    while record.ref_allele.len() > 1
        && alts.iter().all(|alt| alt.len() > 1)
        && alts
            .iter()
            .all(|alt| alt.as_bytes().first() == record.ref_allele.as_bytes().first())
    {
        record.ref_allele.remove(0);
        for alt in &mut alts {
            alt.remove(0);
        }
        record.pos += 1;
    }
    record.alt_allele = alts.join(",");
}

fn merged_alts(
    record: &RawVcfRecord,
    merged_ref: &str,
    mut merged: Vec<String>,
) -> Result<(Vec<String>, Vec<usize>)> {
    let Some(suffix) = merged_ref.strip_prefix(&record.ref_allele) else {
        bail!(
            "incompatible REF alleles at {}:{}",
            record.chrom,
            record.pos
        );
    };
    let mut mapping = vec![0];
    for alt in record.alt_allele.split(',') {
        let expanded = if alt.starts_with('<') || alt == "*" || alt == "." {
            alt.to_string()
        } else {
            format!("{alt}{suffix}")
        };
        let index = match merged.iter().position(|value| value == &expanded) {
            Some(index) => index + 1,
            None => {
                merged.push(expanded);
                merged.len()
            }
        };
        mapping.push(index);
    }
    Ok((merged, mapping))
}

fn remap_sample(
    record: &RawVcfRecord,
    merged_keys: &[String],
    allele_map: &[usize],
    new_alt_count: usize,
    numbers: &BTreeMap<String, FieldNumber>,
) -> Result<String> {
    let source = record.sample_map(0);
    let gt = source.get("GT").map(String::as_str).unwrap_or("./.");
    merged_keys
        .iter()
        .map(|key| {
            let Some(value) = source.get(key) else {
                return Ok(if key == "GT" {
                    "./.".to_string()
                } else {
                    ".".to_string()
                });
            };
            if key == "GT" {
                remap_gt(value, allele_map)
            } else {
                Ok(remap_numbered(
                    value,
                    numbers.get(key).copied().unwrap_or(FieldNumber::Other),
                    allele_map,
                    new_alt_count,
                    gt,
                ))
            }
        })
        .collect::<Result<Vec<_>>>()
        .map(|values| values.join(":"))
}

fn remap_gt(gt: &str, allele_map: &[usize]) -> Result<String> {
    let mut output = String::new();
    let mut token = String::new();
    for ch in gt.chars().chain(std::iter::once('/')) {
        if ch == '/' || ch == '|' {
            if token == "." || token.is_empty() {
                output.push('.');
            } else {
                let old = token.parse::<usize>()?;
                let Some(new) = allele_map.get(old) else {
                    bail!("GT allele index {old} exceeds ALT count");
                };
                output.push_str(&new.to_string());
            }
            token.clear();
            output.push(ch);
        } else {
            token.push(ch);
        }
    }
    output.pop();
    Ok(output)
}

fn remap_numbered(
    value: &str,
    number: FieldNumber,
    allele_map: &[usize],
    new_alt_count: usize,
    gt: &str,
) -> String {
    if value == "." {
        return value.to_string();
    }
    let old: Vec<&str> = value.split(',').collect();
    match number {
        FieldNumber::A => {
            let mut new = vec!["."; new_alt_count];
            for (old_index, value) in old.iter().enumerate() {
                if let Some(&new_index) = allele_map.get(old_index + 1) {
                    new[new_index - 1] = value;
                }
            }
            new.join(",")
        }
        FieldNumber::R => {
            let mut new = vec!["."; new_alt_count + 1];
            for (old_index, value) in old.iter().enumerate() {
                if let Some(&new_index) = allele_map.get(old_index) {
                    new[new_index] = value;
                }
            }
            new.join(",")
        }
        FieldNumber::G => remap_genotypes(&old, allele_map, new_alt_count + 1, gt),
        FieldNumber::Other => value.to_string(),
    }
}

fn remap_genotypes(old: &[&str], allele_map: &[usize], new_alleles: usize, gt: &str) -> String {
    let ploidy = gt.split(['/', '|']).count();
    if ploidy == 1 && old.len() == allele_map.len() {
        let mut new = vec!["."; new_alleles];
        for (old_index, value) in old.iter().enumerate() {
            new[allele_map[old_index]] = value;
        }
        return new.join(",");
    }
    let old_cardinality = allele_map.len() * (allele_map.len() + 1) / 2;
    if ploidy != 2 || old.len() < old_cardinality {
        return old.join(",");
    }
    let old = &old[..old_cardinality];
    let mut new = vec!["."; new_alleles * (new_alleles + 1) / 2];
    for high in 0..allele_map.len() {
        for low in 0..=high {
            let old_index = high * (high + 1) / 2 + low;
            let (mapped_low, mapped_high) = if allele_map[low] <= allele_map[high] {
                (allele_map[low], allele_map[high])
            } else {
                (allele_map[high], allele_map[low])
            };
            let new_index = mapped_high * (mapped_high + 1) / 2 + mapped_low;
            new[new_index] = old[old_index];
        }
    }
    new.join(",")
}

fn merge_scalar(left: &str, right: &str, separator: &str) -> String {
    match (left, right) {
        (".", value) | ("", value) => value.to_string(),
        (value, ".") | (value, "") => value.to_string(),
        (left, right) if left == right => left.to_string(),
        _ => format!("{left}{separator}{right}"),
    }
}

fn max_qual(left: &str, right: &str) -> String {
    match (left.parse::<f32>(), right.parse::<f32>()) {
        (Ok(left_value), Ok(right_value)) if right_value > left_value => right.to_string(),
        (Ok(_), _) => left.to_string(),
        (_, Ok(_)) => right.to_string(),
        _ => ".".to_string(),
    }
}

fn merge_filter(left: &str, right: &str) -> String {
    let mut filters = Vec::new();
    for filter in [left, right] {
        for value in filter.split(';') {
            if value != "." && value != "PASS" && !filters.contains(&value) {
                filters.push(value);
            }
        }
    }
    if filters.is_empty() {
        if left == "." && right == "." {
            ".".to_string()
        } else {
            "PASS".to_string()
        }
    } else {
        filters.join(";")
    }
}

fn merge_info(
    truth: &RawVcfRecord,
    query: &RawVcfRecord,
    truth_map: &[usize],
    query_map: &[usize],
    new_alt_count: usize,
    truth_numbers: &BTreeMap<String, FieldNumber>,
    query_numbers: &BTreeMap<String, FieldNumber>,
) -> String {
    let mut values: BTreeMap<String, Option<String>> = BTreeMap::new();
    let mut order = Vec::new();
    for (record, allele_map, numbers, side_name) in [
        (truth, truth_map, truth_numbers, "TRUTH"),
        (query, query_map, query_numbers, "QUERY"),
    ] {
        for entry in record.info.split(';') {
            if entry == "." || entry.is_empty() {
                continue;
            }
            let (raw_key, raw_value) = match entry.split_once('=') {
                Some((key, value)) => (key, Some(value)),
                None => (entry, None),
            };
            let key = renamed_sample_info_key(raw_key, side_name);
            if !values.contains_key(key.as_ref()) {
                order.push(key.to_string());
            }
            let remapped = raw_value.map(|value| {
                remap_numbered(
                    value,
                    numbers.get(raw_key).copied().unwrap_or(FieldNumber::Other),
                    allele_map,
                    new_alt_count,
                    "0/0",
                )
            });
            match (values.get_mut(key.as_ref()), remapped) {
                (None, value) => {
                    values.insert(key.to_string(), value);
                }
                (Some(Some(existing)), Some(incoming)) => {
                    let mut existing_fields: Vec<&str> = existing.split(',').collect();
                    let incoming_fields: Vec<&str> = incoming.split(',').collect();
                    if existing_fields.len() == incoming_fields.len() {
                        for (slot, incoming) in existing_fields.iter_mut().zip(incoming_fields) {
                            if *slot == "." && incoming != "." {
                                *slot = incoming;
                            }
                        }
                        *existing = existing_fields.join(",");
                    }
                }
                _ => {}
            }
        }
    }
    if order.is_empty() {
        return ".".to_string();
    }
    order
        .into_iter()
        .map(|key| match values.remove(&key).flatten() {
            Some(value) => format!("{key}={value}"),
            None => key,
        })
        .collect::<Vec<_>>()
        .join(";")
}

fn renamed_sample_info_key<'a>(key: &'a str, side_name: &str) -> std::borrow::Cow<'a, str> {
    match key.strip_prefix("SAMPLE_") {
        Some(suffix) => std::borrow::Cow::Owned(format!("{side_name}_{suffix}")),
        None => std::borrow::Cow::Borrowed(key),
    }
}

fn rename_sample_info_values(info: &str, side_name: &str) -> String {
    info.split(';')
        .map(|entry| match entry.split_once('=') {
            Some((key, value)) => format!("{}={value}", renamed_sample_info_key(key, side_name)),
            None => renamed_sample_info_key(entry, side_name).into_owned(),
        })
        .collect::<Vec<_>>()
        .join(";")
}

pub fn annotate_merged_records(
    records: &mut [RawVcfRecord],
    headers: &[String],
    references: &BTreeMap<String, String>,
    mode: ScmpMode,
    qq_field: &str,
) -> Result<()> {
    let qq_is_info = headers
        .iter()
        .any(|line| line.starts_with(&format!("##INFO=<ID={qq_field},")));
    let mut truth = Vec::new();
    let mut query = Vec::new();
    let mut current_chrom: Option<String> = None;
    let mut current_block_end: Option<i64> = None;

    for record_index in 0..records.len() {
        let chrom_changed = current_chrom
            .as_deref()
            .is_some_and(|chrom| chrom != records[record_index].chrom);
        if chrom_changed {
            compare_group(records, references, mode, &truth, &query)?;
            truth.clear();
            query.clear();
            current_block_end = None;
        }
        current_chrom = Some(records[record_index].chrom.clone());

        truth.extend(extract_occurrences(
            &records[record_index],
            record_index,
            0,
        )?);
        query.extend(extract_occurrences(
            &records[record_index],
            record_index,
            1,
        )?);
        let qq = qq_values(&records[record_index], qq_field, qq_is_info)?;
        upsert_format(&mut records[record_index], "QQ", &qq);

        let pos0 = i64::try_from(records[record_index].pos.saturating_sub(1))?;
        let gap_flush = current_block_end.is_some_and(|end| pos0 > end.saturating_add(100_000));
        if gap_flush || truth.len() + query.len() > 1024 {
            // Deliberately after occurrence extraction: this reproduces the
            // surprising legacy behavior where the boundary-pushing record is
            // part of the group being flushed.
            compare_group(records, references, mode, &truth, &query)?;
            truth.clear();
            query.clear();
            current_block_end = None;
        }
        let effective_end =
            records[record_index].effective_end_pos(Path::new("<merged SCMP record>"))?;
        let effective_end0 = i64::try_from(effective_end.saturating_sub(1))?;
        current_block_end = Some(
            current_block_end
                .unwrap_or(effective_end0)
                .max(pos0)
                .max(effective_end0),
        );
    }
    if !truth.is_empty() || !query.is_empty() {
        compare_group(records, references, mode, &truth, &query)?;
    }
    Ok(())
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RefVar {
    pub start: i64,
    pub end: i64,
    pub alt: String,
}

impl RefVar {
    pub fn legacy_key(&self) -> String {
        format!("{}-{}:{}", self.start, self.end, self.alt)
    }
}

#[derive(Clone, Debug)]
struct Occurrence {
    record_index: usize,
    var: RefVar,
}

fn extract_occurrences(
    record: &RawVcfRecord,
    record_index: usize,
    sample_index: usize,
) -> Result<Vec<Occurrence>> {
    let sample = record.sample_map(sample_index);
    let Some(gt) = sample.get("GT") else {
        return Ok(Vec::new());
    };
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    let start = i64::try_from(record.pos)? - 1;
    let mut occurrences = Vec::new();
    for allele in gt.split(['/', '|']) {
        let Ok(index) = allele.parse::<usize>() else {
            continue;
        };
        if index == 0 {
            continue;
        }
        let Some(raw_alt) = alts.get(index - 1) else {
            bail!(
                "GT allele index {index} exceeds ALT count at {}:{}",
                record.chrom,
                record.pos
            );
        };
        let alt = if *raw_alt == "<DEL>" {
            String::new()
        } else if raw_alt.starts_with('<') {
            bail!(
                "unsupported symbolic ALT at {}:{}: {}",
                record.chrom,
                record.pos,
                raw_alt
            );
        } else {
            (*raw_alt).to_string()
        };
        // Pinned SCMP constructs RefVar.end from ALT length. Keep that governed
        // emulation at the VCF-to-RefVar adapter instead of teaching the
        // normalization and matching algorithms the malformed span rule.
        let end = crate::compatibility::scmp_refvar_end(
            crate::compatibility::ScmpRefVarSpanPolicy::LegacyAltLength,
            start,
            record.ref_allele.len(),
            alt.len(),
        )?;
        occurrences.push(Occurrence {
            record_index,
            var: RefVar { start, end, alt },
        });
    }
    Ok(occurrences)
}

fn compare_group(
    records: &mut [RawVcfRecord],
    references: &BTreeMap<String, String>,
    mode: ScmpMode,
    truth: &[Occurrence],
    query: &[Occurrence],
) -> Result<()> {
    if truth.is_empty() && query.is_empty() {
        return Ok(());
    }
    let record_index = truth
        .first()
        .or_else(|| query.first())
        .map(|occurrence| occurrence.record_index)
        .expect("non-empty group");
    let chrom = &records[record_index].chrom;
    let reference = references
        .get(chrom)
        .ok_or_else(|| anyhow::anyhow!("reference has no contig {chrom}"))?;
    // Legacy FastaFile normalizes sequence case before RefVar trimming. The
    // chr21 fixture contains lowercase reference runs, so bytewise matching
    // without this conversion misses representation-equivalent alleles.
    let uppercase_reference;
    let reference = if reference.bytes().any(|base| base.is_ascii_lowercase()) {
        uppercase_reference = reference.to_ascii_uppercase();
        uppercase_reference.as_str()
    } else {
        reference.as_str()
    };
    let (truth_matched, query_matched) = match mode {
        ScmpMode::Alleles => allele_matches(reference, truth, query),
        ScmpMode::Distance { max_distance } => {
            if max_distance < 0 {
                bail!("SCMP maximum distance must be non-negative");
            }
            distance_matches(max_distance, truth, query)
        }
    }?;

    let mut updates: BTreeMap<usize, [String; 2]> = BTreeMap::new();
    for (occurrence, matched) in truth.iter().zip(truth_matched) {
        updates
            .entry(occurrence.record_index)
            .or_insert_with(|| [".".to_string(), ".".to_string()])[0] =
            if matched { "TP" } else { "FN" }.to_string();
    }
    for (occurrence, matched) in query.iter().zip(query_matched) {
        updates
            .entry(occurrence.record_index)
            .or_insert_with(|| [".".to_string(), ".".to_string()])[1] =
            if matched { "TP" } else { "FP" }.to_string();
    }
    for (record_index, values) in updates {
        upsert_format(&mut records[record_index], "BD", &values);
    }
    Ok(())
}

fn allele_matches(
    reference: &str,
    truth: &[Occurrence],
    query: &[Occurrence],
) -> Result<(Vec<bool>, Vec<bool>)> {
    use std::collections::VecDeque;
    let mut truth_by_key: BTreeMap<String, VecDeque<usize>> = BTreeMap::new();
    for (index, occurrence) in truth.iter().enumerate() {
        let mut var = occurrence.var.clone();
        normalize_legacy(reference, &mut var)?;
        truth_by_key
            .entry(var.legacy_key())
            .or_default()
            .push_back(index);
    }
    let mut truth_matched = vec![false; truth.len()];
    let mut query_matched = vec![false; query.len()];
    for (query_index, occurrence) in query.iter().enumerate() {
        let mut var = occurrence.var.clone();
        normalize_legacy(reference, &mut var)?;
        if let Some(truth_index) = truth_by_key
            .get_mut(&var.legacy_key())
            .and_then(VecDeque::pop_front)
        {
            truth_matched[truth_index] = true;
            query_matched[query_index] = true;
        }
    }
    Ok((truth_matched, query_matched))
}

fn distance_matches(
    max_distance: i64,
    truth: &[Occurrence],
    query: &[Occurrence],
) -> Result<(Vec<bool>, Vec<bool>)> {
    let overlaps = |left: &RefVar, right: &RefVar| {
        let left_start = left.start.min(left.end).saturating_sub(max_distance);
        let left_end = left.start.max(left.end).saturating_add(max_distance);
        let right_start = right.start.min(right.end);
        let right_end = right.start.max(right.end);
        left_start <= right_end && right_start <= left_end
    };
    Ok((
        truth
            .iter()
            .map(|left| query.iter().any(|right| overlaps(&left.var, &right.var)))
            .collect(),
        query
            .iter()
            .map(|right| truth.iter().any(|left| overlaps(&right.var, &left.var)))
            .collect(),
    ))
}

fn qq_values(record: &RawVcfRecord, field: &str, is_info: bool) -> Result<[String; 2]> {
    if is_info {
        let value = info_value(&record.info, field)
            .and_then(|value| value.split(',').next())
            .map(render_f32)
            .transpose()?
            .unwrap_or_else(|| ".".to_string());
        return Ok([value.clone(), value]);
    }
    if field == "QUAL" {
        let value = render_f32(&record.qual)?;
        return Ok([value.clone(), value]);
    }
    let values: [String; 2] = std::array::from_fn(|sample_index| {
        record
            .sample_map(sample_index)
            .get(field)
            .cloned()
            .unwrap_or_else(|| ".".to_string())
    });
    if values.iter().any(|value| value.split(',').count() > 1) {
        bail!("FORMAT/{field} must be scalar for SCMP QQ annotation");
    }
    Ok([render_f32(&values[0])?, render_f32(&values[1])?])
}

fn info_value<'a>(info: &'a str, field: &str) -> Option<&'a str> {
    info.split(';').find_map(|entry| {
        let (key, value) = entry.split_once('=')?;
        (key == field).then_some(value)
    })
}

fn render_f32(value: &str) -> Result<String> {
    if value == "." || value.is_empty() {
        return Ok(".".to_string());
    }
    let parsed = value.parse::<f32>()?;
    if parsed == 0.0 {
        Ok("0".to_string())
    } else if parsed.is_nan() {
        Ok(".".to_string())
    } else {
        Ok(parsed.to_string())
    }
}

fn upsert_format<const N: usize>(record: &mut RawVcfRecord, key: &str, values: &[String; N]) {
    debug_assert_eq!(N, record.samples.len());
    let mut keys: Vec<String> = record
        .format_keys()
        .into_iter()
        .map(str::to_string)
        .collect();
    let field_index = match keys.iter().position(|existing| existing == key) {
        Some(index) => index,
        None => {
            keys.push(key.to_string());
            keys.len() - 1
        }
    };
    for (sample_index, sample) in record.samples.iter_mut().enumerate() {
        let mut fields: Vec<String> = sample.split(':').map(str::to_string).collect();
        fields.resize(keys.len(), ".".to_string());
        fields[field_index] = values[sample_index].clone();
        *sample = fields.join(":");
    }
    record.format = Some(keys.join(":"));
}

fn normalize_legacy(reference: &str, var: &mut RefVar) -> Result<()> {
    left_shift_legacy(reference, var)?;
    trim_left_legacy(reference, var, false)?;
    trim_right_legacy(reference, var, false)
}

fn reference_slice(reference: &str, start: i64, end: i64) -> Result<&str> {
    if end < start {
        return Ok("");
    }
    if start < 0 || end < 0 {
        bail!("variant normalization moved before the reference start");
    }
    let start = usize::try_from(start)?;
    let end = usize::try_from(end)?.saturating_add(1);
    reference
        .get(start..end.min(reference.len()))
        .ok_or_else(|| anyhow::anyhow!("reference coordinates are not UTF-8 boundaries"))
}

fn trim_left_legacy(reference: &str, var: &mut RefVar, refpadding: bool) -> Result<()> {
    let reference_allele = reference_slice(reference, var.start, var.end)?.as_bytes();
    let alt = var.alt.as_bytes();
    let minimum = usize::from(refpadding);
    let mut trim = 0;
    while reference_allele.len().saturating_sub(trim) > minimum
        && alt.len().saturating_sub(trim) > minimum
        && reference_allele[trim] == alt[trim]
    {
        trim += 1;
        var.start += 1;
    }
    if trim > 0 {
        var.alt.drain(..trim);
    }
    Ok(())
}

fn trim_right_legacy(reference: &str, var: &mut RefVar, refpadding: bool) -> Result<()> {
    let mut ref_len = var.end - var.start + 1;
    let mut alt_len = i64::try_from(var.alt.len())?;
    let minimum = i64::from(refpadding);
    if ref_len <= minimum || alt_len <= minimum {
        return Ok(());
    }
    let reference_allele = reference_slice(reference, var.start, var.end)?.as_bytes();
    let alt = var.alt.as_bytes();
    while ref_len > minimum
        && alt_len > minimum
        && reference_allele[usize::try_from(ref_len - 1)?] == alt[usize::try_from(alt_len - 1)?]
    {
        ref_len -= 1;
        alt_len -= 1;
    }
    var.end = var.start + ref_len - 1;
    var.alt.truncate(usize::try_from(alt_len)?);
    Ok(())
}

fn left_shift_legacy(reference: &str, var: &mut RefVar) -> Result<()> {
    // C++ `leftShift(..., refpadding=false)` still calls trimLeft/Right
    // with their default `refpadding=true` internally.
    trim_left_legacy(reference, var, true)?;
    trim_right_legacy(reference, var, true)?;
    let initial_ref_len = var.end - var.start + 1;
    if initial_ref_len >= 0
        && initial_ref_len == i64::try_from(var.alt.len())?
        && reference_slice(reference, var.start, var.end)? == var.alt
    {
        return Ok(());
    }
    loop {
        if var.start <= 0 {
            break;
        }
        let mut ref_len = var.end - var.start + 1;
        let previous = reference.as_bytes()[usize::try_from(var.start - 1)?];
        if previous == b'N' {
            break;
        }
        let mut changed = false;
        if ref_len > 0 && !var.alt.is_empty() {
            let reference_last = reference.as_bytes()[usize::try_from(var.start + ref_len - 1)?];
            if reference_last == *var.alt.as_bytes().last().expect("non-empty ALT") {
                ref_len -= 1;
                var.end -= 1;
                var.alt.pop();
                changed = true;
            }
        }
        if ref_len == 0 || var.alt.is_empty() {
            var.start -= 1;
            var.alt.insert(0, char::from(previous));
            changed = true;
        }
        if !changed {
            break;
        }
    }
    trim_left_legacy(reference, var, true)?;
    trim_right_legacy(reference, var, true)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(
        pos: usize,
        reference: &str,
        alt: &str,
        truth_gt: &str,
        query_gt: &str,
    ) -> RawVcfRecord {
        RawVcfRecord {
            chrom: "chr1".to_string(),
            pos,
            id: ".".to_string(),
            ref_allele: reference.to_string(),
            alt_allele: alt.to_string(),
            qual: "12.3456789".to_string(),
            filter: "PASS".to_string(),
            info: ".".to_string(),
            format: Some("GT".to_string()),
            samples: vec![truth_gt.to_string(), query_gt.to_string()],
        }
    }

    fn references(length: usize) -> BTreeMap<String, String> {
        BTreeMap::from([("chr1".to_string(), "A".repeat(length))])
    }

    fn value(record: &RawVcfRecord, sample: usize, key: &str) -> String {
        record
            .sample_map(sample)
            .get(key)
            .cloned()
            .unwrap_or_default()
    }

    #[test]
    fn allele_mode_is_one_to_one_and_last_occurrence_wins() -> Result<()> {
        let headers =
            vec!["##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string()];

        let mut truth_homalt = vec![record(10, "A", "C", "1/1", "0/1")];
        annotate_merged_records(
            &mut truth_homalt,
            &headers,
            &references(100),
            ScmpMode::Alleles,
            "QUAL",
        )?;
        assert_eq!(value(&truth_homalt[0], 0, "BD"), "FN");
        assert_eq!(value(&truth_homalt[0], 1, "BD"), "TP");

        let mut query_homalt = vec![record(10, "A", "C", "0/1", "1/1")];
        annotate_merged_records(
            &mut query_homalt,
            &headers,
            &references(100),
            ScmpMode::Alleles,
            "QUAL",
        )?;
        assert_eq!(value(&query_homalt[0], 0, "BD"), "TP");
        assert_eq!(value(&query_homalt[0], 1, "BD"), "FP");
        Ok(())
    }

    #[test]
    fn allele_and_distance_modes_have_distinct_same_position_semantics() -> Result<()> {
        let headers = Vec::new();
        let mut allele = vec![record(10, "A", "C,G", "1/1", "2/2")];
        annotate_merged_records(
            &mut allele,
            &headers,
            &references(100),
            ScmpMode::Alleles,
            "QUAL",
        )?;
        assert_eq!(value(&allele[0], 0, "BD"), "FN");
        assert_eq!(value(&allele[0], 1, "BD"), "FP");

        let mut distance = vec![record(10, "A", "C,G", "1/1", "2/2")];
        annotate_merged_records(
            &mut distance,
            &headers,
            &references(100),
            ScmpMode::Distance { max_distance: 0 },
            "QUAL",
        )?;
        assert_eq!(value(&distance[0], 0, "BD"), "TP");
        assert_eq!(value(&distance[0], 1, "BD"), "TP");
        Ok(())
    }

    #[test]
    fn legacy_only_allele_mode_preserves_alt_length_refvar_bug() -> Result<()> {
        // Ordinary VCF normalization makes these homopolymer insertions
        // equivalent. Legacy SCMP instead builds RefVar.end from ALT length,
        // so C>CA at POS 1 and A>AA at POS 2 remain distinct.
        let mut records = vec![
            record(1, "C", "CA", "0/1", "0/0"),
            record(2, "A", "AA", "0/0", "0/1"),
        ];
        let reference = BTreeMap::from([("chr1".to_string(), "CAAAAA".to_string())]);
        annotate_merged_records(&mut records, &[], &reference, ScmpMode::Alleles, "QUAL")?;
        assert_eq!(value(&records[0], 0, "BD"), "FN");
        assert_eq!(value(&records[1], 1, "BD"), "FP");
        Ok(())
    }

    #[test]
    fn allele_mode_normalizes_lowercase_reference_before_matching() -> Result<()> {
        let mut records = vec![
            record(1, "C", "CACACAT", "0/1", "0/0"),
            record(7, "C", "T", "0/0", "0/1"),
        ];
        let reference = BTreeMap::from([("chr1".to_string(), "cacacac".to_string())]);

        annotate_merged_records(&mut records, &[], &reference, ScmpMode::Alleles, "QUAL")?;

        assert_eq!(value(&records[0], 0, "BD"), "TP");
        assert_eq!(value(&records[1], 1, "BD"), "TP");
        Ok(())
    }

    #[test]
    fn distance_boundary_is_inclusive_and_many_to_many() -> Result<()> {
        let mut at_boundary = vec![
            record(10, "A", "C", "0/1", "0/0"),
            record(40, "A", "G", "0/0", "1/1"),
        ];
        annotate_merged_records(
            &mut at_boundary,
            &[],
            &references(100),
            ScmpMode::Distance { max_distance: 30 },
            "QUAL",
        )?;
        assert_eq!(value(&at_boundary[0], 0, "BD"), "TP");
        assert_eq!(value(&at_boundary[1], 1, "BD"), "TP");

        at_boundary[1].pos = 41;
        at_boundary[0].format = Some("GT".to_string());
        at_boundary[0].samples = vec!["0/1".to_string(), "0/0".to_string()];
        at_boundary[1].format = Some("GT".to_string());
        at_boundary[1].samples = vec!["0/0".to_string(), "1/1".to_string()];
        annotate_merged_records(
            &mut at_boundary,
            &[],
            &references(100),
            ScmpMode::Distance { max_distance: 30 },
            "QUAL",
        )?;
        assert_eq!(value(&at_boundary[0], 0, "BD"), "FN");
        assert_eq!(value(&at_boundary[1], 1, "BD"), "FP");
        Ok(())
    }

    #[test]
    fn gap_flush_includes_the_boundary_pushing_record() -> Result<()> {
        let mut records = vec![
            record(1, "A", "C", "0/1", "0/0"),
            record(100_002, "A", "G", "0/0", "0/1"),
            record(100_003, "A", "T", "0/0", "0/1"),
        ];
        annotate_merged_records(
            &mut records,
            &[],
            &references(100_010),
            ScmpMode::Distance {
                max_distance: 200_000,
            },
            "QUAL",
        )?;
        assert_eq!(value(&records[0], 0, "BD"), "TP");
        assert_eq!(value(&records[1], 1, "BD"), "TP");
        assert_eq!(value(&records[2], 1, "BD"), "FP");
        Ok(())
    }

    #[test]
    fn qq_is_written_for_homref_while_bd_remains_absent() -> Result<()> {
        let mut records = vec![record(10, "A", "C", "0/0", "./.")];
        annotate_merged_records(
            &mut records,
            &[],
            &references(100),
            ScmpMode::Alleles,
            "QUAL",
        )?;
        assert_eq!(value(&records[0], 0, "QQ"), "12.345679");
        assert_eq!(value(&records[0], 1, "QQ"), "12.345679");
        assert!(!records[0].format_keys().contains(&"BD"));
        Ok(())
    }

    #[test]
    fn compatible_same_position_records_form_multiallelic_and_remap_fields() -> Result<()> {
        let headers = vec![
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Depth\">".to_string(),
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Likelihood\">".to_string(),
            "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">".to_string(),
            "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Allele number\">".to_string(),
            "##INFO=<ID=SAMPLE_GT,Number=.,Type=Integer,Description=\"Original GT\">".to_string(),
        ];
        let mut truth = record(10, "A", "C", "1/1", "0/0");
        truth.format = Some("GT:AD:PL".to_string());
        truth.samples = vec!["1/1:2,8:0,10,20".to_string()];
        truth.info = "SAMPLE_GT=4,4".to_string();
        let mut query = record(10, "A", "G", "0/0", "0/1");
        query.format = Some("GT:AD:PL".to_string());
        query.samples = vec!["0/1:3,7:0,11,22".to_string()];
        query.info = "SAMPLE_GT=4,4".to_string();

        let merged = merge_two_sample_records(&headers, &[truth], &headers, &[query])?;
        assert_eq!(merged.records.len(), 1);
        assert!(
            merged
                .headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=TRUTH_GT,"))
        );
        assert!(
            merged
                .headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=QUERY_GT,"))
        );
        let merged = &merged.records[0];
        assert_eq!(merged.alt_allele, "C,G");
        assert_eq!(value(merged, 0, "GT"), "1/1");
        assert_eq!(value(merged, 1, "GT"), "0/2");
        assert_eq!(value(merged, 0, "AD"), "2,8,.");
        assert_eq!(value(merged, 1, "AD"), "3,.,7");
        assert_eq!(value(merged, 0, "PL"), "0,10,20,.,.,.");
        assert_eq!(value(merged, 1, "PL"), "0,.,.,11,.,22");
        assert_eq!(merged.info, "TRUTH_GT=4,4;QUERY_GT=4,4;AN=4;AC=2,1");
        Ok(())
    }

    #[test]
    fn biallelic_merge_truncates_number_a_and_g_info_but_keeps_unbounded_values() {
        let headers = vec![
            "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency\">".to_string(),
            "##INFO=<ID=SAMPLE_AD,Number=.,Type=Integer,Description=\"Depth\">".to_string(),
            "##INFO=<ID=SAMPLE_PL,Number=G,Type=Integer,Description=\"Likelihood\">".to_string(),
        ];
        let mut source = record(10, "A", "C", "0/1", "0/0");
        source.samples.truncate(1);
        source.info = "AF=0.5,0.5;SAMPLE_PL=100,20,30,10,0,40;SAMPLE_AD=26,0,11".to_string();

        let merged = single_side_record(&source, 1, &BTreeMap::new(), &field_numbers(&headers));

        assert_eq!(merged.info, "QUERY_AD=26,0,11;AF=0.5;QUERY_PL=100,20,30");
    }

    #[test]
    fn merged_records_omit_undeclared_allele_counts() -> Result<()> {
        let headers = vec![
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "##INFO=<ID=SAMPLE_GT,Number=.,Type=Integer,Description=\"Original GT\">".to_string(),
        ];
        let mut truth = record(5, "G", "T", "0/1", "0/0");
        truth.samples.truncate(1);
        truth.info = "SAMPLE_GT=4,4".to_string();
        let query = truth.clone();

        let merged = merge_two_sample_records(&headers, &[truth], &headers, &[query])?;

        assert_eq!(merged.records[0].info, "TRUTH_GT=4,4;QUERY_GT=4,4");
        Ok(())
    }

    #[test]
    fn stale_biallelic_likelihoods_expand_with_missing_multiallelic_slots() {
        assert_eq!(
            remap_numbered("100,20,30,10,0,40", FieldNumber::G, &[0, 1], 2, "0/1",),
            "100,20,30,.,.,."
        );
    }

    #[test]
    fn duplicate_same_position_input_records_follow_bcftools_first_pairing() -> Result<()> {
        // Without an equivalent counterpart, bcftools 1.17
        // `merge --force-samples` pairs records in input order. Remaining
        // records are emitted with a missing genotype on the other side.
        let headers = vec![
            "##contig=<ID=chr1,length=100>".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        ];
        let mut truth_first = record(10, "A", "C", "0/1", "0/0");
        truth_first.id = "t1".to_string();
        truth_first.qual = "11".to_string();
        truth_first.samples.truncate(1);
        let mut truth_second = record(10, "A", "G", "0/1", "0/0");
        truth_second.id = "t2".to_string();
        truth_second.qual = "12".to_string();
        truth_second.samples.truncate(1);
        let mut query = record(10, "A", "T", "0/1", "0/0");
        query.id = "q1".to_string();
        query.qual = "13".to_string();
        query.samples.truncate(1);

        let merged =
            merge_two_sample_records(&headers, &[truth_first, truth_second], &headers, &[query])?;
        assert_eq!(merged.records.len(), 2);

        let paired = &merged.records[0];
        assert_eq!(paired.id, "t1;q1");
        assert_eq!(paired.alt_allele, "C,T");
        assert_eq!(paired.qual, "13");
        assert_eq!(value(paired, 0, "GT"), "0/1");
        assert_eq!(value(paired, 1, "GT"), "0/2");

        let remainder = &merged.records[1];
        assert_eq!(remainder.id, "t2");
        assert_eq!(remainder.alt_allele, "G");
        assert_eq!(remainder.qual, "12");
        assert_eq!(value(remainder, 0, "GT"), "0/1");
        assert_eq!(value(remainder, 1, "GT"), "./.");
        Ok(())
    }

    #[test]
    fn padded_equivalent_indel_is_paired_before_same_position_indel() -> Result<()> {
        let headers = vec![
            "##contig=<ID=chr1,length=100>".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        ];
        let mut truth_first = record(10, "ATT", "A", "0/1", "0/0");
        truth_first.id = "different".to_string();
        truth_first.samples.truncate(1);
        let mut truth_equivalent = record(10, "ATT", "AT", "0/1", "0/0");
        truth_equivalent.id = "equivalent".to_string();
        truth_equivalent.samples.truncate(1);
        let mut query = record(10, "AT", "A", "0/1", "0/0");
        query.id = "query".to_string();
        query.samples.truncate(1);

        let merged = merge_two_sample_records(
            &headers,
            &[truth_first, truth_equivalent],
            &headers,
            &[query],
        )?;

        assert_eq!(merged.records.len(), 2);
        assert_eq!(merged.records[0].id, "equivalent;query");
        assert_eq!(merged.records[0].ref_allele, "AT");
        assert_eq!(merged.records[0].alt_allele, "A");
        assert_eq!(merged.records[1].id, "different");
        assert_eq!(merged.records[1].ref_allele, "ATT");
        assert_eq!(merged.records[1].alt_allele, "A");
        Ok(())
    }

    #[test]
    fn legacy_only_duplicate_indels_preserve_bcftools_reused_count_order() -> Result<()> {
        let headers = vec![
            "##contig=<ID=chr1,length=100>".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        ];
        let mut truth_first = record(10, "AGT", "A", "0/1", "0/0");
        truth_first.id = "truth-first".to_string();
        truth_first.samples.truncate(1);
        let mut truth_middle = record(10, "AGT", "AGTGT", "0/1", "0/0");
        truth_middle.id = "truth-middle".to_string();
        truth_middle.samples.truncate(1);
        let mut truth_last = record(10, "AGT", "AGTGTGT", "0/1", "0/0");
        truth_last.id = "truth-last".to_string();
        truth_last.samples.truncate(1);
        let mut query_first = record(10, "AGT", "AGTGT", "0/1", "0/0");
        query_first.id = "query-first".to_string();
        query_first.samples.truncate(1);
        let mut query_last = record(10, "AGT", "A", "0/1", "0/0");
        query_last.id = "query-last".to_string();
        query_last.samples.truncate(1);

        let merged = merge_two_sample_records(
            &headers,
            &[truth_first, truth_middle, truth_last],
            &headers,
            &[query_first, query_last],
        )?;

        assert_eq!(merged.records.len(), 3);
        assert_eq!(merged.records[0].id, "truth-first;query-last");
        assert_eq!(merged.records[1].id, "truth-last;query-first");
        assert_eq!(merged.records[1].ref_allele, "A");
        assert_eq!(merged.records[1].alt_allele, "AGT,AGTGT");
        assert_eq!(merged.records[2].id, "truth-middle");
        Ok(())
    }

    #[test]
    fn merged_header_keeps_query_provenance_but_drops_input_identity() {
        let truth = vec![
            "##fileformat=VCFv4.1".to_string(),
            "##reference=truth.fa".to_string(),
            "##PEDIGREE=<Derived=truth-child,Original=truth-parent>".to_string(),
        ];
        let query = vec![
            "##fileformat=VCFv4.1".to_string(),
            "##reference=query.fa".to_string(),
            "##source=GATK".to_string(),
            "##PEDIGREE=<Derived=child,Original=parent>".to_string(),
            "##annotator=VariantAnnotator".to_string(),
        ];

        let merged = merged_headers(&truth, &query);
        assert_eq!(merged[0], "##fileformat=VCFv4.2");
        assert!(merged.iter().any(|line| line == "##reference=truth.fa"));
        assert!(
            merged
                .iter()
                .any(|line| line == "##annotator=VariantAnnotator")
        );
        assert!(!merged.iter().any(|line| line == "##reference=query.fa"));
        assert!(!merged.iter().any(|line| line == "##source=GATK"));
        assert!(!merged.iter().any(|line| line.starts_with("##PEDIGREE=")));
    }

    #[test]
    fn duplicate_alt_occurrences_are_preserved() -> Result<()> {
        let occurrences = extract_occurrences(&record(10, "A", "C", "1/1", "0/0"), 0, 0)?;
        assert_eq!(occurrences.len(), 2);
        assert_eq!(occurrences[0].var, occurrences[1].var);
        Ok(())
    }

    #[test]
    fn symbolic_del_is_supported_but_other_symbolic_alts_are_rejected() -> Result<()> {
        let deletion = extract_occurrences(&record(10, "AA", "<DEL>", "0/1", "0/0"), 0, 0)?;
        assert_eq!(
            deletion[0].var,
            RefVar {
                start: 9,
                end: 8,
                alt: String::new()
            }
        );
        assert!(extract_occurrences(&record(10, "A", "<INS>", "0/1", "0/0"), 0, 0).is_err());
        Ok(())
    }
}
