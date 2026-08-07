//! Cohesive quantify annotations responsibility.

use super::{BenchmarkSamples, RegionLevels};
use crate::adapters::report::{self, suffixed_report_path};
use crate::adapters::vcf;
use crate::domain::RawVcfRecord;
use crate::engines::roc;
use anyhow::{Context, Result, bail};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Write;
use std::path::Path;

#[derive(Debug, Eq, PartialEq)]
pub(super) struct Ga4ghAnnotation {
    pub(super) bi: String,
    pub(super) bvt: &'static str,
    pub(super) blt: &'static str,
}

/// Derive the fields that legacy `GA4GHQuantify` computes from each sample's
/// selected genotype alleles. RTG's GA4GH intermediate supplies BD/BK/QQ but
/// does not supply these count-oriented annotations.
pub(super) fn reannotate_ga4gh_record(record: &mut RawVcfRecord) {
    ensure_format_fields(record, &["BI", "BVT", "BLT", "QQ"]);
    let has_overwide_genotype = (0..record.samples.len()).any(|sample_index| {
        record
            .sample_map(sample_index)
            .get("GT")
            .is_some_and(|gt| gt.split(['/', '|']).count() > 2)
    });
    for sample_index in 0..record.samples.len() {
        let gt = record
            .sample_map(sample_index)
            .get("GT")
            .cloned()
            .unwrap_or_else(|| "./.".to_string());
        if has_overwide_genotype {
            set_format_value(record, sample_index, "BI", ".");
            set_format_value(record, sample_index, "BVT", "UNK");
            set_format_value(record, sample_index, "BLT", "ambi");
            let qq = if gt.split(['/', '|']).count() > 2 {
                "."
            } else {
                "0"
            };
            set_format_value(record, sample_index, "QQ", qq);
            continue;
        }
        let annotation = ga4gh_annotation(record, &gt);
        set_format_value(record, sample_index, "BI", &annotation.bi);
        set_format_value(record, sample_index, "BVT", annotation.bvt);
        set_format_value(record, sample_index, "BLT", annotation.blt);
    }
}

pub(super) fn ga4gh_annotation(record: &RawVcfRecord, gt: &str) -> Ga4ghAnnotation {
    let alleles = gt
        .split(['/', '|'])
        .filter(|part| !part.is_empty())
        .map(|part| part.parse::<usize>().ok())
        .collect::<Vec<_>>();
    let all_missing = alleles.is_empty() || alleles.iter().all(Option::is_none);
    let homref = !alleles.is_empty() && alleles.iter().all(|allele| *allele == Some(0));
    let blt = ga4gh_location_type(&alleles);

    let alts = record.alt_allele.split(',').collect::<Vec<_>>();
    let mut selected = alleles
        .iter()
        .filter_map(|allele| allele.as_ref().copied())
        .filter(|allele| *allele > 0)
        .collect::<BTreeSet<_>>();
    // VariantStatistics counts a homozygous alternate allele once. Using a
    // set is equivalent for BI/BVT and also de-duplicates repeated indexes in
    // polyploid input.
    let mut type_bits = 0u8;
    let mut extras = BTreeSet::<String>::new();
    let mut invalid_allele = false;
    for allele in std::mem::take(&mut selected) {
        let Some(alt) = alts.get(allele - 1) else {
            invalid_allele = true;
            continue;
        };
        let (bits, allele_extras) = ga4gh_allele_statistics(record, alt);
        type_bits |= bits;
        extras.extend(allele_extras);
    }

    const SNP: u8 = 1;
    const INS: u8 = 2;
    const DEL: u8 = 4;
    let bvt = if all_missing {
        "NOCALL"
    } else if homref {
        "HOMREF"
    } else if invalid_allele {
        "UNK"
    } else if type_bits == SNP {
        "SNP"
    } else if type_bits & (SNP | INS | DEL) != 0 {
        "INDEL"
    } else {
        "UNK"
    };
    Ga4ghAnnotation {
        bi: if extras.is_empty() {
            ".".to_string()
        } else {
            extras.into_iter().collect::<Vec<_>>().join(",")
        },
        bvt,
        blt,
    }
}

pub(super) fn ga4gh_location_type(alleles: &[Option<usize>]) -> &'static str {
    if alleles.len() > 2 {
        return "ambi";
    }
    if !alleles.is_empty() && alleles.iter().all(|allele| *allele == Some(0)) {
        return "homref";
    }
    if alleles.len() == 2 {
        match (alleles[0], alleles[1]) {
            (Some(0), Some(right)) | (Some(right), Some(0)) if right > 0 => "het",
            (Some(left), Some(right)) if left > 0 && left == right => "homalt",
            (Some(left), Some(right)) if left > 0 && right > 0 => "hetalt",
            (Some(_), None) | (None, Some(_)) => "halfcall",
            _ => "nocall",
        }
    } else if alleles.iter().all(Option::is_none) || alleles.is_empty() {
        "nocall"
    } else if alleles.len() == 1 {
        // GA4GH's legacy VariantStatistics categorises a one-allele call in
        // the same partial-call bucket as `./1`; it does not emit `hemi`.
        "halfcall"
    } else {
        "unknown"
    }
}

/// Return VariantStatistics' low type bits plus its lexically sorted BI set.
pub(super) fn ga4gh_allele_statistics(record: &RawVcfRecord, alt: &str) -> (u8, BTreeSet<String>) {
    const SNP: u8 = 1;
    const INS: u8 = 2;
    const DEL: u8 = 4;

    if alt.starts_with('<') && alt.ends_with('>') {
        let tag = alt[1..alt.len() - 1].to_ascii_uppercase();
        let size = record.ref_allele.len().max(1);
        let (bits, token) = if tag.starts_with("DEL") {
            (DEL, ga4gh_size_token('d', size))
        } else if tag.starts_with("INS") || tag.starts_with("DUP") {
            (INS, ga4gh_size_token('i', size))
        } else {
            (INS | DEL, ga4gh_size_token('c', size))
        };
        return (bits, BTreeSet::from([token]));
    }

    let primitives = crate::engines::align::realign_ref_var(
        record.pos,
        record.ref_allele.as_bytes(),
        alt.as_bytes(),
    );
    let mut total_snp = 0usize;
    let mut total_ins = 0usize;
    let mut total_del = 0usize;
    let mut ti = 0usize;
    let mut tv = 0usize;
    for primitive in primitives {
        if primitive.end < primitive.start {
            total_ins += primitive.alt.len();
        } else if primitive.alt.is_empty() {
            total_del += primitive.end - primitive.start + 1;
        } else if primitive.end == primitive.start && primitive.alt.len() == 1 {
            total_snp += 1;
            let ref_offset = primitive.start.saturating_sub(record.pos);
            let ref_base = record.ref_allele.as_bytes().get(ref_offset).copied();
            let alt_base = primitive.alt.as_bytes().first().copied();
            if ref_base
                .zip(alt_base)
                .is_some_and(|(reference, alternate)| {
                    is_transition_pair(reference as char, alternate as char)
                })
            {
                ti += 1;
            } else {
                tv += 1;
            }
        } else {
            // The shared aligner normally decomposes every concrete allele to
            // SNP/INS/DEL primitives. Retain an INDEL classification if a
            // future representation reaches this fallback.
            total_del += primitive.end - primitive.start + 1;
            total_ins += primitive.alt.len();
        }
    }

    let mut extras = BTreeSet::new();
    if ti > 0 {
        extras.insert("ti".to_string());
    }
    if tv > 0 {
        extras.insert("tv".to_string());
    }
    if total_snp == 0 {
        if total_ins > 0 {
            extras.insert(ga4gh_size_token('i', total_ins));
        }
        if total_del > 0 {
            extras.insert(ga4gh_size_token('d', total_del));
        }
    } else if total_ins + total_del > 0 {
        extras.insert(ga4gh_size_token('c', total_ins + total_del));
    }

    let mut bits = 0u8;
    if total_snp > 0 {
        bits |= SNP;
    }
    if total_ins > 0 {
        bits |= INS;
    }
    if total_del > 0 {
        bits |= DEL;
    }
    (bits, extras)
}

pub(super) fn ga4gh_size_token(prefix: char, size: usize) -> String {
    let bucket = match size {
        0..=5 => "1_5",
        6..=15 => "6_15",
        _ => "16_plus",
    };
    format!("{prefix}{bucket}")
}

pub(super) fn is_transition_pair(reference: char, alternate: char) -> bool {
    matches!(
        (
            reference.to_ascii_uppercase(),
            alternate.to_ascii_uppercase()
        ),
        ('A', 'G') | ('G', 'A') | ('C', 'T') | ('T', 'C')
    )
}

pub(super) fn ensure_format_fields(record: &mut RawVcfRecord, additions: &[&str]) {
    let mut keys = record
        .format
        .as_deref()
        .filter(|format| !format.is_empty() && *format != ".")
        .map(|format| format.split(':').map(str::to_string).collect::<Vec<_>>())
        .unwrap_or_default();
    for sample in &mut record.samples {
        let mut values = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        values.resize(keys.len(), ".".to_string());
        *sample = values.join(":");
    }
    for addition in additions {
        if keys.iter().any(|key| key == addition) {
            continue;
        }
        keys.push((*addition).to_string());
        for sample in &mut record.samples {
            sample.push(':');
            sample.push('.');
        }
    }
    record.format = (!keys.is_empty()).then(|| keys.join(":"));
}

pub(super) fn ensure_ga4gh_headers(headers: &mut Vec<String>) {
    const REQUIRED: [(&str, &str); 8] = [
        (
            "INFO=<ID=BS,",
            "##INFO=<ID=BS,Number=.,Type=Integer,Description=\"Benchmarking superlocus ID for these variants.\">",
        ),
        (
            "FORMAT=<ID=GT,",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        ),
        (
            "FORMAT=<ID=BD,",
            "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">",
        ),
        (
            "FORMAT=<ID=BK,",
            "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">",
        ),
        (
            "FORMAT=<ID=BI,",
            "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">",
        ),
        (
            "FORMAT=<ID=QQ,",
            "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">",
        ),
        (
            "FORMAT=<ID=BVT,",
            "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"High-level variant type (SNP|INDEL).\">",
        ),
        (
            "FORMAT=<ID=BLT,",
            "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"High-level location type (het|homref|hetalt|homalt|nocall).\">",
        ),
    ];
    let insertion_index = headers
        .iter()
        .position(|header| header.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    let mut offset = 0usize;
    for (needle, declaration) in REQUIRED {
        if headers.iter().any(|header| header.contains(needle)) {
            continue;
        }
        headers.insert(insertion_index + offset, declaration.to_string());
        offset += 1;
    }
}

pub(super) fn validate_ga4gh_qq_fields(headers: &[String], records: &[RawVcfRecord]) -> Result<()> {
    let qq_is_string = headers.iter().any(|header| {
        let Some(body) = header.strip_prefix("##FORMAT=<") else {
            return false;
        };
        let fields = body.split(',').collect::<Vec<_>>();
        fields.contains(&"ID=QQ") && fields.contains(&"Type=String")
    });
    if !qq_is_string {
        return Ok(());
    }

    // Legacy's numeric FORMAT reader sizes its result from the encoded BCF
    // string width. A one-character String QQ is accepted as a missing numeric
    // value, while a wider cell yields multiple values and aborts before any
    // report is published (`BCFHelpers.cpp::getFormatFloat`).
    for record in records {
        let Some(qq_index) = record.format_keys().iter().position(|field| *field == "QQ") else {
            continue;
        };
        let encoded_width = record
            .samples
            .iter()
            .filter_map(|sample| sample.split(':').nth(qq_index))
            .map(str::len)
            .max()
            .unwrap_or(0);
        if encoded_width > 1 {
            bail!(
                "too many QQ fields at {}:{}",
                record.chrom,
                record.pos.saturating_sub(1)
            );
        }
    }
    Ok(())
}

pub(super) fn canonicalize_ga4gh_header_order(headers: &mut Vec<String>) {
    let Some(pass_index) = headers
        .iter()
        .position(|header| header.starts_with("##FILTER=<ID=PASS,"))
    else {
        return;
    };
    let pass = headers.remove(pass_index);
    let insertion = headers
        .iter()
        .position(|header| !header.starts_with("##fileformat="))
        .unwrap_or(headers.len());
    headers.insert(insertion, pass);
}

pub(super) fn ensure_info_header(headers: &mut Vec<String>, id: &str, declaration: &str) {
    if headers
        .iter()
        .any(|header| header.contains(&format!("INFO=<ID={id},")))
    {
        return;
    }
    let index = headers
        .iter()
        .position(|header| header.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    headers.insert(index, declaration.to_string());
}

pub(super) fn propagate_ga4gh_superlocus_for_samples(
    records: &mut [RawVcfRecord],
    samples: BenchmarkSamples,
    preserve_missing_query_qq: bool,
    inherit_same_position_tp_qq: bool,
) {
    let mut same_position_tp_qq = BTreeMap::<usize, String>::new();
    if let Some(query_index) = samples.query {
        for record in records.iter() {
            let query = record.sample_map(query_index);
            let Some(score) = (query.get("BD").map(String::as_str) == Some("TP"))
                .then(|| query.get("QQ").cloned())
                .flatten()
                .filter(|score| {
                    score != "." && score.parse::<f64>().is_ok_and(|value| value.is_finite())
                })
            else {
                continue;
            };
            same_position_tp_qq
                .entry(record.pos)
                .and_modify(|current| {
                    if score.parse::<f64>().unwrap() < current.parse::<f64>().unwrap() {
                        *current = score.clone();
                    }
                })
                .or_insert(score);
        }
    }
    let minimum_tp_qq = records
        .iter()
        .filter_map(|record| {
            let query = record.sample_map(samples.query?);
            (query.get("BD").map(String::as_str) == Some("TP"))
                .then(|| query.get("QQ").cloned())
                .flatten()
        })
        .filter(|score| score != "." && score.parse::<f64>().is_ok_and(|value| value.is_finite()))
        .min_by(|left, right| {
            left.parse::<f64>()
                .unwrap()
                .total_cmp(&right.parse::<f64>().unwrap())
        });

    let block_filters = records
        .iter()
        .flat_map(|record| record.filter.split(';'))
        .filter(|filter| !filter.is_empty() && *filter != "." && *filter != "PASS")
        .map(str::to_string)
        .collect::<BTreeSet<_>>();

    for record in records {
        let truth_tp = samples.truth.is_some_and(|truth| {
            record.sample_map(truth).get("BD").map(String::as_str) == Some("TP")
        });
        if let Some(query) = samples.query {
            let query_sample = record.sample_map(query);
            let query_qq = query_sample.get("QQ").cloned();
            if !preserve_missing_query_qq && query_qq.as_deref().is_none_or(|value| value == ".") {
                let inherited = (truth_tp && inherit_same_position_tp_qq)
                    .then(|| same_position_tp_qq.get(&record.pos))
                    .flatten()
                    .map(String::as_str)
                    .unwrap_or("0");
                set_format_value(record, query, "QQ", inherited);
            }
        }
        let query = samples
            .query
            .map(|query| record.sample_map(query))
            .unwrap_or_default();
        let direct_query_qq = (query.get("BD").map(String::as_str) == Some("TP"))
            .then(|| query.get("QQ").cloned())
            .flatten()
            .filter(|score| {
                score != "." && score.parse::<f64>().is_ok_and(|value| value.is_finite())
            });
        let truth_qq = if truth_tp {
            direct_query_qq.as_ref().or(minimum_tp_qq.as_ref())
        } else {
            None
        };
        if let Some(truth) = samples.truth {
            set_format_value(
                record,
                truth,
                "QQ",
                truth_qq.map(String::as_str).unwrap_or("."),
            );
        }

        if truth_tp
            && query.get("BVT").map(String::as_str) == Some("NOCALL")
            && !block_filters.is_empty()
        {
            let mut filters = record
                .filter
                .split(';')
                .filter(|filter| !filter.is_empty() && *filter != "." && *filter != "PASS")
                .map(str::to_string)
                .collect::<BTreeSet<_>>();
            filters.extend(block_filters.iter().cloned());
            record.filter = filters.into_iter().collect::<Vec<_>>().join(";");
        }
    }
}

pub(super) fn normalize_integer_like_format_values(record: &mut RawVcfRecord, key: &str) {
    let Some(index) = record.format_keys().iter().position(|field| *field == key) else {
        return;
    };
    for sample in &mut record.samples {
        let mut fields = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        let Some(value) = fields.get_mut(index) else {
            continue;
        };
        if let Ok(number) = value.parse::<f64>()
            && number.is_finite()
            && number.fract() == 0.0
        {
            *value = number.to_string();
            *sample = fields.join(":");
        }
    }
}

/// Legacy `fastainfo` reports the contig length after trimming only terminal
/// N-runs. Internal N-runs remain part of `Subset.Size`.
pub(super) fn n_trimmed_length(sequence: &str) -> usize {
    let bytes = sequence.as_bytes();
    let leading = bytes
        .iter()
        .take_while(|byte| matches!(**byte, b'N' | b'n'))
        .count();
    let trailing = bytes
        .iter()
        .rev()
        .take_while(|byte| matches!(**byte, b'N' | b'n'))
        .count();
    bytes.len().saturating_sub(leading + trailing)
}

pub(super) fn legacy_regions_extent(record: &RawVcfRecord) -> String {
    effective_reference_range(record)
        .map(|(start, end, _)| format!("{start}-{end}"))
        .unwrap_or_else(|| format!("{}-{}", record.pos, record.end_pos()))
}

/// Port `QuantifyRegions::annotate`'s per-allele reference span. Coordinates
/// are 1-based inclusive; pure insertions bracket both adjacent reference
/// bases, while the pinned legacy qfy assigns the insertion using its left
/// VCF anchor.
pub(super) fn effective_reference_range(record: &RawVcfRecord) -> Option<(usize, usize, bool)> {
    let ref_bytes = record.ref_allele.as_bytes();
    let pos_0b = record.pos.saturating_sub(1) as i64;
    let mut updated_start = i64::MAX;
    let mut updated_end = i64::MIN;
    let mut pure_insertion = false;
    let mut has_nucleotide_alt = false;

    for alt in record.alt_allele.split(',') {
        // classifyAlleleString maps a missing ALT to an empty allele and
        // processes it as a reference-consuming deletion. Any other
        // non-nucleotide allele terminates the loop, so later nucleotide
        // ALTs must not expand the range.
        let normalized_alt = if alt.is_empty() || alt == "." {
            String::new()
        } else if alt.bytes().all(is_legacy_nucleotide_base) {
            alt.to_ascii_uppercase()
        } else {
            break;
        };
        let alt_bytes = normalized_alt.as_bytes();
        let mut ref_len = ref_bytes.len();
        let mut alt_len = alt_bytes.len();
        while ref_len > 0 && alt_len > 0 && ref_bytes[ref_len - 1] == alt_bytes[alt_len - 1] {
            ref_len -= 1;
            alt_len -= 1;
        }
        let mut prefix = 0usize;
        while prefix < ref_len && prefix < alt_len && ref_bytes[prefix] == alt_bytes[prefix] {
            prefix += 1;
        }
        let start = pos_0b + prefix as i64;
        let end = pos_0b + ref_len as i64 - 1;
        if !has_nucleotide_alt {
            pure_insertion = true;
        }
        has_nucleotide_alt = true;
        if end >= start {
            updated_start = updated_start.min(start);
            updated_end = updated_end.max(end);
            pure_insertion = false;
        } else {
            updated_start = updated_start.min(start - 1);
            updated_end = updated_end.max(start);
        }
    }

    if !has_nucleotide_alt {
        return None;
    }
    Some((
        (updated_start + 1) as usize,
        (updated_end + 1) as usize,
        pure_insertion,
    ))
}

pub(super) fn is_legacy_nucleotide_base(base: u8) -> bool {
    matches!(
        base.to_ascii_uppercase(),
        b'A' | b'C'
            | b'G'
            | b'T'
            | b'U'
            | b'R'
            | b'Y'
            | b'K'
            | b'M'
            | b'S'
            | b'W'
            | b'B'
            | b'D'
            | b'H'
            | b'V'
            | b'N'
            | b'X'
    )
}

pub(super) fn remove_region_tag(info: &mut String, unwanted: &str) {
    let mut fields = Vec::new();
    for field in info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
    {
        if let Some(regions) = field.strip_prefix("Regions=") {
            let retained = regions
                .split(',')
                .filter(|region| !region.is_empty() && *region != unwanted)
                .collect::<Vec<_>>();
            if !retained.is_empty() {
                fields.push(format!("Regions={}", retained.join(",")));
            }
        } else {
            fields.push(field.to_string());
        }
    }
    *info = if fields.is_empty() {
        ".".to_string()
    } else {
        fields.join(";")
    };
}

pub(super) fn move_region_to_front(info: &mut String, wanted: &str) {
    let mut fields = Vec::new();
    for field in info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
    {
        if let Some(regions) = field.strip_prefix("Regions=") {
            let mut reordered = regions
                .split(',')
                .filter(|region| !region.is_empty() && *region != wanted)
                .map(str::to_string)
                .collect::<Vec<_>>();
            reordered.insert(0, wanted.to_string());
            fields.push(format!("Regions={}", reordered.join(",")));
        } else {
            fields.push(field.to_string());
        }
    }
    *info = if fields.is_empty() {
        ".".to_string()
    } else {
        fields.join(";")
    };
}

/// Reproduce the old XCMP quantifier rather than consuming GA4GH `BD` fields.
/// XCMP decisions live in record-level `INFO/type`, `kind`, and `ctype`; a
/// finalized hap.py VCF intentionally lacks them, which is why the legacy qfy
/// lane ignores calls inside CONF and labels calls outside CONF as UNK.
#[cfg(test)]
pub(super) fn reannotate_xcmp_record(
    record: &mut RawVcfRecord,
    has_confidence_regions: bool,
    roc_field: &str,
) {
    reannotate_xcmp_record_for_samples(
        record,
        has_confidence_regions,
        roc_field,
        BenchmarkSamples::POSITIONAL,
    );
}

pub(super) fn reannotate_xcmp_record_for_samples(
    record: &mut RawVcfRecord,
    has_confidence_regions: bool,
    roc_field: &str,
    samples: BenchmarkSamples,
) {
    let mut decision = info_value(&record.info, "type").unwrap_or_default();
    let mismatch_kind = info_value(&record.info, "kind").unwrap_or_default();
    let comparison_type = info_value(&record.info, "ctype").unwrap_or_default();
    let hap_match = info_flag(&record.info, "HapMatch");
    let import_fail = info_flag(&record.info, "IMPORT_FAIL");
    let query_filtered = info_flag(&record.info, "Q_FILTERED");

    if hap_match && decision != "TP" && !query_filtered {
        decision = "TP".to_string();
    }
    if has_confidence_regions && !has_region(&record.info, "CONF") {
        decision = "UNK".to_string();
    }
    if import_fail {
        decision = "N".to_string();
    }

    let match_kind = if decision == "TP" {
        "gm"
    } else if mismatch_kind == "gtmismatch" {
        "am"
    } else if mismatch_kind == "almismatch" || comparison_type == "hap:mismatch" {
        "lm"
    } else {
        "."
    };

    let record_qual = record.qual.clone();
    for sample_index in 0..record.samples.len() {
        let fields = record.sample_map(sample_index);
        let selected_roc = if roc_field == "QUAL" {
            Some(record_qual.clone())
        } else if let Some(value) = info_value(&record.info, roc_field) {
            Some(value.split(',').next().unwrap_or(".").to_string())
        } else {
            fields.get(roc_field).cloned()
        }
        .unwrap_or_else(|| ".".to_string());
        let gt = fields.get("GT").map(String::as_str).unwrap_or("./.");
        let no_call = gt
            .split(['/', '|'])
            .all(|allele| allele.is_empty() || allele == ".");
        let suppressed =
            import_fail || no_call || (samples.query == Some(sample_index) && query_filtered);
        let sample_decision = if import_fail {
            "N"
        } else if suppressed || decision.is_empty() {
            "."
        } else if samples.truth == Some(sample_index) && decision == "FP" {
            "FN"
        } else {
            decision.as_str()
        };
        set_format_value(record, sample_index, "BD", sample_decision);
        set_format_value(
            record,
            sample_index,
            "BK",
            if suppressed { "." } else { match_kind },
        );
        set_format_value(
            record,
            sample_index,
            "QQ",
            if suppressed {
                "0"
            } else {
                selected_roc.as_str()
            },
        );
    }
}

#[cfg(test)]
pub(super) fn decorate_quantified_record(
    record: &mut RawVcfRecord,
    annotation_type: &str,
    preserve_info: bool,
    output_vtc: bool,
    has_confidence_regions: bool,
) {
    decorate_quantified_record_for_samples(
        record,
        annotation_type,
        preserve_info,
        output_vtc,
        has_confidence_regions,
        BenchmarkSamples::POSITIONAL,
    );
}

pub(super) fn decorate_quantified_record_for_samples(
    record: &mut RawVcfRecord,
    annotation_type: &str,
    preserve_info: bool,
    output_vtc: bool,
    has_confidence_regions: bool,
    samples: BenchmarkSamples,
) {
    if output_vtc {
        if annotation_type == "xcmp" {
            let mut decision = info_value(&record.info, "type").unwrap_or_default();
            let mut kind = info_value(&record.info, "kind").unwrap_or_default();
            let ctype = info_value(&record.info, "ctype").unwrap_or_default();
            let query_filtered = info_flag(&record.info, "Q_FILTERED");
            if info_flag(&record.info, "HapMatch") && decision != "TP" && !query_filtered {
                kind = format!("hapmatch__{decision}__{kind}");
                decision = "TP".to_string();
            }
            if has_confidence_regions && !has_region(&record.info, "CONF") {
                decision = "UNK".to_string();
            }
            if info_flag(&record.info, "IMPORT_FAIL") {
                decision = "N".to_string();
                kind = "error".to_string();
            }
            let gtt1 = info_value(&record.info, "gtt1").unwrap_or_else(|| ".".to_string());
            let gtt2 = info_value(&record.info, "gtt2").unwrap_or_else(|| ".".to_string());
            set_info_value(
                &mut record.info,
                "XCMP",
                &format!("{decision}:{kind}:{gtt1}:{gtt2}:{ctype}"),
            );
        }
        let truth = samples
            .truth
            .map(|sample_index| record.sample_map(sample_index))
            .unwrap_or_default();
        let query = samples
            .query
            .map(|sample_index| record.sample_map(sample_index))
            .unwrap_or_default();
        let vtc = legacy_vtc(record, &truth, &query);
        if !vtc.is_empty() {
            set_info_value(&mut record.info, "VTC", &vtc);
        }
        if info_value(&record.info, "Regions").as_deref() == Some("TS_boundary") {
            move_info_field_to_end(&mut record.info, "Regions");
        }
    }

    if annotation_type == "xcmp" && !preserve_info {
        retain_info_fields(
            &mut record.info,
            &["END", "VTC", "Regions", "BS", "XCMP", "IMPORT_FAIL"],
        );
    }
}

pub(super) fn set_info_value(info: &mut String, key: &str, value: &str) {
    let mut fields = info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .filter(|field| field.split_once('=').map_or(*field, |(name, _)| name) != key)
        .map(str::to_string)
        .collect::<Vec<_>>();
    fields.push(format!("{key}={value}"));
    *info = fields.join(";");
}

pub(super) fn retain_info_fields(info: &mut String, retained: &[&str]) {
    let fields = info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .filter(|field| {
            let key = field.split_once('=').map_or(*field, |(key, _)| key);
            retained.contains(&key)
        })
        .collect::<Vec<_>>();
    *info = if fields.is_empty() {
        ".".to_string()
    } else {
        fields.join(";")
    };
}

pub(super) fn move_info_field_to_end(info: &mut String, wanted: &str) {
    let mut fields = info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .map(str::to_string)
        .collect::<Vec<_>>();
    let Some(index) = fields
        .iter()
        .position(|field| field.split_once('=').map_or(field.as_str(), |(key, _)| key) == wanted)
    else {
        return;
    };
    let field = fields.remove(index);
    fields.push(field);
    *info = fields.join(";");
}

pub(super) fn legacy_vtc(
    record: &RawVcfRecord,
    truth: &BTreeMap<String, String>,
    query: &BTreeMap<String, String>,
) -> String {
    let mut types = BTreeMap::<u8, String>::new();
    for sample in [truth, query] {
        if sample
            .get("BVT")
            .is_none_or(|value| matches!(value.as_str(), "" | "." | "NOCALL"))
        {
            types.insert(0x80, "nocall__nc".to_string());
            continue;
        }
        let alleles = sample
            .get("GT")
            .map(String::as_str)
            .unwrap_or(".")
            .split(['/', '|'])
            .filter_map(|allele| allele.parse::<usize>().ok())
            .collect::<Vec<_>>();
        let mut combined = 0u8;
        for allele in alleles.iter().copied().filter(|allele| *allele > 0) {
            let Some(alternate) = record.alt_allele.split(',').nth(allele - 1) else {
                continue;
            };
            let bits = allele_edit_bits(&record.ref_allele, alternate);
            combined |= bits;
            for bit in [1u8, 2, 4] {
                if bits & bit != 0 {
                    types.insert(bit, format!("nuc__{}", legacy_type_bits(bit)));
                }
            }
            if bits != 0 {
                types.insert(0x10 | bits, format!("al__{}", legacy_type_bits(bits)));
            }
        }
        if combined == 0 {
            continue;
        }
        let location = match sample.get("BLT").map(String::as_str).unwrap_or("") {
            "het" => 0x30,
            "hetalt" => 0x40,
            "hemi" => 0x50,
            "homalt" => 0x90,
            _ => 0xa0,
        };
        let ref_bit = u8::from(alleles.contains(&0)) * 8;
        types.insert(
            location | ref_bit | combined,
            format!(
                "{}__{}",
                sample.get("BLT").map(String::as_str).unwrap_or("unknown"),
                legacy_type_bits(ref_bit | combined)
            ),
        );
    }
    types.into_values().collect::<Vec<_>>().join(",")
}

pub(super) fn allele_edit_bits(reference: &str, alternate: &str) -> u8 {
    if alternate.starts_with('<') {
        return if alternate.starts_with("<DEL") { 4 } else { 2 };
    }
    let ref_bytes = reference.as_bytes();
    let alt_bytes = alternate.as_bytes();
    let prefix = ref_bytes
        .iter()
        .zip(alt_bytes)
        .take_while(|(left, right)| left == right)
        .count();
    let suffix_limit = (ref_bytes.len() - prefix).min(alt_bytes.len() - prefix);
    let suffix = (0..suffix_limit)
        .take_while(|offset| {
            ref_bytes[ref_bytes.len() - 1 - offset] == alt_bytes[alt_bytes.len() - 1 - offset]
        })
        .count();
    let ref_remaining = ref_bytes.len() - prefix - suffix;
    let alt_remaining = alt_bytes.len() - prefix - suffix;
    match (ref_remaining, alt_remaining) {
        (0, 0) => 0,
        (0, _) => 2,
        (_, 0) => 4,
        (left, right) if left == right => 1,
        (left, right) if left < right => 1 | 2,
        _ => 1 | 4,
    }
}

pub(super) fn legacy_type_bits(bits: u8) -> &'static str {
    const NAMES: [&str; 16] = [
        "nc", "s", "i", "si", "d", "sd", "id", "sid", "r", "rs", "ri", "rsi", "rd", "rsd", "rid",
        "rsid",
    ];
    NAMES[usize::from(bits & 0x0f)]
}

pub(super) fn info_value(info: &str, key: &str) -> Option<String> {
    info.split(';').find_map(|field| {
        field
            .split_once('=')
            .filter(|(field_key, _)| *field_key == key)
            .map(|(_, value)| value.to_string())
    })
}

pub(super) fn info_flag(info: &str, key: &str) -> bool {
    info.split(';').any(|field| field == key)
}

pub(super) fn has_region(info: &str, wanted: &str) -> bool {
    info.split(';')
        .find_map(|field| field.strip_prefix("Regions="))
        .is_some_and(|regions| regions.split(',').any(|region| region == wanted))
}

pub(super) fn set_format_value(
    record: &mut RawVcfRecord,
    sample_index: usize,
    key: &str,
    value: &str,
) {
    let Some(index) = record.format_keys().iter().position(|field| *field == key) else {
        return;
    };
    let Some(sample) = record.samples.get_mut(sample_index) else {
        return;
    };
    let mut fields = sample.split(':').map(str::to_string).collect::<Vec<_>>();
    if let Some(field) = fields.get_mut(index) {
        *field = value.to_string();
        *sample = fields.join(":");
    }
}

pub(super) fn merge_region_tags(info: &mut String, additions: &[String]) {
    let mut tags = Vec::<String>::new();
    let mut fields = Vec::<String>::new();
    for field in info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
    {
        if let Some(regions) = field.strip_prefix("Regions=") {
            tags.extend(
                regions
                    .split(',')
                    .filter(|tag| !tag.is_empty())
                    .map(str::to_string),
            );
        } else {
            fields.push(field.to_string());
        }
    }
    for addition in additions {
        if !tags.contains(addition) {
            tags.push(addition.clone());
        }
    }
    let regions = format!("Regions={}", tags.join(","));
    let insertion = fields
        .iter()
        .position(|field| field.starts_with("RegionsExtent="))
        .unwrap_or(fields.len());
    fields.insert(insertion, regions);
    *info = fields.join(";");
}

pub(super) fn compact_no_roc_outputs(prefix: &Path) -> Result<()> {
    let all_path = suffixed_report_path(prefix, "roc.all.csv.gz");
    let text = vcf::read_text(&all_path)?;
    let rows = text
        .lines()
        .enumerate()
        .filter(|(index, line)| *index == 0 || line.split(',').nth(6).is_some_and(|qq| qq == "*"))
        .map(|(_, line)| line)
        .collect::<Vec<_>>();
    let file = fs::File::create(&all_path)
        .with_context(|| format!("failed to create {}", all_path.display()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    writeln!(encoder, "{}", rows.join("\n"))?;
    encoder
        .finish()
        .with_context(|| format!("failed to finish ROC artifact {}", all_path.display()))?;

    for suffix in [
        "roc.Locations.SNP.csv.gz",
        "roc.Locations.SNP.PASS.csv.gz",
        "roc.Locations.INDEL.csv.gz",
        "roc.Locations.INDEL.PASS.csv.gz",
    ] {
        let path = suffixed_report_path(prefix, suffix);
        if path.exists() {
            fs::remove_file(&path)
                .with_context(|| format!("failed to remove {}", path.display()))?;
        }
    }
    Ok(())
}

pub(super) fn apply_stratification_levels(
    prefix: &Path,
    levels: &RegionLevels,
    preserve_raw_table: bool,
) -> Result<()> {
    if !levels.values().any(|level| *level > 0) {
        return Ok(());
    }

    let csv_path = suffixed_report_path(prefix, "roc.all.csv.gz");
    let csv = vcf::read_text(&csv_path)?;
    let csv = rewrite_subset_levels(&csv, ',', levels);
    let csv_file = fs::File::create(&csv_path)
        .with_context(|| format!("failed to create {}", csv_path.display()))?;
    let mut encoder = GzEncoder::new(csv_file, Compression::default());
    encoder.write_all(csv.as_bytes())?;
    encoder
        .finish()
        .with_context(|| format!("failed to finish ROC artifact {}", csv_path.display()))?;

    if preserve_raw_table {
        let raw_path = suffixed_report_path(prefix, "roc.tsv");
        if raw_path.exists() {
            let raw = fs::read_to_string(&raw_path)
                .with_context(|| format!("failed to read {}", raw_path.display()))?;
            fs::write(&raw_path, rewrite_subset_levels(&raw, '\t', levels))
                .with_context(|| format!("failed to write {}", raw_path.display()))?;
        }
    }
    Ok(())
}

pub(super) fn rewrite_subset_levels(text: &str, delimiter: char, levels: &RegionLevels) -> String {
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return String::new();
    };
    let header_fields = header.split(delimiter).collect::<Vec<_>>();
    let Some(subset_index) = header_fields.iter().position(|field| *field == "Subset") else {
        return text.to_string();
    };
    let Some(level_index) = header_fields
        .iter()
        .position(|field| *field == "Subset.Level")
    else {
        return text.to_string();
    };

    let separator = delimiter.to_string();
    let mut output = vec![header.to_string()];
    for line in lines {
        let mut fields = line
            .split(delimiter)
            .map(str::to_string)
            .collect::<Vec<_>>();
        if let Some(level) = fields
            .get(subset_index)
            .and_then(|subset| levels.get(subset))
            && let Some(field) = fields.get_mut(level_index)
        {
            *field = format!("{:.6}", *level as f64);
        }
        output.push(fields.join(&separator));
    }
    output.join("\n") + "\n"
}
