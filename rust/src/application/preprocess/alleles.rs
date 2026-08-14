//! Cohesive preprocessing responsibility.

use super::genotype::{bcf_encoded_gt, project_split_ad, project_split_genotype, remap_gt};
use super::options::SymbolicDeletionMaterialization;
use crate::application::{PreprocessArgs, SomaticGtMode};
use crate::domain::RawVcfRecord;
use anyhow::{Result, bail};
use std::collections::BTreeSet;

pub(super) fn materialize_symbolic_deletion(
    record: &mut RawVcfRecord,
    end: usize,
    reference: &[u8],
) -> Result<SymbolicDeletionMaterialization> {
    if end < record.pos {
        bail!(
            "symbolic deletion {}:{} has END={} before POS",
            record.chrom,
            record.pos,
            end
        );
    }

    if record.pos == 1 {
        // VariantWriter cannot add a preceding anchor at the start of a
        // contig. The pinned legacy binary instead expands REF through END,
        // keeps `<DEL>` symbolic, and consumes END from INFO.
        let ref_bases = reference.get(..end).ok_or_else(|| {
            anyhow::anyhow!(
                "symbolic deletion {}:1-{} extends beyond the reference sequence",
                record.chrom,
                end
            )
        })?;
        record.ref_allele = String::from_utf8_lossy(ref_bases).to_ascii_uppercase();
        record.alt_allele = record.alt_allele.to_ascii_uppercase();
        record.info = remove_info_field(&record.info, "END");
        return Ok(SymbolicDeletionMaterialization::ContigStart);
    }

    let (start, stop, alt_index) = {
        // Normal case: use the base immediately before the deleted interval
        // as the shared VCF anchor. `end` is 1-based inclusive, and therefore
        // also the exclusive byte index for this zero-based slice.
        (record.pos - 2, end, record.pos - 2)
    };
    let ref_bases = reference.get(start..stop).ok_or_else(|| {
        anyhow::anyhow!(
            "symbolic deletion {}:{}-{} extends beyond the reference sequence",
            record.chrom,
            record.pos,
            end
        )
    })?;
    let anchor = reference.get(alt_index..alt_index + 1).ok_or_else(|| {
        anyhow::anyhow!(
            "symbolic deletion {}:{}-{} has no reference anchor",
            record.chrom,
            record.pos,
            end
        )
    })?;

    let anchor = String::from_utf8_lossy(anchor).to_ascii_uppercase();
    let materialized_alts = record
        .alt_allele
        .split(',')
        .map(|alt| {
            if alt == "<DEL>" {
                anchor.clone()
            } else if is_symbolic_allele(alt) {
                alt.to_ascii_uppercase()
            } else {
                format!("{anchor}{}", alt.to_ascii_uppercase())
            }
        })
        .collect::<Vec<_>>();

    record.pos -= 1;
    record.ref_allele = String::from_utf8_lossy(ref_bases).to_ascii_uppercase();
    record.alt_allele = materialized_alts.join(",");
    record.info = remove_info_field(&record.info, "END");
    Ok(SymbolicDeletionMaterialization::LeadingAnchor)
}

pub(super) fn is_symbolic_allele(alt: &str) -> bool {
    alt == "*" || (alt.starts_with('<') && alt.ends_with('>'))
}

/// Project a called multi-allelic record into biallelic records. Symbolic
/// deletions need this even without primitive decomposition; phased calls need
/// independent normalization streams because legacy never re-aggregates them.
pub(super) struct MaterializedAlleleRecord {
    pub(super) record: RawVcfRecord,
    pub(super) reverse_hetalt_samples: Vec<bool>,
}

pub(super) fn split_called_alleles(record: &RawVcfRecord) -> Vec<MaterializedAlleleRecord> {
    let alts = record.alt_allele.split(',').collect::<Vec<_>>();
    if alts.len() <= 1 {
        return vec![MaterializedAlleleRecord {
            record: record.clone(),
            reverse_hetalt_samples: Vec::new(),
        }];
    }
    let format_keys = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect::<Vec<_>>())
        .unwrap_or_default();
    let gt_index = format_keys.iter().position(|field| *field == "GT");
    let ad_index = format_keys.iter().position(|field| *field == "AD");

    alts.iter()
        .enumerate()
        .filter_map(|(alt_index, alt)| {
            let target = alt_index + 1;
            let called = gt_index.is_some_and(|gt_index| {
                record.samples.iter().any(|sample| {
                    sample
                        .split(':')
                        .nth(gt_index)
                        .is_some_and(|gt| genotype_calls_allele(gt, target))
                })
            });
            if !called {
                return None;
            }

            let mut split = record.clone();
            split.alt_allele = (*alt).to_string();
            let mut reverse_hetalt_samples = Vec::with_capacity(split.samples.len());
            for sample in &mut split.samples {
                let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
                if let Some(gt_index) = gt_index
                    && let Some(gt) = cells.get_mut(gt_index)
                {
                    reverse_hetalt_samples.push(genotype_is_hetalt_without_ref(gt, target));
                    *gt = project_split_genotype(gt, target);
                } else {
                    reverse_hetalt_samples.push(false);
                }
                if let Some(ad_index) = ad_index
                    && let Some(ad) = cells.get_mut(ad_index)
                {
                    *ad = project_split_ad(ad, target);
                }
                *sample = cells.join(":");
            }
            Some(MaterializedAlleleRecord {
                record: split,
                reverse_hetalt_samples,
            })
        })
        .collect()
}

pub(super) fn genotype_calls_allele(gt: &str, target: usize) -> bool {
    gt.split(['/', '|'])
        .any(|allele| allele.parse::<usize>().ok() == Some(target))
}

pub(super) fn genotype_is_hetalt_without_ref(gt: &str, target: usize) -> bool {
    let alleles = gt
        .split(['/', '|'])
        .filter_map(|allele| allele.parse::<usize>().ok())
        .collect::<Vec<_>>();
    gt.contains('/')
        && target > 1
        && alleles.contains(&target)
        && alleles.iter().all(|allele| *allele > 0)
        && alleles.iter().any(|allele| *allele != target)
}

pub(super) fn restore_legacy_hetalt_orientation(
    record: &mut RawVcfRecord,
    reverse_samples: &[bool],
) {
    let Some(gt_index) = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|field| field == "GT"))
    else {
        return;
    };
    for (sample, reverse) in record.samples.iter_mut().zip(reverse_samples) {
        if !reverse {
            continue;
        }
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        if let Some(gt) = cells.get_mut(gt_index) {
            match gt.as_str() {
                "0/1" => *gt = "1/0".to_string(),
                "0|1" => *gt = "1|0".to_string(),
                _ => {}
            }
        }
        *sample = cells.join(":");
    }
}

pub(super) fn remove_info_field(info: &str, key: &str) -> String {
    let retained = info
        .split(';')
        .filter(|entry| {
            !entry.is_empty()
                && *entry != "."
                && entry
                    .split_once('=')
                    .map_or(*entry != key, |(name, _)| name != key)
        })
        .collect::<Vec<_>>();
    if retained.is_empty() {
        ".".to_string()
    } else {
        retained.join(";")
    }
}

pub(super) fn resolve_somatic_mode(args: &PreprocessArgs) -> Option<SomaticGtMode> {
    args.set_gt.or(args.somatic.then_some(SomaticGtMode::Half))
}

/// Mirror `remove_nonref_gt_variants.py`: only the final ALT is treated as
/// `<NON_REF>`, and a record is dropped when any sample's first cell calls its
/// allele index. The legacy script intentionally assumes GT is the first FORMAT
/// cell, so this helper does too.
pub(crate) fn calls_non_ref_allele(record: &RawVcfRecord) -> bool {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    if alts.last() != Some(&"<NON_REF>") {
        return false;
    }
    let non_ref_index = alts.len();
    record.samples.iter().any(|sample| {
        sample
            .split(':')
            .next()
            .unwrap_or_default()
            .split(['/', '|'])
            .any(|token| token.parse::<usize>().ok() == Some(non_ref_index))
    })
}

pub(super) fn trim_uncalled_non_ref(record: &mut RawVcfRecord) -> bool {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    if alts.last() != Some(&"<NON_REF>") {
        return record.alt_allele != ".";
    }
    let concrete = &alts[..alts.len().saturating_sub(1)];
    if concrete.is_empty() {
        return false;
    }
    let mapping: Vec<usize> = (0..=alts.len())
        .map(|index| if index < alts.len() { index } else { 0 })
        .collect();
    let format_keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect())
        .unwrap_or_default();
    let gt_index = format_keys.iter().position(|key| *key == "GT");
    let ad_index = format_keys.iter().position(|key| *key == "AD");
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        if let Some(index) = gt_index
            && let Some(gt) = cells.get_mut(index)
        {
            *gt = remap_gt(gt, &mapping);
        }
        if let Some(index) = ad_index
            && let Some(ad) = cells.get_mut(index)
        {
            let mut depths: Vec<&str> = ad.split(',').collect();
            if depths.len() == alts.len() + 1 {
                depths.pop();
                *ad = depths.join(",");
            }
        }
        *sample = cells.join(":");
    }
    record.alt_allele = concrete.join(",");
    true
}

/// Retain the union of non-reference alleles called by all samples, remapping
/// allele-indexed GT and AD cells. Returns false when the record has no called
/// ALT and is therefore not a variant after `VariantCallsOnly`.
pub(super) fn retain_called_alternates(record: &mut RawVcfRecord) -> bool {
    let alts = record.alt_allele.split(',').collect::<Vec<_>>();
    if alts.is_empty() || alts == ["."] {
        return false;
    }
    let format_keys = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect::<Vec<_>>())
        .unwrap_or_default();
    let Some(gt_index) = format_keys.iter().position(|field| *field == "GT") else {
        return false;
    };
    let ad_index = format_keys.iter().position(|field| *field == "AD");
    let mut called = BTreeSet::new();
    for sample in &record.samples {
        if let Some(gt) = sample.split(':').nth(gt_index) {
            called.extend(
                gt.split(['/', '|'])
                    .filter_map(|allele| allele.parse::<usize>().ok())
                    // `VariantCallsOnly` does not treat a spanning-deletion
                    // (`*`) as an emitted variant. Its call is remapped to
                    // reference when a concrete ALT remains, or the whole
                    // record is discarded when it is the only ALT.
                    .filter(|allele| {
                        *allele > 0 && *allele <= alts.len() && alts[*allele - 1] != "*"
                    }),
            );
        }
    }
    if called.is_empty() {
        return false;
    }

    let mut mapping = vec![0usize; alts.len() + 1];
    let mut retained_alts = Vec::with_capacity(called.len());
    let mut retained_old_indices = Vec::with_capacity(called.len());
    for (alt_index, alt) in alts.iter().enumerate() {
        let old_index = alt_index + 1;
        if called.contains(&old_index) {
            mapping[old_index] = retained_alts.len() + 1;
            retained_alts.push(*alt);
            retained_old_indices.push(old_index);
        }
    }

    for sample in &mut record.samples {
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        if let Some(gt) = cells.get_mut(gt_index) {
            *gt = remap_gt(gt, &mapping);
            // VariantAlleleRemover canonicalizes an unphased REF/ALT call
            // after projecting removed alleles. The regular writer path only
            // sees this orientation for phased genotypes, so do it here.
            if *gt == "1/0" {
                *gt = "0/1".to_string();
            }
        }
        if let Some(ad_index) = ad_index
            && let Some(ad) = cells.get_mut(ad_index)
        {
            let depths = ad.split(',').collect::<Vec<_>>();
            if depths.len() == alts.len() + 1 {
                let mut projected = Vec::with_capacity(retained_old_indices.len() + 1);
                projected.push(depths[0]);
                projected.extend(
                    retained_old_indices
                        .iter()
                        .map(|old_index| depths[*old_index]),
                );
                *ad = projected.join(",");
            }
        }
        *sample = cells.join(":");
    }
    record.alt_allele = retained_alts.join(",");
    true
}

pub(super) fn mask_genotypes_wider_than_diploid(record: &mut RawVcfRecord) {
    let Some(gt_index) = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|field| field == "GT"))
    else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        let Some(gt) = cells.get_mut(gt_index) else {
            continue;
        };
        if gt.split(['/', '|']).count() > 2 {
            *gt = ".".to_string();
        }
        *sample = cells.join(":");
    }
}

/// Port the standalone gVCF conversion pipeline:
/// `N_ALT >= 2` → keep only GT/DP/GQ → trim uncalled alleles → exclude
/// `<NON_REF>`. Returns false when the record is no longer a variant.
pub(super) fn convert_gvcf_record(record: &mut RawVcfRecord) -> bool {
    let alts: Vec<String> = record.alt_allele.split(',').map(str::to_string).collect();
    if alts.len() < 2 {
        return false;
    }

    let format_keys: Vec<String> = record
        .format
        .as_deref()
        .map(|format| format.split(':').map(str::to_string).collect())
        .unwrap_or_default();
    let gt_index = format_keys.iter().position(|key| key == "GT");
    let mut called = BTreeSet::new();
    if let Some(gt_index) = gt_index {
        for sample in &record.samples {
            if let Some(gt) = sample.split(':').nth(gt_index) {
                called.extend(
                    gt.split(['/', '|'])
                        .filter_map(|token| token.parse::<usize>().ok())
                        .filter(|index| *index > 0),
                );
            }
        }
    }

    // If `<NON_REF>` remains called after allele trimming, the final legacy
    // bcftools exclusion removes the whole record.
    if alts
        .iter()
        .enumerate()
        .any(|(index, alt)| alt == "<NON_REF>" && called.contains(&(index + 1)))
    {
        return false;
    }

    let mut kept = Vec::new();
    let mut mapping = vec![0usize; alts.len() + 1];
    for (index, alt) in alts.iter().enumerate() {
        let old_index = index + 1;
        if alt != "<NON_REF>" && called.contains(&old_index) {
            mapping[old_index] = kept.len() + 1;
            kept.push(alt.clone());
        }
    }
    let retained_indices: Vec<usize> = format_keys
        .iter()
        .enumerate()
        .filter_map(|(index, key)| matches!(key.as_str(), "GT" | "DP" | "GQ").then_some(index))
        .collect();
    let retained_keys: Vec<String> = retained_indices
        .iter()
        .map(|index| format_keys[*index].clone())
        .collect();
    for sample in &mut record.samples {
        let cells: Vec<&str> = sample.split(':').collect();
        let mut retained: Vec<String> = retained_indices
            .iter()
            .map(|index| cells.get(*index).copied().unwrap_or(".").to_string())
            .collect();
        if let Some(new_gt_index) = retained_keys.iter().position(|key| key == "GT") {
            retained[new_gt_index] = remap_gt(&retained[new_gt_index], &mapping);
        }
        *sample = retained.join(":");
    }

    record.alt_allele = if kept.is_empty() {
        ".".to_string()
    } else {
        kept.join(",")
    };
    record.info = ".".to_string();
    record.format = (!retained_keys.is_empty()).then(|| retained_keys.join(":"));
    true
}

pub(super) fn filter_gvcf_headers(headers: &mut Vec<String>) {
    headers.retain(|line| {
        if line.starts_with("##INFO=<ID=") {
            return false;
        }
        if let Some(body) = line.strip_prefix("##FORMAT=<ID=") {
            let id = body.split([',', '>']).next().unwrap_or_default();
            return matches!(id, "GT" | "DP" | "GQ");
        }
        true
    });
}

pub(super) fn ensure_missing_ad(record: &mut RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    let keys: Vec<&str> = format.split(':').collect();
    if keys.contains(&"AD") {
        return;
    }
    let alt_count = record
        .alt_allele
        .split(',')
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .count();
    let missing_ad = vec!["."; alt_count + 1].join(",");
    record.format = Some(format!("{format}:AD"));
    for sample in &mut record.samples {
        sample.push(':');
        sample.push_str(&missing_ad);
    }
}

pub(super) fn ensure_missing_dp(record: &mut RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    if format.split(':').any(|field| field == "DP") {
        return;
    }
    record.format = Some(format!("{format}:DP"));
    for sample in &mut record.samples {
        sample.push_str(":0");
    }
}

/// VariantWriter emits unphased, canonical diploid genotypes. Ref/ALT calls
/// use ascending order (`1|0` → `0/1`), while het-of-ALT calls use the legacy
/// aggregator's later/earlier order (`1|2` → `2/1`).
pub(super) fn canonicalize_legacy_genotypes(record: &mut RawVcfRecord) {
    let Some(gt_index) = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|field| field == "GT"))
    else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        let Some(gt) = cells.get_mut(gt_index) else {
            continue;
        };
        if !gt.contains('|') {
            continue;
        }
        let alleles = gt.split(['/', '|']).collect::<Vec<_>>();
        if alleles.len() == 2
            && let (Ok(left), Ok(right)) =
                (alleles[0].parse::<usize>(), alleles[1].parse::<usize>())
        {
            *gt = if left > 0 && right > 0 && left != right {
                format!("{}/{}", left.max(right), left.min(right))
            } else {
                format!("{}/{}", left.min(right), left.max(right))
            };
        } else {
            *gt = gt.replace('|', "/");
        }
        *sample = cells.join(":");
    }
}

pub(super) fn somatic_info_sample_names(headers: &[String]) -> Vec<String> {
    let input_names: Vec<String> = headers
        .iter()
        .find(|line| line.starts_with("#CHROM"))
        .map(|line| line.split('\t').skip(9).map(str::to_string).collect())
        .unwrap_or_default();
    if input_names.len() <= 1 {
        vec!["SAMPLE".to_string()]
    } else {
        input_names
    }
}

/// Legacy `alleles` declares one INFO field per input FORMAT/sample pair.
/// A single input sample is deliberately renamed to the output sample name
/// `SAMPLE`; multi-sample inputs keep their original names as INFO prefixes.
pub(super) fn append_somatic_info_headers(headers: &mut Vec<String>, sample_names: &[String]) {
    let format_headers: Vec<String> = headers
        .iter()
        .filter(|line| line.starts_with("##FORMAT=<ID="))
        .cloned()
        .collect();
    let mut additions = Vec::new();
    for line in format_headers {
        let Some(body) = line.strip_prefix("##FORMAT=<ID=") else {
            continue;
        };
        let id_end = body.find([',', '>']).unwrap_or(body.len());
        let id = &body[..id_end];
        let suffix = &body[id_end..];
        for sample in sample_names {
            additions.push(format!("##INFO=<ID={sample}_{id}{suffix}"));
        }
    }
    let insertion = headers
        .iter()
        .position(|line| line.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    headers.splice(insertion..insertion, additions);
}

pub(super) fn append_info_value(info: &mut String, key: &str, value: &str) {
    if value.is_empty() || value == "." {
        return;
    }
    if info.is_empty() || info == "." {
        *info = format!("{key}={value}");
    } else {
        info.push(';');
        info.push_str(key);
        info.push('=');
        info.push_str(value);
    }
}

pub(super) fn copy_sample_formats_to_info(record: &mut RawVcfRecord, sample_names: &[String]) {
    let keys: Vec<String> = record
        .format
        .as_deref()
        .map(|format| format.split(':').map(str::to_string).collect())
        .unwrap_or_default();
    let samples = record.samples.clone();
    // `alleles` declares fields for every input sample but, due to its use of
    // the translated output header while mutating the record, only the first
    // sample's values survive into the emitted INFO payload. Preserve this
    // long-standing observable quirk for byte parity.
    for (sample_index, sample) in samples.iter().take(1).enumerate() {
        let prefix = sample_names
            .get(sample_index)
            .or_else(|| sample_names.first())
            .map(String::as_str)
            .unwrap_or("SAMPLE");
        for (key, value) in keys.iter().zip(sample.split(':')) {
            let encoded;
            let value = if key == "GT" {
                encoded = bcf_encoded_gt(value);
                encoded.as_str()
            } else {
                value
            };
            append_info_value(&mut record.info, &format!("{prefix}_{key}"), value);
        }
    }
}

pub(super) fn first_sample_gt(record: &RawVcfRecord) -> String {
    let gt_index = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|key| key == "GT"));
    gt_index
        .and_then(|index| record.samples.first()?.split(':').nth(index))
        .unwrap_or(".")
        .to_string()
}

pub(super) fn converted_somatic_gt(mode: SomaticGtMode, has_alt: bool) -> String {
    let allele = if has_alt { "1" } else { "0" };
    match mode {
        SomaticGtMode::Half => format!("./{allele}"),
        SomaticGtMode::Hemi => allele.to_string(),
        SomaticGtMode::Het => format!("0/{allele}"),
        SomaticGtMode::Hom => format!("{allele}/{allele}"),
        SomaticGtMode::First => unreachable!("first mode preserves the source GT"),
    }
}

/// Port the record-level behavior of legacy's `alleles` helper: move every
/// input FORMAT into sample-prefixed INFO, collapse to one output sample, and
/// split each ALT into its own record unless `first` was selected.
pub(super) fn convert_somatic_record(
    record: &RawVcfRecord,
    mode: SomaticGtMode,
    sample_names: &[String],
) -> Vec<RawVcfRecord> {
    let first_gt = first_sample_gt(record);
    let mut base = record.clone();
    copy_sample_formats_to_info(&mut base, sample_names);
    base.format = Some("GT:AD:DP".to_string());

    let with_missing_depths = |gt: String, alt_count: usize| {
        let ad = vec!["."; alt_count.saturating_add(1)].join(",");
        format!("{gt}:{ad}:0")
    };

    if mode == SomaticGtMode::First {
        let alt_count = base
            .alt_allele
            .split(',')
            .filter(|alt| !alt.is_empty() && *alt != ".")
            .count();
        base.samples = vec![with_missing_depths(first_gt, alt_count)];
        return vec![base];
    }

    let alts: Vec<String> = base
        .alt_allele
        .split(',')
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .map(str::to_string)
        .collect();
    if alts.is_empty() {
        base.samples = vec![with_missing_depths(converted_somatic_gt(mode, false), 0)];
        return vec![base];
    }

    alts.into_iter()
        .map(|alt| {
            let mut split = base.clone();
            split.alt_allele = alt;
            split.samples = vec![with_missing_depths(converted_somatic_gt(mode, true), 1)];
            split
        })
        .collect()
}

/// Reproduce the subsequent partial-credit aggregation that consumes the
/// biallelic records emitted by `alleles`. Non-hom modes merge sibling ALTs
/// back into one location with the legacy reversed het-alt GT; hom-alt calls
/// remain separate records.
pub(super) fn finalize_somatic_records(
    mut records: Vec<RawVcfRecord>,
    mode: SomaticGtMode,
) -> Vec<RawVcfRecord> {
    if records.is_empty() {
        return records;
    }
    if mode == SomaticGtMode::First {
        records.retain_mut(trim_uncalled_non_ref);
        return records;
    }
    records.retain(|record| {
        !record.alt_allele.starts_with('<') && record.alt_allele != "*" && record.alt_allele != "."
    });
    if records.is_empty() {
        return records;
    }
    if mode == SomaticGtMode::Hom {
        records.sort_by(|left, right| {
            left.alt_allele
                .len()
                .cmp(&right.alt_allele.len())
                .then(left.alt_allele.cmp(&right.alt_allele))
        });
        return records;
    }
    if records.len() == 1 {
        records[0].samples = vec!["0/1:.,.:0".to_string()];
        return records;
    }

    let mut merged = records.remove(0);
    let mut alts = vec![merged.alt_allele.clone()];
    alts.extend(records.into_iter().map(|record| record.alt_allele));
    alts.sort_by(|left, right| left.len().cmp(&right.len()).then(left.cmp(right)));
    alts.dedup();
    merged.alt_allele = alts.join(",");
    let ad = vec!["."; alts.len() + 1].join(",");
    merged.samples = vec![format!("2/1:{ad}:0")];
    vec![merged]
}

pub(super) fn finalize_somatic_for_pipeline(
    mut converted: Vec<RawVcfRecord>,
    mode: SomaticGtMode,
    normalization_enabled: bool,
) -> Vec<RawVcfRecord> {
    if normalization_enabled {
        finalize_somatic_records(converted, mode)
    } else {
        // scmp-somatic disables the partial-credit normalizer. In that lane
        // the `alleles` helper's half-call (`./1`) and per-ALT record grain
        // pass straight through to scmp; sibling aggregation belongs to the
        // normalizer and must not silently turn the call back into `0/1`.
        for record in &mut converted {
            let gt = first_sample_gt(record);
            record.format = Some("GT".to_string());
            record.samples = vec![gt];
        }
        converted
    }
}

pub(super) fn rewrite_header_for_single_sample(headers: &mut [String]) {
    if let Some(line) = headers.iter_mut().find(|line| line.starts_with("#CHROM")) {
        *line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE".to_string();
    }
}
