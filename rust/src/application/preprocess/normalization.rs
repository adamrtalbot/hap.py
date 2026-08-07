//! Reference-aware left-shift and trimming.

use super::LEFT_SHIFT_WINDOW;
use super::canonical::STALE_INFO_KEYS;
use super::genotype::expand_haploid_gt;
use crate::domain::RawVcfRecord;
use crate::engines::partial_credit;

#[allow(clippy::too_many_arguments)]
pub(super) fn apply_left_shift(record: &mut RawVcfRecord, reference: &[u8], neighbor_end: usize) {
    let ref_len = record.ref_allele.len();
    if ref_len == 0 {
        return;
    }
    let end = record.pos + ref_len - 1;
    let mut rv = partial_credit::RefVar {
        start: record.pos,
        end,
        alt: record.alt_allele.clone(),
    };
    // neighbor_end is the reference end of the previous record: don't allow
    // left-shifting into that span (mirrors legacy partialcredit.py behaviour).
    let pos_min = record
        .pos
        .saturating_sub(LEFT_SHIFT_WINDOW)
        .max(1)
        .max(neighbor_end);
    partial_credit::left_shift(reference, &mut rv, pos_min, true);

    // Rebuild REF from the reference slice now that start/end may have moved.
    let new_ref_len = rv.end as i64 - rv.start as i64 + 1;
    if new_ref_len < 0 || rv.start == 0 {
        return;
    }
    let new_end_usize = rv.start + (new_ref_len.max(0) as usize).saturating_sub(0);
    // Guard against slicing past the reference (noodles reads give us ASCII bytes).
    if rv.start == 0 || new_end_usize > reference.len() + 1 {
        return;
    }
    let ref_bytes = if new_ref_len <= 0 {
        &[][..]
    } else {
        &reference[rv.start - 1..rv.end]
    };
    // Uppercase soft-masked repeat bases so the rebuilt REF matches legacy's
    // canonical output.
    let new_ref = String::from_utf8_lossy(ref_bytes).to_ascii_uppercase();

    // Don't emit zero-length REF or ALT — VCF requires at least one base on
    // each side (partial-credit keeps a left-anchor via `ref_padding=true`, so
    // this is belt-and-braces against edge cases).
    if new_ref.is_empty() || rv.alt.is_empty() {
        return;
    }

    record.pos = rv.start;
    record.ref_allele = new_ref;
    record.alt_allele = rv.alt;
}

pub(super) fn record_reference_matches(record: &RawVcfRecord, reference: &[u8]) -> bool {
    let start = record.pos.saturating_sub(1);
    let end = start.saturating_add(record.ref_allele.len());
    reference
        .get(start..end)
        .is_some_and(|observed| observed.eq_ignore_ascii_case(record.ref_allele.as_bytes()))
}

/// Equivalent of `bcftools norm -f REF -c x -D` for one non-symbolic record.
/// Alleles are normalized independently, then padded to one common site so a
/// multi-allelic record retains its original allele indexes.
pub(super) fn normalize_bcftools_record(record: &mut RawVcfRecord, reference: &[u8]) {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    if alts.iter().any(|alt| {
        alt.is_empty()
            || *alt == "."
            || alt.starts_with('<')
            || *alt == "*"
            || alt.contains(['[', ']'])
    }) {
        return;
    }
    let end = record.end_pos();
    let mut normalized = Vec::with_capacity(alts.len());
    for alt in alts {
        let mut variant = partial_credit::RefVar {
            start: record.pos,
            end,
            alt: alt.to_string(),
        };
        partial_credit::left_shift(reference, &mut variant, 1, false);
        normalized.push(pad_bcftools_allele(variant, reference));
    }
    let common_start = normalized
        .iter()
        .map(|variant| variant.start)
        .min()
        .unwrap_or(record.pos);
    let common_end = normalized
        .iter()
        .map(|variant| variant.end)
        .max()
        .unwrap_or(end);
    if common_start == 0 || common_end < common_start || common_end > reference.len() {
        return;
    }
    let alts = normalized
        .into_iter()
        .map(|variant| {
            let mut allele = Vec::new();
            allele.extend_from_slice(&reference[common_start - 1..variant.start - 1]);
            allele.extend_from_slice(variant.alt.as_bytes());
            allele.extend_from_slice(&reference[variant.end..common_end]);
            String::from_utf8_lossy(&allele).to_ascii_uppercase()
        })
        .collect::<Vec<_>>();
    record.pos = common_start;
    record.ref_allele =
        String::from_utf8_lossy(&reference[common_start - 1..common_end]).to_ascii_uppercase();
    record.alt_allele = alts.join(",");
}

pub(super) fn uppercase_alleles_preserving_breakends(alts: &str) -> String {
    alts.split(',')
        .map(|alt| {
            let Some(first) = alt.find(['[', ']']) else {
                return alt.to_ascii_uppercase();
            };
            let Some(relative_second) = alt[first + 1..].find(['[', ']']) else {
                return alt.to_ascii_uppercase();
            };
            let second = first + 1 + relative_second;
            format!(
                "{}{}{}",
                alt[..=first].to_ascii_uppercase(),
                &alt[first + 1..second],
                alt[second..].to_ascii_uppercase()
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

pub(super) fn materialize_unsupported_import_failure(record: &mut RawVcfRecord) -> bool {
    let breakend = record.alt_allele.contains(['[', ']']);
    let unsupported = record.alt_allele.split(',').any(|alt| {
        alt.contains(['[', ']']) || (alt.starts_with('<') && !matches!(alt, "<DEL>" | "<NON_REF>"))
    });
    if !unsupported {
        return false;
    }

    record.alt_allele = ".".to_string();
    let mut info = record
        .info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .filter(|field| {
            let key = field.split_once('=').map_or(*field, |(key, _)| key);
            (key != "END" || !breakend) && key != "IMPORT_FAIL"
        })
        .map(str::to_string)
        .collect::<Vec<_>>();
    // Breakends have no local span, whereas unsupported symbolic records keep
    // their declared END when the legacy reader materializes IMPORT_FAIL.
    if breakend || !info.iter().any(|field| field.starts_with("END=")) {
        info.push(format!("END={}", record.pos));
    }
    info.push("IMPORT_FAIL".to_string());
    record.info = info.join(";");

    let Some(gt_index) = record.format_keys().iter().position(|key| *key == "GT") else {
        return true;
    };
    let ad_index = record.format_keys().iter().position(|key| *key == "AD");
    for sample in &mut record.samples {
        let mut fields = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        if let Some(gt) = fields.get_mut(gt_index) {
            *gt = "0/0".to_string();
        }
        // The failed BND import has no alternate allele. VariantWriter
        // consequently writes just the reference AD value; retaining the old
        // alternate depth would also make the generated ADO non-zero.
        if let Some(ad_index) = ad_index
            && let Some(ad) = fields.get_mut(ad_index)
        {
            *ad = ad.split(',').next().unwrap_or(".").to_string();
        }
        *sample = fields.join(":");
    }
    true
}

pub(super) fn pad_bcftools_allele(
    mut variant: partial_credit::RefVar,
    reference: &[u8],
) -> partial_credit::RefVar {
    let reference_len = variant.end as i64 - variant.start as i64 + 1;
    if reference_len <= 0 && !variant.alt.is_empty() {
        if variant.start > 1 {
            let anchor = variant.start - 1;
            variant.start = anchor;
            variant.end = anchor;
            variant
                .alt
                .insert(0, reference[anchor - 1].to_ascii_uppercase() as char);
        }
    } else if reference_len > 0 && variant.alt.is_empty() {
        if variant.start > 1 {
            let anchor = variant.start - 1;
            variant.start = anchor;
            variant
                .alt
                .push(reference[anchor - 1].to_ascii_uppercase() as char);
        } else if variant.end < reference.len() {
            variant.end += 1;
            variant
                .alt
                .push(reference[variant.end - 1].to_ascii_uppercase() as char);
        }
    }
    variant
}

/// Remove stale allele-count INFO entries from `record.info` in place.
pub(super) fn strip_stale_info_keys(record: &mut RawVcfRecord) {
    if record.info.is_empty() || record.info == "." {
        return;
    }
    let kept: Vec<&str> = record
        .info
        .split(';')
        .filter(|entry| {
            let key = entry.split('=').next().unwrap_or(entry);
            !STALE_INFO_KEYS.contains(&key)
        })
        .collect();
    record.info = if kept.is_empty() {
        ".".to_string()
    } else {
        kept.join(";")
    };
}

/// Mirror legacy `VariantAlleleSplitter`'s haploid → diploid treatment:
///
/// * `0`   → `0/0`
/// * autosomal `<n>` (n > 0) → `0/<n>` (legacy creates a het half-call)
/// * sex-chromosome `<n>` (n > 0) → `<n>/<n>`
/// * `.`   → `./.`
/// * `./<n>` → `0/<n>` (legacy completes a diploid half-call with REF)
///
/// Fully called diploid inputs pass through unchanged.
pub(super) fn normalise_haploid_genotypes(record: &mut RawVcfRecord, male: bool) {
    if male {
        expand_male_sex_chromosome_genotypes(record);
    }
    let Some(format) = &record.format else {
        return;
    };
    let Some(gt_index) = format.split(':').position(|f| f == "GT") else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        if let Some(cell) = cells.get_mut(gt_index) {
            *cell = expand_haploid_gt(cell, false);
        }
        *sample = cells.join(":");
    }
}

pub(super) fn expand_male_sex_chromosome_genotypes(record: &mut RawVcfRecord) {
    if !matches!(record.chrom.as_str(), "X" | "Y" | "chrX" | "chrY") {
        return;
    }
    let Some(format) = &record.format else {
        return;
    };
    let Some(gt_index) = format.split(':').position(|field| field == "GT") else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        if let Some(gt) = cells.get_mut(gt_index) {
            *gt = expand_haploid_gt(gt, true);
        }
        *sample = cells.join(":");
    }
}
