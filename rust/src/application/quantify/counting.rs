//! Cohesive quantify counting responsibility.

use super::{ClassifiedVariant, INDEL_SUBTYPES, QuantifyCountMaps, QuantifyTypeCounts};
use crate::domain::{CountsBucket, RawVcfRecord, TypeCounts};
use std::collections::{BTreeMap, BTreeSet};

pub(super) fn classify_side(
    record: &RawVcfRecord,
    sample_index: usize,
) -> Option<ClassifiedVariant> {
    let fields = record.sample_map(sample_index);
    let bd = fields.get("BD")?.to_string();
    if bd == "." || bd == "N" {
        return None;
    }
    // XCMP/GA4GH quantification classifies the active genotype alleles, not
    // the record-wide REF/ALT tuple.  This matters for mixed multi-allelic
    // records such as REF=T ALT=C,TATC where the selected allele can be a SNP
    // even though another unselected allele is an insertion.  The quantifier
    // stores that genotype-aware result in BVT/BI/BLT before accumulating.
    let variant_type = fields.get("BVT")?.to_string();
    if !matches!(variant_type.as_str(), "SNP" | "INDEL") {
        return None;
    }
    let info_tokens = fields
        .get("BI")
        .map(String::as_str)
        .unwrap_or(".")
        .split(',')
        .collect::<BTreeSet<_>>();
    let subtypes = if variant_type == "INDEL" {
        info_tokens
            .iter()
            .filter_map(|token| {
                let subtype = token.to_ascii_uppercase();
                INDEL_SUBTYPES
                    .contains(&subtype.as_str())
                    .then_some(subtype)
            })
            .collect()
    } else {
        Vec::new()
    };
    let location_type = fields.get("BLT").map(String::as_str).unwrap_or(".");
    let subsets = parse_subsets(&record.info);
    let fp_class = fp_class(&bd, fields.get("BK").map(String::as_str));
    Some(ClassifiedVariant {
        variant_type,
        subtypes,
        ti: usize::from(info_tokens.contains("ti")),
        tv: usize::from(info_tokens.contains("tv")),
        het: location_type == "het",
        homalt: location_type == "homalt",
        status: bd,
        passes_filter: record.filter == "PASS" || record.filter == ".",
        subsets,
        fp_class,
    })
}

pub(super) fn fp_class(decision: &str, match_kind: Option<&str>) -> Option<&'static str> {
    if decision != "FP" {
        return None;
    }
    match match_kind {
        Some("am") => Some("gt"),
        Some("lm") => Some("al"),
        _ => None,
    }
}

#[cfg(test)]
pub(super) fn query_fp_class(record: &RawVcfRecord) -> Option<&'static str> {
    query_fp_class_for_sample(record, 1)
}

pub(super) fn query_fp_class_for_sample(
    record: &RawVcfRecord,
    sample_index: usize,
) -> Option<&'static str> {
    let fields = record.sample_map(sample_index);
    fp_class(
        fields.get("BD").map(String::as_str).unwrap_or("."),
        fields.get("BK").map(String::as_str),
    )
}

pub(super) fn parse_subsets(info: &str) -> Vec<String> {
    info.split(';')
        .find_map(|entry| entry.strip_prefix("Regions="))
        .map(|entry| {
            entry
                .split(',')
                .filter(|subset| !subset.is_empty() && *subset != "CONF")
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default()
}

pub(super) fn register_subsets(
    subsets_present: &mut BTreeMap<String, BTreeSet<String>>,
    classified: &ClassifiedVariant,
) {
    for subset in &classified.subsets {
        subsets_present
            .entry(classified.variant_type.clone())
            .or_default()
            .insert(subset.clone());
    }
}

fn walk_count_maps(
    counts: &mut QuantifyCountMaps,
    classified: &ClassifiedVariant,
    mut apply: impl FnMut(&mut QuantifyTypeCounts),
) {
    apply(
        counts
            .by_type
            .entry(classified.variant_type.clone())
            .or_default(),
    );
    for subtype in &classified.subtypes {
        apply(
            counts
                .by_subtype
                .entry(classified.variant_type.clone())
                .or_default()
                .entry(subtype.clone())
                .or_default(),
        );
    }
    for subset in &classified.subsets {
        apply(
            counts
                .by_subset_type
                .entry(subset.clone())
                .or_default()
                .entry(classified.variant_type.clone())
                .or_default(),
        );
        for subtype in &classified.subtypes {
            apply(
                counts
                    .by_subset_subtype
                    .entry(subset.clone())
                    .or_default()
                    .entry(classified.variant_type.clone())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default(),
            );
        }
    }
}

pub(super) fn record_truth(counts: &mut QuantifyCountMaps, classified: &ClassifiedVariant) {
    if !matches!(classified.status.as_str(), "TP" | "FN") {
        return;
    }
    walk_count_maps(counts, classified, |stats| {
        add_variant_stats(&mut stats.truth_total, classified);
    });
    walk_count_maps(counts, classified, |stats| {
        add_variant_stats(truth_bucket(stats, &classified.status), classified);
    });
}

pub(super) fn record_truth_total_only(
    counts: &mut QuantifyCountMaps,
    classified: &ClassifiedVariant,
) {
    if !matches!(classified.status.as_str(), "TP" | "FN") {
        return;
    }
    walk_count_maps(counts, classified, |stats| {
        add_variant_stats(&mut stats.truth_total, classified);
    });
}

pub(super) fn record_truth_filtered(
    counts: &mut QuantifyCountMaps,
    classified: &ClassifiedVariant,
) {
    record_truth_total_only(counts, classified);
    if !matches!(classified.status.as_str(), "TP" | "FN") {
        return;
    }
    walk_count_maps(counts, classified, |stats| {
        add_variant_stats(truth_bucket(stats, &classified.status), classified);
    });
}

pub(super) fn record_query(counts: &mut QuantifyCountMaps, classified: &ClassifiedVariant) {
    if !matches!(classified.status.as_str(), "TP" | "FP" | "UNK" | "AMBI") {
        return;
    }
    walk_count_maps(counts, classified, |stats| {
        record_query_stats(stats, classified);
    });
}

pub(super) fn record_query_stats(stats: &mut QuantifyTypeCounts, classified: &ClassifiedVariant) {
    add_variant_stats(&mut stats.query_total, classified);
    add_variant_stats(query_bucket(stats, &classified.status), classified);
    match classified.fp_class {
        Some("gt") => stats.fp_gt += 1,
        Some("al") => stats.fp_al += 1,
        _ => {}
    }
}

pub(super) fn truth_bucket<'a>(stats: &'a mut TypeCounts, status: &str) -> &'a mut CountsBucket {
    match status {
        "TP" => &mut stats.truth_tp,
        "FN" => &mut stats.truth_fn,
        _ => &mut stats.truth_total,
    }
}

pub(super) fn query_bucket<'a>(stats: &'a mut TypeCounts, status: &str) -> &'a mut CountsBucket {
    match status {
        "TP" => &mut stats.query_tp,
        "FP" => &mut stats.query_fp,
        "UNK" | "AMBI" => &mut stats.query_unk,
        _ => &mut stats.query_total,
    }
}

pub(super) fn add_variant_stats(bucket: &mut CountsBucket, classified: &ClassifiedVariant) {
    bucket.total += 1;
    bucket.ti += classified.ti;
    bucket.tv += classified.tv;
    if classified.het {
        bucket.het += 1;
    }
    if classified.homalt {
        bucket.homalt += 1;
    }
}

pub(super) fn derive_pass_truth_false_negatives(counts: &mut QuantifyCountMaps) {
    let derive = |stats: &mut TypeCounts| {
        stats.truth_fn = subtract_bucket(&stats.truth_total, &stats.truth_tp);
    };
    for stats in counts.by_type.values_mut() {
        derive(stats);
    }
    for by_subtype in counts.by_subtype.values_mut() {
        for stats in by_subtype.values_mut() {
            derive(stats);
        }
    }
    for by_type in counts.by_subset_type.values_mut() {
        for stats in by_type.values_mut() {
            derive(stats);
        }
    }
    for by_type in counts.by_subset_subtype.values_mut() {
        for by_subtype in by_type.values_mut() {
            for stats in by_subtype.values_mut() {
                derive(stats);
            }
        }
    }
}

pub(super) fn subtract_bucket(total: &CountsBucket, matched: &CountsBucket) -> CountsBucket {
    CountsBucket {
        total: total.total.saturating_sub(matched.total),
        ti: total.ti.saturating_sub(matched.ti),
        tv: total.tv.saturating_sub(matched.tv),
        het: total.het.saturating_sub(matched.het),
        homalt: total.homalt.saturating_sub(matched.homalt),
    }
}
