//! Extracted cohesive responsibility from the command façade.

use super::AnnotatedRow;
use crate::adapters::vcf::{self, Variant};
use crate::domain::{CountsBucket, RawVcfRecord, TypeCounts};
use std::collections::{BTreeMap, BTreeSet};

struct ComparisonSamples<'a> {
    format: Vec<&'a str>,
    truth: Vec<&'a str>,
    query: Vec<&'a str>,
}

fn comparison_samples(record: &RawVcfRecord) -> Option<ComparisonSamples<'_>> {
    let format = record.format.as_deref()?.split(':').collect();
    let truth = record.samples.first()?.split(':').collect();
    let query = record.samples.get(1)?.split(':').collect();
    Some(ComparisonSamples {
        format,
        truth,
        query,
    })
}

pub(super) fn collect_contigs(
    truth: &[Variant],
    query: &[Variant],
    locations: Option<&[vcf::LocationFilter]>,
) -> BTreeSet<String> {
    let mut contigs = BTreeSet::new();
    if let Some(filters) = locations {
        for filter in filters {
            match filter {
                vcf::LocationFilter::Contig(chrom) => {
                    contigs.insert(chrom.clone());
                }
                vcf::LocationFilter::Range { chrom, .. } => {
                    contigs.insert(chrom.clone());
                }
            }
        }
    }
    for variant in truth.iter().chain(query.iter()) {
        contigs.insert(variant.key.chrom.clone());
    }
    contigs
}

pub(super) fn report_subset_size(
    contig_non_n_lengths: &BTreeMap<String, usize>,
    contigs_in_play: &BTreeSet<String>,
    explicit_bcf: bool,
    implicit_bcf: bool,
) -> usize {
    if explicit_bcf {
        // The legacy BCF path initializes its aggregate region from the
        // complete FASTA dictionary, while the VCF path restricts it to the
        // contigs participating in the comparison. Preserve that observable
        // reporter quirk even though output encoding does not change calls.
        contig_non_n_lengths.values().sum()
    } else {
        let mut selected = contigs_in_play.clone();
        if implicit_bcf {
            // Paired BCF inputs select BCF intermediates and reports, but do
            // not set argparse's explicit `bcf` value. The pinned wrapper's
            // default chromosome discovery consequently retains both aliases
            // when the FASTA declares, for example, `1` and `chr1`.
            for contig in contigs_in_play {
                let alias = contig
                    .strip_prefix("chr")
                    .map(str::to_string)
                    .unwrap_or_else(|| format!("chr{contig}"));
                if contig_non_n_lengths.contains_key(&alias) {
                    selected.insert(alias);
                }
            }
        }
        selected
            .iter()
            .filter_map(|contig| contig_non_n_lengths.get(contig))
            .sum()
    }
}

/// Return the comma-delimited values for one exact INFO key.
///
/// `Regions` is normally the final field, but `--preserve-info` appends the
/// source annotations after it. Parsing the remainder of the INFO string as
/// region names consequently folds the next `;KEY=value` into the last tag.
pub(super) fn info_list_values<'a>(info: &'a str, key: &str) -> Vec<&'a str> {
    info.split(';')
        .find_map(|entry| entry.split_once('=').filter(|(name, _)| *name == key))
        .map(|(_, value)| value.split(',').filter(|value| !value.is_empty()).collect())
        .unwrap_or_default()
}

pub(super) fn derive_subset_counts(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, TypeCounts>> {
    let mut subsets: BTreeMap<String, BTreeMap<String, TypeCounts>> = BTreeMap::new();
    for row in rows {
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };
        let subset_tags = info_list_values(&row.record.info, "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
        if subset_tags.is_empty() {
            continue;
        }

        let truth_sample = SampleView::new(&samples.format, &samples.truth);
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let filtered_out = pass_only && !row.query_pass;

        for subset in subset_tags {
            let type_map = subsets.entry(subset.to_string()).or_default();
            if let Some(variant_type) = truth_sample.variant_type() {
                let stats = type_map.entry(variant_type.to_string()).or_default();
                truth_sample.add_truth(stats, filtered_out);
            }
            if !filtered_out && let Some(variant_type) = query_sample.variant_type() {
                let stats = type_map.entry(variant_type.to_string()).or_default();
                query_sample.add_query(stats);
            }
        }
    }
    subsets
}

/// Per-variant-type (FP.gt, FP.al) tally derived by walking the annotated
/// rows. Rows whose query failed the filter are excluded when `pass_only` is
/// true. Returns `{variant_type: (fp_gt_count, fp_al_count)}` for every
/// variant type observed.
pub(super) fn derive_fp_classes(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, (usize, usize)> {
    let mut out: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for row in rows {
        if pass_only && !row.query_pass {
            continue;
        }
        let Some(class) = row.fp_class else {
            continue;
        };
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let Some(variant_type) = query_sample.variant_type() else {
            continue;
        };
        let bucket = out.entry(variant_type.to_string()).or_default();
        if class == "gt" {
            bucket.0 += 1;
        } else if class == "al" {
            bucket.1 += 1;
        }
    }
    out
}

/// Per-subset variant of `derive_fp_classes`. Each row's `Regions=` tail
/// contributes to every named subset it carries (CONF is filtered out
/// the same way `derive_subset_counts` does so the keys align with the
/// subset rows in extended.csv).
pub(super) fn derive_subset_fp_classes(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, (usize, usize)>> {
    let mut out: BTreeMap<String, BTreeMap<String, (usize, usize)>> = BTreeMap::new();
    for row in rows {
        if pass_only && !row.query_pass {
            continue;
        }
        let Some(class) = row.fp_class else {
            continue;
        };
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };
        let subset_tags = info_list_values(&row.record.info, "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
        if subset_tags.is_empty() {
            continue;
        }
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let Some(variant_type) = query_sample.variant_type() else {
            continue;
        };
        for subset in subset_tags {
            let bucket = out
                .entry(subset.to_string())
                .or_default()
                .entry(variant_type.to_string())
                .or_default();
            if class == "gt" {
                bucket.0 += 1;
            } else if class == "al" {
                bucket.1 += 1;
            }
        }
    }
    out
}

/// Per-(variant_type, subtype) FP class tally for INDEL subtype rows in
/// extended.csv. Legacy emits FP.gt / FP.al at every (INDEL, subtype, *,
/// filter) row and at every (INDEL, subtype, TS_boundary|TS_contained,
/// filter) row. Each FP query row contributes to every indel-class token
/// in its multi-allelic BI (e.g. a hetalt FP with BI `i1_5,i6_15` adds
/// one to both I1_5 and I6_15 — same fanout rule the truth/query stats
/// use).
pub(super) type SubtypeFpClasses = BTreeMap<String, BTreeMap<String, (usize, usize)>>;
pub(super) type SubsetSubtypeFpClasses = BTreeMap<String, SubtypeFpClasses>;

pub(super) fn derive_subtype_fp_classes(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> SubtypeFpClasses {
    let mut out = SubtypeFpClasses::new();
    for row in rows {
        if pass_only && !row.query_pass {
            continue;
        }
        let Some(class) = row.fp_class else {
            continue;
        };
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let Some((variant_type, subtypes)) = query_sample.variant_type_and_subtypes() else {
            continue;
        };
        for subtype in subtypes {
            let bucket = out
                .entry(variant_type.to_string())
                .or_default()
                .entry(subtype)
                .or_default();
            if class == "gt" {
                bucket.0 += 1;
            } else if class == "al" {
                bucket.1 += 1;
            }
        }
    }
    out
}

/// Per-(subset, variant_type, subtype) FP class tally — the cross-product
/// counterpart of `derive_subtype_fp_classes`, plumbed into the
/// (INDEL, subtype, TS_*, filter) extended-csv rows.
pub(super) fn derive_subset_subtype_fp_classes(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> SubsetSubtypeFpClasses {
    let mut out = SubsetSubtypeFpClasses::new();
    for row in rows {
        if pass_only && !row.query_pass {
            continue;
        }
        let Some(class) = row.fp_class else {
            continue;
        };
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };
        let subset_tags = info_list_values(&row.record.info, "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
        if subset_tags.is_empty() {
            continue;
        }
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let Some((variant_type, subtypes)) = query_sample.variant_type_and_subtypes() else {
            continue;
        };
        for subset in subset_tags {
            for subtype in &subtypes {
                let bucket = out
                    .entry(subset.to_string())
                    .or_default()
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype.clone())
                    .or_default();
                if class == "gt" {
                    bucket.0 += 1;
                } else if class == "al" {
                    bucket.1 += 1;
                }
            }
        }
    }
    out
}

pub(super) fn derive_total_counts(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, TypeCounts> {
    let mut totals: BTreeMap<String, TypeCounts> = BTreeMap::new();
    for row in rows {
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };

        let truth_sample = SampleView::new(&samples.format, &samples.truth);
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let filtered_out = pass_only && !row.query_pass;

        if let Some(variant_type) = truth_sample.variant_type() {
            let stats = totals.entry(variant_type.to_string()).or_default();
            truth_sample.add_truth(stats, filtered_out);
        }
        if !filtered_out && let Some(variant_type) = query_sample.variant_type() {
            let stats = totals.entry(variant_type.to_string()).or_default();
            query_sample.add_query(stats);
        }
    }
    totals
}

/// Triple-nested counts for the extended CSV cross-product rows: subset
/// (TS_boundary / TS_contained) → variant_type (INDEL) → subtype (C1_5,
/// I1_5, etc.) → TypeCounts. Legacy emits 72 non-zero rows for these
/// stratifications per case; rust previously zero-filled them which
/// diverged extended.csv on every INDEL subtype × subset × filter cell.
pub(super) fn derive_subset_subtype_counts(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>> {
    let mut out: BTreeMap<String, BTreeMap<String, BTreeMap<String, TypeCounts>>> = BTreeMap::new();
    for row in rows {
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };
        let subset_tags = info_list_values(&row.record.info, "Regions")
            .into_iter()
            .filter(|tag| *tag != "CONF")
            .collect::<Vec<_>>();
        if subset_tags.is_empty() {
            continue;
        }

        let truth_sample = SampleView::new(&samples.format, &samples.truth);
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let filtered_out = pass_only && !row.query_pass;

        for subset in &subset_tags {
            let by_type = out.entry((*subset).to_string()).or_default();
            if let Some((variant_type, subtypes)) = truth_sample.variant_type_and_subtypes() {
                for subtype in subtypes {
                    let stats = by_type
                        .entry(variant_type.to_string())
                        .or_default()
                        .entry(subtype)
                        .or_default();
                    truth_sample.add_truth(stats, filtered_out);
                }
            }
            if !filtered_out
                && let Some((variant_type, subtypes)) = query_sample.variant_type_and_subtypes()
            {
                for subtype in subtypes {
                    let stats = by_type
                        .entry(variant_type.to_string())
                        .or_default()
                        .entry(subtype)
                        .or_default();
                    query_sample.add_query(stats);
                }
            }
        }
    }
    out
}

pub(super) fn derive_subtype_counts(
    rows: &[AnnotatedRow],
    pass_only: bool,
) -> BTreeMap<String, BTreeMap<String, TypeCounts>> {
    let mut subtypes: BTreeMap<String, BTreeMap<String, TypeCounts>> = BTreeMap::new();
    for row in rows {
        let Some(samples) = comparison_samples(&row.record) else {
            continue;
        };

        let truth_sample = SampleView::new(&samples.format, &samples.truth);
        let query_sample = SampleView::new(&samples.format, &samples.query);
        let filtered_out = pass_only && !row.query_pass;

        if let Some((variant_type, sub_list)) = truth_sample.variant_type_and_subtypes() {
            for subtype in sub_list {
                let stats = subtypes
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype)
                    .or_default();
                truth_sample.add_truth(stats, filtered_out);
            }
        }
        if !filtered_out
            && let Some((variant_type, sub_list)) = query_sample.variant_type_and_subtypes()
        {
            for subtype in sub_list {
                let stats = subtypes
                    .entry(variant_type.to_string())
                    .or_default()
                    .entry(subtype)
                    .or_default();
                query_sample.add_query(stats);
            }
        }
    }
    subtypes
}

pub(super) struct SampleView<'a> {
    gt: Option<&'a str>,
    bd: Option<&'a str>,
    bi: Option<&'a str>,
    bvt: Option<&'a str>,
}

impl<'a> SampleView<'a> {
    fn new(format_keys: &[&str], sample_parts: &'a [&str]) -> Self {
        let lookup = |name: &str| -> Option<&'a str> {
            format_keys
                .iter()
                .position(|key| *key == name)
                .and_then(|index| sample_parts.get(index).copied())
        };
        Self {
            gt: lookup("GT"),
            bd: lookup("BD"),
            bi: lookup("BI"),
            bvt: lookup("BVT"),
        }
    }

    fn variant_type(&self) -> Option<&'a str> {
        self.bvt.filter(|value| *value != "NOCALL")
    }

    /// Per-row subtypes for INDEL aggregation. Multi-allelic INDELs emit
    /// comma-joined BI strings (e.g. `i1_5,i6_15` for a 1|2 hetalt with one
    /// primitive in each size bucket). Legacy quantify fans these out per
    /// token — each indel-class primitive contributes to its own subtype
    /// bucket. The `ti`/`tv` tokens that decorate complex INDELs (BI like
    /// `c6_15,tv`) describe the SNP-side of a single complex primitive and
    /// must NOT spawn an extra bucket here, so we drop them.
    fn variant_type_and_subtypes(&self) -> Option<(&'a str, Vec<String>)> {
        let variant_type = self.variant_type()?;
        if variant_type != "INDEL" {
            return None;
        }
        let bi = self.bi?;
        let subtypes: Vec<String> = bi
            .split(',')
            .filter(|tok| !matches!(*tok, "ti" | "tv"))
            .map(|tok| tok.to_uppercase())
            .collect();
        if subtypes.is_empty() {
            return None;
        }
        Some((variant_type, subtypes))
    }

    fn add_truth(&self, stats: &mut TypeCounts, demote_tp_to_fn: bool) {
        // When deriving PASS-tier counts for a row whose query failed the
        // filter, legacy reclassifies the truth side from TP to FN: the good
        // match doesn't count because the query record wouldn't have been
        // considered in a pass-only run.
        let effective_bd = if demote_tp_to_fn && self.bd == Some("TP") {
            Some("FN")
        } else {
            self.bd
        };
        match effective_bd {
            Some("TP") | Some("FN") => add_sample_stats(&mut stats.truth_total, self),
            _ => {}
        }
        match effective_bd {
            Some("TP") => add_sample_stats(&mut stats.truth_tp, self),
            Some("FN") => add_sample_stats(&mut stats.truth_fn, self),
            _ => {}
        }
    }

    fn add_query(&self, stats: &mut TypeCounts) {
        match self.bd {
            Some("TP") | Some("FP") | Some("UNK") => add_sample_stats(&mut stats.query_total, self),
            _ => {}
        }
        match self.bd {
            Some("TP") => add_sample_stats(&mut stats.query_tp, self),
            Some("FP") => add_sample_stats(&mut stats.query_fp, self),
            Some("UNK") => add_sample_stats(&mut stats.query_unk, self),
            _ => {}
        }
    }
}

pub(super) fn add_sample_stats(bucket: &mut CountsBucket, sample: &SampleView<'_>) {
    bucket.total += 1;
    if let Some("SNP") = sample.bvt {
        // BI on multi-allelic hetalt SNPs is comma-separated (e.g.
        // `ti,tv` for GT=1|2 with one transition and one transversion
        // active alt). Legacy fans these out per primitive — one ti
        // and one tv contribution. Iterate the comma list and count
        // each tag once so single-allelic rows still increment by 1.
        if let Some(bi) = sample.bi {
            for tag in bi.split(',') {
                match tag {
                    "ti" => bucket.ti += 1,
                    "tv" => bucket.tv += 1,
                    _ => {}
                }
            }
        }
    }
    // Legacy's summary counts:
    // * het = exactly one allele is the reference index 0 (covers 0/1,
    //   1/0, 0|1, 1|0 AND 0|2, 2|0, 0|3, 3|0, …) — anything heterozygous
    //   with the reference base.
    // * homalt = both alleles equal AND non-zero (1/1, 1|1, 2|2, 3|3, …).
    // Hetalt (1|2, 2|1, …) lands in NEITHER bucket. Earlier the
    // classifier matched only literal `1/1`/`0/1` to mirror QUERY-side
    // counts (queries are split into per-primitive `0/1`/`1/1` rows),
    // but TRUTH-side rows preserve their original multi-allelic GT
    // through bcftools merge so the literal-only rule under-counted
    // every truth-only multi-allelic record.
    if let Some(gt) = sample.gt {
        let alleles: Vec<usize> = gt
            .split(['/', '|'])
            .map(|part| part.parse::<usize>().unwrap_or(0))
            .collect();
        if alleles.len() == 2 {
            let zero_count = alleles.iter().filter(|a| **a == 0).count();
            if zero_count == 1 {
                bucket.het += 1;
            } else if alleles[0] != 0 && alleles[0] == alleles[1] {
                bucket.homalt += 1;
            }
        }
    }
}

pub(super) fn add_variant_stats(bucket: &mut CountsBucket, variant: &Variant) {
    bucket.total += 1;
    if let Some(kind) = snp_bucket_label(variant) {
        if kind == "ti" {
            bucket.ti += 1;
        } else if kind == "tv" {
            bucket.tv += 1;
        }
    }
    if variant.is_het() {
        bucket.het += 1;
    }
    if variant.is_homalt() {
        bucket.homalt += 1;
    }
}

pub(super) fn add_variant_stats_subtype<F>(
    subtype_counts: &mut BTreeMap<String, BTreeMap<String, TypeCounts>>,
    variant_type: &str,
    variant: &Variant,
    bucket_selector: F,
) where
    F: Fn(&mut TypeCounts) -> &mut CountsBucket,
{
    if let Some(subtype) = subtype_label(variant) {
        let family = subtype_counts
            .entry(variant_type.to_string())
            .or_default()
            .entry(subtype)
            .or_default();
        add_variant_stats(bucket_selector(family), variant);
    }
}
