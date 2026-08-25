//! Decorate comparison rows for output: INFO-field rewriting, ROC-field capture, and requantify handoff sanitizing.

#[cfg(test)]
use super::CLUSTER_GAP_BP;
#[cfg(test)]
use super::Cluster;
use super::genotype::parse_gt_alleles;
#[cfg(test)]
use super::matching::build_clusters_with_gap;
use super::matching::effective_refrange;
use super::metrics::info_list_values;
use super::{AnnotatedRow, Variant, VariantKey};
use crate::adapters::{
    metrics_json,
    report::suffixed_report_path,
    vcf::{self},
};
use crate::application::CompareArgs;
use crate::application::preprocess;
use crate::domain::{RawVcfRecord, allele_edit_bits, legacy_type_bits};
use anyhow::Result;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

pub(super) fn rewrite_compare_metrics(
    args: &CompareArgs,
    prefix: &Path,
    commandline: &str,
    roc_indices: &crate::engines::roc::MetricIndices,
) -> Result<()> {
    if args.no_json {
        return Ok(());
    }
    let mut tables = vec![(
        "summary.metrics",
        "summary.metrics",
        suffixed_report_path(prefix, "summary.csv"),
    )];
    if args.write_counts && !args.no_write_counts {
        tables.push((
            "all.metrics",
            "all.metrics",
            suffixed_report_path(prefix, "extended.csv"),
        ));
    }
    for id in &roc_indices.table_order {
        let path = suffixed_report_path(prefix, &format!("{id}.csv.gz"));
        if path.is_file() {
            tables.push((id.as_str(), id.as_str(), path));
        }
    }
    let refs = tables
        .iter()
        .map(|(id, label, path)| (*id, *label, path.as_path()))
        .collect::<Vec<_>>();
    metrics_json::write_metrics_gz_for_module_with_indices(
        &suffixed_report_path(prefix, "metrics.json.gz"),
        "hap.py.comparison",
        "hap.py",
        commandline,
        &refs,
        Some(&roc_indices.tables),
    )
}

pub(super) type InfoKey = (String, usize, String, String);
pub(super) type SemanticInfoKey = (String, usize, String, Vec<String>);

pub(super) fn semantic_info_key(record: &RawVcfRecord) -> SemanticInfoKey {
    let mut alts = record
        .alt_allele
        .split(',')
        .map(str::to_string)
        .collect::<Vec<_>>();
    alts.sort();
    let reference = if alts
        .iter()
        .any(|alt| alt.starts_with('<') && alt.ends_with('>'))
    {
        "*".to_string()
    } else {
        record.ref_allele.clone()
    };
    (record.chrom.clone(), record.pos, reference, alts)
}

#[cfg(test)]
pub(super) fn decorate_output_rows(
    rows: &mut [AnnotatedRow],
    truth: &[RawVcfRecord],
    query: &[RawVcfRecord],
    preserve_info: bool,
    output_vtc: bool,
    roc_field: &str,
) -> Result<()> {
    if !preserve_info && !output_vtc && matches!(roc_field, "QUAL" | "QQ") {
        return Ok(());
    }
    let mut decorations = DecorationIndex::default();
    for record in truth.iter().chain(query) {
        decorations.observe(record, preserve_info, roc_field);
    }
    decorate_output_rows_with_index(rows, &decorations, preserve_info, output_vtc, roc_field)
}

#[derive(Default)]
pub(super) struct DecorationIndex {
    preserved: BTreeMap<InfoKey, BTreeSet<String>>,
    semantic_preserved: BTreeMap<SemanticInfoKey, BTreeSet<String>>,
    roc_values: BTreeMap<InfoKey, String>,
}

impl DecorationIndex {
    pub(super) fn observe(&mut self, record: &RawVcfRecord, preserve_info: bool, roc_field: &str) {
        let key = (
            record.chrom.clone(),
            record.pos,
            record.ref_allele.clone(),
            record.alt_allele.clone(),
        );
        if preserve_info {
            for field in record
                .info
                .split(';')
                .filter(|field| !matches!(*field, "" | "."))
            {
                self.preserved
                    .entry(key.clone())
                    .or_default()
                    .insert(field.to_string());
                self.semantic_preserved
                    .entry(semantic_info_key(record))
                    .or_default()
                    .insert(field.to_string());
            }
        }
        if !matches!(roc_field, "QUAL" | "QQ" | ".") {
            let value = record
                .info
                .split(';')
                .find_map(|entry| {
                    entry
                        .split_once('=')
                        // Pinned xcmp looks up the --qq argument literally.
                        // Prefixes such as INFO. and FORMAT. are not parsed.
                        .filter(|(key, _)| *key == roc_field)
                        .map(|(_, value)| value.to_string())
                })
                .or_else(|| record.sample_map(0).get(roc_field).cloned());
            if let Some(value) = value {
                self.roc_values.insert(key, value);
            }
        }
    }
}

pub(super) fn decorate_output_rows_with_index(
    rows: &mut [AnnotatedRow],
    decorations: &DecorationIndex,
    preserve_info: bool,
    output_vtc: bool,
    roc_field: &str,
) -> Result<()> {
    for row in rows {
        let mut record = row.record.raw().clone();
        let key = (
            record.chrom.clone(),
            record.pos,
            record.ref_allele.clone(),
            record.alt_allele.clone(),
        );
        let current = info_fields_by_key(&record.info);
        let regions = current.get("Regions").cloned();
        let mut base = BTreeMap::<String, String>::new();
        if preserve_info {
            let fields = decorations.preserved.get(&key).into_iter().chain(
                decorations
                    .semantic_preserved
                    .get(&semantic_info_key(&record)),
            );
            for fields in fields {
                for field in fields {
                    let field_key = field.split_once('=').map_or(field.as_str(), |(key, _)| key);
                    if field_key != "Regions" {
                        base.insert(field_key.to_string(), field.clone());
                    }
                }
            }
        }
        if let Some(bs) = current.get("BS") {
            base.insert("BS".to_string(), bs.clone());
        }
        let truth_fields = record.sample_map(0);
        let query_fields = record.sample_map(1);
        let comparison = legacy_comparison_fields(
            &record,
            &truth_fields,
            &query_fields,
            row.xcmp_ctype.unwrap_or("simple:match"),
        );
        if preserve_info {
            for (name, value) in comparison.preserved_fields(row.xcmp_hap_match) {
                base.insert(name.to_string(), value);
            }
            let iqq = if roc_field == "QUAL" {
                Some(comparison.iqq.as_str())
            } else {
                decorations.roc_values.get(&key).map(String::as_str)
            };
            if let Some(iqq) = iqq {
                base.insert("IQQ".to_string(), format!("IQQ={iqq}"));
            }
        }

        if !matches!(roc_field, "QUAL" | "QQ" | ".") && !decorations.roc_values.contains_key(&key) {
            // When xcmp's literal custom-field lookup misses, quantify reads
            // an absent IQQ: called query samples receive NaN, truth samples
            // remain missing, and no-call queries retain the zero sentinel.
            set_comparison_format_value(&mut record, 0, "QQ", ".");
            let query_called = record
                .sample_map(1)
                .get("BVT")
                .is_some_and(|value| !matches!(value.as_str(), "" | "." | "NOCALL"));
            set_comparison_format_value(
                &mut record,
                1,
                "QQ",
                if query_called { "nan" } else { "." },
            );
        }

        let mut info = base.into_values().collect::<Vec<_>>();
        let append_regions_last = regions.as_deref() == Some("Regions=TS_boundary");
        if !append_regions_last && let Some(regions) = regions.as_ref() {
            info.push(regions.clone());
        }
        if preserve_info {
            info.push(format!("RegionsExtent={}", legacy_regions_extent(&record)));
        }
        if output_vtc {
            let mut xcmp_type = comparison.decision;
            let mut xcmp_kind = comparison.kind.to_string();
            if row.xcmp_hap_match && xcmp_type != "TP" {
                xcmp_kind = format!("hapmatch__{xcmp_type}__{xcmp_kind}");
                xcmp_type = "TP".to_string();
            }
            if !info_list_values(&record.info, "Regions").contains(&"CONF") {
                xcmp_type = "UNK".to_string();
            }
            info.push(format!(
                "XCMP={xcmp_type}:{xcmp_kind}:{}:{}:{}",
                comparison.gtt1, comparison.gtt2, comparison.ctype
            ));
            let vtc = legacy_vtc(&record, &truth_fields, &query_fields);
            if !vtc.is_empty() {
                info.push(format!("VTC={vtc}"));
            }
        }
        if append_regions_last && let Some(regions) = regions {
            info.push(regions);
        }
        record.info = if info.is_empty() {
            ".".to_string()
        } else {
            info.join(";")
        };
        row.record = record.into();
    }
    Ok(())
}

pub(super) fn info_fields_by_key(info: &str) -> BTreeMap<String, String> {
    info.split(';')
        .filter(|field| !matches!(*field, "" | "."))
        .map(|field| {
            let key = field.split_once('=').map_or(field, |(key, _)| key);
            (key.to_string(), field.to_string())
        })
        .collect()
}

/// Internal compare rows have already received synthetic truth-set membership
/// tags. Legacy hands qfy the pre-quantification stream instead, so remove
/// those provisional tags from the private re-quantification handoff and let
/// qfy derive them from the final confidence and stratification inputs.
#[cfg(test)]
pub(super) fn sanitize_requantify_handoff_rows(rows: &[AnnotatedRow]) -> Vec<AnnotatedRow> {
    rows.iter()
        .cloned()
        .map(sanitize_requantify_handoff_row)
        .collect()
}

pub(super) fn sanitize_requantify_handoff_row(mut row: AnnotatedRow) -> AnnotatedRow {
    row.record
        .try_update(|record| {
            let entries = record
                .info
                .split(';')
                .filter_map(|entry| {
                    let Some(regions) = entry.strip_prefix("Regions=") else {
                        return Some(entry.to_string());
                    };
                    let retained = regions
                        .split(',')
                        .filter(|tag| !matches!(*tag, "TS_boundary" | "TS_contained"))
                        .collect::<Vec<_>>();
                    (!retained.is_empty()).then(|| format!("Regions={}", retained.join(",")))
                })
                .collect::<Vec<_>>();
            record.info = if entries.is_empty() {
                ".".to_string()
            } else {
                entries.join(";")
            };
            Ok(())
        })
        .expect("sanitizing INFO preserves checked record invariants");
    row
}

pub(super) fn set_comparison_format_value(
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

pub(super) struct LegacyComparison {
    decision: String,
    kind: String,
    gtt1: String,
    gtt2: String,
    ctype: &'static str,
    iqq: String,
}

impl LegacyComparison {
    fn preserved_fields(&self, hap_match: bool) -> Vec<(&'static str, String)> {
        let mut fields = vec![
            ("ctype", format!("ctype={}", self.ctype)),
            ("kind", format!("kind={}", self.kind)),
            ("type", format!("type={}", self.decision)),
        ];
        if self.gtt1 != "." {
            fields.push(("gtt1", format!("gtt1={}", self.gtt1)));
        }
        if self.gtt2 != "." {
            fields.push(("gtt2", format!("gtt2={}", self.gtt2)));
        }
        if hap_match {
            fields.push(("HapMatch", "HapMatch".to_string()));
        }
        fields
    }
}

pub(super) fn legacy_comparison_fields(
    record: &RawVcfRecord,
    truth: &BTreeMap<String, String>,
    query: &BTreeMap<String, String>,
    ctype: &'static str,
) -> LegacyComparison {
    let truth_called = sample_is_called(truth);
    let query_called = sample_is_called(query);
    let (decision, kind) = match (truth_called, query_called) {
        (true, false) => ("FN", "missing"),
        (false, true) => ("FP", "missing"),
        (false, false) => ("N", "match"),
        (true, true) => {
            let truth_bd = truth.get("BD").map(String::as_str).unwrap_or(".");
            let query_bd = query.get("BD").map(String::as_str).unwrap_or(".");
            let bk = query
                .get("BK")
                .or_else(|| truth.get("BK"))
                .map(String::as_str)
                .unwrap_or(".");
            if truth_bd == "FN" || query_bd == "FP" || bk == "am" || bk == "lm" {
                ("FP", legacy_mismatch_kind(record, truth, query))
            } else {
                ("TP", "match")
            }
        }
    };
    let iqq = query
        .get("QQ")
        .filter(|value| !matches!(value.as_str(), "" | "."))
        .cloned()
        .unwrap_or_else(|| "0".to_string());
    LegacyComparison {
        decision: decision.to_string(),
        kind: kind.to_string(),
        gtt1: legacy_gt_label(truth, truth_called),
        gtt2: legacy_gt_label(query, query_called),
        ctype,
        iqq,
    }
}

pub(super) fn sample_is_called(sample: &BTreeMap<String, String>) -> bool {
    sample
        .get("BVT")
        .is_some_and(|value| !matches!(value.as_str(), "" | "." | "NOCALL" | "HOMREF"))
}

pub(super) fn legacy_gt_label(sample: &BTreeMap<String, String>, called: bool) -> String {
    if !called {
        return ".".to_string();
    }
    sample
        .get("BLT")
        .filter(|value| !matches!(value.as_str(), "" | "." | "nocall"))
        .map(|value| format!("gt_{value}"))
        .unwrap_or_else(|| "gt_unknown".to_string())
}

pub(super) fn legacy_mismatch_kind<'a>(
    _record: &RawVcfRecord,
    truth: &'a BTreeMap<String, String>,
    query: &'a BTreeMap<String, String>,
) -> &'static str {
    let truth_gt = truth.get("GT").map(String::as_str).unwrap_or(".");
    let query_gt = query.get("GT").map(String::as_str).unwrap_or(".");
    let truth_alleles = parse_gt_alleles(truth_gt)
        .into_iter()
        .collect::<BTreeSet<_>>();
    let query_alleles = parse_gt_alleles(query_gt)
        .into_iter()
        .collect::<BTreeSet<_>>();
    if truth_alleles == query_alleles {
        return "gtmismatch";
    }
    let truth_nonref = truth_alleles
        .into_iter()
        .filter(|allele| *allele > 0)
        .collect::<BTreeSet<_>>();
    let query_nonref = query_alleles
        .into_iter()
        .filter(|allele| *allele > 0)
        .collect::<BTreeSet<_>>();
    if truth_nonref == query_nonref {
        "gtmismatch"
    } else if !truth_nonref.is_disjoint(&query_nonref) {
        "alpartial"
    } else {
        "almismatch"
    }
}

pub(super) fn legacy_regions_extent(record: &RawVcfRecord) -> String {
    let variant = Variant {
        key: VariantKey {
            chrom: record.chrom.clone(),
            pos: record.pos,
            ref_allele: record.ref_allele.clone(),
            alt_allele: record.alt_allele.clone(),
        },
        qual: record.qual.clone(),
        filter: record.filter.clone(),
        gt: ".".to_string(),
    };
    effective_refrange(&variant)
        .map(|(start, end, _)| format!("{start}-{end}"))
        .unwrap_or_else(|| {
            format!(
                "{}-{}",
                record.pos,
                record.pos + record.ref_allele.len().saturating_sub(1)
            )
        })
}

pub(super) fn legacy_vtc(
    record: &RawVcfRecord,
    truth: &BTreeMap<String, String>,
    query: &BTreeMap<String, String>,
) -> String {
    let mut types = BTreeMap::<u8, String>::new();
    for sample in [truth, query] {
        if !sample_is_called(sample) {
            types.insert(0x80, "nocall__nc".to_string());
            continue;
        }
        let gt = sample.get("GT").map(String::as_str).unwrap_or(".");
        let gt_alleles = parse_gt_alleles(gt);
        let mut allele_bits = 0u8;
        for allele in gt_alleles.iter().copied().filter(|allele| *allele > 0) {
            let Some(alt) = record.alt_allele.split(',').nth(allele - 1) else {
                continue;
            };
            let bits = allele_edit_bits(&record.ref_allele, alt);
            allele_bits |= bits;
            for bit in [1u8, 2, 4] {
                if bits & bit != 0 {
                    types.insert(bit, format!("nuc__{}", legacy_type_bits(bit)));
                }
            }
            if bits != 0 {
                types.insert(0x10 | bits, format!("al__{}", legacy_type_bits(bits)));
            }
        }
        if allele_bits == 0 {
            continue;
        }
        let location = match sample.get("BLT").map(String::as_str).unwrap_or("") {
            "het" => 0x30,
            "hetalt" => 0x40,
            "homalt" => 0x90,
            "hemi" => 0x50,
            _ => 0xa0,
        };
        let ref_bit = u8::from(gt_alleles.contains(&0)) * 8;
        types.insert(
            location | ref_bit | allele_bits,
            format!(
                "{}__{}",
                sample.get("BLT").map(String::as_str).unwrap_or("unknown"),
                legacy_type_bits(ref_bit | allele_bits)
            ),
        );
    }
    types.into_values().collect::<Vec<_>>().join(",")
}

pub(super) fn decorate_existing_comparison_vcf(
    output_path: &Path,
    truth_path: &Path,
    query_path: &Path,
    preserve_info: bool,
    output_vtc: bool,
) -> Result<()> {
    let reader = vcf::open_validated_vcf(output_path)?;
    let mut headers = reader.headers().to_vec();
    if output_vtc {
        let mut chrom_index = headers
            .iter()
            .position(|line| line.starts_with("#CHROM"))
            .unwrap_or(headers.len());
        for declaration in [
            "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">",
            "##INFO=<ID=XCMP,Number=.,Type=String,Description=\"XCMP extra information.\">",
        ] {
            let identity = preprocess::structured_header_identity(declaration);
            let present = headers
                .iter()
                .any(|line| preprocess::structured_header_identity(line) == identity);
            if !present {
                headers.insert(chrom_index, declaration.to_string());
                chrom_index += 1;
            }
        }
    }
    let mut decorations = DecorationIndex::default();
    for path in [truth_path, query_path] {
        for record in vcf::open_validated_vcf(path)? {
            let record = record?;
            decorations.observe(record.raw(), preserve_info, "QUAL");
        }
    }
    let decorated = reader.map(|record| {
        let record = record?;
        let raw = record.raw();
        let mut row = AnnotatedRow {
            sort_key: (raw.chrom.clone(), raw.pos, 0, 0),
            record: raw.clone().into(),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        };
        decorate_output_rows_with_index(
            std::slice::from_mut(&mut row),
            &decorations,
            preserve_info,
            output_vtc,
            "QUAL",
        )?;
        Ok(row.record.into_validated())
    });
    vcf::write_validated_vcf_iter(output_path, &headers, decorated)
}

pub(super) fn build_vcf_headers(
    truth_headers: &[String],
    query_headers: &[String],
    apply_filters_query: bool,
    output_vtc: bool,
    preserve_info: bool,
    roc_field: &str,
) -> Vec<String> {
    let mut supplied: Vec<String> = truth_headers
        .iter()
        .chain(query_headers)
        .filter(|line| line.starts_with("##"))
        .cloned()
        .collect();
    supplied.extend([
        "##INFO=<ID=gtt1,Number=1,Type=String,Description=\"GT of truth call\">".to_string(),
        "##INFO=<ID=gtt2,Number=1,Type=String,Description=\"GT of query call\">".to_string(),
        "##INFO=<ID=type,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">".to_string(),
        "##INFO=<ID=kind,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">".to_string(),
        "##INFO=<ID=ctype,Number=1,Type=String,Description=\"Type of comparison performed\">".to_string(),
        "##INFO=<ID=HapMatch,Number=0,Type=Flag,Description=\"Variant is in matching haplotype block\">".to_string(),
        "##INFO=<ID=BS,Number=1,Type=Integer,Description=\"Start position of the benchmarking superlocus on current chromosome\">".to_string(),
        format!("##INFO=<ID=IQQ,Number=1,Type=Float,Description=\"Quality value for query variants ({roc_field}).\">")
    ]);
    if apply_filters_query {
        supplied.push(
            "##INFO=<ID=Q_FILTERED,Number=0,Type=Flag,Description=\"Filtered call in query\">"
                .to_string(),
        );
    }
    supplied
        .push("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string());

    // xcmp first writes through VariantWriter (base + sorted merged input
    // headers). Quantify then appends its region/FORMAT declarations in this
    // exact order before writing the final two-sample VCF.
    let mut headers = preprocess::canonicalize_legacy_headers(&supplied);
    let chrom = headers.pop().expect("comparison header has #CHROM line");
    let mut quantified_headers =
        vec!["##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">"];
    if preserve_info {
        quantified_headers.push(
            "##INFO=<ID=RegionsExtent,Number=.,Type=String,Description=\"Trimmed reference coordinates matched to regions for this record.\">",
        );
    }
    quantified_headers.extend([
        "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision for call (TP/FP/FN/N)\">",
        "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Sub-type for decision (match/mismatch type)\">",
        "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">",
        "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">",
        "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"High-level variant type (SNP|INDEL).\">",
        "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"High-level location type (het|homref|hetalt|homalt|nocall).\">",
    ]);
    if output_vtc {
        quantified_headers.extend([
            "##INFO=<ID=VTC,Number=.,Type=String,Description=\"Variant types used for counting.\">",
            "##INFO=<ID=XCMP,Number=.,Type=String,Description=\"XCMP extra information.\">",
        ]);
    }
    for line in quantified_headers {
        let identity = preprocess::structured_header_identity(line);
        let already_present = identity.as_ref().is_some_and(|wanted| {
            headers.iter().any(|existing| {
                preprocess::structured_header_identity(existing).as_ref() == Some(wanted)
            })
        });
        if !already_present {
            headers.push(line.to_string());
        }
    }
    headers.push(chrom);
    headers
}

#[cfg(test)]
pub(super) fn build_clusters(truth: &[Variant], query: &[Variant]) -> Vec<Cluster> {
    build_clusters_with_gap(truth, query, CLUSTER_GAP_BP)
}
