//! Cohesive ROC contributions responsibility.

use super::*;

pub(super) fn emit_contributions_with_options<
    F: FnMut(RowKey, Option<f64>, &Cumul, &[String], Option<&str>, Option<&str>),
>(
    row: &AnnotatedRow,
    options: &RocOptions,
    mut emit: F,
) {
    let fields: Vec<&str> = row.record.split('\t').collect();
    if fields.len() < 11 {
        return;
    }

    let info = fields[7];
    let subsets = extract_subsets(info);

    let format_keys: Vec<&str> = fields[8].split(':').collect();
    let truth_parts: Vec<&str> = fields[9].split(':').collect();
    let query_parts: Vec<&str> = fields[10].split(':').collect();
    let truth = Sample::new(&format_keys, &truth_parts);
    let query = Sample::new(&format_keys, &query_parts);

    // Per-side ROC threshold: read FORMAT.QQ from each sample. The
    // record-level QUAL column would lose multi-allelic indels whose
    // xcmp output sets QUAL=0 while still writing the matched per-side
    // quality into FORMAT.QQ.
    let score_field = options.score_field.as_deref().unwrap_or(&options.qq_field);
    let truth_qq = truth.roc_value(score_field, fields[5], info);
    let query_qq = query.roc_value(score_field, fields[5], info);

    let filter_tags = fields[6]
        .split(';')
        .filter(|tag| !tag.is_empty() && *tag != "." && *tag != "PASS")
        .collect::<Vec<_>>();
    let ignores =
        |tag: &str| options.ignored_filters.contains("*") || options.ignored_filters.contains(tag);
    let selectively_filtered = filter_tags.iter().any(|tag| !ignores(tag));

    // Selective filtering adds the legacy SEL tier: filters named by
    // --roc-filter are ignored for SEL, while PASS continues to apply every
    // record filter and ALL applies none.
    let mut tiers = vec![("ALL", false), ("PASS", !row.query_pass)];
    if !options.ignored_filters.is_empty() {
        tiers.push(("SEL", selectively_filtered));
    }
    for (filter, filtered_out) in tiers {
        // Truth side
        if let Some(truth_bvt) = truth.variant_type() {
            let effective = if filtered_out && truth.bd == Some("TP") {
                Some("FN")
            } else {
                truth.bd
            };
            if matches!(effective, Some("TP") | Some("FN")) {
                let bucket = sample_bucket(&truth);
                let mut counts = Cumul::default();
                if effective == Some("TP") {
                    counts.truth_tp = bucket;
                } else {
                    counts.truth_fn = bucket;
                }
                emit_for_axes(
                    AxisContribution {
                        ty: truth_bvt,
                        bi: truth.bi,
                        blt: truth.blt,
                        subsets: &subsets,
                        filter,
                        qq_field: &options.qq_field,
                        qq: truth_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
        }

        // Query side. When the row passes the filter, contribute the
        // BD-typed bucket. When it fails the filter (PASS tier only),
        // emit a zero-counts phantom at every (subtype, subset) the
        // record's BI / region tags imply — legacy's BlockQuantify
        // demotes filter-failed obs but still adds them to the ROC
        // vector with the obs's flag bits, so the per-subtype level
        // bucket exists at the demoted level even though its TTP/QTP
        // counts come entirely from cumulative-above passing records.
        if let Some(query_bvt) = query.variant_type() {
            if !filtered_out && matches!(query.bd, Some("TP") | Some("FP") | Some("UNK")) {
                let bucket = sample_bucket(&query);
                let mut counts = Cumul::default();
                match query.bd {
                    Some("TP") => counts.query_tp = bucket,
                    Some("FP") => {
                        counts.query_fp = bucket;
                        if row.fp_class == Some("gt") {
                            counts.fp_gt = 1;
                        } else if row.fp_class == Some("al") {
                            counts.fp_al = 1;
                        }
                    }
                    Some("UNK") => counts.query_unk = bucket,
                    _ => {}
                }
                emit_for_axes(
                    AxisContribution {
                        ty: query_bvt,
                        bi: query.bi,
                        blt: query.blt,
                        subsets: &subsets,
                        filter,
                        qq_field: &options.qq_field,
                        qq: query_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            } else if filtered_out && matches!(query.bd, Some("TP") | Some("FP") | Some("UNK")) {
                // Filter-failed query phantom (no count contribution).
                // Pass bi through phantom_bi for substat detection.
                let counts = Cumul::default();
                let phantom_bi = if query.bd == Some("TP") {
                    query.bi
                } else {
                    None
                };
                emit_for_axes(
                    AxisContribution {
                        ty: query_bvt,
                        bi: phantom_bi,
                        blt: (query.bd == Some("TP")).then_some(query.blt).flatten(),
                        subsets: &subsets,
                        filter,
                        qq_field: &options.qq_field,
                        qq: query_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
        }
    }

    // Per-Filter rows: legacy emits one (Type, Subtype, Subset='*',
    // Filter=<tag>, QQ='*') row per non-PASS filter tag. Each query
    // record with FILTER=<tag> contributes its BD bucket; truth-side
    // TPs whose paired query carries that filter contribute too. The
    // rendering for these rows zeros TRUTH.TOTAL/FN, QUERY.TOTAL and
    // METRIC.* (handled in render_row); only the *.TP / *.FP / *.UNK
    // buckets carry real counts. The filter column is fields[6]; an
    // empty/`.`/`PASS` filter contributes nothing here (already covered
    // by ALL/PASS tiers above).
    let star_subset = ["*".to_string()];
    if !filter_tags.is_empty() {
        for tag in filter_tags {
            let output_tag = if ignores(tag) {
                format!("SEL_IGN_{tag}")
            } else {
                tag.to_string()
            };
            // Truth side: TP only (FN truths have no query → no filter).
            if let Some(truth_bvt) = truth.variant_type()
                && truth.bd == Some("TP")
            {
                let bucket = sample_bucket(&truth);
                let counts = Cumul {
                    truth_tp: bucket,
                    ..Cumul::default()
                };
                emit_for_axes(
                    AxisContribution {
                        ty: truth_bvt,
                        bi: truth.bi,
                        blt: truth.blt,
                        subsets: &star_subset,
                        filter: &output_tag,
                        qq_field: &options.qq_field,
                        qq: truth_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
            // Query side: all BD types contribute.
            if let Some(query_bvt) = query.variant_type()
                && matches!(query.bd, Some("TP") | Some("FP") | Some("UNK"))
            {
                let bucket = sample_bucket(&query);
                let mut counts = Cumul::default();
                match query.bd {
                    Some("TP") => counts.query_tp = bucket,
                    Some("FP") => {
                        counts.query_fp = bucket;
                        if row.fp_class == Some("gt") {
                            counts.fp_gt = 1;
                        } else if row.fp_class == Some("al") {
                            counts.fp_al = 1;
                        }
                    }
                    Some("UNK") => counts.query_unk = bucket,
                    _ => {}
                }
                emit_for_axes(
                    AxisContribution {
                        ty: query_bvt,
                        bi: query.bi,
                        blt: query.blt,
                        subsets: &star_subset,
                        filter: &output_tag,
                        qq_field: &options.qq_field,
                        qq: query_qq,
                        counts: &counts,
                    },
                    &mut emit,
                );
            }
        }
    }
}

pub(super) struct AxisContribution<'a> {
    pub(super) ty: &'a str,
    pub(super) bi: Option<&'a str>,
    pub(super) blt: Option<&'a str>,
    pub(super) subsets: &'a [String],
    pub(super) filter: &'a str,
    pub(super) qq_field: &'a str,
    pub(super) qq: Option<f64>,
    pub(super) counts: &'a Cumul,
}

pub(super) fn emit_for_axes<
    F: FnMut(RowKey, Option<f64>, &Cumul, &[String], Option<&str>, Option<&str>),
>(
    axes: AxisContribution<'_>,
    emit: &mut F,
) {
    let subtypes = compute_subtypes(axes.ty, axes.bi);
    for subtype in &subtypes {
        for subset in axes.subsets {
            let key =
                RowKey::new_with_qq_field(axes.ty, subtype, subset, axes.filter, axes.qq_field);
            emit(key, axes.qq, axes.counts, &subtypes, axes.bi, axes.blt);
        }
    }
}

pub(super) fn compute_subtypes(variant_type: &str, bi: Option<&str>) -> Vec<String> {
    let mut out = vec!["*".to_string()];
    if variant_type == "INDEL"
        && let Some(bi) = bi
    {
        // Multi-allelic INDELs emit comma-joined BI (e.g. `i1_5,i6_15`).
        // Each indel-class primitive contributes to its own subtype bucket.
        // The `ti`/`tv` decoration on complex INDELs is SNP-side and does
        // not spawn an extra bucket here — same fanout rule extended.csv
        // uses (compare::SampleView::variant_type_and_subtypes).
        for tok in bi.split(',') {
            if matches!(tok, "ti" | "tv") {
                continue;
            }
            let upper = tok.to_uppercase();
            if INDEL_SUBTYPES.contains(&upper.as_str()) && !out.contains(&upper) {
                out.push(upper);
            }
        }
    }
    out
}

pub(super) fn sample_bucket(sample: &Sample<'_>) -> CountsBucket {
    let mut bucket = CountsBucket {
        total: 1,
        ..CountsBucket::default()
    };
    if let Some("SNP") = sample.bvt
        && let Some(bi) = sample.bi
    {
        // Legacy `roc::makeObservationFlags` splits BI on `,` and OR's
        // the matching flags. A multi-allelic SNP with bi="ti,tv"
        // therefore counts as both a transition AND a transversion in
        // the per-substat sweeps. Mirror that by setting both bits.
        for tok in bi.split(',') {
            match tok {
                "ti" => bucket.ti = 1,
                "tv" => bucket.tv = 1,
                _ => {}
            }
        }
    }
    // Mirror compare::add_sample_stats's general het/homalt rule so ROC
    // counts agree with extended.csv: het = exactly one allele equals the
    // reference index 0 (covers 0|2, 0/3, …); homalt = both alleles equal
    // and non-zero (1/1, 2|2, 3/3, …). Hetalt (1|2, 2/1, …) lands in
    // neither bucket. Literal "0/1" / "1/1" matching under-counted truth-
    // side multi-allelic GTs that bcftools merge preserves verbatim.
    if let Some(gt) = sample.gt {
        let alleles: Vec<Option<usize>> = gt
            .split(['/', '|'])
            .map(|part| part.parse::<usize>().ok())
            .collect();
        if let [Some(left), Some(right)] = alleles.as_slice() {
            let zero_count = usize::from(*left == 0) + usize::from(*right == 0);
            if zero_count == 1 {
                bucket.het = 1;
            } else if *left != 0 && left == right {
                bucket.homalt = 1;
            }
        }
    }
    bucket
}

pub(super) fn extract_subsets(info: &str) -> Vec<String> {
    let mut out = vec!["*".to_string()];
    if let Some(tail) = info.split(";Regions=").nth(1) {
        // The Regions= value runs to end-of-INFO (we're already past the
        // FORMAT column, so there's no trailing field; but be defensive
        // against future INFO tags after Regions= by stopping at ';').
        let tag_str = tail.split(';').next().unwrap_or("");
        for tag in tag_str.split(',') {
            if !tag.is_empty() && tag != "CONF" && !out.iter().any(|value| value == tag) {
                out.push(tag.to_string());
            }
        }
    }
    out
}
