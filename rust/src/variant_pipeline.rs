//! VariantInput pipeline stages 5–10 in Rust.
//!
//! Port of the C++ pipeline at `src/c++/lib/variant/VariantInput.cpp`. After
//! the existing preprocessor stages 1–4 in [`crate::preprocess`] have done
//! filtering, multi-allelic ALT splitting, and basic normalization, this
//! module:
//!
//! * **Stage 5 — VariantPrimitiveSplitter** (`primitive_split`)
//!   For each multi-allelic indel record, decompose every ALT allele into
//!   its primitive SNP / pure-insertion / pure-deletion components via
//!   [`crate::align::realign_ref_var`] and emit one VCF record per
//!   primitive.  AD is projected from the original allele list down to
//!   `[AD[0], AD[i+1]]` per split; GT is canonicalized to `0/1` (het) or
//!   `1/1` (hom-alt); each primitive carries the same QUAL/FILTER/INFO as
//!   the parent.
//!
//! * **Stage 6 — VariantAlleleNormalizer 2nd pass** (folded into the
//!   anchoring step for now)
//!   Pure deletions and pure insertions emitted by stage 5 are left-anchored
//!   against the reference so they have at least one base on each side
//!   (`refpadding = true`). This re-aligns insertions to the trailing edge
//!   convention legacy uses (e.g. `CA → CAA` becomes `A → AA` at `pos + 1`).
//!
//! * **Stages 7–10** (VariantLocationAggregator, VariantAlleleUniq,
//!   VariantLeftPadding, VariantCallsOnly) — staged in for the cases that
//!   actually fire in chr21 prepy. We currently rely on the upstream
//!   filtering already done in `preprocess::run`; the aggregator merges
//!   primitives that land at the same position into a multi-allelic record.
//!
//! The split here always runs when an input record is multi-allelic with at
//! least one indel allele. Same-direction multi-allelic insertions (e.g.
//! `T → TG,TTG`) emit primitives at the same position which the aggregator
//! re-merges; mixed-direction or same-direction deletions of different
//! lengths fan out into separate per-position records. Either way, byte
//! parity with legacy comes from one code path.

use crate::align;
use crate::partial_credit::{self, RefVar};
use crate::vcf::RawVcfRecord;
use std::collections::BTreeMap;

/// Decompose any multi-allelic indel record into primitive per-position
/// records and re-aggregate primitives that land at the same position.
///
/// The input is the post-stage-4 record (already through the existing
/// preprocess loop's filter / region / GT-rewrite passes) but BEFORE any
/// `split_multi_allelic` / `apply_left_shift`. Output is a sequence of
/// records that match the legacy primitive splitter's output, each suitable
/// for the existing `insert_ado_format` + `reorder_format_fields` final
/// canonicalization.
pub fn primitive_split(record: &RawVcfRecord, reference: &[u8]) -> Vec<RawVcfRecord> {
    primitive_split_with_floor(record, reference, 0)
}

/// Decompose a record while carrying the previous record's reference end
/// into the primitive normalizer. Legacy keeps this floor across consecutive
/// records; resetting it for every multi-allelic site can slide a later site
/// behind an already-emitted neighbor and make the preprocessed VCF unsorted.
pub(crate) fn primitive_split_with_floor(
    record: &RawVcfRecord,
    reference: &[u8],
    previous_end: usize,
) -> Vec<RawVcfRecord> {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();

    // SNP-only multi-allelics (`C → A,G`) and trivial single-alt records
    // pass through unchanged. Same-direction multi-allelic insertions like
    // `T → TG,TTG` also fall here unless an allele is itself complex
    // (reflen > 1 && altlen > 1) — same-position primitives re-merge below.
    if !needs_primitive_split(record, &alts) {
        let mut passthrough = record.clone();
        // Even for passthrough records the legacy aggregator's reversed
        // het-of-alts GT order applies — see `merge_records` for the
        // derivation. Apply the swap unconditionally here so SNP-only
        // multi-allelic outputs (`C → A,G`) come out as `2/1`, matching
        // legacy's `addAlleleToVariant` ordering.
        if alts.len() > 1 {
            canonicalize_hetalt_gt(&mut passthrough);
        }
        return vec![passthrough];
    }

    // Per-allele decomposition into primitive RefVars.
    let mut primitives: Vec<(usize, Vec<RefVar>)> = Vec::with_capacity(alts.len());
    for (idx, alt) in alts.iter().enumerate() {
        primitives.push((idx, allele_primitives(record, alt, reference)));
    }

    // Build per-primitive RawVcfRecord values, projecting AD and remapping GT
    // for each split allele.
    let format_keys: Vec<String> = record
        .format
        .as_deref()
        .map(|f| f.split(':').map(|s| s.to_string()).collect())
        .unwrap_or_default();
    let ad_index = format_keys.iter().position(|k| k == "AD");
    let gt_index = format_keys.iter().position(|k| k == "GT");
    let phased = gt_index.is_some_and(|gt_index| {
        record.samples.iter().any(|sample| {
            sample
                .split(':')
                .nth(gt_index)
                .is_some_and(|gt| gt.contains('|'))
        })
    });

    let mut output: Vec<RawVcfRecord> = Vec::new();
    for (allele_idx, prims) in primitives {
        let target = (allele_idx + 1) as u32;
        for prim in prims {
            if let Some(rec) = primitive_to_record(
                record,
                &prim,
                target,
                ad_index,
                gt_index,
                &format_keys,
                reference,
            ) {
                output.push(rec);
            }
        }
    }

    // Stage 6 — leftshift each primitive with per-sample `current_maxpos`
    // tracking so siblings can't overlap after sliding. Sort first by the
    // legacy `HCall_less` order (start asc, ref_len desc, alt asc) so the
    // sample's max-position monotonically increases as primitives are
    // processed; the second primitive at the same input anchor sees the
    // first primitive's end as its leftshift floor and stays put — matching
    // legacy's `VariantAlleleNormalizer.cpp:209-238` behaviour.
    output.sort_by(|a, b| {
        a.pos.cmp(&b.pos).then_with(|| {
            b.ref_allele
                .len()
                .cmp(&a.ref_allele.len())
                .then(a.alt_allele.cmp(&b.alt_allele))
        })
    });
    let mut current_maxpos = previous_end;
    for rec in &mut output {
        shift_primitive_record(rec, reference, current_maxpos);
        // Legacy advances from `nv.pos + nv.len - 1` of the emitted record,
        // not from the wider source span that produced the primitive.
        let new_end = rec.pos + rec.ref_allele.len().saturating_sub(1);
        current_maxpos = current_maxpos.max(new_end);
    }

    // Stage 7 — aggregate primitives at the same (chrom, pos, ref) location
    // back into a single multi-allelic record (e.g. `T → TG,TTG`). We only
    // merge when the leading REF span matches; if positions differ even by
    // one base the records stay separate (matching legacy's per-position
    // primitive emission for mixed-length deletions).
    aggregate_same_position(output, &format_keys, phased)
}

/// Run partial-credit left-shift + trim on a single-ALT primitive emitted
/// by `primitive_split`. Mirrors legacy `VariantAlleleNormalizer.cpp:209-238`:
/// `leftShift → trimLeft(refpadding) → trimRight(refpadding) → `
/// `nv.len==0 ⇒ --pos, len=1` insertion fixup. The slide step is
/// load-bearing for byte parity — without it sibling primitives emitted
/// at the same input anchor (e.g. `C → CACAC,CACAT`) would re-aggregate
/// instead of decomposing into per-anchor records, and ambiguous indels
/// in homopolymer / microsatellite contexts would sit at the wrong pos.
/// `current_maxpos` is the per-batch leftshift floor: this primitive may
/// not slide left of the previous primitive's reference end.
fn shift_primitive_record(record: &mut RawVcfRecord, reference: &[u8], current_maxpos: usize) {
    if record.alt_allele.contains(',')
        || record.alt_allele == "."
        || record.alt_allele.is_empty()
        || record.alt_allele == "<NON_REF>"
    {
        return;
    }
    let ref_len = record.ref_allele.len();
    if ref_len == 0 {
        return;
    }
    let mut rv = partial_credit::RefVar {
        start: record.pos,
        end: record.pos + ref_len - 1,
        alt: record.alt_allele.clone(),
    };
    let pos_min = current_maxpos.max(1);
    partial_credit::left_shift(reference, &mut rv, pos_min, true);

    if rv.start == 0 {
        return;
    }
    let new_ref_len = rv.end as i64 - rv.start as i64 + 1;
    if new_ref_len <= 0 || rv.alt.is_empty() {
        return;
    }
    if rv.end > reference.len() {
        return;
    }
    let new_ref: String = reference[(rv.start - 1)..rv.end]
        .iter()
        .map(|b| b.to_ascii_uppercase() as char)
        .collect();
    record.pos = rv.start;
    record.ref_allele = new_ref;
    record.alt_allele = rv.alt;
}

/// Canonicalise hetalt GT entries (e.g. `1/2`) into legacy's reversed
/// `<later>/<earlier>` form (e.g. `2/1`) for every sample on this record.
/// Mirrors what `VariantLocationAggregator::addAlleleToVariant` produces
/// when it merges sibling alt half-calls into a multi-allelic record under
/// `MAX_GT = 2`. Hom and het-with-ref cases pass through unchanged.
pub fn canonicalize_hetalt_gt_public(record: &mut RawVcfRecord) {
    canonicalize_hetalt_gt(record)
}

fn canonicalize_hetalt_gt(record: &mut RawVcfRecord) {
    let format_keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|f| f.split(':').collect())
        .unwrap_or_default();
    let Some(gt_index) = format_keys.iter().position(|k| *k == "GT") else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        if let Some(cell) = cells.get_mut(gt_index) {
            *cell = swap_hetalt_gt(cell);
        }
        *sample = cells.join(":");
    }
}

/// Swap the two GT components when the call is a het-of-alts (both indices
/// greater than zero and not equal). Other shapes (homref, het-with-ref, hom-alt,
/// no-call) pass through unchanged.
fn swap_hetalt_gt(gt: &str) -> String {
    // Legacy's `VariantLocationAggregator::addAlleleToVariant` rewrites the
    // paired GT into `<later>/<earlier>` form — the aggregator only fires on
    // UNPHASED multi-allelic calls; phased (`|`) genotypes carry their own
    // haplotype assignment and are preserved verbatim by the loader. Rust
    // must mirror that: swapping `1|2` into `2|1` would break truth rows
    // like chr21:17566241 (truth source `C→CA,CAA 1|2` must stay `1|2`).
    if !gt.contains('/') {
        return gt.to_string();
    }
    let tokens: Vec<&str> = gt.split('/').collect();
    if tokens.len() != 2 {
        return gt.to_string();
    }
    let parsed: Vec<Option<u32>> = tokens.iter().map(|t| t.parse::<u32>().ok()).collect();
    if let (Some(a), Some(b)) = (parsed[0], parsed[1])
        && a > 0
        && b > 0
        && a != b
    {
        return format!("{}/{}", b, a);
    }
    gt.to_string()
}

fn needs_primitive_split(record: &RawVcfRecord, alts: &[&str]) -> bool {
    if record.alt_allele == "." || record.alt_allele.is_empty() {
        return false;
    }
    if alts.iter().any(|a| *a == "<NON_REF>" || *a == "*") {
        return false;
    }
    let ref_len = record.ref_allele.len();
    // Trigger when (a) we're multi-allelic with ≥1 indel allele OR (b) any
    // allele is a real MNP/complex variant (both reflen>1 and altlen>1).
    let has_multi_indel = alts.len() > 1 && alts.iter().any(|a| a.len() != ref_len);
    let has_complex = alts.iter().any(|a| ref_len > 1 && a.len() > 1);
    has_multi_indel || has_complex
}

/// Split an allele's REF/ALT pair into its primitive RefVars in the legacy
/// convention. Mirrors the per-allele body of
/// `VariantPrimitiveSplitter::advance` plus its post-trim realign-or-passthrough
/// guard (`src/c++/lib/variant/VariantPrimitiveSplitter.cpp:135-178`).
fn allele_primitives(record: &RawVcfRecord, alt: &str, reference: &[u8]) -> Vec<RefVar> {
    let original = RefVar {
        start: record.pos,
        end: record.pos + record.ref_allele.len().saturating_sub(1),
        alt: alt.to_string(),
    };
    // Pre-trim a COPY with refpadding=false to mirror the realignability
    // check the legacy splitter does upfront. The trimmed copy is
    // discarded if the allele turns out not to be realignable — legacy
    // pushes the original (untrimmed) `v.variation[i]` through to the
    // downstream normalizer in that case (`VariantPrimitiveSplitter.cpp:151`).
    let mut probe = original.clone();
    partial_credit::trim_left(reference, &mut probe, false);
    partial_credit::trim_right(reference, &mut probe, false);
    let reflen = (probe.end as i64) - (probe.start as i64) + 1;
    let altlen = probe.alt.len() as i64;
    if (reflen >= 1 && altlen > 1) || (reflen > 1 && altlen >= 1) {
        // Realignable — fan out via the aligner using the trimmed bounds.
        let ref_bytes: Vec<u8> = if reflen > 0 {
            reference[(probe.start - 1)..(probe.end)].to_vec()
        } else {
            Vec::new()
        };
        align::realign_ref_var(probe.start, &ref_bytes, probe.alt.as_bytes())
    } else {
        // Pass the ORIGINAL untrimmed RefVar through — leftshift in stage 6
        // handles its trim-and-slide differently than this pre-check would.
        vec![original]
    }
}

/// Construct the VCF record for a single primitive RefVar.
///
/// Anchors pure deletions / pure insertions on their LEFT base (consuming
/// one base of reference context) so the emitted VCF row is well-formed.
/// SNP primitives stay at their reported position. AD is projected to
/// `[AD[0], AD[target_index]]` and GT canonicalised to the het / hom-alt
/// shape that matches legacy's per-allele split.
fn primitive_to_record(
    record: &RawVcfRecord,
    prim: &RefVar,
    target: u32,
    ad_index: Option<usize>,
    gt_index: Option<usize>,
    format_keys: &[String],
    reference: &[u8],
) -> Option<RawVcfRecord> {
    let mut out = record.clone();
    let reflen_i = (prim.end as i64) - (prim.start as i64) + 1;
    let altlen = prim.alt.len();

    if reflen_i > 0 && altlen == reflen_i as usize && reflen_i == 1 {
        // SNP primitive — emit as-is.
        out.pos = prim.start;
        out.ref_allele = ref_slice(reference, prim.start, prim.end)?
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();
        out.alt_allele = prim.alt.to_ascii_uppercase();
    } else if reflen_i > 0 && altlen == 0 {
        // Pure deletion — right-anchor convention (legacy trailing-edge):
        // pos -> prim.start; REF = deleted_bases + right_anchor; ALT = right_anchor.
        let anchor_pos = prim.end + 1;
        let anchor = ref_slice(reference, anchor_pos, anchor_pos)?[0].to_ascii_uppercase() as char;
        let deleted: String = ref_slice(reference, prim.start, prim.end)?
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();
        out.pos = prim.start;
        out.ref_allele = format!("{deleted}{anchor}");
        out.alt_allele = anchor.to_string();
    } else if reflen_i <= 0 && altlen > 0 {
        // Pure insertion — right-anchor convention (legacy trailing-edge):
        // pos -> prim.start; REF = right_anchor; ALT = inserted_seq + right_anchor.
        let anchor_pos = prim.start;
        let anchor = ref_slice(reference, anchor_pos, anchor_pos)?[0].to_ascii_uppercase() as char;
        out.pos = anchor_pos;
        out.ref_allele = anchor.to_string();
        out.alt_allele = format!("{}{anchor}", prim.alt.to_ascii_uppercase());
    } else if reflen_i > 0 && altlen > 0 {
        // Generic mixed-length primitive (e.g. the original untrimmed
        // RefVar passed through by the legacy splitter when no allele was
        // realignable). Emit verbatim — the leftshift step (stage 6) and
        // the writer-side anchoring will canonicalise it later.
        out.pos = prim.start;
        out.ref_allele = ref_slice(reference, prim.start, prim.end)?
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();
        out.alt_allele = prim.alt.to_ascii_uppercase();
    } else {
        // Empty primitive — nothing to emit.
        return None;
    }

    // Per-sample GT + AD projection.
    let mut new_samples = Vec::with_capacity(record.samples.len());
    for sample in &record.samples {
        let cells: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        new_samples.push(project_sample(
            &cells,
            target,
            ad_index,
            gt_index,
            format_keys.len(),
        ));
    }
    out.samples = new_samples;
    Some(out)
}

fn ref_slice(reference: &[u8], start: usize, end: usize) -> Option<&[u8]> {
    if start == 0 || end < start || end > reference.len() {
        return None;
    }
    Some(&reference[start - 1..end])
}

/// Build the per-sample cell for a split allele:
///   * GT: keep allele `target` as the only non-ref alt (rewrite to `0/1`,
///     `1/1` etc. depending on input shape).
///   * AD: project to `[AD[0], AD[target]]`.
///   * All other fields: pass through unchanged.
fn project_sample(
    cells: &[String],
    target: u32,
    ad_index: Option<usize>,
    gt_index: Option<usize>,
    expected_len: usize,
) -> String {
    let mut new_cells: Vec<String> = (0..expected_len)
        .map(|i| cells.get(i).cloned().unwrap_or_else(|| ".".to_string()))
        .collect();

    if let Some(gi) = gt_index
        && let Some(cell) = new_cells.get_mut(gi)
    {
        *cell = canonical_split_gt(cell, target);
    }
    if let Some(ai) = ad_index
        && let Some(cell) = new_cells.get_mut(ai)
    {
        *cell = project_ad(cell, target);
    }
    new_cells.join(":")
}

/// Project an AD list from the original multi-allelic shape down to
/// `[ref_depth, this_alt_depth]`. Missing or malformed values default to
/// `0`. AD with only one element (no per-allele depths) is preserved as-is.
fn project_ad(ad: &str, target: u32) -> String {
    if ad == "." || ad.is_empty() {
        return ad.to_string();
    }
    let parts: Vec<&str> = ad.split(',').collect();
    if parts.len() <= 2 {
        return ad.to_string();
    }
    let target_idx = target as usize;
    let ref_depth = parts.first().copied().unwrap_or("0");
    let alt_depth = parts.get(target_idx).copied().unwrap_or("0");
    format!("{ref_depth},{alt_depth}")
}

/// Canonicalise a split-allele GT to the legacy form. Matches what the
/// legacy primitive splitter emits: het calls become `0/1`, hom-alt calls
/// become `1/1`, no-call passes through unchanged.
fn canonical_split_gt(gt: &str, target: u32) -> String {
    if gt == "." || gt == "./." || gt == ".|." || gt.is_empty() {
        return gt.to_string();
    }
    let separator = if gt.contains('|') { '|' } else { '/' };
    let mut tokens: Vec<i32> = Vec::new();
    for tok in gt.split(['/', '|']) {
        match tok.parse::<i32>() {
            Ok(v) => tokens.push(v),
            Err(_) => return gt.to_string(),
        }
    }
    let target_i = target as i32;
    let count_target = tokens.iter().filter(|t| **t == target_i).count();
    let count_ref = tokens.iter().filter(|t| **t == 0).count();
    let total = tokens.len();
    if count_target == 0 {
        // Allele not called for this sample — emit homref so the record is
        // still valid VCF; legacy emits a no-call here but the record gets
        // dropped at write-time anyway.
        return vec!["0"; total].join(&separator.to_string());
    }
    if count_target == total {
        return vec!["1"; total].join(&separator.to_string());
    }
    if count_ref + count_target == total || count_target < total {
        // Het call (target + ref OR target + sibling alt that we fold to
        // ref): emit `0/1` in canonical low-allele-first order to match
        // legacy. The original separator (`/` vs `|`) is preserved.
        return ["0", "1"].join(&separator.to_string());
    }
    // Fall-through (shouldn't trigger given the prior branches) — pass GT
    // back unchanged.
    gt.to_string()
}

/// Stage 7 — aggregate primitives that landed at the same anchor position
/// back into a single multi-allelic record. ALT alleles are collected in
/// (length, alphabetical) order so the output matches legacy's canonical
/// `T → TG,TTG` shape.
fn aggregate_same_position(
    records: Vec<RawVcfRecord>,
    format_keys: &[String],
    phased: bool,
) -> Vec<RawVcfRecord> {
    if records.len() <= 1 {
        return records;
    }
    let ad_index = format_keys.iter().position(|k| k == "AD");
    let gt_index = format_keys.iter().position(|k| k == "GT");

    // Group by (chrom, pos, ref_allele).
    let mut groups: BTreeMap<(String, usize, String), Vec<RawVcfRecord>> = BTreeMap::new();
    let mut order: Vec<(String, usize, String)> = Vec::new();
    for rec in records {
        let key = (rec.chrom.clone(), rec.pos, rec.ref_allele.clone());
        if !groups.contains_key(&key) {
            order.push(key.clone());
        }
        groups.entry(key).or_default().push(rec);
    }

    let mut out = Vec::with_capacity(order.len());
    for key in order {
        let group = groups.remove(&key).unwrap();
        if group.len() == 1 {
            out.push(group.into_iter().next().unwrap());
            continue;
        }
        let has_substitution = group
            .iter()
            .any(|record| record.ref_allele.len() == record.alt_allele.len());
        let has_indel = group
            .iter()
            .any(|record| record.ref_allele.len() != record.alt_allele.len());
        if phased && has_substitution && has_indel {
            let mut split = group;
            split.sort_by(|left, right| {
                let left_indel = left.ref_allele.len() != left.alt_allele.len();
                let right_indel = right.ref_allele.len() != right.alt_allele.len();
                left_indel
                    .cmp(&right_indel)
                    .then(left.alt_allele.cmp(&right.alt_allele))
            });
            out.extend(split);
            continue;
        }
        // Sort alts by (length, alphabetical) — matches legacy's canonical
        // multi-allelic ordering.
        let mut sorted = group;
        sorted.sort_by(|a, b| {
            a.alt_allele
                .len()
                .cmp(&b.alt_allele.len())
                .then_with(|| a.alt_allele.cmp(&b.alt_allele))
        });
        let merged = merge_records(sorted, ad_index, gt_index);
        out.push(merged);
    }
    out
}

fn merge_records(
    sorted: Vec<RawVcfRecord>,
    ad_index: Option<usize>,
    gt_index: Option<usize>,
) -> RawVcfRecord {
    let mut base = sorted[0].clone();
    let alt_list: Vec<String> = sorted.iter().map(|r| r.alt_allele.clone()).collect();
    base.alt_allele = alt_list.join(",");

    // Merge per-sample GT and AD across the alleles. We assume each input
    // record was a single-alt projection (GT=0/1 or 1/1, AD=[ref, this]).
    let n_samples = base.samples.len();
    let mut new_samples = Vec::with_capacity(n_samples);
    for s in 0..n_samples {
        let mut cells: Vec<String> = base
            .samples
            .get(s)
            .cloned()
            .unwrap_or_default()
            .split(':')
            .map(|c| c.to_string())
            .collect();

        // Combined AD = [ref_depth, alt1_depth, alt2_depth, ...] taking
        // ref_depth from the first split and per-alt depths from each.
        if let Some(ai) = ad_index {
            let mut ad_parts: Vec<String> = Vec::new();
            // ref depth from the first record
            let first_ad = sorted[0]
                .samples
                .get(s)
                .map(|c| c.split(':').nth(ai).unwrap_or(".").to_string())
                .unwrap_or_default();
            let first_parts: Vec<&str> = first_ad.split(',').collect();
            ad_parts.push(first_parts.first().copied().unwrap_or("0").to_string());
            for rec in &sorted {
                let ad = rec
                    .samples
                    .get(s)
                    .map(|c| c.split(':').nth(ai).unwrap_or(".").to_string())
                    .unwrap_or_default();
                let parts: Vec<&str> = ad.split(',').collect();
                ad_parts.push(parts.get(1).copied().unwrap_or("0").to_string());
            }
            if let Some(cell) = cells.get_mut(ai) {
                *cell = ad_parts.join(",");
            }
        }

        // Combined GT: choose the highest-positional combo from any sample
        // call. For a het+het pair (0/1, 0/1) the merged form is 1/2.
        if let Some(gi) = gt_index {
            let mut alt_calls: Vec<u32> = Vec::new();
            for (i, rec) in sorted.iter().enumerate() {
                let gt_cell = rec
                    .samples
                    .get(s)
                    .map(|c| c.split(':').nth(gi).unwrap_or(".").to_string())
                    .unwrap_or_default();
                let separator = if gt_cell.contains('|') { '|' } else { '/' };
                for tok in gt_cell.split(['/', '|']) {
                    if tok == "1" {
                        alt_calls.push((i + 1) as u32);
                    }
                }
                let _ = separator;
            }
            let separator = '/';
            let merged_gt = if alt_calls.is_empty() {
                "0/0".to_string()
            } else if alt_calls.len() == 1 {
                format!("0{separator}{}", alt_calls[0])
            } else {
                // Multiple alt calls in this sample (hetalt). Legacy's
                // `VariantLocationAggregator` merges by adding incoming
                // alleles into the LAST zero slot of `gt[]` (with MAX_GT=2
                // — see `Variant.hh`'s `#define MAX_GT 2`). Starting from
                // a het-shaped seed `[0, 1]` for the first allele, the
                // second allele's `addAlleleToVariant` finds `gt[0] == 0`
                // and writes there, producing `[<later>, <earlier>]` —
                // i.e. `2/1`, not `1/2`. Mirror that ordering here so
                // multi-allelic byte output matches.
                alt_calls.sort();
                format!("{}{separator}{}", alt_calls[1], alt_calls[0])
            };
            if let Some(cell) = cells.get_mut(gi) {
                *cell = merged_gt;
            }
        }
        new_samples.push(cells.join(":"));
    }
    base.samples = new_samples;
    base
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_record(
        chrom: &str,
        pos: usize,
        r: &str,
        a: &str,
        format: &str,
        sample: &str,
    ) -> RawVcfRecord {
        RawVcfRecord {
            chrom: chrom.to_string(),
            pos,
            id: ".".into(),
            ref_allele: r.into(),
            alt_allele: a.into(),
            qual: ".".into(),
            filter: ".".into(),
            info: ".".into(),
            format: Some(format.into()),
            samples: vec![sample.into()],
        }
    }

    // chr21 reference window covering position 17809620..=17809650 (placeholder
    // bases A..A so leftshift / anchoring math has stable inputs for the test).
    fn fake_ref() -> Vec<u8> {
        let mut v = Vec::new();
        for i in 0u32..30000 {
            v.push(if i % 2 == 0 { b'A' } else { b'C' });
        }
        v
    }

    #[test]
    fn single_allele_passes_through_unchanged() {
        let rec = make_record("chr1", 100, "A", "G", "GT", "0/1");
        let reference = fake_ref();
        let out = primitive_split(&rec, &reference);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "G");
        assert_eq!(out[0].samples[0], "0/1");
    }

    #[test]
    fn snp_only_multi_allelic_passes_through_unchanged() {
        // Same-length multi-allelic SNPs keep their multi shape — the
        // primitive splitter does not fire because no allele is realignable.
        // GT, however, gets the legacy hetalt swap (1/2 → 2/1) to match
        // the LocationAggregator output convention.
        let rec = make_record("chr21", 100, "C", "A,G", "GT:AD", "1/2:2,71,637");
        let reference = fake_ref();
        let out = primitive_split(&rec, &reference);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "A,G");
        assert_eq!(out[0].samples[0], "2/1:2,71,637");
    }

    #[test]
    fn project_ad_keeps_ref_and_target_only() {
        assert_eq!(project_ad("14,0,21", 1), "14,0");
        assert_eq!(project_ad("14,0,21", 2), "14,21");
        // Two-element AD passes through verbatim (already projected).
        assert_eq!(project_ad("14,21", 1), "14,21");
        // Missing AD passes through.
        assert_eq!(project_ad(".", 1), ".");
    }

    /// Build a 1-based reference where positions `start..start+seq.len()`
    /// contain `seq` and surrounding bases are filler `N`. Used to test
    /// primitive splitting against a known short window without loading a
    /// real FASTA.
    fn windowed_ref(start: usize, seq: &[u8]) -> Vec<u8> {
        let mut v = vec![b'N'; start + seq.len() + 1024];
        for (i, &b) in seq.iter().enumerate() {
            v[start - 1 + i] = b;
        }
        v
    }

    #[test]
    fn chr21_17809630_mixed_indel_primitive_splits_correctly() {
        // Replicates the trace `chr21:17809630 CA → CAA,C` documented in
        // PHASE1_BASELINE.md. Reference bytes at 17809630..=17809631 are
        // "CA" — the same bases the C++ pipeline reads from the real chr21
        // FASTA at this position.
        let reference = windowed_ref(17_809_625, b"NNNNNCANNNN");
        // Sanity: pos 17809630 is 'C', pos 17809631 is 'A' in the window.
        assert_eq!(reference[17_809_629], b'C');
        assert_eq!(reference[17_809_630], b'A');

        let rec = make_record("chr21", 17_809_630, "CA", "CAA,C", "GT:AD", "1/2:14,0,21");
        let mut out = primitive_split(&rec, &reference);
        // Sort by (pos, alt) so the assertions are stable regardless of
        // emission order.
        out.sort_by(|a, b| a.pos.cmp(&b.pos).then(a.alt_allele.cmp(&b.alt_allele)));

        assert_eq!(out.len(), 2);

        // The deletion record (`CA → C`) anchors at the original position.
        assert_eq!(out[0].pos, 17_809_630);
        assert_eq!(out[0].ref_allele, "CA");
        assert_eq!(out[0].alt_allele, "C");
        assert_eq!(out[0].samples[0], "0/1:14,21");

        // The insertion record (`CA → CAA`) is repositioned to the trailing
        // edge as legacy emits it: `A → AA` at position 17_809_631.
        assert_eq!(out[1].pos, 17_809_631);
        assert_eq!(out[1].ref_allele, "A");
        assert_eq!(out[1].alt_allele, "AA");
        assert_eq!(out[1].samples[0], "0/1:14,0");
    }

    #[test]
    fn consecutive_repeat_indels_carry_the_previous_record_floor() {
        // Reduced form of the adjacent chr21:30548464/30548467 records that
        // exposed an unsorted germline query.prep.vcf.gz. Resetting the
        // primitive normalizer at the second site slides it back behind the
        // first site's emitted tail; carrying the first record's end keeps
        // the stream monotonic, as legacy VariantAlleleNormalizer does.
        let reference = b"ctctcctctctctctctctctctctctctctc";
        let first = make_record("chr21", 5, "CCTCT", "CTCTCT,C", "GT:AD", "1/2:19,1,13");
        let second = make_record("chr21", 8, "CTCTC", "C,CCTC", "GT:AD", "1/2:33,0,0");

        let first_out = primitive_split(&first, reference);
        let last_first_position = first_out.iter().map(|record| record.pos).max().unwrap();
        let reset_second = primitive_split(&second, reference);
        assert!(
            reset_second
                .iter()
                .any(|record| record.pos < last_first_position),
            "the reduced fixture must exercise the historical backwards slide"
        );

        let carried_second = primitive_split_with_floor(&second, reference, 9);
        assert!(
            carried_second
                .iter()
                .all(|record| record.pos >= last_first_position),
            "carrying the previous REF end must preserve coordinate order"
        );
    }

    #[test]
    fn snp_only_multi_allelic_passes_with_legacy_gt_order() {
        // SNP-only multi-allelic input: `C → A,G GT=1/2`. Legacy's
        // aggregator emits `GT=2/1` (later allele first) due to MAX_GT=2.
        // Our passthrough path applies the same canonicalisation so the
        // record matches byte-for-byte.
        let reference = windowed_ref(10_716_540, b"NNNNNCNNNN");
        let rec = make_record("chr21", 10_716_541, "C", "A,G", "GT:AD", "1/2:2,71,637");
        let out = primitive_split(&rec, &reference);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].samples[0], "2/1:2,71,637");
    }

    #[test]
    fn swap_hetalt_gt_only_swaps_unphased_distinct_alts() {
        assert_eq!(swap_hetalt_gt("1/2"), "2/1");
        assert_eq!(swap_hetalt_gt("2/1"), "1/2");
        // Phased GTs are never swapped — their haplotype order is meaningful
        // (hap1 on the left of `|`) and legacy preserves it verbatim.
        assert_eq!(swap_hetalt_gt("1|2"), "1|2");
        assert_eq!(swap_hetalt_gt("2|1"), "2|1");
        // Het-with-ref left alone.
        assert_eq!(swap_hetalt_gt("0/1"), "0/1");
        assert_eq!(swap_hetalt_gt("1/0"), "1/0");
        // Hom-alt left alone.
        assert_eq!(swap_hetalt_gt("1/1"), "1/1");
        // No-call left alone.
        assert_eq!(swap_hetalt_gt("./."), "./.");
    }

    #[test]
    fn t_tg_ttg_aggregates_into_canonical_multi_allelic() {
        // Replicates `chr21:11101386 T → TTG,TG` (input order) which legacy
        // emits as `T → TG,TTG` after primitive split + aggregation.
        let reference = windowed_ref(11_101_380, b"NNNNNNTNNNN");
        assert_eq!(reference[11_101_385], b'T');
        let rec = make_record("chr21", 11_101_386, "T", "TTG,TG", "GT:AD", "1/2:2,36,137");
        let out = primitive_split(&rec, &reference);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].pos, 11_101_386);
        assert_eq!(out[0].ref_allele, "T");
        assert_eq!(out[0].alt_allele, "TG,TTG");
        // After re-aggregation: AD becomes [ref, TG_depth, TTG_depth] in the
        // sorted-by-length order. Original AD was [2, 36, 137] in the order
        // [ref, TTG, TG]; sorted = [ref, TG, TTG] = [2, 137, 36].
        // GT in legacy hetalt order: later allele first (`2/1` not `1/2`).
        assert_eq!(out[0].samples[0], "2/1:2,137,36");
    }

    #[test]
    fn phased_mixed_snp_and_indel_remain_separate() {
        let reference = windowed_ref(100, b"NNNNNANNNNN");
        let rec = make_record("chr1", 105, "A", "AT,T", "GT:AD", "1|2:5,7,9");
        let out = primitive_split(&rec, &reference);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "T");
        assert_eq!(out[1].alt_allele, "AT");
    }

    #[test]
    fn primitive_floor_advances_from_emitted_span_not_source_span() {
        // Reduced directly from chr21:16528119 in the preprocess-controls
        // truth fixture. Advancing the sibling floor to the original TCA
        // span moves the second insertion one base too far right.
        let reference = windowed_ref(16_528_100, b"GAAATGAAGGTCAGACCCGTCACACACACACACACACACAC");
        let rec = make_record(
            "chr21",
            16_528_119,
            "TCA",
            "TCACA,TCACACACA",
            "GT:AD",
            "2|1:.,.,.",
        );
        let out = primitive_split(&rec, &reference);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].pos, 16_528_119);
        assert_eq!(out[0].ref_allele, "T");
        assert_eq!(out[0].alt_allele, "TCA");
        assert_eq!(out[1].pos, 16_528_120);
        assert_eq!(out[1].ref_allele, "C");
        assert_eq!(out[1].alt_allele, "CACACAC");
    }

    #[test]
    fn aagag_to_a_aag_splits_into_per_position_records() {
        // `chr21:14728309 AAGAG → A,AAG` - both deletions, but at different
        // primitive positions. Legacy emits two separate records:
        //   chr21 14728309 AAGAG → A    (4-base deletion)
        //   chr21 14728311 GAG   → G    (2-base deletion at +2)
        let reference = windowed_ref(14_728_300, b"NNNNNNNNNAAGAGNN");
        assert_eq!(&reference[14_728_308..14_728_313], b"AAGAG");
        let rec = make_record(
            "chr21",
            14_728_309,
            "AAGAG",
            "A,AAG",
            "GT:AD",
            "1/2:15,12,0",
        );
        let mut out = primitive_split(&rec, &reference);
        out.sort_by(|a, b| a.pos.cmp(&b.pos).then(a.alt_allele.cmp(&b.alt_allele)));

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].pos, 14_728_309);
        assert_eq!(out[0].ref_allele, "AAGAG");
        assert_eq!(out[0].alt_allele, "A");
        assert_eq!(out[0].samples[0], "0/1:15,12");

        assert_eq!(out[1].pos, 14_728_311);
        assert_eq!(out[1].ref_allele, "GAG");
        assert_eq!(out[1].alt_allele, "G");
        assert_eq!(out[1].samples[0], "0/1:15,0");
    }

    #[test]
    fn canonical_split_gt_normalises_het_to_zero_one() {
        // Het with target allele 1: 1/2 -> 0/1
        assert_eq!(canonical_split_gt("1/2", 1), "0/1");
        // Het with target allele 2: 1/2 -> 0/1 (other alt folds to ref)
        assert_eq!(canonical_split_gt("1/2", 2), "0/1");
        // Hom-alt of target stays homozygous.
        assert_eq!(canonical_split_gt("1/1", 1), "1/1");
        // Phased separator preserved.
        assert_eq!(canonical_split_gt("1|2", 1), "0|1");
        // No-call passes through.
        assert_eq!(canonical_split_gt("./.", 1), "./.");
    }
}
