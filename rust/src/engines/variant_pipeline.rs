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
//!   `[AD[0], AD[i+1]]` per split; GT is emitted unphased as `0/1` for a
//!   passthrough het, `1/0` for a realigned het, or `1/1` for hom-alt; each
//!   primitive carries the same QUAL/FILTER/INFO as the parent.
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

use crate::domain::{PrimitiveIdentity, RawVcfRecord};
use crate::engines::align;
use crate::engines::partial_credit::{self, RefVar};
use std::collections::BTreeMap;

type TaggedPrimitive = (RefVar, bool);
type AllelePrimitives = (usize, bool, Vec<TaggedPrimitive>);

/// Decompose any multi-allelic indel record into primitive per-position
/// records and re-aggregate primitives that land at the same position.
///
/// The input is the post-stage-4 record (already through the existing
/// preprocess loop's filter / region / GT-rewrite passes) but BEFORE any
/// `split_multi_allelic` / `apply_left_shift`. Output is a sequence of
/// records that match the legacy primitive splitter's output, each suitable
/// for the existing `insert_ado_format` + `reorder_format_fields` final
/// canonicalization.
#[cfg(test)]
fn primitive_split(record: &RawVcfRecord, reference: &[u8]) -> Vec<RawVcfRecord> {
    primitive_split_with_floor(record, reference, 0)
}

/// Decompose a record while carrying the previous record's reference end
/// into the primitive normalizer. Legacy keeps this floor across consecutive
/// records; resetting it for every multi-allelic site can slide a later site
/// behind an already-emitted neighbor and make the preprocessed VCF unsorted.
#[cfg(test)]
fn primitive_split_with_floor(
    record: &RawVcfRecord,
    reference: &[u8],
    previous_end: usize,
) -> Vec<RawVcfRecord> {
    primitive_split_with_context(record, reference, previous_end, false, false)
}

pub(crate) fn primitive_split_with_context(
    record: &RawVcfRecord,
    reference: &[u8],
    previous_end: usize,
    has_following_spanning_deletion: bool,
    equal_floor_blocked: bool,
) -> Vec<RawVcfRecord> {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();

    // SNP-only multi-allelics (`C → A,G`) and trivial single-alt records
    // pass through unchanged. Same-direction multi-allelic insertions like
    // `T → TG,TTG` also fall here unless an allele is itself complex
    // (reflen > 1 && altlen > 1) — same-position primitives re-merge below.
    if !needs_primitive_split(record, &alts) {
        let mut passthrough = record.clone();
        if alts.len() == 1 {
            passthrough.primitive_identity = Some(PrimitiveIdentity {
                start: record.pos,
                end: record.pos + record.ref_allele.len().saturating_sub(1),
                alt: alts[0].to_string(),
            });
        }
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
    let mut primitives: Vec<AllelePrimitives> = Vec::with_capacity(alts.len());
    for (idx, alt) in alts.iter().enumerate() {
        primitives.push((
            idx,
            allele_has_mixed_edit(&record.ref_allele, alt),
            allele_primitives(record, alt, reference),
        ));
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

    let mut output: Vec<(RawVcfRecord, bool)> = Vec::new();
    for (allele_idx, preserve_mixed_anchor, prims) in primitives {
        let target = (allele_idx + 1) as u32;
        for (prim, was_realigned) in prims {
            if let Some(rec) = primitive_to_record(
                record,
                &prim,
                target,
                ad_index,
                gt_index,
                &format_keys,
                reference,
                preserve_mixed_anchor,
                was_realigned,
                previous_end,
                has_following_spanning_deletion,
                equal_floor_blocked,
            ) {
                let mut rec = rec;
                rec.mixed_edit_primitive = preserve_mixed_anchor;
                output.push((rec, preserve_mixed_anchor));
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
    output.sort_by(|(a, _), (b, _)| {
        a.pos.cmp(&b.pos).then_with(|| {
            b.ref_allele
                .len()
                .cmp(&a.ref_allele.len())
                .then(a.alt_allele.cmp(&b.alt_allele))
        })
    });
    let mut current_maxpos = previous_end;
    for (rec, preserve_mixed_anchor) in &mut output {
        // The legacy primitive splitter emits the substitution and indel
        // pieces of a mixed edit at their shared source anchor. Passing
        // those pieces through the generic left shifter moves the indel to
        // its trailing edge and changes both VCF identity and sibling
        // aggregation (for example CAG>A must become C>A plus CAG>C).
        if !*preserve_mixed_anchor {
            shift_primitive_record(rec, reference, current_maxpos);
        }
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
    let mut aggregated = aggregate_same_position(
        output.into_iter().map(|(record, _)| record).collect(),
        &format_keys,
        phased,
    );
    if alts
        .iter()
        .any(|alt| allele_has_mixed_edit(&record.ref_allele, alt))
    {
        aggregated.sort_by(|left, right| {
            left.pos.cmp(&right.pos).then_with(|| {
                let left_indel = left.ref_allele.len() != left.alt_allele.len();
                let right_indel = right.ref_allele.len() != right.alt_allele.len();
                left_indel
                    .cmp(&right_indel)
                    .then(left.ref_allele.len().cmp(&right.ref_allele.len()))
                    .then(left.alt_allele.cmp(&right.alt_allele))
            })
        });
    }
    aggregated
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
    // allele retains both REF and ALT bases after common prefix/suffix trim.
    // The latter includes asymmetric complex alleles such as AT→G and T→CGT;
    // checking only raw lengths >1 misses exactly those shapes.
    let has_multi_indel = alts.len() > 1 && alts.iter().any(|a| a.len() != ref_len);
    let has_complex = alts
        .iter()
        .any(|alt| allele_has_mixed_edit(&record.ref_allele, alt));
    has_multi_indel || has_complex
}

fn allele_has_mixed_edit(reference: &str, alternate: &str) -> bool {
    if reference.len() <= 1 && alternate.len() <= 1 {
        return false;
    }
    let reference = reference.as_bytes();
    let alternate = alternate.as_bytes();
    let mut prefix = 0usize;
    while prefix < reference.len()
        && prefix < alternate.len()
        && reference[prefix].eq_ignore_ascii_case(&alternate[prefix])
    {
        prefix += 1;
    }
    let mut suffix = 0usize;
    while suffix < reference.len().saturating_sub(prefix)
        && suffix < alternate.len().saturating_sub(prefix)
        && reference[reference.len() - suffix - 1]
            .eq_ignore_ascii_case(&alternate[alternate.len() - suffix - 1])
    {
        suffix += 1;
    }
    prefix + suffix < reference.len() && prefix + suffix < alternate.len()
}

/// Split an allele's REF/ALT pair into its primitive RefVars in the legacy
/// convention. Mirrors the per-allele body of
/// `VariantPrimitiveSplitter::advance` plus its post-trim realign-or-passthrough
/// guard (`src/c++/lib/variant/VariantPrimitiveSplitter.cpp:135-178`).
fn allele_primitives(record: &RawVcfRecord, alt: &str, reference: &[u8]) -> Vec<(RefVar, bool)> {
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
            .into_iter()
            .map(|primitive| (primitive, true))
            .collect()
    } else {
        // Pass the ORIGINAL untrimmed RefVar through — leftshift in stage 6
        // handles its trim-and-slide differently than this pre-check would.
        vec![(original, false)]
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
    preserve_mixed_anchor: bool,
    was_realigned: bool,
    previous_end: usize,
    has_following_spanning_deletion: bool,
    equal_floor_blocked: bool,
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
        let anchor_pos = prim.end + 1;
        let anchor = ref_slice(reference, anchor_pos, anchor_pos)?[0].to_ascii_uppercase() as char;
        let deleted: String = ref_slice(reference, prim.start, prim.end)?
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();
        if preserve_mixed_anchor
            && prim.start > 1
            && !(has_following_spanning_deletion
                && (previous_end > prim.start - 1
                    || (equal_floor_blocked && previous_end == prim.start - 1)))
        {
            let left_anchor_pos = prim.start - 1;
            let left_anchor = ref_slice(reference, left_anchor_pos, left_anchor_pos)?[0]
                .to_ascii_uppercase() as char;
            out.pos = left_anchor_pos;
            out.ref_allele = format!("{left_anchor}{deleted}");
            out.alt_allele = left_anchor.to_string();
        } else {
            out.pos = prim.start;
            out.ref_allele = format!("{deleted}{anchor}");
            out.alt_allele = anchor.to_string();
        }
    } else if reflen_i <= 0 && altlen > 0 {
        if preserve_mixed_anchor && prim.start > 1 {
            // An insertion following a substitution uses their common left
            // anchor: A>TT becomes A>T plus A>AT.
            let anchor_pos = prim.start - 1;
            let anchor =
                ref_slice(reference, anchor_pos, anchor_pos)?[0].to_ascii_uppercase() as char;
            out.pos = anchor_pos;
            out.ref_allele = anchor.to_string();
            out.alt_allele = format!("{anchor}{}", prim.alt.to_ascii_uppercase());
        } else {
            // Standalone insertion — right-anchor convention (legacy
            // trailing-edge).
            let anchor_pos = prim.start;
            let anchor =
                ref_slice(reference, anchor_pos, anchor_pos)?[0].to_ascii_uppercase() as char;
            out.pos = anchor_pos;
            out.ref_allele = anchor.to_string();
            out.alt_allele = format!("{}{anchor}", prim.alt.to_ascii_uppercase());
        }
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
            was_realigned,
        ));
    }
    out.samples = new_samples;
    out.primitive_identity = Some(PrimitiveIdentity {
        start: prim.start,
        end: prim.end,
        alt: prim.alt.clone(),
    });
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
    reverse_unphased_het: bool,
) -> String {
    let mut new_cells: Vec<String> = (0..expected_len)
        .map(|i| cells.get(i).cloned().unwrap_or_else(|| ".".to_string()))
        .collect();

    if let Some(gi) = gt_index
        && let Some(cell) = new_cells.get_mut(gi)
    {
        *cell = canonical_split_gt(cell, target);
        if reverse_unphased_het && *cell == "0/1" {
            *cell = "1/0".to_string();
        }
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

/// Canonicalise a split-allele GT to the unphased legacy form. Het calls become
/// `0/1`, hom-alt calls become `1/1`, and no-call passes through unchanged.
fn canonical_split_gt(gt: &str, target: u32) -> String {
    if gt == "." || gt == "./." || gt == ".|." || gt.is_empty() {
        return gt.to_string();
    }
    let mut tokens: Vec<i32> = Vec::new();
    for tok in gt.split(['/', '|']) {
        match tok.parse::<i32>() {
            Ok(v) => tokens.push(v),
            Err(_) => return gt.to_string(),
        }
    }
    let target_i = target as i32;
    let mut projected = tokens
        .iter()
        .map(|allele| if *allele == target_i { 1 } else { 0 })
        .collect::<Vec<_>>();
    projected.sort_unstable();
    projected
        .iter()
        .map(i32::to_string)
        .collect::<Vec<_>>()
        .join("/")
}

/// Stage 7 — aggregate primitives that landed at the same anchor position
/// back into a single multi-allelic record. ALT alleles are collected in
/// (length, alphabetical) order so the output matches legacy's canonical
/// `T → TG,TTG` shape.
fn aggregate_same_position(
    records: Vec<RawVcfRecord>,
    format_keys: &[String],
    _phased: bool,
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
        // A complex allele such as AT→G decomposes into A→G plus AT→A at
        // the same anchor. Legacy's primitive splitter keeps those as two
        // records for both phased and unphased genotypes; only compatible
        // same-kind primitives are eligible for re-aggregation. Restricting
        // this guard to phased input collapsed 37 real GIAB SNP primitives
        // back into complex INDEL rows.
        if has_substitution && has_indel {
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

        // Combined GT: reconstruct phased haplotypes position-by-position.
        // Unphased calls retain the legacy later/earlier ordering.
        if let Some(gi) = gt_index {
            let gt_cells = sorted
                .iter()
                .map(|record| {
                    record
                        .samples
                        .get(s)
                        .and_then(|sample| sample.split(':').nth(gi))
                        .unwrap_or(".")
                })
                .collect::<Vec<_>>();
            let merged_gt = if gt_cells.iter().any(|gt| gt.contains('|')) {
                merge_phased_genotypes(&gt_cells)
            } else {
                merge_unphased_genotypes(&gt_cells)
            };
            if let Some(cell) = cells.get_mut(gi) {
                *cell = merged_gt;
            }
        }
        new_samples.push(cells.join(":"));
    }
    base.samples = new_samples;
    base.primitive_identity = None;
    base
}

/// Aggregate two compatible biallelic calls emitted at the same VCF location.
///
/// Legacy's `VariantLocationAggregator` operates across input-record
/// boundaries, not only across primitives produced from one parent record.
/// It pads shorter REF alleles to the longest sibling span, extends each ALT
/// by the same suffix, and emits one het-alt record. Records that cannot fit
/// that diploid shape pass through unchanged.
pub(crate) fn aggregate_location_records(records: Vec<RawVcfRecord>) -> Vec<RawVcfRecord> {
    aggregate_location_records_inner(records, false)
}

fn aggregate_location_records_inner(
    mut records: Vec<RawVcfRecord>,
    preserve_pair_order: bool,
) -> Vec<RawVcfRecord> {
    let mixed_insertion_order = |left: &RawVcfRecord, right: &RawVcfRecord| {
        let is_insertion =
            |record: &RawVcfRecord| record.ref_allele.len() == 1 && record.alt_allele.len() > 1;
        if is_insertion(left) && is_insertion(right) {
            left.mixed_edit_primitive.cmp(&right.mixed_edit_primitive)
        } else {
            std::cmp::Ordering::Equal
        }
    };
    if records.len() > 2 {
        // Legacy feeds the incoming stream through a one-record buffer.
        // Each call is compared only with the current buffer tail: a failed
        // pair pushes a new tail that can merge with the following call, while
        // a successful pair becomes het-alt and cannot absorb a third allele.
        let has_insertion = records
            .iter()
            .any(|record| record.ref_allele.len() == 1 && record.alt_allele.len() > 1);
        let has_deletion = records
            .iter()
            .any(|record| record.ref_allele.len() > record.alt_allele.len());
        let has_snp = records
            .iter()
            .any(|record| record.ref_allele.len() == 1 && record.alt_allele.len() == 1);
        let mut snp_keys = std::collections::HashSet::new();
        let has_duplicate_snp = records.iter().any(|record| {
            record.ref_allele.len() == 1
                && record.alt_allele.len() == 1
                && !snp_keys.insert((record.ref_allele.as_str(), record.alt_allele.as_str()))
        });
        let preserve_sorted_pair_order =
            has_insertion && has_deletion && has_snp && has_duplicate_snp;
        if has_insertion && has_deletion && has_snp && !has_duplicate_snp {
            records.sort_by_key(|record| {
                usize::from(!(record.ref_allele.len() == 1 && record.alt_allele.len() == 1))
            });
        } else if !(has_insertion && has_deletion) || preserve_sorted_pair_order {
            records.sort_by(|left, right| {
                let class = |record: &RawVcfRecord| {
                    if record.ref_allele.len() == 1 && record.alt_allele.len() == 1 {
                        0
                    } else if record.ref_allele.len() > record.alt_allele.len() {
                        1
                    } else {
                        2
                    }
                };
                class(left)
                    .cmp(&class(right))
                    .then_with(|| mixed_insertion_order(left, right))
                    .then(left.ref_allele.len().cmp(&right.ref_allele.len()))
                    .then(left.alt_allele.len().cmp(&right.alt_allele.len()))
                    .then(left.alt_allele.cmp(&right.alt_allele))
            });
        }
        let mut aggregated: Vec<RawVcfRecord> = Vec::with_capacity(records.len());
        for record in records {
            let Some(previous) = aggregated.pop() else {
                aggregated.push(record);
                continue;
            };
            let merged = aggregate_location_records_inner(
                vec![previous.clone(), record.clone()],
                preserve_sorted_pair_order,
            );
            if merged.len() == 1 {
                aggregated.extend(merged);
            } else {
                aggregated.push(previous);
                aggregated.push(record);
            }
        }
        return aggregated;
    }
    if records.len() != 2
        || records[0].chrom != records[1].chrom
        || records[0].pos != records[1].pos
        || records.iter().any(|record| record.alt_allele.contains(','))
        || records[0].filter != records[1].filter
        || records[0].format != records[1].format
        || records[0].samples.len() != records[1].samples.len()
    {
        return records;
    }
    let first_is_snp = records[0].ref_allele.len() == 1 && records[0].alt_allele.len() == 1;
    let second_is_snp = records[1].ref_allele.len() == 1 && records[1].alt_allele.len() == 1;
    if first_is_snp != second_is_snp {
        return records;
    }
    let format_keys = records[0]
        .format
        .as_deref()
        .map(|format| format.split(':').collect::<Vec<_>>())
        .unwrap_or_default();
    let Some(gt_index) = format_keys.iter().position(|key| *key == "GT") else {
        return records;
    };
    let ad_index = format_keys.iter().position(|key| *key == "AD");
    let compatible_calls = records.iter().all(|record| {
        record.samples.iter().all(|sample| {
            sample
                .split(':')
                .nth(gt_index)
                .is_some_and(|gt| matches!(gt, "0/1" | "1/0"))
        })
    });
    if !compatible_calls {
        return records;
    }
    let opposite_slots = (0..records[0].samples.len()).any(|sample_index| {
        let first_gt = records[0].samples[sample_index]
            .split(':')
            .nth(gt_index)
            .unwrap_or(".");
        let second_gt = records[1].samples[sample_index]
            .split(':')
            .nth(gt_index)
            .unwrap_or(".");
        matches!((first_gt, second_gt), ("0/1", "1/0") | ("1/0", "0/1"))
    });
    let is_insertion = |record: &RawVcfRecord| record.ref_allele.len() < record.alt_allele.len();
    let is_deletion = |record: &RawVcfRecord| record.ref_allele.len() > record.alt_allele.len();
    let ordinary_insertion_with_mixed_deletion = (is_insertion(&records[0])
        && !records[0].mixed_edit_primitive
        && is_deletion(&records[1])
        && records[1].mixed_edit_primitive)
        || (is_insertion(&records[1])
            && !records[1].mixed_edit_primitive
            && is_deletion(&records[0])
            && records[0].mixed_edit_primitive);
    if opposite_slots && ordinary_insertion_with_mixed_deletion {
        return records;
    }
    let both_deletions = records
        .iter()
        .all(|record| record.ref_allele.len() > record.alt_allele.len());
    if both_deletions {
        let mixed_is_longer = records
            .iter()
            .find(|record| record.mixed_edit_primitive)
            .zip(records.iter().find(|record| !record.mixed_edit_primitive))
            .is_some_and(|(mixed, ordinary)| mixed.ref_allele.len() > ordinary.ref_allele.len());
        if opposite_slots && mixed_is_longer {
            return records;
        }
    }

    // VariantAlleleSplitter orders cross-record half-calls before padding.
    // Preserve the source spans here: sorting the final padded alleles would
    // invert overlapping deletions such as GC>G plus GCC>G.
    let ordinary_insertion_and_deletion = records.iter().all(|record| !record.mixed_edit_primitive)
        && records.iter().any(is_insertion)
        && records.iter().any(is_deletion);
    if !preserve_pair_order || ordinary_insertion_and_deletion {
        records.sort_by(|left, right| {
            mixed_insertion_order(left, right)
                .then_with(|| {
                    usize::from(!is_insertion(left)).cmp(&usize::from(!is_insertion(right)))
                })
                .then(left.ref_allele.len().cmp(&right.ref_allele.len()))
                .then(left.alt_allele.len().cmp(&right.alt_allele.len()))
                .then(left.alt_allele.cmp(&right.alt_allele))
        });
    }
    let shorter_mixed_deletion_is_first = both_deletions
        && records[0].mixed_edit_primitive
        && (records[0].ref_allele.len() < records[1].ref_allele.len()
            || (records[0].ref_allele.len() == records[1].ref_allele.len()
                && records[0].alt_allele.len() < records[1].alt_allele.len()))
        && records[0].samples.iter().all(|sample| {
            sample
                .split(':')
                .nth(gt_index)
                .is_some_and(|gt| gt == "1/0")
        })
        && records[1].samples.iter().all(|sample| {
            sample
                .split(':')
                .nth(gt_index)
                .is_some_and(|gt| gt == "0/1")
        });
    let mixed_primitive_is_first = shorter_mixed_deletion_is_first
        || (records[0].mixed_edit_primitive
            && !records[1].mixed_edit_primitive
            && records[0].ref_allele.len() == 1
            && records[0].alt_allele.len() == 1);

    let longest_ref = records
        .iter()
        .map(|record| record.ref_allele.as_str())
        .max_by_key(|reference| reference.len())
        .unwrap_or_default()
        .to_string();
    for record in &mut records {
        let Some(suffix) = longest_ref.strip_prefix(&record.ref_allele) else {
            return records;
        };
        record.ref_allele = longest_ref.clone();
        record.alt_allele.push_str(suffix);
    }
    if records[0].alt_allele == records[1].alt_allele {
        // VariantAlleleUniq.cpp keys alleles on their internal RefVar before
        // final VCF padding. Exact internal duplicates collapse to one hom-alt.
        // Distinct edits can serialize to the same padded REF/ALT spelling;
        // legacy preserves both ALT entries and their aggregator-derived
        // het-alt genotype. xcmp's VariantReader.cpp then deduplicates those
        // spellings when it reads the VCF.
        let distinct_internal_edits = records[0]
            .primitive_identity
            .as_ref()
            .zip(records[1].primitive_identity.as_ref())
            .is_some_and(|(first, second)| first != second);
        if distinct_internal_edits {
            let merged = merge_records(records, ad_index, Some(gt_index));
            return vec![merged];
        }
        let mut merged = records.remove(0);
        for sample in &mut merged.samples {
            let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
            if let Some(gt) = cells.get_mut(gt_index) {
                *gt = "1/1".to_string();
            }
            *sample = cells.join(":");
        }
        merged.primitive_identity = None;
        return vec![merged];
    }

    let mut merged = merge_records(records, ad_index, Some(gt_index));
    if mixed_primitive_is_first {
        for sample in &mut merged.samples {
            let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
            if let Some(gt) = cells.get_mut(gt_index) {
                *gt = "1/2".to_string();
            }
            *sample = cells.join(":");
        }
    }
    vec![merged]
}

fn merge_phased_genotypes(genotypes: &[&str]) -> String {
    let ploidy = genotypes
        .iter()
        .find(|genotype| genotype.contains('|'))
        .map_or(2, |genotype| genotype.split('|').count());
    let mut haplotypes = vec![0usize; ploidy];
    for (alternate_index, genotype) in genotypes.iter().enumerate() {
        let alleles = genotype.split('|').collect::<Vec<_>>();
        if alleles.len() != ploidy {
            return genotype.to_string();
        }
        for (haplotype, allele) in alleles.iter().enumerate() {
            if *allele == "1" {
                haplotypes[haplotype] = alternate_index + 1;
            }
        }
    }
    haplotypes
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join("|")
}

fn merge_unphased_genotypes(genotypes: &[&str]) -> String {
    let mut alt_calls: Vec<u32> = Vec::new();
    for (alternate_index, genotype) in genotypes.iter().enumerate() {
        for token in genotype.split('/') {
            if token == "1" {
                alt_calls.push((alternate_index + 1) as u32);
            }
        }
    }
    if alt_calls.is_empty() {
        "0/0".to_string()
    } else if alt_calls.len() == 1 {
        format!("0/{}", alt_calls[0])
    } else {
        if genotypes.len() == 2
            && genotypes
                .iter()
                .all(|genotype| matches!(*genotype, "0/1" | "1/0"))
        {
            // VariantLocationMap.cpp fills the buffered call's remaining zero
            // slot. The first call therefore decides whether the later
            // alternate lands before (`0/1` -> `2/1`) or after (`1/0` ->
            // `1/2`) the buffered alternate, independently for every sample.
            return if genotypes[0] == "1/0" {
                "1/2".to_string()
            } else {
                "2/1".to_string()
            };
        }
        // `VariantLocationAggregator` fills the last zero slot first,
        // yielding the later alternate before the earlier one.
        alt_calls.sort_unstable();
        format!("{}/{}", alt_calls[1], alt_calls[0])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

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
            mixed_edit_primitive: false,
            primitive_identity: None,
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
    fn location_aggregator_merges_adjacent_real_world_indels() {
        let shorter = make_record("chr1", 963_700, "GC", "G", "GT:AD:ADO:DP", "0/1:0,1:0:0");
        let longer = make_record("chr1", 963_700, "GCC", "G", "GT:AD:ADO:DP", "0/1:0,1:0:0");

        let out = aggregate_location_records(vec![shorter, longer]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].ref_allele, "GCC");
        assert_eq!(out[0].alt_allele, "GC,G");
        assert_eq!(out[0].samples[0], "2/1:0,1,1:0:0");
    }

    #[test]
    fn location_aggregator_preserves_incompatible_calls() {
        let first = make_record("chr1", 100, "A", "G", "GT", "1/1");
        let second = make_record("chr1", 100, "A", "T", "GT", "0/1");

        let out = aggregate_location_records(vec![first, second]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "G");
        assert_eq!(out[1].alt_allele, "T");
    }

    #[test]
    fn location_aggregator_slides_after_an_incompatible_call() {
        let homalt = make_record("chr1", 100, "A", "AT", "GT", "1/1");
        let first_het = make_record("chr1", 100, "A", "ATT", "GT", "0/1");
        let second_het = make_record("chr1", 100, "A", "ATTT", "GT", "0/1");

        let out = aggregate_location_records(vec![second_het, homalt, first_het]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "AT");
        assert_eq!(out[1].alt_allele, "ATT,ATTT");
        assert_eq!(out[1].samples[0], "2/1");
    }

    #[test]
    fn location_aggregator_keeps_different_filters_separate() {
        let mut first = make_record("chr1", 100, "A", "AT", "GT", "0/1");
        first.filter = "OverlapConflict;Silver".to_string();
        let mut second = make_record("chr1", 100, "A", "ATT", "GT", "0/1");
        second.filter = "OverlapConflict".to_string();

        let out = aggregate_location_records(vec![first, second]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "AT");
        assert_eq!(out[1].alt_allele, "ATT");
    }

    #[test]
    fn location_aggregator_stops_after_forming_hetalt() {
        let shortest = make_record("chr1", 100, "A", "AT", "GT", "0/1");
        let middle = make_record("chr1", 100, "A", "ATT", "GT", "0/1");
        let longest = make_record("chr1", 100, "A", "ATTT", "GT", "0/1");

        let out = aggregate_location_records(vec![longest, shortest, middle]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "AT,ATT");
        assert_eq!(out[1].alt_allele, "ATTT");
    }

    #[test]
    fn location_aggregator_preserves_three_record_arrival_order() {
        let insertion = make_record("chr1", 54_932_993, "G", "GTA", "GT", "0/1");
        let shorter_deletion = make_record("chr1", 54_932_993, "GCAT", "G", "GT", "0/1");
        let longer_deletion = make_record("chr1", 54_932_993, "GCATTT", "G", "GT", "0/1");

        let out = aggregate_location_records(vec![insertion, shorter_deletion, longer_deletion]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].ref_allele, "GCAT");
        assert_eq!(out[0].alt_allele, "GTACAT,G");
        assert_eq!(out[0].samples[0], "2/1");
        assert_eq!(out[1].ref_allele, "GCATTT");
        assert_eq!(out[1].alt_allele, "G");
        assert_eq!(out[1].samples[0], "0/1");
    }

    #[test]
    fn location_aggregator_orders_all_ordinary_mixed_types() {
        let ordinary_snp = make_record("chr2", 145_533_317, "A", "G", "GT", "0/1");
        let mut mixed_snp = make_record("chr2", 145_533_317, "A", "G", "GT", "1/0");
        mixed_snp.mixed_edit_primitive = true;
        let mut insertion = make_record("chr2", 145_533_317, "A", "ATGTGTGTG", "GT", "1/0");
        insertion.mixed_edit_primitive = true;
        let deletion = make_record("chr2", 145_533_317, "ATG", "A", "GT", "0/1");

        let out = aggregate_location_records(vec![ordinary_snp, mixed_snp, insertion, deletion]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "G");
        assert_eq!(out[1].ref_allele, "ATG");
        assert_eq!(out[1].alt_allele, "A,ATGTGTGTGTG");
        assert_eq!(out[1].samples[0], "2/1");
    }

    #[test]
    fn location_aggregator_keeps_opposite_slot_deletions_separate() {
        let ordinary_snp = make_record("chr2", 60_749_238, "A", "G", "GT", "0/1");
        let mut mixed_snp = make_record("chr2", 60_749_238, "A", "G", "GT", "1/0");
        mixed_snp.mixed_edit_primitive = true;
        let ordinary_deletion = make_record("chr2", 60_749_238, "AAT", "A", "GT", "0/1");
        let mut mixed_deletion = make_record("chr2", 60_749_238, "AATT", "A", "GT", "1/0");
        mixed_deletion.mixed_edit_primitive = true;

        let out = aggregate_location_records(vec![
            ordinary_snp,
            ordinary_deletion,
            mixed_snp,
            mixed_deletion,
        ]);

        assert_eq!(out.len(), 3);
        assert_eq!(out[0].ref_allele, "A");
        assert_eq!(out[0].alt_allele, "G");
        assert_eq!(out[0].samples[0], "1/1");
        assert_eq!(out[1].ref_allele, "AAT");
        assert_eq!(out[1].alt_allele, "A");
        assert_eq!(out[1].samples[0], "0/1");
        assert_eq!(out[2].ref_allele, "AATT");
        assert_eq!(out[2].alt_allele, "A");
        assert_eq!(out[2].samples[0], "1/0");
    }

    #[test]
    fn location_aggregator_merges_mixed_and_ordinary_opposite_slot_deletions() {
        let mut mixed = make_record("chr4", 6_943_863, "AACTTTTA", "A", "GT", "1/0");
        mixed.mixed_edit_primitive = true;
        let mut ordinary = make_record("chr4", 6_943_863, "AACT", "AA", "GT", "0/1");
        ordinary.mixed_edit_primitive = true;

        let out = aggregate_location_records(vec![ordinary, mixed]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].ref_allele, "AACTTTTA");
        assert_eq!(out[0].alt_allele, "AATTTA,A");
        assert_eq!(out[0].samples[0], "2/1");
    }

    #[test]
    fn location_aggregator_keeps_ordinary_insertion_and_mixed_deletion_separate() {
        let insertion = make_record("chr11", 95_814_076, "T", "TA", "GT", "0/1");
        let mut deletion = make_record("chr11", 95_814_076, "TCT", "T", "GT", "1/0");
        deletion.mixed_edit_primitive = true;

        let out = aggregate_location_records(vec![insertion, deletion]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].ref_allele, "T");
        assert_eq!(out[0].alt_allele, "TA");
        assert_eq!(out[1].ref_allele, "TCT");
        assert_eq!(out[1].alt_allele, "T");
    }

    #[test]
    fn location_aggregator_orders_ordinary_insertion_before_deletion_with_duplicate_snp() {
        let snp = make_record("chr4", 90_057_486, "T", "C", "GT", "1/0");
        let mut mixed_snp = snp.clone();
        mixed_snp.mixed_edit_primitive = true;
        let insertion = make_record("chr4", 90_057_486, "T", "TGGC", "GT", "1/0");
        let mut mixed_insertion = make_record("chr4", 90_057_486, "T", "TC", "GT", "1/0");
        mixed_insertion.mixed_edit_primitive = true;
        let deletion = make_record("chr4", 90_057_486, "TTA", "T", "GT", "1/0");

        let out =
            aggregate_location_records(vec![snp, mixed_snp, insertion, mixed_insertion, deletion]);

        assert_eq!(out.len(), 3);
        assert_eq!(out[0].alt_allele, "C");
        assert_eq!(out[1].ref_allele, "TTA");
        assert_eq!(out[1].alt_allele, "TGGCTA,T");
        assert_eq!(out[2].alt_allele, "TC");
    }

    #[test]
    fn location_aggregator_merges_mixed_shorter_deletion_first() {
        let mut mixed_shorter = make_record("chr10", 34_496_891, "AAAA", "A", "GT", "1/0");
        mixed_shorter.mixed_edit_primitive = true;
        let ordinary_longer = make_record("chr10", 34_496_891, "AAAATAGTATAC", "A", "GT", "0/1");

        let out = aggregate_location_records(vec![mixed_shorter, ordinary_longer]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].ref_allele, "AAAATAGTATAC");
        assert_eq!(out[0].alt_allele, "ATAGTATAC,A");
        assert_eq!(out[0].samples[0], "1/2");
    }

    #[test]
    fn location_aggregator_preserves_equal_span_deletion_slot_order() {
        let mut shorter_alt = make_record("chr21", 14_970_591, "CTCAACTAG", "C", "GT", "1/0");
        shorter_alt.mixed_edit_primitive = true;
        let mut longer_alt = make_record("chr21", 14_970_591, "CTCAACTAG", "CT", "GT", "0/1");
        longer_alt.mixed_edit_primitive = true;

        let out = aggregate_location_records(vec![shorter_alt, longer_alt]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "C,CT");
        assert_eq!(out[0].samples[0], "1/2");
    }

    #[test]
    fn location_aggregator_preserves_mixed_snp_and_indel_calls() {
        let first = make_record("chr1", 100, "A", "G", "GT", "0/1");
        let second = make_record("chr1", 100, "A", "AT", "GT", "0/1");

        let out = aggregate_location_records(vec![first, second]);

        assert_eq!(out.len(), 2);
    }

    #[test]
    fn location_aggregator_keeps_ordinary_insertions_in_length_order() {
        let longer_alt = make_record("chr10", 18_102_686, "A", "ATAT", "GT:AD", "0/1:0,1");
        let shorter_alt = make_record("chr10", 18_102_686, "A", "ATT", "GT:AD", "0/1:0,1");

        let out = aggregate_location_records(vec![longer_alt, shorter_alt]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "ATT,ATAT");
        assert_eq!(out[0].samples[0], "2/1:0,1,1");
    }

    #[test]
    fn location_aggregator_places_mixed_edit_insertion_after_direct_insertion() {
        let direct = make_record("chr2", 7_654_671, "A", "ATTTGGT", "GT:AD", "0/1:0,1");
        let mut mixed = make_record("chr2", 7_654_671, "A", "ATGGT", "GT:AD", "0/1:0,1");
        mixed.mixed_edit_primitive = true;
        let mut sibling_snp = make_record("chr2", 7_654_671, "A", "T", "GT:AD", "1/1:0,1");
        sibling_snp.mixed_edit_primitive = true;

        let out = aggregate_location_records(vec![mixed, direct, sibling_snp]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "T");
        assert_eq!(out[1].alt_allele, "ATTTGGT,ATGGT");
        assert_eq!(out[1].samples[0], "2/1:0,1,1");
    }

    #[test]
    fn location_aggregator_keeps_mixed_edit_snps_in_allele_order() {
        let direct = make_record("chr1", 35_412_089, "C", "A", "GT", "0/1");
        let mut mixed = make_record("chr1", 35_412_089, "C", "G", "GT", "0/1");
        mixed.mixed_edit_primitive = true;

        let out = aggregate_location_records(vec![mixed, direct]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "A,G");
        assert_eq!(out[0].samples[0], "2/1");
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
        // Reference bytes at 17809630..=17809631 are
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
    fn split_genotypes_are_unphased_before_realignment_orientation() {
        assert_eq!(canonical_split_gt("1|2", 1), "0/1");
        assert_eq!(canonical_split_gt("1|2", 2), "0/1");
        assert_eq!(canonical_split_gt("2|1", 1), "0/1");
        assert_eq!(canonical_split_gt("2|1", 2), "0/1");
    }

    #[test]
    fn passthrough_insertions_are_unphased_after_alt_reordering() {
        let reference = windowed_ref(11_101_380, b"NNNNNNTNNNN");
        for input_gt in ["1|2", "2|1"] {
            let sample = format!("{input_gt}:2,36,137");
            let record = make_record("chr21", 11_101_386, "T", "TTG,TG", "GT:AD", &sample);

            let output = primitive_split(&record, &reference);

            assert_eq!(output.len(), 1);
            assert_eq!(output[0].alt_allele, "TG,TTG");
            assert_eq!(output[0].samples[0].split(':').next(), Some("2/1"));
        }
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
    fn unphased_complex_substitution_and_deletion_remain_separate() {
        // The following reference base is A, so the one-base deletion
        // left-shifts back onto the source anchor exactly as at the GIAB
        // AT→G loci that exposed the re-aggregation bug.
        let reference = windowed_ref(100, b"NNNNNATANNN");
        let rec = make_record("chr1", 105, "AT", "G", "GT:AD", "0/1:5,9");
        let out = primitive_split(&rec, &reference);

        assert_eq!(out.len(), 2);
        assert!(out.iter().any(|record| {
            record.pos == 105
                && record.ref_allele == "A"
                && record.alt_allele == "G"
                && record.samples[0] == "1/0:5,9"
        }));
        assert!(out.iter().any(|record| {
            record.pos == 105
                && record.ref_allele == "AT"
                && record.alt_allele == "A"
                && record.samples[0] == "1/0:5,9"
        }));
    }

    #[test]
    fn mixed_deletion_respects_previous_record_floor() {
        let reference = windowed_ref(34_496_889, b"TAAAAATAGTATAC");
        let rec = make_record("chr10", 34_496_890, "AAAA", "C", "GT:AD", "0|1:.,.");

        let out = primitive_split_with_context(&rec, &reference, 34_496_892, true, false);

        assert_eq!(out.len(), 2);
        assert!(out.iter().any(|record| {
            record.pos == 34_496_890 && record.ref_allele == "A" && record.alt_allele == "C"
        }));
        assert!(out.iter().any(|record| {
            record.pos == 34_496_891 && record.ref_allele == "AAAA" && record.alt_allele == "A"
        }));
    }

    #[test]
    fn mixed_deletion_uses_source_anchor_without_a_following_spanning_deletion() {
        let reference = windowed_ref(42_913_280, b"ATATATACACACACGTATATA");
        let rec = make_record("chr1", 42_913_287, "CACACAC", "T", "GT:AD", "0|1:.,.");

        let out = primitive_split_with_context(&rec, &reference, 42_913_291, false, false);

        assert!(out.iter().any(|record| {
            record.pos == 42_913_287 && record.ref_allele == "CACACAC" && record.alt_allele == "C"
        }));
    }

    #[test]
    fn adjacent_preceding_deletion_blocks_mixed_deletion_at_equal_floor() {
        let reference = windowed_ref(69_737_852, b"TCGGGCGGATCACA");
        let rec = make_record("chr14", 69_737_857, "CG", "A", "GT:AD", "1|0:.,.");

        let out = primitive_split_with_context(&rec, &reference, 69_737_857, true, true);

        assert!(out.iter().any(|record| record.pos == 69_737_858));
        assert!(!out.iter().any(|record| {
            record.pos == 69_737_857 && record.ref_allele == "CG" && record.alt_allele == "C"
        }));
    }

    #[test]
    fn mixed_deletion_shares_source_anchor_when_trailing_base_repeats() {
        let reference = windowed_ref(100, b"NNNNNTAGTAGGNN");
        let rec = make_record("chr12", 105, "TAGTAG", "C", "GT:AD", "1|0:5,9");

        let out = primitive_split_with_floor(&rec, &reference, 109);

        assert!(out.iter().any(|record| {
            record.pos == 105 && record.ref_allele == "T" && record.alt_allele == "C"
        }));
        assert!(out.iter().any(|record| {
            record.pos == 105 && record.ref_allele == "TAGTAG" && record.alt_allele == "T"
        }));
    }

    #[test]
    fn unphased_complex_substitution_and_insertion_share_the_source_anchor() {
        let reference = windowed_ref(100, b"NNNNNANNNNN");
        let rec = make_record("chr1", 105, "A", "TT", "GT:AD", "0/1:5,9");
        let out = primitive_split(&rec, &reference);

        assert_eq!(out.len(), 2);
        assert!(out.iter().any(|record| {
            record.pos == 105
                && record.ref_allele == "A"
                && record.alt_allele == "T"
                && record.samples[0] == "1/0:5,9"
        }));
        assert!(out.iter().any(|record| {
            record.pos == 105
                && record.ref_allele == "A"
                && record.alt_allele == "AT"
                && record.samples[0] == "1/0:5,9"
        }));
    }

    #[test]
    fn normative_non_realignable_mixed_allele_uses_canonical_unphased_het() {
        let reference = windowed_ref(100, b"NNNNNCANNNNN");
        for input_gt in ["0/1", "1|0"] {
            let sample = format!("{input_gt}:5,9");
            let record = make_record("chr1", 105, "CA", "CG", "GT:AD", &sample);

            let output = primitive_split(&record, &reference);

            assert_eq!(output.len(), 1);
            assert_eq!(output[0].ref_allele, "CA");
            assert_eq!(output[0].alt_allele, "CG");
            assert_eq!(output[0].samples[0], "0/1:5,9");
        }
    }

    #[test]
    fn legacy_only_realigned_mixed_allele_uses_reversed_unphased_het() {
        let reference = windowed_ref(100, b"NNNNNANNNNN");
        for input_gt in ["0/1", "0|1"] {
            let sample = format!("{input_gt}:5,9");
            let record = make_record("chr1", 105, "A", "TT", "GT:AD", &sample);

            let output = primitive_split(&record, &reference);

            assert_eq!(output.len(), 2);
            assert!(
                output
                    .iter()
                    .all(|primitive| primitive.samples[0] == "1/0:5,9")
            );
        }
    }

    #[test]
    fn duplicate_snp_primitives_fill_both_genotype_slots() {
        let first = make_record("chr1", 105, "C", "A", "GT:AD", "1/0:5,9");
        let second = make_record("chr1", 105, "C", "A", "GT:AD", "1/0:5,9");

        let out = aggregate_location_records(vec![first, second]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "A");
        assert_eq!(out[0].samples[0], "1/1:5,9");
    }

    #[test]
    fn legacy_only_distinct_internal_edits_keep_duplicate_padded_alleles() {
        let mut direct = make_record("chr1", 105, "A", "AT", "GT", "0/1");
        direct.primitive_identity = Some(PrimitiveIdentity {
            start: 105,
            end: 105,
            alt: "AT".to_string(),
        });
        let mut realigned = make_record("chr1", 105, "A", "AT", "GT", "1/0");
        realigned.mixed_edit_primitive = true;
        realigned.primitive_identity = Some(PrimitiveIdentity {
            start: 106,
            end: 105,
            alt: "T".to_string(),
        });

        let out = aggregate_location_records(vec![direct, realigned]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "AT,AT");
        assert_eq!(out[0].samples[0], "2/1");
    }

    #[test]
    fn normative_exact_internal_edits_collapse_duplicate_padded_alleles() {
        let identity = PrimitiveIdentity {
            start: 105,
            end: 105,
            alt: "AT".to_string(),
        };
        let mut first = make_record("chr1", 105, "A", "AT", "GT", "0/1");
        first.primitive_identity = Some(identity.clone());
        let mut second = make_record("chr1", 105, "A", "AT", "GT", "1/0");
        second.primitive_identity = Some(identity);

        let out = aggregate_location_records(vec![first, second]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "AT");
        assert_eq!(out[0].samples[0], "1/1");
    }

    #[test]
    fn legacy_only_distinct_internal_edits_respect_each_sample_buffered_slot() {
        let mut direct = make_record("chr1", 105, "A", "AT", "GT", "0/1");
        direct.samples.push("1/0".to_string());
        direct.primitive_identity = Some(PrimitiveIdentity {
            start: 105,
            end: 105,
            alt: "AT".to_string(),
        });
        let mut realigned = make_record("chr1", 105, "A", "AT", "GT", "1/0");
        realigned.samples.push("0/1".to_string());
        realigned.mixed_edit_primitive = true;
        realigned.primitive_identity = Some(PrimitiveIdentity {
            start: 106,
            end: 105,
            alt: "T".to_string(),
        });

        let out = aggregate_location_records(vec![direct, realigned]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "AT,AT");
        assert_eq!(out[0].samples, ["2/1", "1/2"]);
    }

    #[test]
    fn normative_exact_internal_edits_fold_every_compatible_sample() {
        let identity = PrimitiveIdentity {
            start: 105,
            end: 105,
            alt: "AT".to_string(),
        };
        let mut first = make_record("chr1", 105, "A", "AT", "GT", "0/1");
        first.samples.push("1/0".to_string());
        first.primitive_identity = Some(identity.clone());
        let mut second = make_record("chr1", 105, "A", "AT", "GT", "1/0");
        second.samples.push("0/1".to_string());
        second.primitive_identity = Some(identity);

        let out = aggregate_location_records(vec![first, second]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "AT");
        assert_eq!(out[0].samples, ["1/1", "1/1"]);
    }

    #[test]
    fn mixed_and_direct_snp_use_location_aggregator_hetalt_order() {
        let mut mixed = make_record("chr1", 105, "C", "A", "GT:AD", "1/0:5,9");
        mixed.mixed_edit_primitive = true;
        let direct = make_record("chr1", 105, "C", "G", "GT:AD", "0/1:5,9");

        let out = aggregate_location_records(vec![mixed, direct]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "A,G");
        assert_eq!(out[0].samples[0], "1/2:5,9,9");
    }

    #[test]
    fn mixed_snp_in_second_alt_slot_keeps_location_aggregator_order() {
        let direct = make_record("chr1", 105, "G", "A", "GT:AD", "0/1:5,9");
        let mut mixed = make_record("chr1", 105, "G", "T", "GT:AD", "1/0:5,9");
        mixed.mixed_edit_primitive = true;

        let out = aggregate_location_records(vec![mixed, direct]);

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].alt_allele, "A,T");
        assert_eq!(out[0].samples[0], "2/1:5,9,9");
    }

    #[test]
    fn mixed_four_record_location_aggregates_each_variant_class() {
        let first_snp = make_record("chr2", 105, "A", "G", "GT:AD", "1/0:5,9");
        let first_insertion = make_record("chr2", 105, "A", "AG", "GT:AD", "1/0:5,9");
        let second_snp = make_record("chr2", 105, "A", "G", "GT:AD", "1/0:5,9");
        let second_insertion = make_record("chr2", 105, "A", "AGG", "GT:AD", "1/0:5,9");

        let out = aggregate_location_records(vec![
            first_snp,
            first_insertion,
            second_snp,
            second_insertion,
        ]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].alt_allele, "G");
        assert_eq!(out[0].samples[0], "1/1:5,9");
        assert_eq!(out[1].alt_allele, "AG,AGG");
        assert_eq!(out[1].samples[0], "1/2:5,9,9");
    }

    #[test]
    fn location_aggregator_partitions_snp_before_inverse_indel_pair() {
        let insertion = make_record("chr2", 105, "C", "CT", "GT", "0/1");
        let snp = make_record("chr2", 105, "C", "T", "GT", "0/1");
        let deletion = make_record("chr2", 105, "CT", "C", "GT", "0/1");

        let out = aggregate_location_records(vec![insertion, snp, deletion]);

        assert_eq!(out.len(), 2);
        assert_eq!(out[0].ref_allele, "C");
        assert_eq!(out[0].alt_allele, "T");
        assert_eq!(out[1].ref_allele, "CT");
        assert_eq!(out[1].alt_allele, "CTT,C");
        assert_eq!(out[1].samples[0], "2/1");
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
    fn canonical_split_gt_uses_ref_then_alt_for_unphased_diploid_het() {
        assert_eq!(canonical_split_gt("1/2", 1), "0/1");
        assert_eq!(canonical_split_gt("1/2", 2), "0/1");
        // Hom-alt of target stays homozygous.
        assert_eq!(canonical_split_gt("1/1", 1), "1/1");
        // Phased source calls are also emitted unphased.
        assert_eq!(canonical_split_gt("1|2", 1), "0/1");
        // No-call passes through.
        assert_eq!(canonical_split_gt("./.", 1), "./.");
    }

    proptest! {
        #[test]
        fn split_genotype_canonicalization_is_idempotent_across_ploidies_and_missing_calls(
            alleles in proptest::collection::vec(proptest::option::of(0u32..4), 1..=6),
            phased in any::<bool>(),
        ) {
            let separator = if phased { "|" } else { "/" };
            let gt = alleles
                .iter()
                .map(|allele| allele.map_or_else(|| ".".to_string(), |value| value.to_string()))
                .collect::<Vec<_>>()
                .join(separator);
            let once = canonical_split_gt(&gt, 1);
            prop_assert_eq!(canonical_split_gt(&once, 1), once);
        }

        #[test]
        fn unphased_projection_preserves_ploidy_and_uses_ref_then_alt(
            alleles in proptest::collection::vec(0u32..4, 2..=2),
            target in 1u32..4,
        ) {
            let gt = alleles.iter().map(u32::to_string).collect::<Vec<_>>().join("/");
            let projected = canonical_split_gt(&gt, target);
            let observed = projected
                .split('/')
                .map(str::parse::<u32>)
                .collect::<Result<Vec<_>, _>>()
                .expect("projected genotype is numeric");
            let mut expected = alleles
                .iter()
                .map(|allele| u32::from(*allele == target))
                .collect::<Vec<_>>();
            expected.sort_unstable();
            prop_assert_eq!(observed, expected);
        }

        #[test]
        fn phased_projection_becomes_unphased_ref_then_alt(
            alleles in proptest::collection::vec(0u32..4, 2..=2),
            target in 1u32..4,
        ) {
            let gt = alleles.iter().map(u32::to_string).collect::<Vec<_>>().join("|");
            let projected = canonical_split_gt(&gt, target);
            let observed = projected.split('/').collect::<Vec<_>>();
            prop_assert_eq!(observed.len(), alleles.len());
            let mut expected = alleles
                .iter()
                .map(|allele| u32::from(*allele == target))
                .collect::<Vec<_>>();
            expected.sort_unstable();
            prop_assert_eq!(observed, expected.iter().map(u32::to_string).collect::<Vec<_>>());
        }
    }
}
