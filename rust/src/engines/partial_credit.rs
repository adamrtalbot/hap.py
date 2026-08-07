//! Partial-credit left-shift and primitive decomposition pass.
//!
//! Port of the C++ `src/c++/lib/align/RefVar.cpp` helpers used by legacy
//! hap.py's `preprocess` binary (the second stage of `pre.py` — see
//! `src/python/Haplo/partialcredit.py`). The legacy pipeline runs this pass
//! on every preprocessed VCF to:
//!
//! 1. Trim common prefix + suffix bases from REF/ALT.
//! 2. Left-shift ambiguous indels as far as possible within a per-sample
//!    window (`leftshift_limit`, 1 kbp in hap.py).
//! 3. Decompose multi-allelic records into primitive SNP + indel pieces.
//!
//! This module is intentionally self-contained: inputs are plain strings and
//! a reference sequence slice so it can be unit-tested without going through
//! the VCF reader.

/// A single REF/ALT variant in 1-based, inclusive `[start, end]` coordinates.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct RefVar {
    pub start: usize,
    pub end: usize,
    pub alt: String,
}

/// Case-insensitive ASCII equality treating `N` as a wildcard.
fn base_equal(a: u8, b: u8) -> bool {
    let au = a.to_ascii_uppercase();
    let bu = b.to_ascii_uppercase();
    au == bu || au == b'N' || bu == b'N'
}

fn reflen(rv: &RefVar) -> i64 {
    rv.end as i64 - rv.start as i64 + 1
}

/// Fetch an inclusive slice of the reference in 1-based coordinates.
fn query(reference: &[u8], start: usize, end: usize) -> &[u8] {
    if start == 0 || end < start || end > reference.len() {
        return &[];
    }
    &reference[start - 1..end]
}

/// Trim common prefix bases between REF and ALT. `ref_padding=true` matches
/// legacy semantics where at least one base is left on each side (VCF anchor).
pub(crate) fn trim_left(reference: &[u8], rv: &mut RefVar, ref_padding: bool) {
    let ref_min = if ref_padding { 1 } else { 0 };
    // Legacy fetches the REF slice once and indexes into it as start moves.
    let ref_slice = query(reference, rv.start, rv.end).to_vec();
    let alt = rv.alt.as_bytes().to_vec();
    let mut rel = 0usize;
    while ref_slice.len().saturating_sub(rel) > ref_min
        && alt.len().saturating_sub(rel) > ref_min
        && rel < ref_slice.len()
        && rel < alt.len()
        && base_equal(ref_slice[rel], alt[rel])
    {
        rel += 1;
        rv.start += 1;
    }
    if rel > 0 {
        rv.alt = rv.alt[rel..].to_string();
    }
}

/// Trim common suffix bases between REF and ALT. Mirrors legacy `trimRight`.
pub(crate) fn trim_right(reference: &[u8], rv: &mut RefVar, ref_padding: bool) {
    let min_len: i64 = if ref_padding { 1 } else { 0 };
    let mut reflen = reflen(rv);
    let mut altlen = rv.alt.len() as i64;
    if reflen <= min_len || altlen <= min_len {
        return;
    }
    let ref_slice = query(reference, rv.start, rv.end);
    if ref_slice.is_empty() {
        return;
    }
    let alt = rv.alt.as_bytes();
    while reflen > min_len
        && altlen > min_len
        && base_equal(ref_slice[(reflen - 1) as usize], alt[(altlen - 1) as usize])
    {
        reflen -= 1;
        altlen -= 1;
    }
    rv.end = rv.start + reflen as usize - 1;
    rv.alt = if altlen > 0 {
        rv.alt[..altlen as usize].to_string()
    } else {
        String::new()
    };
}

/// Left-shift an ambiguous indel within `[pos_min, rv.start]`. Port of
/// `RefVar.cpp:leftShift`. `pos_min` is 1-based and inclusive.
pub(crate) fn left_shift(reference: &[u8], rv: &mut RefVar, mut pos_min: usize, ref_padding: bool) {
    if pos_min == 0 {
        pos_min = 1;
    }

    trim_left(reference, rv, ref_padding);
    trim_right(reference, rv, ref_padding);

    let mut reflen_i = reflen(rv);

    // Reserve an extra anchor position to the left for pure insertions.
    if ref_padding && reflen_i <= 1 && rv.alt.len() as i64 > reflen_i {
        let mut pad_left = true;
        if reflen_i == 1 {
            let ref_fb = query(reference, rv.start, rv.start);
            if let (Some(first_alt), Some(first_ref)) = (rv.alt.as_bytes().first(), ref_fb.first())
            {
                pad_left = base_equal(*first_alt, *first_ref);
            }
        }
        if pad_left {
            pos_min = pos_min.saturating_add(1);
        }
    }

    if reflen_i < 0 && rv.alt.is_empty() {
        return;
    }

    if reflen_i >= 0 && reflen_i == rv.alt.len() as i64 {
        let ref_al = query(reference, rv.start, rv.end);
        if ref_al == rv.alt.as_bytes() {
            return;
        }
    }

    // Maintain a cached window [rstart, rend] of the reference to avoid
    // re-querying every iteration (mirrors the C++ buffering behaviour).
    let mut cached_start: Option<usize> = None;
    let mut cached_end: Option<usize> = None;
    let mut cached: Vec<u8> = Vec::new();

    let mut done = false;
    while !done {
        done = true;
        reflen_i = reflen(rv);

        let need_new = match (cached_start, cached_end) {
            (Some(cs), Some(ce)) => rv.start <= cs || rv.end > ce,
            _ => true,
        };
        if need_new {
            let rstart = rv.start.saturating_sub(20).max(1);
            let rend = if reflen_i <= 0 { rv.start } else { rv.end };
            cached = query(reference, rstart, rend).to_vec();
            cached_start = Some(rstart);
            cached_end = Some(rend);
        }

        if rv.start <= pos_min {
            break;
        }

        let rs = cached_start.expect("cache populated above");
        let rel_start = rv.start as i64 - rs as i64;

        // Don't shift onto 'N', and keep at least one left-anchor base.
        if rel_start < 1
            || cached.is_empty()
            || (cached.len() as i64) < rel_start + reflen_i
            || cached[(rel_start - 1) as usize].eq_ignore_ascii_case(&b'N')
        {
            break;
        }

        let alt_bytes = rv.alt.as_bytes();
        // Slide-enable check: the deleted segment's last base must match the
        // alt's last base for a left-shift step to preserve semantics.
        // `base_equal` treats N as a wildcard for trim_left/trim_right (where
        // legacy permits soft-masked / unknown bases to cancel an anchor),
        // but the slide loop must use STRICT equality: an N at the right
        // edge of the deletion has unknown identity, so legacy is conservative
        // and refuses to slide through it. chr21:44035893 (`CGNNNNNNNN→C`,
        // chr21 case) exposes this — with the wildcard rule rust slid 7 bases
        // left into the surrounding `cacacac` repeat at 44035886 because
        // every right-edge N "matched" the prepended anchor; legacy stayed
        // put because the deleted run ends in N.
        let last_ref = cached[(rel_start + reflen_i - 1) as usize].to_ascii_uppercase();
        let last_alt = alt_bytes.last().copied().unwrap_or(0).to_ascii_uppercase();
        let slide_match =
            !alt_bytes.is_empty() && last_ref != b'N' && last_alt != b'N' && last_ref == last_alt;
        if reflen_i > 0 && slide_match {
            // Right-trim one base to enable the slide.
            let new_alt_len = rv.alt.len() - 1;
            rv.end -= 1;
            rv.alt.truncate(new_alt_len);
            done = false;
            reflen_i -= 1;
        }

        if reflen_i == 0 || rv.alt.is_empty() {
            // Soft-masked lowercase bases get normalised to uppercase so the
            // rebuilt ALT matches legacy's canonical output.
            let prepend = cached[(rel_start - 1) as usize].to_ascii_uppercase() as char;
            rv.alt = format!("{}{}", prepend, rv.alt);
            rv.start -= 1;
            done = false;
        }
    }

    trim_left(reference, rv, ref_padding);
    trim_right(reference, rv, ref_padding);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rv(start: usize, end: usize, alt: &str) -> RefVar {
        RefVar {
            start,
            end,
            alt: alt.to_string(),
        }
    }

    #[test]
    fn trim_left_strips_common_prefix() {
        // REF=GCT ALT=GCG at pos 4..6 over reference "AAGCGCT"
        let reference = b"AAGCGCT";
        let mut v = rv(5, 7, "GCG");
        // ref is GCT, alt is GCG — common prefix "GC", one diff at end.
        trim_left(reference, &mut v, true);
        assert_eq!(v.start, 7);
        assert_eq!(v.end, 7);
        assert_eq!(v.alt, "G");
    }

    #[test]
    fn trim_right_strips_common_suffix() {
        // REF=CAT ALT=CGT over reference "CAT" at pos 1..3
        let reference = b"CAT";
        let mut v = rv(1, 3, "CGT");
        trim_right(reference, &mut v, true);
        assert_eq!(v.start, 1);
        assert_eq!(v.end, 2);
        assert_eq!(v.alt, "CG");
    }

    #[test]
    fn left_shift_slides_indel_through_homopolymer() {
        // reference: AAAAAC, an A insertion at the last A position should
        // left-shift through the homopolymer until it reaches the left anchor
        // (pos_min=1 gets bumped to 2 for pure insertions to preserve the
        // VCF left-anchor base).
        //   pos:     1 2 3 4 5 6
        //   ref:     A A A A A C
        let reference = b"AAAAAC";
        let mut v = rv(5, 5, "AA");
        left_shift(reference, &mut v, 1, true);
        assert_eq!(v.start, 2);
        assert_eq!(v.end, 2);
        assert_eq!(v.alt, "AA");
    }

    #[test]
    fn left_shift_respects_pos_min_boundary() {
        let reference = b"AAAAAC";
        let mut v = rv(5, 5, "AA");
        // Boundary at pos=3 — bumped to pos=4 for pure insertions, so the
        // slide stops at position 4.
        left_shift(reference, &mut v, 3, true);
        assert_eq!(v.start, 4);
        assert_eq!(v.end, 4);
        assert_eq!(v.alt, "AA");
    }

    #[test]
    fn left_shift_stops_at_n_base() {
        // reference: ANAAAC at pos 1..6 — sliding into the N is blocked.
        let reference = b"ANAAAC";
        let mut v = rv(5, 5, "AA");
        left_shift(reference, &mut v, 1, true);
        // We can slide from 5 to 3 (past "AAA") but not onto the N at pos 2.
        assert!(v.start >= 3);
    }

    /// Pin Class F (chr21:44035893): a deletion whose deleted segment
    /// ends in N must NOT slide left through the surrounding repeat.
    /// Pre-fix, `base_equal` treated N as a wildcard inside the slide
    /// loop, letting any prepended anchor "match" the right-edge N and
    /// drag the deletion 7 bases left into the `cacacac` upstream window.
    /// Reference window (1-based): pos 1..20 = `cacacacaccacacgnnnnn`.
    /// Variant: pos=15 ref="CGNNNNNN" alt="C" (an 8-base deletion whose
    /// last base is N). Expected: no shift; rv stays at start=15.
    #[test]
    fn left_shift_refuses_to_slide_when_deletion_ends_in_n() {
        let reference = b"cacacacaccacacgnnnnn";
        // pos=15 -> 1-based byte index 14 = 'g'? Wait, the reference layout:
        //   pos 1..14 = "cacacacaccacac" (14 bp)
        //   pos 15..20 = "gnnnnn" (6 bp)
        // The variant `CGNNNNNN→C` at pos 15 wants ref bytes at 15..22, but
        // we only have 6 N-tail bytes. Trim test: use a smaller window.
        // Variant `CGNN→C` at pos 14 instead (3-base deletion, last base N).
        //   pos 14 = 'c', 15 = 'g', 16 = 'n', 17 = 'n'. ref "CGNN", alt "C".
        let mut v = rv(14, 17, "C");
        left_shift(reference, &mut v, 1, true);
        // Pre-fix: rv would slide left into the cacacac window (start <= 12).
        // Post-fix: the right-edge N blocks the slide → start stays at 14.
        assert_eq!(v.start, 14, "start must not slide when right edge is N");
        assert_eq!(v.end, 17, "end must remain at original end");
        assert_eq!(v.alt, "C", "alt must remain at the original anchor");
    }
}
