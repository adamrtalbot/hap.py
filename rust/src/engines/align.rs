//! Pairwise primitive alignment for preprocessing stage 5.
//!
//! Port of the C++ helper at `src/c++/lib/align/Alignment.cpp:299-377`
//! (`realignRefVar`, list-output form). The legacy
//! `VariantPrimitiveSplitter` relies on this to break complex REF/ALT pairs
//! into primitive SNP / pure-insertion / pure-deletion pieces before the
//! second normalization + aggregation passes.
//!
//! We compute a global pairwise alignment between `ref_allele` and
//! `alt_allele` using Needleman-Wunsch with linear gap penalties and a
//! deterministic tie-breaking rule (prefer diagonal → horizontal → vertical,
//! which steers gaps to the right and keeps the output consistent with the
//! legacy C++ aligner on the samples used by the chr21 fixtures). The
//! resulting CIGAR is then decomposed into [`RefVar`] primitives.
//!
//! Primitive conventions follow legacy semantics (inherited from
//! `partial_credit::RefVar`):
//!
//! - **SNP**: `start == end`, `alt` is a single base.
//! - **Deletion** (k ≥ 1 bases deleted starting at `start`): `end = start + k - 1`,
//!   `alt` is empty.
//! - **Insertion** (k ≥ 1 bases inserted immediately before `start`):
//!   `end = start - 1` (so reflen = 0), `alt` is the inserted sequence.
//!
//! `realign_ref_var` early-exits to [`to_primitives`] when either side has
//! fewer than 2 bases — that matches the legacy guard and keeps simple SNPs
//! / trivial indels on the fast path.

use crate::engines::partial_credit::RefVar;

const MATCH_SCORE: i32 = 2;
const MISMATCH_SCORE: i32 = -4;
const GAP_SCORE: i32 = -4;

/// Decompose a REF/ALT pair into legacy-style primitives. `ref_start` is the
/// 1-based absolute position of the first base of `ref_allele` on the contig.
/// `alt_allele` may be empty only when `ref_allele` is non-empty (pure
/// deletion); both empty yields no primitives.
pub(crate) fn realign_ref_var(
    ref_start: usize,
    ref_allele: &[u8],
    alt_allele: &[u8],
) -> Vec<RefVar> {
    if ref_allele.len() < 2 || alt_allele.len() < 2 {
        return to_primitives(ref_start, ref_allele, alt_allele);
    }
    let cigar = needleman_wunsch(ref_allele, alt_allele);
    cigar_to_primitives(ref_start, ref_allele, alt_allele, &cigar)
}

/// Left-to-right primitive decomposition. Mirrors the C++ `toPrimitives()`
/// (`src/c++/lib/align/RefVar.cpp:514`): walk REF and ALT position by
/// position, emit a single-base SNP for each mismatch, and flush whatever
/// tail remains as either a pure deletion (REF has more bases) or a pure
/// insertion (ALT has more). This is the fast path for `reflen < 2 ||
/// altlen < 2` and also correctly handles simple equal-length substitutions.
pub(crate) fn to_primitives(ref_start: usize, ref_allele: &[u8], alt_allele: &[u8]) -> Vec<RefVar> {
    let mut out = Vec::new();
    let mut rstart = ref_start;
    let rend = ref_start + ref_allele.len().saturating_sub(1);
    let mut pos = 0usize;
    let mut reflen = ref_allele.len() as i64;
    let mut altlen = alt_allele.len() as i64;

    while reflen > 0 && altlen > 0 {
        let r = ref_allele[pos].to_ascii_uppercase();
        let a = alt_allele[pos].to_ascii_uppercase();
        if r != a {
            out.push(RefVar {
                start: rstart,
                end: rstart,
                alt: (a as char).to_string(),
            });
        }
        rstart += 1;
        pos += 1;
        reflen -= 1;
        altlen -= 1;
    }

    if reflen > 0 {
        // Remaining ref bases become a pure deletion.
        out.push(RefVar {
            start: rstart,
            end: rend,
            alt: String::new(),
        });
    } else if altlen > 0 {
        // Remaining alt bases become a pure insertion at `rstart`.
        // RefVar convention: end = start - 1 (reflen == 0).
        let inserted: String = alt_allele[pos..]
            .iter()
            .map(|b| b.to_ascii_uppercase() as char)
            .collect();
        out.push(RefVar {
            start: rstart,
            end: rstart.saturating_sub(1),
            alt: inserted,
        });
    }
    let _ = rend; // only used in the `reflen > 0` branch; tolerate for the insertion path

    out
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Op {
    /// Match or mismatch — consumes one base from ref and alt.
    M,
    /// Insertion — consumes one base from alt only (extra base vs reference).
    I,
    /// Deletion — consumes one base from ref only (base missing in alt).
    D,
}

/// Trace-source flag for backtrace tie-breaking.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Trace {
    Diag,
    Up,   // from (i-1, j): ref advanced without alt → deletion
    Left, // from (i, j-1): alt advanced without ref → insertion
}

fn score_of(a: u8, b: u8) -> i32 {
    if a.eq_ignore_ascii_case(&b) {
        MATCH_SCORE
    } else {
        MISMATCH_SCORE
    }
}

/// Global pairwise alignment returning a CIGAR as a flat `Op` list (one entry
/// per consumed base, not run-length encoded). Tie-breaks prefer diagonal
/// steps, then `Up` (deletion), then `Left` (insertion), which in practice
/// places gaps at the rightmost valid position — matching the orientation
/// legacy `ksw`-style alignment produces on the REF/ALT pairs prepy stresses.
fn needleman_wunsch(ref_allele: &[u8], alt_allele: &[u8]) -> Vec<Op> {
    let m = ref_allele.len();
    let n = alt_allele.len();

    let mut score = vec![vec![0i32; n + 1]; m + 1];
    let mut trace = vec![vec![Trace::Diag; n + 1]; m + 1];

    for i in 1..=m {
        score[i][0] = GAP_SCORE * i as i32;
        trace[i][0] = Trace::Up;
    }
    for j in 1..=n {
        score[0][j] = GAP_SCORE * j as i32;
        trace[0][j] = Trace::Left;
    }

    for i in 1..=m {
        for j in 1..=n {
            let diag = score[i - 1][j - 1] + score_of(ref_allele[i - 1], alt_allele[j - 1]);
            let up = score[i - 1][j] + GAP_SCORE; // delete ref base
            let left = score[i][j - 1] + GAP_SCORE; // insert alt base

            // Tie-break rule: strict preference for Diag only when it is
            // *strictly* better than both gap directions. On a tie between
            // Diag and Left we prefer Left so the insertion ends up on the
            // right-hand side of the alignment — this matches legacy's
            // left-to-right `toPrimitives` walk which always emits tail
            // insertions after all matches have been consumed. Among gap
            // directions, Left (insertion) beats Up (deletion) on ties.
            let (best_score, best_trace) = if diag > up && diag > left {
                (diag, Trace::Diag)
            } else if left >= up {
                (left, Trace::Left)
            } else {
                (up, Trace::Up)
            };

            score[i][j] = best_score;
            trace[i][j] = best_trace;
        }
    }

    let mut ops = Vec::with_capacity(m + n);
    let mut i = m;
    let mut j = n;
    while i > 0 || j > 0 {
        match trace[i][j] {
            Trace::Diag if i > 0 && j > 0 => {
                ops.push(Op::M);
                i -= 1;
                j -= 1;
            }
            Trace::Up if i > 0 => {
                ops.push(Op::D);
                i -= 1;
            }
            Trace::Left if j > 0 => {
                ops.push(Op::I);
                j -= 1;
            }
            // Fall-back for degenerate edges (i==0 or j==0 where trace is
            // forced by the border).
            _ => {
                if i > 0 {
                    ops.push(Op::D);
                    i -= 1;
                } else {
                    ops.push(Op::I);
                    j -= 1;
                }
            }
        }
    }
    ops.reverse();
    ops
}

/// Walk the per-base CIGAR and emit primitives. Runs of identical `M`
/// positions with matching bases are silently consumed; any mismatch under
/// `M` becomes a single-base SNP. `I` runs become pure insertions and `D`
/// runs become pure deletions, each in the legacy `RefVar` convention.
fn cigar_to_primitives(
    ref_start: usize,
    ref_allele: &[u8],
    alt_allele: &[u8],
    ops: &[Op],
) -> Vec<RefVar> {
    let mut primitives = Vec::new();
    let mut ref_index = 0usize;
    let mut alt_index = 0usize;

    let mut k = 0usize;
    while k < ops.len() {
        match ops[k] {
            Op::M => {
                let ref_byte = ref_allele[ref_index];
                let alt_byte = alt_allele[alt_index];
                if !ref_byte.eq_ignore_ascii_case(&alt_byte) {
                    primitives.push(RefVar {
                        start: ref_start + ref_index,
                        end: ref_start + ref_index,
                        alt: std::str::from_utf8(&[alt_byte.to_ascii_uppercase()])
                            .expect("ASCII base")
                            .to_string(),
                    });
                }
                ref_index += 1;
                alt_index += 1;
                k += 1;
            }
            Op::I => {
                // Collapse a run of I ops into a single insertion primitive.
                let mut run_end = k;
                while run_end < ops.len() && ops[run_end] == Op::I {
                    run_end += 1;
                }
                let inserted: String = alt_allele[alt_index..alt_index + (run_end - k)]
                    .iter()
                    .map(|b| b.to_ascii_uppercase() as char)
                    .collect();
                // Insertion before position `ref_start + ref_index` in 1-based
                // ref coordinates. End = start - 1 (reflen = 0).
                let insertion_start = ref_start + ref_index;
                primitives.push(RefVar {
                    start: insertion_start,
                    end: insertion_start.saturating_sub(1),
                    alt: inserted,
                });
                alt_index += run_end - k;
                k = run_end;
            }
            Op::D => {
                // Collapse a run of D ops into a single deletion primitive.
                let mut run_end = k;
                while run_end < ops.len() && ops[run_end] == Op::D {
                    run_end += 1;
                }
                let del_start = ref_start + ref_index;
                let del_end = del_start + (run_end - k) - 1;
                primitives.push(RefVar {
                    start: del_start,
                    end: del_end,
                    alt: String::new(),
                });
                ref_index += run_end - k;
                k = run_end;
            }
        }
    }

    primitives
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn apply_primitives(reference: &[u8], ref_start: usize, primitives: &[RefVar]) -> Vec<u8> {
        let mut reconstructed = reference.to_vec();
        for primitive in primitives.iter().rev() {
            let offset = primitive.start - ref_start;
            let replaced = if primitive.end >= primitive.start {
                primitive.end - primitive.start + 1
            } else {
                0
            };
            reconstructed.splice(offset..offset + replaced, primitive.alt.bytes());
        }
        reconstructed
    }

    proptest! {
        #[test]
        fn primitive_edits_reconstruct_the_alternate_haplotype(
            reference in proptest::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..10),
            alternate in proptest::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 0..10),
            ref_start in 1usize..1000,
        ) {
            let primitives = realign_ref_var(ref_start, &reference, &alternate);
            prop_assert_eq!(apply_primitives(&reference, ref_start, &primitives), alternate);
        }
    }

    #[test]
    fn caa_to_c_yields_two_base_deletion() {
        // REF=CAA at pos 100, ALT=C. Expect a single deletion primitive
        // covering the trailing "AA" at positions 101..=102.
        let primitives = realign_ref_var(100, b"CAA", b"C");
        assert_eq!(primitives.len(), 1);
        assert_eq!(primitives[0].start, 101);
        assert_eq!(primitives[0].end, 102);
        assert_eq!(primitives[0].alt, "");
    }

    #[test]
    fn ca_to_caa_yields_single_base_insertion() {
        // REF=CA at pos 100, ALT=CAA. Expect a pure insertion of "A"
        // immediately after the 2-base match — i.e. before ref position 102.
        // RefVar convention: end = start - 1 (reflen 0).
        let primitives = realign_ref_var(100, b"CA", b"CAA");
        assert_eq!(primitives.len(), 1);
        assert_eq!(primitives[0].start, 102);
        assert_eq!(primitives[0].end, 101);
        assert_eq!(primitives[0].alt, "A");
    }

    #[test]
    fn cgcg_to_ctct_yields_two_snps() {
        // Two mismatches on an otherwise aligned 4-base block.
        let primitives = realign_ref_var(50, b"CGCG", b"CTCT");
        assert_eq!(primitives.len(), 2);
        assert_eq!(primitives[0].start, 51);
        assert_eq!(primitives[0].end, 51);
        assert_eq!(primitives[0].alt, "T");
        assert_eq!(primitives[1].start, 53);
        assert_eq!(primitives[1].end, 53);
        assert_eq!(primitives[1].alt, "T");
    }

    #[test]
    fn cga_to_c_yields_two_base_deletion() {
        // Global NW prefers the aligned "C" to stay matched, putting the
        // deletion at the end. Other alignments score worse under
        // MATCH=2, MISMATCH=-4, GAP=-4.
        let primitives = realign_ref_var(10, b"CGA", b"C");
        assert_eq!(primitives.len(), 1);
        assert_eq!(primitives[0].start, 11);
        assert_eq!(primitives[0].end, 12);
        assert_eq!(primitives[0].alt, "");
    }

    #[test]
    fn cc_to_ccc_yields_single_base_insertion() {
        // All-A homopolymer extension — the insertion anchors after the last
        // matched base (tie-break prefers diagonal, pushing the gap right).
        let primitives = realign_ref_var(200, b"CC", b"CCC");
        assert_eq!(primitives.len(), 1);
        assert_eq!(primitives[0].start, 202);
        assert_eq!(primitives[0].end, 201);
        assert_eq!(primitives[0].alt, "C");
    }

    #[test]
    fn trivial_snp_falls_back_to_to_primitives() {
        // reflen == 1 hits the fast path.
        let primitives = realign_ref_var(5, b"A", b"G");
        assert_eq!(primitives.len(), 1);
        assert_eq!(primitives[0].start, 5);
        assert_eq!(primitives[0].end, 5);
        assert_eq!(primitives[0].alt, "G");
    }

    #[test]
    fn to_primitives_single_base_deletion() {
        let primitives = to_primitives(5, b"A", b"");
        assert_eq!(primitives.len(), 1);
        assert_eq!(primitives[0].start, 5);
        assert_eq!(primitives[0].end, 5);
        assert_eq!(primitives[0].alt, "");
    }

    #[test]
    fn to_primitives_single_base_insertion() {
        let primitives = to_primitives(5, b"", b"G");
        assert_eq!(primitives.len(), 1);
        assert_eq!(primitives[0].start, 5);
        assert_eq!(primitives[0].end, 4);
        assert_eq!(primitives[0].alt, "G");
    }
}
