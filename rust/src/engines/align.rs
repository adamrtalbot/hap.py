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

// Legacy klib (ksw) global alignment scores — the BWA-MEM defaults declared in
// hap.py `AlignmentParameters` (src/c++/include/Alignment.hh): match +1,
// mismatch -4, gap-open 6, gap-extend 1, with an affine cost of
// `GAP_OPEN + len * GAP_EXTEND` for a gap of `len` bases. Affine gaps make one
// contiguous indel cheaper than several spread ones, so complex alleles keep
// their insertions/deletions together (and prefer SNPs over a delete+insert
// pair) exactly as legacy does.
const MATCH_SCORE: i32 = 1;
const MISMATCH_SCORE: i32 = -4;
const GAP_OPEN: i32 = 6;
const GAP_EXTEND: i32 = 1;
const NEG_INF: i32 = i32::MIN / 4;

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

/// Which affine-DP layer a cell's optimum came from, for backtrace.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Layer {
    /// Diagonal step — ref and alt both advance (match or mismatch).
    Match,
    /// Gap in ref — alt advances only (insertion).
    Insert,
    /// Gap in alt — ref advances only (deletion).
    Delete,
}

fn score_of(a: u8, b: u8) -> i32 {
    if a.eq_ignore_ascii_case(&b) {
        MATCH_SCORE
    } else {
        MISMATCH_SCORE
    }
}

/// Affine-gap global alignment (Gotoh) returning a per-base CIGAR as a flat
/// `Op` list. Three layers are tracked: `h` (ends in a diagonal step), `e`
/// (ends in an insertion — a gap in ref), and `f` (ends in a deletion — a gap
/// in alt). A gap of `len` costs `GAP_OPEN + len * GAP_EXTEND`.
///
/// Tie-breaking mirrors legacy `ksw_global`, which right-aligns gaps: on equal
/// score a diagonal step is preferred over opening a gap, and among gap layers
/// the alignment prefers to *keep extending* an already-open gap rather than
/// close and reopen it. This keeps a complex allele's indel contiguous and at
/// the right-hand end of the aligned block.
fn needleman_wunsch(ref_allele: &[u8], alt_allele: &[u8]) -> Vec<Op> {
    let m = ref_allele.len();
    let n = alt_allele.len();

    // h/e/f score layers and their backtrace pointers.
    let mut h = vec![vec![NEG_INF; n + 1]; m + 1];
    let mut e = vec![vec![NEG_INF; n + 1]; m + 1];
    let mut f = vec![vec![NEG_INF; n + 1]; m + 1];
    // For each layer, which layer we stepped from.
    let mut h_from = vec![vec![Layer::Match; n + 1]; m + 1];
    let mut e_from = vec![vec![Layer::Insert; n + 1]; m + 1];
    let mut f_from = vec![vec![Layer::Delete; n + 1]; m + 1];

    h[0][0] = 0;
    for i in 1..=m {
        // Deletion of the first `i` reference bases.
        f[i][0] = -(GAP_OPEN + i as i32 * GAP_EXTEND);
        f_from[i][0] = Layer::Delete;
    }
    for j in 1..=n {
        // Insertion of the first `j` alternate bases.
        e[0][j] = -(GAP_OPEN + j as i32 * GAP_EXTEND);
        e_from[0][j] = Layer::Insert;
    }

    for i in 1..=m {
        for j in 1..=n {
            // Insertion (gap in ref): alt advanced. Prefer extending an open
            // insertion over reopening on a tie.
            let open_e = h[i][j - 1] - (GAP_OPEN + GAP_EXTEND);
            let extend_e = e[i][j - 1] - GAP_EXTEND;
            if extend_e >= open_e {
                e[i][j] = extend_e;
                e_from[i][j] = Layer::Insert;
            } else {
                e[i][j] = open_e;
                e_from[i][j] = Layer::Match;
            }

            // Deletion (gap in alt): ref advanced.
            let open_f = h[i - 1][j] - (GAP_OPEN + GAP_EXTEND);
            let extend_f = f[i - 1][j] - GAP_EXTEND;
            if extend_f >= open_f {
                f[i][j] = extend_f;
                f_from[i][j] = Layer::Delete;
            } else {
                f[i][j] = open_f;
                f_from[i][j] = Layer::Match;
            }

            // Diagonal step: best predecessor across layers, plus sub score.
            let diag_prev = h[i - 1][j - 1].max(e[i - 1][j - 1]).max(f[i - 1][j - 1]);
            h[i][j] = diag_prev + score_of(ref_allele[i - 1], alt_allele[j - 1]);
            h_from[i][j] =
                if h[i - 1][j - 1] >= e[i - 1][j - 1] && h[i - 1][j - 1] >= f[i - 1][j - 1] {
                    Layer::Match
                } else if e[i - 1][j - 1] >= f[i - 1][j - 1] {
                    Layer::Insert
                } else {
                    Layer::Delete
                };
        }
    }

    // Start in the highest-scoring layer at the corner. On a tie prefer a gap
    // layer (Insert, then Delete) over Match so a trailing gap right-aligns.
    let (mut layer, _) = [
        (Layer::Insert, e[m][n]),
        (Layer::Delete, f[m][n]),
        (Layer::Match, h[m][n]),
    ]
    .into_iter()
    .fold(
        (Layer::Match, NEG_INF),
        |(best_layer, best), (layer, score)| {
            if score > best {
                (layer, score)
            } else {
                (best_layer, best)
            }
        },
    );

    let mut ops = Vec::with_capacity(m + n);
    let mut i = m;
    let mut j = n;
    while i > 0 || j > 0 {
        match layer {
            Layer::Match => {
                ops.push(Op::M);
                layer = h_from[i][j];
                i -= 1;
                j -= 1;
            }
            Layer::Insert => {
                ops.push(Op::I);
                layer = e_from[i][j];
                j -= 1;
            }
            Layer::Delete => {
                ops.push(Op::D);
                layer = f_from[i][j];
                i -= 1;
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
    fn affine_gaps_keep_complex_insertions_contiguous() {
        // CC → TCTCT in a CT repeat. Affine gaps prefer one contiguous
        // insertion (SNP C>T then a 3-base insertion) over three spread single
        // insertions, matching legacy klib. Positions are relative to ref_start.
        let primitives = realign_ref_var(216, b"CC", b"TCTCT");
        assert_eq!(primitives.len(), 2, "{primitives:?}");
        // SNP C>T at the first base.
        assert_eq!((primitives[0].start, primitives[0].end), (216, 216));
        assert_eq!(primitives[0].alt, "T");
        // Contiguous insertion of "TCT" after the second base (reflen 0).
        assert_eq!((primitives[1].start, primitives[1].end), (218, 217));
        assert_eq!(primitives[1].alt, "TCT");
    }

    #[test]
    fn affine_prefers_two_snps_over_delete_insert_on_a_swap() {
        // GC → CG. Two mismatches (-8) beat a delete+insert pair (two gap
        // opens), so a transposition decomposes into two SNPs, not an indel
        // pair.
        let primitives = realign_ref_var(26, b"GC", b"CG");
        assert_eq!(primitives.len(), 2, "{primitives:?}");
        assert_eq!((primitives[0].start, primitives[0].alt.as_str()), (26, "C"));
        assert_eq!((primitives[1].start, primitives[1].alt.as_str()), (27, "G"));
    }

    #[test]
    fn affine_finds_a_clean_middle_insertion() {
        // GA → GCA. A single clean insertion beats a mismatch-plus-insertion,
        // so the extra base is inserted between the two matches.
        let primitives = realign_ref_var(17, b"GA", b"GCA");
        assert_eq!(primitives.len(), 1, "{primitives:?}");
        assert_eq!((primitives[0].start, primitives[0].end), (18, 17));
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
