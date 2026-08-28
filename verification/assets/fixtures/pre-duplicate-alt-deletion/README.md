# Duplicate-alt deletion stays two rows

This minimized fixture preserves the preprocessing cause seen on the HG003
PEPPER long-read callset at `chr1:118558359`. The source query has two
colocated records on opposite haplotypes: `TGAGA>C` (`1/0`) and `TGAGA>T`
(`0/1`). Decomposing the complex `TGAGA>C` allele yields a `T>C` SNP plus a
`TGAGA>T` deletion, so the locus produces two `TGAGA>T` deletions whose final
VCF padding gives them the same REF/ALT spelling but distinct internal edits.

Legacy `pre.py --decompose --leftshift` keeps the two identical-spelling
**deletions** as separate rows, each retaining its own genotype (`0/1` for the
direct record, then `1/0` for the one recovered from the complex allele). This
differs from identical-spelling **insertions** (see
`hg001-primitive-identity`), which legacy aggregates into one `GT=2/1` het-alt
row. The measured rule is: duplicate padded insertions aggregate, duplicate
padded deletions stay split.

The fixture reduces the real locus onto a 200 bp synthetic contig copied
verbatim from GRCh38 `chr1:118558260-118558459`, so the left-shift context
around the `GAGA` repeat is preserved and both engines resolve the deletion to
the same anchor. Legacy and hap-rs both emit `T>C 1/0`, `TGAGA>T 0/1`,
`TGAGA>T 1/0`.

The same shape recurs at `chr1:158893037`, `chr1:226545091`, and
`chr1:45831780`. PREPY is the highest command seam that observes this cause
directly: the HAPPY seam reshapes it through xcmp's read-side deduplication of
duplicate spellings.
