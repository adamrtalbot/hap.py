# HG001 primitive identity and genotype orientation

This minimized fixture preserves the preprocessing cause exposed by the public
HG001 Platinum Genomes comparison at `chr6:141113702-141113704`. The source
truth contains insertions at both positions. The source query contains an
insertion with the opposite phase and two overlapping multi-allelic records at
141113704, including `A>TGTGTGTGT`.

The fixture reduces that locus to one direct insertion, one SNP, and one
complex allele at a single local anchor. Decomposition turns the complex allele
into the same SNP and an insertion. The direct and realigned insertions are
distinct internal edits even though final VCF padding gives them the same ALT
spelling. Legacy preprocessing therefore retains both ALT entries with
`GT=2/1`; exact duplicate SNP edits collapse to one `GT=1/1` record. The
realigned heterozygous primitive also uses legacy's reversed unphased
orientation.

Commit `5e308a0` diagnosed and fixed primitive identity through the external
sort rewrite, duplicate-allele projection, and the genotype-orientation rule.
On its parent, the paired PREPY comparator reports `result.vcf.gz` different:
hap-rs emits one insertion ALT with `GT=1/1`. With the fix, both implementations
emit duplicate insertion ALTs with `GT=2/1` and the SNP with `GT=1/1`.

The ordinary HAPPY seam was evaluated first with the five public records and
several smaller matching forms. Its classified output projected or reshaped
the preprocessing difference, leaving pre-fix and post-fix hap-rs artifacts
byte-identical or exposing unrelated matching behavior. PREPY is therefore the
highest existing command seam that directly observes this diagnosed cause.
