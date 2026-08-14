# HG001 graph record ordering

This behavior-based example preserves the record-order relationship that caused
the public HG001 Platinum Genomes comparison to diverge at
`chr11:67884797-67884860`. It is minimized to three records on a local reference:

- truth lists a heterozygous same-length MNP before a heterozygous SNV at the
  same coordinate;
- query carries the same SNV with the opposite phase; and
- the SNV spelling shared across inputs must reach the haplotype graph before
  the truth-only MNP, regardless of truth file order.

With the binary from `31e249f^`, the unchanged paired-output comparator reports
`result.vcf.gz` different at line 52: the truth-only MNP is `FN` with `BK=.`
instead of the legacy implementation's `BK=lm`. Commit `31e249f` fixes graph
record ordering; with that fix the complete comparison has no differences.

The samplesheet uses the supported `--no-json`, `--no-roc`, and
`--no-write-counts` output options because metrics and report rendering are
tracked separately. Both implementations expose the same
artifact inventory, including `result.vcf.gz`, which remains compared in full.
The equal-length substitution shape keeps every retained aggregate ROC row in
the SNP family, where transition and transversion cells are applicable. This
removes the former fixture's unrelated INDEL empty-cell rendering difference
without changing the comparison rules or artifact requirements.
