# HG001 graph record ordering

This behavior-based example preserves the record-order relationship that caused
the public HG001 Platinum Genomes comparison to diverge at
`chr11:67884797-67884860`. It is minimized to three records on a local reference:

- truth lists a homozygous deletion before a heterozygous insertion at the same
  coordinate;
- query carries the same insertion with the opposite phase; and
- the insertion spelling shared across inputs must reach the haplotype graph
  before the truth-only deletion, regardless of truth file order.

With the binary from `31e249f^`, the unchanged paired-output comparator reports
`result.vcf.gz` different at line 52: the truth-only deletion is `FN` with
`BK=.` instead of the legacy implementation's `BK=lm`. Commit `31e249f` fixes
graph record ordering; at `b55bf91` the normalized VCF has no differences.

The samplesheet uses the supported `--no-json`, `--no-roc`, and
`--no-write-counts` output options because metrics and report rendering are
tracked separately. Both implementations expose the same
artifact inventory, including `result.vcf.gz`, which remains compared in full.
At `b55bf91`, the complete comparison remains red only in
`result.roc.all.csv.gz`: non-applicable report cells are `.` rather than empty.
That float/report-rendering dependency belongs to issue #12 and is not changed
or hidden here.
