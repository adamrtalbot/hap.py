# Phase 1 parity baseline (historical)

> This file preserves the useful findings from the original migration log. It
> is historical evidence, not a current release-status report. The live gate
> contract is defined by `verification/nextflow.config`, the samplesheets under
> `verification/assets/`, and `verification/tests/main.nf.test`.

## Original snapshot

The first end-to-end parity work used real chr21 NA12878 data and the legacy
`hap.py_rtg-tools:3d21f155636b42e0` image. That early matrix contained three
FTX rows, three preprocessing rows, and four germline rows; the somatic Rust
lane was still disabled.

The snapshot exposed four broad classes of migration work:

- repeated full-reference allocations made validation and preprocessing
  effectively quadratic on chr21;
- gzip output was not yet BGZF/index compatible;
- preprocessing lacked parts of the legacy allele normalization and primitive
  decomposition pipeline;
- germline comparison still differed in VCF classification, summaries,
  extended tables, and ROC output.

Those observations explain the staged implementation history below. They do
not describe the current product surface or current matrix.

## Historical milestones

### Preprocessing and indexed output

The first phase replaced per-record reference copies with byte-slice access,
made REF validation case-insensitive for soft-masked FASTA, and moved VCF.gz
output to BGZF. The preprocessing port then added the legacy normalization
sequence: allele filtering and splitting, left shifting, primitive
decomposition, location aggregation, genotype/depth projection, INFO cleanup,
and deterministic FORMAT ordering.

The historical three-case prepy subset moved from a 1,013-line body diff to
matching VCF bodies and indexes. This was evidence for that subset only; the
current prepy contract is the 26-row samplesheet.

### FTX feature extraction

FTX work moved preprocessing and caller-specific feature extraction into the
Rust code path. It established deterministic CSV rendering and shared Strelka,
MuTect, VarScan2, and Pisces feature handling, including BAM-derived depth
normalization.

The original three-case FTX subset matched its oracle artifacts by the end of
that phase. The current governed FTX matrix has 11 rows, supplemented by the
10-case `verification/scripts/verify-ftx-options.sh` direct oracle.

### Germline comparison

The germline phase ported the comparison and quantification behavior instead
of retaining legacy-output shims. Work included preprocessing reuse, genotype-
aware variant typing, block/haplotype matching, confidence-region annotation,
FP/FN subclassification, summary and extended-table aggregation, and ROC
serialization.

The final recorded checkpoint closed the VCF-body differences for the four
original chr21 germline rows. That checkpoint predates the current 12-row
happy samplesheet and is not a substitute for a full release-gate run.

## Current product boundary

The migration produced one pure-Rust product binary, `hap`, with these public
subcommands:

- `hap germline`
- `hap somatic`
- `hap pre`
- `hap ftx`
- `hap quantify`
- `hap validate`

The product implements BCF and CSI handling in Rust and has no Python, C, or
C++ runtime dependency. The legacy container and auxiliary command-line tools
belong only to the external parity harness.

The separate `verify-fixtures` helper is not part of a default build. Cargo
requires the `verification` feature for every verifier build or invocation,
for example:

```bash
cargo run --features verification --bin verify-fixtures -- check-rust
```

## Current verification contract

The default matrix contains 63 samplesheet rows:

| Lane | Rows |
|---|---:|
| happy | 12 |
| sompy | 12 |
| prepy | 26 |
| ftxpy | 11 |
| qfy | 1 |
| vcfcheck | 1 |

The harness pins Nextflow 26.04.6, nf-test 0.9.5, and an immutable legacy image.
It checks complete samplesheet publication, symmetric artifact sets, structured
or exact output comparison as appropriate, and record retrieval through TBI or
CSI indexes. The qfy and vcfcheck lanes govern `hap quantify` and `hap validate`
alongside the four original tool families.

Focused direct-oracle coverage is provided by:

- `verification/scripts/verify-ftx-options.sh` — 10 FTX cases;
- `verification/scripts/verify-somatic-options.sh` — 3 somatic cases.

These statements describe what the gate checks. A particular revision's
release result must come from its completed nf-test run, not from this
historical baseline.
