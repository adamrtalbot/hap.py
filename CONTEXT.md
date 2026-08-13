# hap-rs Migration

This context describes the language used to replace the legacy hap.py toolkit
while preserving its supported scientific and artifact contracts.

## Language

**hap-rs**:
The native successor to the legacy hap.py toolkit, distributed as the single
`hap` executable.
_Avoid_: Rust rewrite, replacement binary

**Legacy implementation**:
The pinned hap.py, som.py, pre.py, ftx.py, qfy.py, and vcfcheck behavior against
which migration compatibility is evaluated.
_Avoid_: Old implementation, Python version

**Parity case**:
A named, reproducible comparison of legacy and hap-rs behavior for one defined
input and invocation contract.
_Avoid_: Edge case, test row

**Discovery workflow**:
A full-scale legacy-versus-hap-rs workflow used to discover discrepancies and
confirm fixes in their original real-data context, including on remote compute.
_Avoid_: Full parity case, integration test

**Public-data case**:
One intact real-world truth/query comparison represented by a single row in
the shared public-data samplesheet.
_Avoid_: Discrepancy, regression fixture

**Authoritative regression**:
A durable, focused parity case that reproduces a diagnosed cause from a
discovery workflow and determines whether that behavior remains fixed.
_Avoid_: Discovery workflow, public-data test

**Example case**:
One extremely small legacy-versus-hap-rs scenario represented by a row in an
ordinary verification samplesheet and exercised by nf-test. An example case
may cover multiple defects only when they arise from the same minimal scenario.
_Avoid_: Public-data case, one row per defect

**Discrepancy**:
An observable difference between legacy and hap-rs within a parity case that
has not yet been diagnosed and governed.
_Avoid_: Edge case, mismatch, bug

**Minimized fixture**:
The smallest durable input that reproduces a diagnosed discrepancy without
changing its cause.
_Avoid_: Toy case, synthetic example

**Diagnostic evidence**:
Machine-readable differences and intermediate observations used to localize a
discrepancy without changing the parity gate's acceptance semantics.
_Avoid_: Comparator output, debug log

**Parity gate**:
The acceptance boundary requiring every committed parity case to have no
unapproved differences.
_Avoid_: Test suite, comparator

**Representative corpus**:
A documented set of real datasets selected to exercise supported commands,
engines, callers, assemblies, and materially different variant structures.
_Avoid_: Full coverage, exhaustive dataset

**Supported clean machine**:
A machine with Rust/Cargo, Java, Nextflow, nf-test, a supported container
runtime, Git, network access, and sufficient storage, but no legacy scientific
toolchain or pre-fetched corpus data.
_Avoid_: Blank machine, fresh laptop

**Intentional divergence**:
A maintainer-approved hap-rs behavior that deliberately differs from a
diagnosed legacy defect under an explicit versioned contract.
_Avoid_: Acceptable difference, exception

**Migration complete**:
The state in which all committed parity cases and the representative corpus
pass parity, determinism, resource, and performance gates, with unsupported
input classes explicitly documented.
_Avoid_: Feature complete, mostly done
