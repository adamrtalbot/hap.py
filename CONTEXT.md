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
which migration compatibility is evaluated, as observed by running those tools
unmodified in the reference environment.
_Avoid_: Old implementation, Python version

**Reference environment**:
The single immutable container image digest defining authoritative legacy
behavior for every lane, run on linux/amd64. Legacy runs there unmodified; an
observation taken through a patched interpreter, patched source, or library
option override is not evidence of legacy behavior. A setting that configures the
interpreter or VM legacy runs on, or the order the harness schedules it in, is a
harness constraint, allowed only when named in ADR 0004. Its recipe and conda lock are governed artifacts,
because a digest alone cannot say which library produced a number.
_Avoid_: Legacy container, pinned images, the oracle

**Re-baseline**:
What moving any pin does. Legacy is re-run on the new image and its fresh output
becomes the expectation; no earlier observation is preserved or reconciled,
because the harness stores none. Before hap-rs 1.0 a re-baseline invalidates no
claim.
_Avoid_: Regenerate the golden files, update the snapshots

**Drop-in claim**:
The compatibility claim hap-rs makes: every meaningful observable agrees with
the legacy implementation, rather than scientific results or the invocation
interface alone.
_Avoid_: Full parity, byte parity, bug-for-bug compatibility

**Meaningful observable**:
What the drop-in claim covers — data cells, the categorical labels carrying a
classification verdict, and record identity and order.
_Avoid_: Numerical output, the numbers

**Provenance field**:
An emitted value describing the run rather than its result, such as a version,
timestamp, command line, or generated description. Outside the drop-in claim.
_Avoid_: Metadata, excluded field, noise

**Covered invocation surface**:
The set of invocations the drop-in claim applies to: `hap <subcommand>` with the
options, positionals, and input forms the pinned parsers accept. For those, the
claim reaches exit status and produced artifacts, not message text or stream
choice. Malformed invocations are outside it.
_Avoid_: CLI compatibility, supported flags

**Exemption register**:
The sealed list of ratified intentional divergences. Additions require
maintainer sign-off, so its length is checkable at release.
_Avoid_: Waiver list, known differences

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
