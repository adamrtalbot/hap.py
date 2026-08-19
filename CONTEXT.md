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
harness constraint, allowed only when named in verification/nextflow.config . Its recipe and conda lock are governed artifacts, because a digest alone cannot say which library produced a number.
_Avoid_: Legacy container, pinned images, the oracle

**Re-baseline**:
What moving any pin does. Legacy is re-run on the new image and its fresh output
becomes the expectation; no earlier observation is preserved or reconciled,
because the harness stores none.
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
choice. Where legacy exits 0, hap-rs exits 0 and the observed set matches; where
legacy exits non-zero, hap-rs's status and artifacts are its own. Malformed
invocations are outside it, and hap-rs reports them with the argument parser's
status. An environment-supplied reference is outside it too: the reference is
passed as an argument. No caller sits inside it, because a caller only issues a
command line.
_Avoid_: CLI compatibility, supported flags

**Observed set**:
Every file whose name begins with the output prefix an invocation names, in the
directory that prefix names. What the drop-in claim compares. The two
implementations' sets must be equal name for name; nothing outside is compared,
so scratch, temporary files, lock files, and logs carry no comparison.
_Avoid_: Output glob, artifact list, result files

**Product artifact**:
A member of the observed set, written because the invocation asked for it.
_Avoid_: Output, deliverable, report

**Campaign evidence**:
What the harness writes about a run rather than what a tool produced:
`comparison.json`, `verification.json`, the streamed artifact records,
`.command.*`, and captured output. Never compared as product.
_Avoid_: Results, output, artifacts

**Equivalence contract**:
How the observed set is compared, one artifact class at a time. Byte identity is
the default; a documented content comparison replaces it only where the encoding
carries a provenance field or is not a function of the result.
_Avoid_: Diff rules, comparison logic, tolerance

**Unstable encoding**:
A serialization whose bytes are not determined by the result alone, such as a
gzip container carrying an mtime, a deflate stream, or an index storing offsets
into one. The reason a class is compared by content rather than by bytes.
_Avoid_: Nondeterministic output, binary noise

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

**Confirmation run**:
A discovery workflow issued through a third-party pipeline rather than this
repository's harness, run once the port is complete to demonstrate hap-rs in a
real caller and surface missed edge cases. It gates nothing and settles no value.
_Avoid_: Admitted caller, integration test, caller lane

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
engines, assemblies, and materially different variant structures.
_Avoid_: Full coverage, exhaustive dataset

**Supported clean machine**:
A machine with Rust/Cargo, Java, Nextflow, nf-test, a supported container
runtime, Git, network access, and sufficient storage, but no legacy scientific
toolchain or pre-fetched corpus data.
_Avoid_: Blank machine, fresh laptop

**Intentional divergence**:
A maintainer-approved hap-rs behavior that legacy does not share, stated as hap-rs
behavior rather than as a departure from legacy. The two in force are that
`--engine-vcfeval-path` and `--engine-vcfeval-template` are ignored, because the
native engine reads a FASTA and needs no SDF bundle, and that
`hap validate --help` exits 0.
_Avoid_: Acceptable difference, exception

**Migration complete**:
The state in which all committed parity cases and the representative corpus
pass parity, determinism, resource, and performance gates, with unsupported
input classes explicitly documented.
_Avoid_: Feature complete, mostly done
