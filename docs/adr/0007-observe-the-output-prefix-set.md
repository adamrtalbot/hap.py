# Observe the output prefix set

The inputs and the result files are the parity surface. Everything a run leaves
beside them is auxiliary and carries no comparison: scratch directories,
temporary files, lock files, logs, and both streams.

The observed set is every file whose name begins with the output prefix the
invocation names, in the directory that prefix names. There is no enumerated
list of artifacts. Completeness is symmetric discovery: for a covered invocation
that legacy completes with exit status 0, the two implementations' observed name
sets must be equal, and one name present on one side and absent on the other
falsifies the claim.

## How the prefix is named

It differs per command, so the rule names each one. Measured against the pinned
image:

| command | prefix comes from | example |
|---|---|---|
| `germline` | `-o` / `--report-prefix` | `-o result` observes `result*` |
| `quantify` | `-o` / `--report-prefix` | `-o result` observes `result*` |
| `somatic` | `-o` | `-o result` observes `result*` |
| `ftx` | `-o` | `-o result` observes `result*` |
| `pre` | the output path positional | `result.vcf.gz` observes `result.vcf.gz*` |
| `validate` | `-o` / `--output-json` / `--output-file` | `--output-json result.json` observes `result.json*` |

A prefix that names a directory other than the working directory moves the
observed set with it. The prefix is not a glob the caller writes; it is the
value already on the command line.

## The measured artifact sets

Observed by running the pinned digest
`community.wave.seqera.io/library/happy-0.3.15@sha256:4dda6b77c0b1bd778300ea33c498f9d6267537c05d1b655bf9385a11b8702d1c`
on the small matrix fixtures. This table documents; it does not bind. A missing
condition here is a documentation bug, not a falsified claim, because the rule
above is set equality rather than this list.

| command | default artifacts | measured conditions |
|---|---|---|
| `germline` | `.summary.csv`, `.extended.csv`, `.runinfo.json`, `.metrics.json.gz`, `.roc.all.csv.gz`, `.roc.Locations.{SNP,INDEL}.csv.gz`, `.roc.Locations.{SNP,INDEL}.PASS.csv.gz`, `.vcf.gz`, `.vcf.gz.tbi` | `--no-json` drops `.metrics.json.gz`. `--no-roc` drops the four `.roc.Locations.*` and **keeps** `.roc.all.csv.gz`. The annotated VCF and its index appeared in every measured run, including runs that passed no `-V`, and `--help` advertises no flag that suppresses them. `--roc-regions '*'` and a two-column `--stratification` TSV changed nothing. |
| `quantify` | `.summary.csv`, `.extended.csv`, `.metrics.json.gz`, `.roc.all.csv.gz` | No `.runinfo.json` and no VCF by default. `--write-vcf` adds `.vcf.gz` and `.vcf.gz.tbi`. `--no-write-counts` drops `.extended.csv`. `--no-json` drops `.metrics.json.gz`. `--no-roc` keeps `.roc.all.csv.gz`. |
| `pre` | the output path, plus one index companion | `.vcf.gz` takes `.vcf.gz.tbi`, `.bcf` takes `.bcf.csi`. |
| `somatic` | `.stats.csv`, `.metrics.json`, uncompressed | `--feature-table <name>` adds `.features.csv`. `-P --feature-table <name> --happy-stats` adds `.summary.csv`; no `.extended.csv` appeared in any measured run, including with the typed tables `admix.strelka.snv` and `hcc.strelka.indel`. `--happy-stats` without both exits 1. |
| `ftx` | `.csv` | None found. `--feature-label` adds a column, not a file. |
| `validate` | none | A bare invocation writes no file at all, so its only observable on success is exit status 0. `--output-file <path>` writes that JSON. Legacy has no error-BED option: `--error-bed` is an unrecognised option and exits 1, so hap-rs's `--errors-bed` is a declared extension with nothing to compare against. |

Index companions and compression are inventory facts carried in the artifact
name. The companion extension follows the data extension, and compression state
is visible in the name, `metrics.json.gz` compressed against `metrics.json`
uncompressed. Set equality therefore already binds both. Whether the bytes
inside an artifact must match belongs to
[Choose the release artifact-equivalence contract](https://github.com/adamrtalbot/hap.py/issues/17).

## Product artifact against campaign evidence

A product artifact is what the tool writes inside the observed set because the
invocation asked for it. Campaign evidence is what the harness writes about a
run: `comparison.json`, `verification.json`, the `jsonl.gz` streams,
`.command.*`, and captured output. Evidence is never compared as product, and
capturing a stream does not move its content inside the claim, which
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` already placed
outside.

## Exit status, narrowing the claim

`0006-state-the-claim-at-1-0-0-against-the-pinned-pair.md` claims "the same exit
status and the same artifacts" for any invocation the pinned parsers accept.
That is narrowed here.

For a covered invocation where legacy exits 0, hap-rs exits 0 and the observed
name sets are equal. Where legacy exits non-zero, neither hap-rs's exit status
nor its artifacts are claimed, and hap-rs's own convention governs: 0 for
success, including `--help` and `--version`; 1 for a failure hap-rs detects;
the ordinary Unix status where the cause is external, such as 137 for an
out-of-memory kill.

Legacy's failure statuses are one of the behaviours the port sets out to
correct, so reproducing them was never wanted. A partial artifact set on a
failure path records where execution happened to stop, and a claim over it would
be a contract on internal write ordering. Measured: `germline`
`--scratch-prefix ./scr` with `./scr` absent exits 1 with `result.runinfo.json`
already on disk. argparse accepts that invocation, so it is covered rather than
malformed, which is what makes the narrowing necessary rather than academic.

## Artifact names reproduce legacy

hap-rs emits the legacy names, including where legacy is inconsistent: germline
writes `metrics.json.gz` while somatic writes an uncompressed `metrics.json`.
Set equality is what binds, so renaming or recompressing an artifact would
falsify the claim, and the exemption register is sealed at two entries. Changing
a name later is a documented interface change of the same kind as
`hap.py` becoming `hap germline`, not an exemption.

Where a nicer interface costs nothing it is free to happen: a file whose name
does not begin with the output prefix is outside the observed set, so hap-rs may
add artifacts there without touching the claim.

## Auxiliary writes

Nothing outside the observed set is compared, and the measurements show why an
unbounded reading is unworkable. Every legacy Python tool writes `.pyc` files
into `/opt/conda` on import, which nothing native can reproduce. Legacy
`germline` and `pre` leak an orphan `tmpXXXXXX.vcf.gz.tbi` into `TMPDIR`,
`germline` doing so even when `--scratch-prefix` is supplied; `somatic` leaves a
`TMPDIR` directory behind; `ftx` creates its scratch there and deletes it;
`quantify` and `validate` wrote nothing there. With `--keep-scratch` the
retained trees carry per-run random infixes, `truth.pp0Gt9Oi.vcf.gz` and
`hap.py.result.cEGAD7.vcf.gz`, so nobody can compare their contents.

This corrects a premise in
`0005-keep-the-claim-on-the-command-line.md` without changing its decision. That
ADR reads as though legacy consults no `TMPDIR` and as though `pre`, `quantify`
and `ftx` create no scratch. Measured, legacy consults `TMPDIR` in four of the
six tools, `pre` writes there and leaks, and `ftx` writes there and cleans up.
Not consulting `TMPDIR` is an hap-rs design choice, not agreement with legacy.

hap-rs writes only inside the output directory and inside an explicitly
requested scratch prefix. This is normative hap-rs behaviour rather than
agreement with legacy, so it sits beside the claim rather than inside it. Two
consequences follow, both hap-rs changes outside this map:

- The publication lock is removed. `rust/src/output.rs:638` and
  `rust/src/adapters/vcf.rs:923` place a lock under
  `$TMPDIR/hap-rs-publication-locks`, and `output.rs` uses a single
  `global.lock` shared by every hap process on the host, which serialises
  unrelated runs and fails outright where another user owns the path.
- Scratch stops coming from `TMPDIR`. `rust/src/application/somatic.rs:65` and
  `rust/src/application/ftx/mod.rs:56` derive it there today. ADR 0005 named
  `engines/roc.rs:47` and `roc_publication.rs:73`; neither exists now.

## Considered options

Observing the whole working directory was rejected. It would catch a second
hap-rs artifact that a caller's glob matches while the output prefix does not,
which is the residual risk
`0005-keep-the-claim-on-the-command-line.md` records, but it also drags staged
inputs, retained scratch trees with random names, and any future auxiliary file
into the comparison. That risk stays open and is accepted: the outputs matter,
the auxiliary files do not.

An enumerated closed set per command and flag combination was rejected. It is a
second inventory that drifts from the image, which is the failure
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` rejected for the
option surface and `0004-pin-the-legacy-baseline-to-one-container-identity.md`
rejected for the dependency list. `--no-roc` keeping `.roc.all.csv.gz` is the
kind of exception such a matrix has to encode by hand.

Having no rule at all, on the reading that
`0005-keep-the-claim-on-the-command-line.md` owes none, was also rejected. The
product documentation was wrong about the germline VCF and the somatic reports
at the time this was written, so prose with nothing binding underneath it does
not stay true.

Claiming `TMPDIR` behaviour was rejected on the measurement above.

## Consequences

`modules/diff.nf` already implements the rule: it globs the prefix on both
sides and reports an `artifact_set` difference when the name sets differ. This
decision writes down what the comparator does, so no comparator change follows
from it.

The product documentation is corrected against the measurements: the germline
annotated VCF is not conditional on `--write-vcf`, `--no-roc` does not remove
`.roc.all.csv.gz`, somatic's `--happy-stats` reports need `-P` and a feature
table and produced only `.summary.csv`, and the report-family table in the
reference no longer calls the germline VCF optional.

`CONTEXT.md` gains observed set, product artifact, and campaign evidence.

One statement in that documentation is inferred rather than measured. The
germline page now says the annotated VCF is written by default, which is legacy
behaviour plus the fact that the germline lane passes artifact-set equality
without `-V`. hap-rs was not run in the session that wrote this. If hap-rs gates
the VCF behind `-V`, that is a parity bug to fix, not a documentation error to
revert.

[Choose the release artifact-equivalence contract](https://github.com/adamrtalbot/hap.py/issues/17)
inherits the observed set as the thing whose members it decides equivalence for,
one artifact at a time, and inherits exit status as claimed only where legacy
exits 0.
