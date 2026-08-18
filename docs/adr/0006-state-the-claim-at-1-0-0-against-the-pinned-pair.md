# State the claim at 1.0.0 against the pinned pair

The drop-in claim is a statement about two named things: one hap-rs version and
one legacy identity. It is first made at hap-rs 1.0.0. No 0.x release carries it.

## The claim

> hap-rs 1.0.0 is a drop-in replacement for the Legacy implementation, pinned as
> `community.wave.seqera.io/library/happy-0.3.15@sha256:4dda6b77c0b1bd778300ea33c498f9d6267537c05d1b655bf9385a11b8702d1c`
> together with `verification/containers/happy-0.3.15.conda-lock.txt`.
>
> For any invocation the pinned parsers accept, `hap <subcommand>` produces the
> same exit status and the same artifacts as the corresponding legacy tool,
> agreeing on every data cell, every categorical verdict label, and record
> identity and order. Observed on linux/amd64. Legacy has no build for any other
> platform, so no comparison exists elsewhere, and hap-rs behaviour on other
> platforms is defined by hap-rs and tested against its own expectations.
>
> The only permitted departures are the two entries in the sealed exemption
> register. A single covered invocation where a claimed observable disagrees, and
> which is not a register entry, falsifies this.

The falsifier is cheap to state and cheap to run: one invocation, one disagreeing
cell or label or record position, not on the register.

## What the claim attaches to

Both sides are named, because neither alone identifies a behaviour.
`0002-drop-in-claim-against-one-unpatched-legacy-image.md` spent a whole decision
establishing that `hap.py 0.3.15` named two different behaviours, so "drop-in with
hap.py" is not checkable. The claim therefore pairs a hap-rs version with the
digest and its conda lock, and
`0004-pin-the-legacy-baseline-to-one-container-identity.md` already governs what
happens when either side moves: a pin move re-baselines, and legacy's fresh output
becomes the expectation.

There is no calendar expiry. A date is unenforceable and nothing measures it.
hap.py 0.3.15 is the last legacy release, so the legacy side of the pair is not
expected to move at all.

1.0.0 rather than the first published version. A drop-in claim on a 0.x release is
a promise semver permits the next release to break, and the register's length
being checkable at release only pays off at the release the register is sealed for.

## Not claimed

Provenance fields: version, timestamp, command line, generated descriptions,
runtime VCF headers. Message text and which stream carries it. Malformed
invocations. A reference supplied through `HGREF` or `HG19`. VCF on standard
input. hap-rs's own added options. `--bam` feature extraction. Platforms other
than linux/amd64.

Runtime and resource use is not a compatibility question. Neither is independent
scientific correctness: this is an agreement claim, and hap-rs reproduces legacy
including where legacy is wrong.

The list stays at class level. Naming every provenance field in every artifact
would create a second inventory that can drift from the one
[Inventory product artifacts and observable side effects](https://github.com/adamrtalbot/hap.py/issues/22)
and
[Choose the release artifact-equivalence contract](https://github.com/adamrtalbot/hap.py/issues/17)
are filed to build, which is the failure
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` rejected for the
option surface and `0004-pin-the-legacy-baseline-to-one-container-identity.md`
rejected for the dependency list.

Nothing is said about what a caller may not depend on, beyond the exclusions
above. Standard error is not output. hap-rs is free to emit clearer error messages
than legacy, and doing so is not a divergence.

## The register is permanent, and its first entry is smaller than it reads

Both entries stay, and neither has a removal date.
`--engine-vcfeval-path` and `--engine-vcfeval-template` are supported as they are
today: accepted, ignored, warned about. The scheduled removal at 1.0.0 is
withdrawn.

Measurement withdrew it. Probes run in the pinned image on linux/amd64 with
`RTG_JAVA_OPTS=-Xint`, against the `vcfeval-matrix` fixture:

| invocation | legacy exit | effect on claimed observables |
|---|---|---|
| `--engine-vcfeval-template <SDF matching the reference>` | 0 | none. `summary.csv`, `extended.csv`, every `roc.*.csv.gz` and the VCF body are byte-identical to the run without it |
| no template | 0 | legacy builds its own SDF, warning that a supplied one would be faster |
| `--engine-vcfeval-template /definitely/missing/template.sdf` | 0 | none. Legacy ignores a nonexistent template and falls back to building one |
| `--engine-vcfeval-template <valid SDF for a different reference>` | **1** | "no sequence names in common between the reference and the supplied variant sets" |
| `--engine-vcfeval-path /definitely/missing/rtg` | **1** | "Error running rtg tools. Return code was 127" |

The whole difference between the first two runs is provenance: the `##CL=` header,
`runinfo.json`, and `metrics.json`'s `description`, `runInfo` and `timestamp`. So
where the rtg path works and the SDF matches the reference, hap-rs ignoring both
options agrees with legacy everywhere the claim reaches. The divergence is
confined to a path that does not resolve to a working rtg, or an SDF that
disagrees with the reference, where legacy exits 1 and hap-rs exits 0.

Removing the options would have widened that. Legacy runs
`--engine-vcfeval-template <matching SDF>` to exit 0 with identical results, and
that is the exact spelling `verification/modules/happy.nf` issues for all five
legacy vcfeval rows. hap-rs would reject the option and exit non-zero, turning an
invocation that agrees today into one that fails. The removal traded two narrow
disagreements for one broad one, on the option pair every vcfeval caller passes.

The register's wording stays broad: the options are accepted and ignored. The
measured triggers are recorded here rather than in the register entry.

## Deprecations

The two deprecations dated "before 1.0.0" resolve in opposite directions.

The vcfeval option pair is no longer deprecated. It is supported.

`pre` and `quantify` returning exit 0 on an unknown option is fixed rather than
scheduled. It is legacy behaviour worth not reproducing, it is already outside the
claim because
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` puts malformed
invocations outside the covered surface, so changing it narrows nothing and breaks
no claim. An unknown option returns a non-zero status now. Matching the exit 1
that `germline` already returns for an invalid argument is the natural choice, and
the exact status is a hap-rs implementation question rather than a policy one.

The claimed version therefore schedules no narrowing of its own surface.

## Platforms

The claim names linux/amd64 because that is where it was observed, and because
the conda lock is entirely `linux-64`, so legacy has no other build to observe.
arm64 is expected to work and is not this claim's evidence.

That is the same shape as `--bam`: hap-rs behaviour exists, no legacy behaviour
exists to compare it against, so hap-rs behaviour is normative there. Which
platforms hap-rs supports stays with
[Define supported execution environments and the reproducibility boundary](https://github.com/adamrtalbot/hap.py/issues/29),
now with a sharper question: only one supported platform can ever carry a legacy
comparison.

## Considered options

A per-release restatement of the claim, and a claim spanning a whole 1.x line,
were both rejected in favour of the version pair. The first is not falsifiable
without naming the legacy side; the second breaks the moment a 1.x patch touches a
covered observable.

A sentence in the claim governing how the claim itself may later narrow was
rejected as over-specification. Semver covers it, and
`0004-pin-the-legacy-baseline-to-one-container-identity.md` already covers the
re-baseline half.

`CONTEXT.md`'s definition of intentional divergence said a divergence differs from
a diagnosed legacy defect. That framing is dropped rather than widened. A register
entry states hap-rs behaviour: `--engine-vcfeval-path` and
`--engine-vcfeval-template` are ignored because the native engine reads a FASTA and
needs no SDF bundle, and `hap validate --help` exits 0. Whether legacy's
corresponding behaviour is a defect does not enter it, and the register entries and
the vocabulary now read the same way.

## Consequences

The claim statement lives at the top of the compatibility policy, replacing the
proto-claim already there. Ordinary users do not need to read this ADR.

Three follow-ons fall outside this decision, all in hap-rs rather than policy:
the `--engine-vcfeval-path` and `--engine-vcfeval-template` warning text promises
a 1.0.0 removal that is no longer happening; the unknown-option exit status for
`pre` and `quantify` becomes non-zero; and the compatibility inventory's
regression evidence for the vcfeval pair should exercise the two measured
divergence triggers, not only the ignored-flags case that
`matrix_vcfeval_deprecated_flags` covers today.

`0002-drop-in-claim-against-one-unpatched-legacy-image.md` records the register at
two entries and dates the `pre` and `quantify` behaviour to 1.0.0. Both entries
survive, both are permanent, and neither the vcfeval removal nor the 1.0.0 date
does; a note there points here.
