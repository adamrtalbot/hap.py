# Pin the legacy baseline to one container identity

The legacy reference baseline is a container digest URI plus the conda lock
extracted from that image. Two containers carry an identity: the legacy
reference and the comparator. Nothing else is pinned.

## What the identity is

Wave builds and stores the reference image, so the image is its own archive.
The digest fixes the contents and the lock names them.
`verification/containers/happy-0.3.15.conda-lock.txt` lists 108 packages with
URL and md5, regenerated from the image by
`verification/scripts/dump-conda-lock.sh`, which takes the image as a required
argument so a stale default cannot produce a governed artifact. We store no copy
of the image and keep no separate manifest file.

The lock already pins everything that would otherwise be listed separately:
rtg-tools 3.12.1, openjdk 11.0.8, libstdcxx 16.1.0, bcftools 1.17,
samtools 1.18, python 2.7.15, numpy 1.12.1, scipy 1.2.1, and pandas 0.20.3.
No second list is kept. A second list can drift from the image
without anyone noticing, which is the failure
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` already rejected for
the option surface.

The comparator gets the same treatment because it decides pass and fail. Before
this decision `DIFF_OUTPUTS` inherited `container = null` and ran host `python3`
and host `bcftools`: `bcftools view --no-version -Ov` to decode BCF and
`bcftools index --stats` to check index companions. CI pinned
`bcftools=1.19-1build2`, an Ubuntu archive version string that is superseded and
dropped from the archive over time, and a local run used whatever was on `PATH`.
A bcftools change can flip a verdict without either implementation moving.

## Architecture

linux/amd64. The conda lock is entirely `linux-64`, so the image has no other
build, and Docker on an arm64 host selects the amd64 image and emulates it.
`verification/scripts/dump-conda-lock.sh` already passes
`--platform linux/amd64`; the harness that produces the observations did not, and
nothing recorded which substrate produced a number.

Parity is observed on amd64 Linux. CI runs `ubuntu-24.04`, which is native amd64,
so the authoritative runs already satisfy this. Local runs on an arm64 host are
development diagnostics, not parity evidence.

arm64 support in hap-rs is a property of hap-rs, not a parity surface. Which
platforms hap-rs supports is decided in
[Define supported execution environments and the reproducibility boundary](https://github.com/adamrtalbot/hap.py/issues/29).

## Environment settings

A setting that configures the interpreter or virtual machine legacy runs on, or
the order in which the harness schedules legacy tasks, is a harness constraint
and is allowed. A setting that changes the behaviour of a library the tool calls
is an override and is not.

Allowed under that line: `RTG_JAVA_OPTS=-Xint` on the legacy vcfeval path in
`verification/modules/happy.nf`, which puts the JVM in interpreted mode, and
`maxForks = 1` on `HAPPY_LEGACY`, which constrains scheduling.

Excluded: source patches, and the `sitecustomize.py` shim, which changed how
pandas answered som.py rather than how CPython ran it. The distinction is the
layer, not the mechanism: an environment variable that reaches the runtime is
fine, one that reaches a library is not.

Neither allowed setting is an exemption, so the register stays at two entries.
Every allowed setting is named in this ADR with its reason; one that is not named
here is not allowed.

## Host tooling

The host JVM runs Nextflow and nf-test. Legacy runs on the openjdk 11.0.8 inside
the image. The host JVM is not part of the baseline. CI installs Nextflow 26.04.6
and nf-test 0.9.5 with `nf-core/setup-nextflow` and `nf-core/setup-nf-test`, both
pinned to a commit SHA, and asserts both versions before the gate runs.

Verification lives in the same repository as the code, so the repository commit
is the harness revision and no separate record is kept.

## Reference format conversion

`--engine vcfeval` needs an RTG SDF. Legacy is given a frozen SDF bundle from the
samplesheet's `reference_sdf` column; hap-rs is given the FASTA and reads no SDF.
Removing that input is a hap-rs feature, so the asymmetry is permanent.

Each bundle records its sha256, the generating command, the RTG source tag and
the source commit in a README beside the fixture. No equivalence check between
bundle and FASTA is required: legacy fails if they disagree.

Five samplesheet rows run this path against legacy: `matrix_vcfeval_real`,
`matrix_vcfeval_paths`, `matrix_vcfeval_pass_only`, `matrix_vcfeval_custom_roc`,
and `matrix_vcfeval_bcf`. Two hap-rs contract cases exercise it with no legacy
side.

## Reference data

Downloaded once and used from a local copy. Provenance rules for public data are
decided in
[Set public-data admission and provenance rules](https://github.com/adamrtalbot/hap.py/issues/24).

## The change rule

One legacy identity is authoritative at a time. Moving a pin re-baselines:
legacy is re-run and its fresh output is the expectation. Nothing is preserved
and nothing is reconciled.

Nothing can be preserved, because nothing is stored. The harness keeps no legacy
baselines. `verification/results/` is gitignored with no tracked files, there are
no nf-test snapshots, and `verification/tests/main.nf.test` asserts only the row
count per lane and that every comparison in the run has an empty difference list.
Legacy and hap-rs run together in one pipeline execution and are diffed live.

hap-rs has published no version, so a pin move invalidates no claim. From 1.0 on,
a pin move that changes a covered observable is a versioned change.

### Applied to the pandas re-pin

Rebuild at pandas 0.20.3, re-run every comparison, keep the result. There is no
differential run against the old digest, and no requirement that the 102 rows
observed on the Wave image return byte-identical.
`0002-drop-in-claim-against-one-unpatched-legacy-image.md` set that gate on the
assumption that the harness keeps baselines. It keeps none.

The gate that remains is the ordinary one: hap-rs is frozen while legacy moves,
so a legacy change surfaces immediately as a failed comparison.

## Consequences

The comparator moves into its own container with its own lock.

`verification/containers/sompy-0.3.10.yml` is deleted. It is referenced nowhere,
pins hap.py 0.3.10 rather than 0.3.15, and names pandas 0.20.3 at the wrong tool
version, so it reads as authoritative while being unused.

The 102-row gate text in `verification/README.md` is replaced.

No caller-side environment is pinned. Once the quay image retires, one container
produces every legacy observation.

Both changes to `verification/` belong to
[test(verification): unify the legacy reference environment](https://github.com/adamrtalbot/hap.py/issues/38).
