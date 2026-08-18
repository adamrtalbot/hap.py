# Verification

Read the [website guide](https://adamrtalbot.github.io/hap.py/project/verification/)
for an overview. Use this file for work in `verification/`.

This directory contains the parity gate for the Rust `hap` migration. For every
samplesheet row, Nextflow runs the named legacy hap.py tool in the pinned
reference container and the matching `hap` subcommand from `PATH`, then
compares the emitted `result*` files in the pinned comparator container. nf-test fails when a comparison contains a
difference outside the metadata exclusions below.

## Legacy reference environment

One reference environment governs every lane, and legacy runs there unmodified:
no interpreter shim, no source patch, no environment override. An observation
taken through a patch is not evidence of legacy behavior, so a legacy tool that
will not run means the environment is wrong. The tool version is not the
reference; `hap.py 0.3.15` names more than one behavior. See
`docs/adr/0002-drop-in-claim-against-one-unpatched-legacy-image.md`.

The reference is
`community.wave.seqera.io/library/happy-0.3.15:2c2b5746d6b0da37`, the frozen Wave
build of `containers/happy-0.3.15.yml`. It carries all six legacy tools plus
`rtg`, `bcftools`, `samtools`, and `java`. Because the digest is opaque,
`containers/happy-0.3.15.conda-lock.txt` records the 108 packages it installs,
including the `libstdcxx` and `openjdk` builds the recipe never names.
Regenerate it with `scripts/dump-conda-lock.sh <image>`. The recipe lists what
was requested; the lock lists what was installed, and the lock is what explains a
number.

It pins pandas 0.20.3, which is the version window where one image runs all six
tools unpatched. Measured on the image: index-level `groupby` works,
`display.height` survives as a deprecation so som.py needs no shim, and `to_csv`
renders every float as Python 2 `str()` does, identically to the 0.19.2 build it
replaces. numpy 1.12.1, scipy 1.2.1, and libstdcxx 16.1.0 did not move. The
rebuild did move 17 peripheral packages the earlier solve had resolved
differently; the lock is the record of which.

The comparator is pinned the same way, because it decides pass and fail.
`DIFF_OUTPUTS` runs
`community.wave.seqera.io/library/hap-comparator:3fd84dea2dc8fda1`, the frozen
Wave build of `containers/hap-comparator.yml`, with its contents in
`containers/hap-comparator.conda-lock.txt`. It carries python 3.12.3 and the same
`bcftools` 1.17 build the reference image carries, so a BCF decode difference
cannot be an htslib-version artifact. It previously inherited `container = null`
and ran host `python3` and host `bcftools`, so a host change could flip a verdict
without either implementation moving.

Moving a pin re-baselines: legacy is re-run on the new image and its fresh output
is the expectation. Nothing is preserved, because nothing is stored: `results/`
is gitignored, there are no nf-test snapshots, and the harness runs legacy and
`hap` together in one execution and diffs them live.

The gate that remains is the ordinary one. `hap` is frozen while legacy moves, so
a legacy change surfaces immediately as a failed comparison. See
`docs/adr/0004-pin-the-legacy-baseline-to-one-container-identity.md`.

The four `--bam` rows are out of the truth set: `ftx_bam_depth`, `ftx_multi_bam`,
`somatic_bam_depth`, `matrix_multi_bam`. The original never pinned pandas, so it
permits installations where its own `--bam` paths cannot run, and any observation
of them reflects the pandas version this project picks rather than legacy
behavior. `--bam` is classed **no legacy reference** in the compatibility policy,
which makes hap-rs behavior there normative. Their fixtures stay under
`assets/fixtures/ftx-bam/` for the `normative_` expectations, which are separate
work. The matrix is 154 six-lane comparisons plus two hap-rs vcfeval contract
cases with no legacy side.

The line that rule draws is the layer a setting reaches. Configure the
interpreter or VM legacy runs on, or the order the harness schedules it in, and
that is a harness constraint. Change the behaviour of a library the tool calls
and that is an override. `maxForks = 1` on `HAPPY_LEGACY` constrains scheduling;
`RTG_JAVA_OPTS=-Xint` on the legacy vcfeval path in `modules/happy.nf` puts the
JVM in interpreted mode. The `sitecustomize.py` shim changed how pandas answered
som.py, so it fails the same test. Allowed settings are named in
`docs/adr/0004-pin-the-legacy-baseline-to-one-container-identity.md`; a setting
not named there is not allowed.

Parity is observed on linux/amd64. The conda lock is entirely `linux-64`, so the
image has no other build, and Docker on an arm64 host emulates it. CI runs
`ubuntu-24.04` natively, so authoritative runs already satisfy this. A local run
on arm64 is a development diagnostic, not parity evidence.

## Lanes

| Lane | Legacy command | Rust command |
|---|---|---|
| HAPPY | `hap.py` | `hap germline` |
| SOMPY | `som.py` | `hap somatic` |
| PREPY | `pre.py` | `hap pre` |
| FTXPY | `ftx.py` | `hap ftx` |
| QFY | `qfy.py` | `hap quantify` |
| VCFCHECK | `vcfcheck` | `hap validate` |

The six `assets/samplesheet.*.csv` files list the test matrix. Each row supplies
inputs and arguments for one case. Add a row to reproduce a
discrepancy. Keep existing rows, comparison rules, and expected outputs
unchanged.

## Run

Build the Rust binary, make it available on `PATH`, then run the gate:

```bash
cargo build --release
export PATH="$(pwd)/target/release:$PATH"
cd verification
nf-test test --ci tests/main.nf.test
```

Use `HAP_TEST_CASES` to diagnose one lane:

```bash
export PATH="$(pwd)/target/release:$PATH"
cd verification
HAP_TEST_CASES=sompy nf-test test --ci tests/main.nf.test
```

CI runs all six lanes. Use a narrowed run to inspect one command.

## Run the public germline discovery workflow

Public-data cases stay outside nf-test. Build `hap`, then run the intact public
germline samplesheet directly with Nextflow. The release-binary directory on
`PATH` must be absolute because Nextflow tasks run from their own work
directories:

```bash
cargo build --release --bin hap
export PATH="$(pwd)/target/release:$PATH"
cd verification
nextflow run main.nf \
  --cases happy \
  --happy_samplesheet "$PWD/assets/samplesheet.happy.public.csv" \
  --outdir "$PWD/results/public-germline"
```

Accept the run only when `results/public-germline/verification.json` reports
`comparison_count: 3` and every comparison has `ok: true` with an empty
`differences` list. The three comparisons are the intact
`germline_hg001_platinum` public-data case and the two existing HAPPY contract
cases. Also inspect
`results/public-germline/happy/germline_hg001_platinum/comparison.json`; it must
independently report `ok: true`, an empty `differences` list, and identical
legacy and hap-rs artifact inventories.

`--engine vcfeval` uses separate reference formats for the two implementations.
Nextflow gives the pinned legacy image an RTG Tools 3.12.1 SDF `.tar.gz` bundle
from the samplesheet's `reference_sdf` column. It gives `hap` the corresponding
FASTA. The product runs without RTG or Java. A README beside each fixture
records the bundle checksum and provenance.

## Comparison rules

`modules/diff.nf` contains the comparator. It requires equal artifact sets,
compares ordered text and non-ROC CSV content, compares ROC CSV rows as
multisets with duplicate counts, and compares typed JSON trees after
canonicalizing ROC table row order and generated table indexes. The comparator
records the lane, case, artifact, and location for each difference. It ignores
these runtime and provenance fields:

- JSON version, timestamp, command-line, generated description fields, and the
  `engine_vcfeval_template` argument the native engine ignores.
- CSV columns named `sompyversion` and `sompycmd`.
- VCF runtime headers such as source, date, and bcftools command/version.

Each exclusion applies to all cases. All remaining content must match.
nf-test reads the combined `verification.json` to decide pass or failure. Each
lane also writes a `comparison.json` with its differences.

## Rules for changes

- Keep the named lane processes visible.
- Keep the reference image immutable and pinned by digest, never a mutable tag.
- Never patch legacy to make it run. No interpreter shim, no source patch, no
  library option override. A legacy tool that will not run means the environment
  is wrong. Settings that change how the same code executes or when it is
  scheduled are harness constraints; they are allowed only when named in
  `docs/adr/0004-pin-the-legacy-baseline-to-one-container-identity.md`.
- Put product behavior in Rust and artifact comparison in Nextflow.
- Fix the Rust implementation when a case differs. Keep the case and comparison
  rule intact.
- Add a samplesheet row and run it before editing Rust.

Keep each third-party license notice beside its imported fixture. The gate uses
those fixtures as test inputs.
