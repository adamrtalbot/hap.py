# Verification

Read the [website guide](https://adamrtalbot.github.io/hap.py/project/verification/)
for an overview. Use this file for work in `verification/`.

This directory contains the parity gate for the Rust `hap` migration. For every
samplesheet row, Nextflow runs the named legacy hap.py tool in one of two pinned
container environments and the matching `hap` subcommand from `PATH`, then
compares the emitted `result*` files. nf-test fails when a comparison contains a
difference outside the metadata exclusions below.

## Legacy reference environment

One reference environment governs every lane, and legacy runs there unmodified:
no interpreter shim, no source patch, no environment override. An observation
taken through a patch is not evidence of legacy behavior, so a legacy tool that
will not run means the environment is wrong. The tool version is not the
reference; `hap.py 0.3.15` names more than one behavior. See
`docs/adr/0002-drop-in-claim-against-one-unpatched-legacy-image.md`.

The reference is
`community.wave.seqera.io/library/happy-0.3.15:41c2102638513597`, the frozen Wave
build of `containers/happy-0.3.15.yml`. It carries all six legacy tools plus
`rtg`, `bcftools`, `samtools`, and `java`. Because the digest is opaque,
`containers/happy-0.3.15.conda-lock.txt` records the 108 packages it installs,
including the `libstdcxx` and `openjdk` builds the recipe never names.
Regenerate it with `scripts/dump-legacy-conda-lock.sh`. The recipe lists what was
requested; the lock lists what was installed, and the lock is what explains a
number.

The SOMPY and FTXPY lanes have not moved onto it yet. They still run
`quay.io/biocontainers/hap.py@sha256:d63b963a6cb01b4830393b22369e7b91d298e4156dde353739e74e4cfa4f96d0`,
whose pandas 0.24.2 writes CSV floats as shortest-repr where the reference
image's pandas 0.19.2 truncates to twelve significant digits — `1/3` as
`0.3333333333333333` against `0.333333333333`. That image also dropped pandas'
`display.height`, which som.py sets while rendering ambiguity explanations, so
the SOMPY process injects a `sitecustomize.py` shim to swallow it. The second
image and the shim are both known violations of the rule above, kept only until
those lanes are repointed and their legacy CSVs re-observed at pandas 0.19.2. The
reference image needs no shim: pandas removed `display.height` in 0.20, and the
reference pins 0.19.2.

`HAPPY_LEGACY` runs with `maxForks = 1`. Serialization is a scheduling
constraint on the harness, not a modification of legacy, so it is compatible with
the rule above.

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
  deprecated vcfeval template argument that the Rust engine ignores.
- CSV columns named `sompyversion` and `sompycmd`.
- VCF runtime headers such as source, date, and bcftools command/version.

Each exclusion applies to all cases. All remaining content must match.
nf-test reads the combined `verification.json` to decide pass or failure. Each
lane also writes a `comparison.json` with its differences.

## Rules for changes

- Keep the named lane processes visible.
- Keep the reference image immutable and pinned by digest, never a mutable tag.
- Never patch legacy to make it run. No interpreter shim, no source patch, no
  environment override. A legacy tool that will not run means the environment is
  wrong.
- Put product behavior in Rust and artifact comparison in Nextflow.
- Fix the Rust implementation when a case differs. Keep the case and comparison
  rule intact.
- Add a samplesheet row and run it before editing Rust.

Keep each third-party license notice beside its imported fixture. The gate uses
those fixtures as test inputs.
