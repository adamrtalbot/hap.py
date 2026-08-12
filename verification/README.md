# Verification

Read the [website guide](https://adamrtalbot.github.io/hap.py/project/verification/)
for an overview. Use this file for work in `verification/`.

This directory contains the parity gate for the Rust `hap` migration. For every
samplesheet row, Nextflow runs the named legacy hap.py tool in a pinned Wave
image and the matching `hap` subcommand from `PATH`, then compares the emitted
`result*` files. nf-test fails when a comparison contains a difference outside
the metadata exclusions below.

## Legacy reference image

`containers/happy-0.3.15.yml` defines the legacy hap.py environment. Build it
with Wave CLI and freeze the result:

```bash
wave --conda-file verification/containers/happy-0.3.15.yml \
  --platform linux/amd64 --freeze --await --output json
```

Set the frozen Wave image name in `nextflow.config`; do not use an expiring
request-scoped Wave URL as a parity reference. The current reference is
`community.wave.seqera.io/library/happy-0.3.15:41c2102638513597`.

The SOMPY and FTXPY lanes use the digest-pinned hap.py 0.3.15 image from the
nf-core/variantbenchmarking v1.6.0dev full-size run. Its exact dependency set
affects legacy CSV float serialization, so it must not be replaced by another
hap.py 0.3.15 build without first observing and updating the legacy contract.
The SOMPY legacy process supplies a Python environment shim for pandas' removed
`display.height` option, which som.py sets only while rendering ambiguity
explanations. The shim ignores that display-only setting and does not alter
som.py or its output data.

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
cd verification
PATH="../target/release:$PATH" nf-test test --ci tests/main.nf.test
```

Use `HAP_TEST_CASES` to diagnose one lane:

```bash
PATH="../target/release:$PATH" HAP_TEST_CASES=sompy nf-test test --ci tests/main.nf.test
```

CI runs all six lanes. Use a narrowed run to inspect one command.

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
- Keep legacy tools unmodified inside their pinned images.
- Put product behavior in Rust and artifact comparison in Nextflow.
- Fix the Rust implementation when a case differs. Keep the case and comparison
  rule intact.
- Add a samplesheet row and run it before editing Rust.

Keep each third-party license notice beside its imported fixture. The gate uses
those fixtures as test inputs.
