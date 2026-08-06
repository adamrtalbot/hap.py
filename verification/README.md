# Verification

This directory is the parity gate for the Rust `hap` migration.  For every
samplesheet row, Nextflow runs the named legacy hap.py tool in a pinned Wave
image and the matching `hap` subcommand from `PATH`, then compares the emitted
`result*` files.  nf-test fails when a comparison contains a meaningful
difference.

## Lanes

| Lane | Legacy command | Rust command |
|---|---|---|
| HAPPY | `hap.py` | `hap germline` |
| SOMPY | `som.py` | `hap somatic` |
| PREPY | `pre.py` | `hap pre` |
| FTXPY | `ftx.py` | `hap ftx` |
| QFY | `qfy.py` | `hap quantify` |
| VCFCHECK | `vcfcheck` | `hap validate` |

The six `assets/samplesheet.*.csv` files are the executable test matrix. Each
row supplies inputs and arguments for one case. Add a row to reproduce and
govern a newly discovered discrepancy; do not change a row, comparison rule,
or expected output to conceal one.

## Run

Build the Rust binary and make it available on `PATH`. The full matrix also
requires RTG Tools 3.12.1-1 for the real `--engine vcfeval` case; for example:

```bash
micromamba create -n hap-parity -c conda-forge -c bioconda rtg-tools=3.12.1=hdfd78af_1
micromamba activate hap-parity
```

Then run the gate:

```bash
cargo build --release
cd verification
PATH="../target/release:$PATH" nf-test test --ci tests/main.nf.test
```

To diagnose one lane locally, use `HAP_TEST_CASES`:

```bash
PATH="../target/release:$PATH" HAP_TEST_CASES=sompy nf-test test --ci tests/main.nf.test
```

CI must run all six lanes. A narrowed run is diagnostic only.

## Comparison rules

The comparator is embedded in `modules/diff.nf`. It requires equal artifact
sets, compares ordered text and CSV content, recursively compares typed JSON,
and records a structured difference containing the lane, case, artifact, and
location. The only ignored fields are global runtime/provenance metadata:

- JSON version, timestamp, command-line, and generated description fields.
- CSV columns named `sompyversion` and `sompycmd`.

These exclusions are global rather than case-specific. All remaining content
must match exactly. The resulting `verification.json` is the authoritative
machine-readable gate result; lane-level `comparison.json` files contain the
individual differences.

## Rules for changes

- Keep named lane processes visible; do not replace them with generic runners.
- Keep legacy tools unmodified inside their pinned images.
- Keep Rust-only behavior in Rust; comparison behavior belongs in Nextflow.
- Fix the Rust implementation when a case differs. Never bless, weaken, or
  delete a case merely to make the gate pass.
- Add a samplesheet row and let nf-test pinpoint the discrepancy before adding
  a regression fix.

Third-party fixture license notices remain beside their imported files. The
fixtures themselves are part of the test inputs and must not be removed during
documentation cleanup.
