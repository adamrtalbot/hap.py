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

Build the Rust binary, make it available on `PATH`, then run the gate:

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

`--engine vcfeval` uses two deliberately different reference inputs in this
gate. The pinned legacy image receives a committed RTG Tools 3.12.1 SDF
`.tar.gz` bundle from the samplesheet's `reference_sdf` column. The Rust
process receives only the corresponding FASTA. No host RTG or Java runtime is
used by the product implementation. Bundle checksums and generation provenance
are recorded beside each fixture.

## Comparison rules

The comparator is embedded in `modules/diff.nf`. It requires equal artifact
sets, compares ordered text and CSV content, recursively compares typed JSON,
and records a structured difference containing the lane, case, artifact, and
location. The only ignored fields are global runtime/provenance metadata:

- JSON version, timestamp, command-line, generated description fields, and the
  deprecated vcfeval template argument that native Rust intentionally ignores.
- CSV columns named `sompyversion` and `sompycmd`.
- VCF runtime headers such as source, date, and bcftools command/version.

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
