---
title: Verification
description: Compare hap-rs output with pinned legacy tools.
---

Nextflow runs a pinned legacy tool and the local `hap` command for 158
samplesheet cases. It also runs two HAPPY contract checks against the Rust
command alone. nf-test collects all 160 results after the Rust unit and
integration tests pass.

## Coverage

The default matrix contains 158 paired comparisons and two Rust-only contract
checks:

| Lane | Legacy tool | hap-rs command | Checks |
|---|---|---|---:|
| `happy` | `hap.py` | `hap germline` | 34 paired + 2 Rust-only |
| `sompy` | `som.py` | `hap somatic` | 29 |
| `prepy` | `pre.py` | `hap pre` | 47 |
| `ftxpy` | `ftx.py` | `hap ftx` | 27 |
| `qfy` | `qfy.py` | `hap quantify` | 9 |
| `vcfcheck` | `vcfcheck` | `hap validate` | 12 |

Six `verification/assets/samplesheet.*.csv` files list the 158 paired inputs and
arguments. `verification/main.nf` defines the two additional Rust-only checks.

## Run the complete gate

The repository pins Nextflow 26.04.6 and nf-test 0.9.5.

```bash
cargo build --release --bin hap
export PATH="$(pwd)/target/release:$PATH"
cd verification
nf-test test --ci --coverage
```

CI runs all six lanes. Use a narrowed run to inspect one command.

## Run one lane

Use `HAP_TEST_CASES` to diagnose a discrepancy:

```bash
export PATH="$(pwd)/target/release:$PATH"
cd verification
HAP_TEST_CASES=sompy \
  nf-test test --ci tests/main.nf.test
```

Valid values are `happy`, `sompy`, `prepy`, `ftxpy`, `qfy`, and `vcfcheck`.

## Comparison contract

The comparator requires:

- identical artifact sets;
- exact ordered text and non-ROC CSV content after global exclusions;
- equal ROC CSV headers and row multisets, including duplicate counts;
- equal typed JSON trees after canonicalizing ROC table row order and generated
  table indexes;
- equal VCF records after the comparator removes runtime headers.

A missing or extra file fails the case. `comparison.json` records the lane,
case, artifact, and first location for each difference. nf-test reads the
combined `verification.json` to decide pass or failure.

The comparator excludes these runtime and provenance fields from all cases:

- version, timestamp, command line, and generated description fields;
- deprecated vcfeval template metadata ignored by the native engine;
- `sompyversion` and `sompycmd` CSV columns;
- runtime VCF headers such as source, date, and bcftools command/version.

## vcfeval references

Nextflow gives the legacy lane a committed RTG Tools 3.12.1 SDF bundle and the
Rust lane the matching FASTA. Contributors need Java and a supported container
runtime for verification; RTG runs inside the pinned container, and `hap` uses
neither RTG nor Java at runtime.

## Governing a discrepancy

1. Add the smallest input fixture needed to reproduce the behavior.
2. Add a samplesheet row to the affected named lane.
3. Run that lane and inspect its structured difference.
4. Fix the Rust implementation and add local regression coverage.
5. Run the complete matrix.

Keep exclusions, legacy output, and failing cases intact. Fix the Rust
implementation until the outputs match.
