---
title: Contributing
description: Add a verification case, fix the Rust implementation, and prove parity.
---

Use a verification case for each behavior change. The six samplesheets list
the inputs and arguments that `nf-test` compares with pinned hap.py tools.

## Requirements

| Tool | Required version or role |
|---|---|
| Rust | 1.89 or newer |
| Java | 17 |
| Nextflow | 26.04.6 |
| nf-test | 0.9.5 |
| Docker | Runs the pinned legacy tools |
| `bcftools`, `tabix` | Indexed-VCF verification queries |

## 1. Add your use case to a samplesheet

Choose the samplesheet for the command you are changing:

| Samplesheet | Rust command | Lane |
|---|---|---|
| `verification/assets/samplesheet.happy.csv` | `hap germline` | `happy` |
| `verification/assets/samplesheet.sompy.csv` | `hap somatic` | `sompy` |
| `verification/assets/samplesheet.prepy.csv` | `hap pre` | `prepy` |
| `verification/assets/samplesheet.ftxpy.csv` | `hap ftx` | `ftxpy` |
| `verification/assets/samplesheet.qfy.csv` | `hap quantify` | `qfy` |
| `verification/assets/samplesheet.vcfcheck.csv` | `hap validate` | `vcfcheck` |

Add the smallest row that reproduces the behavior. Put new inputs in
`verification/assets/fixtures/` and reference them with a
`local:assets/fixtures/...` path. Nextflow runs one comparison between the
legacy tool and `hap-rs` for each row.

Add and run the samplesheet case before editing Rust.

## 2. Run nf-test or Nextflow

Build the release binary, enter `verification/`, and run the affected lane. The
example below runs the germline lane:

```bash
cargo build --release --bin hap
export PATH="$(pwd)/target/release:$PATH"
cd verification
HAP_TEST_CASES=happy \
  nf-test test --ci tests/main.nf.test
```

Run Nextflow when you need process logs or work files:

```bash
export PATH="$(pwd)/target/release:$PATH"
cd verification
nextflow run main.nf \
  --cases happy --outdir results
```

Replace `happy` with `sompy`, `prepy`, `ftxpy`, `qfy`, or `vcfcheck`.
Use `nf-test` for the final check.

## 3. Identify the failure or gap

Confirm that the new case fails for the expected reason. Inspect the lane's
`comparison.json` and the `verification.json` that Nextflow writes to the output
directory. Each difference names the lane, case, artifact, and location.

Reduce the fixture and arguments to one behavior. Keep the pinned legacy tools,
comparison rules, exclusions, and existing cases unchanged.

## 4. Fix the Rust implementation

Make the smallest correction in `rust/src/`. Use Rust types and ownership to
model the behavior, return errors with context, and add a regression test for
the code you changed.

Run the Rust checks before returning to the verification pipeline:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --release --all-targets
```

## 5. Re-run the verification pipeline

Rebuild `hap` and repeat the affected nf-test lane:

```bash
cargo build --release --bin hap
export PATH="$(pwd)/target/release:$PATH"
cd verification
HAP_TEST_CASES=happy \
  nf-test test --ci tests/main.nf.test
```

Once the focused case passes, run the complete gate:

```bash
export PATH="$(pwd)/target/release:$PATH"
cd verification
nf-test test --ci --coverage
```

Before opening a pull request, run all six lanes. `nf-test` expects 158 unique
comparisons, each with `ok: true` and an empty `differences` list.

## Rules that protect the signal

- Add verification coverage before changing production behavior.
- Use focused lanes for diagnosis and the complete matrix for the final check.
- Fix product differences in Rust. Change comparison normalization only for a
  proven representation difference, with positive and negative regression cases.
- Keep failing comparisons until the Rust output matches.
- Keep the pinned legacy images unchanged.

See the [verification guide](../verification/) for the lane counts, comparison
rules, output files, and fixture provenance. The root
[`CONTRIBUTING.md`](https://github.com/adamrtalbot/hap.py/blob/master/CONTRIBUTING.md)
is the checkout version of this workflow.
