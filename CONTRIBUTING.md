# Contributing to hap-rs

Use a verification case for each behavior change. The six samplesheets in
[`verification/assets/`](verification/assets/) list the inputs and arguments
that `nf-test` compares with pinned hap.py tools.

## Requirements

- Rust 1.89 or newer
- Java 17
- Nextflow 26.04.6
- nf-test 0.9.5
- Docker
- `bcftools` and `tabix`

## Development workflow

### 1. Add your use case to a samplesheet

Choose the samplesheet for the affected command:

| Samplesheet | Rust command |
|---|---|
| `verification/assets/samplesheet.happy.csv` | `hap germline` |
| `verification/assets/samplesheet.sompy.csv` | `hap somatic` |
| `verification/assets/samplesheet.prepy.csv` | `hap pre` |
| `verification/assets/samplesheet.ftxpy.csv` | `hap ftx` |
| `verification/assets/samplesheet.qfy.csv` | `hap quantify` |
| `verification/assets/samplesheet.vcfcheck.csv` | `hap validate` |

Add the smallest row that reproduces the behavior. Put new inputs in
`verification/assets/fixtures/` and reference them with a
`local:assets/fixtures/...` path. Nextflow runs one legacy-versus-Rust
comparison for each row.

Add and run the samplesheet case before editing Rust.

### 2. Run the verification case

Build `hap`, move into `verification/`, and run the affected lane. For example:

```bash
cargo build --release --bin hap
export PATH="$(pwd)/target/release:$PATH"
cd verification
HAP_TEST_CASES=happy \
  nf-test test --ci tests/main.nf.test
```

The lane names are `happy`, `sompy`, `prepy`, `ftxpy`, `qfy`, and `vcfcheck`.

Run Nextflow with the same lane name when you need process logs or work files:

```bash
export PATH="$(pwd)/target/release:$PATH"
cd verification
nextflow run main.nf \
  --cases happy --outdir results
```

### 3. Identify the failure or gap

Confirm that the new case fails for the expected reason. Inspect the lane's
`comparison.json` and the `verification.json` that Nextflow writes to the output
directory. Each difference names the lane, case, artifact, and location.

Reduce the case to one behavior. Keep the pinned legacy tools, comparison
rules, exclusions, and existing cases unchanged.

### 4. Fix the Rust implementation

Make the smallest correction in `rust/src/`. Use Rust types and ownership to
model the behavior, return errors with context, and add a regression test for
the code you changed.

Before returning to the verification pipeline, run the Rust checks:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test --release --all-targets
```

### 5. Re-run verification

Rebuild the release binary and re-run the affected lane:

```bash
cargo build --release --bin hap
export PATH="$(pwd)/target/release:$PATH"
cd verification
HAP_TEST_CASES=happy \
  nf-test test --ci tests/main.nf.test
```

Run the complete gate after the focused case passes:

```bash
export PATH="$(pwd)/target/release:$PATH"
cd verification
nf-test test --ci --coverage
```

Before opening a pull request, run all six lanes. `nf-test` expects 158 unique
comparisons, each with `ok: true` and an empty `differences` list.

Dependency and attribution policy can be checked locally with:

```bash
cargo audit --deny warnings
cargo deny check
python3 scripts/check-compression-backends.py
python3 scripts/check-notices.py
```

## Verification rules

- Add verification coverage before changing production behavior.
- Use focused lanes for diagnosis and the complete matrix for the final check.
- Fix product differences in Rust. Change comparison normalization only for a
  proven representation difference, with positive and negative regression cases.
- Keep failing comparisons until the Rust output matches.
- Keep the pinned legacy images unchanged.

## Compatibility changes

Read the [compatibility policy](docs/src/content/docs/project/compatibility.md)
before preserving or removing a legacy quirk. New emulation needs all of:

- pinned parity or captured upstream regression evidence;
- a nearby explanatory source-and-rationale comment;
- an inventory entry with a governance class and removal decision; and
- separately named `legacy_only_` and `normative_` tests where both contracts
  exist.

User-facing deprecations must warn on stderr and document the replacement and
removal release. Removing governed behavior follows the versioning and release
note process in the compatibility policy.

See [`verification/README.md`](verification/README.md) for the harness layout,
comparison rules, and fixture details.
