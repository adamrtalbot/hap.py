# hap-rs

Pure-Rust implementation of hap.py, distributed as one `hap` executable. The
binary provides `germline`, `somatic`, `pre`, `ftx`, `quantify`, and `validate`
subcommands, including native BCF and CSI handling. It has no Python, C, or C++
runtime dependency.

## Build and test

```bash
cargo build
cargo test
```

Run `cargo run -- --help` for the command list and per-command help.

## Saved-fixture verification

The verification commands in this section are repository-development
workflows. The `verification/` harness is intentionally not included in the
published Cargo package; use a source checkout to run them. The packaged
`verification` feature includes only the `verify-fixtures` helper and its saved
fixtures.

The `verify-fixtures` helper is excluded from normal product builds. Enable the
`verification` feature to build or run it.

Refresh expected outputs from the immutable legacy container:

```bash
cargo run --features verification --bin verify-fixtures -- bless-legacy \
  --image community.wave.seqera.io/library/hap.py_rtg-tools@sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6
```

Check that the container still reproduces those saved outputs:

```bash
cargo run --features verification --bin verify-fixtures -- check-legacy \
  --image community.wave.seqera.io/library/hap.py_rtg-tools@sha256:1562a301d591cb42e4465465b579a68fc84429b8ed7c7cbda7f728f2e0a57be6
```

Check the Rust implementation against the same fixtures:

```bash
cargo run --features verification --bin verify-fixtures -- check-rust
```

## Release parity gate

This gate is available from a source checkout, not from the published Cargo
package.

`verification/` runs each governed samplesheet row through both the pinned
legacy oracle and the local `hap` binary. The default matrix contains 63 rows:

| Lane | Rust command | Rows |
|---|---|---:|
| happy | `hap germline` | 12 |
| sompy | `hap somatic` | 12 |
| prepy | `hap pre` | 26 |
| ftxpy | `hap ftx` | 11 |
| qfy | `hap quantify` | 1 |
| vcfcheck | `hap validate` | 1 |

The gate checks complete samplesheet publication, symmetric legacy/Rust
artifact sets, and every per-case comparator status. Comparisons preserve
ordering and duplicates: stable text is exact, JSON is recursive with a narrow
runtime-metadata allowlist, and VCF/BCF index checks retrieve source records in
order through TBI or CSI. The workflow produces reports; nf-test turns failed
comparison statuses or incomplete coverage into a failed gate.

The evaluation toolchain is pinned to Nextflow 26.04.6 and nf-test 0.9.5:

```bash
cargo build --release --features verification --bin hap --bin verify-fixtures
cd verification
HAP_BIN="$PWD/../target/release/hap" \
  VERIFY_BIN="$PWD/../target/release/verify-fixtures" \
  NXF_SYNTAX_PARSER=v2 NXF_VER=26.04.6 nf-test test --ci
```

Two focused scripts exercise option combinations directly against the same
pinned oracle image:

- `verification/scripts/verify-ftx-options.sh`: 10 FTX cases, including BCF
  input, normalization, region controls, caller tables, and BAM depth inputs.
- `verification/scripts/verify-somatic-options.sh`: 3 somatic cases covering
  explanation output, false-positive regions, and normalization/filtering.

See `verification/README.md` for prerequisites, samplesheet layout, comparison
rules, focused runs, and direct-oracle commands. `rust/FLAG_PARITY.md` is the
detailed flag-level compatibility contract.

## Repository layout

- `Cargo.toml`: crate manifest and product/verification feature boundary
- `rust/src/`: Rust implementation
- `tests/`: integration tests and static fixtures
- `verification/`: pinned legacy-to-Rust release gate
