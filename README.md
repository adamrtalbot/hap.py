# hap-rs

Pure-Rust implementation of hap.py, distributed as one `hap` executable. The
binary provides `germline`, `somatic`, `pre`, `ftx`, `quantify`, and `validate`
subcommands, including native BCF and CSI handling. Its native engines have no
Python, C, or C++ runtime dependency. The optional `--engine vcfeval` mode
delegates to RTG Tools, matching hap.py's external-engine contract.

## Build and test

```bash
cargo build
cargo test
```

Run `cargo run -- --help` for the command list and per-command help.

## Release parity gate

This gate is available from a source checkout, not from the published Cargo
package.

`verification/` runs each governed samplesheet row through both the pinned
legacy reference and the local `hap` binary. The default matrix contains 147 rows:

| Lane | Rust command | Rows |
|---|---|---:|
| happy | `hap germline` | 26 |
| sompy | `hap somatic` | 28 |
| prepy | `hap pre` | 45 |
| ftxpy | `hap ftx` | 27 |
| qfy | `hap quantify` | 9 |
| vcfcheck | `hap validate` | 12 |

The gate compares the legacy and Rust artifact sets directly. It requires exact
ordered content after only global runtime/provenance metadata is removed; a
shared missing or extra artifact also fails. nf-test turns any structured
difference into a failed gate.

The evaluation toolchain is pinned to Nextflow 26.04.6, nf-test 0.9.5, and RTG
Tools 3.12.1-1. RTG is required by the real vcfeval parity case:

```bash
cargo build --release
cd verification
PATH="$PWD/../target/release:$PATH" nf-test test --ci tests/main.nf.test
```

See `verification/README.md` for the test matrix, comparison rules, and focused
diagnostic runs.

## Repository layout

- `Cargo.toml`: crate manifest and product/verification feature boundary
- `rust/src/`: Rust implementation
- `tests/`: integration tests and static fixtures
- `verification/`: pinned legacy-to-Rust release gate
