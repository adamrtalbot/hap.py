# verification/

Pinned legacy-to-Rust parity gate for the pure-Rust `hap` executable.

Each samplesheet row runs twice:

- **Legacy:** the corresponding hap.py tool in the immutable Wave container.
- **Rust:** the matching `hap` subcommand built from this checkout.

The row supplies the same inputs and arguments to both sides. A staged Rust
verifier compares their output trees and writes one `status.json`; nf-test
checks every unique selected row has a published result and requires every
status to pass.

## Governed matrix

| Lane | Legacy tool | Rust command | Default rows |
|---|---|---|---:|
| happy | `hap.py` | `hap germline` | 19 |
| sompy | `som.py` | `hap somatic` | 23 |
| prepy | `pre.py` | `hap pre` | 37 |
| ftxpy | `ftx.py` | `hap ftx` | 24 |
| qfy | `qfy.py` | `hap quantify` | 4 |
| vcfcheck | `vcfcheck` | `hap validate` | 4 |
| **Total** | | | **111** |

The qfy lane generates bounded chr21 xcmp-annotated VCFs from pinned
truth/query/reference fixtures, then gives the same generated VCF and index to
both quantifiers. The default `cases` setting selects all six lanes;
`HAP_TEST_CASES` narrows a diagnostic run without changing the release matrix
definition. `OPTION_MATRIX.md` defines the constrained-pairwise selection
method, while `assets/option-coverage.csv` gives every samplesheet row a unique
coverage ID and records the options exercised.

## Layout

```text
verification/
├── main.nf                       workflow entrypoint
├── nextflow.config               pins and default matrix parameters
├── nf-test.config
├── assets/
│   ├── samplesheet.happy.csv
│   ├── samplesheet.sompy.csv
│   ├── samplesheet.prepy.csv
│   ├── samplesheet.ftxpy.csv
│   ├── samplesheet.qfy.csv
│   ├── samplesheet.vcfcheck.csv
│   ├── option-coverage.csv
│   └── fixtures/                    compact option-matrix fixtures
├── modules/                      legacy/Rust lane processes and reports
├── scripts/
│   ├── verify-ftx-options.sh
│   └── verify-somatic-options.sh
├── OPTION_MATRIX.md              matrix design and comparison contract
└── tests/main.nf.test            coverage and parity assertions
```

`rust/src/parity_verifier.rs` owns the file-type-aware comparator used by the
workflow. `rust/src/verify_fixtures.rs` is the feature-gated helper binary
staged into comparison tasks.

## Prerequisites

- Rust 1.89 or newer (the crate's declared MSRV in `Cargo.toml`).
- Docker with access to the pinned legacy image.
- Nextflow 26.04.6, enforced by `nextflow.config`.
- nf-test 0.9.5.
- `tabix` for independent TBI/VCF-CSI retrieval checks.

Install the pinned nf-test release with:

```bash
curl -fsSL https://get.nf-test.com | bash -s 0.9.5
```

These are evaluation dependencies, not product runtime dependencies. The
shipped `hap` executable is pure Rust and does not require Python, C, or C++.

## Run the gate

From the repository root, build both binaries with the verifier feature:

```bash
cargo build --release --features verification --bin hap --bin verify-fixtures
```

Then run the clean default matrix from this directory:

```bash
cd verification
HAP_BIN="$PWD/../target/release/hap" \
  VERIFY_BIN="$PWD/../target/release/verify-fixtures" \
  NXF_SYNTAX_PARSER=v2 NXF_VER=26.04.6 nf-test test tests/main.nf.test --ci
```

For a focused or resumable diagnostic:

```bash
HAP_TEST_CASES=happy,prepy \
  NF_TEST_NEXTFLOW_OPTIONS=-resume \
  NXF_SYNTAX_PARSER=v2 NXF_VER=26.04.6 nf-test test tests/main.nf.test
```

To inspect workflow output without the nf-test assertion layer, run from the
repository root:

```bash
NXF_SYNTAX_PARSER=v2 NXF_VER=26.04.6 \
  nextflow run verification/main.nf --cases qfy,vcfcheck
```

The workflow intentionally exits after producing `report.md` and `report.csv`.
The release decision comes from nf-test, which checks the row coverage contract
and each authoritative comparator status.

## Direct option oracles

The matrix is supplemented by two direct, pinned-container scripts:

```bash
cargo build --release --bin hap
HAP_BIN="$PWD/target/release/hap" verification/scripts/verify-ftx-options.sh
HAP_BIN="$PWD/target/release/hap" verification/scripts/verify-somatic-options.sh
```

- `verify-ftx-options.sh` runs 10 cases and compares each emitted feature table
  byte-for-byte. It additionally requires `bcftools` to construct its BCF input.
- `verify-somatic-options.sh` runs 3 cases and compares the symmetric output set,
  canonicalizing only the metrics timestamp. It additionally requires Python 3
  for that verification-only JSON canonicalization.

Neither script adds a runtime dependency to `hap`.

## Samplesheets

Every CSV under `assets/` is append-only. A row becomes one legacy task, one
Rust task, and one comparison status. `sample_id` is the join key and must be
unique within its lane.

Paths are resolved against `params.fixture_base`, which defaults to a pinned
upstream commit. Override it with an absolute local checkout for development:

```bash
NXF_SYNTAX_PARSER=v2 NXF_VER=26.04.6 \
  nextflow run verification/main.nf --fixture_base "$PWD"
```

The free-form `args` column is inserted before positional inputs on both
invocations. Add a new row with a distinct `sample_id` to govern another flag
combination; quote shell-sensitive values carefully.

## Outputs

```text
results/
├── report.md
├── report.csv
└── <lane>/<sample_id>/
    ├── legacy/result*
    ├── rust/result*
    ├── diff.log
    └── status.json
```

The nf-test contract rejects missing or extra samples, duplicate IDs, absent
artifact trees, stale/incomplete status payloads, empty comparisons, and any
file verdict whose `ok` value is false.

## Comparison rules

| Pattern | Rule |
|---|---|
| `*.summary.csv`, `*.extended.csv`, `*.stats.csv`, `*.features.csv` | Compare decoded text in emitted order. |
| `*.roc.*.csv.gz`, `*.csv.gz` | Decompress, then compare text in emitted order. |
| `*.vcf`, `*.vcf.gz` | Compare records and nonvolatile headers in emitted order; only the explicit producer/time/command/reference allowlist is excluded. |
| `*.metrics.json*`, `*.runinfo.json*` | Compare recursively; exclude only governed runtime identity, environment, timestamp, and representational metadata. |
| VCF `*.tbi`, `*.csi` | `tabix` must retrieve every companion VCF record in source order. |
| BCF `*.csi` | The built-in CSI reader follows indexed BGZF chunks and must retrieve every BCF record in source order. |
| `*.fai` | Byte exact. |
| Everything else | Byte exact. |

Artifact comparison is symmetric and preserves ordering and duplicates. A file
present on only one side fails, as does an empty matching artifact set.

## Saved and real-world fixture helper

The separate helper remains available only with the `verification` feature:

```bash
cargo run --features verification --bin verify-fixtures -- check-rust
cargo run --features verification --bin verify-fixtures -- check-realworld
```

`check-realworld` is a manual network/container diagnostic. The samplesheet
workflow and nf-test assertions above define the release gate.
