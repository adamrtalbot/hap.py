# Rust architecture

The Rust crate is a private layered implementation behind the single public
entry point `hap_rs::run`. Modules may depend on another module in the same
layer or on a layer to their left, never to their right:

```text
domain <- adapters <- engines <- application <- crate root
   ^           ^          ^            ^
   +-----------+----------+---- cli_compat
```

`cli_compat` is an input boundary. It owns Clap and the legacy argument
surface; application commands consume those private request types. The crate
root performs process-level parsing and dispatch.

`tests/architecture.rs` parses the crate with `syn`, resolves grouped,
absolute, re-exported, and relative paths to canonical module owners, rejects
outward dependencies, and detects cycles between modules in every layer.

| Layer | Directory | Ownership |
|---|---|---|
| Domain | `rust/src/domain` | Filesystem-free variant/interval models, benchmark rows, count models, and mathematical/statistical primitives |
| Adapters | `rust/src/adapters` | VCF/BCF/FASTA parsing and publication, report and metrics serialization |
| Engines | `rust/src/engines` | In-memory alignment, normalization, comparison, ROC calculation, and caller-specific transforms |
| Application | `rust/src/application` | Command validation, workflow orchestration, scratch lifecycle, adapter/engine wiring, and output publication |
| CLI compatibility | `rust/src/cli_compat` | Clap definitions, aliases, legacy exit behavior, and compatibility argument normalization |

There are no compatibility re-exports at the crate root: the only public crate
item is `run`. Application modules may collaborate with each other, but domain,
adapters, engines, and CLI compatibility code cannot import application
modules. Filesystem entry points for SCMP, vcfeval, and ROC publication live in
the application layer and pass neutral in-memory records to engines.

## Cohesive command and engine modules

The large command implementations remain orchestration façades and delegate
filesystem-free responsibilities to child modules:

- `compare/{matching,rows,metrics,output,genotype}.rs`: matching, row
  construction, aggregation, publication preparation, and genotype rules.
- `somatic/{features,normalization,metrics,reports,allele_frequency}.rs`:
  feature extraction, normalization, aggregation, report rendering, and AF
  stratification.
- `preprocess/{alleles,blocksplit,canonical,normalization,options,genotype}.rs`:
  allele projection, scheduling, header/field canonicalization, left shifting,
  option policy, and genotype encoding.
- `quantify/{annotations,counting,reporting,stratification,regions}.rs`:
  GA4GH annotation, counting, table rendering, subset assignment, and interval
  arithmetic.
- `roc/{accumulation,contributions,legacy,rendering}.rs`: observation
  accumulation, axis expansion, legacy ordering, and in-memory CSV rendering.

Focused unit tests use strings, numbers, or in-memory domain records; command
regression suites are isolated in each façade's `test_suite.rs`. VCF/BCF/FASTA
adapters and application publication modules retain filesystem integration
coverage.
