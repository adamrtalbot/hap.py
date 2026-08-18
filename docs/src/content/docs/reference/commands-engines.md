---
title: Commands & Engines
description: Top-level command aliases and germline comparison engines.
---

## Command index

| Canonical command | Accepted aliases | Complete options |
|---|---|---|
| `hap germline` | `hap compare` | [Germline](../../tools/germline/) |
| `hap somatic` | None | [Somatic](../../tools/somatic/) |
| `hap pre` | `hap preprocess`, `hap prepy` | [Preprocess](../../tools/pre/) |
| `hap ftx` | `hap ftxpy` | [Feature extraction](../../tools/ftx/) |
| `hap quantify` | `hap qfy` | [Quantify](../../tools/quantify/) |
| `hap validate` | `hap vcfcheck` | [Validate](../../tools/validate/) |

Global options are `-h, --help` and `-V, --version`. Run
`hap <command> --help` to inspect the installed version's interface.

## Engine selection

| Engine | Best suited to | Reference input | Notes |
|---|---|---|---|
| `xcmp` | General germline benchmarking | FASTA + `.fai` | Default; evaluates compatible haplotype paths. |
| `vcfeval` | vcfeval-style germline comparison | FASTA + `.fai` | Native implementation; ignores legacy SDF template inputs. |
| `scmp-somatic` | Somatic-like allele comparison through germline reports | FASTA + `.fai` | Uses somatic matching semantics. |
| `scmp-distance` | Calls allowed to match within a distance | FASTA + `.fai` | Configure distance with `--scmp-distance`. |

Engine choice changes comparison semantics. All engines use the same report
prefix and core quantification tables.

## Wrapper compatibility

`hap` accepts several flags that existing wrappers pass. It ignores
`--engine-vcfeval-path` and `--engine-vcfeval-template` because the Rust engine
runs in process with a FASTA. Both options warn on stderr and stay supported, so
wrappers that pass them keep working; new wrappers should supply
`--reference <FASTA>` instead. See the [germline page](../../tools/germline/) and
[compatibility policy](../../project/compatibility/).
