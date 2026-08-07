---
title: Introduction
description: What hap-rs does and which command to use.
---

`hap-rs` is a native toolkit for benchmarking and preparing genomic variant
callsets. One `hap` executable provides comparison, preprocessing,
quantification, feature extraction, and validation workflows.

`hap-rs` reimplements hap.py in Rust and keeps its command contracts, report
formats, and benchmarking terms.

## Choose a command

| Goal | Command |
|---|---|
| Compare germline truth and query callsets | [`hap germline`](../../tools/germline/) |
| Benchmark somatic calls without requiring genotype agreement | [`hap somatic`](../../tools/somatic/) |
| Normalize or transform one callset | [`hap pre`](../../tools/pre/) |
| Extract somatic caller features | [`hap ftx`](../../tools/ftx/) |
| Calculate reports from an annotated comparison VCF | [`hap quantify`](../../tools/quantify/) |
| Detect malformed or reference-inconsistent records | [`hap validate`](../../tools/validate/) |

## Rust implementation

- Haplotype enumeration and exact matching (`xcmp`)
- A FASTA-backed `vcfeval` comparison engine
- Somatic and distance-based allele matching
- VCF preprocessing and report generation
- VCF, BGZF, BCF, Tabix, and CSI I/O

You need the `hap` executable and the input files referenced by each command.
Contributors use Java and pinned legacy tools to run the
[verification gate](../../project/verification/).

## Compatibility

The following aliases help existing workflows move to the single executable:

| Previous entry point | hap-rs command |
|---|---|
| `hap.py` | `hap germline` |
| `som.py` | `hap somatic` |
| `pre.py` | `hap pre` |
| `ftx.py` | `hap ftx` |
| `qfy.py` | `hap quantify` or `hap qfy` |
| `vcfcheck` | `hap validate` or `hap vcfcheck` |

[Install hap-rs](../installation/) or run the [Quick Start](../quick-start/).
