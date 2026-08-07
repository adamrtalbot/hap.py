---
title: Core Concepts
description: Truth, query, regions, engines, and report prefixes.
---

## Truth and query

The **truth** callset supplies the expected variants. `hap` compares the
**query** callset with it.

- `hap` labels an unmatched truth call as a false negative (`FN`).
- `hap` labels an unmatched query call as a false positive (`FP`).
- `hap` labels a match as a true positive (`TP`).

Germline comparison considers haplotypes and genotypes. Somatic comparison uses
alleles because callers can encode somatic genotypes in different ways.

## Regions

Several region controls serve different purposes:

| Control | Meaning |
|---|---|
| Confident/false-positive BED (`-f`) | Defines where unmatched query calls count toward precision. |
| Restrict regions (`-R`) | Selects records by overlap; xcmp requires non-overlapping intervals. |
| Target regions (`-T`) | Selects target intervals and can represent overlapping targets. |
| Location (`-l`) | Limits work to a genomic location expression. |

All region files must use the same assembly and contig naming as the callsets.

## Stratification

A stratification TSV maps a label to a BED file, one pair per line:

```text
LOW_COMPLEXITY  low-complexity.bed
SEGMENTAL_DUP   segmental-duplications.bed
```

`hap` resolves paths from the TSV location. Use `--stratification-region` to select
specific labels and `--stratification-fixchr` when the stratification contigs
need `chr` prefix reconciliation.

## Comparison engines

`hap germline` supports four engines:

- `xcmp`: default haplotype-aware comparison.
- `vcfeval`: native vcfeval-style path matching with the FASTA passed through
  `--reference`.
- `scmp-somatic`: loose somatic allele matching.
- `scmp-distance`: distance-bounded allele matching.

See [Commands & Engines](../../reference/commands-engines/) for selection
guidance.

## Report prefixes

Commands that create a report family take a prefix. With
`--report-prefix results/sample`, `hap` writes
`results/sample.summary.csv`, `results/sample.extended.csv`, and so on.

The parent directory must already exist for germline and quantify reports.
