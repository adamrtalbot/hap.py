---
title: Inputs & Outputs
description: Input formats and report files for hap-rs.
---

## Variant files

`hap-rs` reads VCF, BGZF-compressed VCF, and BCF. Indexed queries use the index
next to the data file:

| Data | Index |
|---|---|
| `sample.vcf.gz` | `sample.vcf.gz.tbi` |
| `sample.bcf` | `sample.bcf.csi` |

Input callsets should declare contigs and sample columns in the VCF header.
Truth, query, reference, and region files must describe the same assembly.

## References

Commands that compare alleles against a reference accept FASTA. Place a
samtools-compatible `.fai` index next to it:

```text
reference.fa
reference.fa.fai
```

The Rust `vcfeval` engine reads this FASTA. Compatibility flags accept an RTG
SDF path but ignore its value.

## Region files

- BED files use zero-based, half-open intervals.
- xcmp rejects overlapping restrict-region BED intervals.
- Stratification TSVs map a label to a BED path.
- Contig names must match; the command's `fixchr` option can reconcile prefixes.

## Report families

| Suffix | Produced by | Meaning |
|---|---|---|
| `.summary.csv` | germline, quantify; optional somatic | Headline aggregate metrics. |
| `.extended.csv` | germline, quantify; optional somatic | Stratified and subtype metrics. |
| `.stats.csv` | somatic | Somatic counts and performance metrics. |
| `.metrics.json` / `.metrics.json.gz` | comparison tools | Structured report data. |
| `.runinfo.json` | germline | Effective options and provenance. |
| `.roc.*.csv.gz` | germline, quantify, somatic | Threshold-specific performance data. |
| `.vcf.gz` / `.bcf` | optional comparison or preprocessing output | Annotated or normalized records. |
| `.csv` | ftx | Extracted feature table. |
| JSON path from `-o` | validate | Validation summary and diagnostics. |
| BED path from `-e` | validate | Record intervals with validation errors. |

Each [tool page](../../tools/germline/) marks default outputs and the flags that
enable other files.
