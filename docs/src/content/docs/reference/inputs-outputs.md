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

Pass the FASTA with `--reference`. `hap` ships no reference data. GRCh37, hg19
and b37 work the same way as any other assembly: supply the FASTA and its `.fai`
companion as arguments.

## Region files

- BED files use zero-based, half-open intervals.
- xcmp rejects overlapping restrict-region BED intervals.
- Stratification TSVs map a label to a BED path.
- Contig names must match; the command's `fixchr` option can reconcile prefixes.

## Report families

| Suffix | Produced by | Meaning |
|---|---|---|
| `.summary.csv` | germline, quantify; somatic with `-P`, a feature table and `--happy-stats` | Headline aggregate metrics. |
| `.extended.csv` | germline; quantify unless `--no-write-counts` | Stratified and subtype metrics. Not observed from somatic, including with `--happy-stats` and a typed feature table. |
| `.stats.csv` | somatic | Somatic counts and performance metrics. |
| `.metrics.json.gz` | germline, quantify | Structured report data, compressed. |
| `.metrics.json` | somatic | Structured report data, uncompressed. |
| `.runinfo.json` | germline | Effective options and provenance. |
| `.roc.all.csv.gz` | germline, quantify | Written even with `--no-roc`. |
| `.roc.Locations.*.csv.gz` | germline, quantify; somatic per `--roc` profile | Threshold-specific performance data. |
| `.vcf.gz` / `.bcf` | germline by default; quantify with `--write-vcf`; pre as its output path | Annotated or normalized records. |
| `.vcf.gz.tbi` / `.bcf.csi` | beside every indexed data file above | Index companion, extension follows the data extension. |
| `.features.csv` | somatic with `--feature-table` | Extracted feature table. |
| `.csv` | ftx | Extracted feature table. |
| JSON path from `-o` | validate | Validation summary and diagnostics. A bare `hap validate` writes no file. |
| BED path from `--errors-bed` | validate | Record intervals with validation errors. An hap-rs extension; legacy `vcfcheck` has no such option. |

Each [tool page](../../tools/germline/) marks default outputs and the flags that
enable other files. What the comparison against legacy observes is every file
whose name begins with the output prefix an invocation names, and nothing else;
see `docs/adr/0007-observe-the-output-prefix-set.md`.
