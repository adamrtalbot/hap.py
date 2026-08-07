---
title: Quantify
description: Complete reference for hap quantify.
---

`hap quantify` turns an annotated comparison VCF into summary, extended,
structured, and ROC reports. Use it when comparison and quantification are
separate pipeline stages.

Alias: `hap qfy`; successor to `qfy.py`.

## Usage

```text
hap quantify [OPTIONS] --report-prefix <PREFIX> --reference <FASTA> <INPUT_VCF>
```

| Argument | Required | Description |
|---|---:|---|
| `<INPUT_VCF>` | Yes | Indexed annotated comparison VCF or BCF. |
| `-o, --report-prefix <PREFIX>` | Yes | Prefix for report files. The parent directory must exist. |
| `-r, --reference <FASTA>` | Yes | Reference FASTA for types and region sizes. |

## Example

```bash
hap quantify comparison.vcf.gz \
  --reference reference.fa \
  --false-positives confident.bed \
  --stratification regions.tsv \
  --report-prefix results/sample
```

## Annotation and regions

| Option | Description |
|---|---|
| `-t, --type <TYPE>` | Input annotation contract: `xcmp` or `ga4gh`; xcmp is the effective default. |
| `-f, --false-positives <BED>` | Confident regions for FP versus unknown classification. |
| `--stratification <TSV>` | Map stratification labels to BED files. |
| `--stratification-region <LABEL>` | Select a stratification label; repeatable. |
| `--stratification-fixchr` | Reconcile `chr` prefixes in stratification inputs. |
| `--adjust-conf-regions <VALUE>` | Preserve the established confident-region adjustment argument. |

## Reports and ROC output

| Option | Description |
|---|---|
| `-V, --write-vcf` | Write the region-annotated VCF or BCF. |
| `-X, --write-counts` | Write detailed counts; default: on. |
| `--no-write-counts` | Suppress the extended count report. |
| `--output-vtc` | Include variant-type count output. |
| `--preserve-info` | Preserve source INFO and region-extent annotations. |
| `--bcf` | Write annotated variant output as BCF. |
| `--roc <FIELD>` | Score ROC output with this field; default `QUAL`. |
| `--no-roc` | Disable ROC tables. |
| `--roc-regions <LABEL>` | Select a stratification label for ROC output; repeatable. |
| `--roc-filter <FILTER>` | Limit ROC observations by filter. |
| `--roc-delta <VALUE>` | ROC threshold spacing; default `0.5`. |
| `--ci-alpha <VALUE>` | Confidence-interval alpha; default `0.0` (off). |
| `--no-json` | Suppress `metrics.json.gz`. |

## Execution and compatibility

| Option | Description |
|---|---|
| `--threads <N>` | Requested worker count. |
| `--logfile <PATH>` | Write command logging to a file. |
| `--verbose` | Increase logging; conflicts with `--quiet`. |
| `--quiet` | Reduce logging; conflicts with `--verbose`. |
| `--force-interactive` | Preserve the compatibility switch for interactive execution. |
| `-v, --version` | Print the legacy subcommand version response after argument validation. |
| `-h, --help` | Print command help. |

## Outputs

- `<prefix>.summary.csv`
- `<prefix>.extended.csv` unless `--no-write-counts`
- `<prefix>.metrics.json.gz` unless `--no-json`
- `<prefix>.roc.*.csv.gz` when the input has benchmark samples and ROC output
  is on
- `<prefix>.vcf.gz` or `<prefix>.bcf` with `--write-vcf`

The input and output variant paths must differ. Quantification requires an index
for compressed VCF or BCF input.
