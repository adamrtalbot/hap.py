---
title: Feature Extraction
description: Complete reference for hap ftx.
---

`hap ftx` extracts generic or caller-specific features from a VCF and writes a
CSV feature table. Pass a BAM to add depth features.

Alias: `hap ftxpy`; successor to `ftx.py`.

## Usage

```text
hap ftx [OPTIONS] --output <OUTPUT> <INPUT>
```

| Argument | Required | Description |
|---|---:|---|
| `<INPUT>` | Yes | Source callset. |
| `-o, --output <PATH>` | Yes | CSV destination. `hap` changes another extension to `.csv`. |

## Example

```bash
hap ftx calls.vcf.gz \
  --feature-table generic \
  --reference reference.fa \
  --output results/features.csv
```

## Options

| Option | Description |
|---|---|
| `-l, --location <REGION>` | Limit extraction to a genomic location. |
| `-R, --restrict-regions <BED>` | Restrict records to intervals. |
| `-T, --target-regions <BED>` | Apply target intervals. |
| `-P, --include-nonpass` | Include filtered records. |
| `--feature-table <NAME>` | Feature profile; default `generic`. |
| `--feature-label <LABEL>` | Add a label to the feature output. |
| `--bam <PATH>` | Add a BAM for depth-based features; repeatable. |
| `-r, --reference <FASTA>` | Reference FASTA for normalization or BAM features. |
| `--normalize` | Normalize records before feature extraction. |
| `--fix-chr` | Repair systematic `chr` prefix differences. |
| `-h, --help` | Print command help. |

## Feature profiles

The generic profile copies common variant annotations. Caller profiles support
the established Strelka SNV/INDEL layouts and generic legacy-caller parsing for
Mutect, VarScan2, and Pisces data.

## Output

The CSV contains one row per variant observation in the output. The feature
profile and BAM input determine the columns.
