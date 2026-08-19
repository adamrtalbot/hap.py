---
title: Somatic Comparison
description: Complete reference for hap somatic.
---

`hap somatic` compares somatic truth and query callsets by allele and location,
then reports true positives, false positives, false negatives, unknown calls,
and ambiguous calls.

Successor to the `som.py` entry point.

## Usage

```text
hap somatic [OPTIONS] --output <OUTPUT> <TRUTH> <QUERY>
```

| Argument | Required | Description |
|---|---:|---|
| `<TRUTH>` | Yes | Trusted somatic callset. |
| `<QUERY>` | Yes | Somatic callset to evaluate. |
| `-o, --output <PREFIX>` | Yes | Output prefix for stats and metrics files. |
| `-r, --reference <FASTA>` | For normalization | Reference FASTA. Needed by `--normalize-*` and by an automatic FP denominator; plain allele comparison runs without one. `hap` reads no environment variable in its place. |

## Example

```bash
hap somatic truth.vcf.gz query.vcf.gz \
  --reference reference.fa \
  --false-positives confident.bed \
  --include-nonpass \
  --output results/somatic
```

## Regions and classification

| Option | Description |
|---|---|
| `-l, --location <REGION>` | Limit comparison to a genomic location. |
| `-R, --restrict-regions <BED>` | Restrict records to intervals. |
| `-T, --target-regions <BED>` | Apply target intervals. |
| `-f, --false-positives <BED>` | Define confident regions for FP versus unknown classification. |
| `-a, --ambiguous <BED>` | Add an ambiguous-region BED; repeatable. |
| `--ambi-fp` | Count ambiguous query calls as false positives. |
| `--no-ambi-fp` | Do not count ambiguous calls as false positives. |
| `--count-unk` | Include unknown calls in reported counts. |
| `--no-count-unk` | Exclude unknown calls from the selected count behavior. |
| `-e, --explain_ambiguous` | Add ambiguity class and reason detail to metrics output. |
| `--fp-region-size <N>` | Override the false-positive region size in rate calculations. |

## Filtering and normalization

| Option | Description |
|---|---|
| `-P, --include-nonpass` | Include filtered records. |
| `--normalize-truth` | Normalize truth alleles. |
| `--normalize-query` | Normalize query alleles. |
| `-N, --normalize-all` | Normalize both truth and query. |
| `--fixchr-truth` | Repair truth contig prefixes. Alias: `--fix-chr-truth`. |
| `--fixchr-query` | Repair query contig prefixes. Alias: `--fix-chr-query`. |
| `--no-fixchr-truth` | Disable truth contig-prefix repair. |
| `--no-fixchr-query` | Disable query contig-prefix repair. |
| `--no-order-check` | Disable input record-order checking. |

## Feature and ROC reports

| Option | Description |
|---|---|
| `--feature-table <NAME>` | Emit a supported generic or caller-specific feature table. |
| `--happy-stats` | Also write hap.py-style summary output; requires `--include-nonpass` and `--feature-table`. |
| `--bam <PATH>` | Add a BAM for feature extraction; repeatable. |
| `--roc <MODE>` | Select a supported caller-specific ROC mode. |
| `--bin-afs` | Stratify by allele frequency; requires `--feature-table`. |
| `--af-binsize <VALUE>` | Allele-frequency bin width; default `0.2`. |
| `--af-truth <FIELD>` | Truth allele-frequency field; default `I.T_ALT_RATE`. |
| `--af-query <FIELD>` | Query allele-frequency field; default `T_AF`. |
| `--count-filtered-fn` | Count filtered false negatives; requires `--include-nonpass` and `--feature-table`. |
| `--ci-level <VALUE>` | Confidence level; default `0.95`. |

Feature tables include `generic`, Strelka, Mutect, VarScan2, and Pisces
SNV/INDEL profiles from the `admix` and `hcc` workflows. `hap` rejects unknown
profile names before comparison.

## Execution

| Option | Description |
|---|---|
| `--scratch-prefix <PATH>` | Parent or prefix for temporary files. |
| `--keep-scratch` | Keep temporary files after completion. |
| `--continue` | Preserve established continuation behavior. |
| `--logfile <PATH>` | Write command logging to a file. |
| `--verbose` | Increase logging; conflicts with `--quiet`. |
| `--quiet` | Reduce logging; conflicts with `--verbose`. |
| `-h, --help` | Print command help. |

## Outputs

- `<prefix>.stats.csv`: primary somatic counts and rates
- `<prefix>.metrics.json`: structured metrics, uncompressed, unlike the
  `.metrics.json.gz` that `hap germline` and `hap quantify` write
- `<prefix>.features.csv`: with `--feature-table <name>`
- `<prefix>.summary.csv`: hap.py-style reports, which need `-P` and a feature
  table alongside `--happy-stats`; `--happy-stats` on its own exits 1
- caller-specific ROC output when `--roc` selects one of its supported profiles
