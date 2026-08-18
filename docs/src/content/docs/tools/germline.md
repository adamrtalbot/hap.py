---
title: Germline Comparison
description: Complete reference for hap germline.
---

`hap germline` compares a truth VCF with a query VCF, classifies matching
variants, and produces summary, detailed, structured, and ROC reports.

Aliases: `hap compare`; successor to the `hap.py` entry point.

## Usage

```text
hap germline [OPTIONS] --report-prefix <REPORT_PREFIX> <TRUTH> <QUERY>
```

| Argument | Required | Description |
|---|---:|---|
| `<TRUTH>` | Yes | Trusted VCF, compressed VCF, or BCF. |
| `<QUERY>` | Yes | Callset to evaluate. |
| `-o, --report-prefix <PATH>` | Yes | Prefix applied to every report file. The parent directory must exist. |
| `-r, --reference <FASTA>` | For most runs | Reference FASTA. Without this option, `hap` checks `HG19`, then `HGREF`. |

## Example

```bash
hap germline truth.vcf.gz query.vcf.gz \
  --reference reference.fa \
  --false-positives confident.bed \
  --stratification regions.tsv \
  --engine xcmp \
  --write-vcf \
  --report-prefix results/sample
```

## Comparison engine

| Option | Description |
|---|---|
| `--engine <ENGINE>` | Comparison engine: `xcmp` (default), `vcfeval`, `scmp-somatic`, or `scmp-distance`. |
| `--unhappy` | Disable haplotype comparison; alias: `--no-haplotype-comparison`. |
| `-w, --window-size <N>` | Haplotype comparison window; default `50`. |
| `--xcmp-enumeration-threshold <N>` | Maximum xcmp path enumeration threshold; default `16768`. |
| `--xcmp-expand-hapblocks <N>` | xcmp haplotype-block expansion; default `30`. |
| `--scmp-distance <N>` | Distance for `scmp-distance`; default `30`. Alias: `--lose-match-distance`. |
| `--engine-vcfeval-path <PATH>` | Wrapper compatibility option. Ignored by the native engine; use `--engine vcfeval --reference <FASTA>`. Warns on stderr and stays supported. |
| `--engine-vcfeval-template <PATH>` | Wrapper compatibility option. Ignored by the native engine; use `--engine vcfeval --reference <FASTA>`. Warns on stderr and stays supported. |

## Preprocessing

| Option | Description |
|---|---|
| `--pass-only` | Discard filtered query records. Truth preprocessing remains pass-only. |
| `--preprocess-truth` | Enable the optional truth preprocessing control. |
| `--convert-gvcf-truth` | Convert truth gVCF records before comparison. |
| `--convert-gvcf-query` | Convert query gVCF records before comparison. |
| `--convert-gvcf-to-vcf` | Enable gVCF-to-VCF conversion for both sides. |
| `--usefiltered-truth` | Retain filtered truth records. |
| `--filters-only <LIST>` | Retain the named filter values and discard the rest. |
| `--preprocessing-window-size <N>` | Preprocessing window size; default `10000`. |
| `--adjust-conf-regions` | Adjust confident regions around normalized variants; default: on. |
| `--no-adjust-conf-regions` | Disable confident-region adjustment. |
| `-L, --leftshift` | Enable left shifting. |
| `--no-leftshift` | Disable left shifting. |
| `--decompose` | Enable decomposition of complex records. |
| `-D, --no-decompose` | Disable decomposition. |
| `--bcftools-norm` | Use bcftools-compatible normalization behavior. |
| `--fixchr` | Add or remove `chr` prefixes when input/reference naming differs. |
| `--no-fixchr` | Disable automatic contig-prefix repair. |
| `--filter-nonref` | Filter `<NON_REF>` and related non-reference records. |
| `--somatic` | Convert genotypes using the established somatic mode. |
| `--set-gt <MODE>` | Set genotype mode: `half`, `hemi`, `het`, `hom`, or `first`. |
| `--gender <VALUE>` | Ploidy handling: `male`, `female`, `auto` (default), or `none`. |
| `--bcf` | Request BCF intermediate/output handling. |

## Regions and stratification

A stratification table can load several BED files in one run. `hap` resolves
paths from the TSV file:

```tsv title="regions.tsv"
SUBSTITUTIONS	substitutions.bed
INDELS	indels.bed
```

The separator is a tab. Pass the table once with `--stratification regions.tsv`;
the extended report contains a subset for every table row.

| Option | Description |
|---|---|
| `-R, --restrict-regions <BED>` | Restrict comparison to non-overlapping intervals. |
| `-T, --target-regions <BED>` | Limit records to target intervals; accepts overlapping targets. |
| `-f, --false-positives <BED>` | Confident regions that distinguish false positives from unknown calls. |
| `-l, --location <REGION>` | Limit work to a genomic location expression. |
| `--stratification <TSV>` | Map stratification labels to BED files. |
| `--stratification-region <LABEL>` | Select a stratification label; repeat for multiple labels. |
| `--stratification-fixchr` | Reconcile `chr` prefixes in stratification inputs. |

## Reports and ROC output

| Option | Description |
|---|---|
| `-V, --write-vcf` | Write the annotated comparison VCF. |
| `-X, --write-counts` | Write detailed count reports; default: on. |
| `--no-write-counts` | Suppress the extended count report. |
| `--output-vtc` | Include variant-type count output. |
| `--preserve-info` | Preserve source INFO fields in comparison output. |
| `--roc <FIELD>` | Score ROC output with this field; default `QUAL`. |
| `--no-roc` | Disable ROC tables. |
| `--roc-regions <LABEL>` | Select a stratification label for ROC output; repeatable. |
| `--roc-filter <FILTER>` | Limit ROC observations by filter. |
| `--roc-delta <VALUE>` | ROC threshold spacing; default `0.5`. |
| `--ci-alpha <VALUE>` | Confidence-interval alpha; default `0` (off). |
| `--no-json` | Suppress structured metrics JSON. |

## Execution and compatibility

| Option | Description |
|---|---|
| `--threads <N>` | Requested worker count. |
| `--scratch-prefix <PATH>` | Parent or prefix for temporary comparison files. |
| `--keep-scratch` | Retain temporary files after completion. |
| `--force-interactive` | Preserve the compatibility switch for interactive execution. |
| `--logfile <PATH>` | Write command logging to a file. |
| `--verbose` | Increase logging; conflicts with `--quiet`. |
| `--quiet` | Reduce logging; conflicts with `--verbose`. |
| `-t, --type <TYPE>` | Compatibility annotation type: `xcmp` or `ga4gh`; the engine selects the effective type. |
| `-v, --version` | Print the legacy subcommand version response. |
| `-h, --help` | Print command help. |

## Outputs

The default prefix produces:

- `.summary.csv`
- `.extended.csv`
- `.runinfo.json`
- `.metrics.json.gz` unless `--no-json`
- `.roc.*.csv.gz` unless `--no-roc`
- `.vcf.gz` or `.bcf` with `--write-vcf`

See [Metrics](../../reference/metrics/) for interpretation and
[Inputs & Outputs](../../reference/inputs-outputs/) for format requirements.
