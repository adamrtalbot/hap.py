---
title: Validate
description: Complete reference for hap validate.
---

`hap validate` checks VCF headers, samples, records, ploidy, allele forms,
ordering, and BCF translation. Pass a FASTA to check reference consistency.

Alias: `hap vcfcheck`; successor to `vcfcheck`.

## Usage

```text
hap validate [OPTIONS] [INPUT]
```

Provide a positional input or use `--input-file`. The command accepts one form
per invocation.

## Example

```bash
hap validate calls.vcf.gz \
  --reference reference.fa \
  --output-json results/check.json \
  --errors-bed results/errors.bed \
  --check-bcf-errors true \
  --all-warnings true
```

## Inputs and outputs

| Option | Description |
|---|---|
| `[INPUT]` | Positional VCF/BCF input. Required if you omit `--input-file`. |
| `--input-file <INPUT>` | Named form of the input argument. |
| `-r, --reference <FASTA>` | Validate coordinates and REF alleles against a FASTA. |
| `-o, --output-json <PATH>` | Write the count summary as JSON. Alias: `--output-file`. |
| `-e, --errors-bed <PATH>` | Write reference mismatches and invalid spans as BED records. |

Without `--output-json` or `--errors-bed`, `hap` writes diagnostics to standard
error.

## Selection

| Option | Description |
|---|---|
| `-l, --location <REGION>` | Limit validation to a genomic location. |
| `-R, --regions <BED>` | Restrict records to intervals. |
| `-T, --targets <BED>` | Apply target intervals. |
| `-f, --apply-filters <BOOL>` | When `true`, discard filtered records; default `false`. |
| `--limit-records <N>` | Stop after N records. |
| `--message-every <N>` | Emit progress every N processed records. |

## Validation controls

| Option | Description |
|---|---|
| `-H, --strict-homref <BOOL>` | Enable strict homozygous-reference checks; default `false`. |
| `--check-bcf-errors <BOOL>` | Fail on records that cannot translate to BCF; default `false`. |
| `-W, --all-warnings <BOOL>` | Emit every supported warning class; default `false`. |
| `-h, --help` | Print command help. |

Boolean options require an explicit `true` or `false` value.

## Summary fields

The JSON summary includes record, reference/non-reference, ploidy, inferred
male, overlap, reference-padding, symbolic-ALT, and uncertain-length counts.
The error BED adds a reason such as `OUT_OF_RANGE`, `REF_MISMATCH`, or
`POLYPLOID_GT` to each affected interval.
