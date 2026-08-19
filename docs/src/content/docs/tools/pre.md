---
title: Preprocess
description: Complete reference for hap pre.
---

`hap pre` filters, normalizes, decomposes, and rewrites a variant callset. Use
it to prepare inputs before a comparison.

Aliases: `hap preprocess`, `hap prepy`; successor to `pre.py`.

## Usage

```text
hap pre [OPTIONS] <INPUT> <OUTPUT>
```

| Argument | Required | Description |
|---|---:|---|
| `<INPUT>` | Yes | Source VCF, compressed VCF, or BCF. |
| `<OUTPUT>` | Yes | Destination VCF/BCF path. |

## Example

```bash
hap pre input.vcf.gz normalized.vcf.gz \
  --reference reference.fa \
  --pass-only \
  --leftshift \
  --decompose
```

## Reference, regions, and filters

| Option | Description |
|---|---|
| `-r, --reference <FASTA>` | Reference FASTA. `hap` reads no environment variable in its place. |
| `-l, --location <REGION>` | Limit processing to a genomic location. |
| `-R, --restrict-regions <BED>` | Restrict records to intervals. |
| `-T, --target-regions <BED>` | Apply target intervals. |
| `--pass-only` | Discard filtered records. |
| `--filters-only <LIST>` | Retain records with selected filter values and discard the rest. |
| `--filter-nonref` | Filter terminal non-reference records; default: on. |

## Normalization and transformation

| Option | Description |
|---|---|
| `--fixchr` | Repair systematic `chr` prefix differences. |
| `--no-fixchr` | Disable contig-prefix repair. |
| `--somatic` | Apply established somatic genotype conversion. |
| `--set-gt <MODE>` | Set `half`, `hemi`, `het`, `hom`, or `first` genotype mode. |
| `--convert-gvcf-to-vcf` | Convert eligible gVCF records to VCF records. |
| `--bcf` | Write BCF output. |
| `--bcftools-norm` | Use bcftools-compatible normalization behavior. |
| `-L, --leftshift` | Left-shift eligible indels; default: on. |
| `--no-leftshift` | Disable left shifting. |
| `--decompose` | Split eligible complex and multiallelic records; default: on. |
| `-D, --no-decompose` | Disable decomposition. |
| `--gender <VALUE>` | Ploidy handling: `male`, `female`, `auto` (default), or `none`. |
| `-w, --window-size <N>` | Normalization window size; default `10000`. |

## Execution

| Option | Description |
|---|---|
| `--threads <N>` | Requested worker count. |
| `--logfile <PATH>` | Write command logging to a file. |
| `--verbose` | Increase logging; conflicts with `--quiet`. |
| `--quiet` | Reduce logging; conflicts with `--verbose`. |
| `--force-interactive` | Preserve the compatibility switch for interactive execution. |
| `-v, --version` | Print the legacy subcommand version response. |
| `-h, --help` | Print command help. |

## Output

The output extension selects VCF, compressed VCF, or BCF. `hap` writes a Tabix
or CSI index beside indexed output and replaces the data and index as a pair.
