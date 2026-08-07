---
title: Quick Start
description: Run a first germline benchmark and inspect its reports.
---

The command below compares a query callset with a truth set inside confident
regions.

## Prepare inputs

You need:

- `truth.vcf.gz`: trusted variants
- `query.vcf.gz`: variants to evaluate
- `reference.fa` and `reference.fa.fai`: reference assembly and index
- `confident.bed`: regions where the truth set is callable

All inputs must use the same reference assembly and compatible contig names.

## Run the comparison

```bash
mkdir -p results

hap germline truth.vcf.gz query.vcf.gz \
  --reference reference.fa \
  --false-positives confident.bed \
  --report-prefix results/sample
```

Pass the confident-region BED through `--false-positives`. The option controls
the precision denominator.

## Read the results

The default run creates:

| File | Contents |
|---|---|
| `sample.summary.csv` | Headline SNP and INDEL counts, recall, and precision. |
| `sample.extended.csv` | Metrics split by type, subtype, filter, genotype, and region. |
| `sample.runinfo.json` | Effective command configuration and run metadata. |
| `sample.metrics.json.gz` | Structured forms of the report tables. |
| `sample.roc.*.csv.gz` | Thresholded ROC tables when you enable ROC output. |

Add `--write-vcf` to retain the annotated comparison VCF.

## Next steps

- Use [`--engine vcfeval`](../../tools/germline/#comparison-engine) for a
  vcfeval-style comparison backed by `reference.fa`.
- Add a [stratification TSV](../concepts/#stratification) to measure difficult
  region sets.
- Read the [metrics reference](../../reference/metrics/) before comparing
  headline values across datasets.
