---
title: Metrics
description: Interpret hap-rs counts and rates.
---

## Primary counts

| Metric | Meaning |
|---|---|
| `TRUTH.TOTAL` | Truth variants considered in the row. |
| `QUERY.TOTAL` | Query variants considered in the row. |
| `TRUTH.TP` | Truth calls that the query matches. |
| `TRUTH.FN` | Truth calls with no query match. |
| `QUERY.TP` | Query calls that match truth. |
| `QUERY.FP` | Query calls with no truth match. |
| `QUERY.UNK` | Query calls outside assessable/confident regions. |

Truth and query true-positive counts can differ when complex records decompose
into different primitive observations.

## Rates

Recall measures the fraction of truth calls that the query matched:

```text
recall = TRUTH.TP / (TRUTH.TP + TRUTH.FN)
```

Precision measures the fraction of query calls within assessment that match
truth:

```text
precision = QUERY.TP / (QUERY.TP + QUERY.FP)
```

`METRIC.Frac_NA` is the fraction of query observations outside assessment,
often because they fall beyond confident regions.

## Rows and subsets

Summary reports contain SNP and INDEL aggregates. Extended reports add
dimensions such as:

- variant subtype and size
- genotype
- filter status
- genomic stratification region
- quality threshold for ROC output

Compare rows with the same subset labels. `*`, `PASS`, and named
difficult-region subsets use different denominators.

## Confidence intervals

`--ci-alpha` on germline/quantify and `--ci-level` on somatic control confidence
interval output. Their conventions differ because each command preserves its
established interface; consult the command page before changing the default.
