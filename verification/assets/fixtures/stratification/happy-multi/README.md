# Multiple stratification regions

`regions.tsv` demonstrates the legacy hap.py stratification-table format: each
tab-separated row gives a report label and a BED path relative to the table.
This example loads two BED files in one run, producing `SUBSTITUTIONS` and
`INDELS` subset rows in the extended report.
