# EBI VCF Validator conformance fixtures

These files are retained byte-for-byte from the official
`EBIvariation/vcf-validator` test corpus at immutable commit
`0aacc4a44430ab9cee87d9925aacb28a2fb0a9fb` (release v0.10.2):

| Local file | Upstream path | Git blob | Purpose |
|---|---|---|---|
| `passed_meta_format_P_1.vcf` | `test/input_files/v4.4/passed/passed_meta_format_P_1.vcf` | `efd0044` | Valid VCF 4.4 `FORMAT Number=P` and symbolic CNV record. |
| `failed_body_samples_ploidy_000.vcf` | `test/input_files/v4.1/failed/failed_body_samples_ploidy_000.vcf` | `75ba520` | Invalid PL cardinality for a diploid biallelic genotype. |

Source: <https://github.com/EBIvariation/vcf-validator/tree/0aacc4a44430ab9cee87d9925aacb28a2fb0a9fb/test/input_files>

License: Apache License 2.0, as recorded at the same immutable commit:
<https://github.com/EBIvariation/vcf-validator/blob/0aacc4a44430ab9cee87d9925aacb28a2fb0a9fb/LICENSE>.
Redistribution and modification are permitted under that license. These files
are unmodified; the repository's existing license notices remain separate.
The complete license text is distributed at
`THIRD_PARTY_LICENSES/ebi-vcf-validator-apache-2.0.txt`.

Full source revisions and SHA-256 checksums are recorded in
`THIRD_PARTY_LICENSES/manifest.json` and checked by
`scripts/check-notices.py`.
