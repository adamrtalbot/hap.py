# External adversarial evidence

The audit searched official repositories and release-pinned primary sources.
Only compact files with clear redistribution terms were eligible for import.
Immutable Git object IDs identify upstream bytes; SHA-256 identifies the local
copy used by the matrix.

## Imported fixtures

| Fixture | Immutable source | License | Git blob | SHA-256 | Matrix purpose |
|---|---|---|---|---|---|
| `passed_meta_format_P_1.vcf` | EBI VCF Validator `0aacc4a44430ab9cee87d9925aacb28a2fb0a9fb` | Apache-2.0 | `efd0044c1d958d6583e3b6c9f4aa624e984452e0` | `7a94657671a660a4d4525f6d01234d9d17ae0414cc367546617863180fcb7aba` | Valid VCF 4.4, `Number=P`, symbolic CNV. |
| `failed_body_samples_ploidy_000.vcf` | EBI VCF Validator `0aacc4a44430ab9cee87d9925aacb28a2fb0a9fb` | Apache-2.0 | `75ba520a92741c6687912b39cf1fbc685fe4c5be` | `92271c32027f458fb3a89fe13b11210d3a219beca50329327f4d243776e2c209` | Invalid genotype-likelihood cardinality. |

Upstream tree: <https://github.com/EBIvariation/vcf-validator/tree/0aacc4a44430ab9cee87d9925aacb28a2fb0a9fb/test/input_files>

License at the same commit: <https://github.com/EBIvariation/vcf-validator/blob/0aacc4a44430ab9cee87d9925aacb28a2fb0a9fb/LICENSE>

The redistributed license text is
`THIRD_PARTY_LICENSES/ebi-vcf-validator-apache-2.0.txt`.

The imported files are unmodified. Their local README records their purpose
and provenance beside the bytes.

## Reviewed but not imported

| Source | Immutable revision and license | Candidate behavior | Decision |
|---|---|---|---|
| RTG Tools | `f72c7991210776631b2ee36b8038a64b45deb6da`, BSD-2-Clause | vcfeval phase-obedience, tetraploid paths, BNDs, indexed VCF/TBI | Retain as residual oracle work; templates/triads need a larger fixture and RTG-specific execution lane. |
| GATK | `76edc75c26504da94bbaee66584e107e76ee15de`, Apache-2.0 | Spanning `*` plus `<NON_REF>` gVCFs; realistic Mutect2/Strelka headers and polyploid calls | Do not import until the matching reference interval and caller sample convention can be minimized without losing semantics. |
| GA4GH benchmarking-tools | `e07e6ca9dc372d7b61148ecd0677f7da7cce53a4`, Apache-2.0 | Dense chr21 superloci and relative stratification TSV | Existing matrix already uses pinned chr21 data; a new slice would duplicate large evidence without a newly isolated trigger. |
| Pisces 5.3.0.0 | `becd35f6c3ebeb852bd76034d0ccf494bdac3e1c`, GPL-3.0 | Real multiallelic/depth feature records | Prefer independently authored synthetic records to avoid expanding distribution obligations. Current `master` is PolyForm Strict and was explicitly rejected. |
| Strelka 2.9.10 test data | `1c8f1be330dc6f58568c82468a9da1c5f32fb54a`, GPL-3.0 | Swapped/stale indexed BED pairs | Prefer synthetic index mutations; current `master` licensing is more restrictive and was explicitly rejected. |

The hts-specs malformed corpus, current Strelka/Pisces branches, VarScan
repository examples, and precisionFDA submissions were not imported because
their redistribution terms were absent, ambiguous, or unsuitable. Public
availability alone was not treated as a license.
