# Verification

Read the [website guide](https://adamrtalbot.github.io/hap.py/project/verification/)
for an overview. Use this file for work in `verification/`.

This directory contains the parity gate for the Rust `hap` migration. For every
samplesheet row, Nextflow runs the named legacy hap.py tool in the pinned
reference container and the matching `hap` subcommand from `PATH`, then
compares the emitted `result*` files in the pinned comparator container. nf-test fails when a comparison contains a
difference outside the metadata exclusions below.

## Legacy reference environment

One reference environment governs every lane, and legacy runs there unmodified:
no interpreter shim, no source patch, no environment override. An observation
taken through a patch is not evidence of legacy behavior, so a legacy tool that
will not run means the environment is wrong. The tool version is not the
reference; `hap.py 0.3.15` names more than one behavior. See
`docs/adr/0002-drop-in-claim-against-one-unpatched-legacy-image.md`.

The reference is
`community.wave.seqera.io/library/happy-0.3.15:2c2b5746d6b0da37`, the frozen Wave
build of `containers/happy-0.3.15.yml`. It carries all six legacy tools plus
`rtg`, `bcftools`, `samtools`, and `java`. Because the digest is opaque,
`containers/happy-0.3.15.conda-lock.txt` records the 108 packages it installs,
including the `libstdcxx` and `openjdk` builds the recipe never names.
Regenerate it with `scripts/dump-conda-lock.sh <image>`. The recipe lists what
was requested; the lock lists what was installed, and the lock is what explains a
number.

It pins pandas 0.20.3, which is the version window where one image runs all six
tools unpatched. Measured on the image: index-level `groupby` works,
`display.height` survives as a deprecation so som.py needs no shim, and `to_csv`
renders every float as Python 2 `str()` does, identically to the 0.19.2 build it
replaces. numpy 1.12.1, scipy 1.2.1, and libstdcxx 16.1.0 did not move. The
rebuild did move 17 peripheral packages the earlier solve had resolved
differently; the lock is the record of which.

The comparator is pinned the same way, because it decides pass and fail.
`DIFF_OUTPUTS` runs
`community.wave.seqera.io/library/hap-comparator:3fd84dea2dc8fda1`, the frozen
Wave build of `containers/hap-comparator.yml`, with its contents in
`containers/hap-comparator.conda-lock.txt`. It carries python 3.12.3 and the same
`bcftools` 1.17 build the reference image carries, so a BCF decode difference
cannot be an htslib-version artifact. It previously inherited `container = null`
and ran host `python3` and host `bcftools`, so a host change could flip a verdict
without either implementation moving.

Moving a pin re-baselines: legacy is re-run on the new image and its fresh output
is the expectation. Nothing is preserved, because nothing is stored: `results/`
is gitignored, there are no nf-test snapshots, and the harness runs legacy and
`hap` together in one execution and diffs them live.

The gate that remains is the ordinary one. `hap` is frozen while legacy moves, so
a legacy change surfaces immediately as a failed comparison. See
`docs/adr/0004-pin-the-legacy-baseline-to-one-container-identity.md`.

The four `--bam` rows are out of the truth set: `ftx_bam_depth`, `ftx_multi_bam`,
`somatic_bam_depth`, `matrix_multi_bam`. The original never pinned pandas, so it
permits installations where its own `--bam` paths cannot run, and any observation
of them reflects the pandas version this project picks rather than legacy
behavior. `--bam` is classed **no legacy reference** in the compatibility policy,
which makes hap-rs behavior there normative. Their fixtures stay under
`assets/fixtures/ftx-bam/` for the `normative_` expectations, which are separate
work. The matrix is 155 six-lane comparisons. Coverage below tables every row.

The line that rule draws is the layer a setting reaches. Configure the
interpreter or VM legacy runs on, or the order the harness schedules it in, and
that is a harness constraint. Change the behaviour of a library the tool calls
and that is an override. `maxForks = 1` on `HAPPY_LEGACY` constrains scheduling;
`RTG_JAVA_OPTS=-Xint` on the legacy vcfeval path in `modules/happy.nf` puts the
JVM in interpreted mode. The `sitecustomize.py` shim changed how pandas answered
som.py, so it fails the same test. Allowed settings are named in
`docs/adr/0004-pin-the-legacy-baseline-to-one-container-identity.md`; a setting
not named there is not allowed.

Parity is observed on linux/amd64. The conda lock is entirely `linux-64`, so the
image has no other build, and Docker on an arm64 host emulates it. CI runs
`ubuntu-24.04` natively, so authoritative runs already satisfy this. A local run
on arm64 is a development diagnostic, not parity evidence.

## Lanes

| Lane | Legacy command | Rust command |
|---|---|---|
| HAPPY | `hap.py` | `hap germline` |
| SOMPY | `som.py` | `hap somatic` |
| PREPY | `pre.py` | `hap pre` |
| FTXPY | `ftx.py` | `hap ftx` |
| QFY | `qfy.py` | `hap quantify` |
| VCFCHECK | `vcfcheck` | `hap validate` |

The six `assets/samplesheet.*.csv` files list the test matrix. Each row supplies
inputs and arguments for one case. Add a row to reproduce a
discrepancy. Keep existing rows, comparison rules, and expected outputs
unchanged.

## Run

Build the Rust binary, make it available on `PATH`, then run the gate:

```bash
cargo build --release
export PATH="$(pwd)/target/release:$PATH"
cd verification
nf-test test --ci tests/main.nf.test
```

Use `HAP_TEST_CASES` to diagnose one lane:

```bash
export PATH="$(pwd)/target/release:$PATH"
cd verification
HAP_TEST_CASES=sompy nf-test test --ci tests/main.nf.test
```

CI runs all six lanes. Use a narrowed run to inspect one command.

## Run the public germline discovery workflow

Public-data cases stay outside nf-test. Build `hap`, then run the intact public
germline samplesheet directly with Nextflow. The release-binary directory on
`PATH` must be absolute because Nextflow tasks run from their own work
directories:

```bash
cargo build --release --bin hap
export PATH="$(pwd)/target/release:$PATH"
cd verification
nextflow run main.nf \
  --cases happy \
  --happy_samplesheet "$PWD/assets/samplesheet.happy.public.csv" \
  --outdir "$PWD/results/public-germline"
```

Accept the run only when `results/public-germline/verification.json` reports
`comparison_count: 1` and that comparison has `ok: true` with an empty
`differences` list. The one comparison is the intact `germline_hg001_platinum`
public-data case; the two HAPPY contract comparisons this paragraph used to
count left the pipeline in `1a1bbca`. Also inspect
`results/public-germline/happy/germline_hg001_platinum/comparison.json`; it must
independently report `ok: true`, an empty `differences` list, and identical
legacy and hap-rs artifact inventories.

`--engine vcfeval` uses separate reference formats for the two implementations.
Nextflow gives the pinned legacy image an RTG Tools 3.12.1 SDF `.tar.gz` bundle
from the samplesheet's `reference_sdf` column. It gives `hap` the corresponding
FASTA. The product runs without RTG or Java. A README beside each fixture
records the bundle checksum and provenance.

## Coverage

Measured on 2026-08-18 at `18548d1` by running `nextflow run main.nf` with the
default `--cases` list, twice. The measurement drives the pipeline directly
rather than through nf-test, because it needs `-with-trace` and two separate
output directories. `nf-test` adds the assertions: its per-lane expected counts
now sum to 155 and match what the pipeline emitted in both runs, while its
`every { it.ok }` assertion still fails on the two cases below.

The six samplesheets hold 155 rows and the gate emits 155 comparisons, of which
153 have an empty difference list in each run.
The pair moved between the runs: run 1 reported HAPPY `chr21` and
`chr21_region`, run 2 reported `chr21_passonly` and `chr21_xcmp_controls`, and
all four differences sit in ROC artifacts. The replay section below records
hap-rs emitting identical bytes for all four in both runs while the legacy ROC
tables moved. Tracked as
[#46](https://github.com/adamrtalbot/hap.py/issues/46).

| Lane | Command | Rows |
|---|---|---:|
| HAPPY | `hap germline` | 34 |
| SOMPY | `hap somatic` | 28 |
| PREPY | `hap pre` | 47 |
| FTXPY | `hap ftx` | 25 |
| QFY | `hap quantify` | 9 |
| VCFCHECK | `hap validate` | 12 |
| | | **155** |

The two hap-rs vcfeval contract cases are gone. `1a1bbca` deleted the
`happy_contracts` channel, and the leftovers it named have since been removed:
`HAPPY_RUST_CONTRACT`, its import in `main.nf`, and
`scripts/validate_vcfeval_contract.py`
([#47](https://github.com/adamrtalbot/hap.py/issues/47)). The gate is 155
comparisons with no contract case.

One row per committed case follows. The option-bindings column lists what the
row's `args` cell binds. Every row also carries the bindings the harness adds in
`modules/`: `--reference` and the output prefix everywhere, `--threads` on
HAPPY, PREPY and QFY, `--false-positives` on HAPPY, SOMPY and QFY,
`--feature-table` on SOMPY and FTXPY, `-R` on every PREPY row, and
`--engine-vcfeval-template` on the legacy side of the five HAPPY rows that carry
an SDF bundle.

#### HAPPY: `hap germline`, 34 rows

| Case | Input format and index | Engine or table | Region or stratification | Filtering | Other option bindings |
|---|---|---|---|---|---|
| `chr21` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | — | — | — |
| `chr21_passonly` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | — | `--pass-only` | — |
| `chr21_region` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--location chr21:15000000-20000000` | — | — |
| `chr21_passonly_region` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--location chr21:15000000-20000000` | `--pass-only` | — |
| `chr21_artifact_controls` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--location chr21:15000000-20000000` | — | `--no-roc`, `--no-write-counts`, `--no-json` |
| `chr21_unhappy` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--location chr21:15000000-20000000` | — | `--unhappy` |
| `chr21_scmp_somatic` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `scmp-somatic` | `--location chr21:15000000-20000000` | — | — |
| `synthetic_scmp_distance` | truth vcf; query vcf; reference fasta + .fai; fp bed bed.gz + .tbi | `scmp-distance` | — | — | `--scmp-distance 30`, `--set-gt half`, `--no-leftshift`, `--no-decompose` |
| `chr21_xcmp_controls` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--location chr21:15000000-20000000` | — | `--window-size 75`, `--xcmp-enumeration-threshold 8192`, `--xcmp-expand-hapblocks 45` |
| `chr21_preprocess_controls` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--no-adjust-conf-regions`, `--location chr21:15000000-20000000` | — | `--preprocess-truth`, `--preprocessing-window-size 4096` |
| `chr21_filtered_truth` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--location chr21:15000000-20000000` | `--usefiltered-truth` | — |
| `synthetic_filtered_truth_roc` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | `--usefiltered-truth` | — |
| `chr21_output_metadata` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | xcmp (default) | `--location chr21:15000000-20000000` | — | `--output-vtc`, `--preserve-info` |
| `matrix_bcf` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--bcf` |
| `matrix_implicit_bcf` | truth bcf + .csi; query bcf + .csi; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | — |
| `matrix_gvcf_conversion` | truth gvcf; query gvcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | `--filter-nonref` | `--convert-gvcf-truth`, `--convert-gvcf-query`, `--preprocess-truth` |
| `matrix_restrict_filters` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | `--restrict-regions confident.bed` | `--filters-only LowQual` | — |
| `matrix_targets` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | `--target-regions confident.bed` | — | — |
| `matrix_roc_stratification` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | `--roc-regions EXTRA`, `--stratification-region EXTRA:confident.bed` | `--roc-filter LowQual` | `--roc INFO.SCORE`, `--roc-delta 5`, `--ci-alpha 0.05`, `--preserve-info` |
| `matrix_gt_fix_normalize` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--fixchr`, `--set-gt hom`, `--preprocess-truth`, `--decompose`, `--bcftools-norm` |
| `matrix_scmp_default_first` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `scmp-distance` | — | — | `--lose-match-distance 7`, `--fixchr`, `--no-leftshift`, `--no-decompose` |
| `matrix_vcfeval_real` | truth vcf; query vcf; reference fasta + .fai; fp bed bed; vcfeval sdf sdf.tar.gz | `vcfeval` | — | — | `--scmp-distance 7`, `--no-leftshift`, `--no-decompose` |
| `matrix_vcfeval_paths` | truth vcf; query vcf; reference fasta + .fai; fp bed bed; vcfeval sdf sdf.tar.gz | `vcfeval` | `--no-adjust-conf-regions` | — | `--scmp-distance 4`, `--no-leftshift`, `--no-decompose` |
| `matrix_vcfeval_pass_only` | truth vcf; query vcf; reference fasta + .fai; fp bed bed; vcfeval sdf sdf.tar.gz | `vcfeval` | `--no-adjust-conf-regions` | `--pass-only` | `--scmp-distance 4`, `--no-leftshift`, `--no-decompose` |
| `matrix_vcfeval_custom_roc` | truth vcf; query vcf; reference fasta + .fai; fp bed bed; vcfeval sdf sdf.tar.gz | `vcfeval` | `--no-adjust-conf-regions` | — | `--scmp-distance 4`, `--no-leftshift`, `--no-decompose`, `--roc INFO.SCORE` |
| `matrix_vcfeval_bcf` | truth vcf; query vcf; reference fasta + .fai; fp bed bed; vcfeval sdf sdf.tar.gz | `vcfeval` | `--no-adjust-conf-regions` | — | `--scmp-distance 4`, `--no-leftshift`, `--no-decompose`, `--bcf` |
| `matrix_stratification_tsv_fixchr` | truth vcf; query vcf; reference fasta + .fai; fp bed bed; stratification bed/tsv | xcmp (default) | `--stratification stratification/regions.tsv`, `--stratification-fixchr` | — | `--fixchr`, `--write-vcf` |
| `matrix_stratification_tsv_multiple_beds` | truth vcf; query vcf; reference fasta + .fai; fp bed bed; stratification bed/tsv | xcmp (default) | `--stratification stratification/regions.tsv`, `--stratification-fixchr` | — | `--fixchr`, `--write-vcf` |
| `matrix_write_vcf_vtc_counts` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--write-vcf`, `--output-vtc`, `--preserve-info`, `--no-write-counts`, `--write-counts` |
| `matrix_gvcf_conversion_both` | truth gvcf; query gvcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--convert-gvcf-to-vcf`, `--preprocess-truth` |
| `matrix_gender_male` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--gender male`, `--preprocess-truth`, `--write-vcf` |
| `matrix_somatic_mode` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--somatic`, `--preprocess-truth`, `--gender none`, `--write-vcf` |
| `hg001_graph_record_order` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--no-json`, `--no-roc`, `--no-write-counts`, `--write-vcf`, `--gender none` |
| `hg001_float_rendering` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | xcmp (default) | — | — | `--no-roc`, `--gender none` |

#### SOMPY: `hap somatic`, 28 rows

| Case | Input format and index | Engine or table | Region or stratification | Filtering | Other option bindings |
|---|---|---|---|---|---|
| `chr21_includenonpass_countunk` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `hcc.strelka.indel` | — | `--include-nonpass` | `--count-unk` |
| `chr21_includenonpass` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `hcc.strelka.indel` | — | `--include-nonpass` | — |
| `chr21_passonly_countunk` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `hcc.strelka.indel` | — | — | `--count-unk` |
| `chr21_minimal` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `generic` | — | — | — |
| `grch38_indels` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `hcc.strelka.indel` | — | — | — |
| `grch38_indels_happy` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `hcc.strelka.indel` | — | `--include-nonpass` | — |
| `grch38_indels_af` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `hcc.strelka.indel` | — | — | `--bin-afs` |
| `grch38_snvs` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `admix.strelka.snv` | — | — | — |
| `grch38_snvs_roc` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `admix.strelka.snv` + `--roc strelka.snv` | — | — | — |
| `normalize_query` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--normalize-query` |
| `normalize_all_filtered_fn` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | `--include-nonpass` | `--normalize-all`, `--count-filtered-fn` |
| `matrix_auto_query_contig` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | — |
| `matrix_count_toggle_on` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--no-count-unk`, `--count-unk`, `--fp-region-size 30` |
| `matrix_count_toggle_off` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--count-unk`, `--no-count-unk`, `--fp-region-size 30` |
| `matrix_ambiguity_last_off` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--ambiguous classification.bed`, `--ambiguous classification.bed`, `--explain_ambiguous`, `--ambi-fp`, `--no-ambi-fp` |
| `matrix_ambiguity_last_on` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--ambiguous classification.bed`, `--explain_ambiguous`, `--no-ambi-fp`, `--ambi-fp` |
| `matrix_location_fp_denominator` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | `--location chr1` | — | `--count-unk` |
| `matrix_restrict_regions` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | `--restrict-regions selector.bed` | — | — |
| `matrix_target_regions` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | `--target-regions selector.bed` | — | — |
| `matrix_normalize_truth_fixchr_ci` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--normalize-truth`, `--no-fixchr-truth`, `--fix-chr-truth`, `--no-fixchr-query`, `--ci-level 0.99`, `--fp-region-size 30` |
| `matrix_snv_roc_af` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` + `--roc strelka.snv.qss` | — | `--include-nonpass` | `--bin-afs`, `--af-binsize 0.5`, `--af-truth T_AF`, `--af-query T_AF` |
| `matrix_fp_region_auto` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--count-unk`, `--fp-region-size auto` |
| `matrix_fixchr_query_alias` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` | — | — | `--no-fixchr-query`, `--fix-chr-query`, `--no-fixchr-truth`, `--fp-region-size 30` |
| `matrix_roc_strelka_indel_evs` | truth vcf.gz; query vcf.gz; reference fasta + .fai; fp bed bed.gz | `hcc.strelka.indel` + `--roc strelka.indel.evs` | — | `--include-nonpass` | — |
| `matrix_roc_mutect_snv` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` + `--roc mutect.snv` | — | `--include-nonpass` | — |
| `matrix_roc_varscan2_snv` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `generic` + `--roc varscan2.snv` | — | `--include-nonpass` | — |
| `matrix_strelka_indel_af_type_compat` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `hcc.strelka.indel` | — | `--include-nonpass` | `--bin-afs` |
| `somatic_output_parity_bin_afs` | truth vcf; query vcf; reference fasta + .fai; fp bed bed | `hcc.strelka.indel` | — | `--include-nonpass` | `--bin-afs` |

#### PREPY: `hap pre`, 47 rows

| Case | Input format and index | Engine or table | Region or stratification | Filtering | Other option bindings |
|---|---|---|---|---|---|
| `chr21` | input vcf.gz + .tbi; reference fasta + .fai; regions bed bed.gz + .tbi | — | — | — | — |
| `chr21_passonly` | input vcf.gz + .tbi; reference fasta + .fai; regions bed bed.gz + .tbi | — | — | `--pass-only` | — |
| `chr21_no_fixchr` | input vcf.gz + .tbi; reference fasta + .fai; regions bed bed.gz + .tbi | — | — | — | `--no-fixchr` |
| `tiny_default` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_somatic` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--somatic` |
| `tiny_set_gt_half` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--somatic`, `--set-gt half` |
| `tiny_set_gt_hemi` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--somatic`, `--set-gt hemi` |
| `tiny_set_gt_het` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--somatic`, `--set-gt het` |
| `tiny_set_gt_hom` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--somatic`, `--set-gt hom` |
| `tiny_set_gt_first` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--somatic`, `--set-gt first` |
| `tiny_set_gt_then_somatic` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--set-gt first`, `--somatic` |
| `tiny_convert_gvcf` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--convert-gvcf-to-vcf` |
| `tiny_no_normalize` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `tiny_filters_only` | input vcf; reference fasta + .fai; regions bed bed | — | — | `--filters-only LowQ` | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `tiny_window_size` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--window-size 4096` |
| `tiny_switch_precedence` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--no-leftshift`, `--leftshift`, `--no-decompose`, `--decompose` |
| `tiny_gender_auto` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--gender auto` |
| `tiny_gender_half_call` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--gender auto` |
| `tiny_gender_male` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--gender male` |
| `tiny_gender_female` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--gender female` |
| `tiny_fixchr_precedence` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--no-fixchr`, `--fixchr`, `--gender female` |
| `tiny_bcftools_norm` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--bcftools-norm`, `--no-leftshift`, `--no-decompose`, `--gender none` |
| `tiny_bcftools_norm_default` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--bcftools-norm`, `--gender none` |
| `tiny_bcf_output` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--bcf` |
| `tiny_verbose_log` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--verbose`, `--logfile pre.log`, `--force-interactive` |
| `tiny_symbolic_del` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_edge_alleles` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_multisample_star_phase` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_symbolic_types` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_breakend_forms` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_gvcf_edge` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--convert-gvcf-to-vcf` |
| `tiny_ploidy_format` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_boundaries` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `tiny_boundary_symbolic` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | — |
| `matrix_location_postfix` | input vcf; reference fasta + .fai; regions bed bed | — | `--location chr1:2-2` | — | `--fixchr`, `--no-leftshift`, `--no-decompose`, `--gender none` |
| `matrix_fixchr_auto` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `matrix_fixchr_disabled` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--no-fixchr`, `--no-leftshift`, `--no-decompose`, `--gender none` |
| `matrix_region_end_overlap` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `matrix_target_start_only` | input vcf; reference fasta + .fai; regions bed bed | — | `--target-regions end-boundary.bed` | — | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `matrix_filters_any_unexcluded` | input vcf; reference fasta + .fai; regions bed bed | — | — | `--filters-only LowQ` | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `matrix_pass_wins_filters` | input vcf; reference fasta + .fai; regions bed bed | — | — | `--filters-only LowQ`, `--pass-only` | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `matrix_window_isolated` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--window-size 1`, `--gender none` |
| `matrix_window_neighbours` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--window-size 10000`, `--gender none` |
| `matrix_overlapping_locations` | input vcf; reference fasta + .fai; regions bed bed (`--threads 2`) | — | `--location chr1:1-10,chr1:3-20` | — | `--window-size 1`, `--gender none` |
| `matrix_bcf_input` | input bcf + .csi; reference fasta + .fai; regions bed bed | — | — | — | `--no-leftshift`, `--no-decompose`, `--gender none` |
| `chr21_parallel_blocksplit` | input vcf.gz + .tbi; reference fasta + .fai; regions bed bed.gz + .tbi (`--threads 2`) | — | — | — | `--window-size 1` |
| `hg001_primitive_identity` | input vcf; reference fasta + .fai; regions bed bed | — | — | — | `--decompose`, `--leftshift`, `--gender none` |

#### FTXPY: `hap ftx`, 25 rows

| Case | Input format and index | Engine or table | Region or stratification | Filtering | Other option bindings |
|---|---|---|---|---|---|
| `chr21` | input vcf.gz; reference fasta + .fai | `generic` | — | — | — |
| `chr21_includenonpass` | input vcf.gz; reference fasta + .fai | `generic` | — | `--include-nonpass` | — |
| `chr21_hcc_strelka_snv` | input vcf.gz; reference fasta + .fai | `hcc.strelka.snv` | — | — | — |
| `chr21_admix_strelka_snv` | input vcf.gz; reference fasta + .fai | `admix.strelka.snv` | — | — | — |
| `grch38_hcc_strelka_snv` | input vcf.gz; reference fasta + .fai | `hcc.strelka.snv` | — | — | — |
| `grch38_hcc_strelka_indel` | input vcf.gz; reference fasta + .fai | `hcc.strelka.indel` | — | — | — |
| `grch38_admix_strelka_indel` | input vcf.gz; reference fasta + .fai | `admix.strelka.indel` | — | — | — |
| `chr21_hcc_strelka_snv_tp` | input vcf.gz; reference fasta + .fai | `hcc.strelka.snv` | — | — | `--feature-label TP` |
| `chr21_hcc_strelka_snv_fn` | input vcf.gz; reference fasta + .fai | `hcc.strelka.snv` | — | — | `--feature-label FN` |
| `ftx_called_nonref` | input vcf; reference fasta + .fai | `generic` | — | — | — |
| `ftx_location` | input vcf; reference fasta + .fai | `generic` | `--location 1:7` | — | `--feature-label reference` |
| `ftx_fix_chr` | input vcf; reference fasta + .fai | `generic` | — | — | `--fix-chr`, `--feature-label reference` |
| `ftx_normalize` | input vcf; reference fasta + .fai | `generic` | — | — | `--normalize`, `--feature-label reference` |
| `ftx_mutect_snv` | input vcf; reference fasta + .fai | `hcc.mutect.snv` | — | — | `--feature-label reference` |
| `ftx_varscan2_snv` | input vcf; reference fasta + .fai | `hcc.varscan2.snv` | — | — | `--feature-label reference` |
| `ftx_pisces_snv` | input vcf; reference fasta + .fai | `hcc.pisces.snv` | — | — | `--feature-label reference` |
| `ftx_empty_label` | input vcf; reference fasta + .fai | `generic` | — | — | `--feature-label ` |
| `ftx_mutect_indel` | input vcf; reference fasta + .fai | `hcc.mutect.indel` | — | — | `--feature-label reference` |
| `ftx_mutect_reversed_samples` | input vcf; reference fasta + .fai | `hcc.mutect.snv` | — | — | `--feature-label reference` |
| `ftx_varscan2_indel` | input vcf; reference fasta + .fai | `hcc.varscan2.indel` | — | — | `--feature-label reference` |
| `ftx_pisces_indel` | input vcf; reference fasta + .fai | `hcc.pisces.indel` | — | — | `--feature-label reference` |
| `ftx_selector_normalize` | input vcf; reference fasta + .fai | `generic` | `--location chr1:5-15` | `--include-nonpass` | `--fix-chr`, `--normalize`, `--feature-label reference` |
| `ftx_empty_input` | input vcf; reference fasta + .fai | `generic` | — | — | `--feature-label empty-reference` |
| `ftx_restrict_regions` | input vcf; reference fasta + .fai; regions bed bed | `generic` | `--restrict-regions end-boundary.bed` | — | `--feature-label reference` |
| `ftx_target_regions` | input vcf; reference fasta + .fai; regions bed bed | `generic` | `--target-regions end-boundary.bed` | — | `--feature-label reference` |

#### QFY: `hap quantify`, 9 rows

| Case | Input format and index | Engine or table | Region or stratification | Filtering | Other option bindings |
|---|---|---|---|---|---|
| `chr21_region` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `--type xcmp` | — | — | — |
| `chr21_ga4gh` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `--type ga4gh` | — | — | — |
| `chr21_ga4gh_roc_controls` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `--type ga4gh` | `--roc-regions CONF` | `--roc-filter LowGQX` | `--roc QQ`, `--roc-delta 10`, `--ci-alpha 0.05` |
| `chr21_ga4gh_artifacts` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `--type ga4gh` | — | — | `--write-vcf`, `--output-vtc`, `--bcf`, `--threads 1`, `--force-interactive`, `--logfile qfy.log`, `--verbose`, `--no-write-counts`, `--no-roc`, `--no-json` |
| `chr21_counts_last_off` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `--type ga4gh` | — | — | `--write-counts`, `--no-write-counts`, `--no-roc`, `--no-json` |
| `chr21_counts_last_on` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `--type ga4gh` | — | — | `--report-prefix ignored`, `--roc QUAL`, `--roc INFO.QQ`, `--no-write-counts`, `--write-counts`, `--no-roc`, `--no-json` |
| `chr21_four_column_strat` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed | `--type ga4gh` | `--stratification-region EXTRA:stratification-four-column.bed` | — | `--no-roc` |
| `chr21_tsv_fixchr` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi; stratification bed/tsv | `--type ga4gh` | `--stratification stratification/regions.tsv`, `--stratification-fixchr` | — | `--no-roc` |
| `chr21_adjust_conf_regions` | truth vcf.gz + .tbi; query vcf.gz + .tbi; reference fasta + .fai; fp bed bed.gz + .tbi | `--type ga4gh` | `--adjust-conf-regions PG_NA12878_chr21.vcf.gz` | — | `--write-vcf` |

Every QFY row quantifies the same artifact: `QFY_ANNOTATE` runs the pinned
`hap.py` once per row over the chr21 fixtures with a fixed
`-l chr21:15000000-20000000 -V --no-json --no-roc`, and both implementations
then read that `seed.vcf.gz`. The input-format column above records what the
seeding step consumes, not what `qfy.py` reads.

#### VCFCHECK: `hap validate`, 12 rows

| Case | Input format and index | Engine or table | Region or stratification | Filtering | Other option bindings |
|---|---|---|---|---|---|
| `formats` | input vcf | — | — | — | — |
| `warning_controls` | input vcf | — | — | `--apply-filters true` | `--limit-records 4`, `--message-every 1`, `--strict-homref true`, `--all-warnings true`, `--check-bcf-errors false` |
| `ploidy` | input vcf | — | — | — | `--limit-records -1` |
| `negative_limit` | input vcf | — | — | — | `--limit-records -2` |
| `missing_contig` | input vcf | — | — | — | `--check-bcf-errors false` |
| `sample_count_mismatch` | input vcf | — | — | — | `--check-bcf-errors false` |
| `lowercase_x` | input vcf | — | — | — | `--check-bcf-errors false` |
| `vcf44_number_p` | input vcf | — | — | — | `--check-bcf-errors false` |
| `invalid_pl_cardinality` | input vcf | — | — | — | `--check-bcf-errors false` |
| `edge_alleles` | input vcf | — | — | — | `--check-bcf-errors false` |
| `location_input_option` | input vcf.gz + .tbi (`--input-file`) | — | `--location chr21:15000000-16000000` | — | — |
| `explicit_boolean_values` | input vcf | — | — | `--apply-filters false` | `--strict-homref false`, `--all-warnings false`, `--check-bcf-errors true` |

#### Option interactions that exist

There is no pairwise scheme. These are the factor combinations the committed
rows bind, counted.

| Lane | Engine or table | Input | Filtering | Region | Stratification | Rows |
|---|---|---|---|---|---|---:|
| HAPPY | `scmp-distance` | vcf | — | — | — | 2 |
| HAPPY | `scmp-somatic` | vcf.gz | — | `--location` | — | 1 |
| HAPPY | `vcfeval` | vcf | `--pass-only` | — | — | 1 |
| HAPPY | `vcfeval` | vcf | — | — | — | 4 |
| HAPPY | xcmp (default) | bcf | — | — | — | 1 |
| HAPPY | xcmp (default) | gvcf | `--filter-nonref` | — | — | 1 |
| HAPPY | xcmp (default) | gvcf | — | — | — | 1 |
| HAPPY | xcmp (default) | vcf | `--filters-only` | `--restrict-regions` | — | 1 |
| HAPPY | xcmp (default) | vcf | `--usefiltered-truth` | — | — | 1 |
| HAPPY | xcmp (default) | vcf | — | `--target-regions` | — | 1 |
| HAPPY | xcmp (default) | vcf | — | — | yes | 3 |
| HAPPY | xcmp (default) | vcf | — | — | — | 7 |
| HAPPY | xcmp (default) | vcf.gz | `--pass-only` | `--location` | — | 1 |
| HAPPY | xcmp (default) | vcf.gz | `--pass-only` | — | — | 1 |
| HAPPY | xcmp (default) | vcf.gz | `--usefiltered-truth` | `--location` | — | 1 |
| HAPPY | xcmp (default) | vcf.gz | — | `--location` | — | 6 |
| HAPPY | xcmp (default) | vcf.gz | — | — | — | 1 |
| SOMPY | `admix.strelka.snv` | vcf.gz | — | — | — | 2 |
| SOMPY | `generic` | vcf | `--include-nonpass` | — | — | 4 |
| SOMPY | `generic` | vcf | — | `--location` | — | 1 |
| SOMPY | `generic` | vcf | — | `--restrict-regions` | — | 1 |
| SOMPY | `generic` | vcf | — | `--target-regions` | — | 1 |
| SOMPY | `generic` | vcf | — | — | — | 9 |
| SOMPY | `generic` | vcf.gz | — | — | — | 1 |
| SOMPY | `hcc.strelka.indel` | vcf | `--include-nonpass` | — | — | 1 |
| SOMPY | `hcc.strelka.indel` | vcf.gz | `--include-nonpass` | — | — | 4 |
| SOMPY | `hcc.strelka.indel` | vcf.gz | — | — | — | 3 |
| PREPY | — | bcf | — | — | — | 1 |
| PREPY | — | vcf | `--filters-only` | — | — | 2 |
| PREPY | — | vcf | `--filters-only`, `--pass-only` | — | — | 1 |
| PREPY | — | vcf | — | `--location` | — | 2 |
| PREPY | — | vcf | — | `--target-regions` | — | 1 |
| PREPY | — | vcf | — | — | — | 36 |
| PREPY | — | vcf.gz | `--pass-only` | — | — | 1 |
| PREPY | — | vcf.gz | — | — | — | 3 |
| FTXPY | `admix.strelka.indel` | vcf.gz | — | — | — | 1 |
| FTXPY | `admix.strelka.snv` | vcf.gz | — | — | — | 1 |
| FTXPY | `generic` | vcf | `--include-nonpass` | `--location` | — | 1 |
| FTXPY | `generic` | vcf | — | `--location` | — | 1 |
| FTXPY | `generic` | vcf | — | `--restrict-regions` | — | 1 |
| FTXPY | `generic` | vcf | — | `--target-regions` | — | 1 |
| FTXPY | `generic` | vcf | — | — | — | 5 |
| FTXPY | `generic` | vcf.gz | `--include-nonpass` | — | — | 1 |
| FTXPY | `generic` | vcf.gz | — | — | — | 1 |
| FTXPY | `hcc.mutect.indel` | vcf | — | — | — | 1 |
| FTXPY | `hcc.mutect.snv` | vcf | — | — | — | 2 |
| FTXPY | `hcc.pisces.indel` | vcf | — | — | — | 1 |
| FTXPY | `hcc.pisces.snv` | vcf | — | — | — | 1 |
| FTXPY | `hcc.strelka.indel` | vcf.gz | — | — | — | 1 |
| FTXPY | `hcc.strelka.snv` | vcf.gz | — | — | — | 4 |
| FTXPY | `hcc.varscan2.indel` | vcf | — | — | — | 1 |
| FTXPY | `hcc.varscan2.snv` | vcf | — | — | — | 1 |
| QFY | `ga4gh` | vcf.gz | — | — | yes | 2 |
| QFY | `ga4gh` | vcf.gz | — | — | — | 6 |
| QFY | `xcmp` | vcf.gz | — | — | — | 1 |
| VCFCHECK | — | vcf | `--apply-filters` | — | — | 2 |
| VCFCHECK | — | vcf | — | — | — | 9 |
| VCFCHECK | — | vcf.gz | — | `--location` | — | 1 |

Eleven rows bind a flag together with its negation, or bind the same option
twice, so the gate observes precedence rather than only the single-binding
case.

| Lane | Case | Bindings in order | Interaction |
|---|---|---|---|
| HAPPY | `matrix_write_vcf_vtc_counts` | `-V --output-vtc --preserve-info --no-write-counts -X` | flag and its negation, last binding wins |
| SOMPY | `matrix_count_toggle_on` | `--no-count-unk --count-unk --fp-region-size 30` | flag and its negation, last binding wins |
| SOMPY | `matrix_count_toggle_off` | `--count-unk --no-count-unk --fp-region-size 30` | flag and its negation, last binding wins |
| SOMPY | `matrix_ambiguity_last_off` | `-a classification.bed -a classification.bed --explain_ambiguous --ambi-fp --no-ambi-fp` | flag and its negation, last binding wins; repeated `--ambiguous` |
| SOMPY | `matrix_ambiguity_last_on` | `-a classification.bed --explain_ambiguous --no-ambi-fp --ambi-fp` | flag and its negation, last binding wins |
| SOMPY | `matrix_normalize_truth_fixchr_ci` | `--normalize-truth --no-fixchr-truth --fix-chr-truth --no-fixchr-query --ci-level 0.99 --fp-region-size 30` | flag and its negation, last binding wins |
| SOMPY | `matrix_fixchr_query_alias` | `--no-fixchr-query --fix-chr-query --no-fixchr-truth --fp-region-size 30` | flag and its negation, last binding wins |
| PREPY | `tiny_switch_precedence` | `--no-leftshift -L --no-decompose --decompose` | flag and its negation, last binding wins |
| PREPY | `tiny_fixchr_precedence` | `--no-fixchr --fixchr --gender female` | flag and its negation, last binding wins |
| QFY | `chr21_counts_last_off` | `--type ga4gh --write-counts --no-write-counts --no-roc --no-json` | flag and its negation, last binding wins |
| QFY | `chr21_counts_last_on` | `--type ga4gh --report-prefix ignored --roc QUAL --roc INFO.QQ --no-write-counts --write-counts --no-roc --no-json` | flag and its negation, last binding wins; repeated `--roc` |

Two PREPY rows bind the same pair of options in opposite order:
`tiny_set_gt_first` binds `--somatic --set-gt first` and
`tiny_set_gt_then_somatic` binds `--set-gt first --somatic`. QFY
`chr21_counts_last_on` also binds `--report-prefix ignored`, which the harness
binds again as `--report-prefix result`.

## Factors and gaps

The denominator is the option surface of the pinned parsers, captured from the
reference image on 2026-08-18 with `--help`: 62 long options for `hap.py`, 41
for `som.py`, 26 for `pre.py`, 12 for `ftx.py`, 27 for `qfy.py`, 11 for
`vcfcheck`, 179 in total. The counts match the ones
`docs/adr/0003-bound-the-invocation-surface-to-the-pinned-parsers.md` records.

An option is **covered** when a committed row or the harness binds it, **absent
by decision** when no row binds it and a recorded decision puts it outside the
claim, and a **gap** otherwise. Where the parser enumerates values, the note
names which values are bound and which are not.

Binding wins over the decision, so an option can be covered on one command and
absent by decision on another: `--logfile` and `--verbose` are bound on `hap pre`
and `hap quantify` and on no other row. The decision still governs what gets
compared. A row binding `--logfile` exercises the option; the log file it writes
is outside the observed `result*` prefix set either way.

| Command | Options | Covered | Absent by decision | Gap |
|---|---:|---:|---:|---:|
| `hap germline` | 62 | 49 | 7 | 6 |
| `hap somatic` | 41 | 29 | 7 | 5 |
| `hap pre` | 26 | 22 | 3 | 1 |
| `hap ftx` | 12 | 10 | 2 | 0 |
| `hap quantify` | 27 | 23 | 3 | 1 |
| `hap validate` | 11 | 9 | 2 | 0 |
| | **179** | **142** | **24** | **13** |

#### `hap germline`, 49 of 62 options bound

| Option | Status | Rows | Note |
|---|---|---:|---|
| `--help` | absent by decision | 0 | standard output content is outside the claim (ADR 0003) |
| `--version` | absent by decision | 0 | prints a version and exits 0 on every subcommand, no comparison possible (ADR 0003 correction) |
| `--reference` | covered | 34 | bound by the harness on every row |
| `--report-prefix` | covered | 34 | bound by the harness on every row |
| `--scratch-prefix` | absent by decision | 0 | scratch placement is a filesystem side effect outside the claim (ADR 0005) |
| `--keep-scratch` | absent by decision | 0 | scratch placement is a filesystem side effect outside the claim (ADR 0005) |
| `--type` | gap | 0 | annotation format; `hap quantify` covers both values |
| `--false-positives` | covered | 34 | bound by the harness on every row |
| `--stratification` | covered | 2 |  |
| `--stratification-region` | covered | 1 |  |
| `--stratification-fixchr` | covered | 2 |  |
| `--write-vcf` | covered | 6 |  |
| `--write-counts` | covered | 1 |  |
| `--no-write-counts` | covered | 3 |  |
| `--output-vtc` | covered | 2 |  |
| `--preserve-info` | covered | 3 |  |
| `--roc` | covered | 2 |  |
| `--no-roc` | covered | 3 |  |
| `--roc-regions` | covered | 1 |  |
| `--roc-filter` | covered | 1 |  |
| `--roc-delta` | covered | 1 |  |
| `--ci-alpha` | covered | 1 |  |
| `--no-json` | covered | 2 |  |
| `--location` | covered | 9 |  |
| `--pass-only` | covered | 3 |  |
| `--filters-only` | covered | 1 |  |
| `--restrict-regions` | covered | 1 |  |
| `--target-regions` | covered | 1 |  |
| `--leftshift` | gap | 0 | only the `--no-leftshift` side is bound here; `hap pre` binds `-L` |
| `--no-leftshift` | covered | 7 |  |
| `--decompose` | covered | 1 |  |
| `--no-decompose` | covered | 7 |  |
| `--bcftools-norm` | covered | 1 |  |
| `--fixchr` | covered | 4 |  |
| `--no-fixchr` | gap | 0 | `hap pre` binds it; germline binds only `--fixchr` |
| `--bcf` | covered | 2 |  |
| `--somatic` | covered | 1 |  |
| `--set-gt` | covered | 2 | values bound: `half`, `hom`; not bound: `hemi`, `het`, `first` |
| `--filter-nonref` | covered | 1 |  |
| `--convert-gvcf-to-vcf` | covered | 1 |  |
| `--gender` | covered | 4 | values bound: `male`, `none`; not bound: `female`, `auto` |
| `--convert-gvcf-truth` | covered | 1 |  |
| `--convert-gvcf-query` | covered | 1 |  |
| `--preprocess-truth` | covered | 6 |  |
| `--usefiltered-truth` | covered | 2 |  |
| `--preprocessing-window-size` | covered | 1 |  |
| `--adjust-conf-regions` | gap | 0 | only `--no-adjust-conf-regions` is bound; `hap quantify` binds the positive form |
| `--no-adjust-conf-regions` | covered | 5 |  |
| `--unhappy` | covered | 1 |  |
| `--no-haplotype-comparison` | gap | 0 | the `--unhappy` spelling of the same flag is bound |
| `--window-size` | covered | 1 |  |
| `--xcmp-enumeration-threshold` | covered | 1 |  |
| `--xcmp-expand-hapblocks` | covered | 1 |  |
| `--threads` | covered | 34 | bound by the harness on every row |
| `--engine` | covered | 8 | values bound: `vcfeval` 5, `scmp-distance` 2, `scmp-somatic` 1; `xcmp` runs as the default on the other 26 rows and is never spelled out |
| `--engine-vcfeval-path` | gap | 0 | was bound by the `matrix_vcfeval_deprecated_flags` contract case, removed in `1a1bbca` |
| `--engine-vcfeval-template` | covered | 5 | legacy side only; hap-rs ignores it, and `/final_args/engine_vcfeval_template` is the exemption register's first entry |
| `--scmp-distance` | covered | 6 |  |
| `--lose-match-distance` | covered | 1 |  |
| `--logfile` | absent by decision | 0 | writes outside the observed `result*` prefix set (ADR 0007) |
| `--verbose` | absent by decision | 0 | log level, stream content is outside the claim (ADR 0003) |
| `--quiet` | absent by decision | 0 | log level, stream content is outside the claim (ADR 0003) |

#### `hap somatic`, 29 of 41 options bound

| Option | Status | Rows | Note |
|---|---|---:|---|
| `--help` | absent by decision | 0 | standard output content is outside the claim (ADR 0003) |
| `--output` | covered | 27 | bound by the harness on every row |
| `--location` | covered | 1 |  |
| `--restrict-regions` | covered | 1 |  |
| `--target-regions` | covered | 1 |  |
| `--false-positives` | covered | 27 | bound by the harness on every row |
| `--ambiguous` | covered | 3 |  |
| `--ambi-fp` | covered | 2 |  |
| `--no-ambi-fp` | covered | 2 |  |
| `--count-unk` | covered | 6 |  |
| `--no-count-unk` | covered | 2 |  |
| `--explain_ambiguous` | covered | 2 |  |
| `--reference` | covered | 27 | bound by the harness on every row |
| `--scratch-prefix` | absent by decision | 0 | scratch placement is a filesystem side effect outside the claim (ADR 0005) |
| `--keep-scratch` | absent by decision | 0 | scratch placement is a filesystem side effect outside the claim (ADR 0005) |
| `--continue` | gap | 0 | it reuses scratch VCFs from a previous run; ADR 0005 records scratch placement, not scratch reuse, so no decision covers it |
| `--include-nonpass` | covered | 9 |  |
| `--feature-table` | covered | 27 | values bound: `admix.strelka.snv`, `generic`, `hcc.strelka.indel`; not bound: `hcc.strelka.snv`, `hcc.pisces.snv`, `hcc.mutect.snv`, `hcc.varscan2.indel`, `hcc.pisces.indel`, `admix.strelka.indel`, `hcc.varscan2.snv`, `hcc.mutect.indel` |
| `--happy-stats` | gap | 0 | writes `summary.csv`, an artifact no committed row observes |
| `--bam` | absent by decision | 0 | no legacy reference, the four `--bam` rows left the truth set (ADR 0002) |
| `--normalize-truth` | covered | 1 |  |
| `--normalize-query` | covered | 1 |  |
| `--normalize-all` | covered | 1 |  |
| `--fixchr-truth` | gap | 0 | the `--fix-chr-truth` alias is bound, this spelling is not |
| `--fixchr-query` | gap | 0 | the `--fix-chr-query` alias is bound, this spelling is not |
| `--fix-chr-truth` | covered | 1 |  |
| `--fix-chr-query` | covered | 1 |  |
| `--no-fixchr-truth` | covered | 2 |  |
| `--no-fixchr-query` | covered | 2 |  |
| `--no-order-check` | gap | 0 | legacy help calls it a dev feature that disables an internal order check |
| `--roc` | covered | 5 | values bound: `strelka.snv.qss`, `mutect.snv`, `strelka.snv`, `strelka.indel.evs`, `varscan2.snv`; not bound: `varscan2.indel`, `mutect.indel`, `strelka.indel`, `strelka.snv.vqsr` |
| `--bin-afs` | covered | 3 |  |
| `--af-binsize` | covered | 1 |  |
| `--af-truth` | covered | 1 |  |
| `--af-query` | covered | 1 |  |
| `--count-filtered-fn` | covered | 1 |  |
| `--fp-region-size` | covered | 5 |  |
| `--ci-level` | covered | 1 |  |
| `--logfile` | absent by decision | 0 | writes outside the observed `result*` prefix set (ADR 0007) |
| `--verbose` | absent by decision | 0 | log level, stream content is outside the claim (ADR 0003) |
| `--quiet` | absent by decision | 0 | log level, stream content is outside the claim (ADR 0003) |

#### `hap pre`, 22 of 26 options bound

| Option | Status | Rows | Note |
|---|---|---:|---|
| `--help` | absent by decision | 0 | standard output content is outside the claim (ADR 0003) |
| `--location` | covered | 2 |  |
| `--pass-only` | covered | 2 |  |
| `--filters-only` | covered | 3 |  |
| `--restrict-regions` | covered | 47 | bound by the harness on every row |
| `--target-regions` | covered | 1 |  |
| `--leftshift` | covered | 2 |  |
| `--no-leftshift` | covered | 12 |  |
| `--decompose` | covered | 2 |  |
| `--no-decompose` | covered | 12 |  |
| `--bcftools-norm` | covered | 2 |  |
| `--fixchr` | covered | 2 |  |
| `--no-fixchr` | covered | 3 |  |
| `--bcf` | covered | 1 |  |
| `--somatic` | covered | 7 |  |
| `--set-gt` | covered | 6 | values bound: `half`, `hemi`, `het`, `hom`, `first` |
| `--filter-nonref` | gap | 0 | `hap germline` binds it |
| `--convert-gvcf-to-vcf` | covered | 2 |  |
| `--gender` | covered | 21 | values bound: `male`, `female`, `auto`, `none` |
| `--version` | absent by decision | 0 | prints a version and exits 0 on every subcommand, no comparison possible (ADR 0003 correction) |
| `--reference` | covered | 47 | bound by the harness on every row |
| `--window-size` | covered | 5 |  |
| `--threads` | covered | 47 | bound by the harness on every row |
| `--logfile` | covered | 1 |  |
| `--verbose` | covered | 1 |  |
| `--quiet` | absent by decision | 0 | log level, stream content is outside the claim (ADR 0003) |

#### `hap ftx`, 10 of 12 options bound

| Option | Status | Rows | Note |
|---|---|---:|---|
| `--help` | absent by decision | 0 | standard output content is outside the claim (ADR 0003) |
| `--output` | covered | 25 | bound by the harness on every row |
| `--location` | covered | 2 |  |
| `--restrict-regions` | covered | 1 |  |
| `--target-regions` | covered | 1 |  |
| `--include-nonpass` | covered | 2 |  |
| `--feature-table` | covered | 25 | values bound: `generic`, `hcc.strelka.snv`, `hcc.strelka.indel`, `admix.strelka.snv`, `admix.strelka.indel`, `hcc.mutect.snv`, `hcc.mutect.indel`, `hcc.varscan2.snv`, `hcc.varscan2.indel`, `hcc.pisces.snv`, `hcc.pisces.indel` |
| `--feature-label` | covered | 17 |  |
| `--bam` | absent by decision | 0 | no legacy reference, the four `--bam` rows left the truth set (ADR 0002) |
| `--reference` | covered | 25 | bound by the harness on every row |
| `--normalize` | covered | 2 |  |
| `--fix-chr` | covered | 2 |  |

#### `hap quantify`, 23 of 27 options bound

| Option | Status | Rows | Note |
|---|---|---:|---|
| `--help` | absent by decision | 0 | standard output content is outside the claim (ADR 0003) |
| `--version` | absent by decision | 0 | prints a version and exits 0 on every subcommand, no comparison possible (ADR 0003 correction) |
| `--adjust-conf-regions` | covered | 1 |  |
| `--type` | covered | 9 | values bound: `xcmp`, `ga4gh` |
| `--false-positives` | covered | 9 | bound by the harness on every row |
| `--stratification` | covered | 1 |  |
| `--stratification-region` | covered | 1 |  |
| `--stratification-fixchr` | covered | 1 |  |
| `--write-vcf` | covered | 2 |  |
| `--write-counts` | covered | 2 |  |
| `--no-write-counts` | covered | 3 |  |
| `--output-vtc` | covered | 1 |  |
| `--preserve-info` | gap | 0 | `hap germline` binds it |
| `--roc` | covered | 3 |  |
| `--no-roc` | covered | 5 |  |
| `--roc-regions` | covered | 1 |  |
| `--roc-filter` | covered | 1 |  |
| `--roc-delta` | covered | 1 |  |
| `--ci-alpha` | covered | 1 |  |
| `--no-json` | covered | 3 |  |
| `--report-prefix` | covered | 10 | bound by the harness on every row, and again by `chr21_counts_last_on` |
| `--reference` | covered | 9 | bound by the harness on every row |
| `--threads` | covered | 10 | bound by the harness on every row, and again by `chr21_ga4gh_artifacts` |
| `--logfile` | covered | 1 |  |
| `--bcf` | covered | 1 |  |
| `--verbose` | covered | 1 |  |
| `--quiet` | absent by decision | 0 | log level, stream content is outside the claim (ADR 0003) |

#### `hap validate`, 9 of 11 options bound

| Option | Status | Rows | Note |
|---|---|---:|---|
| `--help` | absent by decision | 0 | standard output content is outside the claim (ADR 0003) |
| `--version` | absent by decision | 0 | prints a version and exits 0 on every subcommand, no comparison possible (ADR 0003 correction) |
| `--input-file` | covered | 1 |  |
| `--output-file` | covered | 12 | bound by the harness on every row |
| `--location` | covered | 1 |  |
| `--limit-records` | covered | 3 |  |
| `--message-every` | covered | 1 |  |
| `--apply-filters` | covered | 2 |  |
| `--strict-homref` | covered | 2 |  |
| `--check-bcf-errors` | covered | 8 |  |
| `--all-warnings` | covered | 2 |  |

#### Values the pinned parsers enumerate

| Command | Option | Values bound | Values not bound |
|---|---|---|---|
| `hap germline` | `--engine` | `xcmp` as the default on 26 rows, `vcfeval` 5, `scmp-distance` 2, `scmp-somatic` 1 | — |
| `hap germline` | `--set-gt` | `half`, `hom` | `hemi`, `het`, `first`, all five bound on `hap pre` |
| `hap germline` | `--gender` | `male`, `none` | `female`, `auto`, both bound on `hap pre` |
| `hap germline` | `--type` | — | `xcmp`, `ga4gh`, both bound on `hap quantify` |
| `hap pre` | `--set-gt` | `half`, `hemi`, `het`, `hom`, `first` | — |
| `hap pre` | `--gender` | `male`, `female`, `auto`, `none` | — |
| `hap quantify` | `--type` | `xcmp`, `ga4gh` | — |
| `hap ftx` | `--feature-table` | all 11 | — |
| `hap somatic` | `--feature-table` | `generic`, `hcc.strelka.indel`, `admix.strelka.snv` | `hcc.strelka.snv`, `admix.strelka.indel`, `hcc.mutect.snv`, `hcc.mutect.indel`, `hcc.varscan2.snv`, `hcc.varscan2.indel`, `hcc.pisces.snv`, `hcc.pisces.indel` |
| `hap somatic` | `--roc` | `strelka.snv`, `strelka.snv.qss`, `strelka.indel.evs`, `mutect.snv`, `varscan2.snv` | `strelka.indel`, `strelka.snv.vqsr`, `mutect.indel`, `varscan2.indel` |
| `hap somatic` | `--fp-region-size` | a nucleotide count, `30`, and `auto` | none; the default derivation runs on every other row |
| `hap somatic` | `--af-binsize` | one bin size, `0.5` | the comma-separated multi-size form the help documents |
| `hap validate` | `--limit-records` | `4`, `-1`, `-2` | none |

The harness binds `--feature-table` on every SOMPY row, so what is unbound there
is eight of its values. FTXPY binds all eleven tables, because FTXPY extracts
them and SOMPY consumes them.

#### Input format and index kind

| Form | Status | Where |
|---|---|---|
| plain `.vcf` | covered | all six lanes |
| `.vcf.gz` with a staged `.tbi` | covered | HAPPY, PREPY, QFY, VCFCHECK |
| `.vcf.gz` with no staged index | covered | SOMPY, FTXPY, which `main.nf` gives no index companion |
| `.bcf` with a staged `.csi` | covered | HAPPY `matrix_implicit_bcf`, PREPY `matrix_bcf_input` |
| plain `.g.vcf` | covered | HAPPY `matrix_gvcf_conversion`, `matrix_gvcf_conversion_both`, PREPY `tiny_gvcf_edge`, `tiny_convert_gvcf` |
| reference `.fa` with `.fai` | covered | all six lanes |
| reference `.fa.gz` with `.fai` and `.gzi` | gap in the gate | bound only by the public germline row, which runs outside nf-test |
| plain `.bed` | covered | HAPPY, PREPY, FTXPY, SOMPY, QFY |
| `.bed.gz` with `.tbi` | covered | HAPPY, PREPY, QFY |
| stratification `.tsv` | covered | HAPPY `matrix_stratification_tsv_*`, QFY `chr21_tsv_fixchr` |
| RTG SDF `.tar.gz` | covered on the legacy side | the five HAPPY `--engine vcfeval` rows; hap-rs takes the FASTA |
| `.bam` with `.bai` | absent by decision | no legacy reference, ADR 0002 |

#### Filtering mode

| Mode | Status | Where |
|---|---|---|
| default, all filters passed on | covered | every lane |
| `--pass-only` | covered | HAPPY 3, PREPY 2 |
| `-P` / `--include-nonpass` | covered | SOMPY 9, FTXPY 2 |
| `--filters-only` | covered | HAPPY 1, PREPY 3 |
| `--usefiltered-truth` | covered | HAPPY 2 |
| `--filter-nonref` | covered on HAPPY 1, gap on PREPY | |
| `--apply-filters true` and `false` | covered | VCFCHECK `warning_controls`, `explicit_boolean_values` |
| Boost boolean spellings `on`/`off`, `yes`/`no`, `1`/`0` | absent by decision | hap-rs takes `true`/`false` only, ADR 0003 |

#### Region and stratification

| Form | Status | Where |
|---|---|---|
| `-l` / `--location`, single range | covered | HAPPY 9, PREPY 1, SOMPY 1, FTXPY 2, VCFCHECK 1 |
| `-l` with a comma-separated list | covered | PREPY `matrix_overlapping_locations` |
| `-R` / `--restrict-regions` | covered | HAPPY 1, SOMPY 1, FTXPY 1, and every PREPY row through the harness |
| `-T` / `--target-regions` | covered | HAPPY 1, SOMPY 1, PREPY 1, FTXPY 1 |
| `--stratification` TSV | covered | HAPPY 2, QFY 1 |
| `--stratification-region NAME:bed` | covered | HAPPY 1, QFY 1 |
| `--stratification-fixchr` | covered | HAPPY 2, QFY 1 |
| more than one BED in one table | covered | HAPPY `matrix_stratification_tsv_multiple_beds` |
| four-column stratification BED | covered | QFY `chr21_four_column_strat` |
| `--roc-regions` | covered | HAPPY 1, QFY 1 |
| `--adjust-conf-regions` | covered on QFY 1, gap on HAPPY | |

#### Options outside the 179

Two rows bind `--force-interactive`, which the pinned parsers do not register
because the image has no SGE, and which legacy accepts anyway: PREPY
`tiny_verbose_log` and QFY `chr21_ga4gh_artifacts`. ADR 0003 records the same
observation.

## Corpus provenance

Every file the 155 committed rows consume, including the index companions
`main.nf` stages alongside each input: 19 upstream files and 123 repository
fixtures, 142 in total. Sizes and SHA-256 measured on 2026-08-18, the upstream
ones by fetching each pinned URL and the fixtures over the working tree at
`18548d1`.

Upstream files come from `params.fixture_base`, the raw GitHub permalink pinned
to Illumina/hap.py `8401169`. Prefix each path below with
`https://raw.githubusercontent.com/Illumina/hap.py/84011695b2ff2406c16a335106db6831fb67fdfe/`
to get its primary source URL.

| File | Bytes | SHA-256 | Lanes |
|---|---:|---|---|
| `example/chr21.fa` | 49092500 | `33f1357172a6265bb6fb82c7368e7e7efca11c303ca2e6cc25d70f6b44c5821a` | FTXPY HAPPY PREPY QFY SOMPY |
| `example/chr21.fa.fai` | 23 | `f7462bb56b94790653b357950b0f7d1d3267d54fbfc419e42fb2494c3fb5b670` | FTXPY HAPPY PREPY QFY SOMPY |
| `example/happy/NA12878_chr21.vcf.gz` | 7342312 | `93ba1801c0a65ecf20eaf73e9383636d3fda5a3d013caa9b2037aad25ab3f38e` | HAPPY PREPY QFY VCFCHECK |
| `example/happy/NA12878_chr21.vcf.gz.tbi` | 20641 | `1b188241ed2a0c8aca8726b7605b43f68419c072ec18afba4d73c198117c2bdf` | HAPPY PREPY QFY VCFCHECK |
| `example/happy/PG_Conf_chr21.bed.gz` | 801211 | `bc94994a5cac62c2f88c306f2e2c1c8d9d21772438f8494fdd6096c2880ebc87` | HAPPY PREPY QFY |
| `example/happy/PG_Conf_chr21.bed.gz.tbi` | 6196 | `9c4954f01959597797e96ab171b7aafa0eaf031e352064808d39940753ff893f` | HAPPY PREPY QFY |
| `example/happy/PG_NA12878_chr21.vcf.gz` | 484628 | `a69c13262fbe21c690bf8c3ab33977bd0257c70777a50d6a79bcfac2e88298e7` | HAPPY QFY |
| `example/happy/PG_NA12878_chr21.vcf.gz.tbi` | 18953 | `5083646e1315ee3fb3be51970a7af504ba2125d26bb8e589520ed9f3c239c8ab` | HAPPY QFY |
| `example/happy/hg38.chr21.fa` | 47488490 | `634829fe7eba3074f481818a448f7e504e19bdcb992d3388aad4de034739e0b2` | FTXPY SOMPY |
| `example/happy/hg38.chr21.fa.fai` | 23 | `2e601331315e6e04fc6c56b6ccce82676cb829b1025f87111d2932605101a420` | FTXPY SOMPY |
| `example/sompy/FP_admix.bed.gz` | 800141 | `9228d4e6eac52319e90d0399a866759947ac144f4628963f966a9d6b37e1128e` | SOMPY |
| `example/sompy/FP_admix_grch38.bed.gz` | 620976 | `f8b60bba09a19ede10a50be2a10df97b3f185f5d9c4ad3ddc0731683f144c389` | SOMPY |
| `example/sompy/PG_admix_truth_grch38_indels.vcf.gz` | 30352 | `5c74f1a406e15c44d698630f3226203fd41af302032fb086d29b87e678f1d039` | SOMPY |
| `example/sompy/PG_admix_truth_grch38_snvs.vcf.gz` | 139504 | `147463137164e775d017ead605d573572adbe960dc3c9110a2b55c12068bea9f` | SOMPY |
| `example/sompy/PG_admix_truth_snvs.vcf.gz` | 179144 | `87938d32dcbc9cc61cb17339f10a513b98c775762f22838840f128de23852d11` | SOMPY |
| `example/sompy/strelka_admix_snvs.vcf.gz` | 2017897 | `09e696b64dd09c0f8883b24bc3313b860d0137668f273eb766dd4fa40d77520c` | FTXPY SOMPY |
| `example/sompy/strelka_grch38_admix_pass_indels.vcf.gz` | 167545 | `69352869ca46b8f2ba6e6648a54fe317c85426123a7f21d0e05885b36ab1f73e` | FTXPY SOMPY |
| `example/sompy/strelka_grch38_admix_pass_snvs.vcf.gz` | 1042293 | `f33cf3888b6d2f8af58afac8c481a327d9722f82757a91a81d31be895ef126e9` | FTXPY SOMPY |
| `src/data/test_formats.vcf` | 2129 | `bf07a880c7d6ab9e98835ba00ed7a6516e43b740f787df3ecd7a8fe8d81757bb` | VCFCHECK |


Repository fixtures:

| File | Bytes | SHA-256 | Added in | Lanes |
|---|---:|---|---|---|
| `verification/assets/fixtures/external/ebi-vcf-validator/failed_body_samples_ploidy_000.vcf` | 221 | `92271c32027f458fb3a89fe13b11210d3a219beca50329327f4d243776e2c209` | `db87278` | VCFCHECK |
| `verification/assets/fixtures/external/ebi-vcf-validator/passed_meta_format_P_1.vcf` | 1291 | `7a94657671a660a4d4525f6d01234d9d17ae0414cc367546617863180fcb7aba` | `db87278` | VCFCHECK |
| `verification/assets/fixtures/ftx-bam/input.vcf` | 904 | `eecab87b827cf7b289176d3c366c452d4132ba7c86d057421f64e43d35213398` | `43076b7` | SOMPY |
| `verification/assets/fixtures/ftx-bam/mutect.vcf` | 650 | `378b4eab515a8bf8fe0af71d65e6911a636caa41c240db829a80b4f277bc119a` | `43076b7` | FTXPY SOMPY |
| `verification/assets/fixtures/ftx-bam/normalize.vcf` | 114 | `b0562c5ad20d3945c9bfa9e5b91af4ef62a39505220277880356b307f38c3474` | `43076b7` | FTXPY |
| `verification/assets/fixtures/ftx-bam/options.vcf` | 148 | `f90bf5ee4c33bd2e9a86fb17c22126c4e931ad56fde3a9d215e68e42501b8fcf` | `43076b7` | FTXPY |
| `verification/assets/fixtures/ftx-bam/pisces.vcf` | 960 | `aac838d5999c6f7b6bbac5667ca17c526290a3feecd4c9012e5cea5116b1229c` | `43076b7` | FTXPY |
| `verification/assets/fixtures/ftx-bam/query.vcf` | 904 | `eecab87b827cf7b289176d3c366c452d4132ba7c86d057421f64e43d35213398` | `43076b7` | SOMPY |
| `verification/assets/fixtures/ftx-bam/ref.fa` | 107 | `40d53f8f5ba8ef031a24edaa98e02d57c347e483e57631e458389bcacacf973f` | `43076b7` | FTXPY SOMPY |
| `verification/assets/fixtures/ftx-bam/ref.fa.fai` | 19 | `af4a4a72b3699b6ea26cff47b1eb124956b693bcecc3adf0ad4cda4ac808152a` | `43076b7` | FTXPY SOMPY |
| `verification/assets/fixtures/ftx-bam/varscan2.vcf` | 815 | `aad7fc60b18e0201e83d4f6ff1ed063fdf03b8b381186035fae0c0c0b1540901` | `43076b7` | FTXPY SOMPY |
| `verification/assets/fixtures/ftx-matrix/empty.vcf` | 90 | `0243db9494c14f47986d93750709023220ad490ba12e08c3f8910a5e59b3644d` | `eb37118` | FTXPY |
| `verification/assets/fixtures/ftx-matrix/mutect-indel.vcf` | 803 | `9c3f3b3d08e1401a246b14538c58ac7f0d4cc9d0a0b15bcb0641f1cdc22dc8f5` | `eb37118` | FTXPY |
| `verification/assets/fixtures/ftx-matrix/mutect-reversed.vcf` | 788 | `60f0d6a431b2ac9bdaacc4c593d25a217036914062095053cfd2cea5becb2a0a` | `eb37118` | FTXPY |
| `verification/assets/fixtures/ftx-matrix/pisces-indel.vcf` | 962 | `45db77409f4d2c1d4fec64d6aed683455e46b8e0644b6753b718284cb2976661` | `eb37118` | FTXPY |
| `verification/assets/fixtures/ftx-matrix/ref.fa` | 107 | `40d53f8f5ba8ef031a24edaa98e02d57c347e483e57631e458389bcacacf973f` | `eb37118` | FTXPY SOMPY |
| `verification/assets/fixtures/ftx-matrix/ref.fa.fai` | 19 | `af4a4a72b3699b6ea26cff47b1eb124956b693bcecc3adf0ad4cda4ac808152a` | `eb37118` | FTXPY SOMPY |
| `verification/assets/fixtures/ftx-matrix/selector-normalize.vcf` | 238 | `0c2988453b941ce1a231e969a0ee5f04fd51f594ae05e0cb268dff01ca833a26` | `eb37118` | FTXPY |
| `verification/assets/fixtures/ftx-matrix/varscan2-indel.vcf` | 817 | `835c2d234f01e4c19b7e25d8f3ec4ee0ac6d94386ef2b46d5cd89580dbe9b762` | `eb37118` | FTXPY |
| `verification/assets/fixtures/ftx-nonref.vcf` | 284 | `f719fa9d22cde75c18b27671c58f1fa334f572366356b49d9dbc8b714cf3ff8a` | `43076b7` | FTXPY |
| `verification/assets/fixtures/germline-output-parity/confident.bed` | 11 | `904f3dfd384dc3979e5ec96a359b3fb297209d90b3e53bd386a303c64b5bc614` | `4bd5133` | HAPPY |
| `verification/assets/fixtures/germline-output-parity/query.vcf` | 762 | `f81966c5ba3db1597552251a49b9588b9c10f9243c86d9fa596baa81ed431cc0` | `4bd5133` | HAPPY |
| `verification/assets/fixtures/germline-output-parity/ref.fa` | 107 | `f4eb29593a5203a0f403a39160172b7292f9245c8173f54399ab6e0dbbe93176` | `4bd5133` | HAPPY |
| `verification/assets/fixtures/germline-output-parity/ref.fa.fai` | 19 | `af4a4a72b3699b6ea26cff47b1eb124956b693bcecc3adf0ad4cda4ac808152a` | `4bd5133` | HAPPY |
| `verification/assets/fixtures/germline-output-parity/truth.vcf` | 757 | `3607f0729dac22f5a5414538ac303aebe4555e66a3c1769033d965a25854def5` | `4bd5133` | HAPPY |
| `verification/assets/fixtures/happy-matrix/confident.bed` | 30 | `0e274c8294e7290ac23582be501091fecc67b8746f3b83df77d00103c47ed4d4` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/numeric-query.vcf` | 468 | `b404c478e0114c20fd35c3b919fef2ebb76e0eab0743d686c52a63afac47c038` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/numeric-truth.vcf` | 376 | `8a1f4e451af08aa2b1b659bb71e89c9a1c7600ab8356083cd30f75dd6881621c` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/query.g.vcf` | 905 | `b277c9b29c511ab0cfc1f129805078ccb9f7f211ae76cfca803922690b79be52` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/query.vcf` | 482 | `1e7e3f42c8c55f3db9720f6835df6d96810b741831f765287854e710bc36b70e` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/ref.fa` | 154 | `3d3735eb3fda0fd7c05bb610132e5f32627286c108e4976134a3f6613a7393b2` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/ref.fa.fai` | 37 | `12a49d4da0ec7bf5a27501066159d83ffede50fff036094db8bf27bf0245c67a` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/truth.g.vcf` | 736 | `62662a96f0ee715589f1ebb434cfa3d04a9ef3a3d92e9738d5c4f80039d1bc68` | `eb37118` | HAPPY |
| `verification/assets/fixtures/happy-matrix/truth.vcf` | 479 | `fdee10a519beab64d2280198a6e7600ba17949e6e29e10c18d259fbf35801b07` | `eb37118` | HAPPY |
| `verification/assets/fixtures/hg001-graph-record-order/confident.bed` | 10 | `15b3f2eb412a054ce887cfd8f2de7727944152649299167ffa4387e5c76f0b1c` | `67471b5` | HAPPY |
| `verification/assets/fixtures/hg001-graph-record-order/query.vcf` | 255 | `4f65ad05046bc40a48099603568db9ceee7c0d8fdf0c693253c0c0bf2078999b` | `67471b5` | HAPPY |
| `verification/assets/fixtures/hg001-graph-record-order/truth.vcf` | 295 | `993e1948d6d0b4aba51d33aac079c010f56cc350f3f5d7101fbc7d83473907b0` | `67471b5` | HAPPY |
| `verification/assets/fixtures/hg001-primitive-identity/input.vcf` | 319 | `a241aa45ea8c5fb997e7be48bbda1dc7a1d35019f7efcd8eaf39c6f97718e986` | `77e0eb7` | PREPY |
| `verification/assets/fixtures/hg001-primitive-identity/ref.fa` | 27 | `0d9751460b42389148cb92fa989bb6a44e5f558aba33fa07cf9acb4d7d71c367` | `77e0eb7` | PREPY |
| `verification/assets/fixtures/hg001-primitive-identity/ref.fa.fai` | 16 | `622d4f277daffb94dc5110b88c1ef81917f77a3bcbbc4a3f3094af5b25c70dac` | `77e0eb7` | PREPY |
| `verification/assets/fixtures/hg001-primitive-identity/regions.bed` | 10 | `21ee040d280b96cd1c760fe44e987cdc1162875842b415bd206e8a54c538d2aa` | `77e0eb7` | PREPY |
| `verification/assets/fixtures/pre-boundary-symbolic.vcf` | 314 | `805babdd904ab73c8991b4d8064bc605232e06228533ed784eea93bdeb010079` | `50adc1c` | PREPY |
| `verification/assets/fixtures/pre-boundary.vcf` | 382 | `13483b2814031b9aba8be95ea7157d66862646f6829fab7d9994ff07f79179d8` | `50adc1c` | PREPY |
| `verification/assets/fixtures/pre-breakend-forms.vcf` | 533 | `9441e898c56659e5388112c4a714d8326d370909a38c8fc001a327a2721ceaa1` | `50adc1c` | PREPY |
| `verification/assets/fixtures/pre-edge-alleles.vcf` | 498 | `4f64c9fbac18e760870a1ed47e882b3d555d86ae6892273d1ceb7a223bc4a144` | `9e958b8` | PREPY |
| `verification/assets/fixtures/pre-gender-half.vcf` | 274 | `26f4ed96d6eed30fa887340b94f4ca1ad88dae50d130a603e4d89a5b3c7ec08b` | `43076b7` | PREPY |
| `verification/assets/fixtures/pre-gender.bed` | 9 | `c8e769503404c195e4ff1f7b1835dac2d15a2884401cb24acdfa92e86b9a3d13` | `43076b7` | HAPPY PREPY |
| `verification/assets/fixtures/pre-gender.fa` | 12 | `71df3da9ebf0f667eacdbcd1bd64553406785b2454627ddf7c6e1972fd3bbe51` | `43076b7` | HAPPY PREPY |
| `verification/assets/fixtures/pre-gender.fa.fai` | 13 | `9e231aad92943581abc55d7c4ddb0cb9c87a7bc32ba9ee477164c5bc31451228` | `43076b7` | HAPPY PREPY |
| `verification/assets/fixtures/pre-gender.vcf` | 268 | `607cda867e489fbe46765c752508db691e6645439e14b07573971efd81f914eb` | `43076b7` | HAPPY PREPY |
| `verification/assets/fixtures/pre-gvcf-edge.vcf` | 756 | `0b4b56f605eb2ee9d4803610d5be4d94c9935f833c39811dec6f8a118cfa1f82` | `50adc1c` | PREPY |
| `verification/assets/fixtures/pre-matrix/all.bed` | 27 | `a236c3b7eaa515d4805b66c434327852bca2d3da1164aab1f9e869420b48c69e` | `eb37118` | HAPPY PREPY |
| `verification/assets/fixtures/pre-matrix/end-boundary.bed` | 9 | `c22efdb8696e3ac49b4dffc00e2ca277b1e029f1b5bacfff4425bcfc4e92697a` | `eb37118` | FTXPY PREPY |
| `verification/assets/fixtures/pre-matrix/end-span.vcf` | 310 | `19199d96f263a71b0775e19eea7a3cf12a61c517abef2257da790666c8b298ed` | `eb37118` | FTXPY PREPY |
| `verification/assets/fixtures/pre-matrix/filters.vcf` | 340 | `332b334268ab96523dd422e2ac3e531baaedc188cf5276274fb7da0849e20bf8` | `eb37118` | PREPY |
| `verification/assets/fixtures/pre-matrix/fixchr.vcf` | 224 | `d533cd866e97509b7f01a7b952b99df4108f45b0e4bfe635f71a91ee885e9d3a` | `eb37118` | PREPY |
| `verification/assets/fixtures/pre-matrix/input.bcf` | 440 | `54b73ae397f3798779d70ef93b3f5949bdcd04698bfa4031dcfc844b1f34c9a6` | `eb37118` | HAPPY PREPY |
| `verification/assets/fixtures/pre-matrix/input.bcf.csi` | 89 | `771b904f2a089403a1c87412e326b8f250693444a7afbabadc0ef6a745fa8e09` | `eb37118` | HAPPY PREPY |
| `verification/assets/fixtures/pre-matrix/reference.fa` | 78 | `664e138c4ae68ac719916a52ce905b7b9c5d44f7246d9cf35180e755c6aaf555` | `eb37118` | FTXPY HAPPY PREPY |
| `verification/assets/fixtures/pre-matrix/reference.fa.fai` | 47 | `c6d71fd4383c2f67e72cb0b17a6abad81777973aa7faaf7bf4c68b7e82007f9a` | `eb37118` | FTXPY HAPPY PREPY |
| `verification/assets/fixtures/pre-matrix/window.vcf` | 245 | `c6d69ad717cff8323a2d5079ca19cd2c2900fec6861f858a8f6ec0034f3ce6d2` | `eb37118` | PREPY |
| `verification/assets/fixtures/pre-multisample-star-phase.vcf` | 307 | `4ae9ee3ca6a65e2eb90be029735279ce93c1a6fb2ef28d6996a10cb6433566c6` | `50adc1c` | PREPY |
| `verification/assets/fixtures/pre-norm.bed` | 9 | `7c0b71d52436cfeb4bad2ee5baa96636a33ea1e135781a93bb3356f794080a72` | `43076b7` | PREPY |
| `verification/assets/fixtures/pre-norm.fa` | 12 | `978871ae35792640e9a441bc7e972888907b0d956e11cc3f4513504411cdc02b` | `43076b7` | PREPY |
| `verification/assets/fixtures/pre-norm.fa.fai` | 13 | `7336abc97d063a73c595f6dfd26e48f0ad080334f4d47b1a5f7fc8ded90ba7d4` | `43076b7` | PREPY |
| `verification/assets/fixtures/pre-norm.vcf` | 274 | `1a14be96d7d5e33d858e77b2b094e536d66f2900a5003cff8449d9ad5e84a1b9` | `43076b7` | PREPY |
| `verification/assets/fixtures/pre-ploidy-format.vcf` | 498 | `6c54d069fc224fc9748b4bad50fdd5cbcef1c1b6c9e644c8db3a000101459255` | `50adc1c` | PREPY |
| `verification/assets/fixtures/pre-somatic.vcf` | 1051 | `4e75611b4ebeabdb7ae0512d8cae5d8a22bd6d7753f33439a7eea1b6b0f30552` | `43076b7` | PREPY |
| `verification/assets/fixtures/pre-symbolic-del.vcf` | 278 | `0bdb4534870750fd5ab5ea98b30dab1d170478b53f4ffdfee69fe4d654f62fae` | `43076b7` | PREPY |
| `verification/assets/fixtures/pre-symbolic-types.vcf` | 574 | `99e4128c0a6cb23edaa555f0ea2aa1f16e00cd75365839a258e91425bcf01079` | `50adc1c` | PREPY |
| `verification/assets/fixtures/qfy-matrix/stratification-four-column.bed` | 67 | `947859fd785aa5bfb025fb888b84bc80dbb0b257036f9149a02ebd9ec53a7d2d` | `db87278` | QFY |
| `verification/assets/fixtures/roc-filtered-truth/confident.bed` | 21 | `0a7c04b781d393bd5540894b313f3a08af62d5e6f8a84914ccb1984702a16f28` | `df3febe` | HAPPY |
| `verification/assets/fixtures/roc-filtered-truth/query.vcf` | 456 | `7de6dc7ef6790717300697f56735a035b8af01a4eb1dc39bf7f1b034b8423e0c` | `df3febe` | HAPPY |
| `verification/assets/fixtures/roc-filtered-truth/truth.vcf` | 551 | `fa56b9683777e5c3f83de414802d5661edb04fa4fc9f0bcdb53d359fec9de232` | `df3febe` | HAPPY |
| `verification/assets/fixtures/somatic-af-type-parity/confident.bed` | 11 | `904f3dfd384dc3979e5ec96a359b3fb297209d90b3e53bd386a303c64b5bc614` | `f662d6a` | SOMPY |
| `verification/assets/fixtures/somatic-af-type-parity/query.vcf` | 852 | `638db51c883e6cf8339a24cc41f5d2fda18231cba5660cadffcab71eedc29db6` | `f662d6a` | SOMPY |
| `verification/assets/fixtures/somatic-af-type-parity/truth.vcf` | 308 | `8760bed4e65b664c60b1b6078043ff19d03960124e211fde67e5ab7144a697c5` | `f662d6a` | SOMPY |
| `verification/assets/fixtures/somatic-output-parity/confident.bed` | 11 | `904f3dfd384dc3979e5ec96a359b3fb297209d90b3e53bd386a303c64b5bc614` | `ac539dc` | SOMPY |
| `verification/assets/fixtures/somatic-output-parity/query.vcf` | 967 | `58db62ee13a62625194b96fc891837da86a422b35367a8e7d82935cde68cd804` | `ac539dc` | SOMPY |
| `verification/assets/fixtures/somatic-output-parity/truth.vcf` | 469 | `45e800f340989d53e8452c5db166887159a680a1e7fe48469e23b45fe083b7c1` | `ac539dc` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/classification.bed` | 88 | `2ba083bcfe2521f466ade7b036229582be8c86a95e6f26793e1510c30d0a3096` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/query-classification.vcf` | 285 | `7c493e151c9f8d1f060fd43fda1281b67cbe3f7ab803f6667ccc3846e4862f8f` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/query-denominator.vcf` | 253 | `949f82118333214de9706d5e107d893ea8ba5c87cddb1c5d42d8eb87bf548fb1` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/query-normalize.vcf` | 195 | `3ad73e9bf18d37e6b83a44f38968cbf633ff986a0c868308f0c1f2045c184f0f` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/query-regions.vcf` | 226 | `c133698a6b5f128427a620bc6694efffde575d731039b2a930ae1e738adb2d62` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/reference.fa` | 64 | `43b5e784ab11935cd3716f2760ad5309c5b9637840921fc94907637acd6b23f9` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/reference.fa.fai` | 33 | `c0e7542eaf7997b579a05495f337f97bc81f5f6258f10db8eeb0cc748699c554` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/selector.bed` | 9 | `293199aeba3ae70c271079ecbe44c073a10d2b728dbc42666e8441c97d6c9d26` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/truth-normalize.vcf` | 218 | `8ab042353a4af513e375be22c43e2c4a25fcf21251303e6f540f8ec0193891b5` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/truth-regions.vcf` | 226 | `c133698a6b5f128427a620bc6694efffde575d731039b2a930ae1e738adb2d62` | `eb37118` | SOMPY |
| `verification/assets/fixtures/somatic-matrix/truth.vcf` | 223 | `77e07a806774293ef56db852a6a66de4a866b29e7f6ca361a0c45a7c76ca5ee9` | `eb37118` | SOMPY |
| `verification/assets/fixtures/stratification/happy-multi/indels.bed` | 8 | `e5aeb0d80043d0e344fe25ad57aaf21b3874135778ee19923f7870107a0a9b78` | `305b900` | HAPPY |
| `verification/assets/fixtures/stratification/happy-multi/regions.tsv` | 50 | `2ff8fd9aef44c779a8f4756b1c0ad3afde560c8d58d109841fafe66cab8b5f56` | `305b900` | HAPPY |
| `verification/assets/fixtures/stratification/happy-multi/substitutions.bed` | 8 | `1ba37212caf699d14ef7756803d2b172449ca7efaae0ef0eef5d178376198f74` | `305b900` | HAPPY |
| `verification/assets/fixtures/stratification/happy/focus.bed` | 7 | `bc0512e1d820bf96cbaf6e40e43b5105c794880ff577e4d2a68f67ff4809f91e` | `95fb889` | HAPPY |
| `verification/assets/fixtures/stratification/happy/regions.tsv` | 16 | `7bb873df79c2be8fd83bb0ba74ff6e68423a3a302bd65e9063242a376487dab7` | `95fb889` | HAPPY |
| `verification/assets/fixtures/stratification/qfy/focus.bed` | 21 | `6a64b61dc00a544c3a70721f2b5050b7e82676cb2af4c200b67ddb57d2a0a896` | `95fb889` | QFY |
| `verification/assets/fixtures/stratification/qfy/regions.tsv` | 16 | `7bb873df79c2be8fd83bb0ba74ff6e68423a3a302bd65e9063242a376487dab7` | `95fb889` | QFY |
| `verification/assets/fixtures/tiny.all.bed` | 10 | `12f42d9f4d53e4e819dc41a7e082e07c980aa8c9eab91a4cee074b3c91a92225` | `43076b7` | PREPY |
| `verification/assets/fixtures/tiny.fa` | 87 | `794e1b939c12d293c591f4cf42537cb5179537061461348944fc157fcf6bc3f9` | `43076b7` | FTXPY PREPY |
| `verification/assets/fixtures/tiny.fa.fai` | 16 | `477662cce672ee3a959427f9343c90d7459a7357bcc6085691015707d7d54b07` | `43076b7` | FTXPY PREPY |
| `verification/assets/fixtures/vcfcheck-matrix/edge-alleles.vcf` | 604 | `540ffe3f71292ab5fd589574586f33cd0f5b8df2935780727648065fa017aaec` | `9e958b8` | VCFCHECK |
| `verification/assets/fixtures/vcfcheck-matrix/lowercase-x.vcf` | 186 | `f89ce741ff00ec22a3e83c3bec40f4a66b8cc2eec5e1c99f9dab91fde8842e28` | `db87278` | VCFCHECK |
| `verification/assets/fixtures/vcfcheck-matrix/missing-contig.vcf` | 165 | `df82df23b35c87ab940c9429b4d141f55f720fc0585620922bb7554ff8672fe1` | `db87278` | VCFCHECK |
| `verification/assets/fixtures/vcfcheck-matrix/ploidy.vcf` | 365 | `2e9855f94d3d8e26a29bf504b554bec12b162654ed9415b4dd7ba545bf14fb63` | `eb37118` | VCFCHECK |
| `verification/assets/fixtures/vcfcheck-matrix/sample-count-mismatch.vcf` | 198 | `ba1afd5a1c9ee28d943c842f346a34671179de14f10ff4612bd60312772ab3f1` | `db87278` | VCFCHECK |
| `verification/assets/fixtures/vcfcheck-matrix/warning-controls.vcf` | 568 | `ef4fc4ba664ec5b649d457cd16080f8443976a77029410ae0ad80428dba70a5a` | `eb37118` | VCFCHECK |
| `verification/assets/fixtures/vcfeval-matrix/confident.bed` | 11 | `6f986b89f783ae885ae507a593af7ec3990551479e2088cec82b5d46391c88b2` | `bcc4bc6` | HAPPY |
| `verification/assets/fixtures/vcfeval-matrix/query.vcf` | 808 | `51eacdaf7088b06712aa92d65efdce8bcbe95ea5ab8257be5fd62c9e040020b1` | `bcc4bc6` | HAPPY |
| `verification/assets/fixtures/vcfeval-matrix/ref.fa` | 128 | `9db66b90f6796e5e87ac9a46fb1755ddf51e2c94d9f1105302b60e30ad9ffd4e` | `bcc4bc6` | HAPPY |
| `verification/assets/fixtures/vcfeval-matrix/ref.fa.fai` | 17 | `adfa0c9e6693c65144c2a983849d30d9359756d0edb6f7ae3556f8d12db2334c` | `bcc4bc6` | HAPPY |
| `verification/assets/fixtures/vcfeval-matrix/ref.sdf.tar.gz` | 2909 | `aacd186572c8d3726eb7ca2b731bb8d9e5df1f995ebd87b262de6045e4734b89` | `bcc4bc6` | HAPPY |
| `verification/assets/fixtures/vcfeval-matrix/truth.vcf` | 575 | `7277195763fd9319c8c0d9c24f8ace550e5028471ffdf638eb373db9ccd20609` | `bcc4bc6` | HAPPY |
| `verification/assets/scmp_distance/confident.bed` | 10 | `e543ebd97816356369b7740fafc2e547c3ad15749c32f1c7d87e13f8e9b46846` | `43076b7` | HAPPY |
| `verification/assets/scmp_distance/confident.bed.gz` | 69 | `d0777e61aadde2c0c1158536cafc6c424d952ce8f9816f62063c673491d832e0` | `43076b7` | HAPPY |
| `verification/assets/scmp_distance/confident.bed.gz.tbi` | 109 | `83a7624bd7abc04f26ccbedef277aefe8d05655e505f00dfceba40adc4738225` | `43076b7` | HAPPY |
| `verification/assets/scmp_distance/query.vcf` | 193 | `feee9eda135f15c4098d1390c6178a50bbee51c59443c6cfa091bb9d74083f34` | `43076b7` | HAPPY |
| `verification/assets/scmp_distance/ref.fa` | 23 | `463ef094b5e5ca7d8499ca2aee48a27a4709e01261f70707faa210c6ffa2c163` | `43076b7` | HAPPY |
| `verification/assets/scmp_distance/ref.fa.fai` | 16 | `18ede293961da55592b03dcd5694483c4f07ace512f000cd50020038ea9deafe` | `43076b7` | HAPPY |
| `verification/assets/scmp_distance/ref.sdf.tar.gz` | 2935 | `e2b82260576fffaf4b7e576d14418c41329bf34a37ccaccfffdd2dee81b14935` | `bcc4bc6` | HAPPY |
| `verification/assets/scmp_distance/truth.vcf` | 193 | `ea8c7b5e1f60540ab415643da8f7d534736d958180c689951aecfab925e05ca0` | `43076b7` | HAPPY |
| `verification/assets/somatic_options/fp.bed` | 0 | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `43076b7` | SOMPY |
| `verification/assets/somatic_options/query.vcf` | 195 | `2cdf4edd0f3a881b0178b9c6c98d6e8dca5db892bd028e757fde7cdaac50a7ab` | `43076b7` | SOMPY |
| `verification/assets/somatic_options/reference.fa` | 17 | `041042472e03d4c6e733672898fbfdefaa218660ae919803ff5d2c2d2e149a0d` | `43076b7` | SOMPY |
| `verification/assets/somatic_options/reference.fa.fai` | 16 | `6a8b0f12b5c8b286031adada3639537ee65ed38ca9bb58ade510b0178f45f522` | `43076b7` | SOMPY |
| `verification/assets/somatic_options/truth.vcf` | 195 | `133009f83214ba4f74d42580e999337f4a41f42ee1d35e5429e13ad7621adfcf` | `43076b7` | SOMPY |

Four groups of fixtures carry provenance beyond the commit that added them, each
recorded in a README beside the files:

- `fixtures/external/ebi-vcf-validator/`: retained byte for byte from
  EBIvariation/vcf-validator at `0aacc4a`, Apache-2.0. The two checksums above
  match `THIRD_PARTY_LICENSES/manifest.json`.
- `scmp_distance/ref.sdf.tar.gz` and `fixtures/vcfeval-matrix/ref.sdf.tar.gz`:
  RTG Tools 3.12.1 SDF bundles. The checksums above match the two fixture
  READMEs.
- `fixtures/germline-output-parity/`, `fixtures/hg001-graph-record-order/` and
  `fixtures/hg001-primitive-identity/`: minimized from the public GIAB HG001
  v4.2.1 and Illumina Platinum Genomes NA12878 comparison, with the source URLs
  in each README.
- `fixtures/somatic-af-type-parity/` and `fixtures/somatic-output-parity/`:
  synthesized from a full-size legacy som.py run of 2026-08-12. Both are
  consumed by a committed SOMPY row.

### Public data

`assets/samplesheet.happy.public.csv` holds one row, `germline_hg001_platinum`,
which runs outside nf-test. Its inputs are fetched from the NCBI GIAB archive at
run time, so they are named by URL rather than checked in. Sizes are the
`Content-Length` the archive reported on 2026-08-18; the checksums are the ones
`assets/fixtures/public-hg001-platinum/README.md` records.

| File | Bytes | SHA-256 |
|---|---:|---|
| [`HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz) | 125932193 | `93bc4c2c696eaf13515ab058caecc064bfed704f85bac7482330ca91bc730daa` |
| [`HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi) | 1589208 | `68a91b78042e66ed9a53ae4069ac7686abf02ea9a1f88da00e0054a7db06dccc` |
| [`HG001_GRCh38_1_22_v4.2.1_benchmark.bed`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed) | 15479939 | `dc3485e60447a3e863c34dfc063c3710ea5ed5efae50855c5197ba62538263b9` |
| [`NA12878.vcf.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/analysis/Illumina_PlatinumGenomes_NA12877_NA12878_09162015/hg38/2.0.1/NA12878/NA12878.vcf.gz) | 45763575 | `b0608c877e4efdc61fc56927100f81194d2119dc4b47f52e0d0ef7ec8701dead` |
| [`NA12878.vcf.gz.tbi`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/analysis/Illumina_PlatinumGenomes_NA12877_NA12878_09162015/hg38/2.0.1/NA12878/NA12878.vcf.gz.tbi) | 1608320 | `81cb0e3b5eb80b88e004985085d0f3457f52fc300b4c9e10749c10701894dc94` |
| [`GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz) | 886344618 | `3c8def6d325c5d1e934b2dd530c4d1709f027677a15859467071fddf3fff2026` |
| [`GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz.fai`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz.fai) | 7804 | `502a1b8fb73ccd53285c28a0f12df90c818b4fe3de1e862ef47c593ef1a0a4b4` |
| [`GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz.gzi`](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/references/GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz.gzi) | 770664 | `aea5e0202577678b72ab195aa66138a7cdcb8408b4a36a7b3740b6edc2253745` |

The truth data is covered by the NIST data-use policy; the Platinum Genomes
callset is public Illumina data redistributed by GIAB. This row is the only
place the gate reads a bgzip-compressed reference, and the only place it reads
whole-genome inputs.

## Replay and resource baseline

The replay and the resource numbers below come from two full gate runs on one
host on 2026-08-18, at `18548d1`, with `hap` built by `cargo build --release`.
The host is an Apple M4 Pro, 12 cores, macOS 26.6.1, running Nextflow 26.04.6.
`hap` runs natively as an arm64 binary, and the legacy and comparator containers
run in a Colima VM that reports `x86_64` on Ubuntu 24.04.4 with 12 CPUs and 31.3
GiB. The top of this file classes that arrangement as a development diagnostic
and `ubuntu-24.04` in CI as the authoritative one. The emulation sits on the
legacy side.

### Replay

Two full runs, `results/replay-1` and `results/replay-2`, fresh each time with no
`-resume`, the same binary and the same pinned digests. The comparison covers
every published file, 1708 in each run, carrying the same names in both: 622
legacy artifacts, 622 hap-rs artifacts, and 464 harness records. The table below
groups the 622 on each side; the harness records follow it.

| Side | Artifact class | Files | Outcome |
|---|---|---:|---|
| hap-rs | `.csv`, `.tsv` | 174 | byte identical |
| hap-rs | `.csv.gz` | 171 | byte identical |
| hap-rs | `.vcf.gz` | 78 | 77 byte identical; 1 identical after dropping the `##bcftools_*` headers it inherits from the QFY seed |
| hap-rs | `.json`, `.json.gz` | 111 | identical after dropping `runInfo`, `timestamp`, `uname` and `environment` |
| hap-rs | `.bcf`, `.csi`, `.tbi` | 88 | 85 byte identical; the QFY `.bcf` decodes to the same records apart from 16 inherited header lines, and both affected indexes report the same `bcftools index --stats` |
| legacy | `.csv`, `.tsv` | 174 | byte identical |
| legacy | `.csv.gz` | 171 | 157 identical after decompression; **14 differ in content** |
| legacy | `.vcf.gz` | 78 | identical after dropping runtime headers |
| legacy | `.json`, `.json.gz` | 111 | 107 identical after dropping the same provenance keys; **4 differ beyond them** |
| legacy | `.bcf`, `.csi`, `.tbi` | 88 | all 5 `.bcf` decode identically after dropping runtime headers; 66 of 83 indexes differ in bytes, downstream of a stream whose gzip header carries the wall clock |

The harness's own files move with the result: `verification.json`, `report.csv`
and four `comparison.json` differ because the failing case set moved, and the
`.jsonl.gz` artifact records differ wherever they carry the sha256 of a legacy
`.gz`.

hap-rs reproduced its own output in every class above. The fields that vary on
both sides are the
ones already classed as provenance: the JSON `runInfo`, `timestamp`, `uname` and
`environment`, where legacy records the container hostname and the task working
directory; the `Date=` and `/tmp/...` values inside `##bcftools_*Command`
headers; the gzip header mtime on legacy's `.gz`; and the index bytes downstream
of a compressed stream that carries one.

Legacy's 18 differences are content. They sit in the HAPPY ROC tables of
four cases, `chr21`, `chr21_region`, `chr21_passonly` and `chr21_xcmp_controls`:

| Artifact | Run 1 rows | Run 2 rows | Rows only in run 1 | Rows only in run 2 |
|---|---:|---:|---:|---:|
| `chr21` `result.roc.Locations.SNP.csv.gz` | 10106 | 10107 | 0 | 1 |
| `chr21_region` `result.roc.all.csv.gz` | 22639 | 22651 | 1 | 13 |
| `chr21_passonly` `result.roc.Locations.INDEL.PASS.csv.gz` | 3700 | 3699 | 1 | 0 |

Legacy retains a different set of QQ threshold rows from one run to the next.
Rows present in both agree cell for cell. One row legacy emitted in run 1 and not
in run 2 begins `2943.000000,*,*,ALL,*,QUAL,...`, with a number where the variant
type belongs. The four `metrics.json.gz` differences are the same rows seen
through `/metrics/2/data/0/values`.

`maxForks = 1` on `HAPPY_LEGACY` was in force for both runs, so serializing the
Happy matrix does not make the legacy ROC report repeatable here.

So the failing case set moves between runs. Run 1 reported `chr21` and
`chr21_region` different; run 2 reported `chr21_passonly` and
`chr21_xcmp_controls`. hap-rs emitted identical bytes for all four in both runs,
so what moved is the reference. Tracked as
[#46](https://github.com/adamrtalbot/hap.py/issues/46). Whether a native amd64
host shows the same is not measured here.

`docs/adr/0008-compare-artifact-bytes-except-where-the-encoding-carries-provenance.md`
carries a correction for this. Its ROC counts were taken on the image `65476ee`
replaced, and its finding that tightening ROC comparison costs nothing does not
hold while the reference moves.

### Resources

From the Nextflow trace of run 1. `realtime` is per task; the lane column sums
the tasks in that lane.

| Lane | Tasks | Legacy wall | hap-rs wall | hap-rs / legacy | Legacy peak RSS, lane maximum | Longest legacy task | Longest hap-rs task |
|---|---:|---:|---:|---:|---:|---:|---:|
| HAPPY | 34 | 2611 s | 111.3 s | 0.043 | 627 MB | 859 s | 41.7 s |
| SOMPY | 27 | 1478 s | 27.1 s | 0.018 | 321 MB | 276 s | 10.3 s |
| PREPY | 47 | 898 s | 105.0 s | 0.117 | 241 MB | 204 s | 29.9 s |
| FTXPY | 25 | 279 s | 2.5 s | 0.009 | 115 MB | 45 s | 0.3 s |
| QFY | 9 | 457 s | 45.0 s | 0.098 | 264 MB | 122 s | 19.6 s |
| VCFCHECK | 12 | 4 s | 0.8 s | 0.214 | 10 MB | 1 s | 0.5 s |
| **total** | **154** | **5727 s** | **292 s** | **0.051** | | | |

| Lane | Legacy wall run 1 | run 2 | hap-rs wall run 1 | run 2 |
|---|---:|---:|---:|---:|
| HAPPY | 2611 s | 2746 s | 111.3 s | 117.2 s |
| SOMPY | 1478 s | 1388 s | 27.1 s | 28.0 s |
| PREPY | 898 s | 958 s | 105.0 s | 111.8 s |
| FTXPY | 279 s | 266 s | 2.5 s | 2.6 s |
| QFY | 457 s | 535 s | 45.0 s | 42.1 s |
| VCFCHECK | 4 s | 3 s | 0.8 s | 0.8 s |

The second table is the same measurement from run 2. Lane wall time moved in
both directions, by 2 to 17 percent on the five substantial lanes and by a
second on VCFCHECK, where the whole lane runs in four. That is host scheduling
rather than either implementation.

Nextflow samples `peak_rss` from `/proc`, so it has a value for the containerized
legacy and comparator tasks and none for the hap-rs tasks, which run on the macOS
host. The hap-rs figures below come from re-running each lane's longest task
under `/usr/bin/time -l`. Both columns name the same case, and both count in
binary megabytes.

| Lane | Case | hap-rs peak RSS | Legacy peak RSS |
|---|---|---:|---:|
| HAPPY | `chr21` | 829 | 627 |
| QFY | `chr21_adjust_conf_regions` | 209 | 264 |
| PREPY | `chr21_no_fixchr` | 109 | 114 |
| SOMPY | `chr21_includenonpass_countunk` | 64 | 270 |
| FTXPY | `grch38_hcc_strelka_snv` | 41 | 115 |
| VCFCHECK | `location_input_option` | 4 | 10 |

hap-rs holds less memory than legacy on five of the six cases. On HAPPY `chr21`
it holds about a third more, and that case is also the lane maximum for legacy.

No threshold is set here. These are the numbers this host produced on
2026-08-18, and a later run is compared against them. The wall-time ratios are
not portable: legacy runs under emulation and hap-rs runs native, so the ratio
carries an emulation factor of unmeasured size and only a same-host comparison
means anything.

## Comparison rules

`modules/diff.nf` contains the comparator. It requires equal artifact sets,
compares ordered text and non-ROC CSV content, compares ROC CSV rows as
multisets with duplicate counts, and compares typed JSON trees after
canonicalizing ROC table row order and generated table indexes. The comparator
records the lane, case, artifact, and location for each difference. It ignores
these runtime and provenance fields:

- JSON version, timestamp, command-line, generated description fields, and the
  `engine_vcfeval_template` argument the native engine ignores.
- CSV columns named `sompyversion` and `sompycmd`.
- VCF runtime headers such as source, date, and bcftools command/version.

Each exclusion applies to all cases. All remaining content must match.
nf-test reads the combined `verification.json` to decide pass or failure. Each
lane also writes a `comparison.json` with its differences.

## Rules for changes

- Keep the named lane processes visible.
- Keep the reference image immutable and pinned by digest, never a mutable tag.
- Never patch legacy to make it run. No interpreter shim, no source patch, no
  library option override. A legacy tool that will not run means the environment
  is wrong. Settings that change how the same code executes or when it is
  scheduled are harness constraints; they are allowed only when named in
  `docs/adr/0004-pin-the-legacy-baseline-to-one-container-identity.md`.
- Put product behavior in Rust and artifact comparison in Nextflow.
- Fix the Rust implementation when a case differs. Keep the case and comparison
  rule intact.
- Add a samplesheet row and run it before editing Rust.

Keep each third-party license notice beside its imported fixture. The gate uses
those fixtures as test inputs.
