---
title: Compatibility Policy
description: Governed hap.py compatibility behavior, deprecations, and removal policy.
---

## The claim

hap-rs is a drop-in replacement for the Legacy implementation.

For any invocation the pinned parsers accept where legacy exits 0,
`hap <subcommand>` exits 0 and produces the same artifacts as the corresponding
legacy tool, agreeing on every data cell, every categorical verdict label, and
record identity and order. The artifacts compared are every file whose name
begins with the output prefix the invocation names, and nothing else. Where
legacy exits non-zero, neither hap-rs's exit status nor its artifacts are
claimed: hap-rs exits 0 for success, 1 for a failure it detects, the argument
parser's status for a usage error, and the ordinary Unix status where the cause
is external. Observed on linux/amd64. Legacy has no build for any other platform,
so no comparison exists elsewhere, and hap-rs behavior on other platforms is
defined by hap-rs and tested against its own expectations.

## Not claimed

Provenance fields: version, timestamp, command line, generated descriptions,
runtime VCF headers. Message text and which stream carries it. Malformed
invocations. A reference supplied through the environment: `hap` takes the
reference as an argument and reads no variable in its place. VCF on standard
input. hap-rs's own added options. `--bam` feature extraction. The somatic
`--continue` flag, which hap-rs does not implement. Platforms other than
linux/amd64.

Runtime and resource use is not a compatibility question. Neither is independent
scientific correctness: this is an agreement claim, and hap-rs reproduces legacy
including where legacy is wrong.

Standard error is not output. hap-rs may emit clearer error messages than legacy,
and doing so is not a divergence.

## Governance classes

| Class | Meaning | Change policy |
|---|---|---|
| **Indefinite** | Observable behavior on which existing results or automation may depend. | Preserve and regression-test. Removal requires the process below. |
| **Legacy-only** | Emulation activated by a historical command, option, input shape, or output contract. | Keep out of domain algorithms; test with a `legacy_only_` name and test the normative path separately. |
| **Deprecated** | Accepted temporarily to support migration. | Emit a warning with the replacement and removal release. |
| **Accidental/fixable** | An implementation difference with no supported compatibility rationale. | Fix with a regression case; do not turn it into a compatibility promise. |
| **No legacy reference** | Legacy accepts the option, but no reliable legacy behavior exists for it, because the original's own unpinned dependency range cannot run it. | hap-rs behavior is normative: define it, test it against `normative_` expectations, and keep the case out of the truth set. Never invent a legacy observation for it. |

## Compatibility inventory

| Behavior | Class and decision | Rationale or source | Regression evidence |
|---|---|---|---|
| Command aliases (`compare`, `preprocess`, `prepy`, `ftxpy`, `qfy`, `vcfcheck`) and historical option spellings | Indefinite, CLI adapter | Existing wrappers invoke the Python entry-point names. | CLI parser tests in `rust/src/cli_compat/cli.rs` |
| argparse-style repeated options and unique long-option abbreviations | Legacy-only, CLI adapter | hap.py uses argparse and last-value-wins parsing. | `repeated_legacy_options_match_each_wrapper_parser`, abbreviation tests in `rust/src/cli_compat/cli.rs` |
| Unknown options and missing required arguments for `germline`, `somatic`, `pre`, `quantify`, `ftx`, `validate` and their aliases use clap's non-zero usage exit | **Accidental/fixable**, fixed; hap-rs behavior is normative | pre.py and qfy.py return success for unknown-option parser errors, and hap-rs used to exit 1 for `germline` after writing a help page to standard output. Malformed invocations sit outside the covered invocation surface, so routing all six through the argument parser narrows nothing. One path for six subcommands also leaves no legacy success code for a new option to inherit. | `normative_unknown_option_uses_clap_usage_exit_for_every_subcommand`, `normative_missing_required_argument_uses_clap_usage_exit_for_every_subcommand`, `normative_legacy_aliases_use_clap_usage_exit` |
| `--engine-vcfeval-path` and `--engine-vcfeval-template` are accepted and ignored | **Exemption register entry 1**; no removal date | The native engine reads the FASTA supplied with `--reference`, so it needs neither an external RTG install nor an SDF bundle and has nothing to do with either path. Use `--engine vcfeval --reference <FASTA>`; the options stay supported for wrappers that pass them. | `vcfeval_ignores_external_runtime_flags` covers accepting and ignoring both options, `germline_warns_on_stderr_for_each_ignored_vcfeval_option` and `germline_without_the_vcfeval_options_emits_no_warning` cover the warning reaching standard error only when an option is passed, and `vcfeval_ignored_options_warning_names_the_reference_replacement` pins its text; the two measured divergence triggers are not yet exercised |
| Overlapping comma-separated locations on the same contig duplicate selected records when parallel preprocessing is active | **Indefinite**; preserve record multiplicity | Legacy pre.py runs each requested location as an independent blocksplit stream. Deduplicating would silently alter outputs. `LocationStreamPolicy` owns identity and per-contig geometry for independent streams versus normative set union before generic block expansion. | Pinned PREPY case `matrix_overlapping_locations`, `legacy_only_parallel_overlapping_locations_duplicate_same_contig_records`, `normative_set_union_collapses_overlapping_locations`, and `normative_set_union_multiblock_keeps_records_selected_only_by_later_range` |
| SCMP allele matching constructs `RefVar.end` from ALT length rather than REF length | Legacy-only and indefinite within SCMP allele mode | Pinned legacy SCMP's VCF-to-`RefVar` adapter uses ALT length, keeping some representation-equivalent homopolymer insertions distinct. `ScmpRefVarSpanPolicy` isolates this malformed-span rule from normalization and matching. | `legacy_only_allele_mode_preserves_alt_length_refvar_bug`, `legacy_only_scmp_refvar_span_uses_alt_length`, and `normative_scmp_refvar_span_uses_reference_length` |
| SCMP's bcftools-style duplicate-record merge reuses stale allele-count slots | Legacy-only and indefinite within the merge adapter | bcftools 1.17 `merge --force-samples` rebuilds allele keys but does not clear its reused count array, changing duplicate-indel pairing order. `AlleleCountArrayPolicy` governs storage lifecycle before generic selection. | `legacy_only_duplicate_indels_preserve_bcftools_reused_count_order`, `legacy_only_bcftools_count_table_reuses_stale_slots`, and `normative_allele_count_table_clears_all_slots` |
| Legacy VCF header and classified-record order, genotype canonicalization, primitive uniqueness before padding, duplicate-ALT classified-row projection, symbolic-allele materialization, and numeric rendering | Legacy-only, VCF codec/preprocessing adapters | Pinned hap.py artifact parity requires byte- or value-equivalent serialization. Graph rows put filtered truth matches first, then retain upstream side/type ranks before allele and same-key ranks. VariantPrimitiveSplitter orientation and VariantAlleleUniq's internal-edit identity are preserved before the final VCF writer; a SNP and an indel sharing a position never merge and the location aggregator emits the SNP first whatever order the input listed them in; xcmp's reader then deduplicates equal padded ALT spellings for classified output. A haplotype-graph variant set holds one line per file, so a spelling both sides carry pairs `min(truth copies, query copies)` times and any excess copy orders with the single-file class. | `legacy_graph_order_puts_records_paired_across_inputs_first`, `duplicate_deletion_spelling_beside_snp_keeps_local_mismatch_block_kind` (fast guard for HAPPY case `duplicate_alt_deletion_aggregate`), `production_spool_preserves_truth_before_query_row_rank`, `legacy_only_realigned_mixed_allele_uses_reversed_unphased_het`, `legacy_only_distinct_internal_edits_keep_duplicate_padded_alleles`, `location_aggregator_emits_snp_before_indel_at_a_shared_position` (PREPY case `same_position_snp_order`), `legacy_only_duplicate_alt_query_projects_for_paired_classified_rows`, sort-spool identity regressions, and all PREPY/HAPPY parity rows |
| A shared multi-allelic insertion plus extra query SNPs is classified as `hapfail`, leaving unmatched SNP rows without `BK=lm` | Legacy-only, xcmp classified-output adapter | Pinned HAPPY output for the reduced chr21 insertion-conflict shape reports `hapfail:mismatch` even though both graph signatures are evaluable. The trigger remains limited to a same-anchor Insert+Subst conflict with a matching truth insertion and a multi-allelic query insertion, and only while that conflict is still live — a conflict position must retain an unmatched query record after exact-match draining. When both members of the conflict were exact-matched to TP, the genuine disagreement lies at another anchor, the block is `hap:mismatch`, and the leftover rows keep `BK=lm`. The hapfail suppression also requires the multi-allelic query insertion to be a genuine het-alt (two distinct alts); a duplicate-alt aggregate (`ACCCT,ACCCT`) is a persisted homozygous insertion, not a reciprocal het pair, and stays `hap:mismatch` with `BK=lm`. | `legacy_only_shared_insertion_with_extra_query_snps_keeps_missing_block_kind`, `normative_single_alt_shared_insertion_keeps_local_mismatch_block_kind`, `matched_insert_subst_conflict_does_not_suppress_adjacent_compound_het_mismatch`, and `homalt_insertion_vs_compound_het_query_does_not_drain_as_tp` |
| A query duplicate-alt DELETION aggregate opposite a truth homalt splits into an `am` pair plus an `lm` copy, while an insertion aggregate stays one collapsed homalt row | Legacy-only, xcmp classified-output adapter | When a duplicate-alt aggregate (`X,X` selecting both indices) faces a truth homalt for the same edit in a block that declined, legacy renders the deletion case as two het rows: one deletion copy allele-matches the truth (combined `FN:am`/`FP:am`) and the other is an excess `FP:lm`. Insertion aggregates (`ACCCT,ACCCT` at chr1:150042104, chr3:30297039, chr7:68324852, chr13:74256688) stay collapsed as one `FP:lm` homalt row, so the split fires for deletions only. hap-rs decomposes the deletion aggregate into `[X 0/1, X 1/0]` before the same-locus pairing loop so the `0/1` copy am-matches and the `1/0` copy falls through as `lm`, reproducing the pinned oracle at chr6:91567856 (`TC>T,T 1/2` beside a `T>A 1/0` SNP). | `duplicate_alt_deletion_aggregate_splits_into_am_pair_and_lm_copy`, `homalt_insertion_vs_compound_het_query_does_not_drain_as_tp`, and HAPPY case `duplicate_alt_deletion_aggregate` |
| An ordinary insertion and a mixed-edit deletion sharing an anchor merge into one het-alt record unless a variant sits at the deletion's first deleted base (anchor + 1) | Indefinite, preprocessing location aggregator | Legacy's `VariantLocationAggregator` folds the pair to `ALT=<ins>,<del>` with the opposite-slot genotype (`T>TTTC 0/1` + `TT>C 1/0` → `TT>TTTCT,T 2/1`), matching the pinned hap.py 0.3.15 aggregator at chr1:152194725. A record at anchor + 1 blocks the fold and the pair stays split, as at chr11:95814076 where a variant at 95814077 sits in the deletion span; a record at any later deleted base (95814078) does not block. The lookahead is stream-local: only a record from the same block/location stream at anchor + 1 suppresses the merge. | `location_aggregator_merges_isolated_ordinary_insertion_and_mixed_deletion`, `location_aggregator_keeps_pair_split_when_successor_present`, and PREPY cases `multiallelic_aggregate` (merge) and `multiallelic_aggregate_split` (anchor+1 keeps it split) |
| A multi-allelic record is split into one record per called ALT before the location aggregator, under both `--decompose` and `--no-decompose` | Indefinite, preprocessing allele splitter; single-sample only (multi-sample gate is a known ceiling) | Legacy `VariantAlleleSplitter` fans every multi-allelic record out per allele regardless of `--decompose`; the flag suppresses only the per-allele primitive decomposition that follows. The split lets a colocated single merge into a het-alt (`T>A,G` beside `T>C` → `T>A,C 2/1` + `T>G`) and lets a standalone `TGA>T,TG` fan out to its two per-position deletions. hap-rs runs the split only for single-sample records: the aggregator's het-alt re-merge is not lossless for multi-sample rows (a 2-sample somatic `C>T,G` drops an allele), so those stay whole until a multi-sample-correct re-merge lands. Under `--decompose` only the pure-SNP multi-allelics the primitive splitter leaves whole are pre-split; indel/complex alleles still fan out through the splitter. | PREPY cases `no_decompose_multiallelic_split` and `decompose_multiallelic_split` |
| FTX drops a record when the final `<NON_REF>` allele is selected by a sample genotype | Legacy-only, FTX preprocessing adapter | Legacy ftx.py applies the shared pre.py `calls_non_ref_allele` filter before feature extraction. Uncalled final `<NON_REF>` alleles and `<NON_REF>` alleles that are not final remain eligible. | `legacy_only_called_final_non_ref_alleles_are_removed_by_ftx` and `normative_uncalled_or_nonfinal_non_ref_alleles_are_retained_by_ftx` |
| Legacy report schemas, table names, ROC row multiplicity, floating-point rendering, and runtime-field exclusions | Indefinite, report adapters | Downstream consumers depend on report structure; pinned SciPy/C++ behavior establishes numeric parity. Runtime-only provenance is excluded explicitly. | Report/ROC/metrics tests and the six-lane verification matrix |
| Engine matching, haplotype enumeration, normalization, and metric calculations | Normative domain behavior unless a row above names an emulation | These are scientific algorithms, not a compatibility switch. Differences without documented legacy evidence are bugs. | Engine unit tests plus verification fixtures |
| `--bam` feature extraction in `ftx` and `somatic` | **No legacy reference** | The original never pinned pandas: `happy.requirements.txt` and the bioconda package both take an unconstrained `pandas`. Below 0.20, `ftx.py --bam` and `som.py --bam` raise `KeyError: 'CHROM'` on an index-level `groupby`; from 0.23, som.py raises `OptionError` on `display.height`. Only 0.20.x through 0.22.x runs both, and nothing upstream requires it, so any observation reflects this project's pandas choice rather than legacy behavior. | `ftx_bam_depth`, `ftx_multi_bam`, `somatic_bam_depth`, and `matrix_multi_bam` leave the truth set; hap-rs `--bam` behavior carries `normative_` expectations instead |
| Unsupported or malformed inputs with no pinned legacy evidence | Accidental/fixable | Guessing an emulation would create a new undocumented contract. | Add positive and negative regression evidence before a fix |

The pinned reference image, executable fixtures, and comparison results in
`verification/` are the behavioral evidence, and legacy must run there
unmodified.

That identity is the digest URI plus the conda lock extracted from the image,
observed on linux/amd64, and the comparator carries the same pair because it
decides pass and fail. Running unmodified means no interpreter shim, no source
patch, and no library option override; a setting that configures the interpreter
or VM legacy runs on, or the order the harness schedules it in, is a harness
constraint. Moving any pin re-baselines: legacy is re-run and its fresh output becomes the
expectation, because the harness stores no legacy baselines.

## Deprecation transition

There are no scheduled deprecations. The claimed version narrows nothing on its
own surface.

Future versions may introduce new breaking changes from the legacy code once parity is
established, proven and used.

## Adding or removing emulation

New emulation is accepted only when a pinned parity case or captured upstream
regression demonstrates the behavior. The implementation must include a nearby
comment explaining the source and reason, an inventory row here, and tests
whose names distinguish `legacy_only_` behavior from `normative_` behavior.

Removing an indefinite or legacy-only quirk requires all of the following:

1. Classify it as deprecated in this inventory and emit a user-visible warning
   for at least one minor release when the behavior has a detectable trigger.
2. Document the replacement and the first removal release.
3. Remove it only in a semver-compatible release; observable result changes are
   reserved for a major release after 1.0. Before 1.0, announce them in the
   release notes and version the verification baseline explicitly.
4. Add release notes under a **Compatibility** or **Breaking changes** heading,
   including affected commands, old and new behavior, and migration steps.
5. Replace legacy-only regression expectations with normative expectations;
   never delete the evidence without recording why the policy changed.

Accidental/fixable behavior does not require a deprecation window, but its
release note must still identify any user-visible correction.
