---
title: Compatibility Policy
description: Governed hap.py compatibility behavior, deprecations, and removal policy.
---

## The claim

hap-rs 1.0.0 is a drop-in replacement for the Legacy implementation, pinned as
`community.wave.seqera.io/library/happy-0.3.15@sha256:4dda6b77c0b1bd778300ea33c498f9d6267537c05d1b655bf9385a11b8702d1c`
together with `verification/containers/happy-0.3.15.conda-lock.txt`.

For any invocation the pinned parsers accept, `hap <subcommand>` produces the
same exit status and the same artifacts as the corresponding legacy tool,
agreeing on every data cell, every categorical verdict label, and record identity
and order. Observed on linux/amd64. Legacy has no build for any other platform,
so no comparison exists elsewhere, and hap-rs behavior on other platforms is
defined by hap-rs and tested against its own expectations.

The only permitted departures are the two entries in the sealed exemption
register below. A single covered invocation where a claimed observable disagrees,
and which is not a register entry, falsifies this.

No 0.x release carries the claim. See
`docs/adr/0006-state-the-claim-at-1-0-0-against-the-pinned-pair.md`.

## Not claimed

Provenance fields: version, timestamp, command line, generated descriptions,
runtime VCF headers. Message text and which stream carries it. Malformed
invocations. A reference supplied through `HGREF` or `HG19`. VCF on standard
input. hap-rs's own added options. `--bam` feature extraction. Platforms other
than linux/amd64.

Runtime and resource use is not a compatibility question. Neither is independent
scientific correctness: this is an agreement claim, and hap-rs reproduces legacy
including where legacy is wrong.

Standard error is not output. hap-rs may emit clearer error messages than legacy,
and doing so is not a divergence.

Recorded quirks belong to one of the governed classes below and are kept at an
input, CLI, codec, or report adapter boundary where practical.

## Exemption register

Sealed at two entries. Within the covered invocation surface, these are the only
behaviors that disagree with the pinned legacy implementation:

1. `--engine-vcfeval-path` and `--engine-vcfeval-template` are accepted and
   ignored. The native engine reads the FASTA supplied with `--reference`, so it
   needs neither an external RTG install nor an SDF bundle.
2. `hap validate --help` returns exit 0.

Both are permanent and neither has a removal date. Adding a third entry requires
maintainer sign-off, so the register's length is checkable at release.

Entry 1 is narrow. Supply a working RTG path and an SDF matching `--reference`, and
hap-rs produces byte-identical `summary.csv`, `extended.csv`, `roc.*.csv.gz`, and
VCF body. Only a path with no working RTG behind it, or an SDF that disagrees with
the reference, falls outside the claim.

`pre` and `quantify` returning exit 0 on an unknown option previously held the
second slot. Malformed invocations now sit outside the supported invocation
surface, so that behavior needs no exemption and is fixed rather than emulated.
See `docs/adr/0003-bound-the-invocation-surface-to-the-pinned-parsers.md`.

## Governance classes

| Class | Meaning | Change policy |
|---|---|---|
| **Indefinite** | Observable behavior on which existing results or automation may depend. | Preserve and regression-test. Removal requires the process below. |
| **Legacy-only** | Emulation activated by a historical command, option, input shape, or output contract. | Keep out of domain algorithms; test with a `legacy_only_` name and test the normative path separately. |
| **Deprecated** | Accepted temporarily to support migration. | Emit a warning with the replacement and removal release. |
| **Accidental/fixable** | An implementation difference with no supported compatibility rationale. | Fix with a regression case; do not turn it into a compatibility promise. |
| **No legacy reference** | Legacy accepts the option, but no reliable legacy behavior exists for it, because the original's own unpinned dependency range cannot run it. | hap-rs behavior is normative: define it, test it against `normative_` expectations, and keep the case out of the truth set. Never invent a legacy observation for it. |

“Legacy-only” describes an explicit legacy trigger, not a global mode. There is
no switch that silently changes all algorithms into a second implementation.

## Compatibility inventory

| Behavior | Class and decision | Rationale or source | Regression evidence |
|---|---|---|---|
| Command aliases (`compare`, `preprocess`, `prepy`, `ftxpy`, `qfy`, `vcfcheck`) and historical option spellings | Indefinite, CLI adapter | Existing wrappers invoke the Python entry-point names. | CLI parser tests in `rust/src/cli_compat/cli.rs` |
| argparse-style repeated options and unique long-option abbreviations | Legacy-only, CLI adapter | hap.py uses argparse and last-value-wins parsing. | `repeated_legacy_options_match_each_wrapper_parser`, abbreviation tests in `rust/src/cli_compat/cli.rs` |
| Invalid `germline` arguments exit 1 and render historical help | Documented hap-rs behavior, CLI adapter | Malformed invocations are outside the covered invocation surface, so this is no longer a compatibility promise. It stays because hap.py automation may distinguish its exit 1 from clap's exit 2. | `legacy_only_germline_invalid_option_retains_failure_exit_one` |
| Unknown `pre` and `quantify` options | **Accidental/fixable**; return a non-zero status | pre.py/qfy.py return success for unknown-option parser errors. That is legacy behavior worth not reproducing, and malformed invocations sit outside the covered invocation surface, so changing it narrows nothing. Missing required arguments already exit non-zero. | Regression evidence moves from the legacy-only expectation to a `normative_` non-zero expectation |
| Missing `pre`/`quantify` arguments and invalid commands without a governed override use clap's non-zero usage exit | Normative | New interfaces must not inherit legacy success codes accidentally. | `normative_pre_and_quantify_missing_arguments_use_standard_nonzero_exit`, `normative_validate_invalid_option_uses_standard_nonzero_exit` |
| `--engine-vcfeval-path` and `--engine-vcfeval-template` are accepted and ignored | **Exemption register entry 1**; no removal date | The native engine reads the FASTA supplied with `--reference`, so it needs neither an external RTG install nor an SDF bundle and has nothing to do with either path. Use `--engine vcfeval --reference <FASTA>`; the options stay supported for wrappers that pass them. | `matrix_vcfeval_deprecated_flags` verification contract; the two measured divergence triggers are not yet exercised |
| Overlapping comma-separated locations on the same contig duplicate selected records when parallel preprocessing is active | **Indefinite**; preserve record multiplicity | Legacy pre.py runs each requested location as an independent blocksplit stream. Deduplicating would silently alter outputs. `LocationStreamPolicy` owns identity and per-contig geometry for independent streams versus normative set union before generic block expansion. | Pinned PREPY case `matrix_overlapping_locations`, `legacy_only_parallel_overlapping_locations_duplicate_same_contig_records`, `normative_set_union_collapses_overlapping_locations`, and `normative_set_union_multiblock_keeps_records_selected_only_by_later_range` |
| SCMP allele matching constructs `RefVar.end` from ALT length rather than REF length | Legacy-only and indefinite within SCMP allele mode | Pinned legacy SCMP's VCF-to-`RefVar` adapter uses ALT length, keeping some representation-equivalent homopolymer insertions distinct. `ScmpRefVarSpanPolicy` isolates this malformed-span rule from normalization and matching. | `legacy_only_allele_mode_preserves_alt_length_refvar_bug`, `legacy_only_scmp_refvar_span_uses_alt_length`, and `normative_scmp_refvar_span_uses_reference_length` |
| SCMP's bcftools-style duplicate-record merge reuses stale allele-count slots | Legacy-only and indefinite within the merge adapter | bcftools 1.17 `merge --force-samples` rebuilds allele keys but does not clear its reused count array, changing duplicate-indel pairing order. `AlleleCountArrayPolicy` governs storage lifecycle before generic selection. | `legacy_only_duplicate_indels_preserve_bcftools_reused_count_order`, `legacy_only_bcftools_count_table_reuses_stale_slots`, and `normative_allele_count_table_clears_all_slots` |
| Legacy VCF header and classified-record order, genotype canonicalization, primitive uniqueness before padding, duplicate-ALT classified-row projection, symbolic-allele materialization, and numeric rendering | Legacy-only, VCF codec/preprocessing adapters | Pinned hap.py artifact parity requires byte- or value-equivalent serialization. Graph rows put filtered truth matches first, then retain upstream side/type ranks before allele and same-key ranks. VariantPrimitiveSplitter orientation and VariantAlleleUniq's internal-edit identity are preserved before the final VCF writer; xcmp's reader then deduplicates equal padded ALT spellings for classified output. | `legacy_graph_order_puts_records_paired_across_inputs_first`, `production_spool_preserves_truth_before_query_row_rank`, `legacy_only_realigned_mixed_allele_uses_reversed_unphased_het`, `legacy_only_distinct_internal_edits_keep_duplicate_padded_alleles`, `legacy_only_duplicate_alt_query_projects_for_paired_classified_rows`, sort-spool identity regressions, and all PREPY/HAPPY parity rows |
| A shared multi-allelic insertion plus extra query SNPs is classified as `hapfail`, leaving unmatched SNP rows without `BK=lm` | Legacy-only, xcmp classified-output adapter | Pinned HAPPY output for the reduced chr21 insertion-conflict shape reports `hapfail:mismatch` even though both graph signatures are evaluable. The trigger remains limited to a same-anchor Insert+Subst conflict with a matching truth insertion and a multi-allelic query insertion. | `legacy_only_shared_insertion_with_extra_query_snps_keeps_missing_block_kind` and `normative_single_alt_shared_insertion_keeps_local_mismatch_block_kind` |
| FTX drops a record when the final `<NON_REF>` allele is selected by a sample genotype | Legacy-only, FTX preprocessing adapter | Legacy ftx.py applies the shared pre.py `calls_non_ref_allele` filter before feature extraction. Uncalled final `<NON_REF>` alleles and `<NON_REF>` alleles that are not final remain eligible. | `legacy_only_called_final_non_ref_alleles_are_removed_by_ftx` and `normative_uncalled_or_nonfinal_non_ref_alleles_are_retained_by_ftx` |
| Legacy report schemas, table names, ROC row multiplicity, floating-point rendering, and runtime-field exclusions | Indefinite, report adapters | Downstream consumers depend on report structure; pinned SciPy/C++ behavior establishes numeric parity. Runtime-only provenance is excluded explicitly. | Report/ROC/metrics tests and the six-lane verification matrix |
| Engine matching, haplotype enumeration, normalization, and metric calculations | Normative domain behavior unless a row above names an emulation | These are scientific algorithms, not a compatibility switch. Differences without documented legacy evidence are bugs. | Engine unit tests plus verification fixtures |
| `--bam` feature extraction in `ftx` and `somatic` | **No legacy reference** | The original never pinned pandas: `happy.requirements.txt` and the bioconda package both take an unconstrained `pandas`. Below 0.20, `ftx.py --bam` and `som.py --bam` raise `KeyError: 'CHROM'` on an index-level `groupby`; from 0.23, som.py raises `OptionError` on `display.height`. Only 0.20.x through 0.22.x runs both, and nothing upstream requires it, so any observation reflects this project's pandas choice rather than legacy behavior. | `ftx_bam_depth`, `ftx_multi_bam`, `somatic_bam_depth`, and `matrix_multi_bam` leave the truth set; hap-rs `--bam` behavior carries `normative_` expectations instead |
| Unsupported or malformed inputs with no pinned legacy evidence | Accidental/fixable | Guessing an emulation would create a new undocumented contract. | Add positive and negative regression evidence before a fix |

The pinned reference image, executable fixtures, and comparison results in
`verification/` are the behavioral evidence, and legacy must run there
unmodified. A tool version alone is insufficient: `hap.py 0.3.15` named two
environments with different dependency sets, which is why authority attaches to
an image digest instead. Where a comparison cannot explain intent, an emulation
must name the upstream hap.py/Python/C++ behavior in a nearby comment. See
`docs/adr/0002-drop-in-claim-against-one-unpatched-legacy-image.md`.

That identity is the digest URI plus the conda lock extracted from the image,
observed on linux/amd64, and the comparator carries the same pair because it
decides pass and fail. Running unmodified means no interpreter shim, no source
patch, and no library option override; a setting that configures the interpreter
or VM legacy runs on, or the order the harness schedules it in, is a harness
constraint, allowed only when named in
`docs/adr/0004-pin-the-legacy-baseline-to-one-container-identity.md`. Moving
any pin re-baselines: legacy is re-run and its fresh output becomes the
expectation, because the harness stores no legacy baselines.

The single-environment rule is in force for all six lanes. The reference is
`community.wave.seqera.io/library/happy-0.3.15:2c2b5746d6b0da37`, rebuilt at
pandas 0.20.3, which runs all six tools unpatched and renders `to_csv` floats
exactly as the 0.19.2 build it replaces. The second image and the som.py
interpreter shim are gone. See the corrections in
`docs/adr/0002-drop-in-claim-against-one-unpatched-legacy-image.md`.

## Deprecation transition

There are no scheduled deprecations. The claimed version narrows nothing on its
own surface.

`--engine-vcfeval-path` and `--engine-vcfeval-template` were previously scheduled
for removal in 1.0.0. That schedule is withdrawn and both options stay supported.
Wrappers that pass them keep working; new wrappers should supply
`--engine vcfeval --reference <FASTA>` instead.

`pre` and `quantify` returning exit 0 on an unknown option was previously
scheduled to change in 1.0.0. It is fixed rather than scheduled: an unknown option
returns a non-zero status.

Warnings go to stderr so report artifacts and machine-readable stdout remain
unchanged.

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
