# CLI compatibility ledger

This ledger deliberately separates three questions that were previously
collapsed into one "implemented" label:

1. **Spelling** — does clap accept the legacy token and value shape?
2. **Behaviour** — does the Rust runtime use the parsed value?
3. **Oracle coverage** — has that behaviour been compared with the pinned
   legacy implementation, rather than merely parser-tested?

The legacy source of truth is commit `8401169`:

```text
git show 8401169:src/python/hap.py
git show 8401169:src/python/som.py
git show 8401169:src/python/pre.py
git show 8401169:src/python/ftx.py
git show 8401169:src/python/qfy.py
git show 8401169:src/c++/main/vcfcheck.cpp
```

The Rust spelling source is `rust/src/cli.rs`. Behaviour claims below require
a read in the corresponding backing module (`compare.rs`, `somatic.rs`,
`preprocess.rs`, `ftx/mod.rs`, `quantify.rs`, or `validate.rs`); a field that is
only serialized into run information is **inert**, not wired.

Coverage labels:

- **nf-test** — exercised by the legacy-vs-Rust Nextflow evaluation pipeline.
- **live fixture** — reproduced through the pinned legacy container.
- **saved fixture** — exercised against saved expected output by Rust; this is
  not a fresh legacy execution unless also labelled **live fixture**.
- **parser** — spelling/default/exit-code coverage only.
- **none** — no focused compatibility evidence.

## Executable and subcommand names

The legacy tools are separate executables. The Rust binary exposes them as
subcommands while retaining discoverable legacy aliases.

| Legacy entry point | Rust entry point | Spelling | Coverage |
|---|---|---|---|
| `hap.py` | `hap germline` (hidden `compare` alias) | compatible wrapper form | parser + nf-test |
| `som.py` | `hap somatic` | compatible wrapper form | nf-test + live/saved fixture |
| `pre.py` | `hap pre` (hidden `preprocess`, `prepy` aliases) | compatible wrapper form | nf-test + live/saved fixture |
| `ftx.py` | `hap ftx` (hidden `ftxpy` alias) | compatible wrapper form | nf-test |
| `qfy.py` | `hap quantify` (`qfy` visible alias) | compatible wrapper form | parser + saved fixture |
| `vcfcheck` | `hap validate` (`vcfcheck` visible alias) | compatible wrapper form | parser + live/saved fixture |

`hap --version` is supported by clap. Germline also accepts the legacy
`-v/--version` spelling before required-input validation; `pre` and `qfy`
retain the legacy requirement that their normal required arguments parse
first. The wrapper spellings reproduce the pinned oracle's empty-version
lines (`Hap.py `, `pre.py `, and `qfy.py `).

## Germline: `hap.py` → `hap germline`

`hap.py` also imports the `pre.py` and `qfy.py` option groups. Those inherited
options are included here when the Rust germline surface accepts them.

| Legacy spelling | Rust spelling | Spelling | Behaviour | Oracle coverage |
|---|---|---|---|---|
| truth/query positionals | same | exact | wired | nf-test + fixture |
| `-r`, `--reference` | same | exact; Rust requires it | wired | nf-test + fixture |
| `-o`, `--report-prefix` | same | exact | wired | nf-test + fixture |
| `--pass-only` | same | exact | wired | nf-test + fixture |
| `-R`, `--restrict-regions` | same | exact | wired | fixture |
| `-T`, `--target-regions` | same | exact | wired with legacy start-position semantics (`-R` uses record-span overlap) | unit test |
| `-f`, `--false-positives` | same | exact | wired | nf-test + fixture |
| `-l`, `--location` | same | exact | wired | nf-test |
| `--scratch-prefix` | same | exact | wired | unit test |
| `--keep-scratch` | same | exact | wired | unit test |
| `--threads` | same | exact | propagated to preprocessing and external vcfeval; the native comparator remains deterministic | nf-test + unit test |
| `--stratification`, repeated `--stratification-region`, `--stratification-fixchr` | same | exact | requantifies comparison output with the requested region tables | unit test |
| `--engine`, `--engine-vcfeval-path`, `--engine-vcfeval-template`, `--scmp-distance`, `--lose-match-distance` | same | exact | native xcmp/SCMP lanes and preprocessed RTG vcfeval handoff | parser + nf-test for xcmp/SCMP; fake-RTG contract test + saved focused legacy vcfeval summary/VCF oracle |
| `--preprocess-truth` | same | exact | enables the query preprocessing policy for truth; truth defaults to no left-shift/decomposition | unit test |
| inherited `-L/--leftshift`, `--no-leftshift` | same | exact | query defaults on; truth follows only with `--preprocess-truth`; last override wins | parser + unit test |
| inherited `--decompose`, `-D/--no-decompose` | same | exact | query defaults on; truth follows only with `--preprocess-truth`; last override wins | parser + unit test |
| `--convert-gvcf-truth`, `--convert-gvcf-query`, `--usefiltered-truth` | same | exact | wired independently for truth/query preprocessing | unit test |
| `--preprocessing-window-size`, `-w/--window-size`, `--xcmp-enumeration-threshold`, `--xcmp-expand-hapblocks` | same | exact | preprocessing interaction reset, clustering, enumeration and block expansion are configurable | unit test |
| `--adjust-conf-regions`, `--no-adjust-conf-regions` | same | exact | legacy truth-variant confidence padding can be enabled or disabled | nf-test + unit test |
| `--unhappy`, `--no-haplotype-comparison` | same | exact aliases | selects record-level comparison semantics | nf-test + unit test |
| inherited `--filters-only`, `--bcftools-norm`, `--fixchr`, `--no-fixchr`, `--bcf`, `--somatic`, `--set-gt`, `--filter-nonref`, `--convert-gvcf-to-vcf`, `--gender` | same | exact | forwarded into side-specific preprocessing; argparse last-token overrides are preserved | parser + unit test |
| inherited `-V/--write-vcf`, `-X/--write-counts`, `--no-write-counts`, `--no-json` | same | exact | controls the comparison artifact set | unit test |
| inherited `-t/--type`, `--output-vtc`, `--preserve-info`, `--roc`, `--no-roc`, `--roc-regions`, `--roc-filter`, `--roc-delta`, `--ci-alpha` | same | exact | wired through direct reporting or GA4GH requantification | unit test + focused oracle |
| `--logfile`, `--verbose`, `--quiet`, `--force-interactive` | same | exact | logging controls are wired; force-interactive is accepted for standalone compatibility | unit test |

The governed `scmp-distance` row uses local synthetic inputs with
`--set-gt half --no-leftshift -D`. This produces all nine legacy artifacts and
avoids a pinned-upstream `hap.py` failure in its implicit `first` preprocessing
path (`FORMAT/AD` cardinality becomes invalid before SCMP runs). The Rust
implementation retains focused source-golden tests for the default `first`
semantics; the matrix invocation itself passes all nine artifacts against the
pinned container.

The saved vcfeval oracle is deliberately narrow: one exact SNP exercises Rust
preprocessing, RTG command construction, the GA4GH handoff, quantification, and
summary/VCF output. Normal `cargo test` stubs the external RTG process with the
saved GA4GH VCF, so it does not claim to test the RTG executable itself or
whole-chromosome vcfeval parity.

## Somatic: `som.py` → `hap somatic`

| Legacy spelling | Rust spelling | Spelling | Behaviour | Oracle coverage |
|---|---|---|---|---|
| truth/query positionals, `-o/--output`, `-r/--reference` | same | exact; Rust requires reference | wired | nf-test + fixture |
| `-l/--location` | same | exact | wired with literal post-fixchr contig matching; also restricts automatic FP-region size | unit test |
| `-R/--restrict-regions`, `-T/--target-regions` | same | exact | wired with legacy span (`-R`) versus start-position (`-T`) semantics | unit test |
| `-f/--false-positives` | same | exact | wired | nf-test + fixture |
| `-a/--ambiguous` | same | exact, repeatable | wired, including truth-fixchr BED rewriting and fifth-column labels | unit test |
| `--ambi-fp`, `--no-ambi-fp` | same | exact | wired; last token wins like argparse | unit test |
| `--count-unk`, `--no-count-unk` | same | exact | wired; last token wins like argparse | nf-test + fixture + unit test |
| `-e/--explain_ambiguous` | same legacy underscore | exact | writes legacy class/reason count tables | unit test |
| `-P/--include-nonpass` | same | exact | wired | nf-test + fixture |
| `--fp-region-size` | same | exact | numeric override and automatic FP/reference sizing, restricted by `-l` | unit test |
| `--feature-table` | same | exact legacy allowlist | generic, Strelka SNV/indel (HCC + admix), MuTect, VarScan2, and Pisces tables | nf-test + pinned oracle + parser + unit |
| `--happy-stats` | same | exact | writes legacy summary/extended compatibility reports with prerequisite validation | nf-test + unit test |
| `--roc` | same | exact legacy choices | feature-driven legacy somatic ROC tables, including AF subsets | nf-test + unit test |
| `--bin-afs`, `--af-binsize`, `--af-truth`, `--af-query` | same | exact | wired for stats, extended summaries, and AF ROC subsets | nf-test + unit test |
| repeatable `--bam` | same | exact | shared pure-Rust BAM depth scanner; per-contig coverage is averaged and tripled like legacy before caller feature normalization | pinned direct oracle + nf-test + unit + parser |
| `--normalize-truth`, `--normalize-query`, `-N/--normalize-all` | same | exact | conditional left-alignment, REF mismatch exclusion, and exact deduplication | nf-test + unit test |
| `--fixchr-truth`, `--fixchr-query`, `--fix-chr-truth`, `--fix-chr-query`, `--no-fixchr-truth`, `--no-fixchr-query` | same | exact | independent legacy prefix rewriting; truth setting also controls FP/ambiguous BED rewriting | unit test |
| `--count-filtered-fn` (`-FN`) | same | exact | filtered TP/FP/FN-derived columns with legacy prerequisites; argv rewrite preserves the multi-character short option | unit test |
| `--no-order-check` | same | exact | bypasses the caller TP CHROM/POS order safety check | unit test |
| `--ci-level` | same | exact | wired and range-validated | unit test |
| `--scratch-prefix`, `--keep-scratch`, `--continue` | same | exact | explicit/kept scratch retention, default cleanup, and normalized-input reuse | unit test |
| `--logfile`, `--verbose`, `--quiet` | same | exact | logfile routing, informational verbose mode, and summary suppression | unit test |

The somatic surface now routes every accepted legacy control into comparison,
reporting, normalization, or operational behavior. The governed matrix contains
23 byte-parity cases, including a deterministic BAM-depth row.

## Preprocess: `pre.py` → `hap pre`

| Legacy spelling | Rust spelling | Spelling | Behaviour | Oracle coverage |
|---|---|---|---|---|
| input/output positionals, `-r/--reference` | same | exact | explicit path or legacy `HG19` → `HGREF` → `/opt/hap.py-data/hg19.fa` lookup | nf-test + unit |
| `-l/--location`, `--pass-only` | same | exact | wired | nf-test/fixture |
| `-R/--restrict-regions`, `-T/--target-regions` | same | exact | wired through one effective BED lane | fixture covers `-R` |
| `--fixchr`, `--no-fixchr` | same | compatible; Rust additionally accepts an explicit bool after `--fixchr` | legacy auto-detection and last-option-wins precedence | nf-test + unit |
| `--somatic`, `--set-gt` | same | exact choices | wired | fixture |
| `--filter-nonref`, `--convert-gvcf-to-vcf` | same | exact | wired | fixture construction |
| `--decompose`, `--no-decompose`, `-D`, `--leftshift`, `--no-leftshift`, `-L` | same | exact | default-on with legacy last-option-wins overrides; both disabled preserves the bcftools-view record shape | nf-test + oracle + unit |
| `--filters-only` | same | exact | selected FILTER labels are excluded, with `--pass-only` taking precedence | nf-test + oracle + unit |
| `--gender male/female/auto/none` | same | exact | vcfcheck-compatible chrX inference and male chrX/chrY haploid expansion | nf-test + oracle + unit |
| `-w/--window-size` | same | exact | controls legacy-compatible block-boundary resets when parallel splitting is eligible; inputs below the legacy 100-variant minimum remain unsplit | live option matrix + unit |
| `--bcftools-norm` | same | exact | built-in `norm -f REF -c x -D` equivalent: left-align, exclude REF mismatches, remove exact duplicates | nf-test + oracle + unit |
| `--bcf` and `.bcf` input/output | same | exact semantic BCF2 | dependency-free BCF2 reader/writer with BGZF and CSI output; `--bcf` appends `.bcf` like legacy | nf-test + oracle + round-trip unit |
| `--threads` | same | exact | controls whether eligible inputs use legacy-compatible preprocessing block boundaries | live option matrix + unit |
| `--logfile` | same | compatible diagnostic output | creates the requested path and redirects diagnostics; dynamic Python log text is outside the result-prefix byte contract | live acceptance + unit |
| `--verbose`, `--quiet` | same | exact control flow | diagnostic routing and mutually-exclusive verbosity are wired; diagnostic stream bytes are not report artifacts | live acceptance + unit |
| `--force-interactive` | same | exact | accepted no-op; the standalone Rust process has no SGE dispatch path to override | parser + nf-test |

## Feature extraction: `ftx.py` → `hap ftx`

| Legacy spelling | Rust spelling | Spelling | Behaviour | Oracle coverage |
|---|---|---|---|---|
| VCF, VCF.gz, or BCF input | same | exact | shared dependency-free VCF/BCF loader; FTX preprocessing and feature dispatch are format-independent | pinned FTX BCF oracle + shared BCF oracle + unit |
| input positional, `-o/--output`, `-r/--reference` | same | exact | explicit path or legacy `HG19` → `HGREF` → `/opt/hap.py-data/hg19.fa` lookup; matching legacy, even an explicitly missing reference is ignored unless `--normalize` needs it | pinned oracle + nf-test + unit + parser |
| `-l/--location`, `-R/--restrict-regions`, `-T/--target-regions` | same | exact | wired | pinned option oracles + unit |
| `-P/--include-nonpass` | same | exact | wired | live option matrix |
| `--feature-table`, `--feature-label` | same | exact | wired for generic TP/FN schemas, Strelka SNV/indel (HCC + admix), MuTect, VarScan2, and Pisces tables | nf-test + unit |
| `--fix-chr` | same | exact | wired; opt-in legacy prefix rewrite | pinned option oracle + unit |
| repeatable `--bam` | same | exact | pure-Rust BGZF/BAM depth extraction; per-reference mapped count × sampled mean query length, averaged across BAMs and tripled like legacy | pinned single- and multi-BAM byte oracles + nf-test + unit + parser |
| `--normalize` | same | exact | wired; lazy reference loading, conditional left-alignment, REF mismatch exclusion, and exact deduplication | pinned option oracle + unit |

Legacy `ftx.py` defines no logfile, verbosity, or quiet CLI flags; its use of
Python logging is internal, so there is no missing FTX logging surface.

## Quantification: `qfy.py` → `hap quantify` / `hap qfy`

| Legacy spelling | Rust spelling | Spelling | Behaviour | Oracle coverage |
|---|---|---|---|---|
| input positional, `-o/--report-prefix`, `-r/--reference` | same | exact; Rust requires reference where legacy allowed omission | wired | live xcmp + GA4GH matrix |
| `-V/--write-vcf` | same | exact | wired, including indexed VCF publication | live option matrix + unit |
| `-X/--write-counts` | same | exact | wired; default remains true | parser + live option matrix |
| `--no-write-counts` | same | exact; conflicts with positive switch | wired | parser + live option matrix |
| `--no-json` | same | exact | wired | live option matrix |
| `-t/--type` | same | exact; accepts legacy `xcmp`/`ga4gh` values | both xcmp and GA4GH annotation schemas wired | parser + live xcmp/GA4GH matrix + unit |
| `-f/--false-positives` | same | exact | wired; records outside confidence are excluded on truth and counted as query UNK; confidence size is reported | live option matrix + unit |
| `--stratification`, `--stratification-region`, `--stratification-fixchr` | same | exact; direct regions repeatable | wired for TSV-relative and `NAME:BED` inputs, overlap counts, region sizes, and optional chr normalization | live option matrix + unit |
| `--output-vtc`, `--preserve-info` | same | exact | VTC/XCMP publication and XCMP clean-INFO/preserve split wired | direct live oracle + unit |
| `--roc`, `--no-roc` | same | exact | default `QUAL`, custom INFO/FORMAT fields, and disable switch wired | parser + live option matrix + unit |
| `--roc-regions`, `--roc-filter`, `--roc-delta`, `--ci-alpha` | same | exact | wired, including named ROC regions, filter tiers, custom deltas, and legacy bit-exact Jeffreys intervals | parser + live option matrix + unit |
| `--adjust-conf-regions` | same | exact | truth-derived confidence padding wired | unit |
| `--bcf` | same | exact | BCF/CSI report publication wired | direct live oracle + unit |
| `--threads` | same | exact | accepted; the dependency-free Rust quantifier remains deterministic and serial | parser |
| `--logfile` | same | compatible | validates and creates the requested file; Python logging text is not reproduced | unit |
| `--verbose` | same | exact report artifacts | retains the private Python-compatible `.roc.tsv` intermediate | live option matrix + unit |
| `--quiet`, `--force-interactive` | same | compatible | accepted no-ops because the standalone Rust process has no default interactive summary or SGE dispatch path | parser |

## Validation: `vcfcheck` → `hap validate` / `hap vcfcheck`

| Legacy spelling | Rust spelling | Spelling | Behaviour | Oracle coverage |
|---|---|---|---|---|
| input positional (`input-file` in Boost) | positional or `--input-file` | exact; forms are mutually exclusive | wired | parser + fixture |
| `-o/--output-file` | `-o/--output-json`, visible `--output-file` alias | exact legacy spelling accepted | wired | parser + fixture |
| `-l/--location` | same | accepts legacy `chr` and `chr:start-end` | wired | parser only |
| `-f/--apply-filters <bool>` | same | exact, including required boolean value | wired | parser + unit |
| `--limit-records` | same | exact signed value shape | wired; `-1` is unlimited and other negative limits process zero records as legacy does | parser + unit |
| `--message-every` | same | exact signed value shape | wired; progress goes to stderr instead of legacy stdout so JSON/stdout consumers remain clean | parser + unit |
| `-H/--strict-homref <bool>` | same | exact | wired into legacy-style adjacent overlap accounting | parser + unit |
| `-W/--all-warnings <bool>` | same | exact | wired for reference-padding, overlap, symbolic-ALT, and uncertain-length diagnostics | parser + unit |
| `--check-bcf-errors <bool>` | same | exact | `false` accepted; `true` fails explicitly because the dependency-free parser cannot reproduce htslib BCF translation validation | parser + unit |
| Rust-only `-r/--reference`, `-e/--errors-bed`, `-R/--regions`, `-T/--targets` | no legacy equivalent | extension | wired | fixture covers reference/errors |

## Current compatibility priorities

1. Extend vcfcheck BCF translation checks only if htslib-compatible validation
   is a required compatibility boundary; the boolean and record-limit controls
   now have live matrix coverage.
2. Reproduce qfy's Python logging text only if diagnostic log contents become
   a compatibility requirement; verbose report artifacts are byte-compared.
3. Add focused oracle cases for every currently wired flag labelled `none`;
   parser acceptance alone is not byte-for-byte evidence.
