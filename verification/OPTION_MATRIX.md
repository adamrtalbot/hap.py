# Option coverage contract

This document defines the finite compatibility matrix for the six governed
legacy wrappers. It does not claim that every member of the command-line power
set is meaningful or executable.

## Scope

| Legacy wrapper | Rust command | Matrix lane |
|---|---|---|
| `hap.py` | `hap germline` | `happy` |
| `som.py` | `hap somatic` | `sompy` |
| `pre.py` | `hap pre` | `prepy` |
| `ftx.py` | `hap ftx` | `ftxpy` |
| `qfy.py` | `hap quantify` | `qfy` |
| `vcfcheck` | `hap validate` | `vcfcheck` |

The historical installation also exposed `bamstats.py`, `cnx.py`, `ovc.py`,
and the C++ utilities `multimerge`, `hapenum`, `dipenum`, `blocksplit`,
`hapcmp`, `xcmp`, `scmp`, `alleles`, `quantify`, `vcfhdr2json`, `roc`,
`validatevcf`, `preprocess`, `fastainfo`, and `gvcf2bed`. Those standalone
interfaces are outside this wrapper-compatibility contract. In particular,
the Rust `preprocess` and `quantify` names refer to the `pre.py` and `qfy.py`
wrappers, not the same-named C++ command lines.

Legacy option spellings are catalogued from commit `8401169`; accepted Rust
spellings are catalogued from `rust/src/cli.rs`.

## Finite methodology

The live matrix is a constrained-pairwise sample of behaviorally distinct
classes, not an exhaustive claim over every class listed below. Rows prioritize
successful wrapper paths evaluated under the governed comparator contract;
stable report bytes are compared exactly:

1. **Individual classes:** representative false/true switches, enum values,
   numeric default/interior/boundary values, and valid path classes.
2. **Constrained pairwise interactions:** valid pairs that share a processing
   stage or prerequisite. Invalid or meaningless pairs are excluded explicitly.
3. **Ordering contracts:** both argument orders for last-token-wins controls;
   zero/one/two/duplicate/order classes for repeatable controls.
4. **Format matrix:** representative plain VCF, indexed compressed VCF, and
   BCF/CSI paths, including mixed-format multi-input rows and BCF publication.
5. **Failure contracts:** parser and unit tests cover required input, invalid
   values, missing files/indexes/references, malformed data, unmet
   prerequisites, overwrite protection, and conflicts. These are not currently
   live nf-test rows because the output-tree comparator governs successful
   artifact publication.

The release matrix uses live legacy-to-Rust output comparison. Parser, unit,
saved-fixture, and direct-oracle tests remain useful evidence, but they do not
replace a live matrix coverage ID. `assets/option-coverage.csv` maps each live
coverage ID to one samplesheet row. `factor_family` and `interaction_class` are
human-maintained review labels, not enforced enum values. The nf-test contract
rejects unknown lanes, unknown sample IDs, duplicate coverage IDs, and a
default release run that omits a governed lane.

## Factor catalogue

This catalogue is the compatibility inventory used to select and extend rows.
An entry here is not, by itself, a claim of live coverage; the authoritative
live set is the 111 IDs in `assets/option-coverage.csv`.

### Germline (`happy`)

| Family | Options and classes |
|---|---|
| Invocation | truth/query; `-r`; `-o`; `-v`; hidden `compare` alias |
| Input and output formats | VCF/VCF.gz/BCF per truth/query role; `--bcf`; `-V`; `-X`; `--no-write-counts`; `--no-json` |
| Selection | `--pass-only`; `-R`; `-T`; `-f`; `-l`; `--filters-only`; `--usefiltered-truth` |
| Preprocessing | `--preprocess-truth`; gVCF truth/query conversion; `-L`/`--no-leftshift`; decomposition pair; `--bcftools-norm`; `--filter-nonref` |
| Genotype and contig | `--fixchr`/`--no-fixchr`; `--somatic`/`--set-gt`; every gender value |
| Comparison engine | xcmp; vcfeval path/template; scmp-somatic; scmp-distance and distance alias; unhappy mode |
| Tuning | preprocessing window; comparison window; enumeration threshold; haplotype-block expansion; threads |
| Confidence and stratification | adjust/no-adjust; stratification TSV; repeated direct region; stratification fixchr |
| ROC and reporting | annotation type; output VTC; preserve INFO; ROC field/disable/regions/filter/delta; CI alpha |
| Operations and failures | scratch/keep; logfile; verbose/quiet conflict; force-interactive; missing/malformed/prerequisite/overwrite failures |

### Somatic (`sompy`)

| Family | Options and classes |
|---|---|
| Invocation and formats | truth/query; output/reference; VCF/VCF.gz/BCF per role |
| Selection | location; restrict/target regions; FP BED; include non-PASS; automatic and explicit FP size |
| Ambiguity and unknowns | repeated ambiguous BED; explain; ambi/no-ambi and count/no-count precedence |
| Feature extraction | every accepted generic/Strelka/Mutect/VarScan2/Pisces SNV/indel table; zero/one/two BAMs |
| Normalization and contigs | truth/query/all normalization; independent truth/query fixchr spellings and negative precedence |
| Reporting | happy stats; every ROC choice; AF enable/bin size/truth/query fields; filtered-FN prerequisites; CI level |
| Operations and failures | order check; scratch/keep/continue; logfile; verbose/quiet; missing reference/index and unmet prerequisite failures |

### Preprocess (`prepy`)

| Family | Options and classes |
|---|---|
| Invocation and formats | input/output; reference; VCF/VCF.gz/BCF input crossed with VCF/VCF.gz/BCF output; version; hidden `preprocess`/`prepy` aliases |
| Selection | location; PASS-only; filters-only; restrict/target span semantics |
| Normalization | leftshift pair; decomposition pair; bcftools norm; filter non-reference; gVCF conversion; symbolic alleles |
| Genotype and contig | fixchr precedence and auto mode; somatic/set-gt precedence and every mode; every gender value |
| Tuning and operations | window equivalence classes; threads; logfile; verbose/quiet; force-interactive |
| Failures | missing reference/index; malformed VCF/BCF/BED/FASTA; output overwrite; invalid values and conflicts |

### Feature extraction (`ftxpy`)

| Family | Options and classes |
|---|---|
| Invocation and formats | input/output/reference; VCF/VCF.gz/BCF input; empty input; hidden `ftxpy` alias |
| Selection | location; restrict/target semantics; include non-PASS |
| Feature tables | generic; every Strelka/Mutect/VarScan2/Pisces SNV/indel table; empty and nonempty labels |
| Transformation | normalize; fix-chr; selector-plus-normalize interaction |
| BAM | zero/one/two BAMs; sample-order detection; missing data/index failures |
| Failures | unsupported table; missing normalization reference; malformed input; output errors |

### Quantification (`qfy`)

| Family | Options and classes |
|---|---|
| Invocation and formats | input/report/reference; xcmp and GA4GH annotations; VCF/VCF.gz/BCF input; visible `qfy` alias |
| Confidence and stratification | FP BED; TSV regions; repeated direct regions; stratification fixchr |
| Artifact set | write VCF; write/no-write counts conflict; no JSON |
| ROC | custom field; disable; repeated regions; filter; delta boundaries; CI alpha boundaries |
| Failures | invalid controls; duplicate/reserved regions; missing inputs; input/output overwrite |

### Validation (`vcfcheck`)

| Family | Options and classes |
|---|---|
| Invocation and formats | positional/`--input-file` conflict; VCF/VCF.gz/BCF input; JSON/stdout output; visible `vcfcheck` alias |
| Selection | location; regions; targets; apply-filters true/false |
| Limits and diagnostics | record limit classes; message interval classes; strict hom-ref; all warnings; errors BED |
| BCF validation | check-BCF-errors false and explicit unsupported-true failure |
| Failures | missing/malformed input; invalid booleans/numerics; ambiguous input forms; output errors |

## Maintenance rules

- Give every live behavior class or interaction a stable coverage ID.
- A coverage ID may move to a better row, but must not disappear silently.
- Samplesheet row deletion requires removing or remapping its manifest entries.
- Focused local runs may select fewer lanes. An unoverridden release run must
  select exactly all six governed lanes.
- New accepted options must be added to this catalogue before release.
