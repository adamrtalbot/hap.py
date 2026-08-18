# Compare artifact bytes except where the encoding carries provenance

Byte identity is the default for every member of the observed set. A documented
content comparison replaces it only where byte identity is unachievable because
the encoding carries a provenance field or is not a function of the result. The
comparator implements this contract, and a change to the comparator is governed
in the direction it moves the contract.

## Why byte identity is the default

`0002-drop-in-claim-against-one-unpatched-legacy-image.md` makes data cells,
verdict labels, and record identity and order the meaningful observables, and
excludes provenance. Neither statement licenses a comparison weaker than bytes
for anything else. Where the bytes of an artifact are a function of the result
alone, comparing them is the cheapest faithful reading of the claim, and any
normalization applied before comparison is a hole nobody can see into: it
absorbs a whole class of difference and no case can ever surface it.

## Where byte identity is unachievable, and why

Measured across 131 paired cases from a published harness run, covering the
`ftxpy`, `happy`, `prepy`, `sompy` and `vcfcheck` lanes.

| artifact class | raw-byte matches | cause |
|---|---|---|
| `.gz` | **0 of 179** | legacy's report `.gz` are plain gzip whose header carries the wall-clock second the file was written: `mtime=1786990993`, `XFL=2`, against hap-rs's `mtime=0`, `XFL=0` |
| `.bcf` | **0 of 4** | BGZF-framed binary; deflate output is not a function of the record set |
| `.tbi`, `.csi` | **0 of 68** | an index stores virtual offsets into the BGZF stream, so its bytes are downstream of the compressed container's bytes |

The gzip mtime settles the compression question by itself. It is a timestamp,
and `0006-state-the-claim-at-1-0-0-against-the-pinned-pair.md` already excludes
timestamps, so requiring hap-rs to reproduce it would require hap-rs to forge
legacy's clock. Decompressing before comparison is therefore the faithful
implementation of that exclusion rather than a relaxation of it.

The index case follows from the compression case rather than standing on its
own. Requiring index bytes to match would require the container's bytes to
match, which would require reproducing zlib's output exactly. Compression
metadata, Tabix and CSI, and BCF are one forced decision, not three open ones.

`.vcf.gz` and `.vcf.gz.tbi` are BGZF on both sides, `mtime=0` with the EOF block
present; their residual difference is deflate parameters and blocking only.

## The contract, one artifact class at a time

### Uncompressed text, and non-ROC CSV

Raw bytes, after removing the named provenance columns at byte level. No line
splitting, no CSV re-serialization, no row padding.

This tightens the comparator. It previously compared through `splitlines()` and
a `csv.writer` round trip, which masked seven classes of difference. Each was
probed directly against the comparator; all seven passed silently while a
control cell change was caught:

CRLF against LF, a missing trailing newline, U+000B/U+000C/U+0085/U+2028 read as
line breaks, CSV quoting style, a short row against a padded row, a CSV trailing
newline, and a VCF body line terminated with CRLF.

None of the seven is active. Measured across 374 paired text artifacts,
uncompressed and decompressed alike: zero CRLF on either side, zero exotic line
breaks, zero cases where the trailing newline differs, zero row-width
differences, and quote characters present on both sides in the 17 `ftxpy`
`.csv` where they appear at all. Seventy-one artifacts lack a final newline on
both sides, so hap-rs already reproduces legacy's missing terminator.

None of the seven was chosen, either. They arrived together in
`0bee33c test(verification): consolidate parity gate`, whose stated purpose was
moving the comparison out of shell scripts into Python; `splitlines()` and
`csv.writer` are the idiomatic Python calls and their masking is a side effect
of that translation. The compression and BCF relaxations, by contrast, have a
stated reason in `fbe03d0 test(verification): compare parity artifacts
semantically`: "without requiring unstable encoded bytes." That reason is
correct and those relaxations stay.

Byte identity here costs nothing to adopt, and this is measured against the
contract itself rather than a proxy. Comparing decompressed content byte for byte
after byte-level column removal matches on **227 of 227** artifacts in the
population this clause governs: 93 ROC `.csv.gz`, 8 ROC plain, 100 non-ROC plain,
and the 26 carrying `sompyversion` and `sompycmd`. Removing those two columns
while preserving every kept field's original bytes is what makes the last group
match, and it is also the direct evidence that the `csv.writer` round trip
absorbed nothing: the original quoting is already identical on both sides.

### ROC tables

Ordered text, like every other CSV. There is no multiset comparison and no
per-case strength setting.

`0002-drop-in-claim-against-one-unpatched-legacy-image.md` names record order as
a meaningful observable, so a multiset comparison declined to check something the
claim covers. Ordered comparison is a pure tightening: it implies multiset
equality.

The relaxation is also spent. Comparing every ROC CSV and every ROC table in
JSON as ordered passes all 131 cases. It was built for a real defect —
`8666456 fix(verification): compare ROC rows strictly` reads "ignore row order
without accepting malformed, missing, extra, or duplicate ROC rows" — and that
defect is fixed. This is a retirement, not the removal of an accident.

Two further canonicalizations go with it, on the same evidence. Sorting ROC
table rows inside `metrics.json` and renumbering a table's generated index
column both absorb nothing: removing the table canonicalization entirely passes
131 of 131.

Retiring the relaxation also removes an inconsistency nobody chose. ROC
classification tested `prefix + '.roc.'`, so `somatic`'s per-allele-fraction
files, `result.SNVs.0.500000-1.000000.roc.csv` among them, were already compared
as ordered text while `result.roc.all.csv.gz` was not.

### Typed JSON

The parsed tree, with provenance pointers excluded, and every compared number
required to be spelled identically.

The pointer list stays governed as comparator code with cases in
`verification/tests/diff.nf.test`, and is deliberately not reproduced here.
`0006-state-the-claim-at-1-0-0-against-the-pinned-pair.md` keeps its exclusions
at class level precisely to avoid a second inventory that can drift, which is the
failure `0003-bound-the-invocation-surface-to-the-pinned-parsers.md` rejected for
the option surface and `0004-pin-the-legacy-baseline-to-one-container-identity.md`
rejected for the dependency list. The rule is that provenance is excluded; the
pointers are how the comparator says so.

Fifteen of the comparator's nineteen named exclusions absorb a real difference.
The largest are the bcftools header pattern at 63 cases, and `/runInfo`,
`/timestamp` and `/metadata/required/description` at 49 each. Four absorb
nothing anywhere: `/version`, `/metadata/required/version`, `##fileDate=` and
the `sompyversion` CSV column. Three of them stay — `/version`,
`/metadata/required/version` and `sompyversion` — because they implement the
provenance classes ADR 0006 states rather than patching an observed difference;
`/version` absorbs nothing only because hap-rs currently reproduces legacy's
empty version string. `##fileDate=` does not stay as a named pattern: the header
rule below subsumes it.

Requiring identical number spelling is new and costs nothing to adopt: across 83
JSON artifacts and 11,689 float tokens in the compared subtrees, zero differ.
It closes a hole this project has already been bitten by. Tree comparison parses
before comparing, so two spellings of the same double are equal, and ADR 0002
spent a correction on exactly that class of difference —
`123456789.123456791` against `123456789.12345679`, `full_repr_float` against
`python_repr_float`.

Object key order stays uncompared, and this is the one place the contract
knowingly declines a difference that is live: key order differs in 23 of 83 JSON
artifacts, all `happy` `runinfo.json`. Key order is a serialization property that
no parser shows a consumer, so it is not record order in the sense ADR 0002
means. Requiring it would convert a passing lane into port work for no
observable gain.

One exclusion is not provenance and is named here so it cannot be mistaken for a
third register entry. `/final_args/engine_vcfeval_template` is the artifact
footprint of the exemption register's first entry: legacy echoes the supplied
SDF path, `vcfeval-template.sdf`, while hap-rs ignores the option and echoes
`None`, in the 5 vcfeval cases. The footprint is complete, because
`engine_vcfeval_path` never appears in `final_args` at all. The register stays at
two entries.

### VCF records against runtime headers

A header line is excluded when the observed invocation wrote it, unless the line
belongs to the named structural set, which is always compared: `##fileformat`,
`##INFO`, `##FORMAT`, `##FILTER`, `##ALT`, and `##contig`. Every other header
line, and every record, is compared as ordered text.

`##fileformat` is in that set because measurement put it there. It is
invocation-written in 48 cases and agrees in all 68, so excluding it would stop
checking that hap-rs writes a well-formed VCF at all — and unlike the rest of the
set it is not a structured declaration, which is why the set is named rather than
described.

`##reference=` is the boundary case and it is classed **provenance**, so it is
excluded. Its value came off the command line, and the rule takes no carve-out
beyond the named set. The cost is accepted and stated below: 48 agreeing lines
stop being compared, and a regression writing the wrong reference name would
pass.

This replaces four literal patterns with a rule, so a provenance header nobody
anticipated is excluded without a comparator change. The four patterns were
narrow but not wrong: three earn their place, `##bcftools_*Command/Version` at
63 cases and `##source=` and `##CL=` at 5 each, while `##fileDate=` absorbs
nothing.

The structural set is what makes the rule safe, and it is not a detail. The
germline output carries 27 self-written `##INFO`, `##FORMAT` and `##FILTER`
declarations, including `BD`, `BK`, `BI`, `BVT`, `BLT` and `QQ` — the very fields
whose values ADR 0002 names as the verdict labels. An authorship rule without the
set would stop comparing the schema of the product.

Authorship is computed per case against that case's own input headers, over the
66 of 68 paired VCF and BCF cases whose inputs resolve locally; the 2 remaining
`prepy` cases name `example/happy/NA12878_chr21.vcf.gz`, which is not present.
Measured that way, the rule sorts the header population as follows:

| group | keys | today |
|---|---|---|
| written by the invocation, disagrees | `##bcftools_{view,concat,annotate,norm,merge}{Command,Version}`, `##CL=` | legacy emits, hap-rs emits none, so line counts differ |
| written by the invocation, native-engine substitution | `##source=` in the 5 vcfeval cases | legacy `RTG Tools 3.12.1 / Core 1581c65779`, hap-rs `hap-rs native vcfeval (RTG Tools 3.12.1 compatible)` |
| written by the invocation, agrees | `##reference=` in 48 lines, value `hg19`; `##fileDate=` in 5 | agrees, and stops being compared |
| inherited from the input files | 2 `##bcftools_viewCommand` and 2 `##bcftools_viewVersion` | agrees, and **starts** being compared, because the retired pattern rule discarded them by key |
| named structural set | 27 `##INFO`/`##FORMAT`/`##FILTER` per germline output, plus `##fileformat`, `##ALT`, `##contig` | agrees, and stays compared |

The cost is named rather than hidden: 48 `##reference=` lines and 5
`##fileDate=` lines currently agree and stop being compared. The gain measured on
the resolvable corpus is narrower than the rule's motivation suggests: 4 inherited
bcftools header lines start being compared, where the retired pattern rule
discarded them by key regardless of who wrote them. The inherited provenance the
rule is really for — `##gvcftools_cmdline`, `##annotator`,
`##annotationservice{uri,version}`, `##UnifiedGenotyper`, `##CombineVariants`,
`##ApplyRecalibration`, `##MaxDepth_chr*`, and an inherited `##source=GATK 1.6` —
sits entirely in the 2 cases whose input file is absent, so it is unmeasured
here.

Authorship is computed, not judged: a line the invocation wrote is a line absent
from the input files' headers. The comparator therefore needs those headers,
which is a new input to `DIFF_OUTPUTS`, and a case whose inputs are unavailable
cannot be classified at all.

### BCF

Decoded through the pinned comparator's `bcftools view --no-version -Ov` and
compared as VCF text under the header rule above. Raw bytes are never sufficient
and never a shortcut.

Decoding is legitimate because the instrument is pinned.
`0004-pin-the-legacy-baseline-to-one-container-identity.md` moved the comparator
into a container with its own digest and lock for this reason: before that,
`DIFF_OUTPUTS` ran whatever `bcftools` was on `PATH`, and a bcftools change could
flip a verdict without either implementation moving.

### Tabix and CSI

A basic check only. Each index must be valid against its own data file, read
through the pinned comparator's `bcftools index --stats`. The contents of the two
sides' indexes are not compared.

An index is a derived access structure, not a result, and its bytes are
downstream of a container the contract does not compare. The name is still
observed: `0007-observe-the-output-prefix-set.md` makes name-set equality
binding, so an index that fails to appear on one side falsifies the claim. Only
its content is unclaimed.

The stronger option was measured and declined. Comparing the `index --stats`
output between the two sides — contig, length, and record count — passes 131 of
131, so it was available at no cost. It is not adopted, because it would put a
derived structure inside the contract for a guarantee the compared `.vcf.gz`
content already carries.

### Compression metadata

Unobserved. No claim is made on deflate parameters, compression level, BGZF
block size, the gzip header's mtime, XFL or OS bytes, or the presence of a BGZF
EOF block.

Compression state is nonetheless bound, because it is visible in the artifact
name: `metrics.json.gz` against `metrics.json`, which ADR 0007 already settles
through set equality. Recompressing an artifact changes its name and falsifies
the claim.

### Binary, and any unknown extension

sha256 over the raw bytes, which fails loudly.

This is the designed default rather than a gap. Every extension in the observed
set today is in a known class — `.csv`, `.csv.gz`, `.json`, `.json.gz`,
`.vcf.gz`, `.bcf`, `.tbi`, `.csi` — so the binary branch has no instances. An
artifact class nobody has decided about must therefore fail until somebody
decides, which is the correct forcing function. A future compressed artifact
with a new extension failing on sha256 is the mechanism working.

## Governing a comparator change

The comparator decides pass and fail, so a change to it is a change to the
claim's evidence. Governance is asymmetric in the direction the change moves the
contract.

A change that **narrows** what is compared — a new exclusion, a new
canonicalization, a new relaxation, a new artifact class admitted to a weaker
comparison — requires an amendment to this ADR, maintainer sign-off, and a case
in `verification/tests/diff.nf.test` demonstrating the new behaviour.

A change that **widens** what is compared — a tightening, a new artifact class
admitted to a stricter comparison, or a bug fix — requires only the test case.

Both tiers are checkable, because `verification/tests/diff.nf.test` is the
record and it lives beside the comparator. It carries 14 cases today, covering
text line reporting, VCF provenance against preserved header comparison, runinfo
container identity, ROC order and multiplicity, and BCF and CSI semantics.

The comparator's container identity is pinned exactly as the legacy reference is,
so moving that digest re-baselines under the rule in ADR 0004.

Nothing in the comparator may weaken the contract for one case. `strict_roc_order`
was a samplesheet column defaulting to the looser comparison and set on exactly
1 of 34 `happy` rows, which made weakening the contract for any case a CSV cell
edit with no ADR, no review, and no record. Retiring ordered ROC comparison
removes the column, and no per-case strength setting replaces it.

## Considered options

Requiring raw bytes everywhere was rejected on measurement: 0 of 179 `.gz`, 0 of
4 `.bcf`, and 0 of 68 indexes match, and the gzip mtime makes the first of those
a provenance field ADR 0006 already excludes.

Recording the compression relaxation as an exemption was rejected. The register
is sealed at two entries, and there is nothing to exempt: the claim never reached
compression framing, because framing is not a data cell, a verdict label, or a
record position.

Enumerating the excluded JSON pointers and header patterns in this ADR was
rejected as the second inventory ADRs 0003, 0004 and 0006 each rejected in turn.
The rule belongs here; the pointers belong in the comparator with its tests.

Comparing index contents through `bcftools index --stats` was measured, found
free, and declined, on the grounds above.

Requiring object key order in JSON was rejected. It fails 23 cases today for a
serialization property no parser exposes.

Stating the header rule by value shape — excluding any header whose value is a
path, date, version or command line — was rejected on measured cost. It would
stop comparing `##reference=` in 50 files, `##fileDate=` in 5, the 2 inherited
`##source=`, and the whole inherited provenance tail across 2 files, all of which
agree today. A namespace rule was rejected in the other direction as barely
broader than the four literal patterns it replaced.

Symmetric governance, requiring an ADR for every comparator change, was rejected
for taxing the tightenings and bug fixes that should be cheap. Freezing the
comparator at 1.0.0 was rejected because a comparator bug found after release
would then be unfixable without re-baselining the claim.

## Consequences

The comparator changes, in `verification/modules/diff.nf`, and every change here
is a widening except the header rule, which narrows and is therefore recorded by
this ADR:

- Uncompressed text and decompressed `.gz` text compare as raw bytes. The line
  splitter becomes byte-level so a `/line/N` location can still be reported, and
  the named provenance columns are removed without re-serializing the row.
- The ROC multiset path, the ROC table row sort, and the generated-table-index
  renumbering are deleted, along with the `strict_roc_order` parameter.
- Compared JSON numbers must be spelled identically.
- Header exclusion becomes the authorship rule with the schema exemption, and
  `DIFF_OUTPUTS` gains the input VCF headers as an input.

`strict_roc_order` leaves `verification/assets/samplesheet.happy.csv` and
`verification/assets/samplesheet.happy.public.csv`.

`verification/tests/diff.nf.test` gains cases for the four changes above, and for
three behaviours it does not currently cover: an `artifact_set` difference, the
sha256 fallback on an unknown extension, and `.gz` decompression.

The exclusion lists in `verification/README.md` and
`docs/src/content/docs/project/verification.md` are rewritten against this ADR.
Both currently describe ROC CSVs as row multisets and JSON as canonicalized for
ROC row order and generated table indexes, which this decision retires, and both
describe the VCF exclusion as named runtime fields rather than as a rule.

`CONTEXT.md` gains **equivalence contract** and **unstable encoding**.

The exemption register is unchanged at two entries.

## Residual risks

The evidence corpus has no `qfy` lane. All 131 paired cases come from `ftxpy`,
`happy`, `prepy`, `sompy` and `vcfcheck`, so quantify's `.metrics.json.gz` and
`.roc.all.csv.gz` are unmeasured here, and the tightenings this ADR adopts are
free on five of six commands rather than on six. Which cases the campaign runs
belongs to
[Define the coverage model and case selection rule](https://github.com/adamrtalbot/hap.py/issues/23).

Four things stay outside the contract by decision, each recoverable only by
amending this ADR: JSON object key order, index content, compression framing, and
the `##reference=` header. The first is live today at 23 artifacts, and the last
costs 48 agreeing lines; the middle two are not observable at all without
reproducing zlib's output.

The header rule's benefit is measured on 66 of 68 paired VCF and BCF cases. The
inherited GATK and gvcftools provenance that motivates the rule sits in the 2
cases this environment cannot resolve, so the rule is adopted on its reasoning
and its cost, with its principal benefit unmeasured here.

An artifact whose bytes are a function of the result but whose extension is
unknown will fail on sha256 rather than being compared as the text or JSON it
is. That is intended, and the fix is a decision under the widening tier, not a
comparator patch.

## Correction: the ROC evidence was taken on the previous reference image

This correction covers the ROC clause. The ROC clause above rests on three
counts: ordered comparison of every ROC CSV
and every ROC table in JSON "passes all 131 cases", removing the table
canonicalization entirely "passes 131 of 131", and 11,689 float tokens across 83
JSON artifacts agree. All three were measured on
`community.wave.seqera.io/library/happy-0.3.15:41c2102638513597`. `65476ee` moved
`params.legacy_image` to `...:2c2b5746d6b0da37`, and under
`0004-pin-the-legacy-baseline-to-one-container-identity.md` that re-baselines,
so the counts describe an image the harness no longer runs.

Measured on the current pin at `18548d1`, by running the gate twice on one host
and comparing the two runs against each other: legacy's own ROC output moves
between runs. 14 of its 171 `.csv.gz` and 4 of its 111 JSON artifacts differ in
content between the two runs, all of them HAPPY ROC tables belonging to `chr21`,
`chr21_region`, `chr21_passonly` and `chr21_xcmp_controls`. Legacy retains a
different set of QQ threshold rows each time; rows present in both agree cell for
cell, and one row emitted in the first run and not the second carries a number
where the variant type belongs. hap-rs reproduced every artifact byte for byte
across the same pair. The tables are in `verification/README.md` under "Replay
and resource baseline", and the defect is
[#46](https://github.com/adamrtalbot/hap.py/issues/46).

The decision is unchanged, and this ADR's consequences are not yet implemented in
`verification/modules/diff.nf`. What changes is the claim that tightening ROC
comparison costs nothing: on the current image the multiset comparison this ADR
retires already fails two cases per run, a different two each time, and the
ordered comparison it adopts would fail them harder. The tightening waits on #46
rather than on a fresh count.

The header-authorship counts in the same ADR come from the same image and are
superseded too, without threatening the rule they support. That table reads 68
paired VCF and BCF cases with 66 resolvable, `##reference=` in 48 lines and
`##fileDate=` in 5. The current corpus publishes 78 `.vcf.gz` and 5 `.bcf` on
each side, and legacy reproduced all 83 across the two runs once runtime headers
are dropped, so nothing there needs a decision, only a fresh count when the
comparator change lands.

Whether the earlier image was stable and this one is not, or whether the movement
was always present and the old expectation happened to sit on the stable side, is
unmeasured. The measurement above was taken on an arm64 host with the legacy
container emulated, which `verification/README.md` classes as a development
diagnostic, so a native amd64 run is owed before anything is concluded about the
image itself.
