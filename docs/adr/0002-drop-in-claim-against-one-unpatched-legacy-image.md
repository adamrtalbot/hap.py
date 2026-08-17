# Drop-in claim against one unpatched Legacy image

hap-rs claims drop-in compatibility with the Legacy implementation: every
meaningful observable agrees, not a scientific-result or interface subset.
"Meaningful" is data cells, categorical result labels, and record identity and
order; provenance fields such as version, timestamp, and command line are
excluded. Authority for what Legacy does rests on one immutable container image
digest running all six tools unmodified. Every observation taken from that image
is authoritative, including one that looks wrong.

## Considered options

A scientific-result or artifact-only claim was rejected because the classified
labels and record ordering are the product, and several ordering contracts
already pass. An interface-only claim was rejected as too weak to mean anything
for a benchmarking tool.

Anchoring authority on the tool version was rejected because `hap.py 0.3.15`
named two behaviours: the frozen Wave image pins pandas 0.19.2 while the
biocontainers image used by the som.py and ftx.py lanes shipped pandas 0.24.2,
and the two render CSV floats differently — `1/3` as `0.333333333333` against
`0.3333333333333333`. Mutable tags were rejected for the same reason: with every
observation authoritative there is no mechanism to decline to follow a reference
that has moved.

Patching Legacy to make it run was rejected outright. The som.py lane had
acquired a `sitecustomize.py` shim that swallowed pandas' removed
`display.height` option, which is how the second image entered the gate and how
`full_repr_float` came to exist alongside `python_repr_float` in the report
adapter. An observation taken through a patch cannot be graded lower than an
unpatched one when all observations are authoritative, so the patch has to go
rather than be documented. This is satisfiable: pandas dropped `display.height`
in 0.20, so pinning 0.19.2 lets som.py run untouched, and the frozen Wave image
already does exactly that.

## Consequences

Adopting one image moves the reference for the som.py and ftx.py lanes from
pandas 0.24.2 to 0.19.2, so their legacy CSVs must be re-observed and
`full_repr_float` retired. Their numpy also differs, 1.16.5 against 1.12.1, and
that effect is unmeasured.

Two deliberate divergences already shipped and stay, as a sealed register rather
than a general escape hatch: `--engine-vcfeval-path` and
`--engine-vcfeval-template` are accepted and ignored, and `pre` and `quantify`
return zero on an unknown option until 1.0.0. Nothing joins that register
without maintainer sign-off, and its length is checkable.

The register's second entry was later replaced. Once malformed invocations moved
outside the claim, unknown-option exit codes stopped needing an exemption, and
`hap validate --help` returning 0 took that slot. The register is still two
entries; see `0003-bound-the-invocation-surface-to-the-pinned-parsers.md`.

## Correction: the named digest cannot run all six tools

The claim above that the Wave image serves as the single authority for all six
lanes was wrong, and measurement settled it. pandas 0.19.2 cannot group by an
index level name. `Tools/bamstats.py` returns a frame indexed on `CHROM`, and both
`ftx.py` and `som.py` then call `pandas.concat(bams).groupby("CHROM")` when
`--bam` is supplied, which raises `KeyError: 'CHROM'` at 0.19.2. Grouping by index
level arrived in pandas 0.20.0. Four samplesheet rows depend on it:
`ftx_bam_depth`, `ftx_multi_bam`, `somatic_bam_depth`, `matrix_multi_bam`.

The related claim that pandas removed `display.height` in 0.20 was also wrong.
It survives as a deprecation at 0.20.3 and 0.22.0, and the removal that breaks
som.py landed later, which is why 0.24.2 raises `OptionError`. So the 0.19.2 pin
was never the only way to run som.py unpatched.

There is a window where one image runs everything unpatched. At pandas 0.20.3
with numpy 1.12.1, index-level grouping works, `display.height` still exists, and
both `to_csv` and `to_json` render `1/3`, `0.1+0.2`, and `123456789.123456789`
exactly as 0.19.2 does, `123456789.123456791` in JSON included. som.py's seven
`.ix` accessors, which pandas 0.20 deprecated, print nothing: CPython 2.7 silences
`DeprecationWarning`, measured at zero bytes on stderr. The
`hap.py-0.3.15-py27hcb73b3d_0` package
pins neither pandas nor numpy, and `pandas-0.20.3-np112py27_0` exists on
conda-forge with the same np112 tagging as the lock's current
`pandas-0.19.2-np112py27_1`, so the rebuild moves one line of the lock rather
than re-solving the environment.

The reference is therefore re-pinned to a rebuilt image at pandas 0.20.3, with
numpy, scipy, and libstdcxx unchanged. The rule this ADR establishes is
unaffected: one immutable digest, six tools, no patching. Only the digest changes.

Adoption is gated on a differential run, and the comparison is three-way rather
than two, because the four `--bam` rows cannot run on the Wave image at all:

| set | currently observed on | disposition at 0.20.3 |
|---|---|---|
| 102 rows: happy 34, prepy 47, qfy 9, vcfcheck 12 | Wave 0.19.2 | **byte-identical.** This is the gate. |
| 52 sompy and ftxpy rows without `--bam` | quay 0.24.2 | move on float rendering, which this ADR already accepted |
| 4 `--bam` rows | quay 0.24.2 | **leave the truth set.** See below. |

If any of the 102 moves, the bump is not free and the decision returns for a
second look. The matrix goes from 158 six-lane comparisons to 154.

Patching legacy to call `groupby(level="CHROM")` was rejected under the
no-patching rule.

## The `--bam` paths have no legacy reference

The original never pinned pandas. `happy.requirements.txt` lists a bare `pandas`,
and `hap.py-0.3.15-py27hcb73b3d_0` depends on an unconstrained `pandas` too. So
the original permits installations where its own features cannot run: below 0.20
`ftx.py --bam` and `som.py --bam` raise `KeyError: 'CHROM'`, and from 0.23 som.py
raises `OptionError` on `display.height`. Only 0.20.x through 0.22.x runs both, and
nothing upstream requires that window. The 0.19.2 pin was a permitted choice, so
this is a defect in the original rather than an error in the recipe.

That makes any `--bam` observation a property of the pandas version this project
selects, not of the legacy implementation. There is no legacy behaviour there to
be authoritative about. The four rows therefore leave the truth set:
`ftx_bam_depth`, `ftx_multi_bam`, `somatic_bam_depth`, `matrix_multi_bam`.

`--bam` is still an accepted option, so it does not fall under the covered
invocation surface as
`0003-bound-the-invocation-surface-to-the-pinned-parsers.md` defines it, nor is it
an hap-rs addition. It gets a third classification, **no legacy reference**:
legacy accepts the option, no reliable legacy behaviour exists for it, and hap-rs
behaviour is therefore normative, defined by hap-rs and tested against its own
expectations. The exemption register stays at two entries, because there is no
observed legacy behaviour to deliberately diverge from.

Removing the rows also removes the one place the numpy 1.16.5 to 1.12.1 residual
could not be isolated. It stays visible in the 52.

Because the digest is opaque, its recipe and
`verification/containers/happy-0.3.15.conda-lock.txt` are governed artifacts:
they are the only readable account of which pandas, scipy, numpy, and libstdc++
produced a given number.
