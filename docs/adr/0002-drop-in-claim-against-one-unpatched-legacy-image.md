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

Because the digest is opaque, its recipe and
`verification/containers/happy-0.3.15.conda-lock.txt` are governed artifacts:
they are the only readable account of which pandas, scipy, numpy, and libstdc++
produced a given number.
