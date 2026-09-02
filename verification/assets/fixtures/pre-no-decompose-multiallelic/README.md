# Multi-allelic split before location aggregation

Legacy `pre.py` runs `VariantAlleleSplitter` on every multi-allelic record
before the location aggregator, regardless of `--decompose`. `--no-decompose`
suppresses only the primitive decomposition of each allele, not the split. So a
multi-allelic record that shares a position with a colocated single-allele
record still splits per allele and re-aggregates: the two engines must agree on
record content, not just order.

The fixture is a synthetic 118 bp contig with three anchors.

`chr1:41` is a colocated pure-SNP case: `T>A,G 1/2` alongside `T>C 0/1`. Legacy
splits `A,G`, aggregates `A` with `C` into one het-alt `T>A,C 2/1`, and leaves
`T>G 0/1` on its own.

`chr1:68` is a colocated deletion case: `TGA>T,TG 1/2` alongside `T>C 0/1`. The
two deletion alleles normalise to different positions (`TGA>T` at 68, `GA>G` at
69), so they cannot re-form a het-alt; legacy emits three separate `0/1` rows
plus the `T>C` SNP.

`chr1:95` is a standalone deletion case: `TGA>T,TG 1/2` with no neighbour. The
split still runs — legacy fans it out to `TGA>T` at 95 and `GA>G` at 96 — so the
split is not conditional on a colocated partner; it is the allele splitter
firing on every single-sample multi-allelic.

Both `--no-decompose` and `--decompose` reach the same output here, because no
anchor needs primitive decomposition beyond the allele split. Observed in the
reference environment `happy-0.3.15:7c701db9f05e454a`.
