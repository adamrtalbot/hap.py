# SNP before indel at a shared position

A SNP and an indel that share a position never merge into one record, and
`pre.py --decompose --leftshift` emits the SNP first whatever order the input
listed them in. Observed in the reference environment: the ordering holds for a
deletion and an insertion alike, for matching and opposite haplotypes, for
`1/1` calls, and when the two rows carry different FILTER values. It is not an
allele sort — the SNP leads even where its ALT sorts after the indel's
(`A>T` before `A>ACC`).

The fixture is a synthetic 118 bp contig with no repeat around any anchor, so
both engines resolve every deletion and insertion to the input position and
only the row order is under test. Each of the four positions lists the indel
first: `chr1:20` a `LowQual` deletion against a `PASS` SNP, `chr1:41` a
deletion on the same haplotype, `chr1:68` an insertion on the opposite
haplotype, and `chr1:100` a hom-alt pair.

`--no-leftshift --no-decompose` skips normalisation entirely and both engines
then preserve the input order, so the ordering belongs to the location
aggregator rather than the writer.
