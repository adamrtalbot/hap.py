# Filtered-truth ROC retention

This fixture keeps three filtered truth calls and mixes matched, truth-only,
and query-only SNP observations across distinct QUAL thresholds. It exercises
the legacy C++ raw ROC table's temporary genotype rows, unordered-map rehashes,
threshold retention, and terminal cumulative row under `--usefiltered-truth`.
The Rust behavior predates this fixture; this lane pins that existing contract
and enables ordered ROC comparison through its samplesheet metadata.
