---
title: I/O resource limits
description: Streaming behavior, decoder limits, and memory benchmarks.
---

VCF and BCF inputs are decoded record by record and compressed inputs are
decompressed incrementally. Validation retains only the current record.
Preprocessing uses 65,536-record external-sort chunks and 32-way merge passes.
Comparison retains at most one 10,000-variant cluster and returns a contextual
resource error instead of splitting a larger connected cluster and changing
its semantics. It externally orders 65,536-row chunks with the same bounded
fan-in. Completed preprocess and comparison chunks are closed paths, so only
the active merge's at most 32 readers are open. Comparison metadata is spooled
per contig and rescanned over the actual emitted span, including left-shifted
rows. Quantify uses a disk-backed transformed spool. Germline and somatic ROC
observations use 16,384-record sorted chunks and bounded 32-way merge passes;
there is no whole-genome observation or threshold-count cap. Somatic comparison
closes each per-contig spool before opening the next and retains records for
only the active contig; feature rows are written to five classification spools
and renumbered with streaming passes. Validation writes `--errors-bed` records
directly to a buffered output.

The parser rejects inputs that exceed these hard limits before allocating the
declared payload:

- VCF logical line: 64 MiB
- BCF header: 64 MiB
- BCF shared plus individual record payload: 64 MiB
- BCF samples per record: 1,000,000
- BCF INFO or FORMAT fields per record: 65,535
- Quantify records per active superlocus: 1,000,000
- Somatic records per active contig and side: 10,000,000
- Comparison variants per active cluster: 10,000

These limits are intentionally far above normal variant records but low enough
to prevent corrupt length and count fields from causing unbounded allocations.
Errors include the path and, for BCF records, the one-based record number.

Run `benches/stream_memory.sh` to capture repeated peak RSS measurements for
validate with `--errors-bed`, preprocess, compare with filtered truth and ROC
disabled, quantify with ROC both disabled and enabled, and somatic generic
feature/happy reports. It compares a chromosome-scale 100,000-record input to
a whole-genome-scale 2,400,000-record input and uses the median of three runs
by default. Complete `time` reports and the median table are preserved under
`target/stream-memory-reports` (or `HAP_MEMORY_REPORT_DIR`) after input cleanup.
The generated whole-genome case distributes records across 24 representative
contigs. The script fails if a workflow or peak-RSS measurement fails, rejects
missing or non-positive RSS values, and proves only that these named scenarios
grow by no more than 64 MiB. Dense connected clusters, superloci, and somatic
contigs are instead protected by the explicit limits above; ROC observation
counts are handled by external sorting rather than rejected.
