---
title: I/O resource limits
description: Streaming behavior, decoder limits, and memory benchmarks.
---

VCF and BCF inputs are decoded record by record and compressed inputs are
decompressed incrementally. Validation retains only the current record.
Preprocessing uses 65,536-record external-sort chunks, comparison retains at
most one 10,000-variant cluster, quantify uses a disk-backed transformed spool,
and somatic comparison uses per-contig disk spools rather than whole inputs.

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

Run `benches/stream_memory.sh` to capture peak RSS for a chromosome-scale
100,000-record input and a whole-genome-scale 2,400,000-record input. The
larger case contains 24 times as many equal-width records. The script fails if
peak RSS grows by more than 64 MiB, retaining the complete `time` reports as
benchmark evidence.
