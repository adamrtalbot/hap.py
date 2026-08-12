# Somatic allele-frequency type parity fixture

This synthetic query-only indel reproduces the allele-frequency type-label
behavior observed in the full-size Seqera som.py run from 2026-08-12. Pinned
som.py repeats the same binned counts under `records`, `SNVs`, and `indels`.
The zero-truth case keeps unrelated confidence-interval values exactly aligned.
