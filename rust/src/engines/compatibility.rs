//! Explicit policies for pinned engine quirks.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum AlleleCountArrayPolicy {
    LegacyReuse,
    ResetPerTable,
}

/// Pinned SCMP constructs `RefVar.end` from ALT length, not reference length.
pub(crate) fn scmp_refvar_end(start: i64, alt_len: usize) -> anyhow::Result<i64> {
    Ok(start + i64::try_from(alt_len)? - 1)
}

pub(crate) fn begin_allele_count_table(policy: AlleleCountArrayPolicy, counts: &mut Vec<usize>) {
    if policy == AlleleCountArrayPolicy::ResetPerTable {
        counts.clear();
    }
}

pub(crate) fn record_allele_count(counts: &mut Vec<usize>, slot: usize, first_record: bool) {
    if counts.len() <= slot {
        counts.resize(slot + 1, 0);
    }
    if first_record {
        counts[slot] = 1;
    } else {
        counts[slot] += 1;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_only_scmp_refvar_span_uses_alt_length() {
        assert_eq!(scmp_refvar_end(9, 3).unwrap(), 11);
    }

    #[test]
    fn legacy_only_bcftools_count_table_reuses_stale_slots() {
        let mut counts = vec![2, 7];
        begin_allele_count_table(AlleleCountArrayPolicy::LegacyReuse, &mut counts);
        record_allele_count(&mut counts, 0, true);
        assert_eq!(counts, [1, 7]);
    }

    #[test]
    fn normative_allele_count_table_clears_all_slots() {
        let mut counts = vec![2, 7];
        begin_allele_count_table(AlleleCountArrayPolicy::ResetPerTable, &mut counts);
        record_allele_count(&mut counts, 0, true);
        assert_eq!(counts, [1]);
    }
}
