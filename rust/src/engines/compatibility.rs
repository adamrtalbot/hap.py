//! Explicit policies for pinned engine quirks.

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ScmpRefVarSpanPolicy {
    LegacyAltLength,
    ReferenceLength,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum AlleleCountArrayPolicy {
    LegacyReuse,
    ResetPerTable,
}

pub(crate) fn scmp_refvar_end(
    policy: ScmpRefVarSpanPolicy,
    start: i64,
    reference_len: usize,
    alt_len: usize,
) -> anyhow::Result<i64> {
    let span = match policy {
        ScmpRefVarSpanPolicy::LegacyAltLength => alt_len,
        ScmpRefVarSpanPolicy::ReferenceLength => reference_len,
    };
    Ok(start + i64::try_from(span)? - 1)
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
        assert_eq!(
            scmp_refvar_end(ScmpRefVarSpanPolicy::LegacyAltLength, 9, 1, 3).unwrap(),
            11
        );
    }

    #[test]
    fn normative_scmp_refvar_span_uses_reference_length() {
        assert_eq!(
            scmp_refvar_end(ScmpRefVarSpanPolicy::ReferenceLength, 9, 1, 3).unwrap(),
            9
        );
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
