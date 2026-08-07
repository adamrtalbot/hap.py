//! Filesystem-free genotype parsing and equivalence rules.

pub(super) fn parse_gt_alleles(gt: &str) -> Vec<usize> {
    gt.split(['/', '|'])
        .map(|part| part.parse::<usize>().unwrap_or(0))
        .collect()
}

pub(super) fn equivalent_gt(left: &str, right: &str) -> bool {
    let mut left_alleles = parse_gt_alleles(left);
    let mut right_alleles = parse_gt_alleles(right);
    left_alleles.sort_unstable();
    right_alleles.sort_unstable();
    left_alleles == right_alleles
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genotype_equivalence_is_unphased_and_tolerates_missing_tokens() {
        assert!(equivalent_gt("1|0", "0/1"));
        assert_eq!(parse_gt_alleles("./2"), vec![0, 2]);
    }
}
