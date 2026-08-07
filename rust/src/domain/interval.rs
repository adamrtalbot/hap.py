//! Coordinate interval independent of BED parsing and filesystem adapters.

/// Zero-based half-open interval on a named reference sequence.
#[derive(Clone, Debug)]
pub(crate) struct Interval {
    pub(crate) chrom: String,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

impl Interval {
    pub(crate) fn matches(&self, chrom: &str, pos: usize) -> bool {
        self.chrom == chrom
            && pos.saturating_sub(1) >= self.start
            && pos.saturating_sub(1) < self.end
    }

    pub(crate) fn overlaps(&self, chrom: &str, start: usize, end: usize) -> bool {
        self.chrom == chrom && start.saturating_sub(1) < self.end && end > self.start
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uses_zero_based_half_open_coordinates() {
        let interval = Interval {
            chrom: "chr1".to_string(),
            start: 9,
            end: 12,
        };
        assert!(interval.matches("chr1", 10));
        assert!(interval.matches("chr1", 12));
        assert!(!interval.matches("chr1", 13));
        assert!(interval.overlaps("chr1", 12, 14));
    }
}
