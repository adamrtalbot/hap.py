//! Pure interval union and intersection calculations.

use crate::domain::Interval;

pub(super) fn region_size(intervals: &[Interval]) -> usize {
    intervals
        .iter()
        .map(|interval| interval.end.saturating_sub(interval.start))
        .sum()
}

fn merged_regions(intervals: &[Interval]) -> Vec<Interval> {
    let mut sorted = intervals.to_vec();
    sorted.sort_by(|left, right| {
        left.chrom
            .cmp(&right.chrom)
            .then(left.start.cmp(&right.start))
            .then(left.end.cmp(&right.end))
    });
    let mut merged: Vec<Interval> = Vec::with_capacity(sorted.len());
    for interval in sorted {
        if let Some(last) = merged.last_mut()
            && last.chrom == interval.chrom
            && interval.start <= last.end
        {
            last.end = last.end.max(interval.end);
        } else {
            merged.push(interval);
        }
    }
    merged
}

pub(super) fn region_intersection_size(left: &[Interval], right: &[Interval]) -> usize {
    let left = merged_regions(left);
    let right = merged_regions(right);
    let (mut left_index, mut right_index, mut size) = (0, 0, 0usize);
    while left_index < left.len() && right_index < right.len() {
        let left_interval = &left[left_index];
        let right_interval = &right[right_index];
        match left_interval.chrom.cmp(&right_interval.chrom) {
            std::cmp::Ordering::Less => {
                left_index += 1;
                continue;
            }
            std::cmp::Ordering::Greater => {
                right_index += 1;
                continue;
            }
            std::cmp::Ordering::Equal => {}
        }
        size += left_interval
            .end
            .min(right_interval.end)
            .saturating_sub(left_interval.start.max(right_interval.start));
        if left_interval.end <= right_interval.end {
            left_index += 1;
        } else {
            right_index += 1;
        }
    }
    size
}

#[cfg(test)]
mod tests {
    use super::*;

    fn interval(start: usize, end: usize) -> Interval {
        Interval {
            chrom: "chr1".to_string(),
            start,
            end,
        }
    }

    #[test]
    fn unions_overlaps_before_intersection() {
        let left = [interval(0, 5), interval(4, 10)];
        let right = [interval(3, 7)];
        assert_eq!(region_size(&left), 11);
        assert_eq!(region_intersection_size(&left, &right), 4);
    }
}
