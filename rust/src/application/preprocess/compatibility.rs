//! Governed location-stream behavior for preprocessing.

use crate::adapters::vcf::LocationFilter;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum LocationStreamPolicy {
    IndependentLegacyStreams,
    #[cfg(test)]
    SetUnion,
}

pub(super) fn location_stream_groups(
    policy: LocationStreamPolicy,
    filters: &[LocationFilter],
    chrom: &str,
    pos: usize,
) -> Vec<usize> {
    let matching = filters
        .iter()
        .enumerate()
        .filter_map(|(index, filter)| filter.matches(chrom, pos).then_some(index))
        .collect::<Vec<_>>();
    match policy {
        LocationStreamPolicy::IndependentLegacyStreams => matching,
        #[cfg(test)]
        LocationStreamPolicy::SetUnion if matching.is_empty() => Vec::new(),
        #[cfg(test)]
        LocationStreamPolicy::SetUnion => vec![0],
    }
}

pub(super) fn location_stream_final_end(
    policy: LocationStreamPolicy,
    filters: &[LocationFilter],
    stream_id: usize,
    chrom: &str,
) -> Option<usize> {
    match policy {
        LocationStreamPolicy::IndependentLegacyStreams => filters
            .get(stream_id)
            .and_then(|filter| location_filter_end(filter, chrom)),
        #[cfg(test)]
        LocationStreamPolicy::SetUnion => {
            let mut maximum = None;
            for filter in filters {
                match filter {
                    LocationFilter::Contig(expected) if expected == chrom => return None,
                    LocationFilter::Range {
                        chrom: expected,
                        end,
                        ..
                    } if expected == chrom => {
                        maximum = Some(maximum.map_or(*end, |current: usize| current.max(*end)));
                    }
                    _ => {}
                }
            }
            maximum
        }
    }
}

fn location_filter_end(filter: &LocationFilter, chrom: &str) -> Option<usize> {
    match filter {
        LocationFilter::Range {
            chrom: expected,
            end,
            ..
        } if expected == chrom => Some(*end),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_only_overlapping_locations_expand_to_independent_streams() {
        let filters = [
            LocationFilter::Range {
                chrom: "chr1".into(),
                start: 1,
                end: 100,
            },
            LocationFilter::Range {
                chrom: "chr1".into(),
                start: 51,
                end: 120,
            },
        ];
        assert_eq!(
            location_stream_groups(
                LocationStreamPolicy::IndependentLegacyStreams,
                &filters,
                "chr1",
                75,
            ),
            [0, 1]
        );
    }

    #[test]
    fn normative_set_union_collapses_overlapping_locations() {
        let filters = [
            LocationFilter::Contig("chr1".into()),
            LocationFilter::Contig("chr1".into()),
        ];
        assert_eq!(
            location_stream_groups(LocationStreamPolicy::SetUnion, &filters, "chr1", 75),
            [0]
        );
        assert!(
            location_stream_groups(LocationStreamPolicy::SetUnion, &filters, "chr2", 75).is_empty()
        );
    }

    #[test]
    fn normative_set_union_geometry_includes_later_range_end() {
        let filters = [
            LocationFilter::Range {
                chrom: "chr1".into(),
                start: 1,
                end: 100,
            },
            LocationFilter::Range {
                chrom: "chr1".into(),
                start: 1_000,
                end: 2_000,
            },
        ];
        assert_eq!(
            location_stream_final_end(LocationStreamPolicy::SetUnion, &filters, 0, "chr1"),
            Some(2_000)
        );
    }
}
