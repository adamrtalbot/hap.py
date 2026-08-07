//! Cohesive preprocessing responsibility.

use super::alleles::{
    calls_non_ref_allele, convert_gvcf_record, convert_somatic_record,
    finalize_somatic_for_pipeline, materialize_symbolic_deletion, trim_uncalled_non_ref,
};
use super::canonical::validate_record_reference;
use super::normalization::{normalize_bcftools_record, record_reference_matches};
use super::options::{add_legacy_chr_prefix, has_non_reference_genotype, passes_filters_only};
use super::{
    BlocksplitContigState, BlocksplitJob, BlocksplitObservation, BlocksplitSelection,
    LEGACY_MIN_BLOCK_VARIANTS,
};
use crate::adapters::vcf;
use crate::cli_compat::cli::{PreprocessArgs, SomaticGtMode};
use crate::domain::{Interval, RawVcfRecord};
use anyhow::Result;
use std::collections::HashSet;
use std::path::Path;

pub(super) fn collect_blocksplit_observations(
    records: &[RawVcfRecord],
    args: &PreprocessArgs,
    fixchr: bool,
    normalization_enabled: bool,
    somatic_mode: Option<SomaticGtMode>,
    somatic_sample_names: Option<&[String]>,
    reference_sequences: &std::collections::BTreeMap<String, String>,
    regions: Option<&[Interval]>,
    targets: Option<&[Interval]>,
    locations: Option<&[vcf::LocationFilter]>,
) -> Result<Vec<BlocksplitObservation>> {
    let input_path = Path::new(&args.input);
    let mut observations = Vec::new();
    let mut normalized_seen = HashSet::new();

    for source in records {
        let mut record = source.clone();
        if fixchr {
            record.chrom = add_legacy_chr_prefix(&record.chrom);
        }
        if (args.pass_only && !record.is_pass())
            || (!args.pass_only
                && !passes_filters_only(&record.filter, args.filters_only.as_deref()))
        {
            continue;
        }
        let effective_end = record.effective_end_pos(input_path)?;
        if !vcf::matches_interval_filters(
            &record.chrom,
            record.pos,
            effective_end,
            regions,
            targets,
            locations,
        ) {
            continue;
        }
        if args.convert_gvcf_to_vcf {
            if !convert_gvcf_record(&mut record) {
                continue;
            }
        } else {
            let drop_record = args.filter_nonref
                && (calls_non_ref_allele(&record)
                    || (normalization_enabled
                        && somatic_mode.is_none()
                        && !trim_uncalled_non_ref(&mut record)));
            if drop_record {
                continue;
            }
        }

        let blocksplit_pos = record.pos;
        let blocksplit_end = effective_end;
        let symbolic_deletion =
            normalization_enabled && record.alt_allele.split(',').any(|alt| alt == "<DEL>");
        if symbolic_deletion {
            let reference = reference_sequences
                .get(&record.chrom)
                .ok_or_else(|| anyhow::anyhow!("reference contig {} not found", record.chrom))?;
            materialize_symbolic_deletion(&mut record, effective_end, reference.as_bytes())?;
        }

        if args.bcftools_norm {
            let Some(reference) = reference_sequences.get(&record.chrom) else {
                continue;
            };
            if !record_reference_matches(&record, reference.as_bytes()) {
                continue;
            }
        } else if normalization_enabled {
            validate_record_reference(&record, reference_sequences)?;
        }

        let converted =
            if let (Some(mode), Some(sample_names)) = (somatic_mode, somatic_sample_names) {
                finalize_somatic_for_pipeline(
                    convert_somatic_record(&record, mode, sample_names),
                    mode,
                    normalization_enabled,
                )
            } else {
                vec![record]
            };

        for mut record in converted {
            if args.bcftools_norm {
                let Some(reference) = reference_sequences.get(&record.chrom) else {
                    continue;
                };
                normalize_bcftools_record(&mut record, reference.as_bytes());
                let key = (
                    record.chrom.clone(),
                    record.pos,
                    record.ref_allele.clone(),
                    record.alt_allele.clone(),
                );
                if !normalized_seen.insert(key) {
                    continue;
                }
            }
            let (pos, end) = if args.bcftools_norm && !symbolic_deletion {
                (record.pos, record.end_pos())
            } else {
                (blocksplit_pos, blocksplit_end)
            };
            let called = has_non_reference_genotype(&record);
            let location_groups = locations.map_or_else(
                || vec![0],
                |filters| {
                    super::compatibility::location_stream_groups(
                        super::LOCATION_STREAM_POLICY,
                        filters,
                        &record.chrom,
                        pos,
                    )
                },
            );
            observations.push(BlocksplitObservation {
                chrom: record.chrom,
                pos,
                end,
                called,
                location_groups,
            });
        }
    }
    Ok(observations)
}

pub(super) fn select_blocksplit_resets(
    observations: &[BlocksplitObservation],
    window_size: i64,
    block_count: usize,
    locations: Option<&[vcf::LocationFilter]>,
) -> BlocksplitSelection {
    select_blocksplit_resets_with_policy(
        observations,
        window_size,
        block_count,
        locations,
        super::LOCATION_STREAM_POLICY,
    )
}

pub(super) fn select_blocksplit_resets_with_policy(
    observations: &[BlocksplitObservation],
    window_size: i64,
    block_count: usize,
    locations: Option<&[vcf::LocationFilter]>,
    location_stream_policy: super::compatibility::LocationStreamPolicy,
) -> BlocksplitSelection {
    if block_count == 0 {
        return BlocksplitSelection::default();
    }
    let mut states = std::collections::BTreeMap::<(usize, String), BlocksplitContigState>::new();
    for (index, observation) in observations.iter().enumerate() {
        if !observation.called {
            continue;
        }
        for &location_group in &observation.location_groups {
            let state = states
                .entry((location_group, observation.chrom.clone()))
                .or_default();
            state.total_variants += 1;
            state.candidate_variants += 1;
            if state.last_called_end.is_some_and(|last_end| {
                observation.pos as i128 > last_end as i128 + window_size as i128
            }) && state.candidate_variants > LEGACY_MIN_BLOCK_VARIANTS
            {
                state.candidates.push((index, state.candidate_variants));
                state.candidate_variants = 0;
            }
            state.last_called_end = Some(
                state
                    .last_called_end
                    .map_or(observation.end, |last_end| last_end.max(observation.end)),
            );
        }
    }

    if states.is_empty() {
        return BlocksplitSelection::default();
    }

    let mut partition_resets = std::collections::BTreeMap::new();
    for (partition, state) in &states {
        let target_variants = LEGACY_MIN_BLOCK_VARIANTS.max(state.total_variants / block_count);
        let mut accumulated = 0usize;
        let mut selected = Vec::new();
        for &(index, candidate_variants) in &state.candidates {
            accumulated += candidate_variants;
            if accumulated > target_variants {
                selected.push(index);
                accumulated = 0;
            }
        }
        partition_resets.insert(partition.clone(), selected);
    }
    let jobs = build_blocksplit_jobs(
        observations,
        &states,
        &partition_resets,
        locations,
        location_stream_policy,
    );
    BlocksplitSelection { jobs: Some(jobs) }
}

pub(super) fn build_blocksplit_jobs(
    observations: &[BlocksplitObservation],
    states: &std::collections::BTreeMap<(usize, String), BlocksplitContigState>,
    partition_resets: &std::collections::BTreeMap<(usize, String), Vec<usize>>,
    locations: Option<&[vcf::LocationFilter]>,
    location_stream_policy: super::compatibility::LocationStreamPolicy,
) -> Vec<BlocksplitJob> {
    let mut jobs = Vec::new();
    for ((location_group, chrom), state) in states {
        let chrom_indices = observations
            .iter()
            .enumerate()
            .filter_map(|(index, observation)| (observation.chrom == *chrom).then_some(index))
            .collect::<Vec<_>>();
        let boundaries = partition_resets
            .get(&(*location_group, chrom.clone()))
            .map(Vec::as_slice)
            .unwrap_or_default();
        let location_start = locations.and_then(|_| {
            chrom_indices
                .iter()
                .copied()
                .find(|&index| observations[index].location_groups.contains(location_group))
        });

        // blocksplit's no-breakpoint fallback returns the complete contig,
        // even when the requested location covers only one part of it.
        if state.candidates.is_empty() || boundaries.is_empty() {
            let included_indices = chrom_indices.into_iter().collect::<HashSet<_>>();
            let reset_before_indices = location_start
                .filter(|index| *index > 0 && included_indices.contains(index))
                .into_iter()
                .collect();
            jobs.push(BlocksplitJob {
                included_indices,
                reset_before_indices,
            });
            continue;
        }

        let final_end = locations.and_then(|filters| {
            super::compatibility::location_stream_final_end(
                location_stream_policy,
                filters,
                *location_group,
                chrom,
            )
        });
        let mut block_start = None;
        let mut block_start_index = None;
        for &boundary_index in boundaries {
            let block_end = observations[boundary_index].pos;
            let included_indices = chrom_indices
                .iter()
                .copied()
                .filter(|&index| {
                    let pos = observations[index].pos;
                    block_start.is_none_or(|start| pos >= start) && pos < block_end
                })
                .collect::<HashSet<_>>();
            if !included_indices.is_empty() {
                let reset_before_indices = [block_start_index, location_start]
                    .into_iter()
                    .flatten()
                    .filter(|index| *index > 0 && included_indices.contains(index))
                    .collect();
                jobs.push(BlocksplitJob {
                    included_indices,
                    reset_before_indices,
                });
            }
            block_start = Some(block_end);
            block_start_index = Some(boundary_index);
        }
        let included_indices = chrom_indices
            .into_iter()
            .filter(|&index| {
                let pos = observations[index].pos;
                block_start.is_none_or(|start| pos >= start)
                    && final_end.is_none_or(|end| pos < end)
            })
            .collect::<HashSet<_>>();
        if !included_indices.is_empty() {
            let reset_before_indices = [block_start_index, location_start]
                .into_iter()
                .flatten()
                .filter(|index| *index > 0 && included_indices.contains(index))
                .collect();
            jobs.push(BlocksplitJob {
                included_indices,
                reset_before_indices,
            });
        }
    }
    jobs
}
