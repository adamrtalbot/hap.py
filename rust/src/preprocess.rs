use crate::cli::{PreprocessArgs, PreprocessGender, SomaticGtMode};
use crate::{fasta, partial_credit, variant_pipeline, vcf};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeSet, HashSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};

/// Per-sample left-shift boundary used by legacy hap.py's C++ preprocess
/// binary — indels can slide at most 1 kbp leftward before the pass stops.
const LEFT_SHIFT_WINDOW: usize = 1024;

/// `blocksplit` does not create a preprocessing boundary until the candidate
/// block contains more than its default minimum of 100 called variants.
const LEGACY_MIN_BLOCK_VARIANTS: usize = 100;
const LEGACY_MAX_BLOCKS: usize = 40;

#[derive(Clone, Debug)]
struct BlocksplitObservation {
    chrom: String,
    pos: usize,
    end: usize,
    called: bool,
    location_groups: Vec<usize>,
}

#[derive(Default)]
struct BlocksplitContigState {
    total_variants: usize,
    candidate_variants: usize,
    last_called_end: Option<usize>,
    candidates: Vec<(usize, usize)>,
}

#[derive(Default)]
struct BlocksplitSelection {
    /// `None` represents legacy's all-empty fallback, which processes the
    /// complete prepared stream instead of scheduling no jobs.
    jobs: Option<Vec<BlocksplitJob>>,
}

struct BlocksplitJob {
    included_indices: HashSet<usize>,
    /// Resets can land on records discarded later by VariantCallsOnly. They
    /// must therefore be applied immediately after selecting the prepared
    /// input record, before genotype-based filtering.
    reset_before_indices: HashSet<usize>,
}

#[cfg(test)]
impl BlocksplitSelection {
    fn reset_indices(&self) -> HashSet<usize> {
        self.jobs
            .iter()
            .flatten()
            .flat_map(|job| job.reset_before_indices.iter().copied())
            .collect()
    }

    fn included_indices(&self) -> Option<HashSet<usize>> {
        self.jobs.as_ref().map(|jobs| {
            jobs.iter()
                .flat_map(|job| job.included_indices.iter().copied())
                .collect()
        })
    }
}

pub fn run(args: PreprocessArgs) -> Result<()> {
    if args.version {
        println!("{}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    let mut logger = PreprocessLogger::new(&args)?;
    logger.info(&format!("Preprocessing {}", args.input))?;
    let output_path = if args.bcf && !args.output.ends_with(".bcf") {
        PathBuf::from(format!("{}.bcf", args.output))
    } else {
        PathBuf::from(&args.output)
    };
    require_output_parent(&output_path)?;
    let reference_path = resolve_reference(args.reference.as_deref())?;
    let (mut headers, records) = vcf::load_raw_vcf(Path::new(&args.input))?;
    require_vcf_sample(&headers)?;
    let reference_index = fasta::read_index(&reference_path)?;
    let reference_sequences = fasta::read_sequences(&reference_path)?;
    let reference_contigs: BTreeSet<String> = reference_index.keys().cloned().collect();
    // Legacy passes selectors to bcftools after its optional CHROM rewrite;
    // selector names themselves are never normalized against the reference.
    let literal_contigs = BTreeSet::new();
    let regions = args
        .regions_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &literal_contigs))
        .transpose()?;
    let targets = args
        .targets_bedfile
        .as_ref()
        .map(|path| vcf::load_bed(Path::new(path), &literal_contigs))
        .transpose()?;
    let locations = args
        .locations
        .as_deref()
        .map(|text| vcf::parse_locations(text, &literal_contigs))
        .transpose()?;

    if args.convert_gvcf_to_vcf {
        filter_gvcf_headers(&mut headers);
    }
    let string_format_fields = string_format_fields(&headers);
    let somatic_mode = resolve_somatic_mode(&args);
    let somatic_sample_names = somatic_mode.map(|_| somatic_info_sample_names(&headers));
    if let Some(sample_names) = somatic_sample_names.as_deref() {
        append_somatic_info_headers(&mut headers, sample_names);
    }
    let input_contigs = records.iter().map(|record| record.chrom.clone()).collect();
    let requested_fixchr = if args.no_fixchr {
        Some(false)
    } else {
        args.fixchr
    };
    let fixchr = resolve_fixchr(requested_fixchr, &reference_contigs, &input_contigs);
    let gender = resolve_gender(args.gender, &records);
    let leftshift = args.leftshift && !args.no_leftshift;
    let decompose = args.decompose && !args.no_decompose;
    let normalization_enabled = leftshift || decompose || gender == PreprocessGender::Male;
    let effective_threads = effective_thread_count(args.threads);
    let blocksplit_selection = if normalization_enabled && effective_threads > 1 {
        let observations = collect_blocksplit_observations(
            &records,
            &args,
            fixchr,
            normalization_enabled,
            somatic_mode,
            somatic_sample_names.as_deref(),
            &reference_sequences,
            regions.as_deref(),
            targets.as_deref(),
            locations.as_deref(),
        )?;
        select_blocksplit_resets(
            &observations,
            args.window_size,
            LEGACY_MAX_BLOCKS.min(effective_threads.saturating_mul(4)),
            locations.as_deref(),
        )
    } else {
        BlocksplitSelection::default()
    };

    let mut output = Vec::new();
    let job_count = blocksplit_selection.jobs.as_ref().map_or(1, Vec::len);
    for job_index in 0..job_count {
        let job = blocksplit_selection
            .jobs
            .as_ref()
            .map(|jobs| &jobs[job_index]);
        let mut normalized_seen = HashSet::new();
        // Each legacy partial-credit job owns an independent left-shift
        // boundary, so overlapping jobs must process their copies separately.
        let mut prev_end_by_chrom: std::collections::HashMap<String, usize> =
            std::collections::HashMap::new();
        let mut prepared_record_index = 0usize;
        for mut record in records.iter().cloned() {
            if fixchr {
                record.chrom = add_legacy_chr_prefix(&record.chrom);
            }
            if (args.pass_only && !record.is_pass())
                || (!args.pass_only
                    && !passes_filters_only(&record.filter, args.filters_only.as_deref()))
            {
                continue;
            }
            let effective_end = record.effective_end_pos(Path::new(&args.input))?;
            if !vcf::matches_interval_filters(
                &record.chrom,
                record.pos,
                effective_end,
                regions.as_deref(),
                targets.as_deref(),
                locations.as_deref(),
            ) {
                continue;
            }

            // Standalone pre defaults this field to true because legacy pre.py
            // accidentally leaves preprocess()'s default in force. hap.py,
            // however, explicitly forwards its own default=false. Honor the
            // resolved value here so germline keeps called <NON_REF> alleles unless
            // the user requests filtering, without changing standalone pre parity.
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

            // VariantReader.cpp treats a symbolic deletion as an empty ALT over
            // the complete INFO/END span. VariantWriter then adds the preceding
            // reference base so the emitted VCF contains an ordinary, anchored
            // deletion. Do this before REF validation: when END is present the
            // legacy reader intentionally ignores the input REF and rebuilds it
            // from the reference sequence.
            let symbolic_deletion = if normalization_enabled
                && record.alt_allele.split(',').any(|alt| alt == "<DEL>")
            {
                let reference = reference_sequences.get(&record.chrom).ok_or_else(|| {
                    anyhow::anyhow!("reference contig {} not found", record.chrom)
                })?;
                Some(materialize_symbolic_deletion(
                    &mut record,
                    effective_end,
                    reference.as_bytes(),
                )?)
            } else {
                None
            };

            if args.bcftools_norm {
                let Some(reference) = reference_sequences.get(&record.chrom) else {
                    continue;
                };
                if !record_reference_matches(&record, reference.as_bytes()) {
                    continue;
                }
            } else if normalization_enabled {
                validate_record_reference(&record, &reference_sequences)?;
            }

            let converted = if let (Some(mode), Some(sample_names)) =
                (somatic_mode, somatic_sample_names.as_deref())
            {
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
                let record_index = prepared_record_index;
                prepared_record_index += 1;
                if job.is_some_and(|job| !job.included_indices.contains(&record_index)) {
                    continue;
                }
                if job.is_some_and(|job| job.reset_before_indices.contains(&record_index)) {
                    prev_end_by_chrom.remove(&record.chrom);
                }
                if !normalization_enabled {
                    output.push(record);
                    continue;
                }
                if args.convert_gvcf_to_vcf {
                    ensure_missing_ad(&mut record);
                }
                // Allele-count INFO fields become stale after preprocessing splits
                // multi-allelics or decomposes complex variants. Legacy hap.py's C++
                // preprocess binary (via `VariantAlleleRemover` → `VariantCallsOnly`)
                // drops them so downstream tools recompute. We mirror that here on the
                // main code path so output stays deterministic for every caller, not
                // just the parity reference.
                strip_stale_info_keys(&mut record);

                // Legacy `VariantReader`/`VariantWriter` drops variant IDs (rsIDs) on
                // output — `preprocess` rebuilds records from CHROM/POS/REF/ALT only,
                // so the ID column always becomes `.`.
                record.id = ".".to_string();

                // Upper-case REF/ALT so soft-masked lowercase bases from the reference
                // come out as the canonical uppercase form legacy emits. We already
                // tolerate case when validating; now we normalise on output.
                record.ref_allele = record.ref_allele.to_ascii_uppercase();
                record.alt_allele = uppercase_alleles_preserving_breakends(&record.alt_allele);
                let import_failed = materialize_unsupported_import_failure(&mut record);

                // Legacy `VariantWriter` emits FILTER=`.` for PASS records (htslib's
                // `bcf_update_filter` treats an empty filter vector as `.`). Normalise
                // PASS to `.` so downstream byte-diff matches.
                if record.filter == "PASS" {
                    record.filter = ".".to_string();
                }

                // Legacy stores INFO as an ordered map keyed alphabetically by tag
                // (std::map<std::string,...>) before serialising — htslib then writes
                // entries in that order. Sort our INFO tags the same way so output
                // byte-matches.
                sort_info_keys(&mut record);

                // Legacy represents PL internally as a single integer per sample (the
                // `v.asInt()` path in `VariantWriter.cpp` line 563). When bcftools
                // emits the record the array is truncated to the last stored value,
                // which for biallelic diploid sites is the HOM_ALT likelihood. We
                // reproduce that truncation here so SAMPLE cells byte-match.
                collapse_pl_to_last_value(&mut record);

                // Multi-allelic indel decomposition + primitive splitting.
                //
                // Replaces the previous `should_split_multi_allelic_indel +
                // split_multi_allelic` branch: `variant_pipeline::primitive_split`
                // implements the legacy `VariantPrimitiveSplitter` (stage 5 of
                // `VariantInput.cpp`) plus the `aggregate_hetalt` re-merge that
                // brings same-position primitives back to a multi-allelic shape.
                // Same-direction insertions (`T → TG,TTG`) re-merge into a single
                // record; mixed-direction or different-length deletions fan out
                // into per-primitive records anchored at their canonical position.
                // Insert the ADO field BEFORE primitive_split *and* before any GT
                // normalisation. Legacy computes `ad_other` (= ADO) from the
                // *original* GT+AD inside VariantReader, then preserves that value
                // across both the half-call AlleleSplitter and the per-allele
                // PrimitiveSplitter. Computing ADO after splitting OR after the
                // haploid → homalt expansion below would widen the "called" set
                // if computed after expansion (input `GT=1` ad=[1,23] becomes
                // `1/1` ad=[1,23] → ADO=AD[0]=1, the unused ref depth, which
                // matches legacy). Computing ADO from the original GT=1 first
                // preserves that ref-depth signal correctly.
                ensure_missing_ad(&mut record);
                insert_ado_format(&mut record);
                ensure_missing_dp(&mut record);

                // The legacy C++ Variant representation has MAX_GT=2. During
                // active preprocessing, wider calls are converted to no-calls
                // before VariantCallsOnly removes their now-uncalled record.
                // This is observable in hap.py's vcfeval handoff for triploid
                // and tetraploid query records.
                if somatic_mode.is_none() {
                    mask_genotypes_wider_than_diploid(&mut record);
                }

                // VariantCallsOnly removes ALT alleles that no sample calls before
                // primitive decomposition. Besides reducing ordinary multi-allelic
                // records, this prevents the primitive splitter from emitting a
                // hom-ref sibling for every uncalled component of a complex allele.
                // Capture ADO first: legacy derives it from the original GT/AD and
                // then projects GT and AD onto the retained alleles.
                if somatic_mode.is_none()
                    && !args.convert_gvcf_to_vcf
                    && !import_failed
                    && !retain_called_alternates(&mut record)
                {
                    continue;
                }

                // Normalise haploid GTs to the legacy het / hom shape. The C++
                // VariantAlleleSplitter treats a single haploid alt call (ngt == 1
                // with gt[0] > 0) as a het half-call and emits `0/1` after the
                // half-call merge — see `VariantAlleleSplitter.cpp:180-227`. Mirror
                // that here so that haploid input lines (`GT=1` on autosomes) come
                // out byte-identical to legacy without a chrX/Y-specific shim.
                if somatic_mode.is_none() {
                    normalise_haploid_genotypes(&mut record, gender == PreprocessGender::Male);
                }

                // Capture before record is potentially consumed by vec![record].
                let record_chrom = record.chrom.clone();
                let orig_end = record.pos + record.ref_allele.len().max(1) - 1;
                let previous_end = prev_end_by_chrom.get(&record_chrom).copied().unwrap_or(0);
                let prev_end = previous_end;

                let source_records = if matches!(
                    symbolic_deletion,
                    Some(SymbolicDeletionMaterialization::LeadingAnchor)
                ) {
                    split_called_alleles(&record)
                } else {
                    vec![MaterializedAlleleRecord {
                        record,
                        reverse_hetalt_samples: Vec::new(),
                    }]
                };
                let mut emitted_groups = Vec::with_capacity(source_records.len());
                for source in source_records {
                    let reverse_hetalt_samples = source.reverse_hetalt_samples;
                    let source = source.record;
                    let mut emitted =
                        if decompose && !source.alt_allele.split(',').any(is_symbolic_allele) {
                            if let Some(reference) = reference_sequences.get(&source.chrom) {
                                variant_pipeline::primitive_split_with_floor(
                                    &source,
                                    reference.as_bytes(),
                                    prev_end,
                                )
                            } else {
                                vec![source]
                            }
                        } else {
                            vec![source]
                        };
                    if decompose {
                        for record in &mut emitted {
                            restore_legacy_hetalt_orientation(record, &reverse_hetalt_samples);
                        }
                    }
                    emitted_groups.push(emitted);
                }
                for emitted in emitted_groups {
                    // Only leftshift records that came through primitive_split
                    // unchanged. Fanned-out primitives are already canonical.
                    let leftshift_eligible = emitted.len() == 1;
                    for mut split in emitted {
                        if leftshift
                            && leftshift_eligible
                            && !split.alt_allele.contains(',')
                            && split.alt_allele != "."
                            && !split.alt_allele.is_empty()
                            && !is_symbolic_allele(&split.alt_allele)
                            && let Some(reference) = reference_sequences.get(&split.chrom)
                        {
                            apply_left_shift(&mut split, reference.as_bytes(), prev_end);
                        }
                        canonicalize_multi_allelic_order(&mut split);
                        canonicalize_legacy_genotypes(&mut split);
                        if split.qual.is_empty() || split.qual == "." {
                            split.qual = "0".to_string();
                        }
                        blank_secondary_sample_annotations(
                            &mut split,
                            args.bcf,
                            &string_format_fields,
                        );
                        // Canonical FORMAT ordering: GT → AD/ADO/DP (fixed) → other
                        // integer-typed fields alphabetical → float-typed alphabetical →
                        // string-typed alphabetical. Mirrors the `int_fmts/float_fmts/
                        // string_fmts` loop order in `VariantWriter.cpp` combined with
                        // dynamic per-value type detection.
                        reorder_format_fields(&mut split);
                        output.push(split);
                    }
                }
                // Advance the per-chromosome boundary so the next variant cannot
                // left-shift into this record's reference span.
                prev_end_by_chrom
                    .entry(record_chrom)
                    .and_modify(|e| *e = (*e).max(orig_end))
                    .or_insert(orig_end);
            }
        }
    }

    if somatic_mode.is_some() {
        rewrite_header_for_single_sample(&mut headers);
    }

    if fixchr {
        ensure_emitted_contig_headers(&mut headers, &output);
    }

    if normalization_enabled {
        sort_normalized_records(&mut output);
        headers = canonicalize_legacy_headers(&headers);
    } else {
        ensure_pass_filter_header(&mut headers);
    }

    vcf::write_raw_vcf(&output_path, &headers, &output)?;
    if output_path
        .extension()
        .and_then(|extension| extension.to_str())
        == Some("vcf")
    {
        bail!(
            "plain VCF output {} cannot be indexed; legacy pre.py exits unsuccessfully",
            output_path.display()
        );
    }
    logger.info(&format!(
        "Wrote {} records to {}",
        output.len(),
        output_path.display()
    ))?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn collect_blocksplit_observations(
    records: &[vcf::RawVcfRecord],
    args: &PreprocessArgs,
    fixchr: bool,
    normalization_enabled: bool,
    somatic_mode: Option<SomaticGtMode>,
    somatic_sample_names: Option<&[String]>,
    reference_sequences: &std::collections::BTreeMap<String, String>,
    regions: Option<&[vcf::BedInterval]>,
    targets: Option<&[vcf::BedInterval]>,
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
                    filters
                        .iter()
                        .enumerate()
                        .filter_map(|(index, filter)| {
                            filter.matches(&record.chrom, pos).then_some(index)
                        })
                        .collect()
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

fn select_blocksplit_resets(
    observations: &[BlocksplitObservation],
    window_size: i64,
    block_count: usize,
    locations: Option<&[vcf::LocationFilter]>,
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
    let jobs = build_blocksplit_jobs(observations, &states, &partition_resets, locations);
    BlocksplitSelection { jobs: Some(jobs) }
}

fn build_blocksplit_jobs(
    observations: &[BlocksplitObservation],
    states: &std::collections::BTreeMap<(usize, String), BlocksplitContigState>,
    partition_resets: &std::collections::BTreeMap<(usize, String), Vec<usize>>,
    locations: Option<&[vcf::LocationFilter]>,
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

        let final_end = locations
            .and_then(|filters| filters.get(*location_group))
            .and_then(|filter| match filter {
                vcf::LocationFilter::Range {
                    chrom: expected,
                    end,
                    ..
                } if expected == chrom => Some(*end),
                _ => None,
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

/// Primitive splitting can move one allele to the trailing edge of an
/// overlapping input record. The following input record may start inside
/// that span, so emission order is no longer guaranteed to be coordinate
/// order even when the source VCF was indexed. Legacy's location aggregator
/// restores position order before writing; do the same while preserving the
/// source contig order and stable order among records at the same position.
fn sort_normalized_records(records: &mut [vcf::RawVcfRecord]) {
    let mut contig_ranks = std::collections::HashMap::new();
    let mut next_rank = 0usize;
    for record in records.iter() {
        contig_ranks.entry(record.chrom.clone()).or_insert_with(|| {
            let rank = next_rank;
            next_rank += 1;
            rank
        });
    }
    records.sort_by(|left, right| {
        contig_ranks[&left.chrom]
            .cmp(&contig_ranks[&right.chrom])
            .then(left.pos.cmp(&right.pos))
    });
}

fn effective_thread_count(threads: Option<usize>) -> usize {
    let available_threads = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1);
    effective_thread_count_with_available(threads, available_threads)
}

fn effective_thread_count_with_available(
    threads: Option<usize>,
    available_threads: usize,
) -> usize {
    threads.unwrap_or_else(|| available_threads.max(1))
}

fn has_non_reference_genotype(record: &vcf::RawVcfRecord) -> bool {
    let Some(format) = record.format.as_deref() else {
        return false;
    };
    let Some(gt_index) = format.split(':').position(|field| field == "GT") else {
        return false;
    };
    record.samples.iter().any(|sample| {
        sample
            .split(':')
            .nth(gt_index)
            .unwrap_or(".")
            .split(['/', '|'])
            .any(|allele| allele.parse::<usize>().is_ok_and(|allele| allele > 0))
    })
}

struct PreprocessLogger {
    verbose: bool,
    quiet: bool,
    file: Option<File>,
}

impl PreprocessLogger {
    fn new(args: &PreprocessArgs) -> Result<Self> {
        let file = args
            .logfile
            .as_deref()
            .map(|path| {
                if let Some(parent) = Path::new(path)
                    .parent()
                    .filter(|parent| !parent.as_os_str().is_empty())
                {
                    fs::create_dir_all(parent)?;
                }
                File::create(path).with_context(|| format!("failed to create logfile {path}"))
            })
            .transpose()?;
        Ok(Self {
            verbose: args.verbose,
            quiet: args.quiet,
            file,
        })
    }

    fn info(&mut self, message: &str) -> Result<()> {
        if !self.verbose || self.quiet {
            return Ok(());
        }
        if let Some(file) = &mut self.file {
            writeln!(file, "INFO {message}")?;
            file.flush()?;
        } else {
            eprintln!("[I] {message}");
        }
        Ok(())
    }
}

/// The Perl prefixing stage in legacy `pre.py` changes record CHROM values
/// without rewriting the input declarations. The following `bcftools view`
/// pass notices each newly used sequence and appends a length-less contig
/// declaration. Mirror that repair while leaving existing declarations and
/// their order untouched.
fn ensure_emitted_contig_headers(headers: &mut Vec<String>, records: &[vcf::RawVcfRecord]) {
    let mut declared: BTreeSet<String> = headers
        .iter()
        .filter_map(|line| {
            line.strip_prefix("##contig=<ID=")
                .and_then(|body| body.split([',', '>']).next())
                .map(str::to_string)
        })
        .collect();
    let mut additions = Vec::new();
    for record in records {
        if declared.insert(record.chrom.clone()) {
            additions.push(format!("##contig=<ID={}>", record.chrom));
        }
    }
    let insert_at = headers
        .iter()
        .position(|line| line.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    headers.splice(insert_at..insert_at, additions);
}

fn ensure_pass_filter_header(headers: &mut Vec<String>) {
    if headers
        .iter()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"))
    {
        return;
    }
    let index = headers
        .iter()
        .position(|line| line.starts_with("##fileformat="))
        .map_or(0, |index| index + 1);
    headers.insert(
        index,
        "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
    );
}

fn resolve_reference(explicit: Option<&str>) -> Result<PathBuf> {
    let hg19 = std::env::var_os("HG19").map(PathBuf::from);
    let hgref = std::env::var_os("HGREF").map(PathBuf::from);
    resolve_reference_candidates(
        explicit.map(Path::new),
        hg19.as_deref(),
        hgref.as_deref(),
        Path::new("/opt/hap.py-data/hg19.fa"),
    )
}

fn require_output_parent(output: &Path) -> Result<()> {
    if let Some(parent) = output
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        && !parent.is_dir()
    {
        bail!("output parent does not exist: {}", parent.display());
    }
    Ok(())
}

fn require_vcf_sample(headers: &[String]) -> Result<()> {
    let has_sample = headers
        .iter()
        .find(|line| line.starts_with("#CHROM"))
        .is_some_and(|line| {
            line.split('\t')
                .nth(9)
                .is_some_and(|sample| !sample.is_empty())
        });
    if !has_sample {
        bail!("input VCF has no samples");
    }
    Ok(())
}

fn resolve_reference_candidates(
    explicit: Option<&Path>,
    hg19: Option<&Path>,
    hgref: Option<&Path>,
    fallback: &Path,
) -> Result<PathBuf> {
    if let Some(path) = explicit {
        return Ok(path.to_path_buf());
    }
    [hg19, hgref, Some(fallback)]
        .into_iter()
        .flatten()
        .find(|path| path.is_file())
        .map(Path::to_path_buf)
        .context("no reference file found; pass --reference or set HG19/HGREF")
}

fn has_chr_prefix(contigs: &BTreeSet<String>) -> Option<bool> {
    let plain = (0..23)
        .map(|value| value.to_string())
        .chain(["X".into(), "Y".into(), "MT".into()])
        .filter(|name| contigs.contains(name))
        .count();
    let prefixed = (0..23)
        .map(|value| format!("chr{value}"))
        .chain(["chrX".into(), "chrY".into(), "chrM".into()])
        .filter(|name| contigs.contains(name))
        .count();
    match prefixed.cmp(&plain) {
        std::cmp::Ordering::Greater => Some(true),
        std::cmp::Ordering::Less => Some(false),
        std::cmp::Ordering::Equal => None,
    }
}

fn resolve_fixchr(
    requested: Option<bool>,
    reference_contigs: &BTreeSet<String>,
    input_contigs: &BTreeSet<String>,
) -> bool {
    requested.unwrap_or_else(|| {
        has_chr_prefix(reference_contigs) == Some(true)
            && has_chr_prefix(input_contigs) == Some(false)
    })
}

fn add_legacy_chr_prefix(chrom: &str) -> String {
    if chrom == "chrMT" {
        return "chrM".to_string();
    }
    if chrom.starts_with("chr") {
        return chrom.to_string();
    }
    let first = chrom.as_bytes().first().copied();
    if !matches!(first, Some(b'0'..=b'9' | b'X' | b'Y' | b'M')) {
        return chrom.to_string();
    }
    if chrom == "MT" || chrom == "M" {
        "chrM".to_string()
    } else {
        format!("chr{chrom}")
    }
}

fn passes_filters_only(filter: &str, filters_only: Option<&str>) -> bool {
    let Some(filters_only) = filters_only.filter(|value| !value.is_empty()) else {
        return true;
    };
    if filter.is_empty() || matches!(filter, "." | "PASS") {
        return true;
    }
    let excluded: BTreeSet<&str> = filters_only.split(',').collect();
    filter.split(';').any(|name| !excluded.contains(name))
}

fn resolve_gender(requested: PreprocessGender, records: &[vcf::RawVcfRecord]) -> PreprocessGender {
    if requested != PreprocessGender::Auto {
        return requested;
    }
    let mut haploid_x = false;
    let mut diploid_x = false;
    for record in records
        .iter()
        .filter(|record| matches!(record.chrom.as_str(), "X" | "chrX" | "chrx"))
    {
        let Some(format) = record.format.as_deref() else {
            continue;
        };
        let Some(gt_index) = format.split(':').position(|field| field == "GT") else {
            continue;
        };
        for sample in &record.samples {
            let gt = sample.split(':').nth(gt_index).unwrap_or(".");
            // vcfcheck classifies ploidy from the encoded GT vector length,
            // including missing slots. Thus `./1` has ngt == 2 and unequal
            // alleles (-1 and 1), so it is diploid rather than haploid.
            let alleles: Vec<&str> = gt.split(['/', '|']).collect();
            if alleles.len() == 1 {
                haploid_x = true;
            } else if alleles.len() > 2 || (alleles.len() == 2 && alleles[0] != alleles[1]) {
                diploid_x = true;
            }
        }
    }
    if haploid_x && !diploid_x {
        PreprocessGender::Male
    } else {
        PreprocessGender::Female
    }
}

pub(crate) fn infer_gender(path: &Path) -> Result<PreprocessGender> {
    let (_, records) = vcf::load_raw_vcf(path)?;
    Ok(resolve_gender(PreprocessGender::Auto, &records))
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SymbolicDeletionMaterialization {
    LeadingAnchor,
    ContigStart,
}

fn materialize_symbolic_deletion(
    record: &mut vcf::RawVcfRecord,
    end: usize,
    reference: &[u8],
) -> Result<SymbolicDeletionMaterialization> {
    if end < record.pos {
        bail!(
            "symbolic deletion {}:{} has END={} before POS",
            record.chrom,
            record.pos,
            end
        );
    }

    if record.pos == 1 {
        // VariantWriter cannot add a preceding anchor at the start of a
        // contig. The pinned legacy binary instead expands REF through END,
        // keeps `<DEL>` symbolic, and consumes END from INFO.
        let ref_bases = reference.get(..end).ok_or_else(|| {
            anyhow::anyhow!(
                "symbolic deletion {}:1-{} extends beyond the reference sequence",
                record.chrom,
                end
            )
        })?;
        record.ref_allele = String::from_utf8_lossy(ref_bases).to_ascii_uppercase();
        record.alt_allele = record.alt_allele.to_ascii_uppercase();
        record.info = remove_info_field(&record.info, "END");
        return Ok(SymbolicDeletionMaterialization::ContigStart);
    }

    let (start, stop, alt_index) = {
        // Normal case: use the base immediately before the deleted interval
        // as the shared VCF anchor. `end` is 1-based inclusive, and therefore
        // also the exclusive byte index for this zero-based slice.
        (record.pos - 2, end, record.pos - 2)
    };
    let ref_bases = reference.get(start..stop).ok_or_else(|| {
        anyhow::anyhow!(
            "symbolic deletion {}:{}-{} extends beyond the reference sequence",
            record.chrom,
            record.pos,
            end
        )
    })?;
    let anchor = reference.get(alt_index..alt_index + 1).ok_or_else(|| {
        anyhow::anyhow!(
            "symbolic deletion {}:{}-{} has no reference anchor",
            record.chrom,
            record.pos,
            end
        )
    })?;

    let anchor = String::from_utf8_lossy(anchor).to_ascii_uppercase();
    let materialized_alts = record
        .alt_allele
        .split(',')
        .map(|alt| {
            if alt == "<DEL>" {
                anchor.clone()
            } else if is_symbolic_allele(alt) {
                alt.to_ascii_uppercase()
            } else {
                format!("{anchor}{}", alt.to_ascii_uppercase())
            }
        })
        .collect::<Vec<_>>();

    record.pos -= 1;
    record.ref_allele = String::from_utf8_lossy(ref_bases).to_ascii_uppercase();
    record.alt_allele = materialized_alts.join(",");
    record.info = remove_info_field(&record.info, "END");
    Ok(SymbolicDeletionMaterialization::LeadingAnchor)
}

fn is_symbolic_allele(alt: &str) -> bool {
    alt == "*" || (alt.starts_with('<') && alt.ends_with('>'))
}

/// Project a called multi-allelic record into biallelic records. Symbolic
/// deletions need this even without primitive decomposition; phased calls need
/// independent normalization streams because legacy never re-aggregates them.
struct MaterializedAlleleRecord {
    record: vcf::RawVcfRecord,
    reverse_hetalt_samples: Vec<bool>,
}

fn split_called_alleles(record: &vcf::RawVcfRecord) -> Vec<MaterializedAlleleRecord> {
    let alts = record.alt_allele.split(',').collect::<Vec<_>>();
    if alts.len() <= 1 {
        return vec![MaterializedAlleleRecord {
            record: record.clone(),
            reverse_hetalt_samples: Vec::new(),
        }];
    }
    let format_keys = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect::<Vec<_>>())
        .unwrap_or_default();
    let gt_index = format_keys.iter().position(|field| *field == "GT");
    let ad_index = format_keys.iter().position(|field| *field == "AD");

    alts.iter()
        .enumerate()
        .filter_map(|(alt_index, alt)| {
            let target = alt_index + 1;
            let called = gt_index.is_some_and(|gt_index| {
                record.samples.iter().any(|sample| {
                    sample
                        .split(':')
                        .nth(gt_index)
                        .is_some_and(|gt| genotype_calls_allele(gt, target))
                })
            });
            if !called {
                return None;
            }

            let mut split = record.clone();
            split.alt_allele = (*alt).to_string();
            let mut reverse_hetalt_samples = Vec::with_capacity(split.samples.len());
            for sample in &mut split.samples {
                let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
                if let Some(gt_index) = gt_index
                    && let Some(gt) = cells.get_mut(gt_index)
                {
                    reverse_hetalt_samples.push(genotype_is_hetalt_without_ref(gt, target));
                    *gt = project_split_genotype(gt, target);
                } else {
                    reverse_hetalt_samples.push(false);
                }
                if let Some(ad_index) = ad_index
                    && let Some(ad) = cells.get_mut(ad_index)
                {
                    *ad = project_split_ad(ad, target);
                }
                *sample = cells.join(":");
            }
            Some(MaterializedAlleleRecord {
                record: split,
                reverse_hetalt_samples,
            })
        })
        .collect()
}

fn genotype_calls_allele(gt: &str, target: usize) -> bool {
    gt.split(['/', '|'])
        .any(|allele| allele.parse::<usize>().ok() == Some(target))
}

fn genotype_is_hetalt_without_ref(gt: &str, target: usize) -> bool {
    let alleles = gt
        .split(['/', '|'])
        .filter_map(|allele| allele.parse::<usize>().ok())
        .collect::<Vec<_>>();
    gt.contains('/')
        && target > 1
        && alleles.contains(&target)
        && alleles.iter().all(|allele| *allele > 0)
        && alleles.iter().any(|allele| *allele != target)
}

fn restore_legacy_hetalt_orientation(record: &mut vcf::RawVcfRecord, reverse_samples: &[bool]) {
    let Some(gt_index) = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|field| field == "GT"))
    else {
        return;
    };
    for (sample, reverse) in record.samples.iter_mut().zip(reverse_samples) {
        if !reverse {
            continue;
        }
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        if let Some(gt) = cells.get_mut(gt_index) {
            match gt.as_str() {
                "0/1" => *gt = "1/0".to_string(),
                "0|1" => *gt = "1|0".to_string(),
                _ => {}
            }
        }
        *sample = cells.join(":");
    }
}

fn project_split_genotype(gt: &str, target: usize) -> String {
    let separator = if gt.contains('|') { '|' } else { '/' };
    let alleles = gt.split(['/', '|']).collect::<Vec<_>>();
    if alleles.contains(&".") {
        return gt.to_string();
    }
    let called = alleles
        .iter()
        .filter(|allele| allele.parse::<usize>().ok() == Some(target))
        .count();
    match called {
        0 => vec!["0"; alleles.len()].join(&separator.to_string()),
        count if count == alleles.len() => vec!["1"; alleles.len()].join(&separator.to_string()),
        _ => ["0", "1"].join(&separator.to_string()),
    }
}

fn project_split_ad(ad: &str, target: usize) -> String {
    if ad == "." || ad.is_empty() {
        return ad.to_string();
    }
    let values = ad.split(',').collect::<Vec<_>>();
    if values.len() <= 2 {
        return ad.to_string();
    }
    format!(
        "{},{}",
        values.first().copied().unwrap_or("0"),
        values.get(target).copied().unwrap_or("0")
    )
}

fn remove_info_field(info: &str, key: &str) -> String {
    let retained = info
        .split(';')
        .filter(|entry| {
            !entry.is_empty()
                && *entry != "."
                && entry
                    .split_once('=')
                    .map_or(*entry != key, |(name, _)| name != key)
        })
        .collect::<Vec<_>>();
    if retained.is_empty() {
        ".".to_string()
    } else {
        retained.join(";")
    }
}

fn resolve_somatic_mode(args: &PreprocessArgs) -> Option<SomaticGtMode> {
    args.set_gt.or(args.somatic.then_some(SomaticGtMode::Half))
}

/// Mirror `remove_nonref_gt_variants.py`: only the final ALT is treated as
/// `<NON_REF>`, and a record is dropped when any sample's first cell calls its
/// allele index. The legacy script intentionally assumes GT is the first FORMAT
/// cell, so this helper does too.
fn calls_non_ref_allele(record: &vcf::RawVcfRecord) -> bool {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    if alts.last() != Some(&"<NON_REF>") {
        return false;
    }
    let non_ref_index = alts.len();
    record.samples.iter().any(|sample| {
        sample
            .split(':')
            .next()
            .unwrap_or_default()
            .split(['/', '|'])
            .any(|token| token.parse::<usize>().ok() == Some(non_ref_index))
    })
}

fn trim_uncalled_non_ref(record: &mut vcf::RawVcfRecord) -> bool {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    if alts.last() != Some(&"<NON_REF>") {
        return record.alt_allele != ".";
    }
    let concrete = &alts[..alts.len().saturating_sub(1)];
    if concrete.is_empty() {
        return false;
    }
    let mapping: Vec<usize> = (0..=alts.len())
        .map(|index| if index < alts.len() { index } else { 0 })
        .collect();
    let format_keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect())
        .unwrap_or_default();
    let gt_index = format_keys.iter().position(|key| *key == "GT");
    let ad_index = format_keys.iter().position(|key| *key == "AD");
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        if let Some(index) = gt_index
            && let Some(gt) = cells.get_mut(index)
        {
            *gt = remap_gt(gt, &mapping);
        }
        if let Some(index) = ad_index
            && let Some(ad) = cells.get_mut(index)
        {
            let mut depths: Vec<&str> = ad.split(',').collect();
            if depths.len() == alts.len() + 1 {
                depths.pop();
                *ad = depths.join(",");
            }
        }
        *sample = cells.join(":");
    }
    record.alt_allele = concrete.join(",");
    true
}

/// Retain the union of non-reference alleles called by all samples, remapping
/// allele-indexed GT and AD cells. Returns false when the record has no called
/// ALT and is therefore not a variant after `VariantCallsOnly`.
fn retain_called_alternates(record: &mut vcf::RawVcfRecord) -> bool {
    let alts = record.alt_allele.split(',').collect::<Vec<_>>();
    if alts.is_empty() || alts == ["."] {
        return false;
    }
    let format_keys = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect::<Vec<_>>())
        .unwrap_or_default();
    let Some(gt_index) = format_keys.iter().position(|field| *field == "GT") else {
        return false;
    };
    let ad_index = format_keys.iter().position(|field| *field == "AD");
    let mut called = BTreeSet::new();
    for sample in &record.samples {
        if let Some(gt) = sample.split(':').nth(gt_index) {
            called.extend(
                gt.split(['/', '|'])
                    .filter_map(|allele| allele.parse::<usize>().ok())
                    // `VariantCallsOnly` does not treat a spanning-deletion
                    // (`*`) as an emitted variant. Its call is remapped to
                    // reference when a concrete ALT remains, or the whole
                    // record is discarded when it is the only ALT.
                    .filter(|allele| {
                        *allele > 0 && *allele <= alts.len() && alts[*allele - 1] != "*"
                    }),
            );
        }
    }
    if called.is_empty() {
        return false;
    }

    let mut mapping = vec![0usize; alts.len() + 1];
    let mut retained_alts = Vec::with_capacity(called.len());
    let mut retained_old_indices = Vec::with_capacity(called.len());
    for (alt_index, alt) in alts.iter().enumerate() {
        let old_index = alt_index + 1;
        if called.contains(&old_index) {
            mapping[old_index] = retained_alts.len() + 1;
            retained_alts.push(*alt);
            retained_old_indices.push(old_index);
        }
    }

    for sample in &mut record.samples {
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        if let Some(gt) = cells.get_mut(gt_index) {
            *gt = remap_gt(gt, &mapping);
            // VariantAlleleRemover canonicalizes an unphased REF/ALT call
            // after projecting removed alleles. The regular writer path only
            // sees this orientation for phased genotypes, so do it here.
            if *gt == "1/0" {
                *gt = "0/1".to_string();
            }
        }
        if let Some(ad_index) = ad_index
            && let Some(ad) = cells.get_mut(ad_index)
        {
            let depths = ad.split(',').collect::<Vec<_>>();
            if depths.len() == alts.len() + 1 {
                let mut projected = Vec::with_capacity(retained_old_indices.len() + 1);
                projected.push(depths[0]);
                projected.extend(
                    retained_old_indices
                        .iter()
                        .map(|old_index| depths[*old_index]),
                );
                *ad = projected.join(",");
            }
        }
        *sample = cells.join(":");
    }
    record.alt_allele = retained_alts.join(",");
    true
}

fn mask_genotypes_wider_than_diploid(record: &mut vcf::RawVcfRecord) {
    let Some(gt_index) = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|field| field == "GT"))
    else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        let Some(gt) = cells.get_mut(gt_index) else {
            continue;
        };
        if gt.split(['/', '|']).count() > 2 {
            *gt = ".".to_string();
        }
        *sample = cells.join(":");
    }
}

/// Port the standalone gVCF conversion pipeline:
/// `N_ALT >= 2` → keep only GT/DP/GQ → trim uncalled alleles → exclude
/// `<NON_REF>`. Returns false when the record is no longer a variant.
fn convert_gvcf_record(record: &mut vcf::RawVcfRecord) -> bool {
    let alts: Vec<String> = record.alt_allele.split(',').map(str::to_string).collect();
    if alts.len() < 2 {
        return false;
    }

    let format_keys: Vec<String> = record
        .format
        .as_deref()
        .map(|format| format.split(':').map(str::to_string).collect())
        .unwrap_or_default();
    let gt_index = format_keys.iter().position(|key| key == "GT");
    let mut called = BTreeSet::new();
    if let Some(gt_index) = gt_index {
        for sample in &record.samples {
            if let Some(gt) = sample.split(':').nth(gt_index) {
                called.extend(
                    gt.split(['/', '|'])
                        .filter_map(|token| token.parse::<usize>().ok())
                        .filter(|index| *index > 0),
                );
            }
        }
    }

    // If `<NON_REF>` remains called after allele trimming, the final legacy
    // bcftools exclusion removes the whole record.
    if alts
        .iter()
        .enumerate()
        .any(|(index, alt)| alt == "<NON_REF>" && called.contains(&(index + 1)))
    {
        return false;
    }

    let mut kept = Vec::new();
    let mut mapping = vec![0usize; alts.len() + 1];
    for (index, alt) in alts.iter().enumerate() {
        let old_index = index + 1;
        if alt != "<NON_REF>" && called.contains(&old_index) {
            mapping[old_index] = kept.len() + 1;
            kept.push(alt.clone());
        }
    }
    let retained_indices: Vec<usize> = format_keys
        .iter()
        .enumerate()
        .filter_map(|(index, key)| matches!(key.as_str(), "GT" | "DP" | "GQ").then_some(index))
        .collect();
    let retained_keys: Vec<String> = retained_indices
        .iter()
        .map(|index| format_keys[*index].clone())
        .collect();
    for sample in &mut record.samples {
        let cells: Vec<&str> = sample.split(':').collect();
        let mut retained: Vec<String> = retained_indices
            .iter()
            .map(|index| cells.get(*index).copied().unwrap_or(".").to_string())
            .collect();
        if let Some(new_gt_index) = retained_keys.iter().position(|key| key == "GT") {
            retained[new_gt_index] = remap_gt(&retained[new_gt_index], &mapping);
        }
        *sample = retained.join(":");
    }

    record.alt_allele = if kept.is_empty() {
        ".".to_string()
    } else {
        kept.join(",")
    };
    record.info = ".".to_string();
    record.format = (!retained_keys.is_empty()).then(|| retained_keys.join(":"));
    true
}

fn filter_gvcf_headers(headers: &mut Vec<String>) {
    headers.retain(|line| {
        if line.starts_with("##INFO=<ID=") {
            return false;
        }
        if let Some(body) = line.strip_prefix("##FORMAT=<ID=") {
            let id = body.split([',', '>']).next().unwrap_or_default();
            return matches!(id, "GT" | "DP" | "GQ");
        }
        true
    });
}

fn ensure_missing_ad(record: &mut vcf::RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    let keys: Vec<&str> = format.split(':').collect();
    if keys.contains(&"AD") {
        return;
    }
    let alt_count = record
        .alt_allele
        .split(',')
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .count();
    let missing_ad = vec!["."; alt_count + 1].join(",");
    record.format = Some(format!("{format}:AD"));
    for sample in &mut record.samples {
        sample.push(':');
        sample.push_str(&missing_ad);
    }
}

fn ensure_missing_dp(record: &mut vcf::RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    if format.split(':').any(|field| field == "DP") {
        return;
    }
    record.format = Some(format!("{format}:DP"));
    for sample in &mut record.samples {
        sample.push_str(":0");
    }
}

fn remap_gt(gt: &str, mapping: &[usize]) -> String {
    let separator = if gt.contains('|') { "|" } else { "/" };
    gt.split(['/', '|'])
        .map(|token| {
            token.parse::<usize>().ok().map_or_else(
                || token.to_string(),
                |allele| mapping.get(allele).copied().unwrap_or(0).to_string(),
            )
        })
        .collect::<Vec<_>>()
        .join(separator)
}

/// VariantWriter emits unphased, canonical diploid genotypes. Ref/ALT calls
/// use ascending order (`1|0` → `0/1`), while het-of-ALT calls use the legacy
/// aggregator's later/earlier order (`1|2` → `2/1`).
fn canonicalize_legacy_genotypes(record: &mut vcf::RawVcfRecord) {
    let Some(gt_index) = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|field| field == "GT"))
    else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        let Some(gt) = cells.get_mut(gt_index) else {
            continue;
        };
        if !gt.contains('|') {
            continue;
        }
        let alleles = gt.split(['/', '|']).collect::<Vec<_>>();
        if alleles.len() == 2
            && let (Ok(left), Ok(right)) =
                (alleles[0].parse::<usize>(), alleles[1].parse::<usize>())
        {
            *gt = if left > 0 && right > 0 && left != right {
                format!("{}/{}", left.max(right), left.min(right))
            } else {
                format!("{}/{}", left.min(right), left.max(right))
            };
        } else {
            *gt = gt.replace('|', "/");
        }
        *sample = cells.join(":");
    }
}

fn somatic_info_sample_names(headers: &[String]) -> Vec<String> {
    let input_names: Vec<String> = headers
        .iter()
        .find(|line| line.starts_with("#CHROM"))
        .map(|line| line.split('\t').skip(9).map(str::to_string).collect())
        .unwrap_or_default();
    if input_names.len() <= 1 {
        vec!["SAMPLE".to_string()]
    } else {
        input_names
    }
}

/// Legacy `alleles` declares one INFO field per input FORMAT/sample pair.
/// A single input sample is deliberately renamed to the output sample name
/// `SAMPLE`; multi-sample inputs keep their original names as INFO prefixes.
fn append_somatic_info_headers(headers: &mut Vec<String>, sample_names: &[String]) {
    let format_headers: Vec<String> = headers
        .iter()
        .filter(|line| line.starts_with("##FORMAT=<ID="))
        .cloned()
        .collect();
    let mut additions = Vec::new();
    for line in format_headers {
        let Some(body) = line.strip_prefix("##FORMAT=<ID=") else {
            continue;
        };
        let id_end = body.find([',', '>']).unwrap_or(body.len());
        let id = &body[..id_end];
        let suffix = &body[id_end..];
        for sample in sample_names {
            additions.push(format!("##INFO=<ID={sample}_{id}{suffix}"));
        }
    }
    let insertion = headers
        .iter()
        .position(|line| line.starts_with("#CHROM"))
        .unwrap_or(headers.len());
    headers.splice(insertion..insertion, additions);
}

fn append_info_value(info: &mut String, key: &str, value: &str) {
    if value.is_empty() || value == "." {
        return;
    }
    if info.is_empty() || info == "." {
        *info = format!("{key}={value}");
    } else {
        info.push(';');
        info.push_str(key);
        info.push('=');
        info.push_str(value);
    }
}

fn bcf_encoded_gt(gt: &str) -> String {
    let mut phased_next = false;
    let mut values = Vec::new();
    let mut token = String::new();
    let flush = |token: &mut String, phased: bool, values: &mut Vec<String>| {
        if token.is_empty() {
            return;
        }
        let encoded = if token == "." {
            0
        } else {
            token
                .parse::<i32>()
                .map(|allele| ((allele + 1) << 1) | i32::from(phased))
                .unwrap_or(0)
        };
        values.push(encoded.to_string());
        token.clear();
    };
    for character in gt.chars() {
        match character {
            '/' | '|' => {
                flush(&mut token, phased_next, &mut values);
                phased_next = character == '|';
            }
            _ => token.push(character),
        }
    }
    flush(&mut token, phased_next, &mut values);
    values.join(",")
}

fn copy_sample_formats_to_info(record: &mut vcf::RawVcfRecord, sample_names: &[String]) {
    let keys: Vec<String> = record
        .format
        .as_deref()
        .map(|format| format.split(':').map(str::to_string).collect())
        .unwrap_or_default();
    let samples = record.samples.clone();
    // `alleles` declares fields for every input sample but, due to its use of
    // the translated output header while mutating the record, only the first
    // sample's values survive into the emitted INFO payload. Preserve this
    // long-standing observable quirk for byte parity.
    for (sample_index, sample) in samples.iter().take(1).enumerate() {
        let prefix = sample_names
            .get(sample_index)
            .or_else(|| sample_names.first())
            .map(String::as_str)
            .unwrap_or("SAMPLE");
        for (key, value) in keys.iter().zip(sample.split(':')) {
            let encoded;
            let value = if key == "GT" {
                encoded = bcf_encoded_gt(value);
                encoded.as_str()
            } else {
                value
            };
            append_info_value(&mut record.info, &format!("{prefix}_{key}"), value);
        }
    }
}

fn first_sample_gt(record: &vcf::RawVcfRecord) -> String {
    let gt_index = record
        .format
        .as_deref()
        .and_then(|format| format.split(':').position(|key| key == "GT"));
    gt_index
        .and_then(|index| record.samples.first()?.split(':').nth(index))
        .unwrap_or(".")
        .to_string()
}

fn converted_somatic_gt(mode: SomaticGtMode, has_alt: bool) -> String {
    let allele = if has_alt { "1" } else { "0" };
    match mode {
        SomaticGtMode::Half => format!("./{allele}"),
        SomaticGtMode::Hemi => allele.to_string(),
        SomaticGtMode::Het => format!("0/{allele}"),
        SomaticGtMode::Hom => format!("{allele}/{allele}"),
        SomaticGtMode::First => unreachable!("first mode preserves the source GT"),
    }
}

/// Port the record-level behavior of legacy's `alleles` helper: move every
/// input FORMAT into sample-prefixed INFO, collapse to one output sample, and
/// split each ALT into its own record unless `first` was selected.
fn convert_somatic_record(
    record: &vcf::RawVcfRecord,
    mode: SomaticGtMode,
    sample_names: &[String],
) -> Vec<vcf::RawVcfRecord> {
    let first_gt = first_sample_gt(record);
    let mut base = record.clone();
    copy_sample_formats_to_info(&mut base, sample_names);
    base.format = Some("GT:AD:DP".to_string());

    let with_missing_depths = |gt: String, alt_count: usize| {
        let ad = vec!["."; alt_count.saturating_add(1)].join(",");
        format!("{gt}:{ad}:0")
    };

    if mode == SomaticGtMode::First {
        let alt_count = base
            .alt_allele
            .split(',')
            .filter(|alt| !alt.is_empty() && *alt != ".")
            .count();
        base.samples = vec![with_missing_depths(first_gt, alt_count)];
        return vec![base];
    }

    let alts: Vec<String> = base
        .alt_allele
        .split(',')
        .filter(|alt| !alt.is_empty() && *alt != ".")
        .map(str::to_string)
        .collect();
    if alts.is_empty() {
        base.samples = vec![with_missing_depths(converted_somatic_gt(mode, false), 0)];
        return vec![base];
    }

    alts.into_iter()
        .map(|alt| {
            let mut split = base.clone();
            split.alt_allele = alt;
            split.samples = vec![with_missing_depths(converted_somatic_gt(mode, true), 1)];
            split
        })
        .collect()
}

/// Reproduce the subsequent partial-credit aggregation that consumes the
/// biallelic records emitted by `alleles`. Non-hom modes merge sibling ALTs
/// back into one location with the legacy reversed het-alt GT; hom-alt calls
/// remain separate records.
fn finalize_somatic_records(
    mut records: Vec<vcf::RawVcfRecord>,
    mode: SomaticGtMode,
) -> Vec<vcf::RawVcfRecord> {
    if records.is_empty() {
        return records;
    }
    if mode == SomaticGtMode::First {
        records.retain_mut(trim_uncalled_non_ref);
        return records;
    }
    records.retain(|record| {
        !record.alt_allele.starts_with('<') && record.alt_allele != "*" && record.alt_allele != "."
    });
    if records.is_empty() {
        return records;
    }
    if mode == SomaticGtMode::Hom {
        records.sort_by(|left, right| {
            left.alt_allele
                .len()
                .cmp(&right.alt_allele.len())
                .then(left.alt_allele.cmp(&right.alt_allele))
        });
        return records;
    }
    if records.len() == 1 {
        records[0].samples = vec!["0/1:.,.:0".to_string()];
        return records;
    }

    let mut merged = records.remove(0);
    let mut alts = vec![merged.alt_allele.clone()];
    alts.extend(records.into_iter().map(|record| record.alt_allele));
    alts.sort_by(|left, right| left.len().cmp(&right.len()).then(left.cmp(right)));
    alts.dedup();
    merged.alt_allele = alts.join(",");
    let ad = vec!["."; alts.len() + 1].join(",");
    merged.samples = vec![format!("2/1:{ad}:0")];
    vec![merged]
}

fn finalize_somatic_for_pipeline(
    mut converted: Vec<vcf::RawVcfRecord>,
    mode: SomaticGtMode,
    normalization_enabled: bool,
) -> Vec<vcf::RawVcfRecord> {
    if normalization_enabled {
        finalize_somatic_records(converted, mode)
    } else {
        // scmp-somatic disables the partial-credit normalizer. In that lane
        // the `alleles` helper's half-call (`./1`) and per-ALT record grain
        // pass straight through to scmp; sibling aggregation belongs to the
        // normalizer and must not silently turn the call back into `0/1`.
        for record in &mut converted {
            let gt = first_sample_gt(record);
            record.format = Some("GT".to_string());
            record.samples = vec![gt];
        }
        converted
    }
}

fn rewrite_header_for_single_sample(headers: &mut [String]) {
    if let Some(line) = headers.iter_mut().find(|line| line.starts_with("#CHROM")) {
        *line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE".to_string();
    }
}

fn validate_record_reference(
    record: &vcf::RawVcfRecord,
    reference_sequences: &std::collections::BTreeMap<String, String>,
) -> Result<()> {
    let reference = reference_sequences
        .get(&record.chrom)
        .ok_or_else(|| anyhow::anyhow!("reference contig {} not found", record.chrom))?;
    // DNA references are ASCII, so byte slicing is safe and avoids allocating a
    // 48 M-entry Vec<char> per record on chr-scale references.
    let bases = reference.as_bytes();
    let end = record.end_pos();
    if record.pos == 0 || end > bases.len() {
        bail!(
            "record {}:{} extends beyond the reference sequence",
            record.chrom,
            record.pos
        );
    }
    let observed = &bases[record.pos - 1..end];
    // Reference FASTA may contain soft-masked lowercase bases (repeat regions).
    // Legacy hap.py tolerates case differences between the VCF REF and the
    // reference sequence — any N in either side also matches anything.
    if !ref_bytes_equal(observed, record.ref_allele.as_bytes()) {
        bail!(
            "record {}:{} REF={} does not match reference {}",
            record.chrom,
            record.pos,
            record.ref_allele,
            String::from_utf8_lossy(observed)
        );
    }
    Ok(())
}

/// INFO keys whose semantics assume a specific allele layout. They go stale as
/// soon as preprocessing splits or decomposes variants, so the legacy C++
/// `preprocess` binary (via `VariantCallsOnly`) drops them on output. Keep this
/// list tight — any tag we strip permanently is one a consumer can't rely on.
const STALE_INFO_KEYS: &[&str] = &["AC", "AN", "MLEAC", "MLEAF"];

/// Reorder `record.info` so its tags appear alphabetically, matching legacy
/// `std::map`-ordered serialisation. Flag-only entries (no `=`) keep their
/// position in the sort by bare tag name. Numeric values that round-trip
/// through htslib's float/int parsing are also collapsed to canonical form
/// (e.g. `-0` → `0`, matching legacy's `bcf_update_info_float` write path
/// which loses the sign on negative zero).
fn sort_info_keys(record: &mut vcf::RawVcfRecord) {
    if record.info.is_empty() || record.info == "." {
        return;
    }
    let mut entries: Vec<String> = record
        .info
        .split(';')
        .map(canonicalise_info_entry)
        .collect();
    entries.sort_by_key(|entry| {
        let key = entry.split('=').next().unwrap_or(entry);
        key.to_string()
    });
    record.info = entries.join(";");
}

/// Collapse `-0`, `-0.0`, etc. in INFO values to `0`. Legacy hap.py reads
/// every INFO value through htslib's typed parsers (`bcf_update_info_float`
/// / `bcf_update_info_int32`) which lose the sign on negative zero before
/// re-emitting the record. Mirroring that here on the main code path keeps
/// our output byte-identical without a diff-only shim.
fn canonicalise_info_entry(entry: &str) -> String {
    let Some((key, value)) = entry.split_once('=') else {
        return entry.to_string();
    };
    let normalised: Vec<String> = value
        .split(',')
        .map(|part| {
            // Only collapse a value that *parses* as a number whose float
            // representation is exactly zero. Leaves non-numeric strings
            // (e.g. `set=variant2`, `culprit=FS`) untouched.
            if let Ok(f) = part.parse::<f64>()
                && f == 0.0
                && (part.starts_with('-') || part.contains('-'))
            {
                // Preserve the original integer / float visual shape
                // (e.g. `-0` → `0`, `-0.0` → `0`) so downstream byte
                // comparison sees what htslib would emit.
                return if part.contains('.') || part.contains('e') || part.contains('E') {
                    // Keep float shape — bcftools emits "0" for any
                    // signed-zero float regardless of original
                    // precision; mirror that.
                    "0".to_string()
                } else {
                    "0".to_string()
                };
            }
            part.to_string()
        })
        .collect();
    format!("{key}={}", normalised.join(","))
}

/// Reduce every sample's PL cell to its last comma-separated value. Legacy
/// hap.py's `preprocess` stores PL as a scalar int per sample (see
/// `VariantWriter.cpp:563`) and htslib serialises that as the final value only.
fn collapse_pl_to_last_value(record: &mut vcf::RawVcfRecord) {
    let Some(format) = &record.format else {
        return;
    };
    let Some(pl_index) = format.split(':').position(|field| field == "PL") else {
        return;
    };
    for sample in &mut record.samples {
        let mut fields: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        if let Some(cell) = fields.get_mut(pl_index)
            && let Some(last) = cell.rsplit(',').next()
        {
            *cell = last.to_string();
        }
        *sample = fields.join(":");
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum ScalarType {
    Integer = 0,
    Float = 1,
    String = 2,
}

/// Classify a per-sample value as integer, float, or string based on the
/// first non-missing scalar. Mirrors `VariantWriter.cpp`'s runtime dispatch
/// (`v.asInt()` / `v.isNumeric()` / `v.asString()`).
fn classify_value(value: &str) -> ScalarType {
    if value.is_empty() || value == "." {
        return ScalarType::Integer;
    }
    // Multi-valued entries (`3279,0,716`, `32,126`) classify on the first
    // numeric element — mirrors legacy's int/float scalar decision.
    let first = value.split(',').next().unwrap_or(value);
    if first == "." || first.is_empty() {
        return ScalarType::Integer;
    }
    if first.parse::<i64>().is_ok() {
        ScalarType::Integer
    } else if first.parse::<f64>().is_ok() {
        ScalarType::Float
    } else {
        ScalarType::String
    }
}

/// Reorder FORMAT fields and sample columns so that legacy's type-bucketed
/// canonical order is used: GT, AD, ADO, DP come first (in that order);
/// everything else is partitioned into integer / float / string buckets by
/// runtime-inferred type and each bucket is sorted alphabetically.
fn reorder_format_fields(record: &mut vcf::RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    let field_names: Vec<String> = format.split(':').map(|s| s.to_string()).collect();

    // Parse each sample into its fields and pad to FORMAT width with ".".
    let original_samples: Vec<Vec<String>> = record
        .samples
        .iter()
        .map(|s| {
            let mut cells: Vec<String> = s.split(':').map(|x| x.to_string()).collect();
            while cells.len() < field_names.len() {
                cells.push(".".to_string());
            }
            cells
        })
        .collect();

    const FIXED: &[&str] = &["GT", "AD", "ADO", "DP"];

    let (fixed_fields, rest_fields): (Vec<_>, Vec<_>) = field_names
        .iter()
        .enumerate()
        .partition(|(_, name)| FIXED.contains(&name.as_str()));

    // Preserve GT/AD/ADO/DP order exactly as listed in FIXED.
    let mut ordered: Vec<(usize, String)> = FIXED
        .iter()
        .filter_map(|target| {
            fixed_fields
                .iter()
                .find(|(_, name)| name.as_str() == *target)
                .map(|(idx, name)| (*idx, name.to_string()))
        })
        .collect();

    // Type-bucket the remaining fields. Use the first sample's value to pick.
    let probe_sample = original_samples.first();
    let mut tagged: Vec<(ScalarType, String, usize)> = rest_fields
        .iter()
        .map(|(idx, name)| {
            let value = probe_sample
                .and_then(|s| s.get(*idx))
                .map(String::as_str)
                .unwrap_or(".");
            (classify_value(value), name.to_string(), *idx)
        })
        .collect();
    tagged.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    ordered.extend(tagged.into_iter().map(|(_, name, idx)| (idx, name)));

    let new_format = ordered
        .iter()
        .map(|(_, name)| name.clone())
        .collect::<Vec<_>>()
        .join(":");
    let new_samples: Vec<String> = original_samples
        .iter()
        .map(|sample| {
            ordered
                .iter()
                .map(|(original_idx, _)| {
                    sample
                        .get(*original_idx)
                        .cloned()
                        .unwrap_or_else(|| ".".to_string())
                })
                .collect::<Vec<_>>()
                .join(":")
        })
        .collect();

    record.format = Some(new_format);
    record.samples = new_samples;
}

fn canonicalize_multi_allelic_order(record: &mut vcf::RawVcfRecord) {
    let alts: Vec<String> = record.alt_allele.split(',').map(str::to_string).collect();
    if alts.len() < 2 || alts.iter().any(|alt| alt.starts_with('<') || alt == "*") {
        return;
    }
    let mut ordered: Vec<(usize, String)> = alts.iter().cloned().enumerate().collect();
    ordered.sort_by(|left, right| left.1.len().cmp(&right.1.len()).then(left.1.cmp(&right.1)));
    let mut mapping = vec![0usize; alts.len() + 1];
    for (new_offset, (old_offset, _)) in ordered.iter().enumerate() {
        mapping[old_offset + 1] = new_offset + 1;
    }
    let keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect())
        .unwrap_or_default();
    let gt_index = keys.iter().position(|key| *key == "GT");
    let ad_index = keys.iter().position(|key| *key == "AD");
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        if let Some(index) = gt_index
            && let Some(gt) = cells.get_mut(index)
        {
            let remapped = remap_gt(gt, &mapping);
            if remapped.contains('/') {
                let mut alleles: Vec<&str> = remapped.split('/').collect();
                if alleles.len() == 2
                    && alleles[0] != alleles[1]
                    && alleles
                        .iter()
                        .all(|allele| *allele != "0" && *allele != ".")
                {
                    alleles.sort_by(|left, right| right.cmp(left));
                    *gt = alleles.join("/");
                } else {
                    *gt = remapped;
                }
            } else {
                *gt = remapped;
            }
        }
        if let Some(index) = ad_index
            && let Some(ad) = cells.get_mut(index)
        {
            let depths: Vec<&str> = ad.split(',').collect();
            if depths.len() == alts.len() + 1 {
                let mut reordered = vec![depths[0]];
                reordered.extend(ordered.iter().map(|(old_offset, _)| depths[old_offset + 1]));
                *ad = reordered.join(",");
            }
        }
        *sample = cells.join(":");
    }
    record.alt_allele = ordered
        .into_iter()
        .map(|(_, alt)| alt)
        .collect::<Vec<_>>()
        .join(",");
}

fn string_format_fields(headers: &[String]) -> BTreeSet<String> {
    headers
        .iter()
        .filter(|header| header.starts_with("##FORMAT=<") && header.contains("Type=String"))
        .filter_map(|header| {
            header
                .strip_prefix("##FORMAT=<ID=")
                .and_then(|tail| tail.split_once(',').map(|(id, _)| id.to_string()))
        })
        .collect()
}

fn blank_secondary_sample_annotations(
    record: &mut vcf::RawVcfRecord,
    bcf_output: bool,
    string_fields: &BTreeSet<String>,
) {
    let keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect())
        .unwrap_or_default();
    for sample in record.samples.iter_mut().skip(1) {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        let homref = keys
            .iter()
            .position(|key| *key == "GT")
            .and_then(|index| cells.get(index))
            .is_some_and(|gt| gt == "0/0" || gt == "0|0");
        let missing_gt = keys
            .iter()
            .position(|key| *key == "GT")
            .and_then(|index| cells.get(index))
            .is_some_and(|gt| {
                gt.chars()
                    .all(|character| matches!(character, '.' | '/' | '|'))
            });
        for (index, key) in keys.iter().enumerate() {
            if missing_gt
                && *key != "GT"
                && let Some(value) = cells.get_mut(index)
            {
                *value = if *key == "AD" {
                    ".,.".to_string()
                } else {
                    ".".to_string()
                };
            } else if !matches!(*key, "GT" | "AD" | "ADO" | "DP" | "PL")
                && let Some(value) = cells.get_mut(index)
            {
                *value = if bcf_output && string_fields.contains(*key) {
                    String::new()
                } else {
                    ".".to_string()
                };
            } else if homref
                && *key == "AD"
                && let Some(value) = cells.get_mut(index)
            {
                let width = value.split(',').count().max(1);
                *value = vec!["."; width].join(",");
            } else if homref
                && *key == "ADO"
                && let Some(value) = cells.get_mut(index)
            {
                *value = ".".to_string();
            }
        }
        *sample = cells.join(":");
    }
}

/// Header records installed by the legacy C++ `VariantWriter` constructor.
/// Definitions here win over conflicting input definitions during merging.
const LEGACY_BASE_HEADERS: &[&str] = &[
    "##fileformat=VCFv4.1",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##reference=hg19",
    "##contig=<ID=chr1,length=249250621>",
    "##contig=<ID=chr2,length=243199373>",
    "##contig=<ID=chr3,length=198022430>",
    "##contig=<ID=chr4,length=191154276>",
    "##contig=<ID=chr5,length=180915260>",
    "##contig=<ID=chr6,length=171115067>",
    "##contig=<ID=chr7,length=159138663>",
    "##contig=<ID=chr8,length=146364022>",
    "##contig=<ID=chr9,length=141213431>",
    "##contig=<ID=chr10,length=135534747>",
    "##contig=<ID=chr11,length=135006516>",
    "##contig=<ID=chr12,length=133851895>",
    "##contig=<ID=chr13,length=115169878>",
    "##contig=<ID=chr14,length=107349540>",
    "##contig=<ID=chr15,length=102531392>",
    "##contig=<ID=chr16,length=90354753>",
    "##contig=<ID=chr17,length=81195210>",
    "##contig=<ID=chr18,length=78077248>",
    "##contig=<ID=chr19,length=59128983>",
    "##contig=<ID=chr20,length=63025520>",
    "##contig=<ID=chr21,length=48129895>",
    "##contig=<ID=chr22,length=51304566>",
    "##contig=<ID=chrX,length=155270560>",
    "##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">",
    "##INFO=<ID=IMPORT_FAIL,Number=.,Type=Flag,Description=\"Flag to identify variants that could not be imported.\">",
    "##FORMAT=<ID=AGT,Number=1,Type=String,Description=\"Genotypes at ambiguous locations\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
    "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depths\">",
    "##FORMAT=<ID=ADO,Number=.,Type=Integer,Description=\"Summed depth of non-called alleles.\">",
];

/// Reproduce `VariantWriterImpl::writeHeader`: keep constructor-installed
/// records first, sort supplied meta-lines lexicographically, and retain only
/// the first structured declaration for each `(header key, ID)` pair.
pub(crate) fn canonicalize_legacy_headers(headers: &[String]) -> Vec<String> {
    let mut output: Vec<String> = LEGACY_BASE_HEADERS
        .iter()
        .map(|line| (*line).to_string())
        .collect();
    let mut identities: BTreeSet<String> = output
        .iter()
        .filter_map(|line| structured_header_identity(line))
        .collect();
    let mut exact: BTreeSet<String> = output.iter().cloned().collect();

    let mut supplied: Vec<&String> = headers
        .iter()
        .filter(|line| line.starts_with("##") && !line.starts_with("##fileformat="))
        .collect();
    supplied.sort_unstable();
    for line in supplied {
        if let Some(identity) = structured_header_identity(line) {
            if identities.insert(identity) {
                exact.insert(line.clone());
                output.push(line.clone());
            }
        } else if exact.insert(line.clone()) {
            output.push(line.clone());
        }
    }

    if let Some(chrom) = headers.iter().rev().find(|line| line.starts_with("#CHROM")) {
        output.push(chrom.clone());
    }
    output
}

/// Return htslib's merge identity for `<ID=...>` header records.
pub(crate) fn structured_header_identity(line: &str) -> Option<String> {
    let body = line.strip_prefix("##")?;
    let (key, value) = body.split_once('=')?;
    let id_value = value.strip_prefix("<ID=")?;
    let id = id_value.split([',', '>']).next()?;
    Some(format!("{key}:{id}"))
}

/// Insert `ADO` into a record's FORMAT immediately after `AD`. For each sample
/// compute ADO = sum of AD values for allele indices that are NOT present in
/// the sample's GT call (i.e. the depth of "other" alleles that were not the
/// genotype's choice). Matches `VariantWriter.cpp:603` semantics.
fn insert_ado_format(record: &mut vcf::RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    let fields: Vec<&str> = format.split(':').collect();
    if fields.contains(&"ADO") {
        return;
    }
    let Some(ad_index) = fields.iter().position(|field| *field == "AD") else {
        return;
    };
    let gt_index = fields.iter().position(|field| *field == "GT");

    let mut new_format_fields: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
    new_format_fields.insert(ad_index + 1, "ADO".to_string());
    record.format = Some(new_format_fields.join(":"));

    for sample in &mut record.samples {
        let mut sample_fields: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        // Pad with missing cells if the sample is shorter than FORMAT (legal
        // VCF trailing-missing shorthand); we only need up to `ad_index + 1`.
        while sample_fields.len() <= ad_index {
            sample_fields.push(".".to_string());
        }
        let ado = compute_ado(
            sample_fields
                .get(ad_index)
                .map(String::as_str)
                .unwrap_or("."),
            gt_index.and_then(|i| sample_fields.get(i).map(String::as_str)),
        );
        sample_fields.insert(ad_index + 1, ado.to_string());
        *sample = sample_fields.join(":");
    }
}

/// Compute `ADO` = sum of allele depths for allele indices that don't appear
/// in the genotype call. Returns 0 when either AD or GT is missing, matching
/// legacy's `bcf_int32_missing` → 0 rendering for unknown.
fn compute_ado(ad_cell: &str, gt_cell: Option<&str>) -> i64 {
    if ad_cell == "." || ad_cell.is_empty() {
        return 0;
    }
    let ad_values: Vec<i64> = ad_cell
        .split(',')
        .map(|v| v.parse::<i64>().unwrap_or(0))
        .collect();
    if ad_values.is_empty() {
        return 0;
    }
    let gt = match gt_cell {
        Some(gt) if !gt.is_empty() && gt != "." => gt,
        _ => return 0,
    };
    let called: std::collections::BTreeSet<usize> = gt
        .split(['/', '|'])
        .filter_map(|a| a.parse::<usize>().ok())
        .collect();
    if called.is_empty() {
        return 0;
    }
    ad_values
        .iter()
        .enumerate()
        .filter_map(|(index, depth)| (!called.contains(&index)).then_some(*depth))
        .sum()
}

/// Run partial-credit left-shift + trim on a single-ALT record, mutating
/// `record.pos` / `record.ref_allele` / `record.alt_allele` in place.
/// Multi-allelic records are left untouched by the caller — they need the
/// split/decomposition pass first.
fn apply_left_shift(record: &mut vcf::RawVcfRecord, reference: &[u8], neighbor_end: usize) {
    let ref_len = record.ref_allele.len();
    if ref_len == 0 {
        return;
    }
    let end = record.pos + ref_len - 1;
    let mut rv = partial_credit::RefVar {
        start: record.pos,
        end,
        alt: record.alt_allele.clone(),
    };
    // neighbor_end is the reference end of the previous record: don't allow
    // left-shifting into that span (mirrors legacy partialcredit.py behaviour).
    let pos_min = record
        .pos
        .saturating_sub(LEFT_SHIFT_WINDOW)
        .max(1)
        .max(neighbor_end);
    partial_credit::left_shift(reference, &mut rv, pos_min, true);

    // Rebuild REF from the reference slice now that start/end may have moved.
    let new_ref_len = rv.end as i64 - rv.start as i64 + 1;
    if new_ref_len < 0 || rv.start == 0 {
        return;
    }
    let new_end_usize = rv.start + (new_ref_len.max(0) as usize).saturating_sub(0);
    // Guard against slicing past the reference (noodles reads give us ASCII bytes).
    if rv.start == 0 || new_end_usize > reference.len() + 1 {
        return;
    }
    let ref_bytes = if new_ref_len <= 0 {
        &[][..]
    } else {
        &reference[rv.start - 1..rv.end]
    };
    // Uppercase soft-masked repeat bases so the rebuilt REF matches legacy's
    // canonical output.
    let new_ref = String::from_utf8_lossy(ref_bytes).to_ascii_uppercase();

    // Don't emit zero-length REF or ALT — VCF requires at least one base on
    // each side (partial-credit keeps a left-anchor via `ref_padding=true`, so
    // this is belt-and-braces against edge cases).
    if new_ref.is_empty() || rv.alt.is_empty() {
        return;
    }

    record.pos = rv.start;
    record.ref_allele = new_ref;
    record.alt_allele = rv.alt;
}

fn record_reference_matches(record: &vcf::RawVcfRecord, reference: &[u8]) -> bool {
    let start = record.pos.saturating_sub(1);
    let end = start.saturating_add(record.ref_allele.len());
    reference
        .get(start..end)
        .is_some_and(|observed| observed.eq_ignore_ascii_case(record.ref_allele.as_bytes()))
}

/// Equivalent of `bcftools norm -f REF -c x -D` for one non-symbolic record.
/// Alleles are normalized independently, then padded to one common site so a
/// multi-allelic record retains its original allele indexes.
fn normalize_bcftools_record(record: &mut vcf::RawVcfRecord, reference: &[u8]) {
    let alts: Vec<&str> = record.alt_allele.split(',').collect();
    if alts.iter().any(|alt| {
        alt.is_empty()
            || *alt == "."
            || alt.starts_with('<')
            || *alt == "*"
            || alt.contains(['[', ']'])
    }) {
        return;
    }
    let end = record.end_pos();
    let mut normalized = Vec::with_capacity(alts.len());
    for alt in alts {
        let mut variant = partial_credit::RefVar {
            start: record.pos,
            end,
            alt: alt.to_string(),
        };
        partial_credit::left_shift(reference, &mut variant, 1, false);
        normalized.push(pad_bcftools_allele(variant, reference));
    }
    let common_start = normalized
        .iter()
        .map(|variant| variant.start)
        .min()
        .unwrap_or(record.pos);
    let common_end = normalized
        .iter()
        .map(|variant| variant.end)
        .max()
        .unwrap_or(end);
    if common_start == 0 || common_end < common_start || common_end > reference.len() {
        return;
    }
    let alts = normalized
        .into_iter()
        .map(|variant| {
            let mut allele = Vec::new();
            allele.extend_from_slice(&reference[common_start - 1..variant.start - 1]);
            allele.extend_from_slice(variant.alt.as_bytes());
            allele.extend_from_slice(&reference[variant.end..common_end]);
            String::from_utf8_lossy(&allele).to_ascii_uppercase()
        })
        .collect::<Vec<_>>();
    record.pos = common_start;
    record.ref_allele =
        String::from_utf8_lossy(&reference[common_start - 1..common_end]).to_ascii_uppercase();
    record.alt_allele = alts.join(",");
}

fn uppercase_alleles_preserving_breakends(alts: &str) -> String {
    alts.split(',')
        .map(|alt| {
            let Some(first) = alt.find(['[', ']']) else {
                return alt.to_ascii_uppercase();
            };
            let Some(relative_second) = alt[first + 1..].find(['[', ']']) else {
                return alt.to_ascii_uppercase();
            };
            let second = first + 1 + relative_second;
            format!(
                "{}{}{}",
                alt[..=first].to_ascii_uppercase(),
                &alt[first + 1..second],
                alt[second..].to_ascii_uppercase()
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn materialize_unsupported_import_failure(record: &mut vcf::RawVcfRecord) -> bool {
    let breakend = record.alt_allele.contains(['[', ']']);
    let unsupported = record.alt_allele.split(',').any(|alt| {
        alt.contains(['[', ']']) || (alt.starts_with('<') && !matches!(alt, "<DEL>" | "<NON_REF>"))
    });
    if !unsupported {
        return false;
    }

    record.alt_allele = ".".to_string();
    let mut info = record
        .info
        .split(';')
        .filter(|field| !field.is_empty() && *field != ".")
        .filter(|field| {
            let key = field.split_once('=').map_or(*field, |(key, _)| key);
            (key != "END" || !breakend) && key != "IMPORT_FAIL"
        })
        .map(str::to_string)
        .collect::<Vec<_>>();
    // Breakends have no local span, whereas unsupported symbolic records keep
    // their declared END when the legacy reader materializes IMPORT_FAIL.
    if breakend || !info.iter().any(|field| field.starts_with("END=")) {
        info.push(format!("END={}", record.pos));
    }
    info.push("IMPORT_FAIL".to_string());
    record.info = info.join(";");

    let Some(gt_index) = record.format_keys().iter().position(|key| *key == "GT") else {
        return true;
    };
    let ad_index = record.format_keys().iter().position(|key| *key == "AD");
    for sample in &mut record.samples {
        let mut fields = sample.split(':').map(str::to_string).collect::<Vec<_>>();
        if let Some(gt) = fields.get_mut(gt_index) {
            *gt = "0/0".to_string();
        }
        // The failed BND import has no alternate allele. VariantWriter
        // consequently writes just the reference AD value; retaining the old
        // alternate depth would also make the generated ADO non-zero.
        if let Some(ad_index) = ad_index
            && let Some(ad) = fields.get_mut(ad_index)
        {
            *ad = ad.split(',').next().unwrap_or(".").to_string();
        }
        *sample = fields.join(":");
    }
    true
}

fn pad_bcftools_allele(
    mut variant: partial_credit::RefVar,
    reference: &[u8],
) -> partial_credit::RefVar {
    let reference_len = variant.end as i64 - variant.start as i64 + 1;
    if reference_len <= 0 && !variant.alt.is_empty() {
        if variant.start > 1 {
            let anchor = variant.start - 1;
            variant.start = anchor;
            variant.end = anchor;
            variant
                .alt
                .insert(0, reference[anchor - 1].to_ascii_uppercase() as char);
        }
    } else if reference_len > 0 && variant.alt.is_empty() {
        if variant.start > 1 {
            let anchor = variant.start - 1;
            variant.start = anchor;
            variant
                .alt
                .push(reference[anchor - 1].to_ascii_uppercase() as char);
        } else if variant.end < reference.len() {
            variant.end += 1;
            variant
                .alt
                .push(reference[variant.end - 1].to_ascii_uppercase() as char);
        }
    }
    variant
}

/// Remove stale allele-count INFO entries from `record.info` in place.
fn strip_stale_info_keys(record: &mut vcf::RawVcfRecord) {
    if record.info.is_empty() || record.info == "." {
        return;
    }
    let kept: Vec<&str> = record
        .info
        .split(';')
        .filter(|entry| {
            let key = entry.split('=').next().unwrap_or(entry);
            !STALE_INFO_KEYS.contains(&key)
        })
        .collect();
    record.info = if kept.is_empty() {
        ".".to_string()
    } else {
        kept.join(";")
    };
}

/// Mirror legacy `VariantAlleleSplitter`'s haploid → diploid treatment:
///
/// * `0`   → `0/0`
/// * autosomal `<n>` (n > 0) → `0/<n>` (legacy creates a het half-call)
/// * sex-chromosome `<n>` (n > 0) → `<n>/<n>`
/// * `.`   → `./.`
/// * `./<n>` → `0/<n>` (legacy completes a diploid half-call with REF)
///
/// Fully called diploid inputs pass through unchanged.
fn normalise_haploid_genotypes(record: &mut vcf::RawVcfRecord, male: bool) {
    if male {
        expand_male_sex_chromosome_genotypes(record);
    }
    let Some(format) = &record.format else {
        return;
    };
    let Some(gt_index) = format.split(':').position(|f| f == "GT") else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        if let Some(cell) = cells.get_mut(gt_index) {
            *cell = expand_haploid_gt(cell, false);
        }
        *sample = cells.join(":");
    }
}

fn expand_male_sex_chromosome_genotypes(record: &mut vcf::RawVcfRecord) {
    if !matches!(record.chrom.as_str(), "X" | "Y" | "chrX" | "chrY") {
        return;
    }
    let Some(format) = &record.format else {
        return;
    };
    let Some(gt_index) = format.split(':').position(|field| field == "GT") else {
        return;
    };
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        if let Some(gt) = cells.get_mut(gt_index) {
            *gt = expand_haploid_gt(gt, true);
        }
        *sample = cells.join(":");
    }
}

/// Expand a haploid GT token to its diploid legacy form.
fn expand_haploid_gt(gt: &str, sex_chromosome: bool) -> String {
    if let Some(separator) = gt.chars().find(|separator| matches!(separator, '/' | '|')) {
        let alleles = gt.split(separator).collect::<Vec<_>>();
        if alleles.len() == 2
            && alleles[0] == "."
            && alleles[1].parse::<u32>().is_ok_and(|allele| allele > 0)
        {
            return format!("0{separator}{}", alleles[1]);
        }
        return gt.to_string();
    }
    match gt {
        "." => "./.".to_string(),
        "0" => "0/0".to_string(),
        // A haploid alt call (single non-zero allele index) expands to a
        // heterozygous half-call — the C++ VariantAlleleSplitter treats
        // ngt==1 with gt[0]>0 as a het half-call and emits `0/N` after the
        // merge (VariantAlleleSplitter.cpp:180-227). Legacy result.vcf.gz
        // output confirms this: query GT=1 on autosomes appears as `0/1:het`
        // in the comparison output, not `1/1:homalt`.
        other => match other.parse::<u32>() {
            Ok(n) if n > 0 && sex_chromosome => format!("{n}/{n}"),
            Ok(n) if n > 0 => format!("0/{n}"),
            _ => gt.to_string(),
        },
    }
}

/// Case-insensitive DNA byte-slice equality with IUPAC `N` treated as a wildcard.
/// Mirrors the tolerance of the legacy hap.py / bcftools preprocessing pass.
fn ref_bytes_equal(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter().zip(b.iter()).all(|(x, y)| {
        let xu = x.to_ascii_uppercase();
        let yu = y.to_ascii_uppercase();
        xu == yu || xu == b'N' || yu == b'N'
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn filters_only_removes_records_carrying_only_selected_filters() {
        assert!(passes_filters_only("PASS", Some("LowQual,q10")));
        assert!(passes_filters_only(".", Some("LowQual,q10")));
        assert!(!passes_filters_only("LowQual", Some("LowQual,q10")));
        assert!(!passes_filters_only("LowQual;q10", Some("LowQual,q10")));
        assert!(passes_filters_only("LowQual;s50", Some("LowQual,q10")));
        assert!(passes_filters_only("LowQual", None));
    }

    #[test]
    fn gender_auto_matches_vcfcheck_haploid_x_heuristic() {
        let mut haploid = make_record(".");
        haploid.chrom = "chrX".to_string();
        haploid.format = Some("GT".to_string());
        haploid.samples = vec!["1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[haploid.clone()]),
            PreprocessGender::Male
        );

        let mut heterozygous = haploid.clone();
        heterozygous.samples = vec!["0/1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[haploid, heterozygous]),
            PreprocessGender::Female
        );
        assert_eq!(
            resolve_gender(PreprocessGender::None, &[]),
            PreprocessGender::None
        );
    }

    #[test]
    fn gender_auto_treats_half_called_x_genotype_as_diploid() {
        let mut half_called = make_record(".");
        half_called.chrom = "chrX".to_string();
        half_called.format = Some("GT".to_string());
        half_called.samples = vec!["./1".to_string()];

        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[half_called]),
            PreprocessGender::Female,
            "legacy vcfcheck uses ngt=2 and compares the missing and called slots"
        );
    }

    #[test]
    fn gender_auto_does_not_treat_lowercase_x_as_x_chromosome() {
        let mut lowercase_x = make_record(".");
        lowercase_x.chrom = "x".to_string();
        lowercase_x.format = Some("GT".to_string());
        lowercase_x.samples = vec!["1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[lowercase_x]),
            PreprocessGender::Female,
            "pinned vcfcheck compares its location variable to lowercase x"
        );

        let mut lowercase_chrx = make_record(".");
        lowercase_chrx.chrom = "chrx".to_string();
        lowercase_chrx.format = Some("GT".to_string());
        lowercase_chrx.samples = vec!["1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[lowercase_chrx]),
            PreprocessGender::Male
        );
    }

    #[test]
    fn region_and_auto_fixchr_complete_half_called_x_genotype_like_legacy() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let regions = directory.path().join("regions.bed");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=X,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "X\t2\tx_half_called\tA\tC\t20\tPASS\t.\tGT:AD\t./1:1,9\n",
            ),
        )?;
        fs::write(&reference, ">chrX\nAAAAA\n")?;
        fs::write(&regions, "chrX\t0\t5\n")?;

        let mut args = interval_args(&input, &output, &reference, Some(&regions), None);
        args.fixchr = None;
        args.gender = PreprocessGender::Auto;
        args.decompose = true;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].to_line(),
            "chrX\t2\t.\tA\tC\t20\t.\t.\tGT:AD:ADO:DP\t0/1:1,9:1:0"
        );
        Ok(())
    }

    #[test]
    fn male_gender_duplicates_haploid_x_and_y_only() {
        let mut x = make_record(".");
        x.chrom = "X".to_string();
        x.format = Some("GT:DP".to_string());
        x.samples = vec!["1:9".to_string()];
        expand_male_sex_chromosome_genotypes(&mut x);
        assert_eq!(x.samples, vec!["1/1:9"]);

        let mut autosome = x.clone();
        autosome.chrom = "1".to_string();
        autosome.samples = vec!["1:9".to_string()];
        expand_male_sex_chromosome_genotypes(&mut autosome);
        assert_eq!(autosome.samples, vec!["1:9"]);
    }

    #[test]
    fn reference_candidates_follow_legacy_hg19_then_hgref_precedence() -> Result<()> {
        let directory = tempdir()?;
        let explicit = directory.path().join("explicit.fa");
        let hg19 = directory.path().join("hg19.fa");
        let hgref = directory.path().join("hgref.fa");
        let fallback = directory.path().join("fallback.fa");
        for path in [&explicit, &hg19, &hgref, &fallback] {
            fs::write(path, ">chr1\nA\n")?;
        }
        assert_eq!(
            resolve_reference_candidates(Some(&explicit), Some(&hg19), Some(&hgref), &fallback)?,
            explicit
        );
        assert_eq!(
            resolve_reference_candidates(None, Some(&hg19), Some(&hgref), &fallback)?,
            hg19
        );
        fs::remove_file(&hg19)?;
        assert_eq!(
            resolve_reference_candidates(None, Some(&hg19), Some(&hgref), &fallback)?,
            hgref
        );
        Ok(())
    }

    #[test]
    fn legacy_fixchr_auto_only_adds_a_missing_prefix() {
        let prefixed = BTreeSet::from(["chr1".to_string(), "chrX".to_string()]);
        let plain = BTreeSet::from(["1".to_string(), "X".to_string()]);
        assert!(resolve_fixchr(None, &prefixed, &plain));
        assert!(!resolve_fixchr(None, &plain, &prefixed));
        assert!(!resolve_fixchr(None, &prefixed, &prefixed));

        assert_eq!(add_legacy_chr_prefix("1"), "chr1");
        assert_eq!(add_legacy_chr_prefix("MT"), "chrM");
        assert_eq!(add_legacy_chr_prefix("chrMT"), "chrM");
        assert_eq!(add_legacy_chr_prefix("GL000207.1"), "GL000207.1");
        assert_eq!(add_legacy_chr_prefix("chr1"), "chr1");
    }

    #[test]
    fn fixchr_adds_lengthless_headers_for_rewritten_contigs() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">1\nAAAAA\n>chr1\nAAAAA\n>chrX\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "1\t2\trs1\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.fixchr = None;
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        run(args)?;

        let (headers, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records[0].chrom, "chr1");
        assert!(
            headers
                .iter()
                .any(|line| line == "##contig=<ID=1,length=5>")
        );
        assert!(headers.iter().any(|line| line == "##contig=<ID=chr1>"));
        Ok(())
    }

    #[test]
    fn disabling_leftshift_and_decomposition_preserves_bcftools_view_shape() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FILTER=<ID=LowQual,Description=\"low\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\trs1\tA\tC\t10\tLowQual\tAC=1\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        run(args)?;
        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "rs1");
        assert_eq!(records[0].filter, "LowQual");
        assert_eq!(records[0].info, "AC=1");
        assert_eq!(records[0].format.as_deref(), Some("GT"));
        assert_eq!(records[0].samples, vec!["0/1"]);
        Ok(())
    }

    #[test]
    fn sites_only_vcf_is_rejected_like_vcfcheck() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("sites.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\n",
            ),
        )?;

        let error = run(interval_args(&input, &output, &reference, None, None)).unwrap_err();
        assert!(error.to_string().contains("no samples"));
        assert!(!output.exists());
        Ok(())
    }

    #[test]
    fn disabled_normalization_passes_through_ref_mismatch() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\tmismatch\tC\tT\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;

        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].ref_allele, "C");
        assert_eq!(records[0].id, "mismatch");
        Ok(())
    }

    #[test]
    fn missing_reference_index_is_rejected() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let args = interval_args(&input, &output, &reference, None, None);
        fs::remove_file(format!("{}.fai", reference.display()))?;
        let error = run(args).unwrap_err();
        assert!(error.to_string().contains("is not indexed"));
        assert!(!output.exists());
        Ok(())
    }

    #[test]
    fn empty_reference_index_is_accepted_but_malformed_lengths_are_rejected() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let index = format!("{}.fai", reference.display());

        fs::write(&index, "")?;
        run(interval_args(&input, &output, &reference, None, None))?;
        assert!(output.exists());
        fs::remove_file(&output)?;
        let output_index = PathBuf::from(format!("{}.tbi", output.display()));
        if output_index.exists() {
            fs::remove_file(output_index)?;
        }

        fs::write(&index, "chr1\tbad\t6\t5\t6\n")?;
        let error = run(interval_args(&input, &output, &reference, None, None)).unwrap_err();
        assert!(
            error.to_string().contains("invalid FASTA index length"),
            "{error:#}"
        );
        assert!(!output.exists());
        Ok(())
    }

    #[test]
    fn plain_vcf_output_failure_leaves_unindexed_output() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        args.threads = Some(1);
        let error = run(args).unwrap_err();
        assert!(error.to_string().contains("plain VCF output"));
        assert!(output.is_file());
        assert_eq!(vcf::load_raw_vcf(&output)?.1.len(), 1);
        assert!(!PathBuf::from(format!("{}.tbi", output.display())).exists());
        assert!(!PathBuf::from(format!("{}.csi", output.display())).exists());
        Ok(())
    }

    #[test]
    fn missing_output_parent_is_rejected_without_artifacts() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output_parent = directory.path().join("missing");
        let output = output_parent.join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let error = run(interval_args(&input, &output, &reference, None, None)).unwrap_err();
        assert!(error.to_string().contains("output parent does not exist"));
        assert!(!output_parent.exists());
        Ok(())
    }

    #[test]
    fn negative_window_is_accepted_when_blocksplit_is_unused() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        args.threads = Some(1);
        args.window_size = -1;
        run(args)?;

        assert_eq!(vcf::load_raw_vcf(&output)?.1.len(), 1);
        Ok(())
    }

    #[test]
    fn explicit_gender_controls_haploid_x_expansion() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let male_output = directory.path().join("male.vcf.gz");
        let female_output = directory.path().join("female.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chrX\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chrX,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chrX\t2\t.\tA\tC\t10\tPASS\t.\tGT:AD\t1:1,9\n",
            ),
        )?;
        let mut male = interval_args(&input, &male_output, &reference, None, None);
        male.gender = PreprocessGender::Male;
        run(male)?;
        let mut female = interval_args(&input, &female_output, &reference, None, None);
        female.gender = PreprocessGender::Female;
        run(female)?;
        let (_, male_records) = vcf::load_raw_vcf(&male_output)?;
        let (_, female_records) = vcf::load_raw_vcf(&female_output)?;
        assert!(male_records[0].samples[0].starts_with("1/1:"));
        assert!(female_records[0].samples[0].starts_with("0/1:"));
        Ok(())
    }

    #[test]
    fn bcftools_norm_excludes_ref_mismatches_deduplicates_and_left_aligns() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t3\tfirst\tA\tAA\t10\tPASS\t.\tGT\t0/1\n",
                "chr1\t3\tduplicate\tA\tAA\t20\tPASS\t.\tGT\t0/1\n",
                "chr1\t4\tmismatch\tC\tT\t30\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        args.bcftools_norm = true;
        run(args)?;
        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "first");
        assert_eq!(records[0].pos, 1);
        assert_eq!(records[0].ref_allele, "A");
        assert_eq!(records[0].alt_allele, "AA");
        Ok(())
    }

    #[test]
    fn logfile_and_verbose_emit_operational_messages() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        let logfile = directory.path().join("pre.log");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.logfile = Some(logfile.display().to_string());
        args.verbose = true;
        run(args)?;
        let log = fs::read_to_string(logfile)?;
        assert!(log.contains("INFO Preprocessing"));
        assert!(log.contains("INFO Wrote 0 records"));
        Ok(())
    }

    fn interval_args(
        input: &Path,
        output: &Path,
        reference: &Path,
        regions: Option<&Path>,
        targets: Option<&Path>,
    ) -> PreprocessArgs {
        let mut index = reference.as_os_str().to_os_string();
        index.push(".fai");
        if !Path::new(&index).is_file() {
            write_test_fai(reference).expect("test reference index should be writable");
        }
        PreprocessArgs {
            input: input.display().to_string(),
            output: output.display().to_string(),
            version: false,
            reference: Some(reference.display().to_string()),
            locations: None,
            pass_only: false,
            filters_only: None,
            regions_bedfile: regions.map(|path| path.display().to_string()),
            targets_bedfile: targets.map(|path| path.display().to_string()),
            fixchr: Some(false),
            no_fixchr: false,
            somatic: false,
            set_gt: None,
            filter_nonref: false,
            convert_gvcf_to_vcf: false,
            bcf: false,
            bcftools_norm: false,
            leftshift: true,
            no_leftshift: false,
            decompose: false,
            no_decompose: false,
            gender: PreprocessGender::Auto,
            window_size: 10_000,
            threads: None,
            logfile: None,
            verbose: false,
            quiet: false,
            force_interactive: false,
        }
    }

    fn write_test_fai(reference: &Path) -> Result<()> {
        let sequences = fasta::read_sequences(reference)?;
        let mut offset = 0u64;
        let mut index = String::new();
        for (name, sequence) in sequences {
            let line_bases = sequence.len().max(1);
            index.push_str(&format!(
                "{name}\t{}\t{offset}\t{line_bases}\t{}\n",
                sequence.len(),
                line_bases + 1
            ));
            offset += sequence.len() as u64 + name.len() as u64 + 3;
        }
        fs::write(format!("{}.fai", reference.display()), index)?;
        Ok(())
    }

    #[test]
    fn region_uses_gvcf_end_while_target_uses_start_position() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let boundary = directory.path().join("boundary.bed");
        let region_output = directory.path().join("region.vcf.gz");
        let target_output = directory.path().join("target.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "chr1\t3\t.\tC\tT\t.\tPASS\tEND=5\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;
        fs::write(&boundary, "chr1\t4\t5\n")?;

        run(interval_args(
            &input,
            &region_output,
            &reference,
            Some(&boundary),
            None,
        ))?;
        run(interval_args(
            &input,
            &target_output,
            &reference,
            None,
            Some(&boundary),
        ))?;

        assert_eq!(vcf::load_raw_vcf(&region_output)?.1.len(), 1);
        assert!(vcf::load_raw_vcf(&target_output)?.1.is_empty());
        Ok(())
    }

    #[test]
    fn selectors_remain_literal_after_fixchr_rewrites_records() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let selector = directory.path().join("selector.bed");
        let location_output = directory.path().join("location.vcf.gz");
        let region_output = directory.path().join("region.vcf.gz");
        let target_output = directory.path().join("target.vcf.gz");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(&selector, "1\t0\t5\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let mut location_args = interval_args(&input, &location_output, &reference, None, None);
        location_args.fixchr = Some(true);
        location_args.locations = Some("1:1-5".to_string());
        run(location_args)?;

        let mut region_args =
            interval_args(&input, &region_output, &reference, Some(&selector), None);
        region_args.fixchr = Some(true);
        run(region_args)?;

        let mut target_args =
            interval_args(&input, &target_output, &reference, None, Some(&selector));
        target_args.fixchr = Some(true);
        run(target_args)?;

        for output in [location_output, region_output, target_output] {
            assert!(
                vcf::load_raw_vcf(&output)?.1.is_empty(),
                "selector contig 1 must not be rewritten to match emitted chr1 records"
            );
        }
        Ok(())
    }

    #[test]
    fn leftshift_switch_controls_repeat_indel_normalization() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let shifted_output = directory.path().join("shifted.vcf.gz");
        let unchanged_output = directory.path().join("unchanged.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "chr1\t4\t.\tAA\tA\t30\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAAAAA\n")?;

        let mut shifted = interval_args(&input, &shifted_output, &reference, None, None);
        shifted.leftshift = true;
        run(shifted)?;
        let mut unchanged = interval_args(&input, &unchanged_output, &reference, None, None);
        unchanged.leftshift = false;
        run(unchanged)?;

        let shifted_record = vcf::load_raw_vcf(&shifted_output)?.1.remove(0);
        let unchanged_record = vcf::load_raw_vcf(&unchanged_output)?.1.remove(0);
        assert_eq!(shifted_record.pos, 1);
        assert_eq!(unchanged_record.pos, 4);
        Ok(())
    }

    #[test]
    fn tiny_parallel_inputs_do_not_create_window_block_boundaries() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAAAAAAAAAAAAAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=20>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t3\tbarrier\tA\tC\t30\tPASS\t.\tGT\t0/1\n",
                "chr1\t5\trepeat\tAA\tA\t30\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.threads = Some(2);
        args.window_size = 1;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(
            records.iter().map(|record| record.pos).collect::<Vec<_>>(),
            [3, 3]
        );
        Ok(())
    }

    #[test]
    fn location_start_reset_survives_a_dropped_homref_record() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(120)))?;
        fs::write(
            &input,
            format!(
                concat!(
                    "##fileformat=VCFv4.2\n",
                    "##contig=<ID=chr1,length=120>\n",
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                    "chr1\t2\tbarrier\t{}\tA\t30\tPASS\t.\tGT\t0/1\n",
                    "chr1\t52\thomref\tA\tC\t30\tPASS\t.\tGT\t0/0\n",
                    "chr1\t60\trepeat\tAA\tA\t30\tPASS\t.\tGT\t0/1\n",
                ),
                "A".repeat(50)
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.threads = Some(2);
        args.locations = Some("chr1:2-2,chr1:52-100".into());
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        let repeat = records
            .iter()
            .find(|record| record.ref_allele == "AA")
            .expect("repeat deletion must survive preprocessing");
        assert_eq!(repeat.pos, 1);
        Ok(())
    }

    #[test]
    fn symbolic_deletion_uses_end_and_reference_like_variant_reader() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=10>\n",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                // The deliberately incorrect one-base REF is accepted by the
                // legacy reader when END is present; it rebuilds the complete
                // deletion allele from the reference instead.
                "chr1\t2\trs1\tT\t<DEL>\t30\tPASS\tEND=5\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].to_line(),
            "chr1\t1\t.\tAACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0"
        );
        Ok(())
    }

    #[test]
    fn multi_allelic_symbolic_deletion_splits_like_variant_allele_splitter() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        let decomposed_output = directory.path().join("decomposed.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=10>\n",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\trs1\tT\t<DEL>,T\t30\tPASS\tEND=5\tGT\t1/2\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.decompose = false;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(
            records
                .iter()
                .map(vcf::RawVcfRecord::to_line)
                .collect::<Vec<_>>(),
            [
                "chr1\t1\t.\tAACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0",
                "chr1\t2\t.\tACCG\tT\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0",
            ]
        );

        let mut decomposed = interval_args(&input, &decomposed_output, &reference, None, None);
        decomposed.gender = PreprocessGender::None;
        decomposed.decompose = true;
        run(decomposed)?;
        let (_, records) = vcf::load_raw_vcf(&decomposed_output)?;
        assert_eq!(
            records
                .iter()
                .map(vcf::RawVcfRecord::to_line)
                .collect::<Vec<_>>(),
            [
                "chr1\t1\t.\tAACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0",
                "chr1\t2\t.\tA\tT\t30\t.\t.\tGT:AD:ADO:DP\t1/0:.,.:0:0",
                "chr1\t2\t.\tACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t1/0:.,.:0:0",
            ]
        );
        Ok(())
    }

    #[test]
    fn contig_start_symbolic_deletion_keeps_symbolic_alt() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=10>\n",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t1\trs1\tT\t<DEL>\t30\tPASS\tEND=4\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].to_line(),
            "chr1\t1\t.\tAACC\t<DEL>\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0"
        );
        Ok(())
    }

    #[test]
    fn legacy_headers_replace_conflicting_core_definitions() {
        let input = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"input AD\">".to_string(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"input DP\">".to_string(),
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"input END\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tINPUT".to_string(),
        ];
        let output = canonicalize_legacy_headers(&input);

        assert_eq!(output[0], "##fileformat=VCFv4.1");
        assert!(output.contains(
            &"##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depths\">".to_string()
        ));
        assert!(output.contains(&"##FORMAT=<ID=ADO,Number=.,Type=Integer,Description=\"Summed depth of non-called alleles.\">".to_string()));
        assert!(output.contains(
            &"##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">".to_string()
        ));
        assert!(!output.iter().any(|line| line.contains("input AD")
            || line.contains("input DP")
            || line.contains("input END")));
        assert_eq!(output.last(), input.last());
    }

    #[test]
    fn legacy_headers_sort_and_merge_non_core_input_declarations() {
        let input = vec![
            "##INFO=<ID=ZZ,Number=1,Type=String,Description=\"second\">".to_string(),
            "##FILTER=<ID=LowQ,Description=\"low quality\">".to_string(),
            "##INFO=<ID=ZZ,Number=1,Type=String,Description=\"first\">".to_string(),
            "##custom=z".to_string(),
            "##custom=a".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let output = canonicalize_legacy_headers(&input);
        let appended = &output[LEGACY_BASE_HEADERS.len()..output.len() - 1];

        assert_eq!(
            appended,
            [
                "##FILTER=<ID=LowQ,Description=\"low quality\">",
                "##INFO=<ID=ZZ,Number=1,Type=String,Description=\"first\">",
                "##custom=a",
                "##custom=z",
            ]
        );
    }

    #[test]
    fn ref_bytes_equal_handles_soft_masked_reference() {
        // Reference soft-masks repeat regions as lowercase; the VCF REF stays
        // uppercase. Legacy hap.py accepts this; we must too.
        assert!(ref_bytes_equal(b"t", b"T"));
        assert!(ref_bytes_equal(b"ACGT", b"acgt"));
    }

    #[test]
    fn breakend_normalization_preserves_remote_contig_spelling() {
        assert_eq!(
            uppercase_alleles_preserving_breakends("t]chr1:70],a"),
            "T]chr1:70],A"
        );

        let mut record = make_record(".");
        record.pos = 40;
        record.ref_allele = "T".to_string();
        record.alt_allele = "T]chr1:70]".to_string();
        record.info = "SVTYPE=BND".to_string();
        record.format = Some("GT:AD:GQ".to_string());
        record.samples = vec!["0/1:9,7:45".to_string()];
        normalize_bcftools_record(&mut record, b"ACGTACGTACGT");
        assert_eq!(record.alt_allele, "T]chr1:70]");
        assert!(materialize_unsupported_import_failure(&mut record));
        sort_info_keys(&mut record);
        assert_eq!(record.alt_allele, ".");
        assert_eq!(record.info, "END=40;IMPORT_FAIL;SVTYPE=BND");
        assert_eq!(record.samples, ["0/0:9:45"]);

        let mut symbolic = make_record("END=31");
        symbolic.pos = 31;
        symbolic.ref_allele = "T".to_string();
        symbolic.alt_allele = "<INS>".to_string();
        symbolic.format = Some("GT:AD".to_string());
        symbolic.samples = vec!["0/1:11,4".to_string()];
        assert!(materialize_unsupported_import_failure(&mut symbolic));
        assert_eq!(symbolic.alt_allele, ".");
        assert_eq!(symbolic.info, "END=31;IMPORT_FAIL");
        assert_eq!(symbolic.samples, ["0/0:11"]);
    }

    #[test]
    fn ref_bytes_equal_respects_n_wildcards() {
        assert!(ref_bytes_equal(b"N", b"A"));
        assert!(ref_bytes_equal(b"ACN", b"acg"));
    }

    #[test]
    fn ref_bytes_equal_rejects_real_differences() {
        assert!(!ref_bytes_equal(b"A", b"C"));
        assert!(!ref_bytes_equal(b"ACGT", b"ACGA"));
        assert!(!ref_bytes_equal(b"AC", b"ACG"));
    }

    fn make_record(info: &str) -> vcf::RawVcfRecord {
        vcf::RawVcfRecord {
            chrom: "chr1".into(),
            pos: 1,
            id: ".".into(),
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            qual: ".".into(),
            filter: "PASS".into(),
            info: info.into(),
            format: None,
            samples: vec![],
        }
    }

    proptest! {
        #[test]
        fn bcftools_normalization_is_idempotent(
            reference in proptest::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 8..40),
            alternate in proptest::collection::vec(prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')], 1..8),
            position_seed in 0usize..64,
            ref_len_seed in 1usize..8,
        ) {
            let position = position_seed % reference.len() + 1;
            let ref_len = ref_len_seed.min(reference.len() - position + 1);
            let mut record = make_record(".");
            record.pos = position;
            record.ref_allele = String::from_utf8(reference[position - 1..position - 1 + ref_len].to_vec()).unwrap();
            record.alt_allele = String::from_utf8(alternate).unwrap();
            normalize_bcftools_record(&mut record, &reference);
            let once = record.to_line();
            normalize_bcftools_record(&mut record, &reference);
            prop_assert_eq!(record.to_line(), once);
        }

        #[test]
        fn legacy_genotype_canonicalization_is_idempotent_across_ploidies_and_missing_calls(
            alleles in proptest::collection::vec(proptest::option::of(0usize..4), 1..=6),
            depth in 0usize..1000,
        ) {
            let mut record = make_record(".");
            record.format = Some("GT:DP".into());
            let genotype = alleles
                .iter()
                .map(|allele| allele.map_or_else(|| ".".to_string(), |value| value.to_string()))
                .collect::<Vec<_>>()
                .join("|");
            record.samples = vec![format!("{genotype}:{depth}")];
            canonicalize_legacy_genotypes(&mut record);
            let once = record.samples.clone();
            canonicalize_legacy_genotypes(&mut record);
            prop_assert_eq!(record.samples, once);
        }
    }

    #[test]
    fn normalized_records_restore_position_order_stably() {
        // Reduced ordering shape from chr21:44049615: a primitive emitted at
        // the trailing edge (26) preceded overlapping source records at 17.
        let mut trailing = make_record(".");
        trailing.chrom = "chr21".into();
        trailing.pos = 26;
        trailing.id = "trailing".into();
        let mut leading = trailing.clone();
        leading.pos = 15;
        leading.id = "leading".into();
        let mut overlap = trailing.clone();
        overlap.pos = 17;
        overlap.id = "overlap".into();
        let mut same_position = overlap.clone();
        same_position.id = "same-position".into();
        let mut next_contig = trailing.clone();
        next_contig.chrom = "chr1".into();
        next_contig.pos = 1;

        let mut records = vec![trailing, leading, overlap, same_position, next_contig];
        sort_normalized_records(&mut records);

        assert_eq!(
            records
                .iter()
                .map(|record| (record.chrom.as_str(), record.pos, record.id.as_str()))
                .collect::<Vec<_>>(),
            vec![
                ("chr21", 15, "leading"),
                ("chr21", 17, "overlap"),
                ("chr21", 17, "same-position"),
                ("chr21", 26, "trailing"),
                ("chr1", 1, "trailing"),
            ]
        );
    }

    #[test]
    fn blocksplit_aggregates_candidate_gaps_to_the_target_size() {
        let observations = (0..404)
            .map(|index| {
                let pos = (index / 101) * 10_000 + (index % 101) + 1;
                BlocksplitObservation {
                    chrom: "chr1".into(),
                    pos,
                    end: pos,
                    called: true,
                    location_groups: vec![0],
                }
            })
            .collect::<Vec<_>>();

        // Four groups produce three candidate gaps, but the second pass only
        // emits the middle one after cumulative candidate counts exceed the
        // 404 / 2 target. A naive every-gap reset would emit all three.
        assert_eq!(
            select_blocksplit_resets(&observations, 1, 2, None).reset_indices(),
            HashSet::from([202])
        );
    }

    #[test]
    fn blocksplit_preserves_negative_window_gap_arithmetic() {
        let observations = (0..101)
            .map(|_| BlocksplitObservation {
                chrom: "chr1".into(),
                pos: 1,
                end: 1,
                called: true,
                location_groups: vec![0],
            })
            .collect::<Vec<_>>();

        assert_eq!(
            select_blocksplit_resets(&observations, -1, 2, None).reset_indices(),
            HashSet::from([100])
        );
        assert!(
            select_blocksplit_resets(&observations, 0, 2, None)
                .reset_indices()
                .is_empty()
        );
    }

    #[test]
    fn blocksplit_ignores_uncalled_records_when_tracking_gaps() {
        let mut observations = (1..=101)
            .map(|pos| BlocksplitObservation {
                chrom: "chr1".into(),
                pos,
                end: pos,
                called: true,
                location_groups: vec![0],
            })
            .collect::<Vec<_>>();
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 10_000,
            end: 29_999,
            called: false,
            location_groups: vec![0],
        });
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 30_000,
            end: 30_000,
            called: true,
            location_groups: vec![0],
        });

        assert_eq!(
            select_blocksplit_resets(&observations, 1, 40, None).reset_indices(),
            HashSet::from([102])
        );
    }

    #[test]
    fn blocksplit_uses_effective_end_when_tracking_called_spans() {
        let mut observations = (1..=101)
            .map(|pos| BlocksplitObservation {
                chrom: "chr1".into(),
                pos,
                end: pos,
                called: true,
                location_groups: vec![0],
            })
            .collect::<Vec<_>>();
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 102,
            end: 10_000,
            called: true,
            location_groups: vec![0],
        });
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 10_001,
            end: 10_001,
            called: true,
            location_groups: vec![0],
        });

        assert!(
            select_blocksplit_resets(&observations, 1, 40, None)
                .reset_indices()
                .is_empty()
        );
    }

    #[test]
    fn blocksplit_partitions_explicit_locations_on_the_same_contig() {
        let observations = (0..404)
            .map(|index| {
                let location_group = index / 202;
                let index_within_location = index % 202;
                let cluster = index_within_location / 101;
                let pos = location_group * 1_000_000
                    + cluster * 10_000
                    + (index_within_location % 101)
                    + 1;
                BlocksplitObservation {
                    chrom: "chr1".into(),
                    pos,
                    end: pos,
                    called: true,
                    location_groups: vec![location_group],
                }
            })
            .collect::<Vec<_>>();

        // Legacy invokes blocksplit once per comma-separated location, so
        // each location computes its own total and target even on one contig.
        let locations = [
            vcf::LocationFilter::Contig("chr1".to_string()),
            vcf::LocationFilter::Contig("chr1".to_string()),
        ];
        assert_eq!(
            select_blocksplit_resets(&observations, 1, 2, Some(&locations)).reset_indices(),
            HashSet::from([101, 202, 303])
        );
    }

    #[test]
    fn parallel_comma_locations_duplicate_the_selected_stream_in_position_order() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(120)))?;
        let mut vcf_text = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=120>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        )
        .to_string();
        for pos in 1..=120 {
            vcf_text.push_str(&format!(
                "chr1\t{pos}\tv{pos}\tA\tC\t30\tPASS\t.\tGT\t0/1\n"
            ));
        }
        fs::write(&input, vcf_text)?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.locations = Some("chr1:1-100,chr1:51-120".to_string());
        args.threads = Some(2);
        run(args)?;

        let positions = vcf::load_raw_vcf(&output)?
            .1
            .into_iter()
            .map(|record| record.pos)
            .collect::<Vec<_>>();
        assert_eq!(positions.len(), 240);
        for (offset, pair) in positions.chunks_exact(2).enumerate() {
            assert_eq!(pair, [offset + 1, offset + 1]);
        }
        Ok(())
    }

    #[test]
    fn parallel_comma_locations_preserve_independent_multiblock_jobs() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(3_200)))?;
        let mut vcf_text = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=3200>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        )
        .to_string();
        for cluster_start in [1, 1_001, 2_001, 3_001] {
            for pos in cluster_start..=cluster_start + 100 {
                vcf_text.push_str(&format!(
                    "chr1\t{pos}\tv{pos}\tA\tC\t30\tPASS\t.\tGT\t0/1\n"
                ));
            }
        }
        fs::write(&input, vcf_text)?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.locations = Some("chr1:1-2101,chr1:1001-3101".to_string());
        args.threads = Some(2);
        args.window_size = 10;
        run(args)?;

        let positions = vcf::load_raw_vcf(&output)?
            .1
            .into_iter()
            .map(|record| record.pos)
            .collect::<Vec<_>>();
        assert_eq!(positions.len(), 705);
        assert!(positions.windows(2).all(|pair| pair[0] <= pair[1]));

        let mut multiplicities = std::collections::BTreeMap::new();
        for pos in positions {
            *multiplicities.entry(pos).or_insert(0usize) += 1;
        }
        for pos in 1..=101 {
            assert_eq!(multiplicities.get(&pos), Some(&2));
        }
        for pos in 1_001..=1_101 {
            assert_eq!(multiplicities.get(&pos), Some(&2));
        }
        for pos in 2_001..=2_100 {
            assert_eq!(multiplicities.get(&pos), Some(&2));
        }
        assert_eq!(multiplicities.get(&2_101), Some(&1));
        for pos in 3_001..=3_100 {
            assert_eq!(multiplicities.get(&pos), Some(&1));
        }
        assert_eq!(multiplicities.get(&3_101), None);
        Ok(())
    }

    #[test]
    fn blocksplit_filters_inactive_partitions_and_preserves_all_empty_fallback() {
        let mut observations = vec![
            BlocksplitObservation {
                chrom: "chr1".into(),
                pos: 1,
                end: 1,
                called: true,
                location_groups: vec![0],
            },
            BlocksplitObservation {
                chrom: "chr1".into(),
                pos: 2,
                end: 2,
                called: false,
                location_groups: vec![0],
            },
            BlocksplitObservation {
                chrom: "chr2".into(),
                pos: 1,
                end: 1,
                called: false,
                location_groups: vec![0],
            },
        ];

        let mixed = select_blocksplit_resets(&observations, 1, 2, None);
        assert_eq!(mixed.included_indices(), Some(HashSet::from([0, 1])));

        observations[0].called = false;
        let all_empty = select_blocksplit_resets(&observations, 1, 2, None);
        assert_eq!(all_empty.included_indices(), None);
    }

    #[test]
    fn blocksplit_parallelism_resolves_explicit_and_default_thread_counts() {
        assert_eq!(effective_thread_count_with_available(Some(0), 8), 0);
        assert_eq!(effective_thread_count_with_available(Some(1), 8), 1);
        assert_eq!(effective_thread_count_with_available(Some(2), 1), 2);
        assert_eq!(effective_thread_count_with_available(None, 0), 1);
        assert_eq!(effective_thread_count_with_available(None, 1), 1);
        assert_eq!(effective_thread_count_with_available(None, 2), 2);
    }

    #[test]
    fn strip_stale_info_keys_removes_allele_counts() {
        let mut record = make_record("AC=1;AF=0.5;AN=2;DP=158");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, "AF=0.5;DP=158");
    }

    #[test]
    fn strip_stale_info_keys_handles_missing_info() {
        let mut record = make_record(".");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, ".");
    }

    #[test]
    fn strip_stale_info_keys_collapses_to_dot_when_all_stripped() {
        let mut record = make_record("AC=1;AN=2;MLEAC=1;MLEAF=0.5");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, ".");
    }

    #[test]
    fn strip_stale_info_keys_preserves_flags_without_values() {
        let mut record = make_record("SOMATIC;AC=1;DP=158");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, "SOMATIC;DP=158");
    }

    #[test]
    fn sort_info_keys_orders_alphabetically() {
        let mut record = make_record("HRun=0;AF=0.5;DP=158;Dels=0");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "AF=0.5;DP=158;Dels=0;HRun=0");
    }

    #[test]
    fn sort_info_keys_handles_flag_entries() {
        let mut record = make_record("SOMATIC;DP=158;AF=0.5");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "AF=0.5;DP=158;SOMATIC");
    }

    #[test]
    fn sort_info_keys_collapses_negative_zero_to_zero() {
        // Legacy hap.py reads INFO floats through htslib which loses the
        // sign on negative zero before re-emitting. Match that here.
        let mut record = make_record("SB=-0;FS=12.24;MQRankSum=0");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "FS=12.24;MQRankSum=0;SB=0");
    }

    #[test]
    fn sort_info_keys_preserves_non_zero_negatives() {
        let mut record = make_record("SB=-1.5;BaseQRankSum=-0.123");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "BaseQRankSum=-0.123;SB=-1.5");
    }

    #[test]
    fn expand_haploid_gt_matches_legacy_output() {
        // Haploid alt calls expand to heterozygous half-call, mirroring
        // VariantAlleleSplitter.cpp:180-227. Empirically confirmed: legacy
        // result.vcf.gz shows query GT=1 on autosomes as 0/1:het, not 1/1.
        assert_eq!(expand_haploid_gt("1", false), "0/1");
        assert_eq!(expand_haploid_gt("2", false), "0/2");
        assert_eq!(expand_haploid_gt("1", true), "1/1");
        assert_eq!(expand_haploid_gt("2", true), "2/2");
        assert_eq!(expand_haploid_gt("0", true), "0/0");
        assert_eq!(expand_haploid_gt(".", true), "./.");
        // Diploid inputs pass through.
        assert_eq!(expand_haploid_gt("0/1", true), "0/1");
        assert_eq!(expand_haploid_gt("1|2", true), "1|2");
        assert_eq!(expand_haploid_gt("./.", true), "./.");
    }

    #[test]
    fn active_preprocessing_masks_genotypes_wider_than_diploid() {
        let mut record = make_record(".");
        record.format = Some("GT:GQ".to_string());
        record.samples = vec![
            "0/1:40".to_string(),
            "0/1/1:50".to_string(),
            "0|0|1|1:60".to_string(),
        ];

        mask_genotypes_wider_than_diploid(&mut record);

        assert_eq!(record.samples, ["0/1:40", ".:50", ".:60"]);
    }

    #[test]
    fn sort_info_keys_collapses_negative_zero_in_lists() {
        let mut record = make_record("AF=-0,0.5");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "AF=0,0.5");
    }

    #[test]
    fn collapse_pl_reduces_to_last_value() {
        let mut record = make_record(".");
        record.format = Some("GT:AD:PL".to_string());
        record.samples = vec!["0/1:32,126:3279,0,716".to_string()];
        collapse_pl_to_last_value(&mut record);
        assert_eq!(record.samples[0], "0/1:32,126:716");
    }

    #[test]
    fn collapse_pl_leaves_scalar_alone() {
        let mut record = make_record(".");
        record.format = Some("GT:PL".to_string());
        record.samples = vec!["0/1:716".to_string()];
        collapse_pl_to_last_value(&mut record);
        assert_eq!(record.samples[0], "0/1:716");
    }

    #[test]
    fn ado_zero_when_genotype_covers_all_alleles() {
        assert_eq!(compute_ado("32,126", Some("0/1")), 0);
        assert_eq!(compute_ado("32,126", Some("0|1")), 0);
    }

    #[test]
    fn ado_captures_ref_depth_for_hom_alt() {
        // GT=1/1 with AD=22,413 → ADO = AD[0] = 22 (the unused reference depth)
        assert_eq!(compute_ado("22,413", Some("1/1")), 22);
        assert_eq!(compute_ado("1,17", Some("1/1")), 1);
    }

    #[test]
    fn ado_captures_alt_depth_for_hom_ref() {
        assert_eq!(compute_ado("100,25", Some("0/0")), 25);
    }

    #[test]
    fn ado_handles_multi_allelic() {
        // GT=1/2 across original alleles — both indices 1 and 2 are called,
        // so AD[0] (ref) is the only "other" depth.
        assert_eq!(compute_ado("10,50,70", Some("1/2")), 10);
        // GT=0/2 — index 1 is the "other" allele.
        assert_eq!(compute_ado("10,50,70", Some("0/2")), 50);
    }

    #[test]
    fn ado_defaults_to_zero_on_missing_gt_or_ad() {
        assert_eq!(compute_ado(".", Some("0/1")), 0);
        assert_eq!(compute_ado("10,20", Some("./.")), 0);
        assert_eq!(compute_ado("10,20", None), 0);
    }

    #[test]
    fn classify_value_distinguishes_int_float_string() {
        assert_eq!(classify_value("99"), ScalarType::Integer);
        assert_eq!(classify_value("3279,0,716"), ScalarType::Integer);
        assert_eq!(classify_value("95.77"), ScalarType::Float);
        assert_eq!(classify_value("0.934"), ScalarType::Float);
        assert_eq!(classify_value("PASS"), ScalarType::String);
        assert_eq!(classify_value("."), ScalarType::Integer);
    }

    #[test]
    fn reorder_format_keeps_gt_ad_ado_dp_first_then_buckets_by_type() {
        // Record with GQ=95.77 (float) should push GQ after the integer
        // bucket — legacy canonical shape `GT:AD:ADO:DP:GQX:MQ:PL:GQ:VF`.
        let mut record = make_record(".");
        record.format = Some("GT:AD:ADO:DP:GQ:GQX:MQ:PL:VF".to_string());
        record.samples = vec!["0/1:22,306:0:328:95.77:96:36:96:0.933".to_string()];
        reorder_format_fields(&mut record);
        assert_eq!(
            record.format.as_deref().unwrap(),
            "GT:AD:ADO:DP:GQX:MQ:PL:GQ:VF"
        );
        assert_eq!(record.samples[0], "0/1:22,306:0:328:96:36:96:95.77:0.933");
    }

    #[test]
    fn reorder_format_treats_integer_gq_as_integer() {
        // When the sample's GQ happens to parse as an int (e.g. 99), legacy
        // places it with the ints — alphabetical means GQ < GQX.
        let mut record = make_record(".");
        record.format = Some("GT:AD:ADO:DP:GQ:GQX:MQ:PL:VF".to_string());
        record.samples = vec!["1/1:0,101:0:101:99:99:51:0:1".to_string()];
        reorder_format_fields(&mut record);
        assert_eq!(
            record.format.as_deref().unwrap(),
            "GT:AD:ADO:DP:GQ:GQX:MQ:PL:VF"
        );
    }

    #[test]
    fn somatic_conversion_uses_legacy_gt_modes_and_preserves_sample_formats_in_info() {
        let mut record = make_record("SOMATIC");
        record.alt_allele = "C,G".to_string();
        record.format = Some("GT:DP:VF".to_string());
        record.samples = vec!["0/1:12:0.25".to_string(), "1/2:30:0.75".to_string()];
        let sample_names = vec!["NORMAL".to_string(), "TUMOR".to_string()];

        let half = convert_somatic_record(&record, SomaticGtMode::Half, &sample_names);
        assert_eq!(half.len(), 2);
        assert_eq!(half[0].alt_allele, "C");
        assert_eq!(half[1].alt_allele, "G");
        assert_eq!(half[0].format.as_deref(), Some("GT:AD:DP"));
        assert_eq!(half[0].samples, ["./1:.,.:0"]);
        assert!(half[0].info.contains("NORMAL_GT=2,4"));
        assert!(half[0].info.contains("NORMAL_DP=12"));
        assert!(half[0].info.contains("NORMAL_VF=0.25"));
        assert!(!half[0].info.contains("TUMOR_GT="));
        assert!(!half[0].info.contains("TUMOR_DP=30"));
        assert!(!half[0].info.contains("TUMOR_VF=0.75"));

        let hemi = convert_somatic_record(&record, SomaticGtMode::Hemi, &sample_names);
        assert_eq!(hemi[0].samples, ["1:.,.:0"]);
        let het = convert_somatic_record(&record, SomaticGtMode::Het, &sample_names);
        assert_eq!(het[0].samples, ["0/1:.,.:0"]);
        let hom = convert_somatic_record(&record, SomaticGtMode::Hom, &sample_names);
        assert_eq!(hom[0].samples, ["1/1:.,.:0"]);
    }

    #[test]
    fn somatic_first_mode_preserves_first_gt_and_does_not_split_alts() {
        let mut record = make_record(".");
        record.alt_allele = "C,G".to_string();
        record.format = Some("GT:DP".to_string());
        record.samples = vec!["2|1:12".to_string(), "0/2:30".to_string()];
        let sample_names = vec!["NORMAL".to_string(), "TUMOR".to_string()];

        let converted = convert_somatic_record(&record, SomaticGtMode::First, &sample_names);
        assert_eq!(converted.len(), 1);
        assert_eq!(converted[0].alt_allele, "C,G");
        assert_eq!(converted[0].format.as_deref(), Some("GT:AD:DP"));
        assert_eq!(converted[0].samples, ["2|1:.,.,.:0"]);
        assert!(converted[0].info.contains("NORMAL_GT=6,5"));
        assert!(!converted[0].info.contains("TUMOR_GT="));
    }

    #[test]
    fn somatic_partial_credit_merges_nonhom_siblings_but_not_hom_calls() {
        let mut record = make_record("TAG=multi");
        record.alt_allele = "T,G".to_string();
        record.format = Some("GT:DP".to_string());
        record.samples = vec!["1/2:30".to_string()];
        let names = vec!["TUMOR".to_string()];

        let half = finalize_somatic_records(
            convert_somatic_record(&record, SomaticGtMode::Half, &names),
            SomaticGtMode::Half,
        );
        assert_eq!(half.len(), 1);
        assert_eq!(half[0].alt_allele, "G,T");
        assert_eq!(half[0].samples, ["2/1:.,.,.:0"]);

        let hom = finalize_somatic_records(
            convert_somatic_record(&record, SomaticGtMode::Hom, &names),
            SomaticGtMode::Hom,
        );
        assert_eq!(hom.len(), 2);
        assert_eq!(hom[0].samples, ["1/1:.,.:0"]);
    }

    #[test]
    fn somatic_half_calls_bypass_partial_credit_for_scmp() {
        let mut record = make_record(".");
        record.alt_allele = "C,G".to_string();
        record.format = Some("GT".to_string());
        record.samples = vec!["1/2".to_string()];
        let names = vec!["TUMOR".to_string()];
        let converted = convert_somatic_record(&record, SomaticGtMode::Half, &names);

        let scmp = finalize_somatic_for_pipeline(converted.clone(), SomaticGtMode::Half, false);
        assert_eq!(scmp.len(), 2);
        assert_eq!(scmp[0].format.as_deref(), Some("GT"));
        assert_eq!(scmp[0].samples, ["./1"]);
        assert_eq!(scmp[1].samples, ["./1"]);

        let normalized = finalize_somatic_for_pipeline(converted, SomaticGtMode::Half, true);
        assert_eq!(normalized.len(), 1);
        assert_eq!(normalized[0].samples, ["2/1:.,.,.:0"]);
    }

    #[test]
    fn somatic_headers_declare_sample_prefixed_format_info_fields() {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL".to_string(),
        ];
        let names = somatic_info_sample_names(&headers);
        assert_eq!(names, ["TUMOR", "NORMAL"]);
        append_somatic_info_headers(&mut headers, &names);

        assert!(headers.contains(
            &"##INFO=<ID=TUMOR_GT,Number=1,Type=String,Description=\"Genotype\">".to_string()
        ));
        assert!(headers.contains(
            &"##INFO=<ID=NORMAL_DP,Number=1,Type=Integer,Description=\"Depth\">".to_string()
        ));
    }

    #[test]
    fn somatic_single_input_sample_uses_output_sample_prefix() {
        let headers =
            vec!["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tORIGINAL".to_string()];
        assert_eq!(somatic_info_sample_names(&headers), ["SAMPLE"]);
    }

    #[test]
    fn non_ref_filter_drops_only_records_whose_gt_calls_the_non_ref_allele() {
        let mut record = make_record(".");
        record.alt_allele = "C,<NON_REF>".to_string();
        record.format = Some("GT:DP".to_string());
        record.samples = vec!["0/1:12".to_string(), "0/0:30".to_string()];
        assert!(!calls_non_ref_allele(&record));

        record.samples[1] = "0/2:30".to_string();
        assert!(calls_non_ref_allele(&record));
    }

    #[test]
    fn uncalled_non_ref_is_trimmed_and_reference_blocks_are_dropped() {
        let mut variant = make_record(".");
        variant.alt_allele = "C,<NON_REF>".to_string();
        variant.format = Some("GT:AD".to_string());
        variant.samples = vec!["0/1:7,8,0".to_string()];
        assert!(trim_uncalled_non_ref(&mut variant));
        assert_eq!(variant.alt_allele, "C");
        assert_eq!(variant.samples, ["0/1:7,8"]);

        let mut block = make_record("END=10");
        block.alt_allele = "<NON_REF>".to_string();
        block.format = Some("GT".to_string());
        block.samples = vec!["0/0".to_string()];
        assert!(!trim_uncalled_non_ref(&mut block));
    }

    #[test]
    fn called_alt_projection_remaps_gt_and_ad_and_drops_homref_records() {
        let mut called = make_record(".");
        called.alt_allele = "C,G,T".to_string();
        called.format = Some("GT:AD:ADO".to_string());
        called.samples = vec!["0|2:10,1,8,2:3".to_string()];
        assert!(retain_called_alternates(&mut called));
        assert_eq!(called.alt_allele, "G");
        assert_eq!(called.samples, ["0|1:10,8:3"]);

        let mut homref = make_record(".");
        homref.alt_allele = "C,G".to_string();
        homref.format = Some("GT:AD".to_string());
        homref.samples = vec!["0/0:12,0,0".to_string()];
        assert!(!retain_called_alternates(&mut homref));

        let mut spanning_deletion = make_record(".");
        spanning_deletion.alt_allele = "T,*".to_string();
        spanning_deletion.format = Some("GT:AD".to_string());
        spanning_deletion.samples = vec!["1/2:8,5,6".to_string()];
        assert!(retain_called_alternates(&mut spanning_deletion));
        assert_eq!(spanning_deletion.alt_allele, "T");
        assert_eq!(spanning_deletion.samples, ["0/1:8,5"]);
    }

    #[test]
    fn preprocessing_removes_uncalled_multi_alleles_before_primitive_split() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=20>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "chr1\t2\t.\tA\tATC,ATCTC\t.\tPASS\t.\tGT:AD\t2|0:5,0,7\n",
                "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t0/0:9,0,0\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAAAAAAAAAAAAAAAAAAAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.decompose = true;
        args.window_size = 4096;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].pos, 2);
        assert_eq!(records[0].ref_allele, "A");
        assert_eq!(records[0].alt_allele, "ATCTC");
        assert_eq!(records[0].qual, "0");
        assert_eq!(records[0].samples, ["0/1:5,7:0:0"]);
        Ok(())
    }

    #[test]
    fn legacy_writer_unphases_and_canonicalizes_genotypes() {
        let mut record = make_record(".");
        record.format = Some("GT".to_string());
        record.samples = vec![
            "1|1".to_string(),
            "1|0".to_string(),
            "0|1".to_string(),
            "1|2".to_string(),
            "2/1".to_string(),
        ];
        canonicalize_legacy_genotypes(&mut record);
        assert_eq!(record.samples, ["1/1", "0/1", "0/1", "2/1", "2/1"]);
    }

    #[test]
    fn multi_allelic_order_remaps_gt_and_ad_like_legacy_aggregator() {
        let mut record = make_record(".");
        record.alt_allele = "T,G".to_string();
        record.format = Some("GT:AD".to_string());
        record.samples = vec!["2/1:10,11,9".to_string(), "0/1:14,14,0".to_string()];
        canonicalize_multi_allelic_order(&mut record);
        assert_eq!(record.alt_allele, "G,T");
        assert_eq!(record.samples[0], "2/1:10,9,11");
        assert_eq!(record.samples[1], "0/2:14,0,14");
    }

    #[test]
    fn secondary_samples_keep_core_fields_but_not_dynamic_annotations() {
        let mut record = make_record(".");
        record.format = Some("GT:AD:ADO:DP:GQ:TXT".to_string());
        record.samples = vec![
            "0/1:7,8:0:15:20:variant".to_string(),
            "0/0:16,0:0:16:25:normal".to_string(),
        ];
        let string_fields = BTreeSet::from(["GT".to_string(), "TXT".to_string()]);
        blank_secondary_sample_annotations(&mut record, false, &string_fields);
        assert_eq!(record.samples[0], "0/1:7,8:0:15:20:variant");
        assert_eq!(record.samples[1], "0/0:.,.:.:16:.:.");

        record.samples[1] = "0/0:16,0:0:16:25:normal".to_string();
        blank_secondary_sample_annotations(&mut record, true, &string_fields);
        assert_eq!(record.samples[1], "0/0:.,.:.:16:.:");
    }

    #[test]
    fn gvcf_conversion_trims_non_ref_and_unrequested_fields() {
        let mut record = make_record("END=10;AC=1;DP=40");
        record.alt_allele = "C,<NON_REF>".to_string();
        record.format = Some("GT:DP:GQ:AD:PL".to_string());
        record.samples = vec!["0/1:12:40:8,4,0:0,10,100,20,200,300".to_string()];

        assert!(convert_gvcf_record(&mut record));
        assert_eq!(record.alt_allele, "C");
        assert_eq!(record.info, ".");
        assert_eq!(record.format.as_deref(), Some("GT:DP:GQ"));
        assert_eq!(record.samples, ["0/1:12:40"]);
    }

    #[test]
    fn gvcf_conversion_strips_unrequested_info_and_format_headers() {
        let mut headers = vec![
            "##INFO=<ID=TAG,Number=1,Type=String,Description=\"Tag\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">".to_string(),
            "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"GQ\">".to_string(),
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">".to_string(),
            "##FORMAT=<ID=TXT,Number=1,Type=String,Description=\"TXT\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS".to_string(),
        ];
        filter_gvcf_headers(&mut headers);
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=GT,"))
        );
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=DP,"))
        );
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=GQ,"))
        );
        assert!(!headers.iter().any(|line| line.starts_with("##INFO=")));
        assert!(
            !headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=AD,"))
        );
        assert!(
            !headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=TXT,"))
        );
    }

    #[test]
    fn gvcf_conversion_drops_single_alt_blocks_and_non_ref_calls() {
        let mut reference_block = make_record("END=10");
        reference_block.alt_allele = "<NON_REF>".to_string();
        reference_block.format = Some("GT:DP:GQ".to_string());
        reference_block.samples = vec!["0/0:12:40".to_string()];
        assert!(!convert_gvcf_record(&mut reference_block));

        let mut non_ref_call = make_record(".");
        non_ref_call.alt_allele = "C,<NON_REF>".to_string();
        non_ref_call.format = Some("GT:DP:GQ".to_string());
        non_ref_call.samples = vec!["0/2:12:40".to_string()];
        assert!(!convert_gvcf_record(&mut non_ref_call));
    }

    #[test]
    fn gvcf_conversion_retains_initial_multi_alt_homref_as_reference_record() {
        let mut record = make_record("END=10");
        record.alt_allele = "C,<NON_REF>".to_string();
        record.format = Some("GT:DP:GQ".to_string());
        record.samples = vec!["0/0:12:40".to_string()];

        assert!(convert_gvcf_record(&mut record));
        assert_eq!(record.alt_allele, ".");
        assert_eq!(record.info, ".");
        assert_eq!(record.samples, ["0/0:12:40"]);
    }

    #[test]
    fn parallel_blocksplit_omits_empty_partitions_but_falls_back_when_all_are_empty() -> Result<()>
    {
        let directory = tempdir()?;
        let reference = directory.path().join("ref.fa");
        let mixed_input = directory.path().join("mixed.vcf");
        let mixed_output = directory.path().join("mixed.out.vcf.gz");
        let empty_input = directory.path().join("empty.vcf");
        let empty_output = directory.path().join("empty.out.vcf.gz");
        fs::write(&reference, ">chr1\nAAAAA\n>chr2\nAAAAA\n")?;
        let header = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=5>\n",
            "##contig=<ID=chr2,length=5>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        );
        fs::write(
            &mixed_input,
            format!(
                "{header}{}{}",
                "chr1\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/1:12:40\n",
                "chr2\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/0:12:40\n",
            ),
        )?;
        fs::write(
            &empty_input,
            format!(
                "{header}{}{}",
                "chr1\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/0:12:40\n",
                "chr2\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/0:12:40\n",
            ),
        )?;

        let mut mixed_args = interval_args(&mixed_input, &mixed_output, &reference, None, None);
        mixed_args.gender = PreprocessGender::None;
        mixed_args.threads = Some(2);
        mixed_args.convert_gvcf_to_vcf = true;
        run(mixed_args)?;
        let mixed_records = vcf::load_raw_vcf(&mixed_output)?.1;
        assert_eq!(
            mixed_records
                .iter()
                .map(|record| record.chrom.as_str())
                .collect::<Vec<_>>(),
            ["chr1"]
        );

        let mut empty_args = interval_args(&empty_input, &empty_output, &reference, None, None);
        empty_args.gender = PreprocessGender::None;
        empty_args.threads = Some(2);
        empty_args.convert_gvcf_to_vcf = true;
        run(empty_args)?;
        let empty_records = vcf::load_raw_vcf(&empty_output)?.1;
        assert_eq!(
            empty_records
                .iter()
                .map(|record| (record.chrom.as_str(), record.alt_allele.as_str()))
                .collect::<Vec<_>>(),
            [("chr1", "."), ("chr2", ".")]
        );
        Ok(())
    }
}
