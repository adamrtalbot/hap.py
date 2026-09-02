use crate::application::{PreprocessArgs, PreprocessGender, ValidatedPreprocessArgs};
use crate::domain::{PrimitiveIdentity, QueryProvenance, RawVcfRecord};
use crate::{
    adapters::{fasta, vcf},
    engines::variant_pipeline,
    output::OutputTransaction,
};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

mod alleles;
mod blocksplit;
mod canonical;
mod compatibility;
mod genotype;
mod normalization;
mod options;
mod streaming;

use alleles::*;
use blocksplit::*;
use canonical::*;
use normalization::*;
use options::*;
use streaming::{LocationAggregatedRecords, PreparedRecordSpool, PreprocessSpool};

pub(crate) use alleles::calls_non_ref_allele;
pub(crate) use canonical::{canonicalize_legacy_headers, structured_header_identity};
pub(crate) use options::infer_gender;

/// Per-sample left-shift boundary used by legacy hap.py's C++ preprocess
/// binary — indels can slide at most 1 kbp leftward before the pass stops.
const LEFT_SHIFT_WINDOW: usize = 1024;

/// `blocksplit` does not create a preprocessing boundary until the candidate
/// block contains more than its default minimum of 100 called variants.
const LEGACY_MIN_BLOCK_VARIANTS: usize = 100;
const LEGACY_MAX_BLOCKS: usize = 40;
const LOCATION_STREAM_POLICY: compatibility::LocationStreamPolicy =
    compatibility::LocationStreamPolicy::IndependentLegacyStreams;

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
    included_indices: Vec<usize>,
    /// Resets can land on records discarded later by VariantCallsOnly. They
    /// must therefore be applied immediately after selecting the prepared
    /// input record, before genotype-based filtering.
    reset_before_indices: Vec<usize>,
}

type RecordIdentity = (String, usize, String, String);

fn assign_passthrough_primitive_identity(record: &mut RawVcfRecord) {
    if record.alt_allele.contains(',') {
        return;
    }
    record.primitive_identity = Some(PrimitiveIdentity {
        start: record.pos,
        end: record.pos + record.ref_allele.len().saturating_sub(1),
        alt: record.alt_allele.clone(),
    });
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

/// Cross-record deletion bookkeeping gathered during the input pre-scan and
/// consumed during normalization. Legacy `VariantInput` tracks these five
/// identity sets together, so they travel as one cohesive bundle.
#[derive(Default)]
struct DeletionSets {
    following_spanning: HashSet<RecordIdentity>,
    equal_floor_blocked: HashSet<RecordIdentity>,
    released_mixed: HashSet<RecordIdentity>,
    released_following: HashSet<RecordIdentity>,
    mixed_deletion_merge_partners: HashSet<RecordIdentity>,
}

/// Per-chromosome left-shift state carried across records during normalization.
///
/// `prev_end` is the running maximum reference end feeding the primitive
/// splitter (unchanged legacy semantics). The remaining fields compute the
/// left-shift floor — how far left a record may slide — as the maximum of:
///
/// * `barrier`: the end of a *non-insertion* record (deletion / substitution /
///   complex) at a strictly earlier original position. A left-shift may not
///   cross a base another variant deleted or changed.
/// * substitution ends at the *same* original position: a SNP sharing a
///   position with a deletion still blocks it (chr21:9920194). These come from
///   `NormalizeContext::substitution_end_by_position`, pre-scanned so the
///   barrier is order-invariant regardless of which record is listed first.
///
/// A pure insertion never floors anything, and records at the same position do
/// not floor each other unless one is a substitution — so a colocated het
/// insertion + deletion both reach their shared anchor and the location
/// aggregator re-merges them into one het-alt, order-invariantly.
/// `group_pos`/`pending` stage the current position's non-insertion end until a
/// later position commits it into `barrier`.
#[derive(Clone, Copy, Default)]
struct ShiftFloors {
    prev_end: usize,
    barrier: usize,
    group_pos: usize,
    pending: usize,
}

/// A pure insertion extends the REF prefix on every ALT — it adds bases without
/// deleting or changing any, so it never blocks a neighbour's left-shift.
fn record_is_pure_insertion(record: &crate::domain::RawVcfRecord) -> bool {
    let reference = record.ref_allele.as_bytes();
    record
        .alt_allele
        .split(',')
        .all(|alt| alt.len() > reference.len() && alt.as_bytes().starts_with(reference))
}

/// A substitution changes an existing reference base (SNP, MNP, complex): some
/// ALT is neither a REF prefix (deletion) nor a REF-prefixed extension
/// (insertion). These block a left-shift even at the same position.
fn record_is_substitution(record: &crate::domain::RawVcfRecord) -> bool {
    let reference = record.ref_allele.as_bytes();
    record.alt_allele.split(',').any(|alt| {
        let alt = alt.as_bytes();
        !alt.starts_with(reference) && !reference.starts_with(alt)
    })
}

/// Per-run normalization configuration shared by the block-split workers and
/// the single-threaded fallback. Every reference borrows caller-owned state for
/// the duration of the normalization phase.
struct NormalizeContext<'a> {
    args: &'a PreprocessArgs,
    gender: PreprocessGender,
    leftshift: bool,
    decompose: bool,
    somatic_mode: Option<crate::application::SomaticGtMode>,
    string_format_fields: &'a BTreeSet<String>,
    reference_sequences: &'a std::collections::BTreeMap<String, String>,
    deletions: &'a DeletionSets,
    /// Greatest reference end of a substitution at each `(chrom, pos)`. A
    /// deletion sharing a position with a SNP is floored here regardless of the
    /// two records' input order, so the barrier stays order-invariant
    /// (chr21:9920194). Built in the single input pre-scan.
    substitution_end_by_position: &'a HashMap<(String, usize), usize>,
}

fn observe_following_spanning_deletions(
    record: &crate::domain::RawVcfRecord,
    mixed_deletions_by_position: &mut HashMap<(String, usize), Vec<(RecordIdentity, usize)>>,
    deletions_by_position: &mut HashMap<(String, usize), Vec<(usize, usize)>>,
    sets: &mut DeletionSets,
) {
    let DeletionSets {
        following_spanning: following_spanning_deletions,
        equal_floor_blocked: equal_floor_blocked_deletions,
        released_mixed: released_mixed_deletions,
        released_following: released_following_deletions,
        mixed_deletion_merge_partners,
    } = sets;
    if let Some(previous_position) = record.pos.checked_sub(1)
        && let Some(candidates) =
            mixed_deletions_by_position.get(&(record.chrom.clone(), previous_position))
    {
        for (candidate, mixed_deleted_length) in candidates {
            let following_deleted_length = record
                .alt_allele
                .split(',')
                .filter(|alternate| {
                    record.ref_allele.len() > alternate.len()
                        && record.ref_allele.starts_with(alternate)
                })
                .map(|alternate| record.ref_allele.len() - alternate.len())
                .min();
            if !released_mixed_deletions.contains(candidate)
                && (record.ref_allele.len() > *mixed_deleted_length
                    || following_deleted_length.is_some())
            {
                following_spanning_deletions.insert(candidate.clone());
            }
            if following_deleted_length
                .is_some_and(|deleted_length| deleted_length < *mixed_deleted_length)
            {
                following_spanning_deletions.remove(candidate);
                released_mixed_deletions.insert(candidate.clone());
                let following_identity = (
                    record.chrom.clone(),
                    record.pos,
                    record.ref_allele.clone(),
                    record.alt_allele.clone(),
                );
                released_following_deletions.insert(following_identity.clone());
                mixed_deletion_merge_partners.insert(following_identity);
            } else if equal_floor_blocked_deletions.contains(candidate)
                && following_deleted_length.is_some()
            {
                mixed_deletion_merge_partners.insert((
                    record.chrom.clone(),
                    record.pos,
                    record.ref_allele.clone(),
                    record.alt_allele.clone(),
                ));
            }
        }
    }

    if record.ref_allele.len() <= 1 {
        return;
    }
    let has_mixed_deletion = record.alt_allele.split(',').any(|alternate| {
        alternate.len() == 1
            && record
                .ref_allele
                .as_bytes()
                .first()
                .zip(alternate.as_bytes().first())
                .is_some_and(|(reference, alternate)| reference != alternate)
    });
    if has_mixed_deletion {
        let identity = (
            record.chrom.clone(),
            record.pos,
            record.ref_allele.clone(),
            record.alt_allele.clone(),
        );
        if record.pos > 1
            && deletions_by_position
                .get(&(record.chrom.clone(), record.pos - 1))
                .is_some_and(|deletions| {
                    deletions.iter().any(|(end, reference_length)| {
                        *end >= record.pos && *reference_length == record.ref_allele.len()
                    })
                })
        {
            equal_floor_blocked_deletions.insert(identity.clone());
        }
        mixed_deletions_by_position
            .entry((record.chrom.clone(), record.pos))
            .or_default()
            .push((identity, record.ref_allele.len() - 1));
    }

    if record.alt_allele.split(',').any(|alternate| {
        record.ref_allele.len() > alternate.len() && record.ref_allele.starts_with(alternate)
    }) {
        deletions_by_position
            .entry((record.chrom.clone(), record.pos))
            .or_default()
            .push((record.end_pos(), record.ref_allele.len()));
    }
}

pub(crate) fn run(args: ValidatedPreprocessArgs) -> Result<()> {
    run_with_optional_reference(args, None)
}

pub(crate) fn run_with_reference(
    args: ValidatedPreprocessArgs,
    reference_sequences: &std::collections::BTreeMap<String, String>,
) -> Result<()> {
    run_with_optional_reference(args, Some(reference_sequences))
}

fn run_with_optional_reference(
    args: ValidatedPreprocessArgs,
    reference_sequences: Option<&std::collections::BTreeMap<String, String>>,
) -> Result<()> {
    let mut args = args.into_inner();
    if args.version {
        println!("{}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    // The reference is an argument: fail before any output path is staged.
    resolve_reference(args.reference.as_deref())?;
    let output_path = preprocess_output_path(&args);
    let index_path = if output_path.extension().and_then(|value| value.to_str()) == Some("bcf") {
        output_path.with_extension("bcf.csi")
    } else if output_path.extension().and_then(|value| value.to_str()) == Some("gz") {
        PathBuf::from(format!("{}.tbi", output_path.display()))
    } else {
        PathBuf::new()
    };
    let inputs = preprocess_inputs(&args);
    let outputs = std::iter::once(output_path.clone())
        .chain((!index_path.as_os_str().is_empty()).then_some(index_path.clone()))
        .chain(args.logfile.as_ref().map(PathBuf::from))
        .collect::<Vec<_>>();
    let transaction = OutputTransaction::files(&inputs, &outputs)?;
    let staged_output = transaction.staged_file(&output_path)?.to_path_buf();
    args.output = staged_output.to_string_lossy().into_owned();
    if let Some(logfile) = args.logfile.as_mut() {
        *logfile = transaction
            .staged_file(Path::new(logfile))?
            .to_string_lossy()
            .into_owned();
    }
    let outcome = run_inner(args, reference_sequences).map_err(|error| {
        anyhow::anyhow!(
            "failed to produce preprocess output {}: {error:#}",
            output_path.display()
        )
    });
    if let Err(error) = outcome {
        if output_path.extension().and_then(|value| value.to_str()) == Some("vcf")
            && staged_output.is_file()
        {
            transaction.commit().with_context(|| {
                format!(
                    "failed to publish legacy unindexed VCF {}",
                    output_path.display()
                )
            })?;
        }
        return Err(error);
    }
    if !index_path.as_os_str().is_empty() {
        let produced_index =
            if output_path.extension().and_then(|value| value.to_str()) == Some("bcf") {
                staged_output.with_extension("bcf.csi")
            } else {
                PathBuf::from(format!("{}.tbi", staged_output.display()))
            };
        let planned_index = transaction.staged_file(&index_path)?;
        fs::rename(&produced_index, planned_index).with_context(|| {
            format!(
                "failed to stage index destination {} from {}",
                index_path.display(),
                produced_index.display()
            )
        })?;
    }
    transaction.commit()
}

fn preprocess_output_path(args: &PreprocessArgs) -> PathBuf {
    if args.bcf && !args.output.ends_with(".bcf") {
        PathBuf::from(format!("{}.bcf", args.output))
    } else {
        PathBuf::from(&args.output)
    }
}

fn preprocess_inputs(args: &PreprocessArgs) -> Vec<PathBuf> {
    std::iter::once(args.input.as_str())
        .chain(args.reference.as_deref())
        .chain(args.regions_bedfile.as_deref())
        .chain(args.targets_bedfile.as_deref())
        .map(PathBuf::from)
        .collect()
}

fn run_inner(
    args: PreprocessArgs,
    shared_reference_sequences: Option<&std::collections::BTreeMap<String, String>>,
) -> Result<()> {
    let phase_started = std::time::Instant::now();
    let mut logger = PreprocessLogger::new(&args)?;
    logger.info(&format!("Preprocessing {}", args.input))?;
    let output_path = preprocess_output_path(&args);
    require_output_parent(&output_path)?;
    let reference_path = resolve_reference(args.reference.as_deref())?;
    let input_path = Path::new(&args.input);
    let input = vcf::open_validated_vcf(input_path)?;
    let mut headers = input.headers().to_vec();
    require_vcf_sample(&headers)?;
    // The index is part of the legacy reference contract even when an outer
    // comparison has already loaded and shared the sequence bodies.
    let reference_index = fasta::read_index(&reference_path)?;
    let owned_reference_sequences;
    let reference_sequences = if let Some(reference_sequences) = shared_reference_sequences {
        reference_sequences
    } else {
        owned_reference_sequences = fasta::read_sequences(&reference_path)?;
        &owned_reference_sequences
    };
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
    let mut input_contigs = BTreeSet::new();
    let mut deletion_sets = DeletionSets::default();
    let mut mixed_deletions_by_position = HashMap::new();
    let mut deletions_by_position = HashMap::new();
    let mut substitution_end_by_position: HashMap<(String, usize), usize> = HashMap::new();
    let mut haploid_x = false;
    let mut diploid_x = false;
    for record in input {
        let record = record?;
        observe_following_spanning_deletions(
            record.raw(),
            &mut mixed_deletions_by_position,
            &mut deletions_by_position,
            &mut deletion_sets,
        );
        if record_is_substitution(record.raw()) {
            let raw = record.raw();
            let end = raw.pos + raw.ref_allele.len().max(1) - 1;
            substitution_end_by_position
                .entry((raw.chrom.clone(), raw.pos))
                .and_modify(|current| *current = (*current).max(end))
                .or_insert(end);
        }
        input_contigs.insert(record.raw().chrom.clone());
        if args.gender == PreprocessGender::Auto {
            observe_gender(record.raw(), &mut haploid_x, &mut diploid_x);
        }
    }
    report_phase("input_header_inspection", phase_started);
    let requested_fixchr = if args.no_fixchr {
        Some(false)
    } else {
        args.fixchr
    };
    let fixchr = resolve_fixchr(requested_fixchr, &reference_contigs, &input_contigs);
    let gender = if args.gender == PreprocessGender::Auto {
        if haploid_x && !diploid_x {
            PreprocessGender::Male
        } else {
            PreprocessGender::Female
        }
    } else {
        args.gender
    };
    let leftshift = args.leftshift && !args.no_leftshift;
    let decompose = args.decompose && !args.no_decompose;
    let normalization_enabled = leftshift || decompose || gender == PreprocessGender::Male;
    let effective_threads = effective_thread_count(args.threads);
    let phase_started = std::time::Instant::now();
    let (blocksplit_selection, prepared_records) = if normalization_enabled && effective_threads > 1
    {
        let (observations, prepared_records) = collect_blocksplit_observations(
            vcf::open_validated_vcf(input_path)?,
            BlocksplitObservationParams {
                args: &args,
                fixchr,
                normalization_enabled,
                somatic_mode,
                somatic_sample_names: somatic_sample_names.as_deref(),
                reference_sequences: &reference_sequences,
                regions: regions.as_deref(),
                targets: targets.as_deref(),
                locations: locations.as_deref(),
            },
        )?;
        (
            select_blocksplit_resets(
                &observations,
                args.window_size,
                LEGACY_MAX_BLOCKS.min(effective_threads.saturating_mul(4)),
                locations.as_deref(),
            ),
            Some(prepared_records),
        )
    } else {
        (BlocksplitSelection::default(), None)
    };
    report_phase("block_selection", phase_started);

    let job_count = blocksplit_selection.jobs.as_ref().map_or(1, Vec::len);
    let phase_started = std::time::Instant::now();
    let declared_contigs = declared_contig_order(&headers, fixchr);
    let mut output = PreprocessSpool::new(normalization_enabled, job_count, &declared_contigs)?;
    let ctx = NormalizeContext {
        args: &args,
        gender,
        leftshift,
        decompose,
        somatic_mode,
        string_format_fields: &string_format_fields,
        reference_sequences,
        deletions: &deletion_sets,
        substitution_end_by_position: &substitution_end_by_position,
    };
    if let Some(prepared_records) = prepared_records {
        process_blocksplit_jobs(
            prepared_records,
            &blocksplit_selection,
            effective_threads,
            &ctx,
            &mut output,
        )?;
    } else {
        for job_index in 0..job_count {
            let job = blocksplit_selection
                .jobs
                .as_ref()
                .map(|jobs| &jobs[job_index]);
            let mut normalized_seen = HashSet::new();
            // Each legacy partial-credit job owns an independent left-shift
            // boundary, so overlapping jobs must process their copies separately.
            let mut prev_end_by_chrom: std::collections::HashMap<String, ShiftFloors> =
                std::collections::HashMap::new();
            let mut prepared_record_index = 0usize;
            for record in vcf::open_validated_vcf(input_path)? {
                let mut record = record?.into_raw();
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
                    conform_record_reference(&mut record, &reference_sequences)?;
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
                        output.push(
                            vcf::ValidatedVcfRecord::try_from_raw(
                                record,
                                QueryProvenance::Unavailable,
                            )?,
                            job_index,
                        )?;
                        continue;
                    }
                    process_normalized_record(
                        record,
                        symbolic_deletion,
                        &mut prev_end_by_chrom,
                        &ctx,
                        |record| output.push(record, job_index),
                    )?;
                }
            }
        }
    }
    report_phase("normalization", phase_started);
    output.sort_emitted_contigs();

    if somatic_mode.is_some() {
        rewrite_header_for_single_sample(&mut headers);
    }

    if fixchr {
        ensure_emitted_contig_headers(&mut headers, &output.emitted_contigs);
    }

    if normalization_enabled {
        headers = canonicalize_legacy_headers(&headers);
    } else {
        ensure_pass_filter_header(&mut headers);
    }

    let output_count = output.serial;
    let phase_started = std::time::Instant::now();
    let records = LocationAggregatedRecords::new(output.finish()?, normalization_enabled);
    vcf::write_validated_vcf_iter(&output_path, &headers, records)?;
    report_phase("external_sort_and_output_publication", phase_started);
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
        output_count,
        output_path.display()
    ))?;
    Ok(())
}

fn report_phase(name: &str, started: std::time::Instant) {
    if std::env::var_os("HAP_RS_PROFILE").is_some() {
        eprintln!(
            "HAP_RS_PHASE name={name} elapsed_seconds={:.3}",
            started.elapsed().as_secs_f64()
        );
    }
}

fn process_blocksplit_jobs(
    prepared_records: PreparedRecordSpool,
    selection: &BlocksplitSelection,
    threads: usize,
    ctx: &NormalizeContext<'_>,
    output: &mut PreprocessSpool,
) -> Result<()> {
    let job_count = selection.jobs.as_ref().map_or(1, Vec::len);
    let prepared_spool_bytes = prepared_records.len()?;
    let prepared_record_count = prepared_records.record_count();
    for contig in prepared_records.contigs() {
        output.seed_contig_rank(contig);
    }
    let dispatched_record_reads = selection
        .jobs
        .as_ref()
        .map_or(prepared_record_count, |jobs| {
            jobs.iter().map(|job| job.included_indices.len()).sum()
        });

    let worker_count = threads.max(1).min(job_count);
    if std::env::var_os("HAP_RS_PROFILE").is_some() {
        eprintln!(
            "HAP_RS_PROFILE input_decodes=2 block_jobs={job_count} workers={worker_count} \
             prepared_records={prepared_record_count} prepared_spool_bytes={prepared_spool_bytes} \
             dispatched_record_reads={dispatched_record_reads} job_spool_bytes=0"
        );
    }
    let next_job = std::sync::atomic::AtomicUsize::new(0);
    let (sender, receiver) = std::sync::mpsc::sync_channel(1024);
    let mut first_error = None;
    std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(worker_count);
        for _ in 0..worker_count {
            let sender = sender.clone();
            let next_job = &next_job;
            let prepared_records = &prepared_records;
            let jobs = selection.jobs.as_deref();
            handles.push(scope.spawn(move || {
                loop {
                    let job_index = next_job.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if job_index >= job_count {
                        break;
                    }
                    let mut previous_ends: std::collections::HashMap<String, ShiftFloors> =
                        std::collections::HashMap::new();
                    let mut reader = match prepared_records.reader() {
                        Ok(reader) => reader,
                        Err(error) => {
                            let _ = sender.send((job_index, Err(error)));
                            break;
                        }
                    };
                    let indices: Box<dyn Iterator<Item = usize>> = match jobs {
                        Some(jobs) => Box::new(jobs[job_index].included_indices.iter().copied()),
                        None => Box::new(0..prepared_record_count),
                    };
                    for record_index in indices {
                        let result = (|| {
                            let (record, symbolic_deletion) = reader.read_at(record_index)?;
                            if jobs.is_some_and(|jobs| {
                                jobs[job_index].reset_before_indices.contains(&record_index)
                            }) {
                                previous_ends.remove(&record.chrom);
                            }
                            process_normalized_record(
                                record,
                                symbolic_deletion,
                                &mut previous_ends,
                                ctx,
                                |record| {
                                    sender.send((job_index, Ok(record))).map_err(|_| {
                                        anyhow::anyhow!("preprocess output receiver closed")
                                    })
                                },
                            )
                        })();
                        if let Err(error) = result {
                            let _ = sender.send((job_index, Err(error)));
                            break;
                        }
                    }
                }
            }));
        }
        drop(sender);
        while let Ok((job_index, record)) = receiver.recv() {
            match record {
                Ok(record) if first_error.is_none() => {
                    if let Err(error) = output.push(record, job_index) {
                        first_error = Some(error);
                    }
                }
                Ok(_) => {}
                Err(error) if first_error.is_none() => first_error = Some(error),
                Err(_) => {}
            }
        }
        for handle in handles {
            if handle.join().is_err() && first_error.is_none() {
                first_error = Some(anyhow::anyhow!("preprocess worker panicked"));
            }
        }
    });
    if let Some(error) = first_error {
        return Err(error);
    }
    Ok(())
}

fn process_normalized_record(
    mut record: crate::domain::RawVcfRecord,
    symbolic_deletion: Option<SymbolicDeletionMaterialization>,
    prev_end_by_chrom: &mut std::collections::HashMap<String, ShiftFloors>,
    ctx: &NormalizeContext<'_>,
    mut emit: impl FnMut(vcf::ValidatedVcfRecord) -> Result<()>,
) -> Result<()> {
    let &NormalizeContext {
        args,
        gender,
        leftshift,
        decompose,
        somatic_mode,
        string_format_fields,
        reference_sequences,
        deletions,
        substitution_end_by_position,
    } = ctx;
    let DeletionSets {
        following_spanning: following_spanning_deletions,
        equal_floor_blocked: equal_floor_blocked_deletions,
        released_following: released_following_deletions,
        mixed_deletion_merge_partners,
        ..
    } = deletions;
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
    // `variant_pipeline::primitive_split` implements the legacy
    // `VariantPrimitiveSplitter` (stage 5 of `VariantInput.cpp`) plus the
    // `aggregate_hetalt` re-merge that brings same-position primitives back
    // to a multi-allelic shape. Same-direction insertions (`T → TG,TTG`)
    // re-merge into a single record; mixed-direction or different-length
    // deletions fan out into per-primitive records anchored at their
    // canonical position. Insert the ADO field BEFORE primitive_split *and*
    // before any GT normalisation. Legacy computes `ad_other` (= ADO) from
    // the *original* GT+AD inside VariantReader, then preserves that value
    // across both the half-call AlleleSplitter and the per-allele
    // PrimitiveSplitter. Computing ADO after splitting OR after the
    // haploid → homalt expansion below would widen the "called" set if
    // computed after expansion (input `GT=1` ad=[1,23] becomes `1/1`
    // ad=[1,23] → ADO=AD[0]=1, the unused ref depth, which matches legacy).
    // Computing ADO from the original GT=1 first preserves that ref-depth
    // signal correctly.
    ensure_missing_ad(&mut record);
    insert_ado_format(&mut record);
    ensure_missing_dp(&mut record);

    // The legacy C++ Variant representation has MAX_GT=2. During active
    // preprocessing, wider calls are converted to no-calls before
    // VariantCallsOnly removes their now-uncalled record. This is
    // observable in hap.py's vcfeval handoff for triploid and tetraploid
    // query records.
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
        return Ok(());
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
    let record_pos = record.pos;
    let is_pure_insertion = record_is_pure_insertion(&record);
    let is_substitution = record_is_substitution(&record);
    let floors = prev_end_by_chrom.entry(record_chrom.clone()).or_default();
    // Advancing to a strictly later position commits the previous position's
    // non-insertion end into the barrier and clears the staged term.
    if record_pos > floors.group_pos {
        floors.barrier = floors.barrier.max(floors.pending);
        floors.group_pos = record_pos;
        floors.pending = 0;
    }
    let prev_end = floors.prev_end;
    // A substitution at this position floors a shifting record here (a SNP
    // blocks a colocated deletion). Substitutions do not left-shift themselves,
    // so they don't consult it. The lookup is pre-scanned, so it is independent
    // of whether the SNP or the deletion is listed first.
    let same_position_substitution = if is_substitution {
        0
    } else {
        substitution_end_by_position
            .get(&(record_chrom.clone(), record_pos))
            .copied()
            .unwrap_or(0)
    };
    let leftshift_floor = floors.barrier.max(same_position_substitution);
    let release_leftshift_floor = released_following_deletions.contains(&(
        record.chrom.clone(),
        record.pos,
        record.ref_allele.clone(),
        record.alt_allele.clone(),
    ));
    // Legacy always runs VariantAlleleSplitter, fanning a multi-allelic record
    // into one record per called ALT before the location aggregator. The
    // primitive splitter (decompose) already splits the indel/complex alleles
    // it recognizes, so only pre-split the multi-allelics it leaves whole: in
    // no-decompose mode that is every multi-allelic, in decompose mode just the
    // pure-SNP ones `record_needs_primitive_split` skips. Without this a comma
    // ALT reaches `aggregate_location_records_inner`, which refuses to merge it
    // with a colocated single, so a het-alt like `T>A,C` never forms and a
    // standalone `TGA>T,TG` never fans out to its two per-position deletions.
    //
    // ponytail: single-sample only. The split's per-position primitives are
    // re-merged by `aggregate_location_records_inner`, whose het-alt merge is
    // not lossless for multi-sample records (splitting a 2-sample somatic
    // `C>T,G` drops an allele). Legacy splits regardless of sample count; the
    // upgrade path is a multi-sample-correct aggregator re-merge, after which
    // this `samples.len() == 1` guard can drop.
    let is_multi_allelic = record.alt_allele.contains(',')
        && !record.alt_allele.split(',').any(is_symbolic_allele);
    let split_multi_allelic = is_multi_allelic
        && record.samples.len() == 1
        && (!decompose || !variant_pipeline::record_needs_primitive_split(&record));
    let source_records = if matches!(
        symbolic_deletion,
        Some(SymbolicDeletionMaterialization::LeadingAnchor)
    ) {
        split_called_alleles(&record)
    } else if split_multi_allelic {
        // The reversed het-of-alts orientation is a symbolic-deletion detail;
        // the plain allele split keeps each call's projected GT verbatim.
        split_called_alleles(&record)
            .into_iter()
            .map(|mut materialized| {
                materialized.reverse_hetalt_samples = Vec::new();
                materialized
            })
            .collect()
    } else {
        vec![MaterializedAlleleRecord {
            record,
            reverse_hetalt_samples: Vec::new(),
        }]
    };
    let mut emitted_groups = Vec::with_capacity(source_records.len());
    for source in source_records {
        let reverse_hetalt_samples = source.reverse_hetalt_samples;
        let mut source = source.record;
        let is_mixed_deletion_merge_partner = mixed_deletion_merge_partners.contains(&(
            source.chrom.clone(),
            source.pos,
            source.ref_allele.clone(),
            source.alt_allele.clone(),
        ));
        if is_mixed_deletion_merge_partner {
            source.mixed_edit_primitive = true;
            assign_passthrough_primitive_identity(&mut source);
        }
        let mut emitted = if is_mixed_deletion_merge_partner {
            vec![source]
        } else if decompose && !source.alt_allele.split(',').any(is_symbolic_allele) {
            if let Some(reference) = reference_sequences.get(&source.chrom) {
                variant_pipeline::primitive_split_with_context(
                    &source,
                    reference.as_bytes(),
                    prev_end,
                    following_spanning_deletions.contains(&(
                        source.chrom.clone(),
                        source.pos,
                        source.ref_allele.clone(),
                        source.alt_allele.clone(),
                    )),
                    equal_floor_blocked_deletions.contains(&(
                        source.chrom.clone(),
                        source.pos,
                        source.ref_allele.clone(),
                        source.alt_allele.clone(),
                    )),
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
                if release_leftshift_floor {
                    extend_record_left(&mut split, reference.as_bytes());
                } else {
                    apply_left_shift(&mut split, reference.as_bytes(), leftshift_floor);
                }
            }
            canonicalize_multi_allelic_order(&mut split);
            canonicalize_legacy_genotypes(&mut split);
            if split.qual.is_empty() || split.qual == "." {
                split.qual = "0".to_string();
            }
            blank_secondary_sample_annotations(&mut split, args.bcf, string_format_fields);
            // Canonical FORMAT ordering: GT → AD/ADO/DP (fixed) → other
            // integer-typed fields alphabetical → float-typed alphabetical →
            // string-typed alphabetical. Mirrors the `int_fmts/float_fmts/
            // string_fmts` loop order in `VariantWriter.cpp` combined with
            // dynamic per-value type detection.
            reorder_format_fields(&mut split);
            emit(vcf::ValidatedVcfRecord::try_from_raw(
                split,
                QueryProvenance::Unavailable,
            )?)?;
        }
    }
    // Advance the per-chromosome boundaries. `prev_end` tracks every record's
    // span for the primitive splitter. A non-insertion contributes its end to
    // the current position group, committed into `barrier` at the next
    // position; same-position substitution floors come from the pre-scan.
    let floors = prev_end_by_chrom.entry(record_chrom).or_default();
    floors.prev_end = floors.prev_end.max(orig_end);
    if !is_pure_insertion {
        floors.pending = floors.pending.max(orig_end);
    }
    Ok(())
}

fn extend_record_left(record: &mut crate::domain::RawVcfRecord, reference: &[u8]) {
    if record.pos <= 1 || record.pos > reference.len() {
        return;
    }
    let anchor = reference[record.pos - 2].to_ascii_uppercase() as char;
    record.pos -= 1;
    record.ref_allele.insert(0, anchor);
    record.alt_allele = record
        .alt_allele
        .split(',')
        .map(|alternate| format!("{anchor}{alternate}"))
        .collect::<Vec<_>>()
        .join(",");
}

#[cfg(test)]
mod test_suite;
