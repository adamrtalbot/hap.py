use crate::application::{PreprocessArgs, PreprocessGender, ValidatedPreprocessArgs};
use crate::domain::QueryProvenance;
use crate::{
    adapters::{fasta, vcf},
    engines::variant_pipeline,
    output::OutputTransaction,
};
use anyhow::{Context, Result, bail};
use std::collections::{BTreeSet, HashSet};
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
use streaming::PreprocessSpool;

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

pub(crate) fn run(args: ValidatedPreprocessArgs) -> Result<()> {
    let mut args = args.into_inner();
    if args.version {
        println!("{}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    if args.reference.is_none() {
        args.reference = Some(resolve_reference(None)?.to_string_lossy().into_owned());
    }
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
    let outcome = run_inner(args).map_err(|error| {
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

fn run_inner(args: PreprocessArgs) -> Result<()> {
    let mut logger = PreprocessLogger::new(&args)?;
    logger.info(&format!("Preprocessing {}", args.input))?;
    let output_path = preprocess_output_path(&args);
    require_output_parent(&output_path)?;
    let reference_path = resolve_reference(args.reference.as_deref())?;
    let input_path = Path::new(&args.input);
    let input = vcf::open_validated_vcf(input_path)?;
    let mut headers = input.headers().to_vec();
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
    let mut input_contigs = BTreeSet::new();
    let mut haploid_x = false;
    let mut diploid_x = false;
    for record in input {
        let record = record?;
        input_contigs.insert(record.chrom.clone());
        if args.gender == PreprocessGender::Auto {
            observe_gender(&record, &mut haploid_x, &mut diploid_x);
        }
    }
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
    let blocksplit_selection = if normalization_enabled && effective_threads > 1 {
        let observations = collect_blocksplit_observations(
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
        select_blocksplit_resets(
            &observations,
            args.window_size,
            LEGACY_MAX_BLOCKS.min(effective_threads.saturating_mul(4)),
            locations.as_deref(),
        )
    } else {
        BlocksplitSelection::default()
    };

    let mut output = PreprocessSpool::new(normalization_enabled)?;
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
        for record in vcf::open_validated_vcf(input_path)? {
            let mut record = record?.raw().clone();
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
                    output.push(vcf::ValidatedVcfRecord::try_from_raw(
                        record,
                        QueryProvenance::Unavailable,
                    )?)?;
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
                        output.push(vcf::ValidatedVcfRecord::try_from_raw(
                            split,
                            QueryProvenance::Unavailable,
                        )?)?;
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
        ensure_emitted_contig_headers(&mut headers, &output.emitted_contigs);
    }

    if normalization_enabled {
        headers = canonicalize_legacy_headers(&headers);
    } else {
        ensure_pass_filter_header(&mut headers);
    }

    let output_count = output.serial;
    vcf::write_validated_vcf_iter(&output_path, &headers, output.finish()?)?;
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

#[cfg(test)]
mod test_suite;
