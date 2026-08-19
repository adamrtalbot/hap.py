//! Framework-independent application requests.
//!
//! Command-line and file-format adapters translate their inputs into these
//! types before invoking a use case. Validation here is deliberately limited
//! to invariants that do not require opening input files.

use std::error::Error;
use std::fmt;
use std::ops::Deref;
use std::path::Path;

use crate::domain::{OutputPlan, VariantOutputFormat};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct RequestValidationError {
    field: &'static str,
    message: String,
}

impl RequestValidationError {
    fn new(field: &'static str, message: impl Into<String>) -> Self {
        Self {
            field,
            message: message.into(),
        }
    }

    #[cfg(test)]
    pub(crate) fn field(&self) -> &'static str {
        self.field
    }
}

impl fmt::Display for RequestValidationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let field = match self.field {
            "ci_alpha" => "ci-alpha",
            field => field,
        };
        write!(formatter, "invalid {field}: {}", self.message)
    }
}

impl Error for RequestValidationError {}

pub(crate) type ValidationResult = Result<(), RequestValidationError>;

#[derive(Debug, Clone)]
pub(crate) struct Validated<T>(T);

impl<T> Deref for Validated<T> {
    type Target = T;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> Validated<T> {
    /// Transfers a validated request into its owning application use case.
    pub(crate) fn into_inner(self) -> T {
        self.0
    }
}

fn require_text(value: &str, field: &'static str) -> ValidationResult {
    if value.trim().is_empty() {
        return Err(RequestValidationError::new(field, "must not be empty"));
    }
    Ok(())
}

fn require_optional_text(value: Option<&str>, field: &'static str) -> ValidationResult {
    if let Some(value) = value {
        require_text(value, field)?;
    }
    Ok(())
}

fn require_texts(values: &[String], field: &'static str) -> ValidationResult {
    if values.iter().any(|value| value.trim().is_empty()) {
        return Err(RequestValidationError::new(
            field,
            "entries must not be empty",
        ));
    }
    Ok(())
}

fn require_threads(threads: Option<usize>) -> ValidationResult {
    if threads == Some(0) {
        return Err(RequestValidationError::new(
            "threads",
            "must be greater than zero",
        ));
    }
    Ok(())
}

fn require_roc_delta(value: f64) -> ValidationResult {
    if !value.is_finite() || value < 0.0 {
        return Err(RequestValidationError::new(
            "roc_delta",
            "must be finite and nonnegative",
        ));
    }
    Ok(())
}

fn require_ci_alpha(value: f64) -> ValidationResult {
    if !value.is_finite() || (value != 0.0 && !(0.0 < value && value < 1.0)) {
        return Err(RequestValidationError::new(
            "ci_alpha",
            "must be 0 (disabled) or strictly between 0 and 1",
        ));
    }
    Ok(())
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub(crate) enum CompareEngine {
    #[default]
    Xcmp,
    Vcfeval,
    ScmpSomatic,
    ScmpDistance,
}

impl CompareEngine {
    pub(crate) fn legacy_name(self) -> &'static str {
        match self {
            Self::Xcmp => "xcmp",
            Self::Vcfeval => "vcfeval",
            Self::ScmpSomatic => "scmp-somatic",
            Self::ScmpDistance => "scmp-distance",
        }
    }
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub(crate) enum SomaticGtMode {
    Half,
    Hemi,
    Het,
    Hom,
    First,
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
pub(crate) enum PreprocessGender {
    Male,
    Female,
    #[default]
    Auto,
    None,
}

#[derive(Debug, Clone)]
pub(crate) struct CompareArgs {
    pub(crate) truth: String,
    pub(crate) query: String,
    pub(crate) reference: String,
    pub(crate) report_prefix: String,
    pub(crate) version: bool,
    pub(crate) annotation_type: Option<String>,
    pub(crate) pass_only: bool,
    pub(crate) preprocess_truth: bool,
    pub(crate) convert_gvcf_truth: bool,
    pub(crate) convert_gvcf_query: bool,
    pub(crate) convert_gvcf_to_vcf: bool,
    pub(crate) usefiltered_truth: bool,
    pub(crate) filters_only: Option<String>,
    pub(crate) preprocess_window: usize,
    pub(crate) adjust_conf_regions: bool,
    pub(crate) no_adjust_conf_regions: bool,
    pub(crate) leftshift: bool,
    pub(crate) no_leftshift: bool,
    pub(crate) decompose: bool,
    pub(crate) no_decompose: bool,
    pub(crate) bcftools_norm: bool,
    pub(crate) fixchr: Option<bool>,
    pub(crate) no_fixchr: bool,
    pub(crate) filter_nonref: bool,
    pub(crate) somatic: bool,
    pub(crate) set_gt: Option<SomaticGtMode>,
    pub(crate) gender: PreprocessGender,
    pub(crate) bcf: bool,
    pub(crate) regions_bedfile: Option<String>,
    pub(crate) targets_bedfile: Option<String>,
    pub(crate) fp_bedfile: Option<String>,
    pub(crate) locations: Option<String>,
    pub(crate) threads: Option<usize>,
    pub(crate) strat_tsv: Option<String>,
    pub(crate) strat_regions: Vec<String>,
    pub(crate) strat_fixchr: bool,
    #[allow(dead_code, reason = "retained in the raw DTO for legacy CLI parity")]
    pub(crate) write_vcf: bool,
    pub(crate) write_counts: bool,
    pub(crate) no_write_counts: bool,
    pub(crate) output_vtc: bool,
    pub(crate) preserve_info: bool,
    pub(crate) roc: String,
    pub(crate) no_roc: bool,
    pub(crate) roc_regions: Vec<String>,
    pub(crate) roc_filter: Option<String>,
    pub(crate) roc_delta: f64,
    pub(crate) ci_alpha: f64,
    pub(crate) no_json: bool,
    pub(crate) no_hc: bool,
    pub(crate) window: usize,
    pub(crate) max_enum: usize,
    pub(crate) hb_expand: usize,
    pub(crate) engine: CompareEngine,
    pub(crate) engine_vcfeval: Option<String>,
    pub(crate) engine_vcfeval_template: Option<String>,
    pub(crate) engine_scmp_distance: usize,
    #[allow(dead_code, reason = "retained in the raw DTO for legacy CLI parity")]
    pub(crate) force_interactive: bool,
    pub(crate) scratch_prefix: Option<String>,
    pub(crate) keep_scratch: bool,
    pub(crate) logfile: Option<String>,
    pub(crate) verbose: bool,
    pub(crate) quiet: bool,
}

impl CompareArgs {
    pub(crate) fn validate(&self) -> ValidationResult {
        require_text(&self.truth, "truth")?;
        require_text(&self.query, "query")?;
        require_text(&self.report_prefix, "report_prefix")?;
        OutputPlan::new(
            &self.report_prefix,
            self.write_counts && !self.no_write_counts,
            !self.no_json,
            Some(if self.bcf {
                VariantOutputFormat::Bcf
            } else {
                VariantOutputFormat::Vcf
            }),
        )
        .map_err(|error| RequestValidationError::new("report_prefix", error.to_string()))?;
        require_optional_text(self.regions_bedfile.as_deref(), "regions_bedfile")?;
        require_optional_text(self.targets_bedfile.as_deref(), "targets_bedfile")?;
        require_optional_text(self.fp_bedfile.as_deref(), "fp_bedfile")?;
        require_optional_text(self.strat_tsv.as_deref(), "strat_tsv")?;
        require_optional_text(self.scratch_prefix.as_deref(), "scratch_prefix")?;
        require_optional_text(self.logfile.as_deref(), "logfile")?;
        require_texts(&self.strat_regions, "strat_regions")?;
        require_texts(&self.roc_regions, "roc_regions")?;
        require_threads(self.threads)?;
        require_roc_delta(self.roc_delta)?;
        require_ci_alpha(self.ci_alpha)?;
        if let Some(annotation_type) = self.annotation_type.as_deref() {
            match annotation_type {
                "xcmp" | "ga4gh" => {}
                other => {
                    return Err(RequestValidationError::new(
                        "annotation_type",
                        format!("unsupported value '{other}'"),
                    ));
                }
            }
        }
        require_text(&self.roc, "roc")?;
        if self.engine == CompareEngine::ScmpDistance
            && self.engine_scmp_distance > i64::MAX as usize
        {
            return Err(RequestValidationError::new(
                "engine_scmp_distance",
                "exceeds the supported signed distance range",
            ));
        }
        Ok(())
    }

    pub(crate) fn validated(self) -> Result<ValidatedCompareArgs, RequestValidationError> {
        self.validate()?;
        Ok(ValidatedCompareArgs { values: self })
    }

    #[cfg(test)]
    pub(crate) fn with_paths(
        truth: String,
        query: String,
        reference: String,
        report_prefix: String,
    ) -> Self {
        Self {
            truth,
            query,
            reference,
            report_prefix,
            version: false,
            annotation_type: None,
            pass_only: false,
            preprocess_truth: false,
            convert_gvcf_truth: false,
            convert_gvcf_query: false,
            convert_gvcf_to_vcf: false,
            usefiltered_truth: false,
            filters_only: None,
            preprocess_window: 10_000,
            adjust_conf_regions: true,
            no_adjust_conf_regions: false,
            leftshift: false,
            no_leftshift: false,
            decompose: false,
            no_decompose: false,
            bcftools_norm: false,
            fixchr: None,
            no_fixchr: false,
            filter_nonref: false,
            somatic: false,
            set_gt: None,
            gender: PreprocessGender::Auto,
            bcf: false,
            regions_bedfile: None,
            targets_bedfile: None,
            fp_bedfile: None,
            locations: None,
            threads: None,
            strat_tsv: None,
            strat_regions: Vec::new(),
            strat_fixchr: false,
            write_vcf: false,
            write_counts: true,
            no_write_counts: false,
            output_vtc: false,
            preserve_info: false,
            roc: "QUAL".to_string(),
            no_roc: false,
            roc_regions: Vec::new(),
            roc_filter: None,
            roc_delta: 0.5,
            ci_alpha: 0.0,
            no_json: false,
            no_hc: false,
            window: 50,
            max_enum: 16_768,
            hb_expand: 30,
            engine: CompareEngine::Xcmp,
            engine_vcfeval: None,
            engine_vcfeval_template: None,
            engine_scmp_distance: 30,
            force_interactive: false,
            scratch_prefix: None,
            keep_scratch: false,
            logfile: None,
            verbose: false,
            quiet: false,
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct ValidatedCompareArgs {
    values: CompareArgs,
}

impl ValidatedCompareArgs {
    pub(crate) fn into_inner(self) -> CompareArgs {
        self.values
    }
}

impl Deref for ValidatedCompareArgs {
    type Target = CompareArgs;

    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

#[derive(Debug, Clone)]
pub(crate) struct PreprocessArgs {
    pub(crate) input: String,
    pub(crate) output: String,
    pub(crate) version: bool,
    pub(crate) reference: Option<String>,
    pub(crate) locations: Option<String>,
    pub(crate) pass_only: bool,
    pub(crate) filters_only: Option<String>,
    pub(crate) regions_bedfile: Option<String>,
    pub(crate) targets_bedfile: Option<String>,
    pub(crate) fixchr: Option<bool>,
    pub(crate) no_fixchr: bool,
    pub(crate) somatic: bool,
    pub(crate) set_gt: Option<SomaticGtMode>,
    pub(crate) filter_nonref: bool,
    pub(crate) convert_gvcf_to_vcf: bool,
    pub(crate) bcf: bool,
    pub(crate) bcftools_norm: bool,
    pub(crate) leftshift: bool,
    pub(crate) no_leftshift: bool,
    pub(crate) decompose: bool,
    pub(crate) no_decompose: bool,
    pub(crate) gender: PreprocessGender,
    pub(crate) window_size: i64,
    pub(crate) threads: Option<usize>,
    pub(crate) logfile: Option<String>,
    pub(crate) verbose: bool,
    #[allow(dead_code, reason = "retained in the raw DTO for legacy CLI parity")]
    pub(crate) quiet: bool,
    #[allow(dead_code, reason = "retained in the raw DTO for legacy CLI parity")]
    pub(crate) force_interactive: bool,
}

impl PreprocessArgs {
    fn validate_common(&self) -> ValidationResult {
        require_text(&self.input, "input")?;
        require_text(&self.output, "output")?;
        require_optional_text(self.reference.as_deref(), "reference")?;
        require_optional_text(self.regions_bedfile.as_deref(), "regions_bedfile")?;
        require_optional_text(self.targets_bedfile.as_deref(), "targets_bedfile")?;
        require_optional_text(self.logfile.as_deref(), "logfile")?;
        require_threads(self.threads)?;
        Ok(())
    }

    pub(crate) fn validate(&self) -> ValidationResult {
        self.validate_common()?;
        if !self.bcf
            && Path::new(&self.output)
                .extension()
                .and_then(|value| value.to_str())
                == Some("vcf")
        {
            return Err(RequestValidationError::new(
                "output",
                "plain VCF output cannot be indexed; use compressed VCF or BCF output",
            ));
        }
        Ok(())
    }

    pub(crate) fn validated(self) -> Result<ValidatedPreprocessArgs, RequestValidationError> {
        self.validate()?;
        Ok(Validated(self))
    }

    #[cfg(test)]
    pub(crate) fn validated_with_legacy_plain_vcf(
        self,
    ) -> Result<ValidatedPreprocessArgs, RequestValidationError> {
        self.validate_common()?;
        Ok(Validated(self))
    }
}

#[derive(Debug, Clone)]
pub(crate) struct QuantifyArgs {
    pub(crate) input_vcf: String,
    pub(crate) report_prefix: String,
    pub(crate) reference: String,
    pub(crate) annotation_type: Option<String>,
    pub(crate) fp_bedfile: Option<String>,
    pub(crate) strat_tsv: Option<String>,
    pub(crate) strat_regions: Vec<String>,
    pub(crate) strat_fixchr: bool,
    pub(crate) write_vcf: bool,
    pub(crate) write_counts: bool,
    pub(crate) output_vtc: bool,
    pub(crate) preserve_info: bool,
    pub(crate) adjust_conf_regions: Option<String>,
    pub(crate) threads: Option<usize>,
    pub(crate) bcf: bool,
    pub(crate) logfile: Option<String>,
    pub(crate) verbose: bool,
    #[allow(dead_code, reason = "retained in the raw DTO for legacy CLI parity")]
    pub(crate) quiet: bool,
    #[allow(dead_code, reason = "retained in the raw DTO for legacy CLI parity")]
    pub(crate) force_interactive: bool,
    pub(crate) roc: String,
    pub(crate) do_roc: bool,
    pub(crate) roc_regions: Vec<String>,
    pub(crate) roc_filter: Option<String>,
    pub(crate) roc_delta: f64,
    pub(crate) ci_alpha: f64,
    pub(crate) no_json: bool,
}

impl QuantifyArgs {
    pub(crate) fn validate(&self) -> ValidationResult {
        require_text(&self.input_vcf, "input_vcf")?;
        self.validate_options()
    }

    fn validate_options(&self) -> ValidationResult {
        require_text(&self.report_prefix, "report_prefix")?;
        require_text(&self.reference, "reference")?;
        let output_plan = OutputPlan::new(
            &self.report_prefix,
            self.write_counts,
            !self.no_json,
            self.write_vcf.then_some(if self.bcf {
                VariantOutputFormat::Bcf
            } else {
                VariantOutputFormat::Vcf
            }),
        )
        .map_err(|error| RequestValidationError::new("report_prefix", error.to_string()))?;
        if !self.input_vcf.is_empty()
            && output_plan.conflicts_with_input(Path::new(&self.input_vcf))
        {
            return Err(RequestValidationError::new(
                "report_prefix",
                "cannot overwrite input VCF",
            ));
        }
        require_optional_text(self.fp_bedfile.as_deref(), "fp_bedfile")?;
        require_optional_text(self.strat_tsv.as_deref(), "strat_tsv")?;
        require_optional_text(self.logfile.as_deref(), "logfile")?;
        require_texts(&self.strat_regions, "strat_regions")?;
        require_texts(&self.roc_regions, "roc_regions")?;
        require_threads(self.threads)?;
        require_roc_delta(self.roc_delta)?;
        require_ci_alpha(self.ci_alpha)?;
        match self.annotation_type.as_deref().unwrap_or("xcmp") {
            "xcmp" | "ga4gh" => {}
            other => {
                return Err(RequestValidationError::new(
                    "annotation_type",
                    format!("unsupported value '{other}'"),
                ));
            }
        }
        require_text(&self.roc, "roc")?;
        if self.adjust_conf_regions.is_some() && self.fp_bedfile.is_none() {
            return Err(RequestValidationError::new(
                "adjust_conf_regions",
                "requires fp_bedfile",
            ));
        }
        Ok(())
    }

    pub(crate) fn validated(self) -> Result<ValidatedQuantifyArgs, RequestValidationError> {
        self.validate()?;
        Ok(ValidatedQuantifyArgs { values: self })
    }
}

#[derive(Debug, Clone)]
pub(crate) struct SomaticArgs {
    pub(crate) truth: String,
    pub(crate) query: String,
    pub(crate) output: String,
    pub(crate) reference: Option<String>,
    pub(crate) location: Option<String>,
    pub(crate) regions_bedfile: Option<String>,
    pub(crate) targets_bedfile: Option<String>,
    pub(crate) fp_bedfile: Option<String>,
    pub(crate) ambiguous_beds: Vec<String>,
    pub(crate) ambi_fp: bool,
    pub(crate) no_ambi_fp: bool,
    pub(crate) count_unk: bool,
    pub(crate) no_count_unk: bool,
    pub(crate) explain_ambiguous: bool,
    pub(crate) include_nonpass: bool,
    pub(crate) fp_region_size: Option<String>,
    pub(crate) feature_table: Option<String>,
    pub(crate) happy_stats: bool,
    pub(crate) bams: Vec<String>,
    pub(crate) normalize_truth: bool,
    pub(crate) normalize_query: bool,
    pub(crate) normalize_all: bool,
    pub(crate) fixchr_truth: Option<bool>,
    pub(crate) fixchr_query: Option<bool>,
    pub(crate) no_fixchr_truth: bool,
    pub(crate) no_fixchr_query: bool,
    pub(crate) no_order_check: bool,
    pub(crate) roc: Option<String>,
    pub(crate) af_strat: bool,
    pub(crate) af_strat_binsize: String,
    pub(crate) af_strat_truth: String,
    pub(crate) af_strat_query: String,
    pub(crate) count_filtered_fn: bool,
    pub(crate) ci_level: f64,
    pub(crate) scratch_prefix: Option<String>,
    pub(crate) keep_scratch: bool,
    pub(crate) cont: bool,
    pub(crate) logfile: Option<String>,
    pub(crate) verbose: bool,
    pub(crate) quiet: bool,
}

impl SomaticArgs {
    pub(crate) fn validate(&self) -> ValidationResult {
        require_text(&self.truth, "truth")?;
        require_text(&self.query, "query")?;
        require_text(&self.output, "output")?;
        require_optional_text(self.reference.as_deref(), "reference")?;
        require_optional_text(self.regions_bedfile.as_deref(), "regions_bedfile")?;
        require_optional_text(self.targets_bedfile.as_deref(), "targets_bedfile")?;
        require_optional_text(self.fp_bedfile.as_deref(), "fp_bedfile")?;
        require_optional_text(self.feature_table.as_deref(), "feature_table")?;
        require_optional_text(self.scratch_prefix.as_deref(), "scratch_prefix")?;
        require_optional_text(self.logfile.as_deref(), "logfile")?;
        require_texts(&self.ambiguous_beds, "ambiguous_beds")?;
        require_texts(&self.bams, "bams")?;
        if !(self.ci_level.is_finite() && 0.0 < self.ci_level && self.ci_level < 1.0) {
            return Err(RequestValidationError::new(
                "ci_level",
                "must be finite and strictly between 0 and 1",
            ));
        }
        let has_effective_feature_table = self.feature_table.is_some() || self.roc.is_some();
        if self.af_strat && !has_effective_feature_table {
            return Err(RequestValidationError::new(
                "af_strat",
                "requires feature_table",
            ));
        }
        if self.af_strat {
            require_text(&self.af_strat_truth, "af_strat_truth")?;
            require_text(&self.af_strat_query, "af_strat_query")?;
        }
        if self.count_filtered_fn && (!self.include_nonpass || !has_effective_feature_table) {
            return Err(RequestValidationError::new(
                "count_filtered_fn",
                "requires include_nonpass and feature_table",
            ));
        }
        if self.happy_stats && (!self.include_nonpass || !has_effective_feature_table) {
            return Err(RequestValidationError::new(
                "happy_stats",
                "requires include_nonpass and feature_table",
            ));
        }
        validate_af_bins(&self.af_strat_binsize)
    }

    pub(crate) fn validated(self) -> Result<ValidatedSomaticArgs, RequestValidationError> {
        self.validate()?;
        Ok(Validated(self))
    }
}

fn validate_af_bins(raw: &str) -> ValidationResult {
    require_text(raw, "af_strat_binsize")?;
    let values = raw
        .split(',')
        .map(|value| {
            if value.trim().is_empty() {
                return Err(RequestValidationError::new(
                    "af_strat_binsize",
                    "entries must not be empty",
                ));
            }
            value.parse::<f64>().map_err(|_| {
                RequestValidationError::new(
                    "af_strat_binsize",
                    format!("'{value}' is not a number"),
                )
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut start = 0.0_f64;
    let mut index = 0usize;
    for _ in 0..100 {
        if !start.is_finite() || start >= 1.0 {
            return Ok(());
        }
        let mut end = start + values[index];
        if end >= 1.0 {
            end = 1.000_000_01;
        }
        if start >= end {
            return Ok(());
        }
        start = end;
        index = (index + 1) % values.len();
    }
    Err(RequestValidationError::new(
        "af_strat_binsize",
        "values produce more than 100 bins",
    ))
}

#[derive(Debug, Clone)]
pub(crate) struct FtxArgs {
    pub(crate) input: String,
    pub(crate) output: String,
    pub(crate) location: Option<String>,
    pub(crate) regions_bedfile: Option<String>,
    pub(crate) targets_bedfile: Option<String>,
    pub(crate) include_nonpass: bool,
    pub(crate) features: String,
    pub(crate) label: Option<String>,
    pub(crate) bams: Vec<String>,
    pub(crate) reference: Option<String>,
    pub(crate) normalize: bool,
    pub(crate) fixchr: bool,
}

impl FtxArgs {
    pub(crate) fn validate(&self) -> ValidationResult {
        require_text(&self.input, "input")?;
        require_text(&self.output, "output")?;
        require_text(&self.features, "features")?;
        require_optional_text(self.regions_bedfile.as_deref(), "regions_bedfile")?;
        require_optional_text(self.targets_bedfile.as_deref(), "targets_bedfile")?;
        require_optional_text(self.reference.as_deref(), "reference")?;
        require_texts(&self.bams, "bams")?;
        Ok(())
    }

    pub(crate) fn validated(self) -> Result<ValidatedFtxArgs, RequestValidationError> {
        self.validate()?;
        Ok(Validated(self))
    }
}

#[derive(Debug, Clone)]
pub(crate) struct ValidateArgs {
    pub(crate) input: String,
    pub(crate) reference: Option<String>,
    pub(crate) output_json: Option<String>,
    pub(crate) errors_bed: Option<String>,
    pub(crate) locations: Option<String>,
    pub(crate) regions_bedfile: Option<String>,
    pub(crate) targets_bedfile: Option<String>,
    pub(crate) apply_filters: bool,
    pub(crate) limit_records: Option<i64>,
    pub(crate) message_every: Option<i64>,
    pub(crate) strict_homref: bool,
    pub(crate) check_bcf_errors: bool,
    pub(crate) all_warnings: bool,
}

impl ValidateArgs {
    pub(crate) fn validate(&self) -> ValidationResult {
        require_text(&self.input, "input")?;
        require_optional_text(self.reference.as_deref(), "reference")?;
        require_optional_text(self.output_json.as_deref(), "output_json")?;
        require_optional_text(self.errors_bed.as_deref(), "errors_bed")?;
        require_optional_text(self.regions_bedfile.as_deref(), "regions_bedfile")?;
        require_optional_text(self.targets_bedfile.as_deref(), "targets_bedfile")?;
        Ok(())
    }

    pub(crate) fn validated(self) -> Result<ValidatedValidateArgs, RequestValidationError> {
        self.validate()?;
        Ok(Validated(self))
    }
}

pub(crate) type ValidatedPreprocessArgs = Validated<PreprocessArgs>;
pub(crate) type ValidatedSomaticArgs = Validated<SomaticArgs>;
pub(crate) type ValidatedFtxArgs = Validated<FtxArgs>;
pub(crate) type ValidatedValidateArgs = Validated<ValidateArgs>;

#[derive(Debug, Clone)]
pub(crate) struct ValidatedQuantifyArgs {
    values: QuantifyArgs,
}

impl ValidatedQuantifyArgs {
    pub(crate) fn into_inner(self) -> QuantifyArgs {
        self.values
    }
}

impl Deref for ValidatedQuantifyArgs {
    type Target = QuantifyArgs;

    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn compare() -> CompareArgs {
        CompareArgs::with_paths(
            "truth.vcf.gz".into(),
            "query.vcf.gz".into(),
            String::new(),
            "results/report".into(),
        )
    }

    fn validate_request() -> ValidateArgs {
        ValidateArgs {
            input: "input.vcf.gz".into(),
            reference: None,
            output_json: None,
            errors_bed: None,
            locations: None,
            regions_bedfile: None,
            targets_bedfile: None,
            apply_filters: false,
            limit_records: None,
            message_every: None,
            strict_homref: false,
            check_bcf_errors: false,
            all_warnings: false,
        }
    }

    #[test]
    fn compare_accepts_default_reference_sentinel() {
        assert!(compare().validate().is_ok());
    }

    #[test]
    fn compare_rejects_bad_roc_and_ci_values() {
        let mut request = compare();
        request.roc_delta = f64::NAN;
        assert_eq!(request.validate().unwrap_err().field(), "roc_delta");

        request.roc_delta = 0.5;
        request.ci_alpha = 1.0;
        assert_eq!(request.validate().unwrap_err().field(), "ci_alpha");
    }

    #[test]
    fn compare_rejects_zero_threads_and_empty_prefix() {
        let mut request = compare();
        request.threads = Some(0);
        assert_eq!(request.validate().unwrap_err().field(), "threads");

        request.threads = Some(1);
        request.report_prefix = "  ".into();
        assert_eq!(request.validate().unwrap_err().field(), "report_prefix");
    }

    #[test]
    fn preprocess_rejects_plain_vcf_output_before_io() {
        let request = PreprocessArgs {
            input: "input.vcf.gz".into(),
            output: "output.vcf".into(),
            version: false,
            reference: None,
            locations: None,
            pass_only: false,
            filters_only: None,
            regions_bedfile: None,
            targets_bedfile: None,
            fixchr: None,
            no_fixchr: false,
            somatic: false,
            set_gt: None,
            filter_nonref: true,
            convert_gvcf_to_vcf: false,
            bcf: false,
            bcftools_norm: false,
            leftshift: true,
            no_leftshift: false,
            decompose: true,
            no_decompose: false,
            gender: PreprocessGender::Auto,
            window_size: 10_000,
            threads: None,
            logfile: None,
            verbose: false,
            quiet: false,
            force_interactive: false,
        };
        assert_eq!(request.validate().unwrap_err().field(), "output");
    }

    #[test]
    fn validate_preserves_legacy_negative_controls() {
        let mut request = validate_request();
        request.limit_records = Some(-2);
        request.message_every = Some(-1);
        request
            .validate()
            .expect("legacy negative controls are interpreted by the use case");
    }

    #[test]
    fn validation_error_is_actionable() {
        let mut request = compare();
        request.truth.clear();
        assert_eq!(
            request.validate().unwrap_err().to_string(),
            "invalid truth: must not be empty"
        );
    }
}
