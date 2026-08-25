//! Command orchestration: validate options, invoke adapters and engines, publish outputs.

pub(crate) mod compare;
pub(crate) mod comparison_io;
pub(crate) mod ftx;
pub(crate) mod preprocess;
pub(crate) mod quantify;
mod requests;
pub(crate) mod roc_publication;
pub(crate) mod somatic;
pub(crate) mod validate;

pub(crate) use requests::{
    CompareArgs, CompareEngine, EngineOptions, FtxArgs, PreprocessArgs, PreprocessGender,
    PreprocessOptions, QuantifyArgs, RequestValidationError, RocOptions, SomaticArgs,
    SomaticGtMode, ValidateArgs, ValidatedCompareArgs, ValidatedFtxArgs, ValidatedPreprocessArgs,
    ValidatedQuantifyArgs, ValidatedSomaticArgs, ValidatedValidateArgs,
};
