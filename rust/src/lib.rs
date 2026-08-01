pub mod align;
// Indexed BCF inspection is consumed by the feature-gated parity verifier;
// the product BCF reader/writer remains used without that feature.
#[cfg_attr(not(feature = "verification"), allow(dead_code))]
pub(crate) mod bcf;
pub mod cephes;
pub mod cli;
pub mod compare;
pub mod fasta;
#[cfg(feature = "verification")]
pub mod fixtures;
pub mod ftx;
pub mod metrics_json;
#[cfg(feature = "verification")]
pub mod parity_verifier;
pub mod partial_credit;
pub mod preprocess;
pub mod quantify;
pub mod report;
pub mod roc;
pub mod scmp;
pub mod somatic;
pub mod strelka;
pub mod validate;
pub mod variant_pipeline;
pub mod vcf;
#[cfg(feature = "verification")]
pub mod verification;
