use anyhow::Result;
use clap::{Parser, Subcommand};
use hap_rs::parity_verifier::{
    AggregateReport, CompareOutputs, aggregate_report, compare_outputs, validate_roc,
};
use hap_rs::verification::{
    LEGACY_IMAGE, bless_legacy_outputs, verify_legacy_outputs, verify_realworld_outputs,
    verify_rust_outputs,
};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[command(
    name = "verify-fixtures",
    version,
    about = "Run fixture verification procedures"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    BlessLegacy {
        #[arg(long, default_value = LEGACY_IMAGE)]
        image: String,
    },
    CheckLegacy {
        #[arg(long, default_value = LEGACY_IMAGE)]
        image: String,
    },
    CheckRust,
    CheckRealworld {
        #[arg(long, default_value = LEGACY_IMAGE)]
        image: String,
    },
    /// Compare one legacy/Rust artifact pair and always publish a verdict.
    CompareOutputs {
        #[arg(long)]
        legacy_dir: PathBuf,
        #[arg(long)]
        rust_dir: PathBuf,
        #[arg(long)]
        prefix: String,
        #[arg(long = "expected-artifact", value_name = "NAME")]
        expected_artifact: Vec<String>,
        #[arg(long)]
        case: String,
        #[arg(long)]
        sample: String,
        #[arg(long)]
        report: PathBuf,
        #[arg(long)]
        status: PathBuf,
    },
    /// Aggregate per-case status JSON files into Markdown and CSV reports.
    AggregateReport {
        #[arg(long)]
        image: String,
        #[arg(long)]
        hap_bin: String,
        #[arg(long)]
        markdown: PathBuf,
        #[arg(long)]
        csv: PathBuf,
        inputs: Vec<PathBuf>,
    },
    /// Reject structurally corrupt legacy ROC output.
    ValidateRoc {
        roc: PathBuf,
    },
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::BlessLegacy { image } => bless_legacy_outputs(&image),
        Command::CheckLegacy { image } => verify_legacy_outputs(&image),
        Command::CheckRust => verify_rust_outputs(),
        Command::CheckRealworld { image } => verify_realworld_outputs(&image),
        Command::CompareOutputs {
            legacy_dir,
            rust_dir,
            prefix,
            expected_artifact,
            case,
            sample,
            report,
            status,
        } => compare_outputs(&CompareOutputs {
            legacy_dir: &legacy_dir,
            rust_dir: &rust_dir,
            prefix: &prefix,
            expected_artifacts: &expected_artifact,
            case: &case,
            sample: &sample,
            report: &report,
            status: &status,
        }),
        Command::AggregateReport {
            image,
            hap_bin,
            markdown,
            csv,
            inputs,
        } => aggregate_report(&AggregateReport {
            image: &image,
            hap_bin: &hap_bin,
            markdown: &markdown,
            csv: &csv,
            inputs: &inputs,
        }),
        Command::ValidateRoc { roc } => validate_roc(&roc).map(|_| ()),
    }
}
