//! Application core for the `hap` command-line program.
//!
//! The crate intentionally exposes only [`run`]. Parsing, compatibility
//! adapters, comparison engines, codecs, and report writers are implementation
//! details so they can evolve without becoming an accidental public API.

#![forbid(unsafe_code)]
#![warn(missing_docs, rustdoc::broken_intra_doc_links, unreachable_pub)]

mod adapters;
mod application;
mod cli_compat;
mod domain;
mod engines;
mod output;

/// Parser entry points used exclusively by the opt-in fuzz workspace.
///
/// This module is absent from normal builds, keeping the crate's supported
/// API limited to [`run`] while allowing cargo-fuzz to exercise the native
/// parsers without spawning subprocesses.
#[cfg(feature = "fuzzing")]
#[doc(hidden)]
pub mod fuzzing {
    use std::collections::BTreeSet;
    use std::io::Read;
    use std::path::Path;

    /// Exercise VCF record parsing for every non-header input line.
    pub fn vcf(data: &[u8]) {
        if let Ok(text) = std::str::from_utf8(data) {
            for line in text.lines().filter(|line| !line.starts_with('#')) {
                let _ = crate::domain::RawVcfRecord::from_line(line, Path::new("<fuzz-vcf>"));
            }
        }
    }

    /// Exercise native BCF decoding.
    pub fn bcf(data: &[u8]) {
        const VALID_BGZF_BCF: &[u8] =
            include_bytes!("../../verification/assets/fixtures/pre-matrix/input.bcf");
        if data == VALID_BGZF_BCF {
            let mut uncompressed = Vec::new();
            if noodles_bgzf::io::Reader::new(data)
                .read_to_end(&mut uncompressed)
                .is_ok()
            {
                let _ = crate::adapters::bcf::decode(&uncompressed, Path::new("<fuzz-bcf>"));
            }
        } else {
            let _ = crate::adapters::bcf::decode(data, Path::new("<fuzz-bcf>"));
        }
    }

    /// Exercise FASTA text parsing.
    pub fn fasta(data: &[u8]) {
        if let Ok(text) = std::str::from_utf8(data) {
            let _ = crate::adapters::fasta::parse_sequences(text, Path::new("<fuzz-fasta>"));
        }
    }

    /// Exercise BED text parsing.
    pub fn bed(data: &[u8]) {
        if let Ok(text) = std::str::from_utf8(data) {
            let contigs = BTreeSet::from(["1".to_string(), "chr1".to_string()]);
            let _ = crate::adapters::vcf::parse_bed(text, &contigs, Path::new("<fuzz-bed>"));
        }
    }

    /// Exercise comma-separated location parsing and matching.
    pub fn location(data: &[u8]) {
        if let Ok(text) = std::str::from_utf8(data) {
            let contigs = BTreeSet::from(["1".to_string(), "chr1".to_string()]);
            if let Ok(filters) = crate::adapters::vcf::parse_locations(text, &contigs) {
                for filter in filters {
                    let _ = filter.matches("chr1", 1);
                    let _ = filter.matches("1", usize::MAX);
                }
            }
        }
    }
}

use anyhow::Result;
use clap::{CommandFactory, Parser, error::ErrorKind};
use cli_compat::cli::{
    Cli, Command, process_args_with_legacy_somatic_aliases, requests_legacy_subcommand_version,
    requests_quantify_version, validate_legacy_germline_version_arguments,
};
use cli_compat::compatibility::{self, UsageErrorPolicy};

fn exit_with_legacy_help(
    arguments: &[std::ffi::OsString],
    error: clap::Error,
    policy: UsageErrorPolicy,
) -> ! {
    eprint!("{error}");
    if let Some(warning) = policy.warning() {
        eprintln!("{warning}");
    }
    let canonical = match arguments.get(1).and_then(|value| value.to_str()) {
        Some("germline" | "compare") => "germline",
        Some("pre" | "preprocess" | "prepy") => "pre",
        Some("quantify" | "qfy") => "quantify",
        _ => std::process::exit(policy.exit_code()),
    };
    let mut command = Cli::command();
    if let Some(subcommand) = command.find_subcommand_mut(canonical) {
        print!("{}", subcommand.render_help());
    }
    std::process::exit(policy.exit_code())
}

/// Parse the process arguments and run the selected `hap` command.
///
/// # Errors
///
/// Returns an error when command execution, input parsing, comparison, or
/// output publication fails. Clap usage errors retain their command-specific
/// process exit behavior for compatibility with the legacy tools.
pub fn run() -> Result<()> {
    let arguments = process_args_with_legacy_somatic_aliases();
    if requests_legacy_subcommand_version(&arguments) {
        if let Err(error) = validate_legacy_germline_version_arguments(&arguments) {
            if error.kind() == ErrorKind::UnknownArgument {
                exit_with_legacy_help(&arguments, error, UsageErrorPolicy::LegacyFailure);
            }
            error.exit();
        }
        println!("Hap.py ");
        return Ok(());
    }
    let quantify_version = requests_quantify_version(&arguments);
    let cli = match Cli::try_parse_from(arguments.clone()) {
        Ok(cli) => cli,
        Err(error) if error.kind() == ErrorKind::UnknownArgument => {
            if let Some(policy) = compatibility::usage_error_policy(&arguments) {
                exit_with_legacy_help(&arguments, error, policy);
            }
            error.exit();
        }
        Err(error)
            if error.kind() == ErrorKind::MissingRequiredArgument
                && compatibility::usage_error_policy(&arguments)
                    == Some(UsageErrorPolicy::LegacyFailure) =>
        {
            eprint!("{error}");
            std::process::exit(1);
        }
        Err(error) => error.exit(),
    };
    compatibility::emit_deprecation_warnings(&cli.command);
    match cli.command {
        Command::Germline(args) => application::compare::run(args),
        Command::Somatic(args) => application::somatic::run(args),
        Command::Preprocess(args) if args.version => {
            println!("pre.py ");
            Ok(())
        }
        Command::Preprocess(args) => application::preprocess::run(args),
        Command::Ftx(args) => application::ftx::run(args),
        Command::Quantify(_) if quantify_version => {
            println!("qfy.py ");
            Ok(())
        }
        Command::Quantify(args) => application::quantify::run(args),
        Command::Validate(args) => application::validate::run(args),
    }
}
