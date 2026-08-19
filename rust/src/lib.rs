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

use anyhow::Result;
use clap::Parser;
use cli_compat::cli::{
    Cli, Command, process_args_with_legacy_somatic_aliases, requests_legacy_subcommand_version,
    requests_quantify_version, validate_legacy_germline_version_arguments,
};
use cli_compat::compatibility;

/// Parse the process arguments and run the selected `hap` command.
///
/// # Errors
///
/// Returns an error when command execution, input parsing, comparison, or
/// output publication fails. Usage errors do not return: clap reports them and
/// exits with its own status, the same for every subcommand.
pub fn run() -> Result<()> {
    let arguments = process_args_with_legacy_somatic_aliases();
    if requests_legacy_subcommand_version(&arguments) {
        if let Err(error) = validate_legacy_germline_version_arguments(&arguments) {
            error.exit();
        }
        println!("Hap.py ");
        return Ok(());
    }
    let quantify_version = requests_quantify_version(&arguments);
    let cli = match Cli::try_parse_from(arguments) {
        Ok(cli) => cli,
        Err(error) => error.exit(),
    };
    compatibility::emit_compatibility_warnings(&cli.command);
    match cli.command {
        Command::Germline(args) => application::compare::run(args.try_into()?),
        Command::Somatic(args) => application::somatic::run(args.try_into()?),
        Command::Preprocess(args) if args.version => {
            println!("pre.py ");
            Ok(())
        }
        Command::Preprocess(args) => application::preprocess::run(args.try_into()?),
        Command::Ftx(args) => application::ftx::run(args.try_into()?),
        Command::Quantify(_) if quantify_version => {
            println!("qfy.py ");
            Ok(())
        }
        Command::Quantify(args) => application::quantify::run(args.try_into()?),
        Command::Validate(args) => application::validate::run(args.try_into()?),
    }
}
