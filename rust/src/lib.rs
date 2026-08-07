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

use anyhow::Result;
use clap::{CommandFactory, Parser, error::ErrorKind};
use cli_compat::cli::{
    Cli, Command, legacy_unknown_argument_exit_code, process_args_with_legacy_somatic_aliases,
    requests_legacy_subcommand_version, requests_quantify_version,
    validate_legacy_germline_version_arguments,
};

fn exit_with_legacy_help(arguments: &[std::ffi::OsString], error: clap::Error, code: i32) -> ! {
    eprint!("{error}");
    let canonical = match arguments.get(1).and_then(|value| value.to_str()) {
        Some("germline" | "compare") => "germline",
        Some("pre" | "preprocess" | "prepy") => "pre",
        Some("quantify" | "qfy") => "quantify",
        _ => std::process::exit(code),
    };
    let mut command = Cli::command();
    if let Some(subcommand) = command.find_subcommand_mut(canonical) {
        print!("{}", subcommand.render_help());
    }
    std::process::exit(code)
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
                exit_with_legacy_help(&arguments, error, 1);
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
            if let Some(code) = legacy_unknown_argument_exit_code(&arguments) {
                exit_with_legacy_help(&arguments, error, code);
            }
            error.exit();
        }
        Err(error)
            if error.kind() == ErrorKind::MissingRequiredArgument
                && legacy_unknown_argument_exit_code(&arguments) == Some(1) =>
        {
            eprint!("{error}");
            std::process::exit(1);
        }
        Err(error) => error.exit(),
    };
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
