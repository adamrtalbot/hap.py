use anyhow::Result;
use clap::{CommandFactory, Parser, error::ErrorKind};
use hap_rs::cli::{
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

fn main() -> Result<()> {
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
        Command::Germline(args) => hap_rs::compare::run(args),
        Command::Somatic(args) => hap_rs::somatic::run(args),
        Command::Preprocess(args) if args.version => {
            println!("pre.py ");
            Ok(())
        }
        Command::Preprocess(args) => hap_rs::preprocess::run(args),
        Command::Ftx(args) => hap_rs::ftx::run(args),
        Command::Quantify(_) if quantify_version => {
            println!("qfy.py ");
            Ok(())
        }
        Command::Quantify(args) => hap_rs::quantify::run(args),
        Command::Validate(args) => hap_rs::validate::run(args),
    }
}
