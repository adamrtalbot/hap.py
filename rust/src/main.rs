use anyhow::Result;
use clap::Parser;
use hap_rs::cli::{Cli, Command, process_args_with_legacy_somatic_aliases};

fn main() -> Result<()> {
    let arguments = process_args_with_legacy_somatic_aliases();
    if arguments
        .get(1)
        .and_then(|value| value.to_str())
        .is_some_and(|command| matches!(command, "germline" | "compare"))
        && arguments
            .iter()
            .skip(2)
            .any(|argument| argument == "-v" || argument == "--version")
    {
        println!("{}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    let cli = Cli::parse_from(arguments);
    match cli.command {
        Command::Germline(args) => hap_rs::compare::run(args),
        Command::Somatic(args) => hap_rs::somatic::run(args),
        Command::Preprocess(args) => hap_rs::preprocess::run(args),
        Command::Ftx(args) => hap_rs::ftx::run(args),
        Command::Quantify(args) => hap_rs::quantify::run(args),
        Command::Validate(args) => hap_rs::validate::run(args),
    }
}
