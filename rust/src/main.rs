use anyhow::Result;
use clap::Parser;
use hap_rs::cli::{
    Cli, Command, process_args_with_legacy_somatic_aliases, requests_legacy_subcommand_version,
    requests_quantify_version,
};

fn main() -> Result<()> {
    let arguments = process_args_with_legacy_somatic_aliases();
    if requests_legacy_subcommand_version(&arguments) {
        println!("Hap.py ");
        return Ok(());
    }
    let quantify_version = requests_quantify_version(&arguments);
    let cli = Cli::parse_from(arguments);
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
