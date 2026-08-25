use clap::{
    Arg, ArgAction, ArgMatches, Args, Command as ClapCommand, CommandFactory, Error,
    FromArgMatches, Parser, Subcommand, ValueEnum,
};

#[derive(Parser, Debug)]
#[command(
    name = "hap",
    version,
    about = "Single-binary Rust haplotype comparison tool"
)]
pub(crate) struct Cli {
    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
#[allow(clippy::large_enum_variant)]
pub(crate) enum Command {
    #[command(alias = "compare", args_override_self = true)]
    Germline(CompareArgs),
    #[command(args_override_self = true)]
    Somatic(SomaticArgs),
    #[command(
        name = "pre",
        alias = "preprocess",
        alias = "prepy",
        args_override_self = true
    )]
    Preprocess(PreprocessArgs),
    #[command(alias = "ftxpy", args_override_self = true)]
    Ftx(FtxArgs),
    #[command(visible_alias = "qfy", args_override_self = true)]
    Quantify(QuantifyArgs),
    #[command(visible_alias = "vcfcheck")]
    Validate(ValidateArgs),
}

#[derive(Args, Debug, Clone)]
pub(crate) struct CompareArgs {
    pub truth: String,
    pub query: String,

    #[arg(short = 'r', long = "reference", default_value = "")]
    pub reference: String,

    #[arg(short = 'o', long = "report-prefix")]
    pub report_prefix: String,

    #[arg(short = 'v', long = "version", default_value_t = false)]
    pub version: bool,

    /// Accepted for wrapper compatibility. hap.py records this initial value
    /// in runinfo, then selects the actual quantifier type from `--engine`.
    #[arg(short = 't', long = "type", value_parser = ["xcmp", "ga4gh"])]
    pub annotation_type: Option<String>,

    #[arg(long = "pass-only", default_value_t = false)]
    pub pass_only: bool,

    #[arg(long = "preprocess-truth", default_value_t = false)]
    pub preprocess_truth: bool,

    #[arg(long = "convert-gvcf-truth", default_value_t = false)]
    pub convert_gvcf_truth: bool,

    #[arg(long = "convert-gvcf-query", default_value_t = false)]
    pub convert_gvcf_query: bool,

    #[arg(long = "convert-gvcf-to-vcf", default_value_t = false)]
    pub convert_gvcf_to_vcf: bool,

    #[arg(long = "usefiltered-truth", default_value_t = false)]
    pub usefiltered_truth: bool,

    #[arg(long = "filters-only")]
    pub filters_only: Option<String>,

    #[arg(long = "preprocessing-window-size", default_value_t = 10_000)]
    pub preprocess_window: usize,

    #[arg(
        long = "adjust-conf-regions",
        action = ArgAction::SetTrue,
        default_value_t = true,
        overrides_with = "no_adjust_conf_regions"
    )]
    pub adjust_conf_regions: bool,

    #[arg(
        long = "no-adjust-conf-regions",
        action = ArgAction::SetTrue,
        overrides_with = "adjust_conf_regions"
    )]
    pub no_adjust_conf_regions: bool,

    #[arg(
        short = 'L',
        long = "leftshift",
        action = ArgAction::SetTrue,
        overrides_with = "no_leftshift"
    )]
    pub leftshift: bool,

    #[arg(
        long = "no-leftshift",
        action = ArgAction::SetTrue,
        overrides_with = "leftshift"
    )]
    pub no_leftshift: bool,

    #[arg(
        long = "decompose",
        action = ArgAction::SetTrue,
        overrides_with = "no_decompose"
    )]
    pub decompose: bool,

    #[arg(
        short = 'D',
        long = "no-decompose",
        action = ArgAction::SetTrue,
        overrides_with = "decompose"
    )]
    pub no_decompose: bool,

    #[arg(long = "bcftools-norm", default_value_t = false)]
    pub bcftools_norm: bool,

    #[arg(
        long = "fixchr",
        default_missing_value = "true",
        num_args = 0,
        overrides_with = "no_fixchr"
    )]
    pub fixchr: Option<bool>,

    #[arg(long = "no-fixchr", default_value_t = false, overrides_with = "fixchr")]
    pub no_fixchr: bool,

    #[arg(long = "filter-nonref", default_value_t = false)]
    pub filter_nonref: bool,

    #[arg(long = "somatic", default_value_t = false, overrides_with = "set_gt")]
    pub somatic: bool,

    #[arg(long = "set-gt", overrides_with = "somatic")]
    pub set_gt: Option<SomaticGtMode>,

    #[arg(long = "gender", value_enum, default_value_t = PreprocessGender::Auto)]
    pub gender: PreprocessGender,

    #[arg(long = "bcf", default_value_t = false)]
    pub bcf: bool,

    #[arg(short = 'R', long = "restrict-regions")]
    pub regions_bedfile: Option<String>,

    #[arg(short = 'T', long = "target-regions")]
    pub targets_bedfile: Option<String>,

    #[arg(short = 'f', long = "false-positives")]
    pub fp_bedfile: Option<String>,

    #[arg(short = 'l', long = "location")]
    pub locations: Option<String>,

    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long = "stratification")]
    pub strat_tsv: Option<String>,

    #[arg(long = "stratification-region", action = ArgAction::Append)]
    pub strat_regions: Vec<String>,

    #[arg(long = "stratification-fixchr", default_value_t = false)]
    pub strat_fixchr: bool,

    #[arg(short = 'V', long = "write-vcf", default_value_t = false)]
    pub write_vcf: bool,

    #[arg(
        short = 'X',
        long = "write-counts",
        action = ArgAction::SetTrue,
        default_value_t = true,
        overrides_with = "no_write_counts"
    )]
    pub write_counts: bool,

    #[arg(
        long = "no-write-counts",
        action = ArgAction::SetTrue,
        overrides_with = "write_counts"
    )]
    pub no_write_counts: bool,

    #[arg(long = "output-vtc", default_value_t = false)]
    pub output_vtc: bool,

    #[arg(long = "preserve-info", default_value_t = false)]
    pub preserve_info: bool,

    #[arg(long = "roc", default_value = "QUAL")]
    pub roc: String,

    #[arg(long = "no-roc", default_value_t = false)]
    pub no_roc: bool,

    #[arg(long = "roc-regions", action = ArgAction::Append)]
    pub roc_regions: Vec<String>,

    #[arg(long = "roc-filter")]
    pub roc_filter: Option<String>,

    #[arg(long = "roc-delta", default_value_t = 0.5)]
    pub roc_delta: f64,

    #[arg(long = "ci-alpha", default_value_t = 0.0)]
    pub ci_alpha: f64,

    #[arg(long = "no-json", default_value_t = false)]
    pub no_json: bool,

    #[arg(
        long = "unhappy",
        visible_alias = "no-haplotype-comparison",
        default_value_t = false
    )]
    pub no_hc: bool,

    #[arg(short = 'w', long = "window-size", default_value_t = 50)]
    pub window: usize,

    #[arg(long = "xcmp-enumeration-threshold", default_value_t = 16_768)]
    pub max_enum: usize,

    #[arg(long = "xcmp-expand-hapblocks", default_value_t = 30)]
    pub hb_expand: usize,

    #[arg(long = "engine", value_enum, default_value_t = CompareEngine::Xcmp)]
    pub engine: CompareEngine,

    #[arg(long = "engine-vcfeval-path")]
    pub engine_vcfeval: Option<String>,

    #[arg(long = "engine-vcfeval-template")]
    pub engine_vcfeval_template: Option<String>,

    #[arg(
        long = "scmp-distance",
        visible_alias = "lose-match-distance",
        default_value_t = 30
    )]
    pub engine_scmp_distance: usize,

    #[arg(long = "force-interactive", default_value_t = false)]
    pub force_interactive: bool,

    #[arg(long = "scratch-prefix")]
    pub scratch_prefix: Option<String>,

    #[arg(long = "keep-scratch", default_value_t = false)]
    pub keep_scratch: bool,

    #[arg(long = "logfile")]
    pub logfile: Option<String>,

    #[arg(long = "verbose", default_value_t = false, conflicts_with = "quiet")]
    pub verbose: bool,

    #[arg(long = "quiet", default_value_t = false, conflicts_with = "verbose")]
    pub quiet: bool,
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, ValueEnum)]
pub(crate) enum CompareEngine {
    #[default]
    Xcmp,
    Vcfeval,
    ScmpSomatic,
    ScmpDistance,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub(crate) enum SomaticGtMode {
    Half,
    Hemi,
    Het,
    Hom,
    First,
}

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq, ValueEnum)]
pub(crate) enum PreprocessGender {
    Male,
    Female,
    #[default]
    Auto,
    None,
}

#[derive(Args, Debug, Clone)]
pub(crate) struct PreprocessArgs {
    pub input: String,
    pub output: String,

    #[arg(short = 'v', long = "version", default_value_t = false)]
    pub version: bool,

    #[arg(short = 'r', long = "reference")]
    pub reference: Option<String>,

    #[arg(short = 'l', long = "location")]
    pub locations: Option<String>,

    #[arg(long = "pass-only", default_value_t = false)]
    pub pass_only: bool,

    #[arg(long = "filters-only")]
    pub filters_only: Option<String>,

    #[arg(short = 'R', long = "restrict-regions")]
    pub regions_bedfile: Option<String>,

    #[arg(short = 'T', long = "target-regions")]
    pub targets_bedfile: Option<String>,

    #[arg(
        long = "fixchr",
        default_missing_value = "true",
        num_args = 0,
        overrides_with = "no_fixchr"
    )]
    pub fixchr: Option<bool>,

    #[arg(long = "no-fixchr", default_value_t = false, overrides_with = "fixchr")]
    pub no_fixchr: bool,

    #[arg(long = "somatic", default_value_t = false, overrides_with = "set_gt")]
    pub somatic: bool,

    #[arg(long = "set-gt", overrides_with = "somatic")]
    pub set_gt: Option<SomaticGtMode>,

    // Standalone pre.py fails to forward this argparse value and therefore
    // executes preprocess()'s historical default=true. Germline constructs
    // PreprocessArgs directly and can still supply false, matching hap.py.
    #[arg(long = "filter-nonref", default_value_t = true)]
    pub filter_nonref: bool,

    #[arg(long = "convert-gvcf-to-vcf", default_value_t = false)]
    pub convert_gvcf_to_vcf: bool,

    #[arg(long = "bcf", default_value_t = false)]
    pub bcf: bool,

    #[arg(long = "bcftools-norm", default_value_t = false)]
    pub bcftools_norm: bool,

    /// Whether eligible indels are normalized to their left-most
    /// representation. Germline supplies side-specific values because legacy
    /// hap.py preprocesses the query by default but leaves truth unchanged.
    #[arg(
        short = 'L',
        long = "leftshift",
        action = ArgAction::SetTrue,
        default_value_t = true,
        overrides_with = "no_leftshift"
    )]
    pub leftshift: bool,

    #[arg(
        long = "no-leftshift",
        action = ArgAction::SetTrue,
        overrides_with = "leftshift"
    )]
    pub no_leftshift: bool,

    /// Whether to run `variant_pipeline::primitive_split` on each record.
    /// `hap pre` and `hap ftx` keep it enabled (the legacy preprocessing
    /// stack decomposes mixed multi-allelic indels); `hap germline` disables
    /// it so truth/query counts aggregate at the variant-site level the way
    /// legacy xcmp does post-comparison.
    #[arg(
        long = "decompose",
        action = ArgAction::SetTrue,
        default_value_t = true,
        overrides_with = "no_decompose"
    )]
    pub decompose: bool,

    #[arg(
        short = 'D',
        long = "no-decompose",
        action = ArgAction::SetTrue,
        overrides_with = "decompose"
    )]
    pub no_decompose: bool,

    #[arg(long = "gender", value_enum, default_value_t = PreprocessGender::Auto)]
    pub gender: PreprocessGender,

    #[arg(
        short = 'w',
        long = "window-size",
        default_value_t = 10_000,
        allow_hyphen_values = true
    )]
    pub window_size: i64,

    #[arg(long = "threads")]
    pub threads: Option<usize>,

    #[arg(long = "logfile")]
    pub logfile: Option<String>,

    #[arg(long = "verbose", default_value_t = false, conflicts_with = "quiet")]
    pub verbose: bool,

    #[arg(long = "quiet", default_value_t = false, conflicts_with = "verbose")]
    pub quiet: bool,

    #[arg(long = "force-interactive", default_value_t = false)]
    pub force_interactive: bool,
}

#[derive(Debug, Clone)]
pub(crate) struct QuantifyArgs {
    pub input_vcf: String,
    pub report_prefix: String,
    pub reference: String,
    pub annotation_type: Option<String>,
    pub fp_bedfile: Option<String>,
    pub strat_tsv: Option<String>,
    pub strat_regions: Vec<String>,
    pub strat_fixchr: bool,
    pub write_vcf: bool,
    pub write_counts: bool,
    pub output_vtc: bool,
    pub preserve_info: bool,
    pub adjust_conf_regions: Option<String>,
    pub threads: Option<usize>,
    pub bcf: bool,
    pub logfile: Option<String>,
    pub verbose: bool,
    pub quiet: bool,
    pub force_interactive: bool,
    pub roc: String,
    pub do_roc: bool,
    pub roc_regions: Vec<String>,
    pub roc_filter: Option<String>,
    pub roc_delta: f64,
    pub ci_alpha: f64,
    pub no_json: bool,
}

impl FromArgMatches for QuantifyArgs {
    fn from_arg_matches(matches: &ArgMatches) -> Result<Self, Error> {
        let mut matches = matches.clone();
        Self::from_arg_matches_mut(&mut matches)
    }

    fn from_arg_matches_mut(matches: &mut ArgMatches) -> Result<Self, Error> {
        Ok(Self {
            input_vcf: matches
                .remove_one::<String>("input_vcf")
                .expect("required by clap"),
            report_prefix: matches
                .remove_one::<String>("report_prefix")
                .expect("required by clap"),
            reference: matches
                .remove_one::<String>("reference")
                .expect("required by clap"),
            annotation_type: matches.remove_one::<String>("annotation_type"),
            fp_bedfile: matches.remove_one::<String>("fp_bedfile"),
            strat_tsv: matches.remove_one::<String>("strat_tsv"),
            strat_regions: matches
                .remove_many::<String>("strat_regions")
                .map(Iterator::collect)
                .unwrap_or_default(),
            strat_fixchr: matches.get_flag("strat_fixchr"),
            write_vcf: matches.get_flag("write_vcf"),
            write_counts: !matches.get_flag("no_write_counts"),
            output_vtc: matches.get_flag("output_vtc"),
            preserve_info: matches.get_flag("preserve_info"),
            adjust_conf_regions: matches.remove_one::<String>("adjust_conf_regions"),
            threads: matches.remove_one::<usize>("threads"),
            bcf: matches.get_flag("bcf"),
            logfile: matches.remove_one::<String>("logfile"),
            verbose: matches.get_flag("verbose"),
            quiet: matches.get_flag("quiet"),
            force_interactive: matches.get_flag("force_interactive"),
            roc: matches
                .remove_one::<String>("roc")
                .unwrap_or_else(|| "QUAL".to_string()),
            do_roc: !matches.get_flag("no_roc"),
            roc_regions: {
                let mut regions = vec!["*".to_string()];
                regions.extend(
                    matches
                        .remove_many::<String>("roc_regions")
                        .into_iter()
                        .flatten(),
                );
                regions
            },
            roc_filter: matches.remove_one::<String>("roc_filter"),
            roc_delta: matches.remove_one::<f64>("roc_delta").unwrap_or(0.5),
            ci_alpha: matches.remove_one::<f64>("ci_alpha").unwrap_or(0.0),
            no_json: matches.get_flag("no_json"),
        })
    }

    fn update_from_arg_matches(&mut self, matches: &ArgMatches) -> Result<(), Error> {
        let mut matches = matches.clone();
        self.update_from_arg_matches_mut(&mut matches)
    }

    fn update_from_arg_matches_mut(&mut self, matches: &mut ArgMatches) -> Result<(), Error> {
        if let Some(input_vcf) = matches.remove_one::<String>("input_vcf") {
            self.input_vcf = input_vcf;
        }
        if let Some(report_prefix) = matches.remove_one::<String>("report_prefix") {
            self.report_prefix = report_prefix;
        }
        if let Some(reference) = matches.remove_one::<String>("reference") {
            self.reference = reference;
        }
        if let Some(annotation_type) = matches.remove_one::<String>("annotation_type") {
            self.annotation_type = Some(annotation_type);
        }
        if let Some(fp_bedfile) = matches.remove_one::<String>("fp_bedfile") {
            self.fp_bedfile = Some(fp_bedfile);
        }
        if let Some(strat_tsv) = matches.remove_one::<String>("strat_tsv") {
            self.strat_tsv = Some(strat_tsv);
        }
        if let Some(strat_regions) = matches.remove_many::<String>("strat_regions") {
            self.strat_regions.extend(strat_regions);
        }
        self.strat_fixchr |= matches.get_flag("strat_fixchr");
        self.write_vcf |= matches.get_flag("write_vcf");
        if matches.get_flag("write_counts") {
            self.write_counts = true;
        }
        if matches.get_flag("no_write_counts") {
            self.write_counts = false;
        }
        self.output_vtc |= matches.get_flag("output_vtc");
        self.preserve_info |= matches.get_flag("preserve_info");
        if let Some(adjust_conf_regions) = matches.remove_one::<String>("adjust_conf_regions") {
            self.adjust_conf_regions = Some(adjust_conf_regions);
        }
        if let Some(threads) = matches.remove_one::<usize>("threads") {
            self.threads = Some(threads);
        }
        self.bcf |= matches.get_flag("bcf");
        if let Some(logfile) = matches.remove_one::<String>("logfile") {
            self.logfile = Some(logfile);
        }
        self.verbose |= matches.get_flag("verbose");
        self.quiet |= matches.get_flag("quiet");
        self.force_interactive |= matches.get_flag("force_interactive");
        if let Some(roc) = matches.remove_one::<String>("roc") {
            self.roc = roc;
        }
        if matches.get_flag("no_roc") {
            self.do_roc = false;
        }
        if let Some(roc_regions) = matches.remove_many::<String>("roc_regions") {
            self.roc_regions.extend(roc_regions);
        }
        if let Some(roc_filter) = matches.remove_one::<String>("roc_filter") {
            self.roc_filter = Some(roc_filter);
        }
        if let Some(roc_delta) = matches.remove_one::<f64>("roc_delta") {
            self.roc_delta = roc_delta;
        }
        if let Some(ci_alpha) = matches.remove_one::<f64>("ci_alpha") {
            self.ci_alpha = ci_alpha;
        }
        self.no_json |= matches.get_flag("no_json");
        Ok(())
    }
}

impl Args for QuantifyArgs {
    fn augment_args(command: ClapCommand) -> ClapCommand {
        augment_quantify_args(command, true)
    }

    fn augment_args_for_update(command: ClapCommand) -> ClapCommand {
        augment_quantify_args(command, false)
    }
}

fn augment_quantify_args(command: ClapCommand, required: bool) -> ClapCommand {
    command
        .arg(Arg::new("input_vcf").required(required))
        .arg(
            Arg::new("report_prefix")
                .short('o')
                .long("report-prefix")
                .required(required),
        )
        .arg(
            Arg::new("reference")
                .short('r')
                .long("reference")
                .required(required),
        )
        .arg(
            Arg::new("version")
                .short('v')
                .long("version")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("annotation_type")
                .short('t')
                .long("type")
                .value_parser(["xcmp", "ga4gh"]),
        )
        .arg(Arg::new("fp_bedfile").short('f').long("false-positives"))
        .arg(Arg::new("strat_tsv").long("stratification"))
        .arg(
            Arg::new("strat_regions")
                .long("stratification-region")
                .action(ArgAction::Append),
        )
        .arg(
            Arg::new("strat_fixchr")
                .long("stratification-fixchr")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("write_vcf")
                .short('V')
                .long("write-vcf")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("write_counts")
                .short('X')
                .long("write-counts")
                .action(ArgAction::SetTrue)
                .overrides_with("no_write_counts"),
        )
        .arg(
            Arg::new("no_write_counts")
                .long("no-write-counts")
                .action(ArgAction::SetTrue)
                .overrides_with("write_counts"),
        )
        .arg(
            Arg::new("output_vtc")
                .long("output-vtc")
                .action(ArgAction::SetTrue),
        )
        .arg(
            Arg::new("preserve_info")
                .long("preserve-info")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("adjust_conf_regions").long("adjust-conf-regions"))
        .arg(
            Arg::new("threads")
                .long("threads")
                .value_parser(clap::value_parser!(usize)),
        )
        .arg(Arg::new("bcf").long("bcf").action(ArgAction::SetTrue))
        .arg(Arg::new("logfile").long("logfile"))
        .arg(
            Arg::new("verbose")
                .long("verbose")
                .action(ArgAction::SetTrue)
                .conflicts_with("quiet"),
        )
        .arg(
            Arg::new("quiet")
                .long("quiet")
                .action(ArgAction::SetTrue)
                .conflicts_with("verbose"),
        )
        .arg(
            Arg::new("force_interactive")
                .long("force-interactive")
                .action(ArgAction::SetTrue),
        )
        .arg(Arg::new("roc").long("roc").default_value("QUAL"))
        .arg(Arg::new("no_roc").long("no-roc").action(ArgAction::SetTrue))
        .arg(
            Arg::new("roc_regions")
                .long("roc-regions")
                .action(ArgAction::Append),
        )
        .arg(Arg::new("roc_filter").long("roc-filter"))
        .arg(
            Arg::new("roc_delta")
                .long("roc-delta")
                .default_value("0.5")
                .value_parser(clap::value_parser!(f64)),
        )
        .arg(
            Arg::new("ci_alpha")
                .long("ci-alpha")
                .default_value("0.0")
                .value_parser(clap::value_parser!(f64)),
        )
        .arg(
            Arg::new("no_json")
                .long("no-json")
                .action(ArgAction::SetTrue),
        )
}

#[derive(Args, Debug, Clone)]
pub(crate) struct SomaticArgs {
    pub truth: String,
    pub query: String,

    #[arg(short = 'o', long = "output")]
    pub output: String,

    #[arg(short = 'r', long = "reference")]
    pub reference: Option<String>,

    #[arg(short = 'l', long = "location")]
    pub location: Option<String>,

    #[arg(short = 'R', long = "restrict-regions")]
    pub regions_bedfile: Option<String>,

    #[arg(short = 'T', long = "target-regions")]
    pub targets_bedfile: Option<String>,

    #[arg(short = 'f', long = "false-positives")]
    pub fp_bedfile: Option<String>,

    #[arg(short = 'a', long = "ambiguous")]
    pub ambiguous_beds: Vec<String>,

    #[arg(
        long = "ambi-fp",
        default_value_t = false,
        overrides_with = "no_ambi_fp"
    )]
    pub ambi_fp: bool,

    #[arg(
        long = "no-ambi-fp",
        default_value_t = false,
        overrides_with = "ambi_fp"
    )]
    pub no_ambi_fp: bool,

    #[arg(
        long = "count-unk",
        default_value_t = false,
        overrides_with = "no_count_unk"
    )]
    pub count_unk: bool,

    #[arg(
        long = "no-count-unk",
        default_value_t = false,
        overrides_with = "count_unk"
    )]
    pub no_count_unk: bool,

    #[arg(short = 'e', long = "explain_ambiguous", default_value_t = false)]
    pub explain_ambiguous: bool,

    #[arg(short = 'P', long = "include-nonpass", default_value_t = false)]
    pub include_nonpass: bool,

    #[arg(long = "fp-region-size")]
    pub fp_region_size: Option<String>,

    #[arg(long = "feature-table", value_parser = parse_somatic_feature_table)]
    pub feature_table: Option<SomaticFeatureTable>,

    #[arg(long = "happy-stats", default_value_t = false)]
    pub happy_stats: bool,

    #[arg(long = "bam")]
    pub bams: Vec<String>,

    #[arg(long = "normalize-truth", default_value_t = false)]
    pub normalize_truth: bool,

    #[arg(long = "normalize-query", default_value_t = false)]
    pub normalize_query: bool,

    #[arg(short = 'N', long = "normalize-all", default_value_t = false)]
    pub normalize_all: bool,

    #[arg(
        long = "fixchr-truth",
        visible_alias = "fix-chr-truth",
        default_missing_value = "true",
        num_args = 0,
        overrides_with = "no_fixchr_truth"
    )]
    pub fixchr_truth: Option<bool>,

    #[arg(
        long = "fixchr-query",
        visible_alias = "fix-chr-query",
        default_missing_value = "true",
        num_args = 0,
        overrides_with = "no_fixchr_query"
    )]
    pub fixchr_query: Option<bool>,

    #[arg(
        long = "no-fixchr-truth",
        default_value_t = false,
        overrides_with = "fixchr_truth"
    )]
    pub no_fixchr_truth: bool,

    #[arg(
        long = "no-fixchr-query",
        default_value_t = false,
        overrides_with = "fixchr_query"
    )]
    pub no_fixchr_query: bool,

    #[arg(long = "no-order-check", default_value_t = false)]
    pub no_order_check: bool,

    #[arg(long = "roc", value_parser = parse_somatic_roc)]
    pub roc: Option<SomaticRoc>,

    #[arg(long = "bin-afs", default_value_t = false)]
    pub af_strat: bool,

    #[arg(long = "af-binsize", default_value = "0.2", allow_hyphen_values = true)]
    pub af_strat_binsize: String,

    #[arg(long = "af-truth", default_value = "I.T_ALT_RATE")]
    pub af_strat_truth: String,

    #[arg(long = "af-query", default_value = "T_AF")]
    pub af_strat_query: String,

    #[arg(long = "count-filtered-fn", default_value_t = false)]
    pub count_filtered_fn: bool,

    #[arg(long = "ci-level", default_value_t = 0.95)]
    pub ci_level: f64,

    #[arg(long = "scratch-prefix")]
    pub scratch_prefix: Option<String>,

    #[arg(long = "keep-scratch", default_value_t = false)]
    pub keep_scratch: bool,

    #[arg(long = "logfile")]
    pub logfile: Option<String>,

    #[arg(long = "verbose", default_value_t = false, conflicts_with = "quiet")]
    pub verbose: bool,

    #[arg(long = "quiet", default_value_t = false, conflicts_with = "verbose")]
    pub quiet: bool,
}

#[derive(Args, Debug, Clone)]
pub(crate) struct FtxArgs {
    pub input: String,

    #[arg(short = 'o', long = "output")]
    pub output: String,

    #[arg(short = 'l', long = "location")]
    pub location: Option<String>,

    #[arg(short = 'R', long = "restrict-regions")]
    pub regions_bedfile: Option<String>,

    #[arg(short = 'T', long = "target-regions")]
    pub targets_bedfile: Option<String>,

    #[arg(short = 'P', long = "include-nonpass", default_value_t = false)]
    pub include_nonpass: bool,

    #[arg(long = "feature-table", default_value = "generic")]
    pub features: String,

    #[arg(long = "feature-label")]
    pub label: Option<String>,

    #[arg(long = "bam")]
    pub bams: Vec<String>,

    #[arg(short = 'r', long = "reference")]
    pub reference: Option<String>,

    #[arg(long = "normalize", default_value_t = false)]
    pub normalize: bool,

    #[arg(long = "fix-chr", default_value_t = false)]
    pub fixchr: bool,

    // Declared extension beyond legacy ftx.py, mirroring `germline`/`somatic`.
    #[arg(long = "scratch-prefix")]
    pub scratch_prefix: Option<String>,

    #[arg(long = "keep-scratch", default_value_t = false)]
    pub keep_scratch: bool,
}

/// Closed set of somatic feature tables accepted by `--feature-table`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SomaticFeatureTable {
    Generic,
    AdmixStrelkaSnv,
    AdmixStrelkaIndel,
    HccStrelkaSnv,
    HccStrelkaIndel,
    HccMutectSnv,
    HccMutectIndel,
    HccVarscan2Snv,
    HccVarscan2Indel,
    HccPiscesSnv,
    HccPiscesIndel,
}

impl SomaticFeatureTable {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Generic => "generic",
            Self::AdmixStrelkaSnv => "admix.strelka.snv",
            Self::AdmixStrelkaIndel => "admix.strelka.indel",
            Self::HccStrelkaSnv => "hcc.strelka.snv",
            Self::HccStrelkaIndel => "hcc.strelka.indel",
            Self::HccMutectSnv => "hcc.mutect.snv",
            Self::HccMutectIndel => "hcc.mutect.indel",
            Self::HccVarscan2Snv => "hcc.varscan2.snv",
            Self::HccVarscan2Indel => "hcc.varscan2.indel",
            Self::HccPiscesSnv => "hcc.pisces.snv",
            Self::HccPiscesIndel => "hcc.pisces.indel",
        }
    }
}

/// Closed set of somatic ROC modes accepted by `--roc`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SomaticRoc {
    StrelkaSnvQss,
    StrelkaSnvVqsr,
    StrelkaSnv,
    StrelkaIndel,
    StrelkaIndelEvs,
    Varscan2Snv,
    Varscan2Indel,
    MutectSnv,
    MutectIndel,
}

impl SomaticRoc {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::StrelkaSnvQss => "strelka.snv.qss",
            Self::StrelkaSnvVqsr => "strelka.snv.vqsr",
            Self::StrelkaSnv => "strelka.snv",
            Self::StrelkaIndel => "strelka.indel",
            Self::StrelkaIndelEvs => "strelka.indel.evs",
            Self::Varscan2Snv => "varscan2.snv",
            Self::Varscan2Indel => "varscan2.indel",
            Self::MutectSnv => "mutect.snv",
            Self::MutectIndel => "mutect.indel",
        }
    }
}

fn parse_somatic_feature_table(value: &str) -> Result<SomaticFeatureTable, String> {
    match value {
        "generic" => Ok(SomaticFeatureTable::Generic),
        "admix.strelka.snv" => Ok(SomaticFeatureTable::AdmixStrelkaSnv),
        "admix.strelka.indel" => Ok(SomaticFeatureTable::AdmixStrelkaIndel),
        "hcc.strelka.snv" => Ok(SomaticFeatureTable::HccStrelkaSnv),
        "hcc.strelka.indel" => Ok(SomaticFeatureTable::HccStrelkaIndel),
        "hcc.mutect.snv" => Ok(SomaticFeatureTable::HccMutectSnv),
        "hcc.mutect.indel" => Ok(SomaticFeatureTable::HccMutectIndel),
        "hcc.varscan2.snv" => Ok(SomaticFeatureTable::HccVarscan2Snv),
        "hcc.varscan2.indel" => Ok(SomaticFeatureTable::HccVarscan2Indel),
        "hcc.pisces.snv" => Ok(SomaticFeatureTable::HccPiscesSnv),
        "hcc.pisces.indel" => Ok(SomaticFeatureTable::HccPiscesIndel),
        _ => Err(format!("unsupported somatic feature table '{value}'")),
    }
}

fn parse_somatic_roc(value: &str) -> Result<SomaticRoc, String> {
    match value {
        "strelka.snv.qss" => Ok(SomaticRoc::StrelkaSnvQss),
        "strelka.snv.vqsr" => Ok(SomaticRoc::StrelkaSnvVqsr),
        "strelka.snv" => Ok(SomaticRoc::StrelkaSnv),
        "strelka.indel" => Ok(SomaticRoc::StrelkaIndel),
        "strelka.indel.evs" => Ok(SomaticRoc::StrelkaIndelEvs),
        "varscan2.snv" => Ok(SomaticRoc::Varscan2Snv),
        "varscan2.indel" => Ok(SomaticRoc::Varscan2Indel),
        "mutect.snv" => Ok(SomaticRoc::MutectSnv),
        "mutect.indel" => Ok(SomaticRoc::MutectIndel),
        _ => Err(format!("unsupported somatic ROC mode '{value}'")),
    }
}

/// argparse accepted the historical multi-character short option `-FN`.
/// clap treats that spelling as a cluster (`-F -N`), so normalize only that
/// exact somatic token before parsing and leave every other command untouched.
pub(crate) fn process_args_with_legacy_somatic_aliases() -> Vec<std::ffi::OsString> {
    normalize_legacy_arguments(std::env::args_os().collect())
}

fn normalize_legacy_arguments(mut arguments: Vec<std::ffi::OsString>) -> Vec<std::ffi::OsString> {
    if arguments.get(1).and_then(|value| value.to_str()) == Some("somatic") {
        for argument in arguments.iter_mut().skip(2) {
            if argument == "-FN" {
                *argument = "--count-filtered-fn".into();
            }
        }
    }

    let canonical = match arguments.get(1).and_then(|value| value.to_str()) {
        Some("germline" | "compare") => "germline",
        Some("somatic") => "somatic",
        Some("pre" | "preprocess" | "prepy") => "pre",
        Some("ftx" | "ftxpy") => "ftx",
        Some("quantify" | "qfy") => "quantify",
        // vcfcheck uses Boost.Program_options, which does not implement
        // argparse's unique long-option abbreviations.
        _ => return arguments,
    };
    let command = Cli::command();
    let Some(subcommand) = command.find_subcommand(canonical) else {
        return arguments;
    };
    let long_options = subcommand
        .get_arguments()
        .filter_map(Arg::get_long)
        .collect::<Vec<_>>();
    let mut after_terminator = false;
    for argument in arguments.iter_mut().skip(2) {
        if after_terminator {
            continue;
        }
        let Some(value) = argument.to_str() else {
            continue;
        };
        if value == "--" {
            after_terminator = true;
            continue;
        }
        let Some(option) = value.strip_prefix("--") else {
            continue;
        };
        let (prefix, suffix) = option
            .split_once('=')
            .map_or((option, ""), |(prefix, _value)| {
                (prefix, &option[prefix.len()..])
            });
        if prefix.is_empty() || long_options.contains(&prefix) {
            continue;
        }
        let matches = long_options
            .iter()
            .filter(|candidate| candidate.starts_with(prefix))
            .copied()
            .collect::<Vec<_>>();
        if let [expanded] = matches.as_slice() {
            *argument = format!("--{expanded}{suffix}").into();
        }
    }
    arguments
}

/// hap.py handled its version flag before validating required arguments.
/// Stop at the option terminator so a positional file named `--version`
/// remains reachable.
pub(crate) fn requests_legacy_subcommand_version(arguments: &[std::ffi::OsString]) -> bool {
    let supports_version = arguments
        .get(1)
        .and_then(|value| value.to_str())
        .is_some_and(|command| matches!(command, "germline" | "compare"));

    supports_version
        && arguments
            .iter()
            .skip(2)
            .take_while(|argument| *argument != "--")
            .any(|argument| argument == "-v" || argument == "--version")
}

/// Validate all supplied germline tokens while relaxing the positional and
/// output requirements that hap.py checks only after handling `--version`.
pub(crate) fn validate_legacy_germline_version_arguments(
    arguments: &[std::ffi::OsString],
) -> Result<(), Error> {
    CompareArgs::augment_args_for_update(ClapCommand::new("germline").args_override_self(true))
        .try_get_matches_from(
            std::iter::once(std::ffi::OsString::from("germline"))
                .chain(arguments.iter().skip(2).cloned()),
        )
        .map(|_| ())
}

/// qfy.py validates its required arguments before honoring its version flag.
/// This predicate is therefore consumed only after clap parsing succeeds.
pub(crate) fn requests_quantify_version(arguments: &[std::ffi::OsString]) -> bool {
    let supports_version = arguments
        .get(1)
        .and_then(|value| value.to_str())
        .is_some_and(|command| matches!(command, "quantify" | "qfy"));

    supports_version
        && arguments
            .iter()
            .skip(2)
            .take_while(|argument| *argument != "--")
            .any(|argument| argument == "-v" || argument == "--version")
}

#[derive(Debug, Clone)]
pub(crate) struct ValidateArgs {
    pub input: String,
    pub reference: Option<String>,
    pub output_json: Option<String>,
    pub errors_bed: Option<String>,
    pub locations: Option<String>,
    pub regions_bedfile: Option<String>,
    pub targets_bedfile: Option<String>,
    pub apply_filters: bool,
    pub limit_records: Option<i64>,
    pub message_every: Option<i64>,
    pub strict_homref: bool,
    pub check_bcf_errors: bool,
    pub all_warnings: bool,
}

impl FromArgMatches for ValidateArgs {
    fn from_arg_matches(matches: &ArgMatches) -> Result<Self, Error> {
        let mut matches = matches.clone();
        Self::from_arg_matches_mut(&mut matches)
    }

    fn from_arg_matches_mut(matches: &mut ArgMatches) -> Result<Self, Error> {
        let input = matches
            .remove_one::<String>("input")
            .or_else(|| matches.remove_one::<String>("input_file"))
            .expect("one input form is required by clap");
        Ok(Self {
            input,
            reference: matches.remove_one("reference"),
            output_json: matches.remove_one("output_json"),
            errors_bed: matches.remove_one("errors_bed"),
            locations: matches.remove_one("locations"),
            regions_bedfile: matches.remove_one("regions_bedfile"),
            targets_bedfile: matches.remove_one("targets_bedfile"),
            apply_filters: matches.remove_one("apply_filters").unwrap_or(false),
            limit_records: matches.remove_one("limit_records"),
            message_every: matches.remove_one("message_every"),
            strict_homref: matches.remove_one("strict_homref").unwrap_or(false),
            check_bcf_errors: matches.remove_one("check_bcf_errors").unwrap_or(false),
            all_warnings: matches.remove_one("all_warnings").unwrap_or(false),
        })
    }

    fn update_from_arg_matches(&mut self, matches: &ArgMatches) -> Result<(), Error> {
        let mut matches = matches.clone();
        self.update_from_arg_matches_mut(&mut matches)
    }

    fn update_from_arg_matches_mut(&mut self, matches: &mut ArgMatches) -> Result<(), Error> {
        if let Some(input) = matches
            .remove_one::<String>("input")
            .or_else(|| matches.remove_one::<String>("input_file"))
        {
            self.input = input;
        }
        if let Some(value) = matches.remove_one("reference") {
            self.reference = Some(value);
        }
        if let Some(value) = matches.remove_one("output_json") {
            self.output_json = Some(value);
        }
        if let Some(value) = matches.remove_one("errors_bed") {
            self.errors_bed = Some(value);
        }
        if let Some(value) = matches.remove_one("locations") {
            self.locations = Some(value);
        }
        if let Some(value) = matches.remove_one("regions_bedfile") {
            self.regions_bedfile = Some(value);
        }
        if let Some(value) = matches.remove_one("targets_bedfile") {
            self.targets_bedfile = Some(value);
        }
        if let Some(value) = matches.remove_one("apply_filters") {
            self.apply_filters = value;
        }
        if let Some(value) = matches.remove_one("limit_records") {
            self.limit_records = Some(value);
        }
        if let Some(value) = matches.remove_one("message_every") {
            self.message_every = Some(value);
        }
        if let Some(value) = matches.remove_one("strict_homref") {
            self.strict_homref = value;
        }
        if let Some(value) = matches.remove_one("check_bcf_errors") {
            self.check_bcf_errors = value;
        }
        if let Some(value) = matches.remove_one("all_warnings") {
            self.all_warnings = value;
        }
        Ok(())
    }
}

impl Args for ValidateArgs {
    fn augment_args(command: ClapCommand) -> ClapCommand {
        augment_validate_args(command, true)
    }

    fn augment_args_for_update(command: ClapCommand) -> ClapCommand {
        augment_validate_args(command, false)
    }
}

fn augment_validate_args(command: ClapCommand, required: bool) -> ClapCommand {
    let bool_value = clap::value_parser!(bool);
    let positional_input = Arg::new("input")
        .value_name("INPUT")
        .conflicts_with("input_file");
    let positional_input = if required {
        positional_input.required_unless_present("input_file")
    } else {
        positional_input
    };
    let option_input = Arg::new("input_file")
        .long("input-file")
        .value_name("INPUT")
        .conflicts_with("input");
    let option_input = if required {
        option_input.required_unless_present("input")
    } else {
        option_input
    };
    command
        .arg(positional_input)
        .arg(option_input)
        .arg(Arg::new("reference").short('r').long("reference"))
        .arg(
            Arg::new("output_json")
                .short('o')
                .long("output-json")
                .visible_alias("output-file"),
        )
        .arg(Arg::new("errors_bed").short('e').long("errors-bed"))
        .arg(
            Arg::new("locations")
                .short('l')
                .long("location")
                .value_name("REGION"),
        )
        .arg(Arg::new("regions_bedfile").short('R').long("regions"))
        .arg(Arg::new("targets_bedfile").short('T').long("targets"))
        .arg(
            Arg::new("apply_filters")
                .short('f')
                .long("apply-filters")
                .value_parser(bool_value.clone())
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("limit_records")
                .long("limit-records")
                .value_parser(clap::value_parser!(i64)),
        )
        .arg(
            Arg::new("message_every")
                .long("message-every")
                .value_parser(clap::value_parser!(i64)),
        )
        .arg(
            Arg::new("strict_homref")
                .short('H')
                .long("strict-homref")
                .value_parser(bool_value.clone())
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("check_bcf_errors")
                .long("check-bcf-errors")
                .value_parser(bool_value.clone())
                .action(ArgAction::Set),
        )
        .arg(
            Arg::new("all_warnings")
                .short('W')
                .long("all-warnings")
                .value_parser(bool_value)
                .action(ArgAction::Set),
        )
}

impl From<CompareEngine> for crate::application::CompareEngine {
    fn from(value: CompareEngine) -> Self {
        match value {
            CompareEngine::Xcmp => Self::Xcmp,
            CompareEngine::Vcfeval => Self::Vcfeval,
            CompareEngine::ScmpSomatic => Self::ScmpSomatic,
            CompareEngine::ScmpDistance => Self::ScmpDistance,
        }
    }
}

impl From<SomaticGtMode> for crate::application::SomaticGtMode {
    fn from(value: SomaticGtMode) -> Self {
        match value {
            SomaticGtMode::Half => Self::Half,
            SomaticGtMode::Hemi => Self::Hemi,
            SomaticGtMode::Het => Self::Het,
            SomaticGtMode::Hom => Self::Hom,
            SomaticGtMode::First => Self::First,
        }
    }
}

impl From<PreprocessGender> for crate::application::PreprocessGender {
    fn from(value: PreprocessGender) -> Self {
        match value {
            PreprocessGender::Male => Self::Male,
            PreprocessGender::Female => Self::Female,
            PreprocessGender::Auto => Self::Auto,
            PreprocessGender::None => Self::None,
        }
    }
}

impl TryFrom<CompareArgs> for crate::application::ValidatedCompareArgs {
    type Error = crate::application::RequestValidationError;

    fn try_from(args: CompareArgs) -> Result<Self, Self::Error> {
        crate::application::CompareArgs {
            truth: args.truth,
            query: args.query,
            reference: args.reference,
            report_prefix: args.report_prefix,
            version: args.version,
            annotation_type: args.annotation_type,
            pass_only: args.pass_only,
            preprocess_truth: args.preprocess_truth,
            convert_gvcf_truth: args.convert_gvcf_truth,
            convert_gvcf_query: args.convert_gvcf_query,
            convert_gvcf_to_vcf: args.convert_gvcf_to_vcf,
            usefiltered_truth: args.usefiltered_truth,
            filters_only: args.filters_only,
            preprocess: crate::application::PreprocessOptions {
                preprocess_window: args.preprocess_window,
                leftshift: args.leftshift,
                no_leftshift: args.no_leftshift,
                decompose: args.decompose,
                no_decompose: args.no_decompose,
                bcftools_norm: args.bcftools_norm,
                fixchr: args.fixchr,
                no_fixchr: args.no_fixchr,
                filter_nonref: args.filter_nonref,
                set_gt: args.set_gt.map(Into::into),
                gender: args.gender.into(),
            },
            adjust_conf_regions: args.adjust_conf_regions,
            no_adjust_conf_regions: args.no_adjust_conf_regions,
            somatic: args.somatic,
            bcf: args.bcf,
            regions_bedfile: args.regions_bedfile,
            targets_bedfile: args.targets_bedfile,
            fp_bedfile: args.fp_bedfile,
            locations: args.locations,
            threads: args.threads,
            strat_tsv: args.strat_tsv,
            strat_regions: args.strat_regions,
            strat_fixchr: args.strat_fixchr,
            write_vcf: args.write_vcf,
            write_counts: args.write_counts,
            no_write_counts: args.no_write_counts,
            output_vtc: args.output_vtc,
            preserve_info: args.preserve_info,
            roc: crate::application::RocOptions {
                roc: args.roc,
                no_roc: args.no_roc,
                roc_regions: args.roc_regions,
                roc_filter: args.roc_filter,
                roc_delta: args.roc_delta,
            },
            ci_alpha: args.ci_alpha,
            no_json: args.no_json,
            no_hc: args.no_hc,
            engine: crate::application::EngineOptions {
                engine: args.engine.into(),
                engine_vcfeval: args.engine_vcfeval,
                engine_vcfeval_template: args.engine_vcfeval_template,
                engine_scmp_distance: args.engine_scmp_distance,
                window: args.window,
                max_enum: args.max_enum,
                hb_expand: args.hb_expand,
            },
            force_interactive: args.force_interactive,
            scratch_prefix: args.scratch_prefix,
            keep_scratch: args.keep_scratch,
            logfile: args.logfile,
            verbose: args.verbose,
            quiet: args.quiet,
        }
        .validated()
    }
}

impl TryFrom<PreprocessArgs> for crate::application::ValidatedPreprocessArgs {
    type Error = crate::application::RequestValidationError;

    fn try_from(args: PreprocessArgs) -> Result<Self, Self::Error> {
        crate::application::PreprocessArgs {
            input: args.input,
            output: args.output,
            version: args.version,
            reference: args.reference,
            locations: args.locations,
            pass_only: args.pass_only,
            filters_only: args.filters_only,
            regions_bedfile: args.regions_bedfile,
            targets_bedfile: args.targets_bedfile,
            fixchr: args.fixchr,
            no_fixchr: args.no_fixchr,
            somatic: args.somatic,
            set_gt: args.set_gt.map(Into::into),
            filter_nonref: args.filter_nonref,
            convert_gvcf_to_vcf: args.convert_gvcf_to_vcf,
            bcf: args.bcf,
            bcftools_norm: args.bcftools_norm,
            leftshift: args.leftshift,
            no_leftshift: args.no_leftshift,
            decompose: args.decompose,
            no_decompose: args.no_decompose,
            gender: args.gender.into(),
            window_size: args.window_size,
            threads: args.threads,
            logfile: args.logfile,
            verbose: args.verbose,
            quiet: args.quiet,
            force_interactive: args.force_interactive,
        }
        .validated()
    }
}

impl TryFrom<QuantifyArgs> for crate::application::ValidatedQuantifyArgs {
    type Error = crate::application::RequestValidationError;

    fn try_from(args: QuantifyArgs) -> Result<Self, Self::Error> {
        crate::application::QuantifyArgs {
            input_vcf: args.input_vcf,
            report_prefix: args.report_prefix,
            reference: args.reference,
            annotation_type: args.annotation_type,
            fp_bedfile: args.fp_bedfile,
            strat_tsv: args.strat_tsv,
            strat_regions: args.strat_regions,
            strat_fixchr: args.strat_fixchr,
            write_vcf: args.write_vcf,
            write_counts: args.write_counts,
            output_vtc: args.output_vtc,
            preserve_info: args.preserve_info,
            adjust_conf_regions: args.adjust_conf_regions,
            threads: args.threads,
            bcf: args.bcf,
            logfile: args.logfile,
            verbose: args.verbose,
            quiet: args.quiet,
            force_interactive: args.force_interactive,
            roc: args.roc,
            do_roc: args.do_roc,
            roc_regions: args.roc_regions,
            roc_filter: args.roc_filter,
            roc_delta: args.roc_delta,
            ci_alpha: args.ci_alpha,
            no_json: args.no_json,
        }
        .validated()
    }
}

impl TryFrom<SomaticArgs> for crate::application::ValidatedSomaticArgs {
    type Error = crate::application::RequestValidationError;

    fn try_from(args: SomaticArgs) -> Result<Self, Self::Error> {
        crate::application::SomaticArgs::from(args).validated()
    }
}

impl From<SomaticArgs> for crate::application::SomaticArgs {
    fn from(args: SomaticArgs) -> Self {
        Self {
            truth: args.truth,
            query: args.query,
            output: args.output,
            reference: args.reference,
            location: args.location,
            regions_bedfile: args.regions_bedfile,
            targets_bedfile: args.targets_bedfile,
            fp_bedfile: args.fp_bedfile,
            ambiguous_beds: args.ambiguous_beds,
            ambi_fp: args.ambi_fp,
            no_ambi_fp: args.no_ambi_fp,
            count_unk: args.count_unk,
            no_count_unk: args.no_count_unk,
            explain_ambiguous: args.explain_ambiguous,
            include_nonpass: args.include_nonpass,
            fp_region_size: args.fp_region_size,
            feature_table: args.feature_table.map(|f| f.as_str().to_string()),
            happy_stats: args.happy_stats,
            bams: args.bams,
            normalize_truth: args.normalize_truth,
            normalize_query: args.normalize_query,
            normalize_all: args.normalize_all,
            fixchr_truth: args.fixchr_truth,
            fixchr_query: args.fixchr_query,
            no_fixchr_truth: args.no_fixchr_truth,
            no_fixchr_query: args.no_fixchr_query,
            no_order_check: args.no_order_check,
            roc: args.roc.map(|r| r.as_str().to_string()),
            af_strat: args.af_strat,
            af_strat_binsize: args.af_strat_binsize,
            af_strat_truth: args.af_strat_truth,
            af_strat_query: args.af_strat_query,
            count_filtered_fn: args.count_filtered_fn,
            ci_level: args.ci_level,
            scratch_prefix: args.scratch_prefix,
            keep_scratch: args.keep_scratch,
            logfile: args.logfile,
            verbose: args.verbose,
            quiet: args.quiet,
        }
    }
}

impl TryFrom<FtxArgs> for crate::application::ValidatedFtxArgs {
    type Error = crate::application::RequestValidationError;

    fn try_from(args: FtxArgs) -> Result<Self, Self::Error> {
        crate::application::FtxArgs {
            input: args.input,
            output: args.output,
            location: args.location,
            regions_bedfile: args.regions_bedfile,
            targets_bedfile: args.targets_bedfile,
            include_nonpass: args.include_nonpass,
            features: args.features,
            label: args.label,
            bams: args.bams,
            reference: args.reference,
            normalize: args.normalize,
            fixchr: args.fixchr,
            scratch_prefix: args.scratch_prefix,
            keep_scratch: args.keep_scratch,
        }
        .validated()
    }
}

impl TryFrom<ValidateArgs> for crate::application::ValidatedValidateArgs {
    type Error = crate::application::RequestValidationError;

    fn try_from(args: ValidateArgs) -> Result<Self, Self::Error> {
        crate::application::ValidateArgs {
            input: args.input,
            reference: args.reference,
            output_json: args.output_json,
            errors_bed: args.errors_bed,
            locations: args.locations,
            regions_bedfile: args.regions_bedfile,
            targets_bedfile: args.targets_bedfile,
            apply_filters: args.apply_filters,
            limit_records: args.limit_records,
            message_every: args.message_every,
            strict_homref: args.strict_homref,
            check_bcf_errors: args.check_bcf_errors,
            all_warnings: args.all_warnings,
        }
        .validated()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::{CommandFactory, error::ErrorKind};

    #[test]
    fn public_help_lists_quantify_validate_and_legacy_aliases() {
        let help = Cli::command().render_help().to_string();
        assert!(help.contains("quantify"));
        assert!(help.contains("qfy"));
        assert!(help.contains("validate"));
        assert!(help.contains("vcfcheck"));
    }

    #[test]
    fn adapter_rejects_zero_threads() {
        let cli = Cli::try_parse_from([
            "hap",
            "compare",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-o",
            "report",
            "--threads",
            "0",
        ])
        .expect("clap should accept the syntactic thread value");
        let Command::Germline(args) = cli.command else {
            panic!("compare alias should resolve to germline");
        };

        let error = crate::application::ValidatedCompareArgs::try_from(args).unwrap_err();
        assert_eq!(error.field(), "threads");
    }

    #[test]
    fn adapter_rejects_invalid_quantify_ci_and_roc_values() {
        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
        ])
        .expect("legacy qfy alias should parse");
        let Command::Quantify(mut args) = cli.command else {
            panic!("qfy alias should resolve to quantify");
        };

        args.ci_alpha = 1.0;
        let error = crate::application::ValidatedQuantifyArgs::try_from(args.clone()).unwrap_err();
        assert_eq!(error.field(), "ci_alpha");

        args.ci_alpha = 0.0;
        args.roc_delta = f64::NAN;
        let error = crate::application::ValidatedQuantifyArgs::try_from(args).unwrap_err();
        assert_eq!(error.field(), "roc_delta");
    }

    #[test]
    fn adapter_preserves_negative_validate_controls() {
        let cli = Cli::try_parse_from(["hap", "vcfcheck", "input.vcf.gz"])
            .expect("legacy vcfcheck alias should parse");
        let Command::Validate(mut args) = cli.command else {
            panic!("vcfcheck alias should resolve to validate");
        };

        args.limit_records = Some(-2);
        args.message_every = Some(-1);
        let request = crate::application::ValidatedValidateArgs::try_from(args)
            .expect("legacy negative controls should reach the use case");
        assert_eq!(request.limit_records, Some(-2));
        assert_eq!(request.message_every, Some(-1));
    }

    #[test]
    fn legacy_alias_parses_into_validated_application_request() {
        let cli = Cli::try_parse_from([
            "hap",
            "compare",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-o",
            "report",
            "--engine",
            "scmp-distance",
            "--set-gt",
            "first",
            "--gender",
            "female",
        ])
        .expect("legacy compare spelling should parse");
        let Command::Germline(args) = cli.command else {
            panic!("compare alias should resolve to germline");
        };

        let request = crate::application::ValidatedCompareArgs::try_from(args)
            .expect("valid CLI values should produce an application request");
        assert_eq!(
            request.engine.engine,
            crate::application::CompareEngine::ScmpDistance
        );
        assert_eq!(
            request.preprocess.set_gt,
            Some(crate::application::SomaticGtMode::First)
        );
        assert_eq!(
            request.preprocess.gender,
            crate::application::PreprocessGender::Female
        );
    }

    #[test]
    fn repeated_legacy_options_match_each_wrapper_parser() {
        let cli = Cli::try_parse_from([
            "hap",
            "germline",
            "truth.vcf",
            "query.vcf",
            "-o",
            "first",
            "-o",
            "last",
            "-r",
            "first.fa",
            "-r",
            "last.fa",
            "--roc",
            "QUAL",
            "--roc",
            "INFO.QQ",
            "--pass-only",
            "--pass-only",
        ])
        .expect("argparse accepts repeated germline options");
        let Command::Germline(args) = cli.command else {
            panic!("germline should parse");
        };
        assert_eq!(args.report_prefix, "last");
        assert_eq!(args.reference, "last.fa");
        assert_eq!(args.roc, "INFO.QQ");
        assert!(args.pass_only);

        let cli = Cli::try_parse_from([
            "hap",
            "somatic",
            "truth.vcf",
            "query.vcf",
            "-o",
            "first",
            "-o",
            "last",
            "-r",
            "first.fa",
            "-r",
            "last.fa",
            "-l",
            "chr1",
            "-l",
            "chr2",
        ])
        .expect("argparse accepts repeated somatic options");
        let Command::Somatic(args) = cli.command else {
            panic!("somatic should parse");
        };
        assert_eq!(args.output, "last");
        assert_eq!(args.reference.as_deref(), Some("last.fa"));
        assert_eq!(args.location.as_deref(), Some("chr2"));

        let cli = Cli::try_parse_from([
            "hap",
            "pre",
            "input.vcf",
            "output.vcf",
            "-r",
            "first.fa",
            "-r",
            "last.fa",
            "-w",
            "10",
            "-w",
            "-1",
        ])
        .expect("argparse accepts repeated pre options");
        let Command::Preprocess(args) = cli.command else {
            panic!("pre should parse");
        };
        assert_eq!(args.reference.as_deref(), Some("last.fa"));
        assert_eq!(args.window_size, -1);

        let cli = Cli::try_parse_from([
            "hap",
            "ftx",
            "input.vcf",
            "-o",
            "first.csv",
            "-o",
            "last.csv",
            "--feature-table",
            "generic",
            "--feature-table",
            "hcc.mutect.snv",
        ])
        .expect("argparse accepts repeated ftx options");
        let Command::Ftx(args) = cli.command else {
            panic!("ftx should parse");
        };
        assert_eq!(args.output, "last.csv");
        assert_eq!(args.features, "hcc.mutect.snv");

        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "input.vcf",
            "-o",
            "first",
            "-o",
            "last",
            "-r",
            "ref.fa",
            "--roc",
            "QUAL",
            "--roc",
            "FORMAT.GQ",
        ])
        .expect("argparse accepts repeated qfy options");
        let Command::Quantify(args) = cli.command else {
            panic!("qfy should parse");
        };
        assert_eq!(args.report_prefix, "last");
        assert_eq!(args.roc, "FORMAT.GQ");

        let error = Cli::try_parse_from([
            "hap",
            "vcfcheck",
            "input.vcf",
            "-o",
            "first.json",
            "-o",
            "last.json",
            "--apply-filters",
            "true",
            "--apply-filters",
            "false",
        ])
        .expect_err("program_options rejects repeated scalar options");
        assert_eq!(error.kind(), ErrorKind::ArgumentConflict);
    }

    #[test]
    fn argparse_wrappers_expand_only_unique_long_option_prefixes() {
        for (arguments, expanded) in [
            (
                vec!["hap", "germline", "--threa=3"],
                vec!["hap", "germline", "--threads=3"],
            ),
            (
                vec!["hap", "somatic", "--count-filtered-f"],
                vec!["hap", "somatic", "--count-filtered-fn"],
            ),
            (
                vec!["hap", "prepy", "--wind", "12"],
                vec!["hap", "prepy", "--window-size", "12"],
            ),
            (
                vec!["hap", "ftxpy", "--feature-l", "reference"],
                vec!["hap", "ftxpy", "--feature-label", "reference"],
            ),
            (
                vec!["hap", "qfy", "--report-p", "result"],
                vec!["hap", "qfy", "--report-prefix", "result"],
            ),
        ] {
            let normalized = normalize_legacy_arguments(
                arguments
                    .into_iter()
                    .map(std::ffi::OsString::from)
                    .collect(),
            );
            assert_eq!(
                normalized,
                expanded
                    .into_iter()
                    .map(std::ffi::OsString::from)
                    .collect::<Vec<_>>()
            );
        }

        let ambiguous = ["hap", "pre", "--f"]
            .into_iter()
            .map(std::ffi::OsString::from)
            .collect::<Vec<_>>();
        assert_eq!(normalize_legacy_arguments(ambiguous.clone()), ambiguous);

        let boost = ["hap", "vcfcheck", "--output-f", "result.json"]
            .into_iter()
            .map(std::ffi::OsString::from)
            .collect::<Vec<_>>();
        assert_eq!(normalize_legacy_arguments(boost.clone()), boost);
    }

    #[test]
    fn qfy_alias_honours_legacy_count_switches() {
        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
        ])
        .expect("qfy defaults should parse");
        let Command::Quantify(args) = cli.command else {
            panic!("qfy should resolve to quantify");
        };
        assert!(
            args.write_counts,
            "legacy qfy writes count tables by default"
        );

        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--no-write-counts",
        ])
        .expect("qfy legacy spelling should parse");

        let Command::Quantify(args) = cli.command else {
            panic!("qfy should resolve to quantify");
        };
        assert!(!args.write_counts);

        let cli = Cli::try_parse_from([
            "hap",
            "quantify",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--write-counts",
        ])
        .expect("write-counts should remain accepted");
        let Command::Quantify(args) = cli.command else {
            panic!("quantify should parse");
        };
        assert!(args.write_counts);

        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--write-counts",
            "--no-write-counts",
        ])
        .expect("legacy qfy accepts opposite count switches");
        let Command::Quantify(args) = cli.command else {
            panic!("qfy should resolve to quantify");
        };
        assert!(!args.write_counts, "the final negative switch must win");

        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--no-write-counts",
            "--write-counts",
        ])
        .expect("legacy qfy accepts the reverse switch order");
        let Command::Quantify(args) = cli.command else {
            panic!("qfy should resolve to quantify");
        };
        assert!(args.write_counts, "the final positive switch must win");
    }

    #[test]
    fn qfy_parses_supported_region_and_roc_controls() {
        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "-t",
            "xcmp",
            "-f",
            "conf.bed",
            "--stratification",
            "regions.tsv",
            "--stratification-region",
            "LOW_COMPLEXITY:low-complexity.bed",
            "--stratification-fixchr",
            "--roc",
            "QUAL",
            "--roc-regions",
            "*",
            "--roc-delta",
            "0.5",
            "--ci-alpha",
            "0",
        ])
        .expect("legacy qfy controls should parse");
        let Command::Quantify(args) = cli.command else {
            panic!("qfy should resolve to quantify");
        };
        assert_eq!(args.annotation_type.as_deref(), Some("xcmp"));
        assert_eq!(args.fp_bedfile.as_deref(), Some("conf.bed"));
        assert_eq!(args.strat_tsv.as_deref(), Some("regions.tsv"));
        assert_eq!(
            args.strat_regions,
            vec!["LOW_COMPLEXITY:low-complexity.bed"]
        );
        assert!(args.strat_fixchr);
        assert!(args.do_roc);
        assert_eq!(args.roc, "QUAL");
        assert_eq!(args.roc_regions, vec!["*", "*"]);
        assert_eq!(args.roc_delta, 0.5);
        assert_eq!(args.ci_alpha, 0.0);

        let cli = Cli::try_parse_from([
            "hap",
            "quantify",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--no-roc",
        ])
        .expect("--no-roc should parse");
        let Command::Quantify(args) = cli.command else {
            panic!("quantify should parse");
        };
        assert!(!args.do_roc);
    }

    #[test]
    fn qfy_accepts_standalone_legacy_runtime_controls() {
        let cli = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--output-vtc",
            "--preserve-info",
            "--adjust-conf-regions",
            "truth.vcf.gz",
            "--threads",
            "3",
            "--bcf",
            "--logfile",
            "qfy.log",
            "--verbose",
            "--force-interactive",
        ])
        .expect("standalone qfy compatibility controls should parse");
        let Command::Quantify(args) = cli.command else {
            panic!("qfy should resolve to quantify");
        };
        assert!(args.output_vtc);
        assert!(args.preserve_info);
        assert_eq!(args.adjust_conf_regions.as_deref(), Some("truth.vcf.gz"));
        assert_eq!(args.threads, Some(3));
        assert!(args.bcf);
        assert_eq!(args.logfile.as_deref(), Some("qfy.log"));
        assert!(args.verbose);
        assert!(!args.quiet);
        assert!(args.force_interactive);

        let conflict = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--verbose",
            "--quiet",
        ])
        .expect_err("legacy verbosity controls are mutually exclusive");
        assert_eq!(conflict.kind(), ErrorKind::ArgumentConflict);
    }

    #[test]
    fn qfy_rejects_unknown_annotation_type_during_parsing() {
        let error = Cli::try_parse_from([
            "hap",
            "qfy",
            "annotated.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--type",
            "unknown",
        ])
        .expect_err("unknown annotation formats must not be accepted inertly");
        assert_eq!(error.kind(), ErrorKind::InvalidValue);
        assert_eq!(error.exit_code(), 2);
    }

    #[test]
    fn vcfcheck_alias_accepts_legacy_output_and_location_forms() {
        let cli = Cli::try_parse_from([
            "hap",
            "vcfcheck",
            "input.vcf.gz",
            "--output-file",
            "counts.json",
            "--location",
            "chr1:10-20",
        ])
        .expect("vcfcheck legacy flags should parse");

        let Command::Validate(args) = cli.command else {
            panic!("vcfcheck should resolve to validate");
        };
        assert_eq!(args.output_json.as_deref(), Some("counts.json"));
        assert_eq!(args.locations.as_deref(), Some("chr1:10-20"));
    }

    #[test]
    fn validate_accepts_explicit_input_and_legacy_value_boole() {
        let cli = Cli::try_parse_from([
            "hap",
            "validate",
            "--input-file",
            "input.vcf.gz",
            "--apply-filters",
            "false",
            "--limit-records",
            "7",
            "--message-every",
            "3",
            "-H",
            "true",
            "--check-bcf-errors",
            "false",
            "-W",
            "true",
        ])
        .expect("legacy value-bearing options should parse");

        let Command::Validate(args) = cli.command else {
            panic!("validate should parse");
        };
        assert_eq!(args.input, "input.vcf.gz");
        assert!(!args.apply_filters);
        assert_eq!(args.limit_records, Some(7));
        assert_eq!(args.message_every, Some(3));
        assert!(args.strict_homref);
        assert!(!args.check_bcf_errors);
        assert!(args.all_warnings);
    }

    #[test]
    fn validate_rejects_ambiguous_inputs_and_bare_legacy_boole() {
        let conflict = Cli::try_parse_from([
            "hap",
            "validate",
            "positional.vcf",
            "--input-file",
            "option.vcf",
        ])
        .expect_err("only one input form may be supplied");
        assert_eq!(conflict.kind(), ErrorKind::ArgumentConflict);

        for option in [
            "--apply-filters",
            "--strict-homref",
            "--check-bcf-errors",
            "--all-warnings",
        ] {
            let error = Cli::try_parse_from(["hap", "validate", "input.vcf", option])
                .expect_err("legacy booleans require an explicit value");
            assert_ne!(error.kind(), ErrorKind::UnknownArgument, "{option}");
            assert_eq!(error.exit_code(), 2, "{option}");
        }
    }

    #[test]
    fn germline_accepts_legacy_engine_controls_for_runtime_dispatch() {
        let cli = Cli::try_parse_from([
            "hap",
            "germline",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-r",
            "ref.fa",
            "-o",
            "report",
            "--engine",
            "vcfeval",
            "--engine-vcfeval-path",
            "custom-rtg",
            "--engine-vcfeval-template",
            "template.sdf",
            "--lose-match-distance",
            "42",
        ])
        .expect("legacy engine controls should parse before runtime dispatch");
        let Command::Germline(args) = cli.command else {
            panic!("germline should parse");
        };
        assert_eq!(args.engine, CompareEngine::Vcfeval);
        assert_eq!(args.engine_vcfeval.as_deref(), Some("custom-rtg"));
        assert_eq!(
            args.engine_vcfeval_template.as_deref(),
            Some("template.sdf")
        );
        assert_eq!(args.engine_scmp_distance, 42);
    }

    #[test]
    fn germline_accepts_legacy_preprocessing_switches_and_last_override_wins() {
        let cli = Cli::try_parse_from([
            "hap",
            "germline",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-r",
            "ref.fa",
            "-o",
            "report",
            "--preprocess-truth",
            "-L",
            "--no-leftshift",
            "--decompose",
            "-D",
        ])
        .expect("legacy germline preprocessing switches should parse");
        let Command::Germline(args) = cli.command else {
            panic!("germline should parse");
        };
        assert!(args.preprocess_truth);
        assert!(!args.leftshift);
        assert!(args.no_leftshift);
        assert!(!args.decompose);
        assert!(args.no_decompose);

        let cli = Cli::try_parse_from([
            "hap",
            "germline",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-r",
            "ref.fa",
            "-o",
            "report",
            "--no-leftshift",
            "--leftshift",
            "--no-decompose",
            "--decompose",
        ])
        .expect("later positive preprocessing switches should re-enable behavior");
        let Command::Germline(args) = cli.command else {
            panic!("germline should parse");
        };
        assert!(args.leftshift);
        assert!(!args.no_leftshift);
        assert!(args.decompose);
        assert!(!args.no_decompose);
    }

    #[test]
    fn pre_accepts_remaining_legacy_controls_and_last_override_wins() {
        let cli = Cli::try_parse_from([
            "hap",
            "pre",
            "input.vcf.gz",
            "output.vcf.gz",
            "--no-leftshift",
            "--leftshift",
            "--no-decompose",
            "--decompose",
            "--filters-only",
            "LowQual,q10",
            "--gender",
            "male",
            "--window-size",
            "4096",
            "--bcf",
            "--bcftools-norm",
            "--logfile",
            "pre.log",
            "--verbose",
            "--force-interactive",
            "--no-fixchr",
            "--fixchr",
        ])
        .expect("remaining legacy pre controls should parse without an explicit reference");
        let Command::Preprocess(args) = cli.command else {
            panic!("pre should parse");
        };
        assert!(args.leftshift);
        assert!(!args.no_leftshift);
        assert!(args.decompose);
        assert!(!args.no_decompose);
        assert_eq!(args.filters_only.as_deref(), Some("LowQual,q10"));
        assert_eq!(args.gender, PreprocessGender::Male);
        assert_eq!(args.window_size, 4096);
        assert!(args.bcf);
        assert!(args.bcftools_norm);
        assert_eq!(args.logfile.as_deref(), Some("pre.log"));
        assert!(args.verbose);
        assert!(!args.quiet);
        assert!(args.force_interactive);
        assert_eq!(args.fixchr, Some(true));
        assert!(!args.no_fixchr);
        assert!(args.reference.is_none());

        let cli = Cli::try_parse_from([
            "hap",
            "pre",
            "input.vcf.gz",
            "output.vcf.gz",
            "-L",
            "--no-leftshift",
            "--decompose",
            "-D",
            "--fixchr",
            "--no-fixchr",
        ])
        .expect("later negative switches should disable preprocessing");
        let Command::Preprocess(args) = cli.command else {
            panic!("pre should parse");
        };
        // The positive fields retain their default-on value; the explicit
        // negative marker is the last-wins runtime override.
        assert!(args.leftshift);
        assert!(args.no_leftshift);
        assert!(args.decompose);
        assert!(args.no_decompose);
        assert!(args.fixchr.is_none());
        assert!(args.no_fixchr);

        let error = Cli::try_parse_from([
            "hap",
            "pre",
            "input.vcf.gz",
            "output.vcf.gz",
            "--verbose",
            "--quiet",
        ])
        .expect_err("legacy verbosity controls are mutually exclusive");
        assert_eq!(error.kind(), ErrorKind::ArgumentConflict);

        for arguments in [
            vec![
                "hap",
                "germline",
                "truth.vcf",
                "query.vcf",
                "-o",
                "report",
                "--fixchr",
                "false",
            ],
            vec!["hap", "pre", "input.vcf", "output.vcf", "--fixchr", "false"],
            vec![
                "hap",
                "somatic",
                "truth.vcf",
                "query.vcf",
                "-o",
                "report",
                "--fixchr-truth",
                "false",
            ],
        ] {
            Cli::try_parse_from(arguments)
                .expect_err("legacy fixchr switches do not accept Boolean values");
        }
    }

    #[test]
    fn ftx_accepts_repeatable_bams_for_depth_normalization() {
        let cli = Cli::try_parse_from([
            "hap",
            "ftx",
            "input.vcf.gz",
            "-o",
            "features",
            "-r",
            "ref.fa",
            "--bam",
            "normal.bam",
            "--bam",
            "tumor.bam",
        ])
        .expect("FTX BAM inputs should parse");
        let Command::Ftx(args) = cli.command else {
            panic!("ftx should parse");
        };
        assert_eq!(args.bams, ["normal.bam", "tumor.bam"]);
    }

    #[test]
    fn ftx_accepts_scratch_options_as_a_declared_extension() {
        let cli = Cli::try_parse_from([
            "hap",
            "ftx",
            "input.vcf.gz",
            "-o",
            "features",
            "--scratch-prefix",
            ".",
            "--keep-scratch",
        ])
        .expect("ftx scratch options should parse");
        let Command::Ftx(args) = cli.command else {
            panic!("ftx should parse");
        };
        assert_eq!(args.scratch_prefix.as_deref(), Some("."));
        assert!(args.keep_scratch);
    }

    #[test]
    fn ftx_reference_is_optional_at_parse_time_like_legacy() {
        let cli = Cli::try_parse_from(["hap", "ftx", "input.vcf.gz", "-o", "features"])
            .expect("legacy ftx permits omission of --reference");
        let Command::Ftx(args) = cli.command else {
            panic!("ftx should parse");
        };
        assert!(args.reference.is_none());
    }

    #[test]
    fn somatic_accepts_all_implemented_caller_feature_tables() {
        for feature_table in [
            "hcc.strelka.snv",
            "admix.strelka.indel",
            "hcc.mutect.snv",
            "hcc.varscan2.indel",
            "hcc.pisces.snv",
        ] {
            Cli::try_parse_from([
                "hap",
                "somatic",
                "truth.vcf.gz",
                "query.vcf.gz",
                "-o",
                "report",
                "-r",
                "ref.fa",
                "--feature-table",
                feature_table,
            ])
            .unwrap_or_else(|error| panic!("{feature_table} was rejected: {error}"));
        }
    }

    #[test]
    fn somatic_accepts_repeatable_bam_inputs() {
        let cli = Cli::try_parse_from([
            "hap",
            "somatic",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
            "--bam",
            "normal.bam",
            "--bam",
            "tumor.bam",
        ])
        .expect("somatic BAM inputs should parse");
        let Command::Somatic(args) = cli.command else {
            panic!("somatic should parse");
        };
        assert_eq!(args.bams, ["normal.bam", "tumor.bam"]);
    }

    #[test]
    fn ftx_normalize_and_fix_chr_switches_are_observable() {
        let cli = Cli::try_parse_from([
            "hap",
            "ftx",
            "input.vcf.gz",
            "-o",
            "features",
            "-r",
            "ref.fa",
            "--normalize",
            "--fix-chr",
        ])
        .expect("supported ftx preprocessing switches should parse");
        let Command::Ftx(args) = cli.command else {
            panic!("ftx should resolve to feature extraction");
        };
        assert!(args.normalize);
        assert!(args.fixchr);
    }

    #[test]
    fn germline_accepts_legacy_quantifier_type_hint() {
        let cli = Cli::try_parse_from([
            "hap",
            "germline",
            "truth.vcf",
            "query.vcf",
            "-o",
            "result",
            "-t",
            "ga4gh",
        ])
        .expect("legacy inherited --type should parse");
        let Command::Germline(args) = cli.command else {
            panic!("germline should parse");
        };
        assert_eq!(args.annotation_type.as_deref(), Some("ga4gh"));
    }

    #[test]
    fn somatic_and_set_gt_share_argparse_last_token_wins_semantics() {
        for command in ["germline", "pre"] {
            let prefix = if command == "germline" {
                vec!["hap", command, "truth.vcf", "query.vcf", "-o", "result"]
            } else {
                vec!["hap", command, "input.vcf", "output.vcf"]
            };

            let mut somatic_last = prefix.clone();
            somatic_last.extend(["--set-gt", "first", "--somatic"]);
            let parsed = Cli::try_parse_from(somatic_last).expect("somatic-last ordering parses");
            match parsed.command {
                Command::Germline(args) => {
                    assert!(args.somatic);
                    assert_eq!(args.set_gt, None);
                }
                Command::Preprocess(args) => {
                    assert!(args.somatic);
                    assert_eq!(args.set_gt, None);
                }
                _ => unreachable!(),
            }

            let mut set_gt_last = prefix;
            set_gt_last.extend(["--somatic", "--set-gt", "first"]);
            let parsed = Cli::try_parse_from(set_gt_last).expect("set-gt-last ordering parses");
            match parsed.command {
                Command::Germline(args) => {
                    assert!(!args.somatic);
                    assert_eq!(args.set_gt, Some(SomaticGtMode::First));
                }
                Command::Preprocess(args) => {
                    assert!(!args.somatic);
                    assert_eq!(args.set_gt, Some(SomaticGtMode::First));
                }
                _ => unreachable!(),
            }
        }
    }
}
