//! Command-level regression tests.

#[cfg(test)]
mod tests {
    use super::super::*;
    use crate::cli_compat::cli::{Cli, Command};
    use clap::Parser;

    fn parsed_somatic(extra: &[&str]) -> SomaticArgs {
        let mut argv = vec![
            "hap",
            "somatic",
            "truth.vcf.gz",
            "query.vcf.gz",
            "-o",
            "report",
            "-r",
            "ref.fa",
        ];
        argv.extend_from_slice(extra);
        let cli = Cli::try_parse_from(argv).expect("test arguments should parse");
        let Command::Somatic(args) = cli.command else {
            panic!("somatic command expected");
        };
        args
    }

    fn interval(start: usize, end: usize, label: &str) -> AmbiguousInterval {
        AmbiguousInterval {
            interval: Interval {
                chrom: "chr1".to_string(),
                start,
                end,
            },
            label: label.to_string(),
            details: Vec::new(),
        }
    }

    fn raw_record(line: &str) -> RawVcfRecord {
        RawVcfRecord::from_line(line, Path::new("test.vcf")).expect("valid test VCF record")
    }

    fn unique_test_dir(label: &str) -> PathBuf {
        let id = SOMATIC_SCRATCH_RUN_ID.fetch_add(1, Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "hap-somatic-test-{label}-{}-{id}",
            std::process::id()
        ))
    }

    fn write_test_vcf(path: &Path, position: usize) {
        fs::write(
            path,
            format!(
                "##fileformat=VCFv4.1\n##contig=<ID=chr1,length=20>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t{position}\t.\tA\tC\t.\tPASS\t.\n"
            ),
        )
        .expect("write test VCF");
    }

    #[test]
    fn count_unk_without_fp_regions_classifies_unmatched_calls_as_unknown() {
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &[], true, false),
            QueryClass::Unk
        );
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &[], false, false),
            QueryClass::Fp
        );
    }

    #[test]
    fn explicit_negative_toggle_wins() {
        assert!(resolve_toggle(true, false));
        assert!(!resolve_toggle(true, true));
        assert!(!resolve_toggle(false, false));
    }

    #[test]
    fn paired_somatic_toggles_use_last_token_wins_precedence() {
        let args = parsed_somatic(&["--count-unk", "--no-count-unk"]);
        assert!(!resolve_toggle(args.count_unk, args.no_count_unk));
        let args = parsed_somatic(&["--no-count-unk", "--count-unk"]);
        assert!(resolve_toggle(args.count_unk, args.no_count_unk));

        let args = parsed_somatic(&["--ambi-fp", "--no-ambi-fp"]);
        assert!(!resolve_toggle(args.ambi_fp, args.no_ambi_fp));
        let args = parsed_somatic(&["--no-ambi-fp", "--ambi-fp"]);
        assert!(resolve_toggle(args.ambi_fp, args.no_ambi_fp));
    }

    #[test]
    fn fixchr_pairs_use_last_token_wins_precedence() {
        let args = parsed_somatic(&["--no-fixchr-truth", "--fixchr-truth"]);
        assert!(args.fixchr_truth.unwrap_or(true) && !args.no_fixchr_truth);
        let args = parsed_somatic(&["--fixchr-truth", "--no-fixchr-truth"]);
        assert!(!args.fixchr_truth.unwrap_or(true) || args.no_fixchr_truth);

        let args = parsed_somatic(&["--no-fixchr-query", "--fixchr-query"]);
        assert!(args.fixchr_query.unwrap_or(true) && !args.no_fixchr_query);
        let args = parsed_somatic(&["--fixchr-query", "--no-fixchr-query"]);
        assert!(!args.fixchr_query.unwrap_or(true) || args.no_fixchr_query);
    }

    #[test]
    fn normalize_truth_query_and_all_select_the_legacy_inputs() {
        let truth = parsed_somatic(&["--normalize-truth"]);
        assert_eq!(selected_normalizations(&truth), (true, false));
        let query = parsed_somatic(&["--normalize-query"]);
        assert_eq!(selected_normalizations(&query), (false, true));
        let mut all = parsed_somatic(&["--normalize-all"]);
        all.count_filtered_fn = true;
        assert_eq!(
            selected_normalizations(&all),
            (true, true),
            "-FN reporting must remain independent of -N normalization"
        );
    }

    #[test]
    fn governed_somatic_defaults_remain_supported() {
        validate_args(&parsed_somatic(&[])).expect("default comparison must remain supported");
        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "hcc.strelka.indel",
        ]))
        .expect("nf-test feature-table comparison must remain supported");
        validate_args(&parsed_somatic(&["--feature-table", "admix.strelka.snv"]))
            .expect("caller-specific SNV comparison must remain supported");
    }

    #[test]
    fn transformation_and_reporting_controls_are_accepted() {
        let cases: &[&[&str]] = &[
            &["--normalize-truth"],
            &["--normalize-query"],
            &["--normalize-all"],
            &["--no-fixchr-truth"],
            &["--no-fixchr-query"],
            &["--roc", "strelka.snv"],
        ];
        for case in cases {
            validate_args(&parsed_somatic(case))
                .unwrap_or_else(|error| panic!("legacy control {case:?} was rejected: {error}"));
        }
    }

    #[test]
    fn scratch_lifecycle_logging_and_verbosity_are_operational() {
        let root = unique_test_dir("controls");
        let scratch = root.join("explicit-scratch");
        let logfile = root.join("som.log");
        fs::create_dir_all(&root).expect("create operational test root");

        let mut args = parsed_somatic(&[]);
        args.scratch_prefix = Some(scratch.display().to_string());
        args.logfile = Some(logfile.display().to_string());
        args.verbose = true;
        let mut controls = SomaticOperationalControls::prepare(&args).expect("prepare controls");
        assert!(controls.scratch.path.is_dir());
        assert!(!controls.should_print_summary());
        controls.info("operational log marker").expect("write log");
        controls.scratch.cleanup().expect("retain explicit scratch");
        assert!(scratch.is_dir());
        assert!(
            fs::read_to_string(&logfile)
                .expect("read logfile")
                .contains("operational log marker")
        );

        let mut quiet_args = parsed_somatic(&[]);
        quiet_args.quiet = true;
        let quiet = SomaticOperationalControls::prepare(&quiet_args).expect("prepare quiet mode");
        assert!(!quiet.should_print_summary());
        let quiet_path = quiet.scratch.path.clone();
        quiet.scratch.cleanup().expect("clean quiet scratch");
        assert!(!quiet_path.exists());

        let mut keep_args = parsed_somatic(&[]);
        keep_args.keep_scratch = true;
        let keep = SomaticOperationalControls::prepare(&keep_args).expect("prepare kept scratch");
        let keep_path = keep.scratch.path.clone();
        keep.scratch.cleanup().expect("keep scratch");
        assert!(keep_path.is_dir());

        fs::remove_dir_all(&keep_path).expect("remove retained generated scratch");
        fs::remove_dir_all(&root).expect("remove operational test root");
    }

    #[test]
    fn continue_reuses_cached_normalized_inputs() {
        let root = unique_test_dir("continue");
        let scratch = root.join("scratch");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let reference = root.join("reference.fa");
        fs::create_dir_all(&root).expect("create continue test root");
        fs::write(&reference, ">chr1\nAAAAAAAAAAAAAAAAAAAA\n").expect("write reference");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 7);

        let comparison = |output: &Path, cont: bool| {
            let mut args = parsed_somatic(&[]);
            args.truth = truth.display().to_string();
            args.query = query.display().to_string();
            args.reference = reference.display().to_string();
            args.output = output.display().to_string();
            args.scratch_prefix = Some(scratch.display().to_string());
            args.cont = cont;
            args.quiet = true;
            run(args).expect("run somatic comparison");
        };

        comparison(&root.join("first"), false);
        assert!(scratch.join("normalized_truth.vcf.gz").is_file());
        assert!(scratch.join("normalized_query.vcf.gz").is_file());

        write_test_vcf(&query, 8);
        comparison(&root.join("continued"), true);
        let stats =
            fs::read_to_string(root.join("continued.stats.csv")).expect("read continued stats");
        let snv = stats
            .lines()
            .find(|line| line.starts_with("1,SNVs,"))
            .expect("continued SNV row");
        assert!(snv.starts_with("1,SNVs,1,1,1,0,0,0,0,"));

        fs::remove_dir_all(&root).expect("remove continue test root");
    }

    #[test]
    fn explain_ambiguous_without_features_writes_csv_and_metrics_tables() {
        let root = unique_test_dir("explain-no-features");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let ambiguous = root.join("ambiguous.bed");
        fs::create_dir_all(&root).expect("create explanation test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 8);
        fs::write(
            &ambiguous,
            "chr1\t7\t8\tignored\tunk\t2\tignored\tlow-vaf\n",
        )
        .expect("write ambiguous BED");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.ambiguous_beds = vec![ambiguous.display().to_string()];
        args.explain_ambiguous = true;
        args.fp_region_size = Some("10".to_string());
        args.quiet = true;
        run(args).expect("explanation run must not require a reference or feature table");

        assert!(!root.join("result.features.csv").exists());
        assert!(
            fs::read_to_string(root.join("result.ambiclasses.csv"))
                .expect("read ambiguity classes")
                .contains("ambi-unk,1")
        );
        assert!(
            fs::read_to_string(root.join("result.ambireasons.csv"))
                .expect("read ambiguity reasons")
                .contains("ambi-unk: low-vaf,1")
        );
        let metrics =
            fs::read_to_string(root.join("result.metrics.json")).expect("read explanation metrics");
        assert!(metrics.contains("\"id\": \"ambiclasses\""));
        assert!(metrics.contains("\"id\": \"ambireasons\""));
        assert!(metrics.contains("ambi-unk: low-vaf"));

        fs::remove_dir_all(&root).expect("remove explanation test root");
    }

    #[test]
    fn empty_ambiguity_explanations_do_not_create_detail_tables() {
        let root = unique_test_dir("empty-explanation");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let ambiguous = root.join("ambiguous.bed");
        fs::create_dir_all(&root).expect("create empty explanation test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 7);
        fs::write(
            &ambiguous,
            "chr1\t7\t8\tignored\tunk\t2\tignored\tlow-vaf\n",
        )
        .expect("write ambiguous BED");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.ambiguous_beds = vec![ambiguous.display().to_string()];
        args.explain_ambiguous = true;
        args.fp_region_size = Some("10".to_string());
        args.quiet = true;
        run(args).expect("empty explanation comparison must succeed");

        assert!(!root.join("result.ambiclasses.csv").exists());
        assert!(!root.join("result.ambireasons.csv").exists());
        let metrics =
            fs::read_to_string(root.join("result.metrics.json")).expect("read explanation metrics");
        assert!(!metrics.contains("\"id\": \"ambiclasses\""));
        assert!(!metrics.contains("\"id\": \"ambireasons\""));

        fs::remove_dir_all(&root).expect("remove empty explanation test root");
    }

    #[test]
    fn explicit_fp_bed_participates_in_explanations_and_uppercase_ambiguous_denominator() {
        let root = unique_test_dir("fp-explanation-denominator");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let fp = root.join("fp.bed");
        let ambiguous = root.join("ambiguous.bed");
        fs::create_dir_all(&root).expect("create FP explanation test root");
        write_test_vcf(&truth, 2);
        fs::write(
            &query,
            concat!(
                "##fileformat=VCFv4.1\n",
                "##contig=<ID=chr1,length=20>\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "chr1\t2\t.\tA\tC\t.\tPASS\t.\n",
                "chr1\t5\t.\tA\tG\t.\tPASS\t.\n",
                "chr1\t8\t.\tA\tT\t.\tPASS\t.\n",
            ),
        )
        .expect("write FP explanation query VCF");
        fs::write(&fp, "chr1\t4\t5\t2\tFP\tsource=upper\n").expect("write FP BED");
        fs::write(
            &ambiguous,
            concat!(
                "chr1\t4\t5\t2\tFP\tsource=upper\n",
                "chr1\t7\t8\t3\tfp\tsource=lower\n",
            ),
        )
        .expect("write ambiguity BED");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.fp_bedfile = Some(fp.display().to_string());
        args.ambiguous_beds = vec![ambiguous.display().to_string()];
        args.explain_ambiguous = true;
        args.quiet = true;
        run(args).expect("labeled FP regions must supply the automatic denominator");

        let stats = fs::read_to_string(root.join("result.stats.csv")).expect("read stats");
        assert!(
            stats
                .lines()
                .any(|line| line.starts_with("5,records,") && line.contains(",2,500000.0,")),
            "explicit and uppercase ambiguous FP intervals both count toward the denominator"
        );
        let classes =
            fs::read_to_string(root.join("result.ambiclasses.csv")).expect("read classes");
        assert!(classes.contains(",FP,1\n"));
        assert!(classes.contains(",ambi-fp,1\n"));
        let reasons =
            fs::read_to_string(root.join("result.ambireasons.csv")).expect("read reasons");
        assert!(
            reasons.contains(",FP: rep. count 2,2\n"),
            "the explicit and ambiguous FP entries both contribute reasons"
        );
        assert!(reasons.contains(",ambi-fp: rep. count 3,1\n"));

        fs::remove_dir_all(&root).expect("remove FP explanation test root");
    }

    #[test]
    fn automatic_reference_denominator_uses_truth_contigs_only() {
        let root = unique_test_dir("truth-only-denominator");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let reference = root.join("reference.fa");
        fs::create_dir_all(&root).expect("create truth denominator test root");
        write_test_vcf(&truth, 2);
        fs::write(
            &query,
            concat!(
                "##fileformat=VCFv4.1\n",
                "##contig=<ID=chr1,length=10>\n",
                "##contig=<ID=chr2,length=20>\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "chr1\t2\t.\tA\tC\t.\tPASS\t.\n",
                "chr2\t5\t.\tC\tG\t.\tPASS\t.\n",
            ),
        )
        .expect("write query-only contig VCF");
        fs::write(
            &reference,
            ">chr1\nAAAAAAAAAA\n>chr2\nCCCCCCCCCCCCCCCCCCCC\n",
        )
        .expect("write reference");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = reference.display().to_string();
        args.quiet = true;
        run(args).expect("run truth-only automatic denominator comparison");

        let stats = fs::read_to_string(root.join("result.stats.csv")).expect("read stats");
        assert!(
            stats
                .lines()
                .any(|line| line.starts_with("5,records,") && line.contains(",10,100000.0,")),
            "the query-only chr2 contig must not enlarge the reference denominator"
        );

        fs::remove_dir_all(&root).expect("remove truth denominator test root");
    }

    #[test]
    fn usable_fp_bed_avoids_loading_a_missing_reference() {
        let root = unique_test_dir("lazy-reference-fp-bed");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let fp = root.join("fp.bed");
        fs::create_dir_all(&root).expect("create lazy-reference test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 8);
        fs::write(&fp, "chr1\t0\t20\n").expect("write FP BED");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.fp_bedfile = Some(fp.display().to_string());
        args.quiet = true;
        run(args).expect("usable FP BED must avoid loading the missing reference");

        let stats = fs::read_to_string(root.join("result.stats.csv")).expect("read stats");
        let snv = stats
            .lines()
            .find(|line| line.starts_with("1,SNVs,"))
            .expect("SNV row");
        assert!(snv.contains(",20,50000.0,"));

        fs::remove_dir_all(&root).expect("remove lazy-reference test root");
    }

    #[test]
    fn normalize_all_with_filtered_fn_still_requires_and_uses_reference() {
        let root = unique_test_dir("normalize-all-filtered-fn");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        fs::create_dir_all(&root).expect("create normalization test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 7);

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.normalize_all = true;
        args.count_filtered_fn = true;
        args.include_nonpass = true;
        args.feature_table = Some("generic".to_string());
        args.fp_region_size = Some("10".to_string());
        args.quiet = true;
        let error = run(args).expect_err("normalization must load the missing reference");
        assert!(error.to_string().contains("failed to read"));

        fs::remove_dir_all(&root).expect("remove normalization test root");
    }

    #[test]
    fn af_controls_match_legacy_validation_independently_of_happy_stats() {
        validate_args(&parsed_somatic(&[
            "--feature-table",
            "hcc.strelka.indel",
            "--bin-afs",
        ]))
        .expect("AF stats rows do not require --happy-stats");

        validate_args(&parsed_somatic(&[
            "--af-binsize",
            "0.1",
            "--af-truth",
            "TRUTH_AF",
            "--af-query",
            "QUERY_AF",
        ]))
        .expect("AF controls are inert unless --bin-afs is selected");

        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "generic",
            "--happy-stats",
            "--bin-afs",
        ]))
        .expect("legacy validates selected AF columns after extracting the feature table");

        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "hcc.strelka.indel",
            "--happy-stats",
            "--bin-afs",
            "--af-query",
            "MISSING_AF",
        ]))
        .expect("legacy validates missing AF selectors after input I/O");

        validate_args(&parsed_somatic(&[
            "-P",
            "--feature-table",
            "hcc.strelka.indel",
            "--happy-stats",
            "--bin-afs",
            "--af-binsize",
            "0.25",
            "--af-truth",
            "I.T_ALT_RATE",
            "--af-query",
            "T_AF",
        ]))
        .expect("implemented AF extended-summary controls should remain supported");
    }

    #[test]
    fn af_bin_edge_values_follow_the_pinned_python_loop() {
        for raw in ["0", "-0.1", "nan", "inf"] {
            let args = parsed_somatic(&["--af-binsize", raw]);
            validate_args(&args)
                .unwrap_or_else(|error| panic!("legacy AF bin {raw:?} was rejected: {error}"));
        }
        for raw in ["", "bogus"] {
            let mut args = parsed_somatic(&[]);
            args.af_strat_binsize = raw.to_string();
            assert!(validate_args(&args).is_err());
        }
        let tiny = parsed_somatic(&["--af-binsize", "1e-12"]);
        assert!(
            validate_args(&tiny)
                .unwrap_err()
                .to_string()
                .contains("more than 100 bins")
        );

        assert!(parse_af_bins("0").unwrap().is_empty());
        assert!(parse_af_bins("-0.1").unwrap().is_empty());
        let nan = parse_af_bins("nan").unwrap();
        assert_eq!(nan.len(), 1);
        assert_eq!(nan[0].0, 0.0);
        assert!(nan[0].1.is_nan());
        assert_eq!(format_af_interval(nan[0].0, nan[0].1), "0.000000-nan");
        assert_eq!(parse_af_bins("inf").unwrap(), vec![(0.0, 1.000_000_01)]);
        assert_eq!(format_af_interval(0.0, 1.000_000_01), "0.000000-1.000000");
    }

    #[test]
    fn af_roc_artifacts_are_declared_in_the_transaction_plan() {
        let mut args = parsed_somatic(&["--bin-afs", "--af-binsize", "0.5"]);
        args.roc = Some("generic".to_string());
        let artifacts = somatic_artifacts(&args).unwrap();
        for prefix in ["records", "SNVs", "indels"] {
            for interval in ["0.000000-0.500000", "0.500000-1.000000"] {
                assert!(artifacts.contains(&format!("{prefix}.{interval}.roc.csv")));
            }
        }
    }

    #[test]
    fn af_bin_edge_values_emit_the_legacy_stats_rows() {
        let root = unique_test_dir("af-bin-edges");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        fs::create_dir_all(&root).expect("create AF edge test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 8);

        for (label, value, expected_interval) in [
            ("zero", "0", None),
            ("negative", "-0.1", None),
            ("nan", "nan", Some("records.0.000000-nan")),
            ("infinity", "inf", Some("records.0.000000-1.000000")),
        ] {
            let mut args = parsed_somatic(&[]);
            args.truth = truth.display().to_string();
            args.query = query.display().to_string();
            args.output = root.join(label).display().to_string();
            args.reference = root.join("missing.fa").display().to_string();
            args.feature_table = Some("generic".to_string());
            args.af_strat = true;
            args.af_strat_binsize = value.to_string();
            args.af_strat_truth = "QUAL.truth".to_string();
            args.af_strat_query = "QUAL".to_string();
            args.fp_region_size = Some("10".to_string());
            args.quiet = true;
            run(args).unwrap_or_else(|error| panic!("AF edge {value:?} failed: {error}"));

            let stats = fs::read_to_string(root.join(format!("{label}.stats.csv")))
                .expect("read AF edge stats");
            let intervals = stats
                .lines()
                .filter(|line| line.contains("records."))
                .collect::<Vec<_>>();
            match expected_interval {
                Some(expected) => {
                    assert_eq!(intervals.len(), 1);
                    assert!(intervals[0].contains(expected));
                }
                None => assert!(intervals.is_empty()),
            }
        }

        fs::remove_dir_all(&root).expect("remove AF edge test root");
    }

    #[test]
    fn ambiguous_fp_toggle_matches_legacy_label_semantics() {
        let ambiguous = vec![
            interval(10, 20, "fp"),
            interval(30, 40, "Fp"),
            interval(50, 60, "FP"),
        ];
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &ambiguous, true, false),
            QueryClass::Ambi
        );
        assert_eq!(
            classify_query("chr1", 11, 11, &[], &ambiguous, true, true),
            QueryClass::Fp
        );
        assert_eq!(
            classify_query("chr1", 31, 31, &[], &ambiguous, true, true),
            QueryClass::Ambi,
            "legacy labels are case-sensitive"
        );
        assert_eq!(
            classify_query("chr1", 51, 51, &[], &ambiguous, true, false),
            QueryClass::Fp,
            "uppercase FP is unconditional"
        );
    }

    #[test]
    fn ambiguous_bed_uses_the_fifth_column_as_the_label() {
        let mut bed = tempfile::NamedTempFile::new().expect("temporary BED");
        writeln!(bed, "chr1\t0\t10\tignored\tFP\textra").expect("write BED");
        let contigs = BTreeSet::from(["chr1".to_string()]);
        let intervals =
            load_ambiguous_beds(&[bed.path().to_string_lossy().into_owned()], &contigs, true)
                .expect("load ambiguous BED");
        assert_eq!(intervals.len(), 1);
        assert_eq!(intervals[0].label, "FP");
        assert_eq!(intervals[0].details, ["ignored", "FP", "extra"]);
    }

    #[test]
    fn classification_beds_follow_the_truth_fixchr_switch() {
        let root = unique_test_dir("classification-bed");
        let fp = root.join("fp.bed");
        let ambiguous = root.join("ambiguous.bed");
        fs::create_dir_all(&root).expect("create BED test root");
        fs::write(&fp, "1\t0\t10\n").expect("write FP BED");
        fs::write(&ambiguous, "1\t0\t10\tignored\tfp\n").expect("write ambiguous BED");
        let contigs = BTreeSet::from(["chr1".to_string()]);

        assert_eq!(
            load_classification_bed(&fp, &contigs, true).unwrap()[0].chrom,
            "chr1"
        );
        assert_eq!(
            load_classification_bed(&fp, &contigs, false).unwrap()[0].chrom,
            "1"
        );
        assert_eq!(
            load_ambiguous_beds(&[ambiguous.display().to_string()], &contigs, true).unwrap()[0]
                .interval
                .chrom,
            "chr1"
        );
        assert_eq!(
            load_ambiguous_beds(&[ambiguous.display().to_string()], &contigs, false).unwrap()[0]
                .interval
                .chrom,
            "1"
        );

        fs::remove_dir_all(&root).expect("remove BED test root");
    }

    #[test]
    fn ambiguous_explanation_records_classes_and_reasons() {
        let ambiguous = vec![AmbiguousInterval {
            interval: Interval {
                chrom: "chr1".to_string(),
                start: 10,
                end: 20,
            },
            label: "unk".to_string(),
            details: vec![
                "2".to_string(),
                "ignored".to_string(),
                "low-vaf".to_string(),
            ],
        }];
        let mut classes = BTreeMap::new();
        let mut reasons = BTreeMap::new();
        record_ambiguous_explanation(
            "chr1",
            11,
            11,
            &ambiguous,
            false,
            &mut classes,
            &mut reasons,
        );
        assert_eq!(classes.get("ambi-unk"), Some(&1));
        assert_eq!(reasons.get("ambi-unk: rep. count 2"), Some(&1));
        assert_eq!(reasons.get("ambi-unk: low-vaf"), Some(&1));
    }

    #[test]
    fn ambiguity_reason_csv_preserves_python2_counter_indices() {
        let counts = BTreeMap::from([
            ("ambi-unk: 2".to_string(), 1),
            ("ambi-unk: ignored".to_string(), 1),
            ("ambi-unk: low-vaf".to_string(), 1),
            ("ambi-unk: rep. count ignored".to_string(), 1),
        ]);
        let output = tempfile::NamedTempFile::new().expect("temporary ambiguity CSV");
        write_legacy_count_table(output.path(), "reason", &counts)
            .expect("write ambiguity reasons");
        assert_eq!(
            fs::read_to_string(output.path()).expect("read ambiguity reasons"),
            concat!(
                ",reason,count\n",
                "2,ambi-unk: 2,1\n",
                "1,ambi-unk: ignored,1\n",
                "3,ambi-unk: low-vaf,1\n",
                "0,ambi-unk: rep. count ignored,1\n",
            )
        );
        let metric = count_metric_json("ambireasons", "reason", &counts);
        assert!(metric.contains("\"values\": [2, 1, 3, 0]"));
    }

    #[test]
    fn explicit_fp_regions_take_priority_over_ambiguous_regions() {
        let fp = vec![Interval {
            chrom: "chr1".to_string(),
            start: 10,
            end: 20,
        }];
        let ambiguous = vec![interval(10, 20, "unk")];
        assert_eq!(
            classify_query("chr1", 11, 11, &fp, &ambiguous, true, false),
            QueryClass::Fp
        );
    }

    #[test]
    fn raw_filtering_uses_symbolic_end_for_regions_and_start_for_targets() {
        let path = Path::new("query.vcf");
        let record =
            RawVcfRecord::from_line("chr1\t3\t.\tC\t<DEL>\t.\tPASS\tEND=5\tGT\t0/1", path).unwrap();
        let contigs = BTreeSet::from(["chr1".to_string()]);
        let boundary = vec![Interval {
            chrom: "chr1".to_string(),
            start: 4,
            end: 5,
        }];

        let by_region = filter_raw_records(
            vec![record.clone()],
            path,
            &RawFilterOptions {
                reference_contigs: &contigs,
                fixchr: true,
                pass_only: false,
                regions: Some(&boundary),
                targets: None,
                locations: None,
            },
        )
        .unwrap();
        let by_target = filter_raw_records(
            vec![record.clone()],
            path,
            &RawFilterOptions {
                reference_contigs: &contigs,
                fixchr: true,
                pass_only: false,
                regions: None,
                targets: Some(&boundary),
                locations: None,
            },
        )
        .unwrap();
        let by_location = filter_raw_records(
            vec![record],
            path,
            &RawFilterOptions {
                reference_contigs: &contigs,
                fixchr: true,
                pass_only: false,
                regions: None,
                targets: None,
                locations: Some(&[vcf::LocationFilter::Range {
                    chrom: "chr1".to_string(),
                    start: 3,
                    end: 3,
                }]),
            },
        )
        .unwrap();

        assert_eq!(by_region.len(), 1, "-R must use INFO/END");
        assert!(by_target.is_empty(), "-T must use POS");
        assert_eq!(by_location.len(), 1, "-l must use POS");
    }

    #[test]
    fn fp_region_size_with_regions_and_location_preserves_legacy_zero_bug() {
        let fp_regions = vec![Interval {
            chrom: "chr1".to_string(),
            start: 0,
            end: 50,
        }];
        let references = BTreeMap::from([("chr1".to_string(), "A".repeat(100))]);
        let range = [vcf::LocationFilter::Range {
            chrom: "chr1".to_string(),
            start: 21,
            end: 30,
        }];
        let contig = [vcf::LocationFilter::Contig("chr1".to_string())];

        assert_eq!(
            calculate_fp_region_size(None, &fp_regions, &[], Some(&range), &references, &[]),
            0
        );
        assert_eq!(
            calculate_fp_region_size(None, &fp_regions, &[], Some(&contig), &references, &[]),
            0
        );
        assert_eq!(
            calculate_fp_region_size(None, &[], &[], Some(&range), &references, &[]),
            10
        );
        assert_eq!(
            calculate_fp_region_size(Some("7"), &fp_regions, &[], Some(&range), &references, &[],),
            7
        );
        assert!(!fp_region_size_requires_reference(Some("7"), &[], &[]));
        assert!(!fp_region_size_requires_reference(None, &fp_regions, &[]));
        assert!(fp_region_size_requires_reference(None, &[], &[]));
        assert!(fp_region_size_requires_reference(Some("auto"), &[], &[]));
        assert_eq!(fp_rate_or_blank(1, 0), "inf");
        assert_eq!(fp_rate_or_blank(0, 0), "");
        let rates = vec!["inf".to_string(), String::new()];
        assert_eq!(infer_type(&rates), "double");
        assert_eq!(
            column_json("fp.rate", "fp.rate", "double", &rates, false),
            "{\"values\": [null, null], \"type\": \"double\", \"id\": \"fp.rate\", \"label\": \"fp.rate\"}"
        );
    }

    #[test]
    fn fp_bed_with_dash_range_preserves_legacy_late_failure_artifacts() {
        let root = unique_test_dir("fp-range-failure");
        let truth = root.join("truth.vcf");
        let query = root.join("query.vcf");
        let fp = root.join("fp.bed");
        fs::create_dir_all(&root).expect("create FP range failure test root");
        write_test_vcf(&truth, 7);
        write_test_vcf(&query, 8);
        fs::write(&fp, "chr1\t0\t20\t2\tFP\tsource=upper\n").expect("write FP BED");

        let mut args = parsed_somatic(&[]);
        args.truth = truth.display().to_string();
        args.query = query.display().to_string();
        args.output = root.join("result").display().to_string();
        args.reference = root.join("missing.fa").display().to_string();
        args.fp_bedfile = Some(fp.display().to_string());
        args.location = Some("chr1:1-10".to_string());
        args.feature_table = Some("generic".to_string());
        args.quiet = true;

        let error = run(args).expect_err("legacy denominator parsing must reject dash ranges");
        assert!(error.to_string().contains("invalid literal for int()"));
        assert!(
            root.join("result.features.csv").is_file(),
            "som.py writes the feature table before the denominator failure"
        );
        assert!(!root.join("result.stats.csv").exists());
        assert!(!root.join("result.metrics.json").exists());
        validate_legacy_fp_location_denominator(Some("10"), Some("chr1:1-10"), true)
            .expect("an explicit integer denominator bypasses the legacy range bug");

        fs::remove_dir_all(&root).expect("remove FP range failure test root");
    }

    #[test]
    fn fp_classification_uses_the_reference_span_not_symbolic_end() {
        let record = raw_record("chr1\t10\t.\tA\t<DEL>\t.\tPASS\tEND=100");
        let fp = vec![Interval {
            chrom: "chr1".to_string(),
            start: 49,
            end: 50,
        }];
        assert_eq!(
            classify_query(
                &record.chrom,
                record.pos,
                record.end_pos(),
                &fp,
                &[],
                true,
                false,
            ),
            QueryClass::Unk
        );
        assert_eq!(
            classify_query(
                &record.chrom,
                record.pos,
                record.effective_end_pos(Path::new("test.vcf")).unwrap(),
                &fp,
                &[],
                true,
                false,
            ),
            QueryClass::Fp,
            "the symbolic span would have produced the wrong class"
        );
    }

    #[test]
    fn terminal_non_ref_is_removed_only_when_its_allele_is_called() {
        assert!(calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\tC,<NON_REF>\t.\tPASS\t.\tGT\t0/2"
        )));
        assert!(!calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\tC,<NON_REF>\t.\tPASS\t.\tGT\t0/1"
        )));
        assert!(!calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\t<NON_REF>,C\t.\tPASS\t.\tGT\t0/1"
        )));
        assert!(!calls_terminal_non_ref(&raw_record(
            "chr1\t10\t.\tA\tC,<NON_REF>\t.\tPASS\t."
        )));
    }

    #[test]
    fn sites_only_records_survive_raw_filtering() {
        let record = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\t.");
        let filtered = filter_raw_records(
            vec![record],
            Path::new("test.vcf"),
            &RawFilterOptions {
                reference_contigs: &BTreeSet::from(["chr1".to_string()]),
                fixchr: true,
                pass_only: true,
                regions: None,
                targets: None,
                locations: None,
            },
        )
        .expect("filter sites-only VCF");
        assert_eq!(filtered.len(), 1);
        assert!(filtered[0].record.samples.is_empty());
    }

    #[test]
    fn exact_pairing_preserves_duplicate_occurrences() {
        let record = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\t.");
        let filtered = |record: RawVcfRecord| FilteredRawRecord {
            key: vcf::VariantKey {
                chrom: record.chrom.clone(),
                pos: record.pos,
                ref_allele: record.ref_allele.clone(),
                alt_allele: record.alt_allele.clone(),
            },
            record,
        };
        let truth = vec![filtered(record.clone()), filtered(record.clone())];
        let query = vec![filtered(record.clone()), filtered(record)];
        let (truth_matches, query_matches) = pair_exact_records(&truth, &query);
        assert_eq!(truth_matches, vec![Some(0), Some(1)]);
        assert_eq!(query_matches, vec![Some(0), Some(1)]);
    }

    #[test]
    fn caller_feature_merge_suffixes_truth_and_preserves_numeric_csv_types() {
        let truth = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\teditDistance=1\tGT\t0/0\t0/1");
        let query = raw_record(
            "chr1\t10\t.\tA\tC\t60\tPASS\tNT=ref;QSS_NT=10;VQSR=2;SomaticEVS=3;MQ=50;MQ0=1;SNVSB=0.2;ReadPosRankSum=0.3\tGT:DP:FDP:SDP:AU:CU:GU:TU\t0/0:20:2:1:10,0:0,0:0,0:0,0\t0/1:30:3:2:15,0:5,0:0,0:0,0",
        );
        let table = build_caller_feature_table(
            "admix.strelka.snv",
            &[],
            &[],
            None,
            true,
            &CallerRecordGroups {
                tp_truth: &[truth],
                tp_query: &[query],
                fn_truth: &[],
                fp_query: &[],
                ambi_query: &[],
                unk_query: &[],
            },
        )
        .expect("caller feature table");
        let headers = parse_csv_line(&table.header);
        let row = parse_csv_line(&table.tp[0]);
        let field = |name: &str| {
            let index = headers.iter().position(|header| header == name).unwrap();
            row[index].as_str()
        };
        assert_eq!(
            &headers[..8],
            [
                "",
                "CHROM",
                "POS",
                "tag",
                "REF",
                "REF.truth",
                "ALT",
                "ALT.truth"
            ]
        );
        assert_eq!(field("tag"), "TP");
        assert_eq!(field("REF.truth"), "A");
        assert_eq!(field("I.editDistance"), "1.00000000");
        assert_eq!(field("QSS_NT"), "10.00000000");
        assert_eq!(field("S.2.GT"), "0/1");
        assert_eq!(field("POS"), "10");
    }

    #[test]
    fn no_order_check_controls_caller_tp_order_validation() {
        let truth = raw_record("chr1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/0\t0/1");
        let query = raw_record(
            "chr1\t11\t.\tA\tC\t60\tPASS\tNT=ref;QSS_NT=10\tGT:DP:FDP:SDP:AU:CU:GU:TU\t0/0:20:0:0:20,0:0,0:0,0:0,0\t0/1:30:0:0:15,0:15,0:0,0:0,0",
        );
        let groups = CallerRecordGroups {
            tp_truth: std::slice::from_ref(&truth),
            tp_query: std::slice::from_ref(&query),
            fn_truth: &[],
            fp_query: &[],
            ambi_query: &[],
            unk_query: &[],
        };

        let error =
            match build_caller_feature_table("hcc.strelka.snv", &[], &[], None, true, &groups) {
                Ok(_) => panic!("default order check must reject mismatched TP rows"),
                Err(error) => error,
            };
        assert!(error.to_string().contains("out of order"));
        build_caller_feature_table("hcc.strelka.snv", &[], &[], None, false, &groups)
            .expect("--no-order-check must bypass the developer safety check");
    }

    #[test]
    fn bam_depths_flow_into_somatic_caller_features() {
        let fixture_dir =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("verification/assets/fixtures/ftx-bam");
        let vcf_path = fixture_dir.join("input.vcf");
        let (headers, records) = vcf::load_raw_vcf(&vcf_path).expect("load BAM parity VCF");
        let depths =
            ftx::bam_normalization_depths(&[fixture_dir.join("reads.bam").display().to_string()])
                .expect("scan BAM normalization depths");

        let table = build_caller_feature_table(
            "hcc.strelka.snv",
            &headers,
            &headers,
            Some(&depths),
            true,
            &CallerRecordGroups {
                tp_truth: &records,
                tp_query: &records,
                fn_truth: &[],
                fp_query: &[],
                ambi_query: &[],
                unk_query: &[],
            },
        )
        .expect("caller feature table with BAM depths");
        let columns = parse_csv_line(&table.header);
        let row = parse_csv_line(&table.tp[0]);
        let field = |name: &str| {
            let index = columns.iter().position(|column| column == name).unwrap();
            row[index].as_str()
        };

        assert!((depths["chr1"] - 0.9).abs() < f64::EPSILON);
        assert_eq!(field("N_DP_RATE"), "10.00000000");
        assert_eq!(field("T_DP_RATE"), "20.00000000");
    }

    #[test]
    fn raw_stats_types_match_legacy_labels_and_indexes() {
        assert_eq!(STATS_TYPE_ROWS[0], (0, "indels"));
        assert_eq!(STATS_TYPE_ROWS[1], (1, "SNVs"));
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\tC\t.\tPASS\t.")),
            Some("SNVs")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\tAC\t.\tPASS\t.")),
            Some("indels")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tAC\tGT\t.\tPASS\t.")),
            Some("MNPs")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\t<DEL>\t.\tPASS\t.")),
            Some("others")
        );
        assert_eq!(
            raw_type_label(&raw_record("chr1\t1\t.\tA\t.\t.\tPASS\t.")),
            None
        );
    }

    #[test]
    fn confidence_level_changes_interval_without_changing_point_estimate() {
        let counts = SomaticCounts {
            truth_total: 10,
            query_total: 10,
            tp: 8,
            fp: 2,
            fn_count: 2,
            ..SomaticCounts::default()
        };
        let ci_95 = render_row(
            1,
            "SNVs",
            counts,
            &StatsRowContext {
                fp_region_size: 100,
                ci_alpha: 0.05,
                filtered: None,
                include_filtered_columns: false,
                commandline: "som.py",
            },
        );
        let ci_80 = render_row(
            1,
            "SNVs",
            counts,
            &StatsRowContext {
                fp_region_size: 100,
                ci_alpha: 0.20,
                filtered: None,
                include_filtered_columns: false,
                commandline: "som.py",
            },
        );
        let ci_95 = ci_95.split(',').collect::<Vec<_>>();
        let ci_80 = ci_80.split(',').collect::<Vec<_>>();
        assert_eq!(ci_95[9], ci_80[9]);
        assert_ne!(ci_95[10], ci_80[10]);
        assert_ne!(ci_95[11], ci_80[11]);
    }

    #[test]
    fn filtered_fn_columns_align_with_their_header() {
        let counts = SomaticCounts {
            truth_total: 10,
            query_total: 11,
            tp: 8,
            fp: 3,
            fn_count: 2,
            unk: 0,
            ambi: 0,
        };
        let filtered = FilteredCounts {
            tp: 1,
            fp: 2,
            ..FilteredCounts::default()
        };
        let header = stats_header(true).split(',').count();
        let row = render_row(
            1,
            "SNVs",
            counts,
            &StatsRowContext {
                fp_region_size: 100,
                ci_alpha: 0.05,
                filtered: Some(filtered),
                include_filtered_columns: true,
                commandline: "som.py",
            },
        );
        assert_eq!(header, row.split(',').count());
        assert!(row.contains(",2.0,1.0,0.0,0.0,"));
    }

    #[test]
    fn generic_filtered_fn_counts_default_each_reported_variant_type_to_zero() {
        let filtered = filtered_counts_for_type(true, Some("generic"), "indels", &BTreeMap::new())
            .expect("generic feature tables report filtered counts for every variant type");

        assert_eq!(
            (filtered.fp, filtered.tp, filtered.unk, filtered.ambi),
            (0, 0, 0, 0)
        );

        let row = render_row(
            0,
            "indels",
            SomaticCounts {
                truth_total: 1,
                query_total: 1,
                tp: 1,
                ..SomaticCounts::default()
            },
            &StatsRowContext {
                fp_region_size: 10,
                ci_alpha: 0.05,
                filtered: Some(filtered),
                include_filtered_columns: true,
                commandline: "som.py",
            },
        );
        assert_eq!(
            row,
            concat!(
                "0,indels,1,1,1,0,0,0,0,0.0,0.0,0.0,0.0,",
                "1.0,0.025,1.0,1.0,1.0,",
                "0.025,1.0,0.0,0.0,10,0.0,",
                "1.0,1.0,0.0,0.0,0.0,som.py-,som.py"
            )
        );

        let filtered_values = vec!["0.0".to_string(), "0.0".to_string()];
        assert_eq!(infer_type(&filtered_values), "double");
        assert_eq!(
            column_json(
                "fp.filtered",
                "fp.filtered",
                infer_type(&filtered_values),
                &filtered_values,
                false,
            ),
            "{\"values\": [0.0, 0.0], \"type\": \"double\", \"id\": \"fp.filtered\", \"label\": \"fp.filtered\"}"
        );
    }

    #[test]
    fn af_extended_uses_the_selected_truth_and_query_fields() {
        let output = tempfile::NamedTempFile::new().expect("temporary output");
        let header = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,FILTER,TRUTH_AF,QUERY_AF";
        let rows = vec![
            "0,chr1,10,TP,A,A,C,C,,0.1,0.9".to_string(),
            "0,chr1,20,FN,,A,,C,,0.8,".to_string(),
            "0,chr1,30,FP,A,,C,,,0.7,0.1".to_string(),
        ];
        write_happy_style_extended(
            output.path(),
            header,
            &rows,
            "hcc.strelka.snv",
            "0.5",
            "TRUTH_AF",
            "QUERY_AF",
        )
        .expect("extended output");
        let text = fs::read_to_string(output.path()).expect("read extended output");
        let lines = text.lines().collect::<Vec<_>>();
        assert!(lines[0].starts_with("Type,Subtype,Subset,Filter,"));
        assert!(lines[1].starts_with("SNP,*,\"[0.00,0.50)\",PASS,1,1,0,2,1,0"));
        assert!(lines[2].starts_with("SNP,*,\"[0.00,0.50)\",ALL,1,1,0,2,1,0"));
        assert!(lines[3].contains("SNP,*,\"[0.50,1.00]\",PASS,1,0,1,0,0,0,NA,0.0,NA,NA,NA"));
        assert!(lines[4].contains("SNP,*,\"[0.50,1.00]\",ALL,1,0,1,0,0,0,NA,0.0,NA,NA,NA"));
    }

    #[test]
    fn af_stats_use_truth_af_for_tp_fn_and_query_af_for_other_tags() {
        let header = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,FILTER,TRUTH_AF,QUERY_AF";
        let rows = vec![
            "0,chr1,10,TP,A,A,C,C,,0.1,0.9".to_string(),
            "0,chr1,20,FN,,A,,C,,0.8,".to_string(),
            "0,chr1,30,FP,A,,C,,LowQual,,0.1".to_string(),
            "0,chr1,40,UNK,A,,C,,,,0.7".to_string(),
            "0,chr1,50,AMBI,A,,C,,,,1.0".to_string(),
        ];
        let bins = calculate_af_stats(header, &rows, "0.5", "TRUTH_AF", "QUERY_AF", None)
            .expect("AF stats");
        assert_eq!(bins.len(), 2);
        assert_eq!(bins[0].2.truth_total, 1);
        assert_eq!(bins[0].2.query_total, 2);
        assert_eq!(bins[0].2.tp, 1);
        assert_eq!(bins[0].2.fp, 1);
        assert_eq!(bins[0].3.fp, 1);
        assert_eq!(bins[1].2.truth_total, 1);
        assert_eq!(bins[1].2.fn_count, 1);
        assert_eq!(bins[1].2.query_total, 2);
        assert_eq!(bins[1].2.unk, 1);
        assert_eq!(bins[1].2.ambi, 1);
    }

    #[test]
    fn happy_summary_matches_the_legacy_dataframe_shape_and_counts() {
        let output = tempfile::NamedTempFile::new().expect("temporary output");
        let header = ",CHROM,POS,tag,REF,REF.truth,ALT,ALT.truth,FILTER";
        let rows = vec![
            "0,chr1,10,TP,A,A,C,C,".to_string(),
            "0,chr1,20,FN,,A,,C,".to_string(),
            "0,chr1,30,FP,A,,C,,LowQual".to_string(),
            "0,chr1,40,AMBI,A,,C,,".to_string(),
            "0,chr1,50,UNK,A,,C,,".to_string(),
        ];
        write_happy_style_summary(output.path(), header, &rows, "hcc.strelka.indel")
            .expect("happy summary");
        let text = fs::read_to_string(output.path()).expect("read happy summary");
        let expected = concat!(
            ",Type,Filter,TRUTH.TOTAL,TRUTH.TP,TRUTH.FN,QUERY.TOTAL,QUERY.FP,QUERY.UNK,FP.gt,METRIC.Recall,METRIC.Precision,METRIC.Frac_NA,METRIC.F1_Score,TRUTH.TOTAL.TiTv_ratio,QUERY.TOTAL.TiTv_ratio,TRUTH.TOTAL.het_hom_ratio,QUERY.TOTAL.het_hom_ratio\n",
            "0,INDEL,PASS,2,1,1,3,0,2,NA,0.5,1.0,0.6667,0.6667,NA,NA,NA,NA\n",
            "0,INDEL,ALL,2,1,1,4,1,2,NA,0.5,0.5,0.5,0.5,NA,NA,NA,NA\n"
        );
        assert_eq!(text, expected);
    }

    #[test]
    fn metrics_json_preserves_pandas_nullable_numeric_columns() {
        let values = vec!["0".to_string(), String::new(), ".".to_string()];
        assert_eq!(infer_type(&values), "double");
        assert_eq!(
            column_json("recall2", "recall2", "double", &values, false),
            "{\"values\": [0.0, null, null], \"type\": \"double\", \"id\": \"recall2\", \"label\": \"recall2\"}"
        );

        assert_eq!(infer_type(&["1".to_string(), "2".to_string()]), "int64");
        assert_eq!(infer_type(&["1.5".to_string(), String::new()]), "double");
        assert_eq!(infer_type(&[String::new(), ".".to_string()]), "string");

        assert_eq!(
            column_json(
                "recall",
                "recall",
                "double",
                &["0.9319987680936249".to_string()],
                false,
            ),
            "{\"values\": [0.9319987680936249], \"type\": \"double\", \"id\": \"recall\", \"label\": \"recall\"}"
        );
    }

    #[test]
    fn normalization_left_aligns_trims_and_deduplicates_like_bcftools_norm() {
        let reference = BTreeMap::from([("chr1".to_string(), "AAAAAC".to_string())]);
        let deletion = raw_record("chr1\t3\t.\tAA\tA\t.\tPASS\t.");
        let duplicate = deletion.clone();
        let normalized = normalize_somatic_records(vec![deletion, duplicate], &reference);
        assert_eq!(normalized.len(), 1);
        assert_eq!(normalized[0].pos, 1);
        assert_eq!(normalized[0].ref_allele, "AA");
        assert_eq!(normalized[0].alt_allele, "A");

        let mismatch = raw_record("chr1\t1\t.\tC\tT\t.\tPASS\t.");
        assert!(normalize_somatic_records(vec![mismatch], &reference).is_empty());
    }

    #[test]
    fn chromosome_rewrite_switch_matches_the_legacy_prefix_expression() {
        let prefixed = BTreeSet::from(["chr1".to_string(), "chrM".to_string()]);
        assert_eq!(somatic_chrom("1", &prefixed, true), "chr1");
        assert_eq!(somatic_chrom("MT", &prefixed, true), "chrM");
        assert_eq!(somatic_chrom("1", &prefixed, false), "1");
        assert_eq!(somatic_chrom("GL0001", &prefixed, true), "GL0001");
    }

    #[test]
    fn caller_specific_roc_matches_the_legacy_pandas_csv_shape() {
        let output = tempfile::NamedTempFile::new().expect("temporary ROC");
        let header = ",CHROM,POS,tag,EVS,FILTER,NT";
        let rows = vec![
            "0,chr1,10,TP,10.00000000,,ref".to_string(),
            "0,chr1,20,FP,5.00000000,LowEVS,ref".to_string(),
            "0,chr1,30,FN,,,".to_string(),
        ];
        write_somatic_roc(output.path(), header, &rows, "strelka.snv").expect("write somatic ROC");
        assert_eq!(
            fs::read_to_string(output.path()).expect("read ROC"),
            concat!(
                ",EVS,tp,fp,fn,precision,recall\n",
                "0,0,1,1,1,0.50000000,0.50000000\n",
                "1,5,1,1,1,0.50000000,0.50000000\n",
                "2,10,1,0,1,1.00000000,0.50000000\n"
            )
        );
    }

    #[test]
    fn caller_specific_roc_preserves_integer_rate_columns() {
        let output = tempfile::NamedTempFile::new().expect("temporary ROC");
        let header = ",CHROM,POS,tag,QSS_NT,FILTER,NT";
        let rows = vec!["0,chr1,10,TP,10.00000000,,ref".to_string()];
        write_somatic_roc(output.path(), header, &rows, "strelka.snv.qss")
            .expect("write integer-rate somatic ROC");
        assert_eq!(
            fs::read_to_string(output.path()).expect("read ROC"),
            concat!(",QSS_NT,tp,fp,fn,precision,recall\n", "0,10,1,0,0,1,1\n")
        );
    }

    #[test]
    fn generic_feature_qual_is_serialized_as_pandas_float() {
        let truth = raw_record("chr1\t1\t.\tA\tC\t60\tPASS\t.");
        let row = render_generic_tp_row(0, &truth, &truth);
        assert!(row.ends_with(",60.00000000,60.00000000"));
        assert_eq!(cpp_default_six(0.872_429_31), "0.872429");
        assert_eq!(cpp_default_six(-1.0), "-1");
    }
}
