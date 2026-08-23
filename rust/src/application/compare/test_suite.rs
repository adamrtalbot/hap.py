//! Command-level regression tests.

#[cfg(test)]
mod scratch_tests {
    use super::super::*;
    use std::thread;

    fn comparison_record(line: &str) -> crate::domain::ComparisonRecord {
        RawVcfRecord::from_line(line, Path::new("comparison-test.vcf"))
            .unwrap()
            .into()
    }

    fn run_args(args: CompareArgs) -> Result<()> {
        super::super::run(args.validated()?)
    }

    #[test]
    fn comparison_headers_merge_inputs_and_append_legacy_annotations() {
        let truth = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##INFO=<ID=TRUTH_ONLY,Number=1,Type=String,Description=\"truth\">".to_string(),
            "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"wrong\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tT".to_string(),
        ];
        let query = vec![
            "##INFO=<ID=QUERY_ONLY,Number=1,Type=String,Description=\"query\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQ".to_string(),
        ];
        let headers = build_vcf_headers(&truth, &query, true, false, false, "QUAL");

        assert!(headers.iter().any(|line| line.contains("ID=TRUTH_ONLY,")));
        assert!(headers.iter().any(|line| line.contains("ID=QUERY_ONLY,")));
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=Q_FILTERED,"))
        );
        assert!(headers.iter().any(|line| line == "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Additional comparison information\">"));
        assert!(headers.iter().any(|line| line == "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Variant quality for ROC creation.\">"));
        assert!(
            !headers
                .iter()
                .any(|line| line.contains("Description=\"wrong\""))
        );
        assert_eq!(
            headers.last().unwrap(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY"
        );
    }

    #[test]
    fn comparison_headers_only_declare_filtered_calls_when_enabled() {
        let headers = build_vcf_headers(&[], &[], false, false, false, "QUAL");
        assert!(!headers.iter().any(|line| line.contains("ID=Q_FILTERED,")));
    }

    fn fixture_path(file: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/synth-snp-match")
            .join(file)
    }

    fn test_root(label: &str) -> PathBuf {
        let id = SCRATCH_RUN_ID.fetch_add(1, Ordering::Relaxed);
        let root =
            std::env::temp_dir().join(format!("hap-compare-{label}-{}-{id}", std::process::id()));
        fs::create_dir_all(&root).expect("create test root");
        root
    }

    fn args(report_prefix: PathBuf, scratch_parent: &Path, keep_scratch: bool) -> CompareArgs {
        let mut args = CompareArgs::with_paths(
            fixture_path("truth.vcf").display().to_string(),
            fixture_path("query.vcf").display().to_string(),
            fixture_path("ref.fa").display().to_string(),
            report_prefix.display().to_string(),
        );
        args.scratch_prefix = Some(scratch_parent.display().to_string());
        args.keep_scratch = keep_scratch;
        args
    }

    #[test]
    fn preprocessing_defaults_are_asymmetric_between_truth_and_query() {
        let root = test_root("preprocessing-defaults");
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.bcftools_norm = true;
        let truth = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth.vcf.gz"),
            true,
            options.preprocess_truth,
        );
        let query = build_preprocess_args(
            &options,
            &options.query,
            &root.join("query.vcf.gz"),
            false,
            true,
        );

        assert!(!truth.leftshift);
        assert!(!truth.decompose);
        assert!(!truth.bcftools_norm);
        assert!(query.leftshift);
        assert!(query.decompose);
        assert!(query.bcftools_norm);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn custom_roc_regions_retain_the_implicit_aggregate_region() {
        let mut regions = vec!["CONF".to_string()];
        ensure_aggregate_roc_region(&mut regions);
        ensure_aggregate_roc_region(&mut regions);
        assert_eq!(regions, ["*", "CONF"]);
    }

    #[test]
    fn scmp_engines_apply_legacy_preprocessing_defaults() {
        let root = test_root("scmp-policy");
        let mut somatic = args(root.join("somatic"), &root, false);
        somatic.engine = CompareEngine::ScmpSomatic;
        somatic.preprocess_truth = true;
        somatic.bcftools_norm = true;
        normalize_engine_preprocessing(&mut somatic);
        assert!(somatic.somatic);
        assert_eq!(somatic.set_gt, None);
        assert!(!somatic.preprocess_truth);
        assert!(somatic.no_leftshift);
        assert!(!somatic.bcftools_norm);
        assert!(!effective_decomposition(&somatic));

        let mut distance = args(root.join("distance"), &root, false);
        distance.engine = CompareEngine::ScmpDistance;
        normalize_engine_preprocessing(&mut distance);
        assert_eq!(
            distance.set_gt,
            Some(crate::application::SomaticGtMode::First)
        );
        assert!(!effective_decomposition(&distance));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn comparison_log_honors_verbose_and_quiet_levels() {
        let root = test_root("comparison-log");
        let logfile = root.join("comparison.log");
        let mut options = args(root.join("result"), &root, false);
        options.logfile = Some(logfile.display().to_string());
        options.verbose = true;
        initialize_compare_log(&options).unwrap();
        log_compare_info(&options, "comparison stage").unwrap();
        assert_eq!(
            fs::read_to_string(&logfile).unwrap(),
            "INFO comparison stage\n"
        );

        options.quiet = true;
        log_compare_info(&options, "suppressed").unwrap();
        assert!(!fs::read_to_string(&logfile).unwrap().contains("suppressed"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn preprocess_truth_and_negative_switches_control_both_sides() {
        let root = test_root("preprocessing-overrides");
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.preprocess_truth = true;

        let truth = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth-enabled.vcf.gz"),
            true,
            options.preprocess_truth,
        );
        assert!(truth.leftshift);
        assert!(truth.decompose);

        options.no_leftshift = true;
        options.no_decompose = true;
        let truth_disabled = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth-disabled.vcf.gz"),
            true,
            options.preprocess_truth,
        );
        let query_disabled = build_preprocess_args(
            &options,
            &options.query,
            &root.join("query-disabled.vcf.gz"),
            false,
            true,
        );
        assert!(!truth_disabled.leftshift);
        assert!(!truth_disabled.decompose);
        assert!(!query_disabled.leftshift);
        assert!(!query_disabled.decompose);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn remaining_germline_preprocess_controls_propagate_per_side() {
        let root = test_root("preprocessing-controls");
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.usefiltered_truth = true;
        options.filters_only = Some("LowQual,q10".to_string());
        options.convert_gvcf_truth = true;
        options.convert_gvcf_query = false;
        options.filter_nonref = true;
        options.preprocess_truth = true;
        options.bcftools_norm = true;
        options.fixchr = Some(true);
        options.gender = crate::application::PreprocessGender::Male;
        options.preprocess_window = 4096;

        let truth = build_preprocess_args(
            &options,
            &options.truth,
            &root.join("truth.bcf"),
            !options.usefiltered_truth,
            options.preprocess_truth,
        );
        let query = build_preprocess_args(
            &options,
            &options.query,
            &root.join("query.bcf"),
            options.pass_only,
            true,
        );
        assert!(!truth.pass_only);
        assert_eq!(truth.filters_only, None);
        assert!(truth.convert_gvcf_to_vcf);
        assert!(truth.filter_nonref);
        assert_eq!(query.filters_only.as_deref(), Some("LowQual,q10"));
        assert!(!query.convert_gvcf_to_vcf);
        assert!(query.filter_nonref);
        for side in [&truth, &query] {
            assert!(side.bcftools_norm);
            assert_eq!(side.fixchr, Some(true));
            assert_eq!(side.gender, crate::application::PreprocessGender::Male);
            assert_eq!(side.window_size, 4096);
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn no_roc_no_counts_no_json_emits_the_legacy_artifact_shape() {
        let root = test_root("minimal-artifacts");
        let prefix = root.join("result");
        let mut options = args(prefix.clone(), &root.join("scratch"), false);
        options.no_roc = true;
        options.no_write_counts = true;
        options.no_json = true;
        run_args(options).unwrap();

        for suffix in [
            "summary.csv",
            "vcf.gz",
            "vcf.gz.tbi",
            "roc.all.csv.gz",
            "runinfo.json",
        ] {
            assert!(suffixed_report_path(&prefix, suffix).is_file(), "{suffix}");
        }
        for suffix in [
            "extended.csv",
            "metrics.json.gz",
            "roc.Locations.SNP.csv.gz",
            "roc.Locations.SNP.PASS.csv.gz",
            "roc.Locations.INDEL.csv.gz",
            "roc.Locations.INDEL.PASS.csv.gz",
        ] {
            assert!(!suffixed_report_path(&prefix, suffix).exists(), "{suffix}");
        }
        let roc = vcf::read_text(&suffixed_report_path(&prefix, "roc.all.csv.gz")).unwrap();
        for line in roc.lines().skip(1) {
            let cells = line.split(',').collect::<Vec<_>>();
            assert_eq!(cells[6], "*");
            if cells[0] == "INDEL" {
                for block_start in [16usize, 23, 30, 37, 44, 51, 58] {
                    assert_eq!(cells[block_start + 1], ".");
                    assert_eq!(cells[block_start + 2], ".");
                }
            }
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bcf_mode_keeps_bcf_intermediates_and_publishes_only_the_bcf_report_pair() {
        let root = test_root("bcf-artifacts");
        let prefix = root.join("result");
        let scratch_parent = root.join("scratch");
        let mut options = args(prefix.clone(), &scratch_parent, true);
        options.bcf = true;
        run_args(options).unwrap();

        let bcf_report = suffixed_report_path(&prefix, "bcf");
        assert!(bcf_report.is_file());
        assert!(suffixed_report_path(&prefix, "bcf.csi").is_file());
        assert!(!suffixed_report_path(&prefix, "vcf.gz").exists());
        assert!(!suffixed_report_path(&prefix, "vcf.gz.tbi").exists());
        let (_, records) = vcf::load_raw_vcf(&bcf_report).unwrap();
        assert!(!records.is_empty());

        let scratch_runs = fs::read_dir(&scratch_parent)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect::<Vec<_>>();
        assert_eq!(scratch_runs.len(), 1);
        let scratch = &scratch_runs[0];
        for name in [
            "truth.prep.bcf",
            "truth.prep.bcf.csi",
            "query.prep.bcf",
            "query.prep.bcf.csi",
        ] {
            assert!(scratch.join(name).is_file(), "{name}");
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn paired_bcf_inputs_implicitly_enable_bcf_reports_and_intermediates() {
        let root = test_root("implicit-bcf-artifacts");
        let truth_bcf = root.join("truth.bcf");
        let query_bcf = root.join("query.bcf");
        for (source, destination) in [
            (fixture_path("truth.vcf"), &truth_bcf),
            (fixture_path("query.vcf"), &query_bcf),
        ] {
            let (headers, records) = vcf::load_raw_vcf(&source).unwrap();
            vcf::write_raw_vcf(destination, &headers, &records).unwrap();
        }

        let prefix = root.join("result");
        let scratch_parent = root.join("scratch");
        let mut options = args(prefix.clone(), &scratch_parent, true);
        options.truth = truth_bcf.display().to_string();
        options.query = query_bcf.display().to_string();
        assert!(!options.bcf, "the CLI flag is intentionally absent");
        run_args(options).unwrap();

        assert!(suffixed_report_path(&prefix, "bcf").is_file());
        assert!(suffixed_report_path(&prefix, "bcf.csi").is_file());
        assert!(!suffixed_report_path(&prefix, "vcf.gz").exists());
        let runinfo = fs::read_to_string(suffixed_report_path(&prefix, "runinfo.json")).unwrap();
        assert!(
            runinfo.contains("\"bcf\":false"),
            "implicit output selection must not rewrite the explicit CLI flag: {runinfo}"
        );
        let scratch_runs = child_directories(&scratch_parent);
        assert_eq!(scratch_runs.len(), 1);
        assert!(scratch_runs[0].join("truth.prep.bcf").is_file());
        assert!(scratch_runs[0].join("query.prep.bcf").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn missing_report_parent_is_rejected_without_creating_artifacts() {
        let root = test_root("missing-report-parent");
        let missing_parent = root.join("missing");
        let scratch_parent = root.join("scratch");
        let options = args(missing_parent.join("result"), &scratch_parent, false);

        let error = run_args(options).expect_err("missing report parents must be rejected");
        assert!(error.to_string().contains("output path does not exist"));
        assert!(!missing_parent.exists());
        assert!(!scratch_parent.exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn regions_reject_legacy_overlaps_and_ordering_but_targets_do_not() {
        let root = test_root("region-overlap-check");
        for (label, contents) in [
            ("overlap", "chr1\t0\t10\nchr1\t5\t12\n"),
            ("out-of-order", "chr1\t10\t12\nchr1\t0\t6\n"),
        ] {
            let bed = root.join(format!("{label}.bed"));
            fs::write(&bed, contents).unwrap();
            let mut options = args(
                root.join(format!("{label}-result")),
                &root.join(format!("{label}-scratch")),
                false,
            );
            options.regions_bedfile = Some(bed.display().to_string());
            let error = run_args(options).expect_err("invalid -R BED must fail before comparison");
            assert!(
                error
                    .to_string()
                    .contains("The regions bed file (specified using -R) has overlaps")
            );
        }

        let targets = root.join("targets.bed");
        fs::write(&targets, "chr1\t10\t12\nchr1\t0\t12\n").unwrap();
        let target_prefix = root.join("target-result");
        let mut options = args(target_prefix.clone(), &root.join("target-scratch"), false);
        options.targets_bedfile = Some(targets.display().to_string());
        run_args(options).expect("-T keeps accepting overlapping or out-of-order intervals");
        assert!(suffixed_report_path(&target_prefix, "summary.csv").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn empty_truth_with_default_locations_fails_before_query_comparison() {
        let root = test_root("empty-truth-default-locations");
        let truth = root.join("truth.vcf");
        fs::write(
            &truth,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=16>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n",
            ),
        )
        .unwrap();
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.truth = truth.display().to_string();

        let error = run_args(options).expect_err("legacy derives default contigs from truth calls");
        assert!(
            error
                .to_string()
                .contains("Truth and reference have no chromosomes in common")
        );

        let explicit_prefix = root.join("explicit-result");
        let mut explicit = args(
            explicit_prefix.clone(),
            &root.join("explicit-scratch"),
            false,
        );
        explicit.truth = truth.display().to_string();
        explicit.locations = Some("chr1".to_string());
        run_args(explicit).expect("an explicit contig bypasses legacy default-contig discovery");
        assert!(suffixed_report_path(&explicit_prefix, "summary.csv").is_file());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn default_locations_exclude_query_only_contigs() {
        let root = test_root("query-only-contig");
        let reference = root.join("ref.fa");
        fs::write(
            &reference,
            ">chr1\nAACCGGTTAACCGGTT\n>chrX\nAACCGGTTAACCGGTT\n",
        )
        .unwrap();
        fs::write(
            reference.with_extension("fa.fai"),
            "chr1\t16\t6\t16\t17\nchrX\t16\t29\t16\t17\n",
        )
        .unwrap();
        let truth = root.join("truth.vcf");
        fs::write(
            &truth,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=16>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\n",
                "chr1\t5\t.\tG\tT\t60\tPASS\t.\tGT\t1/1\n",
            ),
        )
        .unwrap();
        let query = root.join("query.vcf");
        fs::write(
            &query,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=16>\n",
                "##contig=<ID=chrX,length=16>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tQUERY\n",
                "chr1\t5\t.\tG\tT\t60\tPASS\t.\tGT\t1/1\n",
                "chrX\t5\t.\tG\tA\t60\tPASS\t.\tGT\t0/1\n",
            ),
        )
        .unwrap();
        let prefix = root.join("result");
        let mut options = CompareArgs::with_paths(
            truth.display().to_string(),
            query.display().to_string(),
            reference.display().to_string(),
            prefix.display().to_string(),
        );
        options.scratch_prefix = Some(root.join("scratch").display().to_string());
        run_args(options).unwrap();

        let (_, records) = vcf::load_raw_vcf(&suffixed_report_path(&prefix, "vcf.gz")).unwrap();
        assert!(records.iter().all(|record| record.chrom == "chr1"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bcf_report_size_spans_the_full_reference() {
        let contig_lengths = BTreeMap::from([("chr1".to_string(), 100), ("chrX".to_string(), 40)]);
        let contigs_in_play = BTreeSet::from(["chr1".to_string()]);

        assert_eq!(
            report_subset_size(&contig_lengths, &contigs_in_play, false, false),
            100,
            "ordinary reports retain their active-contig size"
        );
        assert_eq!(
            report_subset_size(&contig_lengths, &contigs_in_play, true, false),
            140,
            "legacy BCF reports size the complete reference dictionary"
        );

        let aliased_lengths = BTreeMap::from([
            ("1".to_string(), 100),
            ("chr1".to_string(), 100),
            ("chrX".to_string(), 40),
        ]);
        assert_eq!(
            report_subset_size(&aliased_lengths, &contigs_in_play, false, true),
            200,
            "implicit BCF reports retain both declared chromosome aliases"
        );
    }

    #[test]
    fn scmp_bcf_mode_materializes_confidence_padding_without_exposing_vcf() {
        let root = test_root("scmp-bcf-artifacts");
        let prefix = root.join("result");
        let confidence = root.join("confident.bed");
        fs::write(&confidence, "chr1\t0\t16\n").unwrap();
        let mut options = args(prefix.clone(), &root.join("scratch"), false);
        options.engine = CompareEngine::ScmpDistance;
        options.bcf = true;
        options.fp_bedfile = Some(confidence.display().to_string());
        run_args(options).unwrap();

        assert!(suffixed_report_path(&prefix, "bcf").is_file());
        assert!(suffixed_report_path(&prefix, "bcf.csi").is_file());
        assert!(!suffixed_report_path(&prefix, "vcf.gz").exists());
        assert!(!suffixed_report_path(&prefix, "vcf.gz.tbi").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn scmp_output_vtc_is_forwarded_to_ga4gh_quantification() {
        let root = test_root("scmp-output-vtc");
        let prefix = root.join("result");
        let mut options = args(prefix.clone(), &root.join("scratch"), false);
        options.engine = CompareEngine::ScmpDistance;
        options.output_vtc = true;
        run_args(options).unwrap();

        let (headers, records) =
            vcf::load_raw_vcf(&suffixed_report_path(&prefix, "vcf.gz")).unwrap();
        assert!(headers.iter().any(|line| line.contains("##INFO=<ID=VTC,")));
        assert!(records.iter().any(|record| record.info.contains("VTC=")));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn subset_derivation_stops_regions_at_the_next_info_field() {
        let rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            record: comparison_record(concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\t",
                "BS=7;Regions=CONF,TS_boundary,TS_contained;AF=0.5;VTC=nuc__s\t",
                "GT:BD:BK:BVT:BLT:BI\t",
                "0/1:TP:gm:SNP:het:ti\t0/1:TP:gm:SNP:het:ti"
            )),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];

        let counts = derive_subset_counts(&rows, false);
        assert_eq!(
            counts.keys().cloned().collect::<Vec<_>>(),
            vec!["TS_boundary", "TS_contained"]
        );
        for subset in ["TS_boundary", "TS_contained"] {
            let snp = &counts[subset]["SNP"];
            assert_eq!(snp.truth_total.total, 1);
            assert_eq!(snp.query_total.total, 1);
        }
    }

    #[test]
    fn requantify_handoff_drops_only_provisional_truth_set_membership() {
        let rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            record: comparison_record(concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\t",
                "BS=7;Regions=CONF,TS_boundary,EXTRA,TS_contained;RegionsExtent=7-7\t",
                "GT:BD\t0/1:TP\t0/1:TP"
            )),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];

        let sanitized = sanitize_requantify_handoff_rows(&rows);

        assert!(
            rows[0]
                .record
                .info
                .contains("Regions=CONF,TS_boundary,EXTRA,TS_contained")
        );
        assert!(
            sanitized[0]
                .record
                .info
                .contains("BS=7;Regions=CONF,EXTRA;RegionsExtent=7-7")
        );
        assert!(!sanitized[0].record.info.contains("TS_boundary"));
        assert!(!sanitized[0].record.info.contains("TS_contained"));
    }

    #[test]
    fn vcfeval_ignores_external_runtime_flags() {
        let root = test_root("vcfeval-ignored-flags");
        let mut options = CompareArgs::with_paths(
            fixture_path("truth.vcf").display().to_string(),
            fixture_path("query.vcf").display().to_string(),
            fixture_path("ref.fa").display().to_string(),
            root.join("result").display().to_string(),
        );
        options.scratch_prefix = Some(root.join("scratch").display().to_string());
        options.engine = CompareEngine::Vcfeval;
        options.engine_vcfeval = Some("definitely-absent-rtg-for-test".to_string());
        options.engine_vcfeval_template = Some(root.join("absent.sdf").display().to_string());
        run_args(options).unwrap();
        assert!(root.join("result.summary.csv").is_file());
        let (_, records) = vcf::load_raw_vcf(&root.join("result.vcf.gz")).unwrap();
        assert_eq!(records.len(), 1);
        assert!(
            records[0]
                .samples
                .iter()
                .all(|sample| sample.contains(":TP:gm"))
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn native_vcfeval_with_paired_bcf_inputs_publishes_bcf_and_preserves_info() {
        let root = test_root("vcfeval-contract");
        let truth_bcf = root.join("truth.bcf");
        let query_bcf = root.join("query.bcf");
        for (source, destination) in [
            (fixture_path("truth.vcf"), &truth_bcf),
            (fixture_path("query.vcf"), &query_bcf),
        ] {
            let (headers, records) = vcf::load_raw_vcf(&source).unwrap();
            vcf::write_raw_vcf(destination, &headers, &records).unwrap();
        }

        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.truth = truth_bcf.display().to_string();
        options.query = query_bcf.display().to_string();
        options.engine = CompareEngine::Vcfeval;
        options.output_vtc = true;
        run_args(options).unwrap();

        assert!(!root.join("result.vcf.gz").exists());
        let (headers, records) = vcf::load_raw_vcf(&root.join("result.bcf")).unwrap();
        assert!(headers.iter().any(|line| line.contains("ID=VTC,")));
        assert!(headers.iter().any(|line| line.contains("ID=XCMP,")));
        assert!(
            records[0].info.contains("VTC=nuc__s,al__s,homalt__s"),
            "{}",
            records[0].info
        );
        assert!(records[0].info.contains("XCMP="), "{}", records[0].info);

        let preserve_prefix = root.join("preserve-result");
        let preserve_scratch = root.join("preserve-scratch");
        let mut preserve = args(preserve_prefix.clone(), &preserve_scratch, false);
        preserve.engine = CompareEngine::Vcfeval;
        preserve.preserve_info = true;
        run_args(preserve).unwrap();
        assert!(suffixed_report_path(&preserve_prefix, "runinfo.json").is_file());
        for suffix in ["summary.csv", "extended.csv", "vcf.gz", "metrics.json.gz"] {
            assert!(suffixed_report_path(&preserve_prefix, suffix).is_file());
        }
        let (_, preserved_records) =
            vcf::load_raw_vcf(&suffixed_report_path(&preserve_prefix, "vcf.gz")).unwrap();
        assert_eq!(preserved_records.len(), 1);
        assert!(child_directories(&preserve_scratch).is_empty());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn scmp_engines_use_non_haplotype_comparison_semantics() {
        let root = test_root("scmp-engines");
        let fixture = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/synth-homopolymer-insertion");
        let engine_args = |prefix: PathBuf| {
            let mut options = CompareArgs::with_paths(
                fixture.join("truth.vcf").display().to_string(),
                fixture.join("query.vcf").display().to_string(),
                fixture.join("ref.fa").display().to_string(),
                prefix.display().to_string(),
            );
            options.scratch_prefix = Some(root.join("scratch").display().to_string());
            options
        };
        let mut somatic = engine_args(root.join("somatic"));
        somatic.engine = CompareEngine::ScmpSomatic;
        run_args(somatic).unwrap();
        let somatic_summary = fs::read_to_string(root.join("somatic.summary.csv")).unwrap();

        let mut distance = engine_args(root.join("distance"));
        distance.engine = CompareEngine::ScmpDistance;
        distance.engine_scmp_distance = 30;
        run_args(distance).unwrap();
        let distance_summary = fs::read_to_string(root.join("distance.summary.csv")).unwrap();
        // Legacy AlleleMatcher's RefVar constructor uses ALT length for the
        // reference end. Consequently these two ordinary VCF-equivalent
        // homopolymer insertions do not hash alike in scmp-somatic, while
        // distance mode still pairs their overlapping intervals.
        assert_ne!(somatic_summary, distance_summary);

        run_args(engine_args(root.join("xcmp"))).unwrap();
        let xcmp_summary = fs::read_to_string(root.join("xcmp.summary.csv")).unwrap();
        assert_ne!(somatic_summary, xcmp_summary);
        assert!(distance_summary.contains("INDEL,ALL,1,1,0,1,0,0"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn preserve_info_and_vtc_decorate_rows_and_headers() {
        let source = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tC\t30\tPASS\tSCORE=9\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        let mut rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            record: comparison_record(concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\tBS=7;Regions=CONF\t",
                "GT:BD:BK:BVT:BLT:QQ\t",
                "0/1:TP:gm:SNP:het:30\t0/1:TP:gm:SNP:het:30"
            )),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];
        decorate_output_rows(&mut rows, &[source], &[], true, true, "QUAL").unwrap();
        assert!(
            rows[0].record.info.contains(concat!(
                "BS=7;IQQ=30;SCORE=9;ctype=simple:match;gtt1=gt_het;",
                "gtt2=gt_het;kind=match;type=TP;Regions=CONF;RegionsExtent=7-7;",
                "XCMP=TP:match:gt_het:gt_het:simple:match;",
                "VTC=nuc__s,al__s,het__rs"
            )),
            "{}",
            rows[0].record.info
        );
        let headers = build_vcf_headers(&[], &[], false, true, true, "QUAL");
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=VTC,"))
        );
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##INFO=<ID=XCMP,"))
        );
        let regions = headers
            .iter()
            .position(|line| line.starts_with("##INFO=<ID=Regions,"))
            .unwrap();
        let extent = headers
            .iter()
            .position(|line| line.starts_with("##INFO=<ID=RegionsExtent,"))
            .unwrap();
        let vtc = headers
            .iter()
            .position(|line| line.starts_with("##INFO=<ID=VTC,"))
            .unwrap();
        assert_eq!(extent, regions + 1);
        assert!(vtc > extent);
    }

    #[test]
    fn prefixed_custom_roc_field_preserves_xcmp_literal_lookup_miss() {
        let source = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tC\t30\tPASS\tSCORE=9\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        let mut rows = vec![AnnotatedRow {
            sort_key: ("chr1".to_string(), 7, 0, 0),
            record: comparison_record(concat!(
                "chr1\t7\t.\tA\tC\t30\tPASS\tBS=7;Regions=CONF\t",
                "GT:BD:BK:BVT:BLT:QQ\t",
                "0/1:TP:gm:SNP:het:30\t0/1:TP:gm:SNP:het:30"
            )),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }];

        decorate_output_rows(&mut rows, &[source], &[], true, false, "INFO.SCORE").unwrap();

        let record = &rows[0].record;
        assert!(record.info.contains("SCORE=9"));
        assert!(!record.info.contains("INFO.SCORE="));
        assert!(!record.info.contains("IQQ="));
        assert_eq!(
            record.sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("nan")
        );

        let headers = build_vcf_headers(&[], &[], false, false, true, "INFO.SCORE");
        assert!(headers.iter().any(|line| {
            line == "##INFO=<ID=IQQ,Number=1,Type=Float,Description=\"Quality value for query variants (INFO.SCORE).\">"
        }));
    }

    #[test]
    fn metadata_helpers_preserve_legacy_multiallelic_contracts() {
        let record = RawVcfRecord::from_line(
            concat!(
                "chr21\t19323424\t.\tCGTGT\tC,CGTGTGTGTGT,CGT\t0\t.\t.\t",
                "GT:BVT:BLT\t1/2:INDEL:hetalt\t./.:NOCALL:nocall"
            ),
            Path::new("metadata.vcf"),
        )
        .unwrap();
        assert_eq!(legacy_regions_extent(&record), "19323424-19323428");

        let mixed = RawVcfRecord::from_line(
            concat!(
                "chr21\t10\t.\tAGT\tA,AGTGT\t0\t.\t.\t",
                "GT:BVT:BLT\t1/2:INDEL:hetalt\t./.:NOCALL:nocall"
            ),
            Path::new("metadata.vcf"),
        )
        .unwrap();
        assert_eq!(
            legacy_vtc(&mixed, &mixed.sample_map(0), &mixed.sample_map(1)),
            "nuc__i,nuc__d,al__i,al__d,hetalt__id,nocall__nc"
        );
    }

    #[test]
    fn semantic_preserve_key_reorders_alts_without_colliding_deletions() {
        let first = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tAT,ATT\t0\t.\t.\tGT\t1/2",
            Path::new("source.vcf"),
        )
        .unwrap();
        let reordered = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tATT,AT\t0\t.\t.\tGT\t1/2",
            Path::new("source.vcf"),
        )
        .unwrap();
        let other_ref = RawVcfRecord::from_line(
            "chr1\t7\t.\tAA\tAT,ATT\t0\t.\t.\tGT\t1/2",
            Path::new("source.vcf"),
        )
        .unwrap();
        assert_eq!(semantic_info_key(&first), semantic_info_key(&reordered));
        assert_ne!(semantic_info_key(&first), semantic_info_key(&other_ref));

        let symbolic_n = RawVcfRecord::from_line(
            "chr1\t7\t.\tN\t<DEL>\t0\t.\t.\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        let symbolic_a = RawVcfRecord::from_line(
            "chr1\t7\t.\tA\t<DEL>\t0\t.\t.\tGT\t0/1",
            Path::new("source.vcf"),
        )
        .unwrap();
        assert_eq!(
            semantic_info_key(&symbolic_n),
            semantic_info_key(&symbolic_a)
        );
    }

    fn child_directories(parent: &Path) -> Vec<PathBuf> {
        fs::read_dir(parent)
            .expect("read scratch parent")
            .filter_map(|entry| entry.ok())
            .map(|entry| entry.path())
            .filter(|path| path.is_dir())
            .collect()
    }

    #[test]
    fn concurrent_runs_share_parent_without_colliding_and_cleanup() {
        let root = test_root("concurrent");
        let scratch_parent = root.join("scratch");
        fs::create_dir_all(root.join("first")).unwrap();
        fs::create_dir_all(root.join("second")).unwrap();
        let first = args(root.join("first/result"), &scratch_parent, false);
        let second = args(root.join("second/result"), &scratch_parent, false);

        let first_run = thread::spawn(move || run_args(first));
        let second_run = thread::spawn(move || run_args(second));
        first_run.join().expect("first thread panicked").unwrap();
        second_run.join().expect("second thread panicked").unwrap();

        assert_eq!(
            fs::read(root.join("first/result.summary.csv")).unwrap(),
            fs::read(root.join("second/result.summary.csv")).unwrap()
        );
        assert!(
            child_directories(&scratch_parent).is_empty(),
            "successful invocations must delete only their own run directories"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn keep_scratch_retains_run_but_errors_cleanup_by_default() {
        let root = test_root("lifecycle");
        let kept_parent = root.join("kept");
        fs::create_dir(root.join("kept-output")).unwrap();
        run_args(args(root.join("kept-output/result"), &kept_parent, true)).unwrap();

        let kept = child_directories(&kept_parent);
        assert_eq!(kept.len(), 1, "--keep-scratch retains the unique run");
        assert!(kept[0].join("truth.prep.vcf.gz").is_file());
        assert!(kept[0].join("truth.prep.vcf.gz.tbi").is_file());
        assert!(kept[0].join("query.prep.vcf.gz").is_file());
        assert!(kept[0].join("query.prep.vcf.gz.tbi").is_file());

        let error_parent = root.join("error");
        fs::create_dir(root.join("error-output")).unwrap();
        let mut failing = args(root.join("error-output/result"), &error_parent, false);
        failing.truth = root.join("missing.vcf").display().to_string();
        assert!(run_args(failing).is_err());
        assert!(
            child_directories(&error_parent).is_empty(),
            "error paths must delete their invocation directory"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn report_suffixes_preserve_dotted_prefixes() {
        let prefix = Path::new("reports/sample.v1");
        assert_eq!(
            suffixed_report_path(prefix, "summary.csv"),
            Path::new("reports/sample.v1.summary.csv")
        );
    }

    #[test]
    fn stratification_tsv_requantifies_reports_and_vcf() {
        let root = test_root("stratification");
        let bed = root.join("focus.bed");
        let tsv = root.join("regions.tsv");
        fs::write(&bed, "chr1\t4\t5\n").unwrap();
        fs::write(&tsv, "FOCUS\tfocus.bed\n").unwrap();
        let mut options = args(root.join("result"), &root.join("scratch"), false);
        options.strat_tsv = Some(tsv.display().to_string());
        run_args(options).unwrap();

        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        assert!(extended.lines().any(|line| line.contains(",FOCUS,")));
        let (_, records) = vcf::load_raw_vcf(&root.join("result.vcf.gz")).unwrap();
        assert!(
            records
                .iter()
                .any(|record| record.info.contains("Regions=FOCUS"))
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn explicit_scratch_cleanup_propagates_removal_errors() {
        let root = test_root("cleanup-error");
        let scratch = ScratchRun::create(&root, false).unwrap();
        let scratch_path = scratch.path().to_path_buf();
        fs::remove_dir_all(&scratch_path).unwrap();
        fs::write(&scratch_path, "not a directory").unwrap();

        let error = scratch.cleanup().unwrap_err();
        assert!(error.to_string().contains("failed to remove scratch run"));

        fs::remove_file(scratch_path).unwrap();
        fs::remove_dir_all(root).unwrap();
    }
}

#[cfg(test)]
mod memory_guards {
    use super::super::*;

    fn comparison_record(line: &str) -> crate::domain::ComparisonRecord {
        RawVcfRecord::from_line(line, Path::new("comparison-test.vcf"))
            .unwrap()
            .into()
    }

    // Class 1 pinning helper for ergonomic qual overrides in tests.
    impl Variant {
        fn with_qual(mut self, q: &str) -> Self {
            self.qual = q.to_string();
            self
        }
    }

    fn het(gt: &str) -> Variant {
        Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 10,
                ref_allele: "A".to_string(),
                alt_allele: "G".to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: gt.to_string(),
        }
    }

    #[test]
    fn streaming_clusters_follow_sequence_dictionary_across_9_to_10() -> Result<()> {
        let at = |chrom: &str, pos: usize| {
            let mut variant = het("0/1");
            variant.key.chrom = chrom.to_string();
            variant.key.pos = pos;
            variant
        };
        let truth = vec![at("9", 100), at("10", 1)];
        let query = vec![at("9", 200), at("10", 1)];
        let headers = [
            "##contig=<ID=9,length=1000>".to_string(),
            "##contig=<ID=10,length=1000>".to_string(),
        ];
        let ranks = comparison_contig_ranks(&headers, &[]);

        let clusters = StreamingClusters::new(
            truth.into_iter().map(Ok),
            query.into_iter().map(Ok),
            200,
            ranks,
        )
        .collect::<Result<Vec<_>>>()?;

        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].chrom, "9");
        assert_eq!((clusters[0].truth.len(), clusters[0].query.len()), (1, 1));
        assert_eq!(clusters[1].chrom, "10");
        assert_eq!((clusters[1].truth.len(), clusters[1].query.len()), (1, 1));
        Ok(())
    }

    #[test]
    fn estimated_state_count_doubles_per_unphased_het() {
        let variants = vec![het("0/1"), het("0/1"), het("0/1")];
        assert_eq!(estimated_state_count(&variants), 8);
    }

    #[test]
    fn estimated_state_count_phased_het_does_not_double() {
        let variants = vec![het("0|1"), het("0|1"), het("0|1")];
        assert_eq!(estimated_state_count(&variants), 1);
    }

    #[test]
    fn estimated_state_count_homozygous_does_not_double() {
        let variants = vec![het("1/1"), het("1/1"), het("1/1")];
        assert_eq!(estimated_state_count(&variants), 1);
    }

    #[test]
    fn estimated_state_count_saturates_above_threshold() {
        // 40 heterozygous unphased variants would be 2^40 states without the
        // saturating short-circuit; we expect `usize::MAX` as the sentinel so
        // the caller knows enumeration is not affordable.
        let variants: Vec<Variant> = (0..40).map(|_| het("0/1")).collect();
        assert_eq!(estimated_state_count(&variants), usize::MAX);
    }

    #[test]
    fn custom_enumeration_threshold_above_default_is_honored() {
        let variants = (0..15).map(|_| het("0/1")).collect::<Vec<_>>();
        assert_eq!(estimated_state_count(&variants), usize::MAX);
        assert_eq!(estimated_state_count_with_limit(&variants, 32_768), 32_768);
    }

    #[test]
    fn build_clusters_does_not_split_a_connected_cluster_at_variant_cap() {
        // MAX_CLUSTER_VARIANTS + 2 variants packed within 1 bp of each other
        // remains one connected cluster. Production reports an explicit
        // resource error rather than silently changing comparison semantics.
        let mut truth = Vec::new();
        for offset in 0..(MAX_CLUSTER_VARIANTS + 2) {
            truth.push(Variant {
                key: VariantKey {
                    chrom: "chr1".to_string(),
                    pos: 100 + offset,
                    ref_allele: "A".to_string(),
                    alt_allele: "G".to_string(),
                },
                qual: "30".to_string(),
                filter: "PASS".to_string(),
                gt: "0/1".to_string(),
            });
        }
        let clusters = build_clusters(&truth, &[]);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].truth.len(), MAX_CLUSTER_VARIANTS + 2);
    }

    /// Class D pin: chr21:38861935 INDEL hetalt-vs-homalt-of-shared-allele.
    /// Truth `T→TAA,TA 1|1` selects {TAA}, query `T→TAA,TA 1/2` selects
    /// {TAA, TA}. Selected sets differ but overlap on TAA. INDEL → BK=lm.
    #[test]
    fn class_d_indel_overlap_emits_lm() {
        let truth = variant(38861935, "T", "TAA,TA", "1|1");
        let query = variant(38861935, "T", "TAA,TA", "1/2");
        assert_eq!(compute_paired_bk(&truth, &query), "lm");
    }

    /// Class D pin: chr21:9922359 SNP hetalt-vs-het-overlapping. Truth
    /// `T→A,C 1|0` selects {A}, query `T→A,C 1/2` selects {A, C}. SNP →
    /// BK=`.` (legacy quirk).
    #[test]
    fn class_d_snp_overlap_emits_dot() {
        let truth = variant(9922359, "T", "A,C", "1|0");
        let query = variant(9922359, "T", "A,C", "1/2");
        assert_eq!(compute_paired_bk(&truth, &query), ".");
    }

    /// Class D pin: same selected set, different multiset (truth het,
    /// query homalt of same allele) → BK=`am`.
    #[test]
    fn class_d_same_set_diff_multiset_emits_am() {
        let truth = variant(100, "T", "A", "0|1");
        let query = variant(100, "T", "A", "1/1");
        assert_eq!(compute_paired_bk(&truth, &query), "am");
    }

    /// Class C pin (post-#79): chr21:30374431-435 cluster signatures.
    /// Truth has multi-position multi-allelic (insert + multi-allelic
    /// G→GT,T) and query has overlapping insert + subst single-allelic
    /// records. The narrow relaxation (truth's alts cover both query
    /// alleles at the conflict pos) lets the query enumeration produce
    /// the truth-matching haplotype pair.
    #[test]
    fn class_c_cluster_signatures_match_after_relaxation() {
        let mut reference = vec![b'N'; 30374450];
        let window = b"ggccTAATTTGTTTTTTTTTT";
        for (i, b) in window.iter().enumerate() {
            reference[30374425 - 1 + i] = *b;
        }
        let reference = String::from_utf8(reference).unwrap();
        let truth = vec![
            variant(30374431, "A", "AT", "1|0"),
            variant(30374435, "G", "GT,T", "2|1"),
        ];
        let query = vec![
            variant(30374435, "G", "GT", "1/1"),
            variant(30374435, "G", "T", "0/1"),
        ];
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 30374431,
            end: 30374435,
            truth: truth.clone(),
            query: query.clone(),
        };
        let relax = compute_class_c_relaxation_positions(&query, &truth);
        assert_eq!(relax, BTreeSet::from([30374435usize]));
        let truth_sig =
            cluster_signature(&cluster, &truth, &reference, None, &BTreeSet::new()).unwrap();
        let query_sig = cluster_signature(&cluster, &query, &reference, None, &relax).unwrap();
        assert!(truth_sig.is_some() && query_sig.is_some());
        assert!(
            truth_sig
                .as_ref()
                .unwrap()
                .intersection(query_sig.as_ref().unwrap())
                .next()
                .is_some()
        );
    }

    /// Class C negative pin: chr21:16328989 — truth `G→GA 1|1` (homalt
    /// insert) vs query `G→GA 1/1 + G→A 0/1`. Truth's alts {GA} do not
    /// cover the SNP `A` so the relaxation must NOT fire; legacy keeps
    /// the strict drain semantics → BK=`.` on the FP query record.
    #[test]
    fn class_c_relaxation_skips_insert_only_truth_counterpart() {
        let truth = vec![variant(16328989, "G", "GA", "1|1")];
        let query = vec![
            variant(16328989, "G", "GA", "1/1"),
            variant(16328989, "G", "A", "0/1"),
        ];
        let relax = compute_class_c_relaxation_positions(&query, &truth);
        assert!(
            relax.is_empty(),
            "truth must include both Insert and Subst alleles for relaxation"
        );
    }

    fn variant(pos: usize, r: &str, alt: &str, gt: &str) -> Variant {
        Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos,
                ref_allele: r.to_string(),
                alt_allele: alt.to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: gt.to_string(),
        }
    }

    #[test]
    fn combined_tp_uses_max_call_qual_while_qq_stays_query_sourced() {
        let truth = variant(25, "A", "G", "1/1").with_qual("60");
        let query = variant(25, "A", "G", "1/1").with_qual("55");
        let row = tp_combined_row(&truth, &query, "A", 25, "");
        assert_eq!(row.record.qual, "60");
        assert!(row.record.samples[0].ends_with(":55"));
        assert!(row.record.samples[1].ends_with(":55"));

        let higher_query = query.with_qual("65");
        let row = tp_combined_row(&truth, &higher_query, "A", 25, "");
        assert_eq!(row.record.qual, "65");
    }

    #[test]
    fn truth_halfcalls_render_with_legacy_unknown_allele_contract() {
        let truth = variant(984495, "C", ".", "0|.").with_qual("30");

        let n_row = fn_row(&truth, "C", 984494, ";Regions=CONF,TS_contained", ".");
        assert_eq!(n_row.record.ref_allele, "C");
        assert_eq!(n_row.record.alt_allele, ".");
        assert_eq!(
            n_row.record.info,
            "END=984495;BS=984494;Regions=CONF,TS_contained"
        );
        assert_eq!(n_row.record.samples[0], "0|.:N:.:.:UNK:halfcall:.");

        let tp_row = tp_single_side_row(
            &truth,
            "C",
            984494,
            ";Regions=CONF,TS_contained",
            Side::Truth,
            Some("30"),
            ".",
        );
        assert_eq!(
            tp_row.record.info,
            "END=984495;BS=984494;Regions=CONF,TS_contained"
        );
        assert_eq!(tp_row.record.samples[0], "0|.:TP:gm:.:UNK:halfcall:30");
        let mut tallies = FoldedComparisonReports::default();
        tallies.observe(&tp_row);
        assert!(tallies.all_counts.is_empty());
        assert!(tallies.pass_counts.is_empty());

        // A halfcall carries its block's BK like every other record in the
        // block: legacy's quantify derives BK from the block `ctype`, so a
        // `*` allele inside a `hap:mismatch` block is `UNK/lm`, not `UNK/.`.
        let unk_row = unk_truth_row(&truth, "C", 984494, "", "lm");
        assert_eq!(unk_row.record.info, "END=984495;BS=984494");
        assert_eq!(unk_row.record.samples[0], "0|.:UNK:lm:.:UNK:halfcall:.");
        let matched_row = unk_truth_row(&truth, "C", 984494, "", ".");
        assert_eq!(matched_row.record.samples[0], "0|.:UNK:.:.:UNK:halfcall:.");
    }

    #[test]
    fn exact_only_unphased_indel_block_reaches_legacy_hap_match_verdict() {
        let identical = Cluster {
            chrom: "chr21".to_string(),
            start: 15006495,
            end: 15006495,
            truth: vec![variant(15006495, "A", "ATCTC", "0/1")],
            query: vec![variant(15006495, "A", "ATCTC", "0/1")],
        };
        assert_eq!(
            identical_gt_exact_indel_keys(&identical),
            BTreeSet::from([identical.truth[0].key.clone()])
        );

        let mut phased_truth = identical.clone();
        phased_truth.truth[0].gt = "1|0".to_string();
        assert!(
            identical_gt_exact_indel_keys(&phased_truth).is_empty(),
            "the standard phased-truth fixture must retain its legacy lm verdict"
        );
    }

    #[test]
    fn complex_subtype_uses_net_length_change_at_bucket_boundaries() {
        assert_eq!(
            subtype_label(&variant(100, "C", "GAGGTA", "0/1")),
            Some("C1_5,tv".to_string())
        );
        assert_eq!(
            subtype_label(&variant(100, "TTTAGT", "A", "0/1")),
            Some("C1_5,tv".to_string())
        );
        assert_eq!(
            subtype_label(&variant(100, "G", "TAATTTTTAAATTTTT", "0/1")),
            Some("C6_15,tv".to_string())
        );
    }

    #[test]
    fn legacy_graph_order_puts_records_paired_across_inputs_first() {
        let insertion = variant(100, "A", "AT", "0|1");
        let snp = variant(100, "A", "T", "1|1");
        let later = variant(120, "A", "G", "0/1");
        // htslib pairs on allele spelling: the insertion is the only key the
        // other input also carries, so it leads its position regardless of
        // which line came first in this file.
        let paired = BTreeSet::from([insertion.key.clone()]);

        let truth_order =
            legacy_graph_truth_order(&[snp.clone(), insertion.clone(), later.clone()], &paired);
        assert_eq!(
            truth_order
                .iter()
                .map(|v| v.key.clone())
                .collect::<Vec<_>>(),
            vec![insertion.key.clone(), snp.key.clone(), later.key.clone()]
        );

        // Truth keeps raw input order inside a class: nothing paired here, so
        // the insertion stays behind the substitution it followed in the file.
        let raw = legacy_graph_truth_order(
            &[snp.clone(), insertion.clone(), later.clone()],
            &BTreeSet::new(),
        );
        assert_eq!(
            raw.iter().map(|v| v.key.clone()).collect::<Vec<_>>(),
            vec![snp.key.clone(), insertion.key.clone(), later.key.clone()]
        );

        // The query-only remainder instead follows pre.py's trimmed-start
        // order, so the substitution at POS precedes the insertion anchored
        // on it even when the file listed them the other way round.
        let query_order = legacy_graph_query_order(
            &[insertion.clone(), snp.clone(), later.clone()],
            &BTreeSet::new(),
            &[],
        );
        assert_eq!(
            query_order
                .iter()
                .map(|v| v.key.clone())
                .collect::<Vec<_>>(),
            vec![snp.key.clone(), insertion.key.clone(), later.key]
        );

        // Paired query records take the truth stream's order instead.
        let both = BTreeSet::from([snp.key.clone(), insertion.key.clone()]);
        let truth = legacy_graph_truth_order(&[insertion.clone(), snp.clone()], &both);
        let paired_query =
            legacy_graph_query_order(&[snp.clone(), insertion.clone()], &both, &truth);
        assert_eq!(
            paired_query
                .iter()
                .map(|v| v.key.clone())
                .collect::<Vec<_>>(),
            vec![insertion.key, snp.key]
        );
    }

    #[test]
    fn outside_conf_exact_indel_keeps_local_mismatch_with_residual_allele() {
        let reference = BTreeMap::from([("chr21".to_string(), "A".repeat(256))]);
        let exact_truth = variant(100, "A", "AT", "1|1");
        let exact_query = variant(100, "A", "AT", "1/1");
        let residual_truth = variant(105, "A", "AT,AG", "1|1");
        let residual_query = variant(105, "A", "AT,AG", "2/2");
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 105,
            truth: vec![exact_truth, residual_truth],
            query: vec![exact_query, residual_query],
        };
        let mut counts = BTreeMap::new();
        let mut subtype_counts = BTreeMap::new();
        let mut rows = Vec::new();
        process_cluster(
            &cluster,
            &reference,
            Some(&[]),
            ComparisonConfig {
                no_hc: false,
                max_enum: 100_000,
                hb_expand: 0,
            },
            &mut counts,
            &mut subtype_counts,
            &mut rows,
        )
        .unwrap();

        let exact = rows
            .iter()
            .find(|row| row.record.raw().pos == 100)
            .expect("exact pair is emitted");
        assert!(
            exact.record.samples_contain(":UNK:lm:"),
            "rows={:?}",
            rows.iter()
                .map(|row| row.record.raw().to_line())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn repetitive_indel_block_promotes_legacy_hap_matched_gt_mismatch() {
        let earlier_truth = variant(15181523, "A", "AT", "0/1");
        let paired_truth = variant(15181526, "A", "AT", "0/1");
        let shared_truth_snp = variant(15181526, "A", "T", "0/1");
        let paired_query = variant(15181526, "A", "AT", "1/1");
        let shared_query_snp = variant(15181526, "A", "T", "0/1");
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 15181523,
            end: 15181526,
            truth: vec![
                earlier_truth.clone(),
                paired_truth.clone(),
                shared_truth_snp.clone(),
            ],
            query: vec![paired_query.clone(), shared_query_snp.clone()],
        };
        let region_state = RegionState {
            conf_enabled: true,
            any_conf: true,
            any_nonconf: true,
            covered_truth: BTreeSet::from([paired_truth.key.clone(), shared_truth_snp.key.clone()]),
            covered_query: BTreeSet::from([paired_query.key.clone(), shared_query_snp.key.clone()]),
            ..RegionState::default()
        };

        assert_eq!(
            legacy_repetitive_indel_hap_promotions(&cluster, &region_state),
            BTreeSet::from([paired_truth.key.clone()])
        );
        assert_eq!(
            legacy_preprocessed_snp_first_positions(&cluster),
            BTreeSet::from([15181526])
        );

        let reference = "A".repeat(4);
        let mut rows = vec![
            tp_combined_row(&paired_truth, &paired_query, &reference, cluster.start, ""),
            tp_combined_row(
                &shared_truth_snp,
                &shared_query_snp,
                &reference,
                cluster.start,
                "",
            ),
        ];
        apply_legacy_combined_before_truth_only_order(&mut rows);
        apply_legacy_preprocessed_snp_first_order(&mut rows, &cluster);
        rows.sort_by_key(|row| row.sort_key.clone());
        assert_eq!(rows[0].record.raw().alt_allele, "T");
        assert_eq!(rows[1].record.raw().alt_allele, "AT");

        let all_conf = RegionState {
            covered_truth: cluster.truth.iter().map(|v| v.key.clone()).collect(),
            ..region_state.clone()
        };
        assert!(
            legacy_repetitive_indel_hap_promotions(&cluster, &all_conf).is_empty(),
            "promotion requires the balancing truth copy outside CONF"
        );

        let mut phased = cluster.clone();
        phased.truth[1].gt = "0|1".to_string();
        assert!(legacy_preprocessed_snp_first_positions(&phased).is_empty());
    }

    #[test]
    fn xcmp_excludes_filtered_truth_after_preprocessing() {
        let pass = variant(100, "A", "G", "1|1");
        let mut filtered = variant(101, "C", "T", "1|1");
        filtered.filter = "OverlapConflict".to_string();
        let mut variants = vec![pass.clone(), filtered];

        retain_xcmp_truth_calls(&mut variants);

        assert_eq!(variants.len(), 1);
        assert_eq!(variants[0].key, pass.key);
    }

    #[test]
    fn filtered_truth_counterpart_sorts_first_at_shared_locus() {
        let row = |alt: &str, side_rank| AnnotatedRow {
            sort_key: ("chr21".to_string(), 15576177, side_rank, 0),
            record: comparison_record(&format!(
                "chr21\t15576177\t.\tG\t{alt}\t0\t.\tBS=15576177\tGT\t./.\t0/1"
            )),
            query_pass: true,
            fp_class: None,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        };
        let mut rows = vec![row("GAAAGAA", 0), row("A", 2)];
        let keys = BTreeSet::from([VariantKey {
            chrom: "chr21".to_string(),
            pos: 15576177,
            ref_allele: "G".to_string(),
            alt_allele: "A".to_string(),
        }]);

        sort_comparison_rows(&mut rows, &keys);

        assert!(rows[0].record.alt_allele == "A");
    }

    #[test]
    fn cluster_gate_rejects_snp_only_multiallelic_cluster() {
        // chr21:15027483 fixture: truth T→C,TATC (GT 1|0) and query T→C
        // (GT 0/1). The GT-active alt on truth is C (SNP) and on query is
        // C (SNP); TATC is declared but not selected. Legacy's xcmp skips
        // block-level hapcmp in this setup because `n_nonsnp == 0` across
        // all GT-selected alts. The haplotype match that would otherwise
        // rescue this as TP must therefore be suppressed so the unmatched
        // pair stays as FN + FP.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 15027483,
            end: 15027483,
            truth: vec![variant(15027483, "T", "C,TATC", "1|0")],
            query: vec![variant(15027483, "T", "C", "0/1")],
        };
        assert!(!cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn cluster_gate_admits_cluster_with_gt_selected_indel() {
        // chr21:15006495 fixture: truth A→ATCTC,ATC (GT 1|0) selects ATCTC
        // (4-bp insert); query A→ATCTC (GT 0/1) also selects an insert.
        // Hapcmp must run so the haplotype rescue lets this pair through
        // as TP.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 15006495,
            end: 15006499,
            truth: vec![variant(15006495, "A", "ATCTC,ATC", "1|0")],
            query: vec![variant(15006495, "A", "ATCTC", "0/1")],
        };
        assert!(cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn cluster_gate_ignores_non_selected_insertion() {
        // Truth A→G,AT (GT=1|0): GT selects only the SNP G; AT is a
        // non-selected insertion. Non-selected insertions do NOT trigger
        // hapcmp — only non-selected deletions do.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 100,
            truth: vec![variant(100, "A", "G,AT", "1|0")],
            query: vec![variant(100, "A", "G", "0/1")],
        };
        assert!(!cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn cluster_gate_admits_non_selected_deletion() {
        // Truth TA→AA,T (GT=1|0): GT selects only the SNP AA; T is a
        // non-selected deletion allele. Legacy's xcmp counts this record
        // toward n_nonsnp and runs hapcmp, so we must too.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 36476244,
            end: 36476244,
            truth: vec![variant(36476244, "TA", "AA,T", "1|0")],
            query: vec![variant(36476244, "T", "A", "0/1")],
        };
        assert!(cluster_has_gt_selected_nonsnp(&cluster));
    }

    #[test]
    fn trimmed_primitive_lens_collapses_shared_prefix_and_suffix() {
        // TAT→CAT trims to T→C after prefix/suffix removal — still a SNP.
        assert_eq!(trimmed_primitive_lens("TAT", "CAT"), (1, 1));
        // A→AT is a 1-bp insert after trimming the shared A.
        assert_eq!(trimmed_primitive_lens("A", "AT"), (0, 1));
        // AT→A is a 1-bp delete after trimming.
        assert_eq!(trimmed_primitive_lens("AT", "A"), (1, 0));
    }

    #[test]
    fn exact_match_key_rejects_reordered_multiallelic_alts() {
        // chr21:18743964 fixture: truth declares `GAA→GA,G` (GT 2|1) and
        // query declares `GAA→G,GA` (GT 1/2). Under BTreeSet-of-alts
        // semantics both sides share alleles {GA, G}, but legacy's
        // simpleCompare operates on byte-level ALT strings and sees them
        // as distinct (column order differs). The multi-allelic pair
        // must NOT exact-match — when the enclosing block also carries a
        // truth-only C→T 1|1 SNP, the block-wide hap signature is
        // disjoint and everything in the cluster falls to FN/FP.
        let truth = variant(18743964, "GAA", "GA,G", "2|1");
        let query = variant(18743964, "GAA", "G,GA", "1/2");
        assert!(
            !query_matches_truth_key(&query, &truth),
            "multi-allelic ALTs with different column order must not pair"
        );
    }

    #[test]
    fn exact_match_key_pairs_identical_multiallelic_alts_regardless_of_gt_phase() {
        // chr21:16032497 fixture: truth `C→CA,CAA` GT 2|1 (phased) and
        // query `C→CA,CAA` GT 1/2 (unphased). Raw ALT columns are
        // byte-equal; GT multisets are equal. Pairs as TP with
        // `equivalent_gt` doing the multiset check at the call site.
        let truth = variant(16032497, "C", "CA,CAA", "2|1");
        let query = variant(16032497, "C", "CA,CAA", "1/2");
        assert!(query_matches_truth_key(&query, &truth));
        assert!(equivalent_gt(&truth.gt, &query.gt));
    }

    #[test]
    fn bk_path_different_alt_at_same_pos_emits_dot() {
        // Legacy's VariantReader groups records by the full (chrom, pos,
        // ref, alt) tuple. Truth C→T and query C→G at the same position
        // have DIFFERENT alt columns, so they stay as separate Variants
        // objects; compareVariants returns kind=missing for each →
        // BK=. (not almismatch/lm).  The synth-snp-mismatch fixture
        // (generated by legacy) confirms this: G→T truth vs G→A query →
        // BK=. for both FN and FP rows.
        let truth = variant(15181526, "C", "T", "0|1");
        let query_counterpart = variant(15181526, "C", "G", "0/1");
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&query_counterpart), false),
            "."
        );
    }

    #[test]
    fn bk_path_almismatch_emits_lm_on_same_alt_disjoint_gt_selection() {
        // almismatch fires when both sides share the same (chrom, pos,
        // ref, alt) record — i.e. the same multi-allelic VCF entry — but
        // each sample's GT selects a completely disjoint non-ref allele
        // subset. Truth selects {T} from "T,A", query selects {A} →
        // disjoint → BK=lm.
        let truth = variant(15181526, "C", "T,A", "1/1");
        let query_counterpart = variant(15181526, "C", "T,A", "2/2");
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&query_counterpart), false),
            "lm"
        );
    }

    #[test]
    fn bk_path_hap_mismatch_emits_lm_without_same_locus_counterpart() {
        // chr21:15313079 pattern: truth FN at pos P; the cluster's
        // query side has a GT-selected indel at a different locus (not
        // same-locus-to-FN). No almismatch fires, but the block-level
        // haplotype comparator ran and the two signatures disagreed —
        // legacy's ctype="hap:mismatch" → BK=lm. Confidence-region
        // membership of the counterpart is irrelevant.
        let truth = variant(15313088, "A", "G", "0|1");
        let query_indel = variant(15313079, "C", "CA", "0/1");
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&query_indel), true),
            "lm"
        );
    }

    #[test]
    fn bk_path_long_persisted_aggregate_keeps_missing_kind() {
        let long = variant(
            1533412,
            "C",
            &format!("C{},C{}", "A".repeat(600), "A".repeat(700)),
            "2/1",
        );
        assert_eq!(bk_for_row(&long, &[], true), ".");
    }

    #[test]
    fn bk_path_fallthrough_emits_dot_on_snp_only_mismatch() {
        // chr21:15200371 pattern: truth FN at pos P; cluster's only
        // query is a SNP 7bp away at a different locus. No same-locus
        // counterpart means no almismatch; SNP-only cluster means
        // legacy's hap-run gate (n_nonsnp>0) didn't fire, so
        // hap:mismatch is also false. Legacy falls through to BK=.
        // Proximity alone never promotes BK to lm.
        let truth = variant(15200371, "T", "C", "0|1");
        let neighbour_snp = variant(15200378, "T", "C", "0/1");
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&neighbour_snp), false),
            "."
        );
    }

    #[test]
    fn bk_path_filtered_counterpart_not_almismatch() {
        // chr21:15007500 pattern: filtered-out query at the exact
        // same (chrom, pos, ref) as a truth FN. Legacy's simple-compare
        // runs only on post-filter calls, so a filtered counterpart
        // never enters `alleles_seen_2` — almismatch cannot fire.
        // Rust must skip filtered counterparts with the same guard.
        let truth = variant(15007500, "C", "T", "1|0");
        let mut filtered_query = variant(15007500, "C", "G", "0/1");
        filtered_query.filter = "LowQual".to_string();
        assert_eq!(
            bk_for_row(&truth, std::slice::from_ref(&filtered_query), false),
            "."
        );
    }

    // Residual #49 — chr21:17566241 exact-match pair with reordered
    // multi-allelic ALT columns. Truth `C→CA,CAA` (phased 1|2) and query
    // `C→CAA,CA` (unphased 1/2) encode the same diploid allele set.
    // Legacy's simpleCompare matches them via the VariantReader's shared
    // allele-unification table and emits one combined TP:gm row; rust
    // used to fall through to `cluster_signature` and split the pair
    // into truth-only FN + query-only FP. The `simple_compare_pairs_match`
    // predicate plus `canonical_hetalt_gt` close this gap.
    #[test]
    fn simple_compare_matches_reordered_multiallelic_hetalt_indel() {
        let truth = variant(17566241, "C", "CA,CAA", "1|2");
        let query = variant(17566241, "C", "CAA,CA", "1/2");
        let reference = "N".repeat(17566250);
        assert!(simple_compare_pairs_match(
            &truth,
            &query,
            &reference,
            17566241,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
        // Query GT must be remapped into truth's ALT index space —
        // `1` (CAA) → truth idx 2, `2` (CA) → truth idx 1, so `1/2` → `2/1`.
        assert_eq!(canonical_hetalt_gt(&truth.key.alt_allele, &query), "2/1");
    }

    #[test]
    fn reordered_multiallelic_genotype_mismatch_uses_truth_allele_indices() {
        // happy:chr21 at 38861935. Legacy unifies the reordered ALT columns
        // into truth order and emits one combined FN/FP row. Query `2/1`
        // against `TA,TAA` therefore becomes `1/2` against `TAA,TA`.
        let truth = variant(38_861_935, "T", "TAA,TA", "1|1");
        let query = variant(38_861_935, "T", "TA,TAA", "2/1");

        assert!(!query_matches_truth_key(&query, &truth));
        assert!(query_matches_truth_allele_set(&query, &truth));
        assert_ne!(
            selected_alt_sequences(&truth),
            selected_alt_sequences(&query)
        );
        assert_eq!(canonical_hetalt_gt(&truth.key.alt_allele, &query), "1/2");
        assert_eq!(compute_paired_bk(&truth, &query), "lm");
    }

    // Class A pin (post-#79): chr21:27249918 truth-subset match. Truth
    // `CTAAATAAA→C` GT 1|0 selects {C}; query `CTAAATAAA→C,CTAAA` GT 1/2
    // selects {C, CTAAA}. Truth's {C} ⊊ query's selected, the C allele
    // matches between sides, and the multi-allelic query primitive-splits
    // (CTAAA trims to ATAAA→A at pos+4, distinct from the C primitive at
    // pos). Legacy emits a combined TP/gm row at truth's representation
    // plus a residual FP at the orphan primitive.
    #[test]
    fn truth_subset_match_fires_for_chr21_27249918_shape() {
        let truth = variant(27249918, "CTAAATAAA", "C", "1|0");
        let query = variant(27249918, "CTAAATAAA", "C,CTAAA", "1/2");
        let reference = "N".repeat(27249930);
        assert!(truth_subset_match(
            &truth,
            &query,
            &reference,
            27249918,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
        // Query GT remap into truth's index space: query allele 1 (C)
        // matches truth idx 1; query allele 2 (CTAAA) drops to ref 0.
        // Unphased canonicalisation places the smaller index first.
        assert_eq!(remap_query_gt_subset(&truth, &query), "0/1");
    }

    #[test]
    fn truth_subset_match_rejects_homalt_truth_against_hetalt_query() {
        // chr21:10716541 shape: truth `C→G` GT 1|1 (homalt, multiset
        // [G×2]) vs query `C→A,G` GT 2/1 (hetalt, multiset [G×1, A×1]).
        // Set-subset would match {G} ⊆ {A,G} but the multiset check
        // rejects: truth needs G twice, query has it once.
        let truth = variant(10716541, "C", "G", "1|1");
        let query = variant(10716541, "C", "A,G", "2/1");
        let reference = "N".repeat(10716550);
        assert!(!truth_subset_match(
            &truth,
            &query,
            &reference,
            10716541,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
    }

    #[test]
    fn truth_subset_match_rejects_same_anchor_multiallelic_query() {
        // chr21:40096658 shape: truth `T→TAGATAGAG` GT 1|0 vs query
        // `T→TAGATAGAG,TAGATAGAT` GT 1/2. Both query alts trim to the
        // same anchor (T at pos), so `query_primitive_splits` is false
        // and legacy keeps two separate rows rather than emitting a
        // combined TP. The gate must reject this case.
        let truth = variant(40096658, "T", "TAGATAGAG", "1|0");
        let query = variant(40096658, "T", "TAGATAGAG,TAGATAGAT", "1/2");
        let reference = "N".repeat(40096670);
        assert!(!truth_subset_match(
            &truth,
            &query,
            &reference,
            40096658,
            std::slice::from_ref(&truth),
            std::slice::from_ref(&query),
        ));
    }

    #[test]
    fn truth_subset_match_preserves_persisted_location_aggregate() {
        // test_full chr1:963700: preprocessing has already combined the
        // query's two deletion calls into one canonical `2/1` record.
        // Legacy hap.py emits both truth rows plus this one query row; it
        // does not consume one query allele into a combined truth row.
        let truth_short = variant(963700, "GC", "G", "1|0");
        let truth_long = variant(963700, "GCC", "G", "0|1");
        let query = variant(963700, "GCC", "GC,G", "2/1");
        let reference = "N".repeat(963710);
        assert!(!truth_subset_match(
            &truth_long,
            &query,
            &reference,
            963700,
            &[truth_short, truth_long.clone()],
            std::slice::from_ref(&query),
        ));
    }

    // Residual #50 — chr21:15671076 cluster had a truth `T→TATATA` at
    // pos 15671094 plus a truth `T→TA` at pos 15671095. Both are
    // insertions. The previous homopolymer anchor slide pushed both
    // anchors to the same cluster_end-1 position, and `apply_events`
    // then rejected the pair as "two different inserts at the same
    // anchor". Restricting the slide to true homopolymer-extension
    // inserts (seq consists entirely of the anchor base) keeps the
    // distinct events at distinct anchors so the cluster signature
    // resolves and the block fires BK=lm.
    #[test]
    fn normalize_ref_alt_does_not_slide_heterogeneous_insert() {
        // Reference has a T-homopolymer at pos 15671094..15671099.
        // A 5-base `ATATA` insertion must stay at its original anchor
        // 15671094 rather than sliding through the T-homopolymer —
        // the previous slide collapsed distinct cluster inserts to the
        // same anchor, triggering a false conflict in `apply_events`.
        let mut reference = vec![b'N'; 15671110];
        // Lay down "atatatatatatattttttt" starting at pos 15671080.
        let window = b"atatatatatatattttttt";
        for (i, b) in window.iter().enumerate() {
            reference[15671080 - 1 + i] = *b;
        }
        let reference = String::from_utf8(reference).unwrap();
        let events = normalize_ref_alt(15671094, "T", "TATATA", &reference, 15671076, 15671097);
        assert_eq!(events.len(), 1);
        match &events[0] {
            Event::Insert { anchor, seq } => {
                assert_eq!(
                    *anchor, 15671094,
                    "anchor must not slide for non-homopolymer insert"
                );
                assert_eq!(seq, "ATATA");
            }
            other => panic!("expected Insert, got {other:?}"),
        }
    }

    // Residual #50 companion — chr21:16997925 cluster had a query
    // homalt deletion `GCA→G` at pos 16997949 (delete pos 16997950-51)
    // plus a query het insert `A→ACG` at pos 16997951. The insert's
    // anchor (16997951) falls inside the delete's span, so
    // `apply_events` returned `None`. Legacy applies downstream
    // insertions against the shifted (post-delete) sequence without
    // failing; rust now allows insertions at deleted anchors since the
    // output-walk loop emits them unconditionally.
    #[test]
    fn apply_events_allows_insert_anchored_inside_delete() {
        // Reference must span the whole cluster — padding with Ns up to
        // position 16997960 and placing the CACA context inline.
        let mut reference = vec![b'N'; 16997960];
        reference[16997948] = b'G'; // pos 16997949: G
        reference[16997949] = b'C'; // pos 16997950: C (deleted)
        reference[16997950] = b'A'; // pos 16997951: A (deleted, insert anchor)
        reference[16997951] = b'G'; // pos 16997952: G
        let reference = String::from_utf8(reference).unwrap();
        let events = vec![
            Event::Delete {
                start: 16997950,
                end: 16997951,
            },
            Event::Insert {
                anchor: 16997951,
                seq: "CG".to_string(),
            },
        ];
        let result = apply_events(&reference, 16997949, 16997952, &events).unwrap();
        assert!(
            result.is_some(),
            "delete + downstream insert must produce a valid haplotype"
        );
    }

    // Class 4 pin: BI (comparison_info)
    // on multi-allelic SNPs with mixed ti/tv alleles must emit the
    // comma-joined per-allele tokens legacy uses, not a single
    // collapsed ti/tv. Chr21:17562906 `A → G,T GT=2/1`: A→G is ti,
    // A→T is tv; legacy emits `ti,tv`, rust previously collapsed to
    // `tv` via the SNP-single-alt fastpath.
    #[test]
    fn comparison_info_joins_multi_allelic_snp_ti_tv() {
        let var = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 17562906,
                ref_allele: "A".to_string(),
                alt_allele: "G,T".to_string(),
            },
            qual: "1684.9".to_string(),
            filter: "PASS".to_string(),
            gt: "2/1".to_string(),
        };
        assert_eq!(comparison_info(&var, "N"), "ti,tv");
    }

    #[test]
    fn ambiguous_reference_snp_has_no_titv_subtype() {
        let variant = variant(121965037, "N", "T", "1/1");
        assert_eq!(comparison_info(&variant, "N"), ".");
        assert_eq!(snp_bucket_label(&variant), None);
    }

    #[test]
    fn spanning_deletion_halfcall_uses_anchor_for_confidence() {
        let halfcall = variant(100, "ACGT", ".", "0|.");
        let conf = [Interval {
            chrom: "chr21".to_string(),
            start: 99,
            end: 100,
        }];
        assert!(variant_is_conf(&halfcall, "N", 100, 103, &conf));
    }

    #[test]
    fn halfcall_is_tp_only_inside_a_matched_deletion() {
        let deletion = variant(100, "ACGT", "A", "0|1");
        let query = variant(100, "ACGT", "A", "0/1");
        let mut cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 103,
            truth: vec![deletion, variant(102, "G", ".", "0|.")],
            query: vec![query],
        };
        assert!(!halfcall_is_covered_by_matched_deletion(
            &cluster.truth[1],
            &cluster
        ));
        cluster.query.push(variant(103, "T", "TA,TAA", "2/1"));
        assert!(halfcall_is_covered_by_matched_deletion(
            &cluster.truth[1],
            &cluster
        ));
        let outside = variant(104, "T", ".", "0|.");
        assert!(!halfcall_is_covered_by_matched_deletion(&outside, &cluster));
    }

    #[test]
    fn unmatched_covering_deletion_promotes_halfcall_in_hap_match() {
        let deletion = variant(100, "ACGT", "A", "0|1");
        let halfcall = variant(102, "G", ".", "0|.");
        let cluster = Cluster {
            chrom: "chr10".to_string(),
            start: 100,
            end: 103,
            truth: vec![deletion, halfcall.clone()],
            query: vec![variant(102, "G", "T", "1/1")],
        };

        assert!(halfcall_is_covered_by_matched_deletion(&halfcall, &cluster));
    }

    #[test]
    fn overlapping_truth_deletion_makes_gt_discordance_local() {
        let paired_truth = variant(127, "ACGT", "A", "1|0");
        let paired_query = variant(127, "ACGT", "A", "1/1");
        let overlapping_truth = variant(100, "A".repeat(50).as_str(), "A", "0|1");
        let cluster = Cluster {
            chrom: "chr20".to_string(),
            start: 100,
            end: 149,
            truth: vec![overlapping_truth, paired_truth],
            query: vec![paired_query],
        };

        assert!(legacy_overlapping_deletion_mismatch(&cluster));
    }

    #[test]
    fn later_truth_deletions_inside_paired_span_keep_haplotype_match() {
        let paired_truth = variant(100, "GGA", "G", "0/1");
        let paired_query = variant(100, "GGA", "G", "1/1");
        let later_truth = variant(102, "AG", "A", "1/0");
        let cluster = Cluster {
            chrom: "chr13".to_string(),
            start: 100,
            end: 103,
            truth: vec![paired_truth, later_truth],
            query: vec![paired_query],
        };

        assert!(!legacy_overlapping_deletion_mismatch(&cluster));
    }

    #[test]
    fn hapfail_does_not_propagate_row_mismatch_to_exact_rows() {
        assert!(should_propagate_unreconciled_exact_rows(
            "hap:mismatch",
            true,
            false,
        ));
        assert!(!should_propagate_unreconciled_exact_rows(
            "hapfail:mismatch",
            true,
            false,
        ));
        assert!(should_propagate_unreconciled_exact_rows(
            "hapfail:mismatch",
            true,
            true,
        ));
        assert!(!should_propagate_unreconciled_exact_rows(
            "hap:mismatch",
            false,
            true,
        ));
    }

    #[test]
    fn exact_insert_counterpart_does_not_propagate_hapfail_mismatch() {
        let truth_insert = variant(100, "A", "ATTT", "0/1");
        let truth_snp = variant(100, "A", "T", "1/0");
        let query_insert = variant(100, "A", "ATTT", "0/1");
        let query_snp = variant(100, "A", "T", "1/1");
        let counterpart = query_insert_conflict_has_truth_counterpart(
            &[query_insert, query_snp],
            &[truth_insert, truth_snp.clone()],
            &[truth_snp],
        );

        assert_eq!(counterpart, Some(true));
        assert!(!should_propagate_unreconciled_exact_rows(
            "hapfail:mismatch",
            true,
            matches!(counterpart, Some(false)),
        ));
    }

    fn shared_insertion_with_extra_query_snp_rows(
        truth_alt: &str,
        truth_gt: &str,
        query_alt: &str,
        query_gt: &str,
    ) -> Vec<AnnotatedRow> {
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 101,
            truth: vec![variant(100, "T", truth_alt, truth_gt)],
            query: vec![
                variant(100, "T", "G", "0/1"),
                variant(100, "T", query_alt, query_gt),
                variant(101, "T", "G", "0/1"),
            ],
        };
        let reference = BTreeMap::from([("chr21".to_string(), "T".repeat(256))]);
        let mut counts = BTreeMap::new();
        let mut subtype_counts = BTreeMap::new();
        let mut rows = Vec::new();

        process_cluster(
            &cluster,
            &reference,
            None,
            ComparisonConfig {
                no_hc: false,
                max_enum: 100_000,
                hb_expand: 0,
            },
            &mut counts,
            &mut subtype_counts,
            &mut rows,
        )
        .unwrap();

        rows.iter()
            .filter(|row| row.record.raw().alt_allele == "G")
            .cloned()
            .collect()
    }

    #[test]
    fn legacy_only_shared_insertion_with_extra_query_snps_keeps_missing_block_kind() {
        // Reduced chr21:30866581 shape. The multi-allelic insertion is shared
        // with reversed ALT order, while query also carries an FP SNP at the
        // insertion anchor and another adjacent SNP. Both graph signatures can
        // be computed but disagree; legacy classifies the block as hapfail and
        // leaves the SNP rows at BK=`.`.
        let snp_rows =
            shared_insertion_with_extra_query_snp_rows("TTGGG,TTGG", "2|1", "TTGG,TTGGG", "1/2");
        assert_eq!(snp_rows.len(), 2);
        assert!(
            snp_rows
                .iter()
                .all(|row| row.record.samples_contain(":FP:.:"))
        );
    }

    #[test]
    fn normative_single_alt_shared_insertion_keeps_local_mismatch_block_kind() {
        let snp_rows = shared_insertion_with_extra_query_snp_rows("TTGGG", "0|1", "TTGGG", "0/1");
        assert_eq!(snp_rows.len(), 2);
        assert!(
            snp_rows
                .iter()
                .all(|row| row.record.samples_contain(":FP:lm:")),
            "rows={:?}",
            snp_rows
                .iter()
                .map(|row| row.record.raw().to_line())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn matched_insert_subst_conflict_does_not_suppress_adjacent_compound_het_mismatch() {
        // Broad-germline shape (reduced from GRCh37 HG002 7:130838422). Anchor
        // 100 carries a compound-het insertion pair: truth spells it as two
        // separate het records, the query as one merged hetalt. Anchor 101
        // carries an Insert+Subst pair (G>GT insert + G>T SNP) that exact-
        // matches between truth and query and drains to TP before hap-compare.
        //
        // The drained 101 conflict must NOT suppress the block's hap:mismatch
        // verdict: the genuine disagreement is the 100 compound-het pair, and
        // legacy stamps BK=lm on those leftover FN/FP rows. Before the
        // `conflict_positions_still_unmatched` gate, the matched 101 conflict
        // forced hapfail and the 100 rows lost their BK=lm.
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 101,
            truth: vec![
                variant(100, "G", "GT", "1/0"),
                variant(100, "G", "GTT", "0/1"),
                variant(101, "G", "GT", "1/0"),
                variant(101, "G", "T", "0/1"),
            ],
            query: vec![
                variant(100, "G", "GT,GTT", "2/1"),
                variant(101, "G", "GT", "0/1"),
                variant(101, "G", "T", "0/1"),
            ],
        };
        let reference = BTreeMap::from([("chr21".to_string(), "G".repeat(256))]);
        let mut counts = BTreeMap::new();
        let mut subtype_counts = BTreeMap::new();
        let mut rows = Vec::new();
        process_cluster(
            &cluster,
            &reference,
            None,
            ComparisonConfig {
                no_hc: false,
                max_enum: 100_000,
                hb_expand: 0,
            },
            &mut counts,
            &mut subtype_counts,
            &mut rows,
        )
        .unwrap();

        let anchor_rows: Vec<_> = rows
            .iter()
            .filter(|row| row.record.raw().pos == 100)
            .collect();
        assert!(!anchor_rows.is_empty(), "expected leftover rows at anchor 100");
        assert!(
            anchor_rows.iter().all(|row| {
                row.record.samples_contain(":FN:lm:") || row.record.samples_contain(":FP:lm:")
            }),
            "anchor-100 FN/FP rows must carry BK=lm, got {:?}",
            anchor_rows
                .iter()
                .map(|row| row.record.raw().to_line())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn reciprocal_query_aggregate_accepts_truth_insertion_subset() {
        let truth_insert = variant(100, "A", "ATTT", "0/1");
        let query_insert = variant(100, "A", "AT,ATTT", "2/1");
        let query_deletion = variant(100, "ATTT", "A", "0/1");

        let counterpart = query_insert_conflict_has_truth_counterpart(
            &[query_insert, query_deletion],
            &[truth_insert.clone()],
            &[truth_insert],
        );

        assert_eq!(counterpart, Some(true));
    }

    #[test]
    fn insertion_and_snp_without_truth_insertion_remain_a_conflict() {
        let truth_snp = variant(100, "A", "G", "0/1");
        let query_insert = variant(100, "A", "ATTT", "0/1");
        let query_snp = variant(100, "A", "G", "0/1");

        let counterpart = query_insert_conflict_has_truth_counterpart(
            &[query_insert, query_snp],
            &[truth_snp.clone()],
            &[truth_snp],
        );

        assert_eq!(counterpart, Some(false));
    }

    #[test]
    fn one_deletion_covering_multiple_insert_substitution_conflicts_is_local() {
        let query = vec![
            variant(100, "ACGTAC", "A", "1|1"),
            variant(103, "T", "TA", "0|1"),
            variant(103, "T", "G", "0|1"),
            variant(104, "A", "ACC", "0|1"),
            variant(104, "A", "C", "0|1"),
        ];

        assert!(legacy_covered_multi_conflict_mismatch(&query));
    }

    #[test]
    fn heterozygous_deletion_covering_multiple_conflicts_remains_hapfail() {
        let query = vec![
            variant(100, "ACGTAC", "A", "0|1"),
            variant(103, "T", "TA", "0|1"),
            variant(103, "T", "G", "0|1"),
            variant(104, "A", "ACC", "0|1"),
            variant(104, "A", "C", "0|1"),
        ];

        assert!(!legacy_covered_multi_conflict_mismatch(&query));
    }

    #[test]
    fn separate_deletions_do_not_promote_multiple_conflicts() {
        let query = vec![
            variant(100, "ACGT", "A", "0|1"),
            variant(103, "T", "TA", "0|1"),
            variant(103, "T", "G", "0|1"),
            variant(106, "A", "ACC", "0|1"),
            variant(106, "A", "C", "0|1"),
        ];

        assert!(!legacy_covered_multi_conflict_mismatch(&query));
    }

    #[test]
    fn halfcall_order_follows_legacy_companion_kind() {
        let reference = "A".repeat(256);
        let halfcall = variant(102, "A", ".", "0|.");
        let truth = variant(102, "A", "G", "0|1");
        let query = variant(102, "A", "G", "0/1");

        let mut matched = vec![
            fn_row(&halfcall, &reference, 100, "", "."),
            tp_combined_row(&truth, &query, &reference, 100, ""),
        ];
        apply_legacy_halfcall_order(&mut matched);
        assert_eq!(matched[0].sort_key.2, 2);
        assert_eq!(matched[1].sort_key.2, 1);

        let mut truth_only = vec![
            fn_row(&halfcall, &reference, 100, "", "."),
            fn_row(&truth, &reference, 100, "", "."),
        ];
        apply_legacy_halfcall_order(&mut truth_only);
        assert_eq!(truth_only[0].sort_key.2, 2);
        assert_eq!(truth_only[1].sort_key.2, 3);

        let residual_query = variant(102, "A", "AG", "0/1");
        let mut three_grains = vec![
            fn_row(&halfcall, &reference, 100, "", "."),
            fn_row(&truth, &reference, 100, "", "."),
            fp_like_row(&residual_query, &reference, 100, "", "UNK", None, "lm"),
        ];
        apply_legacy_halfcall_order(&mut three_grains);
        assert_eq!(three_grains[0].sort_key.2, 2);
        assert_eq!(three_grains[1].sort_key.2, 3);
        assert_eq!(three_grains[2].sort_key.2, 4);
    }

    #[test]
    fn combined_allele_sorts_before_truth_only_allele_at_same_position() {
        let reference = "A".repeat(256);
        let paired_truth = variant(102, "AAAA", "A", "1|0");
        let paired_query = variant(102, "AAAA", "A", "0/1");
        let truth_only = variant(102, "AAA", "A", "0|1");
        let mut rows = vec![
            fn_row(&truth_only, &reference, 100, "", "lm"),
            tp_combined_row(&paired_truth, &paired_query, &reference, 100, ""),
        ];

        apply_legacy_combined_before_truth_only_order(&mut rows);

        assert_eq!(rows[0].sort_key.2, 1);
        assert_eq!(rows[1].sort_key.2, 0);
    }

    #[test]
    fn discordant_snp_pair_keeps_same_locus_indel_rows_separate() {
        let truth_indel = variant(102, "A", "AT", "0/1");
        let query_indel = variant(102, "A", "AT", "1/1");
        let truth_snp = variant(102, "A", "G", "1/0");
        let query_snp = variant(102, "A", "G", "1/1");
        let cluster = Cluster {
            chrom: "chr6".to_string(),
            start: 100,
            end: 103,
            truth: vec![truth_indel.clone(), truth_snp],
            query: vec![query_indel.clone(), query_snp],
        };

        assert!(mixed_type_same_locus_keeps_indel_rows_separate(
            &cluster,
            &truth_indel,
            &query_indel
        ));

        let indel_only = Cluster {
            chrom: "chr6".to_string(),
            start: 100,
            end: 103,
            truth: vec![truth_indel.clone()],
            query: vec![query_indel.clone()],
        };
        assert!(!mixed_type_same_locus_keeps_indel_rows_separate(
            &indel_only,
            &truth_indel,
            &query_indel
        ));
    }

    #[test]
    fn halfcall_inherits_same_position_truth_local_mismatch_kind() {
        let reference = "A".repeat(256);
        let halfcall = variant(102, "A", ".", "0|.");
        let unmatched_snp = variant(102, "A", "T", "0|1");
        let matched_snp = variant(103, "A", "G", "0|1");
        let matched_query = variant(103, "A", "G", "0/1");

        let mut rows = vec![
            fn_row(&halfcall, &reference, 100, "", "."),
            fn_row(&unmatched_snp, &reference, 100, "", "lm"),
            tp_combined_row(&matched_snp, &matched_query, &reference, 100, ""),
        ];
        apply_legacy_halfcall_block_kind(&mut rows);

        assert!(rows[0].record.samples[0].contains(":N:lm:.:UNK:halfcall:"));
        assert!(rows[1].record.samples[0].contains(":FN:lm:"));
        assert!(rows[2].record.samples[0].contains(":TP:gm:"));
    }

    #[test]
    fn outside_conf_aggregate_mismatch_marks_the_whole_block_local() {
        let cluster = Cluster {
            chrom: "chr9".to_string(),
            start: 113463932,
            end: 113463970,
            truth: vec![
                variant(113463934, "T", "TATTTTTTTTATTGTATTGTATTG", "1|0"),
                variant(113463959, "T", "TA", "1|0"),
                variant(113463959, "T", "TATTTT", "0|1"),
            ],
            query: vec![
                variant(113463934, "T", "TATTTTTTTTATTGTATTGTATTG", "0/1"),
                variant(113463959, "T", "TA,TATTTT", "2/1"),
            ],
        };
        assert_eq!(
            legacy_unknown_aggregate_local_mismatch(&cluster, &RegionState::default(), false),
            Some(true)
        );
    }

    #[test]
    fn matched_outside_conf_reaching_insertion_aggregate_marks_whole_block_local() {
        // Reduced form of test_full chr9:113463932. The first insertion is
        // byte-equal across truth/query and reaches the later two-allele
        // aggregate. The remaining truth primitives and query aggregate
        // reconstruct the same haplotypes, but legacy's graph still stamps
        // the entire outside-CONF block BK=lm.
        let reaching = "TATTTTTTTTATTGTATTGTATTG";
        let aggregate = "TATTTTATTTTATTTTATTTTATTTTATTTTATTTT";
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 3,
            end: 28,
            truth: vec![
                variant(3, "T", reaching, "1|0"),
                variant(28, "T", "TA", "1|0"),
                variant(28, "T", aggregate, "0|1"),
            ],
            query: vec![
                variant(3, "T", reaching, "0/1"),
                variant(28, "T", &format!("TA,{aggregate}"), "2/1"),
            ],
        };
        let reference = BTreeMap::from([("chr21".to_string(), "T".repeat(96))]);
        let mut counts = BTreeMap::new();
        let mut subtype_counts = BTreeMap::new();
        let mut rows = Vec::new();
        process_cluster(
            &cluster,
            &reference,
            Some(&[]),
            ComparisonConfig {
                no_hc: false,
                max_enum: 100_000,
                hb_expand: 0,
            },
            &mut counts,
            &mut subtype_counts,
            &mut rows,
        )
        .unwrap();

        assert!(!rows.is_empty());
        assert!(
            rows.iter().all(|row| !row.record.samples_contain(":UNK:.:")
                && row.record.samples_contain(":UNK:lm:")),
            "rows={:?}",
            rows.iter()
                .map(|row| row.record.raw().to_line())
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn outside_conf_aggregate_with_matching_phases_keeps_block_kind_missing() {
        let cluster = Cluster {
            chrom: "chr2".to_string(),
            start: 1533412,
            end: 1534186,
            truth: vec![
                variant(1533412, "C", "CAAA", "0|1"),
                variant(1533412, "C", "CAAAAA", "1|0"),
                variant(1533413, "CAAAAAA", "C", "1|1"),
            ],
            query: vec![
                variant(1533412, "C", "CAAA,CAAAAA", "2/1"),
                variant(1533413, "CAAAAAA", "C", "1/1"),
            ],
        };
        assert_eq!(
            legacy_unknown_aggregate_local_mismatch(&cluster, &RegionState::default(), true),
            Some(false)
        );
    }

    #[test]
    fn symbolic_output_ref_uses_the_reference_base() {
        let variant = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 2,
                ref_allele: "N".to_string(),
                alt_allele: "<DEL>".to_string(),
            },
            qual: "0".to_string(),
            filter: ".".to_string(),
            gt: "1|0".to_string(),
        };
        assert_eq!(display_ref(&variant, "aTg"), "T");
    }

    #[test]
    fn fully_nonconf_matched_fanout_keeps_fallback_bk() {
        assert_eq!(matched_query_unk_bk(true, false, "."), ".");
        assert_eq!(matched_query_unk_bk(true, true, "."), ".");
        assert_eq!(matched_query_unk_bk(false, false, "."), ".");
        assert_eq!(matched_query_unk_bk(true, false, "lm"), "lm");
    }

    // Class 3 pin: `variant_is_conf` must apply legacy's
    // `!is_pure_insertion || fully_covered` rule using gvcf2bed-style
    // refrange. A pure insertion at a CONF edge (anchor in, anchor+1
    // out) must NOT be classified as covered. Chr21:15859667 `T → TA`
    // anchor 15859667 is last base of CONF `[15859657, 15859667)`,
    // anchor+1 15859668 is not covered — insertion straddles the
    // boundary, legacy emits no CONF tag, rust now agrees.
    #[test]
    fn variant_is_conf_rejects_insertion_at_bed_edge() {
        let var = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15859667,
                ref_allele: "T".to_string(),
                alt_allele: "TA".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1|0".to_string(),
        };
        let intervals = vec![Interval {
            chrom: "chr21".to_string(),
            start: 15859657,
            end: 15859667,
        }];
        // Anchor at 15859667 is in [15859657,15859667)? 15859666 < 15859667 → yes.
        // Anchor+1 at 15859668 is in? 15859667 < 15859667 → NO.
        // Partial → pure insertion → skip CONF.
        assert!(!variant_is_conf(&var, "N", 15859600, 15859700, &intervals));
    }

    // Per-primitive CONF coverage: a multi-allelic
    // query whose deletion primitive sits inside CONF but whose insertion
    // primitive straddles a CONF edge must produce a RegionState where
    // `covered_query` contains the deletion primitive's key but NOT the
    // insertion primitive's, and `any_nonconf` is true. Without this,
    // the cluster collapses to TS_contained on every record and the
    // insertion primitive incorrectly carries a CONF tag.
    //
    // Mirrors chr21:48036437 — query `AGTGTGT → AGTGTGTGT,A` GT=1/2 at
    // pos 37002776 splits into a deletion at pos 37002776 (in CONF) and
    // an insertion T→TGT at pos 37002782 (straddles a CONF gap).
    #[test]
    fn region_state_marks_multi_allelic_primitives_separately() {
        let parent = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 100,
                ref_allele: "AGT".to_string(),
                alt_allele: "AGTGT,A".to_string(),
            },
            qual: "100".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 100,
            end: 102,
            truth: vec![],
            query: vec![parent.clone()],
        };
        // CONF covers 1-based positions 100..101, leaving 102 uncovered.
        // The deletion primitive's effective ref range
        // (refstart=101, refend=102, is_pure_insertion=false) overlaps
        // CONF at position 101 → covered (subst/del has_overlap path).
        // The insertion primitive after right-anchor canonicalisation
        // emits at pos 102 with anchor T → bracket [101, 102]. Position
        // 102 is NOT in CONF → fully_covered=false → primitive must NOT
        // be marked covered.
        let intervals = vec![Interval {
            chrom: "chr21".to_string(),
            start: 99, // 0-based half-open → covers 1-based 100..101
            end: 101,
        }];
        let state = RegionState::from_cluster(&cluster, "N", Some(&intervals));
        assert!(state.any_conf, "deletion primitive must register coverage");
        // The insertion primitive's anchor falls at the CONF edge with
        // anchor+1 outside coverage — primitive must NOT be in
        // covered_query, and any_nonconf must be set.
        let primitives = split_query_primitives_with_neighbors(&parent, "N", 100, &[], &[]);
        let insertion_primitive = primitives
            .iter()
            .find(|p| p.key.alt_allele.len() > p.key.ref_allele.len())
            .expect("split_query_primitives must produce one insertion");
        assert!(
            !state.covered_query.contains(&insertion_primitive.key),
            "insertion primitive at CONF edge must NOT be in covered_query"
        );
        assert!(
            state.any_nonconf,
            "presence of an uncovered insertion primitive must mark cluster as non-CONF"
        );
    }

    #[test]
    fn overlapping_deletion_keeps_persisted_hetalt_deletion_aggregate_final() {
        let aggregate = variant(100, "CTCAACTAG", "C,CT", "1/2");
        let neighbor = variant(101, "TCAACTAGTTAAG", "T", "0/1");

        let split = split_query_primitives_with_neighbors(
            &aggregate,
            "N",
            100,
            &[aggregate.clone(), neighbor],
            &[],
        );

        assert_eq!(split.len(), 1);
        assert_eq!(split[0].key, aggregate.key);
        assert_eq!(split[0].gt, aggregate.gt);
    }

    #[test]
    fn legacy_only_duplicate_alt_query_projects_for_unmatched_classified_rows() {
        for gt in ["2/1", "1/2", "2|1"] {
            let query = variant(104, "A", "AGTGTGTGT,AGTGTGTGT", gt);
            let cluster = Cluster {
                chrom: "chr21".to_string(),
                start: 104,
                end: 104,
                truth: vec![],
                query: vec![query],
            };
            let reference = BTreeMap::from([("chr21".to_string(), "A".repeat(256))]);
            let mut rows = Vec::new();
            process_cluster(
                &cluster,
                &reference,
                None,
                ComparisonConfig {
                    no_hc: false,
                    max_enum: 100_000,
                    hb_expand: 0,
                },
                &mut BTreeMap::new(),
                &mut BTreeMap::new(),
                &mut rows,
            )
            .unwrap();
            assert_eq!(rows.len(), 1);
            let row = &rows[0];

            assert_eq!(row.record.alt_allele, "AGTGTGTGT");
            assert_eq!(row.record.sample_map(1).get("GT").unwrap(), "1/1");
            assert_eq!(row.record.sample_map(1).get("BLT").unwrap(), "homalt");
        }
    }

    #[test]
    fn legacy_only_duplicate_alt_query_projects_for_paired_classified_rows() {
        let truth = variant(104, "A", "AGTGTGTGT", "1|1");
        let query = variant(104, "A", "AGTGTGTGT,AGTGTGTGT", "2/1");
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 104,
            end: 104,
            truth: vec![truth],
            query: vec![query],
        };
        let reference = BTreeMap::from([("chr21".to_string(), "A".repeat(256))]);
        let mut rows = Vec::new();
        process_cluster(
            &cluster,
            &reference,
            None,
            ComparisonConfig {
                no_hc: false,
                max_enum: 100_000,
                hb_expand: 0,
            },
            &mut BTreeMap::new(),
            &mut BTreeMap::new(),
            &mut rows,
        )
        .unwrap();
        assert_eq!(rows.len(), 1);
        let row = &rows[0];

        assert_eq!(row.record.alt_allele, "AGTGTGTGT");
        assert_eq!(row.record.sample_map(0).get("GT").unwrap(), "1|1");
        assert_eq!(row.record.sample_map(1).get("GT").unwrap(), "1/1");
        assert_eq!(row.record.sample_map(1).get("BLT").unwrap(), "homalt");
    }

    #[test]
    fn normative_duplicate_alt_matching_representation_remains_distinct() {
        let aggregate = variant(141113704, "A", "AGTGTGTGT,AGTGTGTGT", "2/1");

        let split = split_query_primitives_with_neighbors(
            &aggregate,
            "A",
            141113702,
            std::slice::from_ref(&aggregate),
            &[],
        );

        assert_eq!(split.len(), 1);
        assert_eq!(split[0].key, aggregate.key);
        assert_eq!(split[0].gt, aggregate.gt);
    }

    #[test]
    fn normative_duplicate_alt_query_projection_leaves_other_shapes_unchanged() {
        for query in [
            variant(141113704, "A", "AGTGTGTGT,ACT", "2/1"),
            variant(141113704, "A", "AGTGTGTGT,AGTGTGTGT", "0/1"),
            variant(141113704, "A", "AGTGTGTGT,AGTGTGTGT", "1/1"),
            variant(141113704, "A", "AGTGTGTGT,AGTGTGTGT,AGTGTGTGT", "2/1"),
            variant(141113704, "A", "<DEL>,<DEL>", "2/1"),
        ] {
            let projected = legacy_duplicate_alt_query_output_projection(&query);
            assert_eq!(projected.key, query.key);
            assert_eq!(projected.gt, query.gt);
        }
    }

    // Class 3 support: SNPs at a single base inside any CONF interval
    // are covered regardless of the insertion-aware fully_covered
    // clause. Sanity-check that the refactored `variant_is_conf`
    // doesn't accidentally reject point-covered SNPs.
    #[test]
    fn variant_is_conf_accepts_snp_inside_interval() {
        let var = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15859632,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1|0".to_string(),
        };
        let intervals = vec![Interval {
            chrom: "chr21".to_string(),
            start: 15859480,
            end: 15859645,
        }];
        assert!(variant_is_conf(&var, "N", 15859600, 15859700, &intervals));
    }

    // Class 3 gvcf2bed-style padding bridges CONF gaps at insertion
    // anchors. Chr21:17562905 `C→CG` pure insertion sits at the edge
    // of CONF interval `[17561431, 17562905)` and the next interval
    // `[17562906, 17564580)` — a 1-base gap at 17562905. Legacy's
    // gvcf2bed emits `chr21 17562904 17562906` which, when merged with
    // the raw CONF bed, closes the gap and makes the adjacent A→G,T
    // SNP at pos 17562906 land inside CONF. Pin the padding function.
    #[test]
    fn gvcf2bed_padding_spans_insertion_anchor_and_next_base() {
        let truth = vec![Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 17562905,
                ref_allele: "C".to_string(),
                alt_allele: "CG".to_string(),
            },
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            gt: "1|1".to_string(),
        }];
        let padding = gvcf2bed_padding(&truth, None);
        assert_eq!(padding.len(), 1);
        assert_eq!(padding[0].chrom, "chr21");
        // 0-based half-open: anchor = 17562904 .. anchor+1+1 = 17562906
        assert_eq!(padding[0].start, 17562904);
        assert_eq!(padding[0].end, 17562906);
    }

    // gvcf2bed `-T <bed>` filter: legacy uses `bcf_sr_set_targets(.., 1, 0)`
    // which gates on the **start position only** — not the full ref span.
    // A record whose 1-based pos (→ 0-based start) lies outside every
    // raw CONF interval is dropped before emission. This is what
    // `IS_CONF.Size` parity hinges on (legacy sums per-file BED lengths
    // without cross-file dedup, so dropping out-of-target records keeps
    // the padding budget honest).
    #[test]
    fn gvcf2bed_padding_target_filter_excludes_out_of_target_record() {
        let truth = vec![
            // pos 100 → pos_0b 99, INSIDE conf [50, 150)
            Variant {
                key: VariantKey {
                    chrom: "chr1".to_string(),
                    pos: 100,
                    ref_allele: "A".to_string(),
                    alt_allele: "G".to_string(),
                },
                qual: ".".to_string(),
                filter: "PASS".to_string(),
                gt: "0/1".to_string(),
            },
            // pos 200 → pos_0b 199, OUTSIDE conf — should be dropped
            Variant {
                key: VariantKey {
                    chrom: "chr1".to_string(),
                    pos: 200,
                    ref_allele: "A".to_string(),
                    alt_allele: "G".to_string(),
                },
                qual: ".".to_string(),
                filter: "PASS".to_string(),
                gt: "0/1".to_string(),
            },
        ];
        let conf = vec![Interval {
            chrom: "chr1".to_string(),
            start: 50,
            end: 150,
        }];
        let padding = gvcf2bed_padding(&truth, Some(&conf));
        assert_eq!(padding.len(), 1, "only in-target record should emit");
        assert_eq!(padding[0].start, 99);
        assert_eq!(padding[0].end, 100);
    }

    // Legacy gvcf2bed emits a BED line **per record** unconditionally,
    // even when every alt is symbolic (`<DEL>`, `<NON_REF>`, etc.). The
    // alt loop's `break` on the first non-NUC alt leaves
    // `nuc_alleles=false`, so refstart/refend stay at the raw
    // [pos, pos+reflen-1] from getLocation — and emission proceeds. Our
    // truth fixtures contain ~14 such records on chr21 (`<DEL>` calls);
    // skipping them under-counts IS_CONF.Size by ~14 bp.
    #[test]
    fn gvcf2bed_padding_emits_symbolic_only_record_with_raw_ref_span() {
        let truth = vec![Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15847471, // 1-based — pos_0b = 15847470
                ref_allele: "N".to_string(),
                alt_allele: "<DEL>".to_string(),
            },
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            gt: "0/1".to_string(),
        }];
        let padding = gvcf2bed_padding(&truth, None);
        assert_eq!(padding.len(), 1, "symbolic-only record must still emit");
        // Raw refrange [pos_0b, pos_0b + reflen - 1] = [15847470, 15847470].
        // Half-open BED: [15847470, 15847471) → 1 bp.
        assert_eq!(padding[0].start, 15847470);
        assert_eq!(padding[0].end, 15847471);
    }

    #[test]
    fn gvcf2bed_padding_preserves_preprocessed_truth_spans() {
        let truth = [
            (10, "C", "A"),
            (20, "T", "A"),
            (29, "ACGTACG", "A"),
            (40, "T", "."),
            (50, "CG", "C"),
        ]
        .into_iter()
        .map(|(pos, reference, alternate)| Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos,
                ref_allele: reference.to_string(),
                alt_allele: alternate.to_string(),
            },
            qual: ".".to_string(),
            filter: "PASS".to_string(),
            gt: String::new(),
        })
        .collect::<Vec<_>>();
        let confidence = [Interval {
            chrom: "chr1".to_string(),
            start: 0,
            end: 120,
        }];
        let padding = gvcf2bed_padding(&truth, Some(&confidence));
        let size = padding
            .iter()
            .map(|interval| interval.end - interval.start)
            .sum::<usize>();

        assert_eq!(size, 10);
    }

    // Class 2 pin: `canonical_hetalt_gt` renders a query hetalt GT in
    // "alpha-later / alpha-earlier" order against the output (truth's)
    // alt list.
    //
    //   * Same-alt canonical case (C→C,G 1/2): alpha-later = G at pos 2
    //     → output `2/1`.
    //   * Same-alt NON-canonical case (T→TAA,TA 1/2): alpha-later =
    //     TAA at pos 1 → output `1/2` (verbatim).
    //   * Reordered case (output `A→ATT,AT` + query `A→AT,ATT 1/2`):
    //     alpha-later = ATT at output pos 1 → output `1/2`.
    //   * Reordered case (output `C→CA,CAA` + query `C→CAA,CA 1/2`):
    //     alpha-later = CAA at output pos 2 → output `2/1`.
    #[test]
    fn canonical_hetalt_gt_snp_canonical_swap() {
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 18280183,
                ref_allele: "T".to_string(),
                alt_allele: "C,G".to_string(),
            },
            qual: "1367.41".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("C,G", &query), "2/1");
    }

    #[test]
    fn canonical_hetalt_gt_indel_noncanonical_verbatim() {
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 21189041,
                ref_allele: "T".to_string(),
                alt_allele: "TAA,TA".to_string(),
            },
            qual: "1067.79".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("TAA,TA", &query), "1/2");
    }

    #[test]
    fn canonical_hetalt_gt_reordered_same_set() {
        // Truth A→ATT,AT vs query A→AT,ATT 1/2 → output 1/2 (ATT later).
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 15712678,
                ref_allele: "A".to_string(),
                alt_allele: "AT,ATT".to_string(),
            },
            qual: "1751.52".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("ATT,AT", &query), "1/2");
    }

    #[test]
    fn canonical_hetalt_gt_reordered_set_swap() {
        // Truth C→CA,CAA vs query C→CAA,CA 1/2 → output 2/1.
        let query = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 16032497,
                ref_allele: "C".to_string(),
                alt_allele: "CAA,CA".to_string(),
            },
            qual: "759.81".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        assert_eq!(canonical_hetalt_gt("CA,CAA", &query), "2/1");
    }

    // Class 1 pin (shared_qq picker): truth-only TP rows on a hap-matched
    // cluster propagate the minimum query QUAL, including zero, across the
    // superlocus. HG001 Platinum Genomes chr1:1876492 has QUAL=0 query
    // alleles; legacy therefore writes QQ=0 on the truth-only TP rows.
    #[test]
    fn shared_qq_picks_zero_query_qual() {
        let queries = [
            variant(15246143, "G", "C", "0/1").with_qual("817.09"),
            variant(15246157, "T", "TA", "0/1").with_qual("174.59"),
            variant(15246160, "A", "AT", "0/1").with_qual("0"),
        ];
        let min_qq: Option<&str> = queries
            .iter()
            .filter_map(|q| {
                q.qual
                    .parse::<f64>()
                    .ok()
                    .filter(|v| v.is_finite() && *v >= 0.0)
                    .map(|v| (v, q.qual.as_str()))
            })
            .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(_, s)| s);
        assert_eq!(min_qq, Some("0"));
    }

    // Class 5 pin: `cluster_signature` returns `Ok(None)` when two homalt
    // deletions on the same side have overlapping ref spans. Var1 (pos=5,
    // GAT→G) claims hap1 and hap2 ref bases 6-7 after placement; Var2
    // (pos=6, AT→A) would start at ref position 6, but h1_end=7 ≥ 6 and
    // h2_end=7 ≥ 6, so the overlap gate in `enumerate_haplotype_assignments`
    // prunes every state for Var2. The resulting empty state list leaves
    // `signatures` empty; the non-empty `variants` slice triggers the
    // `Ok(None)` return so the caller falls back to mismatch.
    #[test]
    fn cluster_signature_overlapping_homalt_deletions_returns_none() {
        let mut reference = vec![b'N'; 10];
        reference[4] = b'G'; // pos 5 (1-based)
        reference[5] = b'A'; // pos 6
        reference[6] = b'T'; // pos 7
        let reference = String::from_utf8(reference).unwrap();
        let var1 = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 5,
                ref_allele: "GAT".to_string(),
                alt_allele: "G".to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let var2 = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 6,
                ref_allele: "AT".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "30".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr1".to_string(),
            start: 5,
            end: 7,
            truth: vec![var1.clone(), var2.clone()],
            query: vec![],
        };
        let variants = vec![var1, var2];
        let result =
            cluster_signature(&cluster, &variants, &reference, None, &BTreeSet::new()).unwrap();
        assert!(
            result.is_none(),
            "overlapping homalt deletions must produce Ok(None)"
        );
    }

    /// Pin Class G: when `cluster_query_filter` aggregates filter tokens
    /// across multiple query records (truth-side TP-row stamping path),
    /// the joined string must be byte-wise sorted to match legacy's
    /// bcftools-merged ordering. chr21:40875336 cluster sources are:
    ///   * pos 40875343 T→A: `TruthSensitivityTranche99.90to100.00;LowGQX`
    ///   * pos 40875344 T→A: `TruthSensitivityTranche99.00to99.90`
    ///   * pos 40875347 A→G: `TruthSensitivityTranche99.90to100.00`
    ///
    /// Source-order union is `T99.90to100.00;LowGQX;T99.00to99.90` (rust
    /// pre-fix). Legacy emits `LowGQX;T99.00to99.90;T99.90to100.00` —
    /// alphabetic.
    #[test]
    fn cluster_query_filter_sorts_aggregated_tokens() {
        let mk = |pos: usize, filter: &str| Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos,
                ref_allele: "T".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: filter.to_string(),
            gt: "0/1".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 40875336,
            end: 40875347,
            truth: vec![],
            query: vec![
                mk(40875343, "TruthSensitivityTranche99.90to100.00;LowGQX"),
                mk(40875344, "TruthSensitivityTranche99.00to99.90"),
                mk(40875347, "TruthSensitivityTranche99.90to100.00"),
            ],
        };
        let got = cluster_query_filter(&cluster);
        assert_eq!(
            got,
            "LowGQX;TruthSensitivityTranche99.00to99.90;TruthSensitivityTranche99.90to100.00"
        );
    }

    /// Empty cluster (no PASS-bearing query) must yield ".".
    #[test]
    fn cluster_query_filter_empty_returns_dot() {
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 0,
            end: 0,
            truth: vec![],
            query: vec![],
        };
        assert_eq!(cluster_query_filter(&cluster), ".");
    }

    /// Pin Class E: chr21:44049606 cluster (chr21_passonly shape) — query
    /// has a homalt deletion AATGATAGATAG→A at 44049606 covering positions
    /// 44049607..44049617, plus a 1/2 multi-allelic at 44049615 whose
    /// second alt registers as an "insert" in `query_insert_conflict_…`.
    /// The deletion blocks the insert anchor on both haplotypes, draining
    /// query enumeration. Truth's only record (TGATA→T at 44049663) is
    /// far outside the deletion's claimed range, so legacy treats this as
    /// an internal query conflict — `BK=.`. Pre-fix the loose
    /// `!truth_remaining.is_empty()` disjunct made this fire `Some(false)`
    /// (→ BK=lm); the tightened gate must return `None`.
    #[test]
    fn deletion_covers_insert_no_proximate_truth_returns_none() {
        let q_del = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 44049606,
                ref_allele: "AATGATAGATAG".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let q_multi = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 44049615,
                ref_allele: "TAGATGATAGAT".to_string(),
                alt_allele: "T,TAGACAGATGATAGAT".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let truth_far = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 44049663,
                ref_allele: "TGATA".to_string(),
                alt_allele: "T".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "0|1".to_string(),
        };
        let query = vec![q_del, q_multi];
        let truth = vec![truth_far.clone()];
        // truth_remaining still contains the unmatched 44049663 record.
        let result = query_insert_conflict_has_truth_counterpart(&query, &truth, &truth);
        assert_eq!(
            result, None,
            "deletion-covers-insert with truth outside the deletion range \
             must return None so hap_mismatch stays false"
        );
    }

    /// Positive pin for Class E gate: when truth has a variant at the
    /// blocked insert anchor, `Some(false)` must still fire so the BK=lm
    /// branch keeps working for genuinely truth-anchored mismatches.
    #[test]
    fn deletion_covers_insert_truth_at_anchor_returns_some_false() {
        let q_del = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 100,
                ref_allele: "ATGATGATGAT".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        // Multi-allelic with a 16-base alt → registers as `insert` at pos 105.
        let q_multi = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 105,
                ref_allele: "GATGAT".to_string(),
                alt_allele: "G,GATCATGAT".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let truth_at_anchor = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 105,
                ref_allele: "G".to_string(),
                alt_allele: "GATC".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "0|1".to_string(),
        };
        let query = vec![q_del, q_multi];
        let truth = vec![truth_at_anchor];
        let result = query_insert_conflict_has_truth_counterpart(&query, &truth, &truth);
        assert_eq!(
            result,
            Some(false),
            "truth at the blocked anchor must keep the BK=lm path firing"
        );
    }

    /// Positive pin for Class E gate: when truth_remaining contains a
    /// variant that overlaps the blocking deletion's claimed range, the
    /// drain represents a genuine mismatch — `Some(false)` must fire.
    #[test]
    fn deletion_covers_insert_truth_in_del_range_returns_some_false() {
        let q_del = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 100,
                ref_allele: "ATGATGATGAT".to_string(),
                alt_allele: "A".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/1".to_string(),
        };
        let q_multi = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 105,
                ref_allele: "GATGAT".to_string(),
                alt_allele: "G,GATCATGAT".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        // Truth variant at pos 107 — inside the deletion's [101, 110] range.
        let truth_in_range = Variant {
            key: VariantKey {
                chrom: "chr1".to_string(),
                pos: 107,
                ref_allele: "T".to_string(),
                alt_allele: "C".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "0|1".to_string(),
        };
        let query = vec![q_del, q_multi];
        let truth = vec![truth_in_range];
        let result = query_insert_conflict_has_truth_counterpart(&query, &truth, &truth);
        assert_eq!(
            result,
            Some(false),
            "unmatched truth inside the blocking deletion range must keep \
             BK=lm firing"
        );
    }

    /// Class B (chr21:21690513). Reproduces the chr21 reference window
    /// `ttatatatatatatatatatatacacacacacacacatacatacatacata` at synthetic
    /// 1-based positions 1..51 (pos 1 ↔ chr21:21690480). Anchor at pos 34
    /// corresponds to chr21:21690513 (`C`). The CA-microsat upstream lets
    /// `partial_credit::left_shift` canonicalize the CACAC primitive at
    /// pos 22 (`T→TACAC`) — distinct from CACAT's stayed-put pos 34.
    /// Truth declares `C→CACAT` at pos 34 only, so the shifted CACAC has
    /// no truth representation → fan out.
    #[test]
    fn class_b_same_anchor_insertion_fans_out_when_truth_at_original_only() {
        let reference = b"ttatatatatatatatatatatacacacacacacacatacatacatacata";
        let trimmed = vec![
            (34, "C".to_string(), "CACAC".to_string()),
            (34, "C".to_string(), "CACAT".to_string()),
        ];
        let cluster_truth = vec![variant(34, "C", "CACAT", "0|1")];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            34,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &[],
        );
        let shifted = result.expect("must fan out — truth at original, orphan at shifted");
        assert_eq!(shifted.len(), 2);
        let cacac = shifted
            .iter()
            .find(|(_, _, alt)| alt.ends_with('C'))
            .expect("CACAC primitive present");
        let cacat = shifted
            .iter()
            .find(|(_, _, alt)| alt.ends_with('T'))
            .expect("CACAT primitive present");
        assert_eq!(
            cacac.0, 22,
            "CACAC must slide through CA-microsat to pos 22"
        );
        assert_eq!(cacac.1, "T", "ref byte at the shifted anchor (pos 22) is T");
        assert_eq!(cacac.2, "TACAC", "alt rotates to T-prefixed canonical form");
        assert_eq!(
            cacat.0, 34,
            "CACAT must stay at the original anchor (pos 34)"
        );
        assert_eq!(cacat.1, "C");
        assert_eq!(cacat.2, "CACAT");
    }

    /// Class B negative — chr21:40096658 shape. The TAGATAGAT primitive
    /// canonicalizes via `left_shift` to a position where truth ALREADY
    /// declares that allele as part of a multi-allelic. The fan-out
    /// gate must suppress the split so the block-level haplotype matcher
    /// can reconcile the multi-allelic record. Reproduces the real chr21
    /// AGAT-microsat upstream of pos 40096658.
    #[test]
    fn class_b_same_anchor_insertion_blocked_when_truth_at_shifted_anchor() {
        // Synthetic positions 1..50: pos 1 ↔ chr21:40096640. Pos 19 ↔
        // chr21:40096658 (anchor `T`). Pos 11 ↔ chr21:40096650 (anchor
        // `C`) where TAGATAGAT shifts to `CAGATAGAT`.
        let reference = b"ctgaagagttcagatagatagatagatagatagatagatagatagacaga";
        let trimmed = vec![
            (19, "T".to_string(), "TAGATAGAG".to_string()),
            (19, "T".to_string(), "TAGATAGAT".to_string()),
        ];
        // Truth declares TAGATAGAG at the original anchor AND
        // CAGATAGAT,CAGATAGATAGAT at the shifted anchor.
        let cluster_truth = vec![
            variant(11, "C", "CAGATAGAT,CAGATAGATAGAT", "0|1"),
            variant(19, "T", "TAGATAGAG", "1|0"),
        ];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            19,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &[],
        );
        assert!(
            result.is_none(),
            "fan-out must be blocked when truth declares the shifted alt"
        );
    }

    /// Class B negative — chr21:32767041 shape. No truth records exist in
    /// the cluster. Even though the primitives shift apart (TCTCTCT slides
    /// to a CTCT-microsat anchor, TCTCACA stays put), the fan-out must be
    /// suppressed because no query alt has truth representation at the
    /// original anchor — there's no truth_subset_match-style emit shape
    /// to license the split.
    #[test]
    fn class_b_same_anchor_insertion_blocked_without_truth_at_original() {
        let reference = b"ttatatatatatatatatatatacacacacacacacatacatacatacata";
        let trimmed = vec![
            (34, "C".to_string(), "CACAC".to_string()),
            (34, "C".to_string(), "CACAT".to_string()),
        ];
        let cluster_truth: Vec<Variant> = vec![];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            34,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &[],
        );
        assert!(
            result.is_none(),
            "fan-out must be blocked when no query alt has truth at the original anchor"
        );
    }

    /// Class B (chr21:21690513 chr21 case). When a neighboring cluster
    /// query record sits at the natural slide target, the slide must
    /// stop one position above so the shifted primitive doesn't share
    /// an anchor with the existing record. Pos 22 holds a SNP neighbor;
    /// the CACAC primitive must canonicalize at pos 23 (`A→ACACA`) — the
    /// same legacy verdict reproduced in the chr21 (no --pass-only) case.
    #[test]
    fn class_b_same_anchor_neighbor_floor_clamps_slide_target() {
        let reference = b"ttatatatatatatatatatatacacacacacacacatacatacatacata";
        let trimmed = vec![
            (34, "C".to_string(), "CACAC".to_string()),
            (34, "C".to_string(), "CACAT".to_string()),
        ];
        let cluster_truth = vec![variant(34, "C", "CACAT", "0|1")];
        // Neighboring SNP at pos 22 (`T→C`) — sliding CACAC onto pos 22
        // would clobber that record's anchor.
        let cluster_neighbors = vec![variant(22, "T", "C", "0/1")];
        let result = try_split_same_anchor_via_shift(
            "chr21",
            34,
            &trimmed,
            reference,
            1,
            &cluster_truth,
            &cluster_neighbors,
        );
        let shifted = result.expect("must still fan out — neighbor only clamps slide depth");
        // The CACAC primitive rotates to (23, A, ACACA) — the slide
        // reduces start by 1 per iteration, alternating the alt's last
        // base. With pos_min clamped to 22 (after pure-insertion bump
        // becomes 23), the slide stops with start = 23 holding 'A' anchor.
        let stayer = shifted
            .iter()
            .find(|(p, _, _)| *p == 34)
            .expect("CACAT primitive must stay at pos 34");
        let shifter = shifted
            .iter()
            .find(|(p, _, _)| *p != 34)
            .expect("shifted primitive must land at a distinct anchor");
        assert_eq!(stayer.1, "C");
        assert_eq!(stayer.2, "CACAT");
        assert_eq!(
            shifter.0, 23,
            "CACAC slide must stop at pos 23 (one above neighbor at pos 22)"
        );
        assert_eq!(shifter.1, "A", "ref byte at pos 23 is A");
        assert_eq!(
            shifter.2, "ACACA",
            "alt rotates to A-prefixed canonical form"
        );
    }

    /// Class F (chr21:47906004). A multi-allelic deletion's parent record
    /// has an `effective_refrange` that reaches into a CONF interval, but
    /// neither fanned-out primitive's range (after the per-primitive
    /// left-shift) touches CONF. Legacy emits no Regions tag on the
    /// per-primitive rows because it operates on post-fan-out records
    /// only — the parent-path's any_conf vote must be suppressed for
    /// fanned-out multi-allelics. Without this gate the cluster picks
    /// up a spurious TS_boundary tag from the parent.
    ///
    /// Synthetic layout (mirrors chr21:47906xxx at smaller positions):
    /// reference `aaaaaaaaaaaaaaaaaaaaagaactaaagt` covers 1-based positions
    /// 1..31. The variant `AGAACTAAA→A,AAAA` at pos 21 produces
    /// per-primitive rows at (17, AAAAGAACT, A) and (21, AGAACT, A) after
    /// the deletion-only slide (range 18..25 and 22..26). The parent's
    /// effective_refrange is 21..28 (alt A reaches further right than
    /// either fanned-out primitive does).
    #[test]
    fn class_f_region_state_skips_parent_path_for_fanned_out_multiallelic() {
        let reference = "aaaaaaaaaaaaaaaaaaaaagaactaaagt".to_string();
        let parent = Variant {
            key: VariantKey {
                chrom: "chr21".to_string(),
                pos: 21,
                ref_allele: "AGAACTAAA".to_string(),
                alt_allele: "A,AAAA".to_string(),
            },
            qual: "0".to_string(),
            filter: "PASS".to_string(),
            gt: "1/2".to_string(),
        };
        let cluster = Cluster {
            chrom: "chr21".to_string(),
            start: 21,
            end: 30,
            truth: vec![],
            query: vec![parent],
        };
        // CONF covers 1-based positions 27..30 — the parent's effective
        // range reaches into 27..28 but the fanned-out primitives' ranges
        // (post-slide 18..25 and 22..26) both stop at or before pos 26.
        let intervals = vec![Interval {
            chrom: "chr21".to_string(),
            start: 26,
            end: 30,
        }];
        let state = RegionState::from_cluster(&cluster, &reference, Some(&intervals));
        assert!(
            !state.any_conf,
            "fanned-out multi-allelic with all primitives outside CONF must NOT \
             register any_conf via the parent path (got any_conf=true)"
        );
        assert!(
            state.any_nonconf,
            "primitives outside CONF must vote any_nonconf"
        );
    }
}
