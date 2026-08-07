//! Command-level regression tests.

#[cfg(test)]
mod tests {
    use super::super::*;
    use flate2::read::GzDecoder;
    use std::io::Read;
    use std::sync::atomic::{AtomicU64, Ordering};

    static TEST_ID: AtomicU64 = AtomicU64::new(0);

    fn fixture(file: &str) -> PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/quantify-simple")
            .join(file)
    }

    fn test_root(label: &str) -> PathBuf {
        let id = TEST_ID.fetch_add(1, Ordering::Relaxed);
        let root =
            std::env::temp_dir().join(format!("hap-quantify-{label}-{}-{id}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        root
    }

    fn write_indexed_vcf_text(path: &Path, text: &str) {
        let source = path.with_extension("source.vcf");
        fs::write(&source, text).unwrap();
        let (headers, records) = vcf::load_raw_vcf(&source).unwrap();
        vcf::write_raw_vcf(path, &headers, &records).unwrap();
        fs::remove_file(source).unwrap();
    }

    fn indexed_fixture(root: &Path) -> PathBuf {
        let input = root.join("input.vcf.gz");
        if !input.exists() {
            // The repository fixture intentionally preserves the historical
            // invalid String-typed, multi-character QQ case used by the parity
            // matrix. Unit tests exercising successful quantify behavior need
            // the valid Float declaration the legacy numeric reader accepts.
            let (mut headers, records) = vcf::load_raw_vcf(&fixture("annotated.vcf")).unwrap();
            for header in &mut headers {
                if header.starts_with("##FORMAT=<ID=QQ,") {
                    *header = header.replace("Type=String", "Type=Float");
                }
            }
            vcf::write_raw_vcf(&input, &headers, &records).unwrap();
        }
        input
    }

    fn args(root: &Path) -> QuantifyArgs {
        QuantifyArgs {
            input_vcf: indexed_fixture(root).display().to_string(),
            report_prefix: root.join("result").display().to_string(),
            reference: fixture("ref.fa").display().to_string(),
            annotation_type: Some("ga4gh".to_string()),
            fp_bedfile: None,
            strat_tsv: None,
            strat_regions: Vec::new(),
            strat_fixchr: false,
            write_vcf: false,
            write_counts: true,
            output_vtc: false,
            preserve_info: false,
            adjust_conf_regions: None,
            threads: None,
            bcf: false,
            logfile: None,
            verbose: false,
            quiet: false,
            force_interactive: false,
            roc: "QUAL".to_string(),
            do_roc: false,
            roc_regions: vec!["*".to_string()],
            roc_filter: None,
            roc_delta: 0.5,
            ci_alpha: 0.0,
            no_json: true,
        }
    }

    fn raw_record(pos: usize, bs: usize, regions: &str, query: &str) -> RawVcfRecord {
        RawVcfRecord::from_line(
            &format!(
                "chr1\t{pos}\t.\tA\tG\t50\tPASS\tBS={bs};Regions={regions}\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:.:.:ti:SNP:het:50\t{query}"
            ),
            Path::new("test.vcf"),
        )
        .unwrap()
    }

    #[test]
    fn superlocus_region_flags_follow_recomputed_confidence_across_the_whole_block() {
        let mut records = vec![
            raw_record(2, 7, "CONF", "0/1:UNK:.:ti:SNP:het:50"),
            raw_record(3, 7, "TS_contained", "0/1:UNK:.:ti:SNP:het:50"),
            raw_record(4, 8, "CONF", "0/1:.:.:ti:SNP:het:50"),
        ];

        propagate_superlocus_annotations(&mut records, "xcmp");

        assert!(has_region(&records[0].info, "TS_boundary"));
        assert!(has_region(&records[1].info, "TS_boundary"));
        assert!(has_region(&records[1].info, "TS_contained"));
        assert!(has_region(&records[2].info, "TS_contained"));
        assert_eq!(
            records[0].sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
    }

    #[test]
    fn boundary_propagation_preserves_existing_contained_membership() {
        let mut records = vec![
            raw_record(2, 7, "CONF,TS_contained", "0/1:TP:gm:ti:SNP:het:50"),
            raw_record(3, 7, "TS_boundary", "0/1:UNK:.:ti:SNP:het:50"),
        ];

        propagate_superlocus_annotations(&mut records, "ga4gh");

        assert!(has_region(&records[0].info, "TS_contained"));
        assert!(has_region(&records[0].info, "TS_boundary"));
    }

    #[test]
    fn ga4gh_reannotation_adds_legacy_fields_from_selected_genotype_alleles() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG,AT\t50\tPASS\tBS=3;Regions=CONF\tGT:BD:BK:QQ\t0/2:TP:gm:40\t1/2:TP:am:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();

        reannotate_ga4gh_record(&mut record);

        assert_eq!(record.format.as_deref(), Some("GT:BD:BK:QQ:BI:BVT:BLT"));
        assert_eq!(
            record.samples[0], "0/2:TP:gm:40:i1_5:INDEL:het",
            "the unused SNP allele must not affect truth BVT/BI"
        );
        assert_eq!(record.samples[1], "1/2:TP:am:35:i1_5,ti:INDEL:hetalt");
        assert_eq!(
            record.sample_map(1).get("BK").map(String::as_str),
            Some("am"),
            "GA4GH BK is supplied by RTG and must be preserved"
        );
    }

    #[test]
    fn ga4gh_reannotation_classifies_homref_nocall_and_halfcall() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3\tGT:BD:BK:QQ\t0/0:TP:gm:40\t./.:N:.:.",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        reannotate_ga4gh_record(&mut record);
        assert_eq!(
            record.sample_map(0).get("BVT").map(String::as_str),
            Some("HOMREF")
        );
        assert_eq!(
            record.sample_map(0).get("BLT").map(String::as_str),
            Some("homref")
        );
        assert_eq!(
            record.sample_map(1).get("BVT").map(String::as_str),
            Some("NOCALL")
        );
        assert_eq!(
            record.sample_map(1).get("BLT").map(String::as_str),
            Some("nocall")
        );

        let halfcall = ga4gh_annotation(&record, "./1");
        assert_eq!(halfcall.bvt, "SNP");
        assert_eq!(halfcall.blt, "halfcall");
        assert_eq!(halfcall.bi, "ti");
        assert_eq!(ga4gh_annotation(&record, "1").blt, "halfcall");
    }

    #[test]
    fn ga4gh_reannotation_marks_overwide_genotype_locations_ambiguous() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t.\tPASS\t.\tGT:BD\t0/1/1:N\t.",
            Path::new("rtg.vcf"),
        )
        .unwrap();

        reannotate_ga4gh_record(&mut record);

        assert_eq!(record.format.as_deref(), Some("GT:BD:BI:BVT:BLT:QQ"));
        assert_eq!(record.samples[0], "0/1/1:N:.:UNK:ambi:.");
        assert_eq!(record.samples[1], ".:.:.:UNK:ambi:0");
    }

    #[test]
    fn ga4gh_count_unk_rewrites_both_samples_outside_confidence() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3\tGT:BD:BK:QQ\t0/1:N:.:40\t0/1:FP:.:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        annotate_regions(&mut record, Some(&[]), &RegionMap::new(), true);
        assert_eq!(
            record.sample_map(0).get("BD").map(String::as_str),
            Some("UNK")
        );
        assert_eq!(
            record.sample_map(1).get("BD").map(String::as_str),
            Some("UNK")
        );

        let mut query_only = RawVcfRecord::from_line(
            "chr1\t4\t.\tA\tG\t50\tPASS\tBS=4\tGT:BD:BK:QQ\t./.:.:.:.\t0/1:FP:.:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        annotate_regions(&mut query_only, Some(&[]), &RegionMap::new(), true);
        assert_eq!(
            query_only.sample_map(0).get("BD").map(String::as_str),
            Some("UNK")
        );
        assert_eq!(
            query_only.sample_map(1).get("BD").map(String::as_str),
            Some("UNK")
        );

        let mut filtered_truth_handoff = RawVcfRecord::from_line(
            "chr1\t5\t.\tA\tG\t50\tPASS\tBS=5\tGT:BD:BK:BI:BVT:BLT:QQ\t./.:.:.:.:NOCALL:nocall:.\t0/1:FP:.:ti:SNP:het:35",
            Path::new("xcmp.vcf"),
        )
        .unwrap();
        annotate_regions_for_samples(
            &mut filtered_truth_handoff,
            Some(&[]),
            &RegionMap::new(),
            true,
            BenchmarkSamples::POSITIONAL,
            true,
        );
        assert_eq!(
            filtered_truth_handoff
                .sample_map(0)
                .get("BD")
                .map(String::as_str),
            Some("."),
            "filtered-truth xcmp handoff preserves a missing truth decision"
        );
        assert_eq!(
            filtered_truth_handoff
                .sample_map(1)
                .get("BD")
                .map(String::as_str),
            Some("UNK")
        );

        let mut filtered_truth_only_handoff = RawVcfRecord::from_line(
            "chr1\t6\t.\tA\tG\t50\tPASS\tBS=6\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t./.:.:.:.:NOCALL:nocall:.",
            Path::new("xcmp.vcf"),
        )
        .unwrap();
        annotate_regions_for_samples(
            &mut filtered_truth_only_handoff,
            Some(&[]),
            &RegionMap::new(),
            true,
            BenchmarkSamples::POSITIONAL,
            true,
        );
        assert_eq!(
            filtered_truth_only_handoff
                .sample_map(0)
                .get("BD")
                .map(String::as_str),
            Some("UNK")
        );
        assert_eq!(
            filtered_truth_only_handoff
                .sample_map(1)
                .get("BD")
                .map(String::as_str),
            Some("."),
            "filtered-truth xcmp handoff preserves either NOCALL side"
        );
    }

    #[test]
    fn ga4gh_fp_match_kinds_feed_count_and_roc_classes() {
        let record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t0/1:FP:am:ti:SNP:het:35",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        let classified = classify_side(&record, 1).unwrap();
        let mut counts = QuantifyCountMaps::default();
        record_query(&mut counts, &classified);

        assert_eq!(counts.by_type["SNP"].fp_gt, 1);
        assert_eq!(counts.by_type["SNP"].fp_al, 0);
        assert_eq!(query_fp_class(&record), Some("gt"));

        let allele_mismatch = RawVcfRecord::from_line(
            "chr1\t4\t.\tA\tG\t50\tPASS\tBS=4\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t0/1:FP:lm:ti:SNP:het:30",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        assert_eq!(query_fp_class(&allele_mismatch), Some("al"));
    }

    #[test]
    fn xcmp_vtc_and_preserve_info_follow_legacy_cleaning_controls() {
        let source = "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3;Regions=CONF;type=TP;kind=match;ctype=simple:match;gtt1=gt_het;gtt2=gt_het;CUSTOM=kept\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:50\t0/1:TP:gm:ti:SNP:het:50";
        let mut cleaned = RawVcfRecord::from_line(source, Path::new("xcmp.vcf")).unwrap();
        decorate_quantified_record(&mut cleaned, "xcmp", false, true, true);
        assert_eq!(
            cleaned.info,
            "BS=3;Regions=CONF;XCMP=TP:match:gt_het:gt_het:simple:match;VTC=nuc__s,al__s,het__rs"
        );

        let mut preserved = RawVcfRecord::from_line(source, Path::new("xcmp.vcf")).unwrap();
        decorate_quantified_record(&mut preserved, "xcmp", true, false, true);
        assert!(preserved.info.contains("CUSTOM=kept"));
        assert!(preserved.info.contains("type=TP"));
    }

    #[test]
    fn adjusted_confidence_padding_uses_truth_records_inside_raw_confidence() {
        let truth = vec![
            RawVcfRecord::from_line(
                "chr1\t3\t.\tA\tAT\t50\tPASS\t.\tGT\t0/1",
                Path::new("truth.vcf"),
            )
            .unwrap(),
        ];
        let confidence = [Interval {
            chrom: "chr1".to_string(),
            start: 2,
            end: 3,
        }];

        let padding = truth_confidence_padding(&truth, &confidence);

        assert_eq!(padding.len(), 1);
        assert_eq!(padding[0].chrom, "chr1");
        assert_eq!((padding[0].start, padding[0].end), (2, 4));
    }

    #[test]
    fn confidence_region_precedes_existing_containment_tag() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t50\tPASS\tBS=3;Regions=TS_contained\tGT:BD\t0/1:TP\t0/1:TP",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        let confidence = [Interval {
            chrom: "chr1".to_string(),
            start: 2,
            end: 4,
        }];

        annotate_regions(&mut record, Some(&confidence), &RegionMap::new(), true);

        assert!(record.info.contains("Regions=CONF,TS_contained"));
    }

    #[test]
    fn region_updates_stay_before_preserved_extent_metadata() {
        let mut info = "BS=5;RegionsExtent=5-5;Regions=EXTRA".to_string();

        merge_region_tags(&mut info, &["CONF".to_string()]);

        assert_eq!(info, "BS=5;Regions=EXTRA,CONF;RegionsExtent=5-5");
    }

    #[test]
    fn partial_custom_stratification_overlap_marks_the_superlocus_boundary() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t8\t.\tTACG\tT\t50\tPASS\tBS=5;Regions=CONF,TS_contained\tGT:BD\t0/1:TP\t0/1:TP",
            Path::new("xcmp.vcf"),
        )
        .unwrap();
        let mut stratifications = RegionMap::new();
        stratifications.insert(
            "EXTRA".to_string(),
            vec![
                Interval {
                    chrom: "chr1".to_string(),
                    start: 0,
                    end: 7,
                },
                Interval {
                    chrom: "chr1".to_string(),
                    start: 9,
                    end: 100,
                },
            ],
        );
        let confidence = stratifications["EXTRA"].clone();

        annotate_regions(&mut record, Some(&confidence), &stratifications, false);

        assert!(has_region(&record.info, "EXTRA"));
        assert!(has_region(&record.info, "TS_boundary"));
        propagate_superlocus_annotations(std::slice::from_mut(&mut record), "ga4gh");
        assert!(has_region(&record.info, "CONF"));
        assert!(has_region(&record.info, "EXTRA"));
        assert!(has_region(&record.info, "TS_boundary"));
        assert!(has_region(&record.info, "TS_contained"));
    }

    #[test]
    fn insertion_requires_both_reference_anchors_in_confidence() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tATC\t50\tPASS\t.\tGT:BD\t0/1:TP\t0/1:TP",
            Path::new("scmp.vcf"),
        )
        .unwrap();
        let confidence = [Interval {
            chrom: "chr1".to_string(),
            start: 2,
            end: 4,
        }];

        annotate_regions(&mut record, Some(&confidence), &RegionMap::new(), true);

        assert!(has_region(&record.info, "CONF"));
        assert_eq!(
            record.sample_map(0).get("BD").map(String::as_str),
            Some("TP")
        );

        for confidence in [
            Interval {
                chrom: "chr1".to_string(),
                start: 2,
                end: 3,
            },
            Interval {
                chrom: "chr1".to_string(),
                start: 3,
                end: 4,
            },
        ] {
            let mut one_anchor = RawVcfRecord::from_line(
                "chr1\t3\t.\tA\tATC\t50\tPASS\t.\tGT:BD\t0/1:TP\t0/1:TP",
                Path::new("scmp.vcf"),
            )
            .unwrap();
            annotate_regions(
                &mut one_anchor,
                Some(std::slice::from_ref(&confidence)),
                &RegionMap::new(),
                true,
            );
            assert!(!has_region(&one_anchor.info, "CONF"));
        }
    }

    #[test]
    fn confidence_prefixed_stratification_is_folded_into_confidence() {
        let root = test_root("conf-vars");
        let base = root.join("base.bed");
        let vars = root.join("vars.bed");
        fs::write(&base, "chr1\t0\t2\n").unwrap();
        fs::write(&vars, "chr1\t2\t4\n").unwrap();
        let mut options = args(&root);
        options.fp_bedfile = Some(base.display().to_string());
        options.strat_regions = vec![format!("CONF_VARS:{}", vars.display())];
        let contigs = ["chr1".to_string()].into_iter().collect();

        let (confidence, regions, _) = load_regions(&options, &contigs).unwrap();

        let confidence = confidence.expect("combined CONF lane");
        assert_eq!(confidence.len(), 2);
        assert_eq!(region_size(&confidence), 4);
        assert!(!regions.contains_key("CONF_VARS"));
    }

    #[test]
    fn pass_truth_false_negatives_include_filtered_query_matches() {
        let mut counts = QuantifyCountMaps::default();
        counts.by_type.insert(
            "SNP".to_string(),
            QuantifyTypeCounts {
                counts: TypeCounts {
                    truth_total: CountsBucket {
                        total: 10,
                        ti: 7,
                        tv: 3,
                        ..CountsBucket::default()
                    },
                    truth_tp: CountsBucket {
                        total: 6,
                        ti: 4,
                        tv: 2,
                        ..CountsBucket::default()
                    },
                    // Explicit FN classification alone misses truth TPs whose
                    // paired query failed FILTER.
                    truth_fn: CountsBucket {
                        total: 2,
                        ti: 1,
                        tv: 1,
                        ..CountsBucket::default()
                    },
                    ..TypeCounts::default()
                },
                ..QuantifyTypeCounts::default()
            },
        );

        derive_pass_truth_false_negatives(&mut counts);

        let stats = &counts.by_type["SNP"].truth_fn;
        assert_eq!(stats.total, 4);
        assert_eq!(stats.ti, 3);
        assert_eq!(stats.tv, 1);
    }

    #[test]
    fn ga4gh_superlocus_propagates_truth_quality_and_query_nocall_filters() {
        let mut records = vec![
            RawVcfRecord::from_line(
                "chr1\t2\t.\tA\tG\t50\tPASS\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:.\t0/1:TP:gm:ti:SNP:het:30",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
            RawVcfRecord::from_line(
                "chr1\t3\t.\tC\tT\t40\tPASS\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:.\t./.:N:.:.:NOCALL:nocall:.",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
            RawVcfRecord::from_line(
                "chr1\t4\t.\tG\tA\t20\tLowQual\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:FN:.:ti:SNP:het:.\t0/1:FP:.:ti:SNP:het:20",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
        ];

        propagate_superlocus_annotations(&mut records, "ga4gh");

        assert_eq!(
            records[0].sample_map(0).get("QQ").map(String::as_str),
            Some("30")
        );
        assert_eq!(
            records[1].sample_map(0).get("QQ").map(String::as_str),
            Some("30"),
            "truth TP without a paired query TP receives the minimum TP score in BS"
        );
        assert_eq!(
            records[2].sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(records[1].filter, "LowQual");
    }

    #[test]
    fn ga4gh_superlocus_does_not_propagate_nan_query_quality_to_truth() {
        let mut records = vec![
            RawVcfRecord::from_line(
                "chr1\t2\t.\tA\tG\t50\tPASS\tBS=9;Regions=CONF\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:.\t0/1:TP:gm:ti:SNP:het:nan",
                Path::new("rtg.vcf"),
            )
            .unwrap(),
        ];

        propagate_superlocus_annotations(&mut records, "ga4gh");

        assert_eq!(
            records[0].sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(
            records[0].sample_map(1).get("QQ").map(String::as_str),
            Some("nan")
        );
    }

    #[test]
    fn ga4gh_output_restores_missing_qualities_and_integer_rendering() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t30\t.\tC\t<DEL>\t.\tPASS\t.\tGT:BD\t0/1:N\t./.:N",
            Path::new("rtg.vcf"),
        )
        .unwrap();
        reannotate_ga4gh_record(&mut record);
        propagate_ga4gh_superlocus_for_samples(
            std::slice::from_mut(&mut record),
            BenchmarkSamples::POSITIONAL,
            false,
            false,
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("0")
        );
        set_format_value(&mut record, 1, "QQ", "60.0");
        normalize_integer_like_format_values(&mut record, "QQ");

        assert_eq!(record.format.as_deref(), Some("GT:BD:BI:BVT:BLT:QQ"));
        assert_eq!(
            record.sample_map(0).get("QQ").map(String::as_str),
            Some(".")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("60")
        );

        let mut somatic_record = RawVcfRecord::from_line(
            "chr1\t30\t.\tC\t<DEL>\t.\tPASS\t.\tGT:BD\t0/1:N\t./.:N",
            Path::new("scmp-somatic.vcf"),
        )
        .unwrap();
        reannotate_ga4gh_record(&mut somatic_record);
        propagate_ga4gh_superlocus_for_samples(
            std::slice::from_mut(&mut somatic_record),
            BenchmarkSamples::POSITIONAL,
            true,
            false,
        );
        assert_eq!(
            somatic_record.sample_map(1).get("QQ").map(String::as_str),
            Some("."),
            "scmp-somatic preserves a missing query score"
        );
    }

    #[test]
    fn ga4gh_output_headers_add_only_missing_legacy_declarations() {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY".to_string(),
        ];
        ensure_ga4gh_headers(&mut headers);
        ensure_ga4gh_headers(&mut headers);
        canonicalize_ga4gh_header_order(&mut headers);
        for id in ["GT", "BD", "BK", "BI", "QQ", "BVT", "BLT"] {
            assert_eq!(
                headers
                    .iter()
                    .filter(|header| header.contains(&format!("FORMAT=<ID={id},")))
                    .count(),
                1,
                "header {id} must be present exactly once"
            );
        }
        assert!(headers.iter().any(|header| header.contains("INFO=<ID=BS,")));
        assert!(headers[1].starts_with("##FILTER=<ID=PASS,"));
        assert!(headers.last().unwrap().starts_with("#CHROM"));
    }

    #[test]
    fn ga4gh_accepts_a_string_qq_declaration() {
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FORMAT=<ID=QQ,Number=1,Type=String,Description=\"Score\">".to_string(),
        ];
        let mut record = RawVcfRecord::from_line(
            "chr1\t2\t.\tA\tG\t60\tPASS\t.\tGT:QQ\t0/1:7\t0/1:.",
            Path::new("fixture.vcf"),
        )
        .unwrap();
        assert!(validate_ga4gh_qq_fields(&headers, std::slice::from_ref(&record)).is_ok());

        record.samples[0] = "0/1:60".to_string();
        let error = validate_ga4gh_qq_fields(&headers, &[record]).unwrap_err();
        assert_eq!(error.to_string(), "too many QQ fields at chr1:1");
    }

    #[test]
    fn region_header_uses_the_legacy_bcf_compatible_declaration() {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        ensure_info_header(
            &mut headers,
            "Regions",
            "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">",
        );
        assert_eq!(
            headers[1],
            "##INFO=<ID=Regions,Number=.,Type=String,Description=\"Tags for regions.\">"
        );
    }

    fn read_gzip(path: &Path) -> String {
        let mut text = String::new();
        GzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    }

    fn write_named_sample_vcf(path: &Path, sample_names: &[&str], samples: &[Vec<&str>]) {
        let mut text = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=10>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision\">\n",
            "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Kind\">\n",
            "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Info\">\n",
            "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"Type\">\n",
            "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"Location\">\n",
            "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Quality\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
        )
        .to_string();
        for name in sample_names {
            text.push('\t');
            text.push_str(name);
        }
        text.push('\n');
        for (record_index, record_samples) in samples.iter().enumerate() {
            let pos = record_index + 2;
            let qual = 30 + record_index * 10;
            text.push_str(&format!(
                "chr1\t{pos}\t.\tA\tG\t{qual}\tPASS\tBS={pos}\tGT:BD:BK:BI:BVT:BLT:QQ"
            ));
            for sample in record_samples {
                text.push('\t');
                text.push_str(sample);
            }
            text.push('\n');
        }
        write_indexed_vcf_text(path, &text);
    }

    fn summary_snp_all(path: &Path) -> Vec<String> {
        fs::read_to_string(path)
            .unwrap()
            .lines()
            .find(|line| line.starts_with("SNP,ALL,"))
            .unwrap()
            .split(',')
            .map(str::to_string)
            .collect()
    }

    #[test]
    fn compressed_quantifier_inputs_require_their_index() {
        let root = test_root("input-index");
        let (mut headers, records) = vcf::load_raw_vcf(&fixture("annotated.vcf")).unwrap();
        let chrom_header = headers
            .iter()
            .position(|header| header.starts_with("#CHROM"))
            .unwrap();
        headers.insert(
            chrom_header,
            "##INFO=<ID=BS,Number=1,Type=Integer,Description=\"Benchmark superlocus\">".to_string(),
        );
        for (name, index_suffix) in [("input.vcf.gz", ".tbi"), ("input.bcf", ".csi")] {
            let input = root.join(name);
            vcf::write_raw_vcf(&input, &headers, &records).unwrap();
            let index = PathBuf::from(format!("{}{}", input.display(), index_suffix));
            assert!(index.is_file());
            fs::remove_file(index).unwrap();

            let case_root = root.join(format!("case-{name}"));
            let mut options = args(&case_root);
            options.input_vcf = input.display().to_string();
            let error = run(options).unwrap_err();
            assert!(
                error.to_string().contains("index"),
                "unexpected missing-index error: {error:#}"
            );
        }

        let plain_root = root.join("plain");
        let plain_input = plain_root.join("input.vcf");
        fs::create_dir_all(&plain_root).unwrap();
        fs::copy(fixture("annotated.vcf"), &plain_input).unwrap();
        let mut options = args(&plain_root);
        options.input_vcf = plain_input.display().to_string();
        let error = run(options).unwrap_err();
        assert!(error.to_string().contains("compressed and indexed"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn named_truth_and_query_samples_drive_counts_and_rocs_in_any_column() {
        let root = test_root("named-samples");
        let cases = [
            (
                "reversed",
                vec!["QUERY", "TRUTH"],
                vec![
                    vec!["0/1:FP:lm:ti:SNP:het:30", "0/1:FN:lm:ti:SNP:het:30"],
                    vec!["0/1:TP:gm:ti:SNP:het:40", "0/1:TP:gm:ti:SNP:het:40"],
                ],
            ),
            (
                "extra",
                vec!["NOISE", "QUERY", "AUX", "TRUTH"],
                vec![
                    vec![
                        "0/1:TP:gm:ti:SNP:het:99",
                        "0/1:FP:lm:ti:SNP:het:30",
                        "0/1:TP:gm:ti:SNP:het:98",
                        "0/1:FN:lm:ti:SNP:het:30",
                    ],
                    vec![
                        "0/1:TP:gm:ti:SNP:het:99",
                        "0/1:TP:gm:ti:SNP:het:40",
                        "0/1:TP:gm:ti:SNP:het:98",
                        "0/1:TP:gm:ti:SNP:het:40",
                    ],
                ],
            ),
        ];

        for (label, names, samples) in cases {
            let case_root = root.join(label);
            fs::create_dir_all(&case_root).unwrap();
            let input = case_root.join("input.vcf.gz");
            write_named_sample_vcf(&input, &names, &samples);
            let mut options = args(&case_root);
            options.input_vcf = input.display().to_string();
            options.do_roc = true;
            run(options).unwrap();

            let fields = summary_snp_all(&case_root.join("result.summary.csv"));
            assert_eq!(&fields[2..10], ["2", "1", "1", "2", "1", "0", "0", "1"]);
            let roc = read_gzip(&case_root.join("result.roc.all.csv.gz"));
            assert!(roc.lines().any(|line| {
                let fields = line.split(',').collect::<Vec<_>>();
                fields.first() == Some(&"SNP")
                    && fields.get(6) == Some(&"*")
                    && fields.get(51) == Some(&"1")
            }));
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn missing_truth_or_query_names_disable_rocs_without_positional_counts() {
        let root = test_root("missing-named-samples");
        let input = root.join("named.vcf.gz");
        write_named_sample_vcf(
            &input,
            &["FIRST", "SECOND"],
            &[vec!["0/1:TP:gm:ti:SNP:het:30", "0/1:FP:lm:ti:SNP:het:30"]],
        );
        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.do_roc = true;
        run(options).unwrap();

        let summary = fs::read_to_string(root.join("result.summary.csv")).unwrap();
        assert_eq!(summary.lines().count(), 3);
        assert!(
            summary
                .lines()
                .skip(1)
                .all(|line| line.contains(",ALL,0,0,0,0,0,0,0,0,"))
        );
        assert!(!root.join("result.roc.Locations.SNP.csv.gz").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn explicit_zero_roc_delta_uses_the_low_level_default() {
        let root = test_root("zero-roc-delta");
        let input = root.join("named.vcf.gz");
        write_named_sample_vcf(
            &input,
            &["TRUTH", "QUERY"],
            &[
                vec!["0/1:TP:gm:ti:SNP:het:10", "0/1:TP:gm:ti:SNP:het:10"],
                vec!["0/1:TP:gm:ti:SNP:het:10.05", "0/1:TP:gm:ti:SNP:het:10.05"],
            ],
        );
        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.do_roc = true;
        options.roc_delta = 0.0;
        run(options).unwrap();

        let roc = read_gzip(&root.join("result.roc.Locations.SNP.csv.gz"));
        let levels = roc
            .lines()
            .skip(1)
            .filter_map(|line| line.split(',').nth(6))
            .filter(|level| *level != "*")
            .collect::<Vec<_>>();
        assert_eq!(levels, ["10.000000"]);
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn preserve_info_output_declares_and_populates_regions_extent() {
        let root = test_root("regions-extent");
        let mut options = args(&root);
        options.annotation_type = Some("xcmp".to_string());
        options.write_vcf = true;
        options.preserve_info = true;
        run(options).unwrap();

        let (headers, records) =
            vcf::load_raw_vcf(&root.join("result.vcf.gz")).expect("quantified VCF");
        assert_eq!(
            headers
                .iter()
                .filter(|line| line.contains("INFO=<ID=RegionsExtent,"))
                .count(),
            1
        );
        assert_eq!(
            records
                .iter()
                .map(|record| info_value(&record.info, "RegionsExtent").unwrap())
                .collect::<Vec<_>>(),
            ["2-2", "8-9"]
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn confidence_and_named_stratifications_drive_counts_and_sizes() {
        let root = test_root("regions");
        let conf = root.join("confidence.bed");
        let early = root.join("early.bed");
        let late = root.join("late.bed");
        let strata = root.join("regions.tsv");
        fs::write(&conf, "1\t0\t5\n").unwrap();
        fs::write(&early, "1\t0\t5\n").unwrap();
        fs::write(&late, "chr1\t5\t10\n").unwrap();
        fs::write(&strata, "EARLY\tearly.bed\n").unwrap();

        let mut options = args(&root);
        options.fp_bedfile = Some(conf.display().to_string());
        options.strat_tsv = Some(strata.display().to_string());
        options.strat_regions = vec![format!("LATE:{}", late.display())];
        options.strat_fixchr = true;
        run(options).unwrap();

        let summary = fs::read_to_string(root.join("result.summary.csv")).unwrap();
        assert!(
            summary.contains("INDEL,ALL,0,0,0,1,0,1,0,0,0.0,,1.0,"),
            "the call outside confidence must be QUERY.UNK: {summary}"
        );
        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        let early_row = extended
            .lines()
            .find(|line| line.starts_with("SNP,*,EARLY,ALL,"))
            .expect("named TSV stratification row");
        let early_fields = early_row.split(',').collect::<Vec<_>>();
        assert_eq!(early_fields[13], "5.000000");
        assert_eq!(early_fields[14], "5.000000");
        assert!(
            extended
                .lines()
                .any(|line| line.starts_with("INDEL,*,LATE,ALL,")),
            "direct NAME:BED stratification must be applied"
        );
        assert!(root.join("result.roc.all.csv.gz").exists());
        assert!(!root.join("result.roc.Locations.SNP.csv.gz").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn four_column_stratifications_add_dynamic_child_lanes_and_parent_membership() {
        let root = test_root("four-column-regions");
        let regions = root.join("hierarchy.bed");
        fs::write(
            &regions,
            "chr1\t1\t2\tcoding_1\nchr1\t4\t5\tunused\nchr1\t7\t9\tcoding_2\n",
        )
        .unwrap();

        let dynamic_root = root.join("dynamic");
        let mut dynamic = args(&dynamic_root);
        dynamic.strat_regions = vec![format!("EXTRA:{}", regions.display())];
        dynamic.write_vcf = true;
        run(dynamic).unwrap();

        let (_, records) = vcf::load_raw_vcf(&dynamic_root.join("result.vcf.gz")).unwrap();
        assert!(has_region(&records[0].info, "EXTRA"));
        assert!(has_region(&records[0].info, "EXTRA_coding_1"));
        assert!(has_region(&records[1].info, "EXTRA"));
        assert!(has_region(&records[1].info, "EXTRA_coding_2"));
        assert!(
            records
                .iter()
                .all(|record| !has_region(&record.info, "EXTRA_unused"))
        );
        assert!(
            !records
                .iter()
                .any(|record| has_region(&record.info, "EXTRA_coding"))
        );

        let extended = fs::read_to_string(dynamic_root.join("result.extended.csv")).unwrap();
        for (variant_type, subset, size, level) in [
            ("SNP", "EXTRA", "4.000000", "0.000000"),
            ("SNP", "EXTRA_coding_1", "1.000000", "1.000000"),
            ("SNP", "EXTRA_coding_2", "2.000000", "1.000000"),
            ("SNP", "EXTRA_unused", "1.000000", "1.000000"),
            ("INDEL", "EXTRA", "4.000000", "0.000000"),
            ("INDEL", "EXTRA_coding_1", "1.000000", "1.000000"),
            ("INDEL", "EXTRA_coding_2", "2.000000", "1.000000"),
            ("INDEL", "EXTRA_unused", "1.000000", "1.000000"),
        ] {
            let row = extended
                .lines()
                .find(|line| line.starts_with(&format!("{variant_type},*,{subset},ALL,")))
                .unwrap_or_else(|| panic!("missing dynamic subset row {variant_type}/{subset}"));
            let fields = row.split(',').collect::<Vec<_>>();
            assert_eq!(fields[13], size, "wrong size for {variant_type}/{subset}");
            assert_eq!(fields[15], level, "wrong level for {variant_type}/{subset}");
        }
        assert!(
            !extended
                .lines()
                .any(|line| { line.split(',').nth(2) == Some("EXTRA_coding") })
        );
        let roc = read_gzip(&dynamic_root.join("result.roc.all.csv.gz"));
        let dynamic_roc = roc
            .lines()
            .find(|line| {
                let fields = line.split(',').collect::<Vec<_>>();
                fields.first() == Some(&"SNP")
                    && fields.get(2) == Some(&"EXTRA_coding_1")
                    && fields.get(6) == Some(&"*")
            })
            .expect("dynamic child ROC row");
        assert_eq!(dynamic_roc.split(',').nth(15), Some("1.000000"));
        for variant_type in ["SNP", "INDEL"] {
            for filter in ["ALL", "PASS"] {
                let empty_roc = roc
                    .lines()
                    .find(|line| {
                        let fields = line.split(',').collect::<Vec<_>>();
                        fields.first() == Some(&variant_type)
                            && fields.get(1) == Some(&"*")
                            && fields.get(2) == Some(&"EXTRA_unused")
                            && fields.get(3) == Some(&filter)
                            && fields.get(6) == Some(&"*")
                    })
                    .unwrap_or_else(|| {
                        panic!("missing empty dynamic ROC row {variant_type}/{filter}")
                    });
                let fields = empty_roc.split(',').collect::<Vec<_>>();
                assert_eq!(fields[13], "1.000000");
                assert_eq!(fields[15], "1.000000");
                assert_eq!(fields[16], "0");
            }
        }

        let fixed_root = root.join("fixed");
        let mut fixed = args(&fixed_root);
        fixed.strat_regions = vec![format!("=EXTRA:{}", regions.display())];
        fixed.write_vcf = true;
        run(fixed).unwrap();
        let (_, fixed_records) = vcf::load_raw_vcf(&fixed_root.join("result.vcf.gz")).unwrap();
        assert!(
            fixed_records
                .iter()
                .all(|record| has_region(&record.info, "EXTRA"))
        );
        assert!(fixed_records.iter().all(|record| {
            !parse_subsets(&record.info)
                .iter()
                .any(|subset| subset.starts_with("EXTRA_"))
        }));

        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn default_qual_roc_preserves_count_tables_and_writes_roc_files() {
        let root = test_root("roc");
        let baseline_root = root.join("baseline");
        let roc_root = root.join("roc");
        let baseline = args(&baseline_root);
        run(baseline).unwrap();
        let mut options = args(&roc_root);
        options.do_roc = true;
        run(options).unwrap();

        assert_eq!(
            fs::read(baseline_root.join("result.summary.csv")).unwrap(),
            fs::read(roc_root.join("result.summary.csv")).unwrap()
        );
        assert_eq!(
            fs::read(baseline_root.join("result.extended.csv")).unwrap(),
            fs::read(roc_root.join("result.extended.csv")).unwrap()
        );
        for suffix in [
            "roc.all.csv.gz",
            "roc.Locations.SNP.csv.gz",
            "roc.Locations.SNP.PASS.csv.gz",
            "roc.Locations.INDEL.csv.gz",
            "roc.Locations.INDEL.PASS.csv.gz",
        ] {
            let path = roc_root.join(format!("result.{suffix}"));
            assert!(
                path.metadata().unwrap().len() > 0,
                "{} is empty",
                path.display()
            );
        }
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn inherited_roc_controls_flow_through_quantify_outputs() {
        let root = test_root("roc-controls");
        let input = root.join("annotated.vcf.gz");
        let source = fs::read_to_string(fixture("annotated.vcf")).unwrap();
        let source = source
            .replace(
                "##FORMAT=<ID=QQ,Number=1,Type=String",
                "##FORMAT=<ID=QQ,Number=1,Type=Float",
            )
            .replace(
                "##FORMAT=<ID=GT",
                "##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"ROC score\">\n##FORMAT=<ID=GT",
            )
            .replace("PASS\tBS=2", "LowQual\tBS=2;SCORE=10.0")
            .replace("PASS\tBS=8", "PASS\tBS=8;SCORE=10.4");
        write_indexed_vcf_text(&input, &source);

        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.roc = "SCORE".to_string();
        options.roc_filter = Some("LowQual".to_string());
        options.roc_delta = 0.0;
        options.ci_alpha = 0.05;
        options.do_roc = true;
        options.no_json = false;
        run(options).unwrap();

        let all = read_gzip(&root.join("result.roc.all.csv.gz"));
        let mut lines = all.lines();
        let header = lines.next().unwrap();
        assert!(header.ends_with("METRIC.Frac_NA.Lower,METRIC.Frac_NA.Upper"));
        assert!(lines.any(|line| {
            let fields = line.split(',').collect::<Vec<_>>();
            fields[5] == "SCORE" && fields[6] == "10.000000"
        }));
        assert!(root.join("result.roc.Locations.SNP.SEL.csv.gz").exists());
        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        let mut extended_lines = extended.lines();
        assert!(
            extended_lines
                .next()
                .unwrap()
                .ends_with("METRIC.Frac_NA.Lower,METRIC.Frac_NA.Upper")
        );
        assert!(extended_lines.all(|line| line.split(',').nth(5) == Some("SCORE")));
        let metrics = read_gzip(&root.join("result.metrics.json.gz"));
        assert!(metrics.contains("\"id\":\"roc.Locations.SNP.SEL\""));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn metrics_json_uses_legacy_table_schema_and_column_types() {
        let root = test_root("metrics-json");
        let mut options = args(&root);
        options.do_roc = true;
        options.no_json = false;
        run(options).unwrap();

        let json = read_gzip(&root.join("result.metrics.json.gz"));
        assert!(json.contains("\"name\":\"qfy.py.comparison\""));
        assert!(json.contains("\"module\":\"qfy.py\""));
        for id in [
            "summary.metrics",
            "all.metrics",
            "roc.all",
            "roc.Locations.INDEL.PASS",
            "roc.Locations.SNP.PASS",
            "roc.Locations.SNP",
            "roc.Locations.INDEL",
        ] {
            assert!(
                json.contains(&format!("\"id\":\"{id}\"")),
                "missing legacy metrics table {id}"
            );
        }
        assert!(
            json.contains("\"type\":\"int64\",\"id\":\"TRUTH.TOTAL\",\"label\":\"TRUTH.TOTAL\"")
        );
        assert!(
            json.contains(
                "\"type\":\"double\",\"id\":\"METRIC.Recall\",\"label\":\"METRIC.Recall\""
            )
        );
        assert!(json.contains("\"type\":\"string\",\"id\":\"Type\",\"label\":\"Type\""));
        assert!(!json.contains("\"tool\":\"hap quantify\""));
        assert!(!json.contains("\"truth_total\""));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn metrics_json_omits_all_table_when_extended_counts_are_disabled() {
        let root = test_root("metrics-json-no-counts");
        let mut options = args(&root);
        options.write_counts = false;
        options.no_json = false;
        run(options).unwrap();

        let json = read_gzip(&root.join("result.metrics.json.gz"));
        assert!(json.contains("\"id\":\"summary.metrics\""));
        assert!(json.contains("\"id\":\"roc.all\""));
        assert!(!json.contains("\"id\":\"all.metrics\""));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn no_write_counts_keeps_summary_and_suppresses_only_extended_table() {
        let root = test_root("no-counts");
        let mut options = args(&root);
        options.write_counts = false;
        run(options).unwrap();

        assert!(root.join("result.summary.csv").is_file());
        assert!(!root.join("result.extended.csv").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn dotted_report_prefix_is_preserved_for_every_output_family() {
        let root = test_root("dotted-prefix");
        let mut options = args(&root);
        options.report_prefix = root.join("sample.v1").display().to_string();
        options.write_vcf = true;
        run(options).unwrap();

        for suffix in [
            "summary.csv",
            "extended.csv",
            "vcf.gz",
            "vcf.gz.tbi",
            "roc.all.csv.gz",
        ] {
            assert!(
                root.join(format!("sample.v1.{suffix}")).is_file(),
                "missing dotted-prefix output {suffix}"
            );
        }
        assert!(!root.join("sample.summary.csv").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn bcf_logfile_and_verbose_controls_publish_legacy_standalone_artifacts() {
        let root = test_root("bcf-logfile");
        let mut options = args(&root);
        options.write_vcf = true;
        options.bcf = true;
        options.logfile = Some(root.join("qfy.log").display().to_string());
        options.verbose = true;
        options.do_roc = false;
        run(options).unwrap();

        assert!(root.join("result.bcf").is_file());
        assert!(root.join("result.bcf.csi").is_file());
        assert!(!root.join("result.vcf.gz").exists());
        assert!(root.join("qfy.log").is_file());
        let raw_roc = fs::read_to_string(root.join("result.roc.tsv")).unwrap();
        let mut lines = raw_roc.lines();
        let header = lines.next().unwrap();
        assert!(header.starts_with("FP.al\tFP.gt\tFilter\tGenotype\tMETRIC.F1_Score"));
        let columns = header.split('\t').collect::<Vec<_>>();
        let qq = columns.iter().position(|column| *column == "QQ").unwrap();
        assert!(lines.all(|line| line.split('\t').nth(qq) == Some("*")));
        assert!(header.contains("QUERY.TOTAL.hetalt"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn write_vcf_refuses_to_overwrite_its_input_before_writing_reports() {
        let root = test_root("input-overwrite");
        let input = root.join("result.vcf.gz");
        fs::copy(fixture("annotated.vcf.gz"), &input).unwrap();
        let before = fs::read(&input).unwrap();
        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.write_vcf = true;

        let error = run(options).unwrap_err();
        assert!(error.to_string().contains("cannot overwrite input VCF"));
        assert_eq!(fs::read(&input).unwrap(), before);
        assert!(!root.join("result.summary.csv").exists());
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn xcmp_mode_rederives_decisions_instead_of_consuming_final_bd_fields() {
        let root = test_root("xcmp-reannotation");
        let input = root.join("annotated.vcf.gz");
        let reference = root.join("ref.fa");
        let confidence = root.join("confidence.bed");
        fs::write(&reference, ">chr1\nNNACGTNN\n").unwrap();
        fs::write(&confidence, "chr1\t0\t3\n").unwrap();
        write_indexed_vcf_text(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "##FORMAT=<ID=BD,Number=1,Type=String,Description=\"Decision\">\n",
                "##FORMAT=<ID=BK,Number=1,Type=String,Description=\"Kind\">\n",
                "##FORMAT=<ID=BI,Number=1,Type=String,Description=\"Info\">\n",
                "##FORMAT=<ID=BVT,Number=1,Type=String,Description=\"Type\">\n",
                "##FORMAT=<ID=BLT,Number=1,Type=String,Description=\"Location\">\n",
                "##FORMAT=<ID=QQ,Number=1,Type=Float,Description=\"Quality\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTRUTH\tQUERY\n",
                // Final BD fields inside CONF are ignored because old XCMP
                // record-level INFO/type is absent.
                "chr1\t2\t.\tA\tG\t60\tPASS\t.\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:60\t0/1:TP:gm:ti:SNP:het:60\n",
                // Outside CONF, XCMP's count_unk rule replaces the absent
                // decision with UNK for both samples; only QUERY contributes.
                "chr1\t5\t.\tG\tA\t50\tPASS\t.\tGT:BD:BK:BI:BVT:BLT:QQ\t0/1:TP:gm:ti:SNP:het:50\t0/1:TP:gm:ti:SNP:het:50\n",
            ),
        );

        let mut options = args(&root);
        options.input_vcf = input.display().to_string();
        options.reference = reference.display().to_string();
        options.fp_bedfile = Some(confidence.display().to_string());
        options.annotation_type = Some("xcmp".to_string());
        run(options).unwrap();

        let summary = fs::read_to_string(root.join("result.summary.csv")).unwrap();
        assert!(summary.contains("SNP,ALL,0,0,0,1,0,1,0,0,0.0,,1.0,"));
        let extended = fs::read_to_string(root.join("result.extended.csv")).unwrap();
        let snp_all = extended
            .lines()
            .find(|line| line.starts_with("SNP,*,*,ALL,"))
            .unwrap()
            .split(',')
            .collect::<Vec<_>>();
        assert_eq!(
            snp_all[13], "4",
            "terminal Ns are excluded from Subset.Size"
        );
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn xcmp_mode_maps_record_level_fp_to_truth_fn_and_query_fp() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t40\tPASS\ttype=FP;kind=gtmismatch;Regions=CONF\tGT:BD:BK:QQ\t0/1:TP:gm:1\t0/1:TP:gm:2",
            Path::new("input.vcf"),
        )
        .unwrap();
        reannotate_xcmp_record(&mut record, true, "QUAL");
        assert_eq!(
            record.sample_map(0).get("BD").map(String::as_str),
            Some("FN")
        );
        assert_eq!(
            record.sample_map(1).get("BD").map(String::as_str),
            Some("FP")
        );
        assert_eq!(
            record.sample_map(0).get("BK").map(String::as_str),
            Some("am")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("40")
        );
    }

    #[test]
    fn xcmp_custom_format_roc_field_is_copied_into_qq() {
        let mut record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\tG\t40\tPASS\ttype=TP;kind=match\tGT:BD:BK:QQ:GQX\t0/1:TP:gm:1:17.5\t0/1:TP:gm:2:23.5",
            Path::new("input.vcf"),
        )
        .unwrap();
        reannotate_xcmp_record(&mut record, false, "GQX");
        assert_eq!(
            record.sample_map(0).get("QQ").map(String::as_str),
            Some("17.5")
        );
        assert_eq!(
            record.sample_map(1).get("QQ").map(String::as_str),
            Some("23.5")
        );
    }

    #[test]
    fn invalid_quantify_controls_fail_before_reading_inputs() {
        let root = test_root("invalid-controls");
        let baseline = args(&root);

        let mut unknown_type = baseline.clone();
        unknown_type.annotation_type = Some("unknown".to_string());
        assert!(
            run(unknown_type)
                .unwrap_err()
                .to_string()
                .contains("unsupported")
        );

        let mut roc_field = baseline.clone();
        roc_field.roc.clear();
        assert!(run(roc_field).unwrap_err().to_string().contains("empty"));

        let mut roc_regions = baseline.clone();
        roc_regions.roc_regions.push(String::new());
        assert!(run(roc_regions).unwrap_err().to_string().contains("empty"));

        let mut roc_delta = baseline.clone();
        roc_delta.roc_delta = -1.0;
        assert!(
            run(roc_delta)
                .unwrap_err()
                .to_string()
                .contains("nonnegative")
        );

        let mut ci = baseline;
        ci.ci_alpha = 1.5;
        assert!(run(ci).unwrap_err().to_string().contains("ci-alpha"));
        fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn reported_bed_size_uses_half_open_bed_span() {
        let intervals = [Interval {
            chrom: "chr1".to_string(),
            start: 0,
            end: 16,
        }];
        assert_eq!(region_size(&intervals), 16);
    }

    #[test]
    fn symbolic_alt_has_no_nucleotide_reference_range() {
        let record = RawVcfRecord::from_line(
            "chr1\t3\t.\tA\t<DEL>\t40\tPASS\tEND=7\tGT\t0/1",
            Path::new("symbolic.vcf"),
        )
        .unwrap();
        assert_eq!(effective_reference_range(&record), None);
    }

    #[test]
    fn reference_range_stops_at_first_symbolic_alt() {
        let parse = |alt: &str| {
            RawVcfRecord::from_line(
                &format!("chr1\t3\t.\tA\t{alt}\t40\tPASS\t.\tGT\t0/1"),
                Path::new("mixed.vcf"),
            )
            .unwrap()
        };

        // A leading symbolic allele prevents the later insertion from being
        // considered at all.
        assert_eq!(effective_reference_range(&parse("<DEL>,AT")), None);
        // A preceding nucleotide allele is retained, while nucleotide alleles
        // after the symbolic one are ignored.
        assert_eq!(
            effective_reference_range(&parse("AT,<DEL>,AG")),
            effective_reference_range(&parse("AT"))
        );
        assert_eq!(effective_reference_range(&parse("AT")), Some((3, 4, true)));
    }

    #[test]
    fn missing_alt_is_processed_as_an_empty_allele() {
        let parse = |alt: &str| {
            RawVcfRecord::from_line(
                &format!("chr1\t3\t.\tAT\t{alt}\t40\tPASS\t.\tGT\t0/1"),
                Path::new("missing.vcf"),
            )
            .unwrap()
        };

        let reference_span = Some((3, 4, false));
        assert_eq!(effective_reference_range(&parse(".")), reference_span);
        assert_eq!(effective_reference_range(&parse("")), reference_span);
        // Missing is retained before the symbolic ALT terminates processing;
        // the later nucleotide ALT cannot narrow that reference span.
        assert_eq!(
            effective_reference_range(&parse(".,<DEL>,A")),
            reference_span
        );
    }

    #[test]
    fn named_subset_report_sizes_distinguish_boundary_and_confidence_intersection() {
        let stratification_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);
        let stratification_confidence_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);

        assert_eq!(
            named_subset_report_sizes(
                "TS_boundary",
                100,
                140,
                Some(141),
                &stratification_sizes,
                &stratification_confidence_sizes,
            ),
            (140, Some(141))
        );
        assert_eq!(
            named_subset_report_sizes(
                "EXTRA",
                100,
                140,
                Some(141),
                &stratification_sizes,
                &stratification_confidence_sizes,
            ),
            (138, Some(138))
        );
    }

    #[test]
    fn confidence_intersection_uses_unique_interval_union() {
        let interval = |chrom: &str, start, end| Interval {
            chrom: chrom.to_string(),
            start,
            end,
        };
        let subset = vec![interval("chr1", 0, 98), interval("chrX", 0, 40)];
        let confidence = vec![
            interval("chr1", 0, 98),
            interval("chrX", 0, 40),
            interval("chr1", 1, 4),
        ];

        assert_eq!(region_size(&confidence), 141);
        assert_eq!(region_intersection_size(&subset, &confidence), 138);
    }
}
