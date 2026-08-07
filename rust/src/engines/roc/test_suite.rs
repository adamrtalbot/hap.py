//! Command-level regression tests.

#[cfg(test)]
mod tests {
    use super::super::*;

    fn comparison_record(line: &str) -> crate::domain::ComparisonRecord {
        crate::domain::RawVcfRecord::from_line(line, std::path::Path::new("roc-test.vcf"))
            .unwrap()
            .into()
    }

    #[test]
    fn half_call_is_not_counted_as_heterozygous_in_roc_stats() {
        let format = ["GT", "BD", "BI", "BVT", "BLT", "QQ"];
        let values = ["./1", "TP", "tv", "SNP", "halfcall", "60"];
        let sample = Sample::new(&format, &values);
        let bucket = sample_bucket(&sample);
        assert_eq!(bucket.het, 0);
        assert_eq!(bucket.homalt, 0);
        assert_eq!(bucket.tv, 1);
    }

    #[test]
    fn legacy_string_hash_matches_pinned_libstdcpp() {
        let key = "SNP\t*\t*\tPASS\tTS_contained\t379.290009";
        assert_eq!(legacy_string_hash(key), 0x1dac_92aa_2553_6c1f);
        assert_eq!(legacy_string_hash(key) % 10_273, 5_747);
    }

    #[test]
    fn filtered_truth_threshold_retention_matches_legacy() {
        // Reduced from the --usefiltered-truth parity lane. The NOCALL/FP
        // and FN/NOCALL pairs are the annotated handoff produced by retained
        // filtered truth calls; together with four TP pairs they exercise the
        // temporary untyped genotype rows and multiple raw-table rehashes.
        let rows = vec![
            annotated(
                "chr1",
                5,
                "60",
                "0/1:TP:gm:ti:SNP:het:58",
                "0/1:TP:gm:ti:SNP:het:58",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                17,
                "55",
                "0/1:TP:gm:tv:SNP:het:53",
                "0/1:TP:gm:tv:SNP:het:53",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                29,
                "48",
                "./.:.:.:.:NOCALL:nocall:.",
                "0/1:FP:.:tv:SNP:het:48",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                53,
                "40",
                "0/1:TP:gm:ti:SNP:het:38",
                "0/1:TP:gm:ti:SNP:het:38",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                65,
                "35",
                "0/1:FN:.:tv:SNP:het:.",
                "./.:.:.:.:NOCALL:nocall:0",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                89,
                "25",
                "0/1:TP:gm:ti:SNP:het:23",
                "0/1:TP:gm:ti:SNP:het:23",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                97,
                "20",
                "./.:.:.:.:NOCALL:nocall:.",
                "0/1:FP:.:tv:SNP:het:20",
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let emitted = groups[&key].emit();

        assert_eq!(
            emitted
                .iter()
                .map(|row| {
                    (
                        row.qq_str.as_str(),
                        row.cum.truth_tp.total,
                        row.cum.truth_fn.total,
                        row.cum.query_tp.total,
                        row.cum.query_fp.total,
                    )
                })
                .collect::<Vec<_>>(),
            [
                ("*", 4, 1, 4, 2),
                ("0.000000", 4, 1, 4, 2),
                ("20.000000", 4, 1, 4, 1),
                ("23.000000", 3, 2, 4, 1),
                ("38.000000", 2, 3, 3, 1),
                ("48.000000", 2, 3, 2, 0),
                ("53.000000", 1, 4, 2, 0),
                ("58.000000", 0, 5, 1, 0),
            ]
        );

        let retained = legacy_metric_raw_order(&groups, 0.5, true);
        assert_eq!(
            retained,
            [
                "SNP\t*\t*\tALL\t*\t48.000000",
                "SNP\t*\t*\tALL\t*\t38.000000",
                "SNP\t*\t*\tALL\t*\t23.000000",
                "SNP\t*\t*\tALL\t*\t20.000000",
                "SNP\t*\t*\tPASS\t*\t58.000000",
                "SNP\t*\t*\tPASS\t*\t53.000000",
                "SNP\t*\t*\tPASS\t*\t20.000000",
                "SNP\t*\t*\tPASS\t*\t*",
                "SNP\t*\t*\tALL\t*\t0.000000",
                "SNP\t*\t*\tALL\t*\t*",
                "SNP\t*\t*\tALL\t*\t53.000000",
                "SNP\t*\t*\tPASS\t*\t38.000000",
                "SNP\t*\t*\tALL\t*\t58.000000",
                "SNP\t*\t*\tPASS\t*\t48.000000",
                "SNP\t*\t*\tPASS\t*\t0.000000",
                "SNP\t*\t*\tPASS\t*\t23.000000",
            ]
            .map(str::to_string)
        );
    }

    #[test]
    fn python27_metrics_dictionary_preserves_pinned_iteration_order() {
        let five_table_insertion = [
            "roc.all",
            "roc.Locations.SNP",
            "roc.Locations.INDEL.PASS",
            "roc.Locations.INDEL",
            "roc.Locations.SNP.PASS",
        ]
        .map(str::to_string);
        assert_eq!(
            python27_dict_iteration_order(&five_table_insertion),
            [
                "roc.Locations.SNP.PASS",
                "roc.all",
                "roc.Locations.INDEL",
                "roc.Locations.SNP",
                "roc.Locations.INDEL.PASS",
            ]
        );

        // The sixth insertion crosses CPython 2.7's two-thirds load factor,
        // so this also pins the resize and rehash path used by SEL reports.
        let seven_table_insertion = [
            "roc.all",
            "roc.Locations.INDEL.SEL",
            "roc.Locations.SNP",
            "roc.Locations.SNP.SEL",
            "roc.Locations.INDEL.PASS",
            "roc.Locations.SNP.PASS",
            "roc.Locations.INDEL",
        ]
        .map(str::to_string);
        assert_eq!(
            python27_dict_iteration_order(&seven_table_insertion),
            [
                "roc.all",
                "roc.Locations.INDEL",
                "roc.Locations.SNP.PASS",
                "roc.Locations.SNP",
                "roc.Locations.INDEL.PASS",
                "roc.Locations.SNP.SEL",
                "roc.Locations.INDEL.SEL",
            ]
        );
    }

    #[test]
    fn introsort_depth_floor_uses_target_pointer_width() {
        assert_eq!(lg_floor(0), 0);
        assert_eq!(lg_floor(1), 0);
        assert_eq!(lg_floor(2), 1);
        assert_eq!(lg_floor(3), 1);
        assert_eq!(lg_floor(16), 4);
        assert_eq!(lg_floor(usize::MAX), usize::BITS as usize - 1);
    }

    #[allow(clippy::too_many_arguments)] // Keeps row fixtures legible at each call site.
    fn annotated(
        chrom: &str,
        pos: usize,
        qual: &str,
        truth_sample: &str,
        query_sample: &str,
        regions: &str,
        query_pass: bool,
        fp_class: Option<&'static str>,
    ) -> AnnotatedRow {
        let regions_tag = if regions.is_empty() {
            String::new()
        } else {
            format!(";Regions={regions}")
        };
        let line = format!(
            "{chrom}\t{pos}\t.\tA\tT\t{qual}\t.\tBS=1{regions_tag}\tGT:BD:BK:BI:BVT:BLT:QQ\t{truth_sample}\t{query_sample}"
        );
        AnnotatedRow {
            sort_key: (chrom.to_string(), pos, 1, 0),
            record: comparison_record(&line),
            query_pass,
            fp_class,
            xcmp_ctype: None,
            xcmp_hap_match: false,
        }
    }

    // Per-side ROC threshold pin: legacy reads each sample's FORMAT.QQ
    // for its own ROC sweep. On chr21 multi-allelic indels the record
    // QUAL is `0` while TRUTH.QQ carries the matched per-side quality
    // (e.g. `829.15`); reading QUAL bins these into the QQ=0 bucket and
    // shifts ~800 truth-TPs out of upper QQ thresholds. Reading
    // FORMAT.QQ keeps each side's records in the right cumulative
    // bucket.
    #[test]
    fn truth_and_query_use_per_side_format_qq() {
        let rows = vec![
            // Multi-allelic-style row: record QUAL=0 but per-side QQ
            // diverge (truth=500, query=0). Truth contributions should
            // land in the 500.000000 bucket, query in 0.000000.
            annotated(
                "chr1",
                100,
                "0",
                "0/1:TP:gm:tv:SNP:het:500",
                "0/1:TP:gm:tv:SNP:het:0",
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");

        // Truth-side bucket at qq=500.0 must hold the TP. (Reading
        // record QUAL=0 would put it in the 0.000000 bucket instead.)
        let bucket_500 = accum
            .numeric_buckets
            .get("500.000000")
            .expect("missing 500.000000 bucket");
        assert_eq!(bucket_500.counts.truth_tp.total, 1);
        assert_eq!(bucket_500.counts.query_tp.total, 0);

        // Query-side bucket at qq=0.0 must hold the matching query TP.
        let bucket_0 = accum
            .numeric_buckets
            .get("0.000000")
            .expect("missing 0.000000 bucket");
        assert_eq!(bucket_0.counts.truth_tp.total, 0);
        assert_eq!(bucket_0.counts.query_tp.total, 1);
    }

    #[test]
    fn cumulates_a_single_snp_group_across_qq_thresholds() {
        // Four rows, all SNP TP/TP matches, in subset="*". QUAL values are
        // 10, 20, 30, and ".". The "." row lands only in the baseline.
        let rows = vec![
            annotated(
                "chr1",
                100,
                "10",
                "0/1:TP:gm:tv:SNP:het:10",
                "0/1:TP:gm:tv:SNP:het:10",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "20",
                "1/1:TP:gm:ti:SNP:homalt:20",
                "1/1:TP:gm:ti:SNP:homalt:20",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                "0/1:TP:gm:ti:SNP:het:30",
                "0/1:TP:gm:ti:SNP:het:30",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                ".",
                "0/1:TP:gm:tv:SNP:het:.",
                "0/1:TP:gm:tv:SNP:het:.",
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        // Baseline sums all 4 rows; the "." row maps to level=0 (mirroring
        // legacy's `if(std::isnan(qq)) qq = 0` in BlockQuantify::observe),
        // producing a fourth numeric threshold at "0.000000". Lex-ASC
        // string order on integer-valued QQs happens to match numeric order
        // here.
        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        assert_eq!(
            qq_strs,
            vec!["*", "0.000000", "10.000000", "20.000000", "30.000000"]
        );

        // Cumulative semantics (strict-above): row at level L reports
        //   tp = total − cum_through_first_at_L
        // The obs vector contains BOTH a truth-side and a query-side
        // observation per row, so 8 records total. The exact tp values
        // at tied levels depend on libstdc++ tie-break ordering between
        // the truth and query siblings, which is implementation-defined
        // and not stable for std::sort. We only assert structural
        // properties: monotonic decreasing across thresholds, baseline
        // holds total, and the highest threshold reaches the boundary.
        let tp_totals: Vec<usize> = emitted.iter().map(|r| r.cum.truth_tp.total).collect();
        assert_eq!(tp_totals[0], 4, "baseline truth_tp should sum all rows");
        for win in tp_totals.windows(2) {
            assert!(win[0] >= win[1], "tp_totals must monotonically decrease");
        }
        let query_totals: Vec<usize> = emitted.iter().map(|r| r.cum.query_tp.total).collect();
        assert_eq!(query_totals[0], 4);
        for win in query_totals.windows(2) {
            assert!(win[0] >= win[1]);
        }
    }

    #[test]
    fn splits_contributions_across_axes_for_pass_snp_with_region() {
        // One SNP PASS row in built-in and named regions contributes to every
        // non-CONF region axis.
        let rows = vec![annotated(
            "chr1",
            100,
            "42",
            "0/1:TP:gm:tv:SNP:het:42",
            "0/1:TP:gm:tv:SNP:het:42",
            "CONF,TS_contained,EXTRA",
            true,
            None,
        )];
        let groups = accumulate(&rows);
        // accumulate now pre-seeds empty baseline entries for every
        // expected (type, subtype, subset, filter) combo so legacy's
        // empty-subtype baseline rows appear in roc.all. Restrict the
        // assertion to SNP groups with the row's observed axes.
        let snp_combos: std::collections::BTreeSet<(String, String, String)> = groups
            .iter()
            .filter(|(k, accum)| {
                k.ty == "SNP" && k.subtype == "*" && accum.baseline.truth_tp.total > 0
            })
            .map(|(k, _)| (k.subset.clone(), k.filter.clone(), k.subtype.clone()))
            .collect();
        let expected: std::collections::BTreeSet<(String, String, String)> = [
            ("*", "ALL", "*"),
            ("*", "PASS", "*"),
            ("TS_contained", "ALL", "*"),
            ("TS_contained", "PASS", "*"),
            ("EXTRA", "ALL", "*"),
            ("EXTRA", "PASS", "*"),
        ]
        .into_iter()
        .map(|(a, b, c)| (a.to_string(), b.to_string(), c.to_string()))
        .collect();
        let combos = snp_combos;
        assert_eq!(combos, expected);
    }

    #[test]
    fn render_emits_expected_cell_formats() {
        // INDEL baseline row with all zero counts: ti/tv cells should be
        // empty, TiTv_ratio should be empty, denom-zero metrics should render
        // "0.0" (matching `metric_ratio`).
        let key = RowKey::new("INDEL", "*", "*", "ALL");
        let emitted = EmittedRow {
            qq_str: "*".to_string(),
            cum: Cumul::default(),
            substats: None,
        };
        let subset_confidence_sizes = BTreeMap::new();
        let subset_sizes = BTreeMap::new();
        let rendered = render_row(
            &key,
            &emitted,
            100,
            140,
            50,
            &subset_sizes,
            &subset_confidence_sizes,
            false,
            0.0,
        );
        let cells: Vec<&str> = rendered.split(',').collect();
        assert_eq!(cells.len(), 65, "expected 65 columns, got {}", cells.len());
        assert_eq!(cells[0], "INDEL");
        assert_eq!(cells[1], "*");
        assert_eq!(cells[6], "*");
        // Zero-denominator metrics render as "0.0" (legacy metric_ratio).
        assert_eq!(cells[7], "0.0");
        // FP.gt / FP.al are raw integers on the base row.
        assert_eq!(cells[11], "0");
        assert_eq!(cells[12], "0");
        // Subset.Size = raw subset_size integer at Subset="*".
        assert_eq!(cells[13], "100");
        // Confidence regions make unsupported INDEL ti/tv cells use `.`.
        assert_eq!(cells[17], ".");
        assert_eq!(cells[18], ".");
        // TiTv_ratio: empty for INDEL.
        assert_eq!(cells[21], "");

        // Het/hom ratio: both zero → empty.
        let het_hom = het_hom_ratio(0, 0);
        assert_eq!(het_hom, "");

        let without_confidence = render_row(
            &key,
            &emitted,
            100,
            140,
            0,
            &subset_sizes,
            &subset_confidence_sizes,
            false,
            0.0,
        );
        let cells = without_confidence.split(',').collect::<Vec<_>>();
        assert_eq!(cells[17], "");
        assert_eq!(cells[18], "");
    }

    #[test]
    fn subset_size_cells_use_boundary_reference_and_named_confidence_intersection() {
        let subset_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);
        let subset_confidence_sizes = BTreeMap::from([("EXTRA".to_string(), 138)]);

        assert_eq!(
            subset_size_cells(
                "TS_boundary",
                "*",
                100,
                140,
                141,
                &subset_sizes,
                &subset_confidence_sizes,
            ),
            ("140.000000".to_string(), "141.000000".to_string())
        );
        assert_eq!(
            subset_size_cells(
                "EXTRA",
                "*",
                100,
                140,
                141,
                &subset_sizes,
                &subset_confidence_sizes,
            ),
            ("138.000000".to_string(), "138.000000".to_string())
        );
        assert_eq!(
            subset_size_cells(
                "*",
                "*",
                100,
                140,
                141,
                &subset_sizes,
                &subset_confidence_sizes,
            ),
            ("100".to_string(), "141.000000".to_string())
        );
    }

    #[test]
    fn absent_variant_type_has_no_roc_rows_or_location_files() {
        let rows = vec![annotated(
            "chr1",
            100,
            "42",
            "0/1:TP:gm:i1_5:INDEL:het:42",
            "0/1:TP:gm:i1_5:INDEL:het:42",
            "",
            true,
            None,
        )];
        let artifacts = calculate(&rows, 100, 0).unwrap();
        let all = artifacts
            .csv
            .iter()
            .find(|artifact| artifact.suffix == "roc.all.csv.gz")
            .unwrap();
        assert!(all.rows.iter().all(|line| line.starts_with("INDEL,")));
        let emitted = artifacts
            .csv
            .iter()
            .filter(|artifact| !artifact.rows.is_empty())
            .map(|artifact| artifact.suffix.as_str())
            .collect::<BTreeSet<_>>();
        assert!(emitted.contains("roc.Locations.INDEL.csv.gz"));
        assert!(emitted.contains("roc.Locations.INDEL.PASS.csv.gz"));
        assert!(!emitted.contains("roc.Locations.SNP.csv.gz"));
        assert!(!emitted.contains("roc.Locations.SNP.PASS.csv.gz"));
    }

    #[test]
    fn configured_unobserved_subset_has_zero_baselines_for_each_active_type() {
        let rows = vec![
            annotated(
                "chr1",
                100,
                "42",
                "0/1:TP:gm:tv:SNP:het:42",
                "0/1:TP:gm:tv:SNP:het:42",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "41",
                "0/1:TP:gm:i1_5:INDEL:het:41",
                "0/1:TP:gm:i1_5:INDEL:het:41",
                "",
                true,
                None,
            ),
        ];
        let options = RocOptions {
            subset_sizes: BTreeMap::from([("EXTRA_unused".to_string(), 7)]),
            ..RocOptions::default()
        };
        let artifacts = calculate_with_options(&rows, 100, 0, &options).unwrap();
        let all = artifacts
            .csv
            .iter()
            .find(|artifact| artifact.suffix == "roc.all.csv.gz")
            .unwrap();
        let unused = all
            .rows
            .iter()
            .map(String::as_str)
            .filter(|line| line.split(',').nth(2) == Some("EXTRA_unused"))
            .map(|line| line.split(',').collect::<Vec<_>>())
            .collect::<Vec<_>>();

        let actual_axes = unused
            .iter()
            .map(|fields| (fields[0], fields[1], fields[3], fields[6]))
            .collect::<BTreeSet<_>>();
        let mut expected_axes = BTreeSet::new();
        for filter in ["ALL", "PASS"] {
            expected_axes.insert(("SNP", "*", filter, "*"));
            expected_axes.insert(("INDEL", "*", filter, "*"));
            for subtype in INDEL_SUBTYPES {
                expected_axes.insert(("INDEL", subtype, filter, "*"));
            }
        }
        assert_eq!(actual_axes, expected_axes);
        assert!(unused.iter().all(|fields| fields[6] == "*"));
        assert!(unused.iter().all(|fields| fields[13] == "7.000000"));
        assert!(unused.iter().all(|fields| fields[16] == "0"));
    }

    #[test]
    fn pass_tier_demotes_tp_when_query_filtered() {
        // Truth/query TP but query_pass=false. Expected: ALL counts truth
        // as TP, PASS counts truth as FN and query contributes nothing.
        let rows = vec![annotated(
            "chr1",
            100,
            "10",
            "0/1:TP:gm:tv:SNP:het:10",
            "0/1:TP:gm:tv:SNP:het:10",
            "",
            false,
            None,
        )];
        let groups = accumulate(&rows);

        let all_key = RowKey::new("SNP", "*", "*", "ALL");
        let pass_key = RowKey::new("SNP", "*", "*", "PASS");
        let all = groups.get(&all_key).expect("ALL group missing").emit();
        let pass = groups.get(&pass_key).expect("PASS group missing").emit();

        // ALL baseline: 1 TP on truth + 1 TP on query.
        assert_eq!(all[0].cum.truth_tp.total, 1);
        assert_eq!(all[0].cum.truth_fn.total, 0);
        assert_eq!(all[0].cum.query_tp.total, 1);

        // PASS baseline: truth demoted to FN, query contributes nothing.
        assert_eq!(pass[0].cum.truth_tp.total, 0);
        assert_eq!(pass[0].cum.truth_fn.total, 1);
        assert_eq!(pass[0].cum.query_tp.total, 0);
    }

    #[test]
    fn pass_tier_ignores_filtered_queries_without_a_roc_decision() {
        let rows = vec![
            annotated(
                "chr1",
                100,
                "10",
                "0/0:.:.:.:NOCALL:homref:.",
                "0/1:UNK:gm:tv:SNP:het:10",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "100",
                "0/0:.:.:.:NOCALL:homref:.",
                "0/1:.:gm:tv:SNP:het:100",
                "",
                false,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let pass_key = RowKey::new("SNP", "*", "*", "PASS");
        let pass = groups.get(&pass_key).expect("PASS group missing");

        // Legacy rocEvaluate calls addROCValue only for TP/FP/UNK query
        // decisions. A filtered record with BD=. must not create an N
        // observation or shift the roc-delta threshold sequence.
        assert_eq!(pass.obs.len(), 1);
        let qq: Vec<_> = pass.emit().into_iter().map(|row| row.qq_str).collect();
        assert_eq!(qq, vec!["*", "10.000000"]);
    }

    #[test]
    fn roc_delta_keeps_well_spaced_rows() {
        // Three SNP contributions: TP/TP matches at QQ=50 and QQ=30, plus a
        // truth-only FN at QQ=40. All three numeric thresholds are >0.5
        // apart, so legacy's roc_delta=0.5 filter keeps all of them — the
        // rust emit must match. This pins the absence of a 7-tuple dedup
        // (which legacy doesn't apply) and preserves the derived-FN
        // semantics across all kept rows.
        let rows = vec![
            annotated(
                "chr1",
                100,
                "50",
                "0/1:TP:gm:tv:SNP:het:50",
                "0/1:TP:gm:tv:SNP:het:50",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "40",
                "0/1:FN:gm:tv:SNP:het:40",
                "0/1:FN:gm:tv:SNP:het:40",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                "0/1:TP:gm:tv:SNP:het:30",
                "0/1:TP:gm:tv:SNP:het:30",
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        // Expected: baseline + three numeric rows (lex-ASC by QQ string).
        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        assert_eq!(qq_strs, vec!["*", "30.000000", "40.000000", "50.000000"]);

        // Baseline cum: 2 truth TP + 1 truth FN + 2 query TP (raw).
        // Numeric rows use legacy's strict-above formula `tp = total −
        // cum_through_first_at_L`. Exact per-row tp values depend on
        // libstdc++ cluster-tie-break ordering between truth and query
        // observations at the same level, which is implementation-
        // defined. Assert structural properties only.
        let tp_totals: Vec<usize> = emitted.iter().map(|r| r.cum.truth_tp.total).collect();
        assert_eq!(tp_totals[0], 2, "baseline truth_tp should sum the two TPs");
        for win in tp_totals.windows(2) {
            assert!(win[0] >= win[1]);
        }
        let fn_totals: Vec<usize> = emitted.iter().map(|r| r.cum.truth_fn.total).collect();
        assert_eq!(fn_totals[0], 1, "baseline truth_fn should sum the one FN");
        for win in fn_totals.windows(2) {
            assert!(win[0] <= win[1]);
        }
        let qtp_totals: Vec<usize> = emitted.iter().map(|r| r.cum.query_tp.total).collect();
        assert_eq!(qtp_totals[0], 2);
        for win in qtp_totals.windows(2) {
            assert!(win[0] >= win[1]);
        }
    }

    #[test]
    fn roc_delta_drops_rows_within_half_unit() {
        // Five SNP query-TP contributions at closely-spaced QQs. Legacy
        // --roc-delta 0.5 keeps only rows where QQ moves >0.5 from the
        // last kept row (ASC walk, first always kept). Expected kept
        // numeric QQs: 30.0, 30.6, 31.2 — 30.3 and 30.9 are within 0.5
        // of the prior kept row and must be dropped. This is the
        // dominant source of rust/legacy roc.all row-count divergence
        // before the fix (~6700 → ~2572 rows per SNP group).
        let rows = vec![
            annotated(
                "chr1",
                100,
                "30.0",
                "0/1:TP:gm:tv:SNP:het:30",
                "0/1:TP:gm:tv:SNP:het:30",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "30.3",
                "0/1:TP:gm:tv:SNP:het:30.3",
                "0/1:TP:gm:tv:SNP:het:30.3",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30.6",
                "0/1:TP:gm:tv:SNP:het:30.6",
                "0/1:TP:gm:tv:SNP:het:30.6",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                "30.9",
                "0/1:TP:gm:tv:SNP:het:30.9",
                "0/1:TP:gm:tv:SNP:het:30.9",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                500,
                "31.2",
                "0/1:TP:gm:tv:SNP:het:31.2",
                "0/1:TP:gm:tv:SNP:het:31.2",
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        // 31.2 rounds through f32 to 31.2000007629... which formats to
        // "31.200001" at 6 decimals — legacy exhibits the same artifact.
        assert_eq!(qq_strs, vec!["*", "30.000000", "30.600000", "31.200001"]);
    }

    #[test]
    fn no_dedup_when_every_qq_changes_the_cumulative_tuple() {
        // Four SNP TP/TP matches at distinct QQs. Every threshold moves both
        // cum.truth_tp.total and cum.query_tp.total by 1, so the 7-tuple is
        // different at every row. Dedup must keep all four numeric rows
        // alongside the baseline — guards against over-dedup regressions.
        let rows = vec![
            annotated(
                "chr1",
                100,
                "10",
                "0/1:TP:gm:tv:SNP:het:10",
                "0/1:TP:gm:tv:SNP:het:10",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                200,
                "20",
                "0/1:TP:gm:ti:SNP:het:20",
                "0/1:TP:gm:ti:SNP:het:20",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                300,
                "30",
                "1/1:TP:gm:ti:SNP:homalt:30",
                "1/1:TP:gm:ti:SNP:homalt:30",
                "",
                true,
                None,
            ),
            annotated(
                "chr1",
                400,
                "40",
                "0/1:TP:gm:tv:SNP:het:40",
                "0/1:TP:gm:tv:SNP:het:40",
                "",
                true,
                None,
            ),
        ];
        let groups = accumulate(&rows);
        let key = RowKey::new("SNP", "*", "*", "ALL");
        let accum = groups.get(&key).expect("SNP/*/*/ALL group missing");
        let emitted = accum.emit();

        assert_eq!(emitted.len(), 5);
        let qq_strs: Vec<String> = emitted.iter().map(|r| r.qq_str.clone()).collect();
        assert_eq!(
            qq_strs,
            vec!["*", "10.000000", "20.000000", "30.000000", "40.000000",]
        );
    }

    #[test]
    fn custom_roc_field_and_delta_control_thresholds() {
        let mut first = annotated(
            "chr1",
            100,
            "90",
            "0/1:TP:gm:tv:SNP:het:90",
            "0/1:TP:gm:tv:SNP:het:90",
            "",
            true,
            None,
        );
        first
            .record
            .try_update(|record| {
                record.info = record.info.replace("BS=1", "BS=1;SCORE=10.0");
                Ok(())
            })
            .unwrap();
        let mut second = annotated(
            "chr1",
            200,
            "80",
            "0/1:TP:gm:tv:SNP:het:80",
            "0/1:TP:gm:tv:SNP:het:80",
            "",
            true,
            None,
        );
        second
            .record
            .try_update(|record| {
                record.info = record.info.replace("BS=1", "BS=1;SCORE=10.4");
                Ok(())
            })
            .unwrap();
        let options = RocOptions {
            qq_field: "SCORE".to_string(),
            delta: 0.0,
            ..RocOptions::default()
        };
        let groups = accumulate_with_options(&[first, second], &options);
        let key = RowKey::new_with_qq_field("SNP", "*", "*", "ALL", "SCORE");
        let rows = groups[&key].emit_with_delta(options.delta);
        assert_eq!(
            rows.iter()
                .map(|row| row.qq_str.as_str())
                .collect::<Vec<_>>(),
            vec!["*", "10.000000", "10.400000"]
        );
        assert!(!groups[&key].numeric_buckets.contains_key("90.000000"));
        let sorted = build_star_sorted(&groups);
        for subtype in ["*", "ti", "tv"] {
            assert!(
                sorted.contains_key(&(
                    "SNP".to_string(),
                    subtype.to_string(),
                    "*".to_string(),
                    "ALL".to_string(),
                )),
                "custom ROC fields need the same subtype sort snapshots as QUAL"
            );
        }
    }

    #[test]
    fn roc_threshold_source_can_differ_from_reported_field() {
        let row = annotated(
            "chr1",
            100,
            "90",
            "0/1:TP:gm:tv:SNP:het:90",
            "0/1:TP:gm:tv:SNP:het:90",
            "",
            true,
            None,
        );
        let options = RocOptions {
            qq_field: "INFO.SCORE".to_string(),
            score_field: Some("QQ".to_string()),
            delta: 0.0,
            ..RocOptions::default()
        };

        let groups = accumulate_with_options(&[row], &options);
        let key = RowKey::new_with_qq_field("SNP", "*", "*", "ALL", "INFO.SCORE");
        let rows = groups[&key].emit_with_delta(options.delta);

        assert_eq!(
            rows.iter()
                .map(|row| row.qq_str.as_str())
                .collect::<Vec<_>>(),
            vec!["*", "90.000000"]
        );
    }

    #[test]
    fn ignored_filter_creates_selective_tier_and_renamed_filter_counts() {
        let mut row = annotated(
            "chr1",
            100,
            "20",
            "0/1:TP:gm:tv:SNP:het:20",
            "0/1:TP:gm:tv:SNP:het:20",
            "",
            false,
            None,
        );
        row.record
            .try_update(|record| {
                record.filter = "LowQual".to_string();
                Ok(())
            })
            .unwrap();
        let options = RocOptions {
            ignored_filters: HashSet::from(["LowQual".to_string()]),
            ..RocOptions::default()
        };
        let groups = accumulate_with_options(&[row], &options);
        let pass = &groups[&RowKey::new_with_qq_field("SNP", "*", "*", "PASS", "QUAL")];
        let selective = &groups[&RowKey::new_with_qq_field("SNP", "*", "*", "SEL", "QUAL")];
        assert_eq!(pass.baseline.truth_fn.total, 1);
        assert_eq!(pass.baseline.query_tp.total, 0);
        assert_eq!(selective.baseline.truth_tp.total, 1);
        assert_eq!(selective.baseline.query_tp.total, 1);
        let ignored =
            &groups[&RowKey::new_with_qq_field("SNP", "*", "*", "SEL_IGN_LowQual", "QUAL")];
        assert_eq!(ignored.baseline.query_tp.total, 1);
    }

    #[test]
    fn roc_regions_preserve_legacy_filter_sweep_switch() {
        let mut row = annotated(
            "chr1",
            100,
            "20",
            "0/1:TP:gm:tv:SNP:het:20",
            "0/1:TP:gm:tv:SNP:het:20",
            "",
            false,
            None,
        );
        row.record
            .try_update(|record| {
                record.filter = "LowQual".to_string();
                Ok(())
            })
            .unwrap();
        let options = RocOptions {
            roc_regions: HashSet::from(["TS_contained".to_string()]),
            ..RocOptions::default()
        };
        let groups = accumulate_with_options(&[row], &options);
        let sorted = build_star_sorted(&groups);
        let rendered = render_rows(
            &groups,
            &sorted,
            RowFilter::All,
            RenderConfig {
                subset_size: 100,
                whole_reference_size: 100,
                conf_size: 0,
                subset_sizes: &BTreeMap::new(),
                subset_confidence_sizes: &BTreeMap::new(),
                delta: 0.5,
                ci_alpha: 0.0,
                filter_counts_only: options.roc_regions.contains("*"),
            },
        );
        assert!(rendered.iter().any(|line| {
            let fields = line.split(',').collect::<Vec<_>>();
            fields[3] == "LowQual" && fields[6] != "*"
        }));
    }

    #[test]
    fn confidence_interval_columns_and_modified_jeffreys_edges() {
        let header = roc_header(0.05);
        assert_eq!(header.split(',').count(), EXTENDED_HEADER.len() + 6);
        let (empty_lower, empty_upper) = jeffreys_interval(0, 0, 0.05);
        assert_eq!((empty_lower, empty_upper), (0.0, 1.0));
        let (lower, upper) = jeffreys_interval(0, 10, 0.05);
        assert_eq!(lower, 0.0);
        assert!((upper - (1.0 - 0.025_f64.powf(0.1))).abs() < 1e-14);
        let (lower, upper) = jeffreys_interval(5, 10, 0.05);
        assert!((lower - 0.223_528_670_252_705_2).abs() < 1e-12);
        assert!((upper - 0.776_471_329_747_294_7).abs() < 1e-12);
        let (lower, upper) = jeffreys_interval(1037, 1156, 0.05);
        assert_eq!(lower.to_bits(), 0x3fec_1d15_8e8d_75af);
        assert_eq!(upper.to_bits(), 0x3fed_3c11_2b31_7194);

        let (lower, _) = jeffreys_interval(1233, 1233, 0.05);
        assert_eq!(lower.to_bits(), 0x3fef_e787_2242_48a7);
        let (_, upper) = jeffreys_interval(0, 136, 0.05);
        assert_eq!(upper.to_bits(), 0x3f9b_66db_9060_1320);
    }

    #[test]
    fn filter_tier_unknown_fraction_ci_uses_zero_query_total() {
        let key = RowKey::new("SNP", "*", "*", "LowMQ");
        let emitted = EmittedRow {
            qq_str: "*".to_string(),
            cum: Cumul {
                truth_tp: CountsBucket {
                    total: 2,
                    ..CountsBucket::default()
                },
                query_tp: CountsBucket {
                    total: 1,
                    ..CountsBucket::default()
                },
                query_fp: CountsBucket {
                    total: 3,
                    ..CountsBucket::default()
                },
                query_unk: CountsBucket {
                    total: 4,
                    ..CountsBucket::default()
                },
                ..Cumul::default()
            },
            substats: None,
        };

        let subset_confidence_sizes = BTreeMap::new();
        let subset_sizes = BTreeMap::new();
        let rendered = render_row(
            &key,
            &emitted,
            100,
            140,
            50,
            &subset_sizes,
            &subset_confidence_sizes,
            true,
            0.05,
        );
        let cells = rendered.split(',').collect::<Vec<_>>();
        assert_eq!(
            &cells[cells.len() - 6..],
            [
                "0.158113883008",
                "1.0",
                "0.0",
                "0.716248320437",
                "0.0",
                "1.0",
            ]
        );
    }
}
