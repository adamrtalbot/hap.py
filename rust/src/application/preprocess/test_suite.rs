//! Command-level regression tests.

#[cfg(test)]
mod tests {
    use super::super::genotype::expand_haploid_gt;
    use super::super::*;
    use crate::application::SomaticGtMode;
    use crate::domain::RawVcfRecord;
    use proptest::prelude::*;
    use std::fs;
    use tempfile::tempdir;

    fn run(args: PreprocessArgs) -> anyhow::Result<()> {
        super::super::run(args.validated_with_legacy_plain_vcf()?)
    }

    #[test]
    fn filters_only_removes_records_carrying_only_selected_filters() {
        assert!(passes_filters_only("PASS", Some("LowQual,q10")));
        assert!(passes_filters_only(".", Some("LowQual,q10")));
        assert!(!passes_filters_only("LowQual", Some("LowQual,q10")));
        assert!(!passes_filters_only("LowQual;q10", Some("LowQual,q10")));
        assert!(passes_filters_only("LowQual;s50", Some("LowQual,q10")));
        assert!(passes_filters_only("LowQual", None));
    }

    #[test]
    fn gender_auto_matches_vcfcheck_haploid_x_heuristic() {
        let mut haploid = make_record(".");
        haploid.chrom = "chrX".to_string();
        haploid.format = Some("GT".to_string());
        haploid.samples = vec!["1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[haploid.clone()]),
            PreprocessGender::Male
        );

        let mut heterozygous = haploid.clone();
        heterozygous.samples = vec!["0/1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[haploid, heterozygous]),
            PreprocessGender::Female
        );
        assert_eq!(
            resolve_gender(PreprocessGender::None, &[]),
            PreprocessGender::None
        );
    }

    #[test]
    fn gender_auto_treats_half_called_x_genotype_as_diploid() {
        let mut half_called = make_record(".");
        half_called.chrom = "chrX".to_string();
        half_called.format = Some("GT".to_string());
        half_called.samples = vec!["./1".to_string()];

        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[half_called]),
            PreprocessGender::Female,
            "legacy vcfcheck uses ngt=2 and compares the missing and called slots"
        );
    }

    #[test]
    fn gender_auto_does_not_treat_lowercase_x_as_x_chromosome() {
        let mut lowercase_x = make_record(".");
        lowercase_x.chrom = "x".to_string();
        lowercase_x.format = Some("GT".to_string());
        lowercase_x.samples = vec!["1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[lowercase_x]),
            PreprocessGender::Female,
            "pinned vcfcheck compares its location variable to lowercase x"
        );

        let mut lowercase_chrx = make_record(".");
        lowercase_chrx.chrom = "chrx".to_string();
        lowercase_chrx.format = Some("GT".to_string());
        lowercase_chrx.samples = vec!["1".to_string()];
        assert_eq!(
            resolve_gender(PreprocessGender::Auto, &[lowercase_chrx]),
            PreprocessGender::Male
        );
    }

    #[test]
    fn region_and_auto_fixchr_complete_half_called_x_genotype_like_legacy() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let regions = directory.path().join("regions.bed");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=X,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "X\t2\tx_half_called\tA\tC\t20\tPASS\t.\tGT:AD\t./1:1,9\n",
            ),
        )?;
        fs::write(&reference, ">chrX\nAAAAA\n")?;
        fs::write(&regions, "chrX\t0\t5\n")?;

        let mut args = interval_args(&input, &output, &reference, Some(&regions), None);
        args.fixchr = None;
        args.gender = PreprocessGender::Auto;
        args.decompose = true;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].to_line(),
            "chrX\t2\t.\tA\tC\t20\t.\t.\tGT:AD:ADO:DP\t0/1:1,9:1:0"
        );
        Ok(())
    }

    #[test]
    fn male_gender_duplicates_haploid_x_and_y_only() {
        let mut x = make_record(".");
        x.chrom = "X".to_string();
        x.format = Some("GT:DP".to_string());
        x.samples = vec!["1:9".to_string()];
        expand_male_sex_chromosome_genotypes(&mut x);
        assert_eq!(x.samples, vec!["1/1:9"]);

        let mut autosome = x.clone();
        autosome.chrom = "1".to_string();
        autosome.samples = vec!["1:9".to_string()];
        expand_male_sex_chromosome_genotypes(&mut autosome);
        assert_eq!(autosome.samples, vec!["1:9"]);
    }

    #[test]
    fn reference_candidates_follow_legacy_hg19_then_hgref_precedence() -> Result<()> {
        let directory = tempdir()?;
        let explicit = directory.path().join("explicit.fa");
        let hg19 = directory.path().join("hg19.fa");
        let hgref = directory.path().join("hgref.fa");
        let fallback = directory.path().join("fallback.fa");
        for path in [&explicit, &hg19, &hgref, &fallback] {
            fs::write(path, ">chr1\nA\n")?;
        }
        assert_eq!(
            resolve_reference_candidates(Some(&explicit), Some(&hg19), Some(&hgref), &fallback)?,
            explicit
        );
        assert_eq!(
            resolve_reference_candidates(None, Some(&hg19), Some(&hgref), &fallback)?,
            hg19
        );
        fs::remove_file(&hg19)?;
        assert_eq!(
            resolve_reference_candidates(None, Some(&hg19), Some(&hgref), &fallback)?,
            hgref
        );
        Ok(())
    }

    #[test]
    fn legacy_fixchr_auto_only_adds_a_missing_prefix() {
        let prefixed = BTreeSet::from(["chr1".to_string(), "chrX".to_string()]);
        let plain = BTreeSet::from(["1".to_string(), "X".to_string()]);
        assert!(resolve_fixchr(None, &prefixed, &plain));
        assert!(!resolve_fixchr(None, &plain, &prefixed));
        assert!(!resolve_fixchr(None, &prefixed, &prefixed));

        assert_eq!(add_legacy_chr_prefix("1"), "chr1");
        assert_eq!(add_legacy_chr_prefix("MT"), "chrM");
        assert_eq!(add_legacy_chr_prefix("chrMT"), "chrM");
        assert_eq!(add_legacy_chr_prefix("GL000207.1"), "GL000207.1");
        assert_eq!(add_legacy_chr_prefix("chr1"), "chr1");
    }

    #[test]
    fn fixchr_adds_lengthless_headers_for_rewritten_contigs() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">1\nAAAAA\n>chr1\nAAAAA\n>chrX\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "1\t2\trs1\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.fixchr = None;
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        run(args)?;

        let (headers, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records[0].chrom, "chr1");
        assert!(
            headers
                .iter()
                .any(|line| line == "##contig=<ID=1,length=5>")
        );
        assert!(headers.iter().any(|line| line == "##contig=<ID=chr1>"));
        Ok(())
    }

    #[test]
    fn disabling_leftshift_and_decomposition_preserves_bcftools_view_shape() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FILTER=<ID=LowQual,Description=\"low\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\trs1\tA\tC\t10\tLowQual\tAC=1\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        run(args)?;
        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "rs1");
        assert_eq!(records[0].filter, "LowQual");
        assert_eq!(records[0].info, "AC=1");
        assert_eq!(records[0].format.as_deref(), Some("GT"));
        assert_eq!(records[0].samples, vec!["0/1"]);
        Ok(())
    }

    #[test]
    fn sites_only_vcf_is_rejected_like_vcfcheck() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("sites.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\n",
            ),
        )?;

        let error = run(interval_args(&input, &output, &reference, None, None)).unwrap_err();
        assert!(error.to_string().contains("no samples"));
        assert!(!output.exists());
        Ok(())
    }

    #[test]
    fn disabled_normalization_passes_through_ref_mismatch() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\tmismatch\tC\tT\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;

        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].ref_allele, "C");
        assert_eq!(records[0].id, "mismatch");
        Ok(())
    }

    #[test]
    fn missing_reference_index_is_rejected() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let args = interval_args(&input, &output, &reference, None, None);
        fs::remove_file(format!("{}.fai", reference.display()))?;
        let error = run(args).unwrap_err();
        assert!(error.to_string().contains("is not indexed"));
        assert!(!output.exists());
        Ok(())
    }

    #[test]
    fn empty_reference_index_is_accepted_but_malformed_lengths_are_rejected() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let index = format!("{}.fai", reference.display());

        fs::write(&index, "")?;
        run(interval_args(&input, &output, &reference, None, None))?;
        assert!(output.exists());
        fs::remove_file(&output)?;
        let output_index = PathBuf::from(format!("{}.tbi", output.display()));
        if output_index.exists() {
            fs::remove_file(output_index)?;
        }

        fs::write(&index, "chr1\tbad\t6\t5\t6\n")?;
        let error = run(interval_args(&input, &output, &reference, None, None)).unwrap_err();
        assert!(
            error.to_string().contains("invalid FASTA index length"),
            "{error:#}"
        );
        assert!(!output.exists());
        Ok(())
    }

    #[test]
    fn plain_vcf_output_failure_leaves_unindexed_output() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        args.threads = Some(1);
        let error = run(args).unwrap_err();
        assert!(error.to_string().contains("plain VCF output"));
        assert!(output.is_file());
        assert_eq!(vcf::load_raw_vcf(&output)?.1.len(), 1);
        assert!(!PathBuf::from(format!("{}.tbi", output.display())).exists());
        assert!(!PathBuf::from(format!("{}.csi", output.display())).exists());
        Ok(())
    }

    #[test]
    fn missing_output_parent_is_rejected_without_artifacts() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output_parent = directory.path().join("missing");
        let output = output_parent.join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let error = run(interval_args(&input, &output, &reference, None, None)).unwrap_err();
        assert!(error.to_string().contains("output parent does not exist"));
        assert!(!output_parent.exists());
        Ok(())
    }

    #[test]
    fn negative_window_is_accepted_when_blocksplit_is_unused() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        write_test_fai(&reference)?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        args.threads = Some(1);
        args.window_size = -1;
        run(args)?;

        assert_eq!(vcf::load_raw_vcf(&output)?.1.len(), 1);
        Ok(())
    }

    #[test]
    fn explicit_gender_controls_haploid_x_expansion() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let male_output = directory.path().join("male.vcf.gz");
        let female_output = directory.path().join("female.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chrX\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chrX,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chrX\t2\t.\tA\tC\t10\tPASS\t.\tGT:AD\t1:1,9\n",
            ),
        )?;
        let mut male = interval_args(&input, &male_output, &reference, None, None);
        male.gender = PreprocessGender::Male;
        run(male)?;
        let mut female = interval_args(&input, &female_output, &reference, None, None);
        female.gender = PreprocessGender::Female;
        run(female)?;
        let (_, male_records) = vcf::load_raw_vcf(&male_output)?;
        let (_, female_records) = vcf::load_raw_vcf(&female_output)?;
        assert!(male_records[0].samples[0].starts_with("1/1:"));
        assert!(female_records[0].samples[0].starts_with("0/1:"));
        Ok(())
    }

    #[test]
    fn bcftools_norm_excludes_ref_mismatches_deduplicates_and_left_aligns() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t3\tfirst\tA\tAA\t10\tPASS\t.\tGT\t0/1\n",
                "chr1\t3\tduplicate\tA\tAA\t20\tPASS\t.\tGT\t0/1\n",
                "chr1\t4\tmismatch\tC\tT\t30\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.leftshift = false;
        args.no_leftshift = true;
        args.decompose = false;
        args.no_decompose = true;
        args.gender = PreprocessGender::None;
        args.bcftools_norm = true;
        run(args)?;
        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].id, "first");
        assert_eq!(records[0].pos, 1);
        assert_eq!(records[0].ref_allele, "A");
        assert_eq!(records[0].alt_allele, "AA");
        Ok(())
    }

    #[test]
    fn logfile_and_verbose_emit_operational_messages() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        let logfile = directory.path().join("pre.log");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.logfile = Some(logfile.display().to_string());
        args.verbose = true;
        run(args)?;
        let log = fs::read_to_string(logfile)?;
        assert!(log.contains("INFO Preprocessing"));
        assert!(log.contains("INFO Wrote 0 records"));
        Ok(())
    }

    fn interval_args(
        input: &Path,
        output: &Path,
        reference: &Path,
        regions: Option<&Path>,
        targets: Option<&Path>,
    ) -> PreprocessArgs {
        let mut index = reference.as_os_str().to_os_string();
        index.push(".fai");
        if !Path::new(&index).is_file() {
            write_test_fai(reference).expect("test reference index should be writable");
        }
        PreprocessArgs {
            input: input.display().to_string(),
            output: output.display().to_string(),
            version: false,
            reference: Some(reference.display().to_string()),
            locations: None,
            pass_only: false,
            filters_only: None,
            regions_bedfile: regions.map(|path| path.display().to_string()),
            targets_bedfile: targets.map(|path| path.display().to_string()),
            fixchr: Some(false),
            no_fixchr: false,
            somatic: false,
            set_gt: None,
            filter_nonref: false,
            convert_gvcf_to_vcf: false,
            bcf: false,
            bcftools_norm: false,
            leftshift: true,
            no_leftshift: false,
            decompose: false,
            no_decompose: false,
            gender: PreprocessGender::Auto,
            window_size: 10_000,
            threads: None,
            logfile: None,
            verbose: false,
            quiet: false,
            force_interactive: false,
        }
    }

    fn write_test_fai(reference: &Path) -> Result<()> {
        let sequences = fasta::read_sequences(reference)?;
        let mut offset = 0u64;
        let mut index = String::new();
        for (name, sequence) in sequences {
            let line_bases = sequence.len().max(1);
            index.push_str(&format!(
                "{name}\t{}\t{offset}\t{line_bases}\t{}\n",
                sequence.len(),
                line_bases + 1
            ));
            offset += sequence.len() as u64 + name.len() as u64 + 3;
        }
        fs::write(format!("{}.fai", reference.display()), index)?;
        Ok(())
    }

    #[test]
    fn region_uses_gvcf_end_while_target_uses_start_position() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let boundary = directory.path().join("boundary.bed");
        let region_output = directory.path().join("region.vcf.gz");
        let target_output = directory.path().join("target.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "chr1\t3\t.\tC\tT\t.\tPASS\tEND=5\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;
        fs::write(&boundary, "chr1\t4\t5\n")?;

        run(interval_args(
            &input,
            &region_output,
            &reference,
            Some(&boundary),
            None,
        ))?;
        run(interval_args(
            &input,
            &target_output,
            &reference,
            None,
            Some(&boundary),
        ))?;

        assert_eq!(vcf::load_raw_vcf(&region_output)?.1.len(), 1);
        assert!(vcf::load_raw_vcf(&target_output)?.1.is_empty());
        Ok(())
    }

    #[test]
    fn selectors_remain_literal_after_fixchr_rewrites_records() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let selector = directory.path().join("selector.bed");
        let location_output = directory.path().join("location.vcf.gz");
        let region_output = directory.path().join("region.vcf.gz");
        let target_output = directory.path().join("target.vcf.gz");
        fs::write(&reference, ">chr1\nAAAAA\n")?;
        fs::write(&selector, "1\t0\t5\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=1,length=5>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "1\t2\t.\tA\tC\t10\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let mut location_args = interval_args(&input, &location_output, &reference, None, None);
        location_args.fixchr = Some(true);
        location_args.locations = Some("1:1-5".to_string());
        run(location_args)?;

        let mut region_args =
            interval_args(&input, &region_output, &reference, Some(&selector), None);
        region_args.fixchr = Some(true);
        run(region_args)?;

        let mut target_args =
            interval_args(&input, &target_output, &reference, None, Some(&selector));
        target_args.fixchr = Some(true);
        run(target_args)?;

        for output in [location_output, region_output, target_output] {
            assert!(
                vcf::load_raw_vcf(&output)?.1.is_empty(),
                "selector contig 1 must not be rewritten to match emitted chr1 records"
            );
        }
        Ok(())
    }

    #[test]
    fn leftshift_switch_controls_repeat_indel_normalization() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let shifted_output = directory.path().join("shifted.vcf.gz");
        let unchanged_output = directory.path().join("unchanged.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "chr1\t4\t.\tAA\tA\t30\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAAAAA\n")?;

        let mut shifted = interval_args(&input, &shifted_output, &reference, None, None);
        shifted.leftshift = true;
        run(shifted)?;
        let mut unchanged = interval_args(&input, &unchanged_output, &reference, None, None);
        unchanged.leftshift = false;
        run(unchanged)?;

        let shifted_record = vcf::load_raw_vcf(&shifted_output)?.1.remove(0);
        let unchanged_record = vcf::load_raw_vcf(&unchanged_output)?.1.remove(0);
        assert_eq!(shifted_record.pos, 1);
        assert_eq!(unchanged_record.pos, 4);
        Ok(())
    }

    #[test]
    fn tiny_parallel_inputs_do_not_create_window_block_boundaries() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, ">chr1\nAAAAAAAAAAAAAAAAAAAA\n")?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=20>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t3\tbarrier\tA\tC\t30\tPASS\t.\tGT\t0/1\n",
                "chr1\t5\trepeat\tAA\tA\t30\tPASS\t.\tGT\t0/1\n",
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.threads = Some(2);
        args.window_size = 1;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(
            records.iter().map(|record| record.pos).collect::<Vec<_>>(),
            [3, 3]
        );
        Ok(())
    }

    #[test]
    fn location_start_reset_survives_a_dropped_homref_record() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(120)))?;
        fs::write(
            &input,
            format!(
                concat!(
                    "##fileformat=VCFv4.2\n",
                    "##contig=<ID=chr1,length=120>\n",
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                    "chr1\t2\tbarrier\t{}\tA\t30\tPASS\t.\tGT\t0/1\n",
                    "chr1\t52\thomref\tA\tC\t30\tPASS\t.\tGT\t0/0\n",
                    "chr1\t60\trepeat\tAA\tA\t30\tPASS\t.\tGT\t0/1\n",
                ),
                "A".repeat(50)
            ),
        )?;
        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.threads = Some(2);
        args.locations = Some("chr1:2-2,chr1:52-100".into());
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        // Location aggregation may pad the shifted deletion to a longer REF
        // when another job emits an allele at the same normalized position.
        // Match the one-base deletion by allele length instead of requiring
        // its pre-aggregation `AA>A` spelling.
        assert!(records.iter().any(|record| {
            record.pos == 1
                && record
                    .alt_allele
                    .split(',')
                    .any(|alt| record.ref_allele.len() == alt.len() + 1)
        }));
        Ok(())
    }

    #[test]
    fn symbolic_deletion_uses_end_and_reference_like_variant_reader() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=10>\n",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                // The deliberately incorrect one-base REF is accepted by the
                // legacy reader when END is present; it rebuilds the complete
                // deletion allele from the reference instead.
                "chr1\t2\trs1\tT\t<DEL>\t30\tPASS\tEND=5\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].to_line(),
            "chr1\t1\t.\tAACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0"
        );
        Ok(())
    }

    #[test]
    fn multi_allelic_symbolic_deletion_splits_like_variant_allele_splitter() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        let decomposed_output = directory.path().join("decomposed.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=10>\n",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t2\trs1\tT\t<DEL>,T\t30\tPASS\tEND=5\tGT\t1/2\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.decompose = false;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(
            records
                .iter()
                .map(RawVcfRecord::to_line)
                .collect::<Vec<_>>(),
            [
                "chr1\t1\t.\tAACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0",
                "chr1\t2\t.\tACCG\tT\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0",
            ]
        );

        let mut decomposed = interval_args(&input, &decomposed_output, &reference, None, None);
        decomposed.gender = PreprocessGender::None;
        decomposed.decompose = true;
        run(decomposed)?;
        let (_, records) = vcf::load_raw_vcf(&decomposed_output)?;
        assert_eq!(
            records
                .iter()
                .map(RawVcfRecord::to_line)
                .collect::<Vec<_>>(),
            [
                "chr1\t1\t.\tAACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0",
                "chr1\t2\t.\tA\tT\t30\t.\t.\tGT:AD:ADO:DP\t1/0:.,.:0:0",
                "chr1\t2\t.\tACCG\tA\t30\t.\t.\tGT:AD:ADO:DP\t1/0:.,.:0:0",
            ]
        );
        Ok(())
    }

    #[test]
    fn contig_start_symbolic_deletion_keeps_symbolic_alt() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=10>\n",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t1\trs1\tT\t<DEL>\t30\tPASS\tEND=4\tGT\t0/1\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAACCGGTTAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(
            records[0].to_line(),
            "chr1\t1\t.\tAACC\t<DEL>\t30\t.\t.\tGT:AD:ADO:DP\t0/1:.,.:0:0"
        );
        Ok(())
    }

    #[test]
    fn legacy_headers_replace_conflicting_core_definitions() {
        let input = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"input AD\">".to_string(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"input DP\">".to_string(),
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"input END\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tINPUT".to_string(),
        ];
        let output = canonicalize_legacy_headers(&input);

        assert_eq!(output[0], "##fileformat=VCFv4.1");
        assert!(output.contains(
            &"##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depths\">".to_string()
        ));
        assert!(output.contains(&"##FORMAT=<ID=ADO,Number=.,Type=Integer,Description=\"Summed depth of non-called alleles.\">".to_string()));
        assert!(output.contains(
            &"##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">".to_string()
        ));
        assert!(!output.iter().any(|line| line.contains("input AD")
            || line.contains("input DP")
            || line.contains("input END")));
        assert_eq!(output.last(), input.last());
    }

    #[test]
    fn legacy_headers_sort_and_merge_non_core_input_declarations() {
        let input = vec![
            "##INFO=<ID=ZZ,Number=1,Type=String,Description=\"second\">".to_string(),
            "##FILTER=<ID=LowQ,Description=\"low quality\">".to_string(),
            "##INFO=<ID=ZZ,Number=1,Type=String,Description=\"first\">".to_string(),
            "##custom=z".to_string(),
            "##custom=a".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let output = canonicalize_legacy_headers(&input);
        let appended = &output[LEGACY_BASE_HEADERS.len()..output.len() - 1];

        assert_eq!(
            appended,
            [
                "##FILTER=<ID=LowQ,Description=\"low quality\">",
                "##INFO=<ID=ZZ,Number=1,Type=String,Description=\"first\">",
                "##custom=a",
                "##custom=z",
            ]
        );
    }

    #[test]
    fn ref_bytes_equal_handles_soft_masked_reference() {
        // Reference soft-masks repeat regions as lowercase; the VCF REF stays
        // uppercase. Legacy hap.py accepts this; we must too.
        assert!(ref_bytes_equal(b"t", b"T"));
        assert!(ref_bytes_equal(b"ACGT", b"acgt"));
    }

    #[test]
    fn breakend_normalization_preserves_remote_contig_spelling() {
        assert_eq!(
            uppercase_alleles_preserving_breakends("t]chr1:70],a"),
            "T]chr1:70],A"
        );

        let mut record = make_record(".");
        record.pos = 40;
        record.ref_allele = "T".to_string();
        record.alt_allele = "T]chr1:70]".to_string();
        record.info = "SVTYPE=BND".to_string();
        record.format = Some("GT:AD:GQ".to_string());
        record.samples = vec!["0/1:9,7:45".to_string()];
        normalize_bcftools_record(&mut record, b"ACGTACGTACGT");
        assert_eq!(record.alt_allele, "T]chr1:70]");
        assert!(materialize_unsupported_import_failure(&mut record));
        sort_info_keys(&mut record);
        assert_eq!(record.alt_allele, ".");
        assert_eq!(record.info, "END=40;IMPORT_FAIL;SVTYPE=BND");
        assert_eq!(record.samples, ["0/0:9:45"]);

        let mut symbolic = make_record("END=31");
        symbolic.pos = 31;
        symbolic.ref_allele = "T".to_string();
        symbolic.alt_allele = "<INS>".to_string();
        symbolic.format = Some("GT:AD".to_string());
        symbolic.samples = vec!["0/1:11,4".to_string()];
        assert!(materialize_unsupported_import_failure(&mut symbolic));
        assert_eq!(symbolic.alt_allele, ".");
        assert_eq!(symbolic.info, "END=31;IMPORT_FAIL");
        assert_eq!(symbolic.samples, ["0/0:11"]);
    }

    #[test]
    fn ref_bytes_equal_respects_n_wildcards() {
        assert!(ref_bytes_equal(b"N", b"A"));
        assert!(ref_bytes_equal(b"ACN", b"acg"));
    }

    #[test]
    fn ref_bytes_equal_rejects_real_differences() {
        assert!(!ref_bytes_equal(b"A", b"C"));
        assert!(!ref_bytes_equal(b"ACGT", b"ACGA"));
        assert!(!ref_bytes_equal(b"AC", b"ACG"));
    }

    fn make_record(info: &str) -> RawVcfRecord {
        RawVcfRecord {
            chrom: "chr1".into(),
            pos: 1,
            id: ".".into(),
            ref_allele: "A".into(),
            alt_allele: "G".into(),
            qual: ".".into(),
            filter: "PASS".into(),
            info: info.into(),
            format: None,
            samples: vec![],
        }
    }

    proptest! {
        #[test]
        fn bcftools_normalization_is_idempotent(
            reference in proptest::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
                8..40,
            ),
            alternate in proptest::collection::vec(
                prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')],
                1..8,
            ),
            position_seed in 0usize..64,
            ref_len_seed in 1usize..8,
        ) {
            let position = position_seed % reference.len() + 1;
            let ref_len = ref_len_seed.min(reference.len() - position + 1);
            let mut record = make_record(".");
            record.pos = position;
            record.ref_allele = String::from_utf8(
                reference[position - 1..position - 1 + ref_len].to_vec(),
            ).expect("generated DNA is UTF-8");
            record.alt_allele = String::from_utf8(alternate)
                .expect("generated DNA is UTF-8");
            normalize_bcftools_record(&mut record, &reference);
            let once = record.to_line();
            normalize_bcftools_record(&mut record, &reference);
            prop_assert_eq!(record.to_line(), once);
        }

        #[test]
        fn legacy_genotype_canonicalization_is_idempotent_across_ploidies_and_missing_calls(
            alleles in proptest::collection::vec(proptest::option::of(0usize..4), 1..=6),
            depth in 0usize..1000,
        ) {
            let mut record = make_record(".");
            record.format = Some("GT:DP".into());
            let genotype = alleles
                .iter()
                .map(|allele| allele.map_or_else(|| ".".to_string(), |value| value.to_string()))
                .collect::<Vec<_>>()
                .join("|");
            record.samples = vec![format!("{genotype}:{depth}")];
            canonicalize_legacy_genotypes(&mut record);
            let once = record.samples.clone();
            canonicalize_legacy_genotypes(&mut record);
            prop_assert_eq!(record.samples, once);
        }
    }

    #[test]
    fn normalized_records_restore_position_order_stably() {
        // Reduced ordering shape from chr21:44049615: a primitive emitted at
        // the trailing edge (26) preceded overlapping source records at 17.
        let mut trailing = make_record(".");
        trailing.chrom = "chr21".into();
        trailing.pos = 26;
        trailing.id = "trailing".into();
        let mut leading = trailing.clone();
        leading.pos = 15;
        leading.id = "leading".into();
        let mut overlap = trailing.clone();
        overlap.pos = 17;
        overlap.id = "overlap".into();
        let mut same_position = overlap.clone();
        same_position.id = "same-position".into();
        let mut next_contig = trailing.clone();
        next_contig.chrom = "chr1".into();
        next_contig.pos = 1;

        let mut records = vec![trailing, leading, overlap, same_position, next_contig];
        sort_normalized_records(&mut records);

        assert_eq!(
            records
                .iter()
                .map(|record| (record.chrom.as_str(), record.pos, record.id.as_str()))
                .collect::<Vec<_>>(),
            vec![
                ("chr21", 15, "leading"),
                ("chr21", 17, "overlap"),
                ("chr21", 17, "same-position"),
                ("chr21", 26, "trailing"),
                ("chr1", 1, "trailing"),
            ]
        );
    }

    #[test]
    fn blocksplit_aggregates_candidate_gaps_to_the_target_size() {
        let observations = (0..404)
            .map(|index| {
                let pos = (index / 101) * 10_000 + (index % 101) + 1;
                BlocksplitObservation {
                    chrom: "chr1".into(),
                    pos,
                    end: pos,
                    called: true,
                    location_groups: vec![0],
                }
            })
            .collect::<Vec<_>>();

        // Four groups produce three candidate gaps, but the second pass only
        // emits the middle one after cumulative candidate counts exceed the
        // 404 / 2 target. A naive every-gap reset would emit all three.
        assert_eq!(
            select_blocksplit_resets(&observations, 1, 2, None).reset_indices(),
            HashSet::from([202])
        );
    }

    #[test]
    fn blocksplit_preserves_negative_window_gap_arithmetic() {
        let observations = (0..101)
            .map(|_| BlocksplitObservation {
                chrom: "chr1".into(),
                pos: 1,
                end: 1,
                called: true,
                location_groups: vec![0],
            })
            .collect::<Vec<_>>();

        assert_eq!(
            select_blocksplit_resets(&observations, -1, 2, None).reset_indices(),
            HashSet::from([100])
        );
        assert!(
            select_blocksplit_resets(&observations, 0, 2, None)
                .reset_indices()
                .is_empty()
        );
    }

    #[test]
    fn blocksplit_ignores_uncalled_records_when_tracking_gaps() {
        let mut observations = (1..=101)
            .map(|pos| BlocksplitObservation {
                chrom: "chr1".into(),
                pos,
                end: pos,
                called: true,
                location_groups: vec![0],
            })
            .collect::<Vec<_>>();
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 10_000,
            end: 29_999,
            called: false,
            location_groups: vec![0],
        });
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 30_000,
            end: 30_000,
            called: true,
            location_groups: vec![0],
        });

        assert_eq!(
            select_blocksplit_resets(&observations, 1, 40, None).reset_indices(),
            HashSet::from([102])
        );
    }

    #[test]
    fn blocksplit_uses_effective_end_when_tracking_called_spans() {
        let mut observations = (1..=101)
            .map(|pos| BlocksplitObservation {
                chrom: "chr1".into(),
                pos,
                end: pos,
                called: true,
                location_groups: vec![0],
            })
            .collect::<Vec<_>>();
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 102,
            end: 10_000,
            called: true,
            location_groups: vec![0],
        });
        observations.push(BlocksplitObservation {
            chrom: "chr1".into(),
            pos: 10_001,
            end: 10_001,
            called: true,
            location_groups: vec![0],
        });

        assert!(
            select_blocksplit_resets(&observations, 1, 40, None)
                .reset_indices()
                .is_empty()
        );
    }

    #[test]
    fn blocksplit_partitions_explicit_locations_on_the_same_contig() {
        let observations = (0..404)
            .map(|index| {
                let location_group = index / 202;
                let index_within_location = index % 202;
                let cluster = index_within_location / 101;
                let pos = location_group * 1_000_000
                    + cluster * 10_000
                    + (index_within_location % 101)
                    + 1;
                BlocksplitObservation {
                    chrom: "chr1".into(),
                    pos,
                    end: pos,
                    called: true,
                    location_groups: vec![location_group],
                }
            })
            .collect::<Vec<_>>();

        // Legacy invokes blocksplit once per comma-separated location, so
        // each location computes its own total and target even on one contig.
        let locations = [
            vcf::LocationFilter::Contig("chr1".to_string()),
            vcf::LocationFilter::Contig("chr1".to_string()),
        ];
        assert_eq!(
            select_blocksplit_resets(&observations, 1, 2, Some(&locations)).reset_indices(),
            HashSet::from([101, 202, 303])
        );
    }

    #[test]
    fn legacy_only_parallel_overlapping_locations_duplicate_same_contig_records() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(120)))?;
        let mut vcf_text = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=120>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        )
        .to_string();
        for pos in 1..=120 {
            vcf_text.push_str(&format!(
                "chr1\t{pos}\tv{pos}\tA\tC\t30\tPASS\t.\tGT\t0/1\n"
            ));
        }
        fs::write(&input, vcf_text)?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.locations = Some("chr1:1-100,chr1:51-120".to_string());
        args.threads = Some(2);
        run(args)?;

        let positions = vcf::load_raw_vcf(&output)?
            .1
            .into_iter()
            .map(|record| record.pos)
            .collect::<Vec<_>>();
        assert_eq!(positions.len(), 240);
        for (offset, pair) in positions.chunks_exact(2).enumerate() {
            assert_eq!(pair, [offset + 1, offset + 1]);
        }
        Ok(())
    }

    #[test]
    fn parallel_location_streams_merge_equal_normalized_positions_roundwise() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(20)))?;
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=20>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
                "chr1\t3\tbarrier\tA\tC\t30\tPASS\t.\tGT\t0/1\n",
                "chr1\t5\trepeat_deletion\tAA\tA\t30\tPASS\t.\tGT\t0/1\n",
            ),
        )?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.locations = Some("chr1:1-10,chr1:3-20".to_string());
        args.threads = Some(2);
        args.window_size = 1;
        run(args)?;

        let alleles = vcf::load_raw_vcf(&output)?
            .1
            .into_iter()
            .map(|record| (record.pos, record.ref_allele, record.alt_allele))
            .collect::<Vec<_>>();
        assert_eq!(
            alleles,
            [
                (3, "A".to_string(), "C".to_string()),
                (3, "A".to_string(), "C".to_string()),
                (3, "AA".to_string(), "A".to_string()),
                (3, "AA".to_string(), "A".to_string()),
            ]
        );
        Ok(())
    }

    #[test]
    fn preprocess_stream_tie_break_state_respects_left_shift_window() -> Result<()> {
        let mut spool = PreprocessSpool::new(true, 1, &[])?;
        for pos in 1..=4_096 {
            let mut record = make_record(".");
            record.pos = pos;
            spool.push(
                vcf::ValidatedVcfRecord::try_from_raw(record, QueryProvenance::Unavailable)?,
                0,
            )?;
        }

        assert!(spool.retained_position_count() <= LEFT_SHIFT_WINDOW + 1);
        Ok(())
    }

    #[test]
    fn preprocess_spool_sorts_numeric_contigs_by_declared_header_order() -> Result<()> {
        let declared = ["1".to_string(), "2".to_string(), "10".to_string()];
        let mut spool = PreprocessSpool::new(true, 1, &declared)?;
        for chrom in ["10", "2", "1"] {
            let mut record = make_record(".");
            record.chrom = chrom.to_string();
            spool.push(
                vcf::ValidatedVcfRecord::try_from_raw(record, QueryProvenance::Unavailable)?,
                0,
            )?;
        }

        let observed = spool
            .finish()?
            .map(|record| record.map(|record| record.raw().chrom.clone()))
            .collect::<Result<Vec<_>>>()?;
        assert_eq!(observed, ["1", "2", "10"]);
        Ok(())
    }

    #[test]
    fn normative_set_union_multiblock_keeps_records_selected_only_by_later_range() {
        let locations = [
            vcf::LocationFilter::Range {
                chrom: "chr1".into(),
                start: 1,
                end: 101,
            },
            vcf::LocationFilter::Range {
                chrom: "chr1".into(),
                start: 10_000,
                end: 20_101,
            },
        ];
        let observations = (0..303)
            .map(|index| {
                let pos = match index {
                    0..101 => index + 1,
                    101..202 => 10_000 + index - 101,
                    _ => 20_000 + index - 202,
                };
                BlocksplitObservation {
                    chrom: "chr1".into(),
                    pos,
                    end: pos,
                    called: true,
                    location_groups: compatibility::location_stream_groups(
                        compatibility::LocationStreamPolicy::SetUnion,
                        &locations,
                        "chr1",
                        pos,
                    ),
                }
            })
            .collect::<Vec<_>>();

        let selection = select_blocksplit_resets_with_policy(
            &observations,
            1,
            40,
            Some(&locations),
            compatibility::LocationStreamPolicy::SetUnion,
        );

        assert_eq!(
            selection.included_indices(),
            Some((0..303).collect::<HashSet<_>>())
        );
        assert!(
            selection
                .included_indices()
                .is_some_and(|indices| indices.contains(&302)),
            "the final record exists only in the later range"
        );
    }

    #[test]
    fn parallel_comma_locations_preserve_independent_multiblock_jobs() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let output = directory.path().join("output.vcf.gz");
        let reference = directory.path().join("ref.fa");
        fs::write(&reference, format!(">chr1\n{}\n", "A".repeat(3_200)))?;
        let mut vcf_text = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=3200>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        )
        .to_string();
        for cluster_start in [1, 1_001, 2_001, 3_001] {
            for pos in cluster_start..=cluster_start + 100 {
                vcf_text.push_str(&format!(
                    "chr1\t{pos}\tv{pos}\tA\tC\t30\tPASS\t.\tGT\t0/1\n"
                ));
            }
        }
        fs::write(&input, vcf_text)?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.gender = PreprocessGender::None;
        args.locations = Some("chr1:1-2101,chr1:1001-3101".to_string());
        args.threads = Some(2);
        args.window_size = 10;
        run(args)?;

        let positions = vcf::load_raw_vcf(&output)?
            .1
            .into_iter()
            .map(|record| record.pos)
            .collect::<Vec<_>>();
        assert_eq!(positions.len(), 705);
        assert!(positions.windows(2).all(|pair| pair[0] <= pair[1]));

        let mut multiplicities = std::collections::BTreeMap::new();
        for pos in positions {
            *multiplicities.entry(pos).or_insert(0usize) += 1;
        }
        for pos in 1..=101 {
            assert_eq!(multiplicities.get(&pos), Some(&2));
        }
        for pos in 1_001..=1_101 {
            assert_eq!(multiplicities.get(&pos), Some(&2));
        }
        for pos in 2_001..=2_100 {
            assert_eq!(multiplicities.get(&pos), Some(&2));
        }
        assert_eq!(multiplicities.get(&2_101), Some(&1));
        for pos in 3_001..=3_100 {
            assert_eq!(multiplicities.get(&pos), Some(&1));
        }
        assert_eq!(multiplicities.get(&3_101), None);
        Ok(())
    }

    #[test]
    fn blocksplit_filters_inactive_partitions_and_preserves_all_empty_fallback() {
        let mut observations = vec![
            BlocksplitObservation {
                chrom: "chr1".into(),
                pos: 1,
                end: 1,
                called: true,
                location_groups: vec![0],
            },
            BlocksplitObservation {
                chrom: "chr1".into(),
                pos: 2,
                end: 2,
                called: false,
                location_groups: vec![0],
            },
            BlocksplitObservation {
                chrom: "chr2".into(),
                pos: 1,
                end: 1,
                called: false,
                location_groups: vec![0],
            },
        ];

        let mixed = select_blocksplit_resets(&observations, 1, 2, None);
        assert_eq!(mixed.included_indices(), Some(HashSet::from([0, 1])));

        observations[0].called = false;
        let all_empty = select_blocksplit_resets(&observations, 1, 2, None);
        assert_eq!(all_empty.included_indices(), None);
    }

    #[test]
    fn blocksplit_parallelism_resolves_explicit_and_default_thread_counts() {
        assert_eq!(effective_thread_count_with_available(Some(0), 8), 0);
        assert_eq!(effective_thread_count_with_available(Some(1), 8), 1);
        assert_eq!(effective_thread_count_with_available(Some(2), 1), 2);
        assert_eq!(effective_thread_count_with_available(None, 0), 1);
        assert_eq!(effective_thread_count_with_available(None, 1), 1);
        assert_eq!(effective_thread_count_with_available(None, 2), 2);
    }

    #[test]
    fn strip_stale_info_keys_removes_allele_counts() {
        let mut record = make_record("AC=1;AF=0.5;AN=2;DP=158");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, "AF=0.5;DP=158");
    }

    #[test]
    fn strip_stale_info_keys_handles_missing_info() {
        let mut record = make_record(".");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, ".");
    }

    #[test]
    fn strip_stale_info_keys_collapses_to_dot_when_all_stripped() {
        let mut record = make_record("AC=1;AN=2;MLEAC=1;MLEAF=0.5");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, ".");
    }

    #[test]
    fn strip_stale_info_keys_preserves_flags_without_values() {
        let mut record = make_record("SOMATIC;AC=1;DP=158");
        strip_stale_info_keys(&mut record);
        assert_eq!(record.info, "SOMATIC;DP=158");
    }

    #[test]
    fn sort_info_keys_orders_alphabetically() {
        let mut record = make_record("HRun=0;AF=0.5;DP=158;Dels=0");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "AF=0.5;DP=158;Dels=0;HRun=0");
    }

    #[test]
    fn sort_info_keys_handles_flag_entries() {
        let mut record = make_record("SOMATIC;DP=158;AF=0.5");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "AF=0.5;DP=158;SOMATIC");
    }

    #[test]
    fn sort_info_keys_collapses_negative_zero_to_zero() {
        // Legacy hap.py reads INFO floats through htslib which loses the
        // sign on negative zero before re-emitting. Match that here.
        let mut record = make_record("SB=-0;FS=12.24;MQRankSum=0");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "FS=12.24;MQRankSum=0;SB=0");
    }

    #[test]
    fn sort_info_keys_preserves_non_zero_negatives() {
        let mut record = make_record("SB=-1.5;BaseQRankSum=-0.123");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "BaseQRankSum=-0.123;SB=-1.5");
    }

    #[test]
    fn expand_haploid_gt_matches_legacy_output() {
        // Haploid alt calls expand to heterozygous half-call, mirroring
        // VariantAlleleSplitter.cpp:180-227. Empirically confirmed: legacy
        // result.vcf.gz shows query GT=1 on autosomes as 0/1:het, not 1/1.
        assert_eq!(expand_haploid_gt("1", false), "0/1");
        assert_eq!(expand_haploid_gt("2", false), "0/2");
        assert_eq!(expand_haploid_gt("1", true), "1/1");
        assert_eq!(expand_haploid_gt("2", true), "2/2");
        assert_eq!(expand_haploid_gt("0", true), "0/0");
        assert_eq!(expand_haploid_gt(".", true), "./.");
        // Diploid inputs pass through.
        assert_eq!(expand_haploid_gt("0/1", true), "0/1");
        assert_eq!(expand_haploid_gt("1|2", true), "1|2");
        assert_eq!(expand_haploid_gt("./.", true), "./.");
    }

    #[test]
    fn active_preprocessing_masks_genotypes_wider_than_diploid() {
        let mut record = make_record(".");
        record.format = Some("GT:GQ".to_string());
        record.samples = vec![
            "0/1:40".to_string(),
            "0/1/1:50".to_string(),
            "0|0|1|1:60".to_string(),
        ];

        mask_genotypes_wider_than_diploid(&mut record);

        assert_eq!(record.samples, ["0/1:40", ".:50", ".:60"]);
    }

    #[test]
    fn sort_info_keys_collapses_negative_zero_in_lists() {
        let mut record = make_record("AF=-0,0.5");
        sort_info_keys(&mut record);
        assert_eq!(record.info, "AF=0,0.5");
    }

    #[test]
    fn collapse_pl_reduces_to_last_value() {
        let mut record = make_record(".");
        record.format = Some("GT:AD:PL".to_string());
        record.samples = vec!["0/1:32,126:3279,0,716".to_string()];
        collapse_pl_to_last_value(&mut record);
        assert_eq!(record.samples[0], "0/1:32,126:716");
    }

    #[test]
    fn collapse_pl_leaves_scalar_alone() {
        let mut record = make_record(".");
        record.format = Some("GT:PL".to_string());
        record.samples = vec!["0/1:716".to_string()];
        collapse_pl_to_last_value(&mut record);
        assert_eq!(record.samples[0], "0/1:716");
    }

    #[test]
    fn ado_zero_when_genotype_covers_all_alleles() {
        assert_eq!(compute_ado("32,126", Some("0/1")), 0);
        assert_eq!(compute_ado("32,126", Some("0|1")), 0);
    }

    #[test]
    fn ado_captures_ref_depth_for_hom_alt() {
        // GT=1/1 with AD=22,413 → ADO = AD[0] = 22 (the unused reference depth)
        assert_eq!(compute_ado("22,413", Some("1/1")), 22);
        assert_eq!(compute_ado("1,17", Some("1/1")), 1);
    }

    #[test]
    fn ado_captures_alt_depth_for_hom_ref() {
        assert_eq!(compute_ado("100,25", Some("0/0")), 25);
    }

    #[test]
    fn ado_handles_multi_allelic() {
        // GT=1/2 across original alleles — both indices 1 and 2 are called,
        // so AD[0] (ref) is the only "other" depth.
        assert_eq!(compute_ado("10,50,70", Some("1/2")), 10);
        // GT=0/2 — index 1 is the "other" allele.
        assert_eq!(compute_ado("10,50,70", Some("0/2")), 50);
    }

    #[test]
    fn ado_defaults_to_zero_on_missing_gt_or_ad() {
        assert_eq!(compute_ado(".", Some("0/1")), 0);
        assert_eq!(compute_ado("10,20", Some("./.")), 0);
        assert_eq!(compute_ado("10,20", None), 0);
    }

    #[test]
    fn classify_value_distinguishes_int_float_string() {
        assert_eq!(classify_value("99"), ScalarType::Integer);
        assert_eq!(classify_value("3279,0,716"), ScalarType::Integer);
        assert_eq!(classify_value("95.77"), ScalarType::Float);
        assert_eq!(classify_value("0.934"), ScalarType::Float);
        assert_eq!(classify_value("PASS"), ScalarType::String);
        assert_eq!(classify_value("."), ScalarType::Integer);
    }

    #[test]
    fn reorder_format_keeps_gt_ad_ado_dp_first_then_buckets_by_type() {
        // Record with GQ=95.77 (float) should push GQ after the integer
        // bucket — legacy canonical shape `GT:AD:ADO:DP:GQX:MQ:PL:GQ:VF`.
        let mut record = make_record(".");
        record.format = Some("GT:AD:ADO:DP:GQ:GQX:MQ:PL:VF".to_string());
        record.samples = vec!["0/1:22,306:0:328:95.77:96:36:96:0.933".to_string()];
        reorder_format_fields(&mut record);
        assert_eq!(
            record.format.as_deref().unwrap(),
            "GT:AD:ADO:DP:GQX:MQ:PL:GQ:VF"
        );
        assert_eq!(record.samples[0], "0/1:22,306:0:328:96:36:96:95.77:0.933");
    }

    #[test]
    fn reorder_format_treats_integer_gq_as_integer() {
        // When the sample's GQ happens to parse as an int (e.g. 99), legacy
        // places it with the ints — alphabetical means GQ < GQX.
        let mut record = make_record(".");
        record.format = Some("GT:AD:ADO:DP:GQ:GQX:MQ:PL:VF".to_string());
        record.samples = vec!["1/1:0,101:0:101:99:99:51:0:1".to_string()];
        reorder_format_fields(&mut record);
        assert_eq!(
            record.format.as_deref().unwrap(),
            "GT:AD:ADO:DP:GQ:GQX:MQ:PL:VF"
        );
    }

    #[test]
    fn somatic_conversion_uses_legacy_gt_modes_and_preserves_sample_formats_in_info() {
        let mut record = make_record("SOMATIC");
        record.alt_allele = "C,G".to_string();
        record.format = Some("GT:DP:VF".to_string());
        record.samples = vec!["0/1:12:0.25".to_string(), "1/2:30:0.75".to_string()];
        let sample_names = vec!["NORMAL".to_string(), "TUMOR".to_string()];

        let half = convert_somatic_record(&record, SomaticGtMode::Half, &sample_names);
        assert_eq!(half.len(), 2);
        assert_eq!(half[0].alt_allele, "C");
        assert_eq!(half[1].alt_allele, "G");
        assert_eq!(half[0].format.as_deref(), Some("GT:AD:DP"));
        assert_eq!(half[0].samples, ["./1:.,.:0"]);
        assert!(half[0].info.contains("NORMAL_GT=2,4"));
        assert!(half[0].info.contains("NORMAL_DP=12"));
        assert!(half[0].info.contains("NORMAL_VF=0.25"));
        assert!(!half[0].info.contains("TUMOR_GT="));
        assert!(!half[0].info.contains("TUMOR_DP=30"));
        assert!(!half[0].info.contains("TUMOR_VF=0.75"));

        let hemi = convert_somatic_record(&record, SomaticGtMode::Hemi, &sample_names);
        assert_eq!(hemi[0].samples, ["1:.,.:0"]);
        let het = convert_somatic_record(&record, SomaticGtMode::Het, &sample_names);
        assert_eq!(het[0].samples, ["0/1:.,.:0"]);
        let hom = convert_somatic_record(&record, SomaticGtMode::Hom, &sample_names);
        assert_eq!(hom[0].samples, ["1/1:.,.:0"]);
    }

    #[test]
    fn somatic_first_mode_preserves_first_gt_and_does_not_split_alts() {
        let mut record = make_record(".");
        record.alt_allele = "C,G".to_string();
        record.format = Some("GT:DP".to_string());
        record.samples = vec!["2|1:12".to_string(), "0/2:30".to_string()];
        let sample_names = vec!["NORMAL".to_string(), "TUMOR".to_string()];

        let converted = convert_somatic_record(&record, SomaticGtMode::First, &sample_names);
        assert_eq!(converted.len(), 1);
        assert_eq!(converted[0].alt_allele, "C,G");
        assert_eq!(converted[0].format.as_deref(), Some("GT:AD:DP"));
        assert_eq!(converted[0].samples, ["2|1:.,.,.:0"]);
        assert!(converted[0].info.contains("NORMAL_GT=6,5"));
        assert!(!converted[0].info.contains("TUMOR_GT="));
    }

    #[test]
    fn somatic_partial_credit_merges_nonhom_siblings_but_not_hom_calls() {
        let mut record = make_record("TAG=multi");
        record.alt_allele = "T,G".to_string();
        record.format = Some("GT:DP".to_string());
        record.samples = vec!["1/2:30".to_string()];
        let names = vec!["TUMOR".to_string()];

        let half = finalize_somatic_records(
            convert_somatic_record(&record, SomaticGtMode::Half, &names),
            SomaticGtMode::Half,
        );
        assert_eq!(half.len(), 1);
        assert_eq!(half[0].alt_allele, "G,T");
        assert_eq!(half[0].samples, ["2/1:.,.,.:0"]);

        let hom = finalize_somatic_records(
            convert_somatic_record(&record, SomaticGtMode::Hom, &names),
            SomaticGtMode::Hom,
        );
        assert_eq!(hom.len(), 2);
        assert_eq!(hom[0].samples, ["1/1:.,.:0"]);
    }

    #[test]
    fn somatic_half_calls_bypass_partial_credit_for_scmp() {
        let mut record = make_record(".");
        record.alt_allele = "C,G".to_string();
        record.format = Some("GT".to_string());
        record.samples = vec!["1/2".to_string()];
        let names = vec!["TUMOR".to_string()];
        let converted = convert_somatic_record(&record, SomaticGtMode::Half, &names);

        let scmp = finalize_somatic_for_pipeline(converted.clone(), SomaticGtMode::Half, false);
        assert_eq!(scmp.len(), 2);
        assert_eq!(scmp[0].format.as_deref(), Some("GT"));
        assert_eq!(scmp[0].samples, ["./1"]);
        assert_eq!(scmp[1].samples, ["./1"]);

        let normalized = finalize_somatic_for_pipeline(converted, SomaticGtMode::Half, true);
        assert_eq!(normalized.len(), 1);
        assert_eq!(normalized[0].samples, ["2/1:.,.,.:0"]);
    }

    #[test]
    fn somatic_headers_declare_sample_prefixed_format_info_fields() {
        let mut headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\tNORMAL".to_string(),
        ];
        let names = somatic_info_sample_names(&headers);
        assert_eq!(names, ["TUMOR", "NORMAL"]);
        append_somatic_info_headers(&mut headers, &names);

        assert!(headers.contains(
            &"##INFO=<ID=TUMOR_GT,Number=1,Type=String,Description=\"Genotype\">".to_string()
        ));
        assert!(headers.contains(
            &"##INFO=<ID=NORMAL_DP,Number=1,Type=Integer,Description=\"Depth\">".to_string()
        ));
    }

    #[test]
    fn somatic_single_input_sample_uses_output_sample_prefix() {
        let headers =
            vec!["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tORIGINAL".to_string()];
        assert_eq!(somatic_info_sample_names(&headers), ["SAMPLE"]);
    }

    #[test]
    fn non_ref_filter_drops_only_records_whose_gt_calls_the_non_ref_allele() {
        let mut record = make_record(".");
        record.alt_allele = "C,<NON_REF>".to_string();
        record.format = Some("GT:DP".to_string());
        record.samples = vec!["0/1:12".to_string(), "0/0:30".to_string()];
        assert!(!calls_non_ref_allele(&record));

        record.samples[1] = "0/2:30".to_string();
        assert!(calls_non_ref_allele(&record));
    }

    #[test]
    fn uncalled_non_ref_is_trimmed_and_reference_blocks_are_dropped() {
        let mut variant = make_record(".");
        variant.alt_allele = "C,<NON_REF>".to_string();
        variant.format = Some("GT:AD".to_string());
        variant.samples = vec!["0/1:7,8,0".to_string()];
        assert!(trim_uncalled_non_ref(&mut variant));
        assert_eq!(variant.alt_allele, "C");
        assert_eq!(variant.samples, ["0/1:7,8"]);

        let mut block = make_record("END=10");
        block.alt_allele = "<NON_REF>".to_string();
        block.format = Some("GT".to_string());
        block.samples = vec!["0/0".to_string()];
        assert!(!trim_uncalled_non_ref(&mut block));
    }

    #[test]
    fn called_alt_projection_remaps_gt_and_ad_and_drops_homref_records() {
        let mut called = make_record(".");
        called.alt_allele = "C,G,T".to_string();
        called.format = Some("GT:AD:ADO".to_string());
        called.samples = vec!["0|2:10,1,8,2:3".to_string()];
        assert!(retain_called_alternates(&mut called));
        assert_eq!(called.alt_allele, "G");
        assert_eq!(called.samples, ["0|1:10,8:3"]);

        let mut homref = make_record(".");
        homref.alt_allele = "C,G".to_string();
        homref.format = Some("GT:AD".to_string());
        homref.samples = vec!["0/0:12,0,0".to_string()];
        assert!(!retain_called_alternates(&mut homref));

        let mut spanning_deletion = make_record(".");
        spanning_deletion.alt_allele = "T,*".to_string();
        spanning_deletion.format = Some("GT:AD".to_string());
        spanning_deletion.samples = vec!["1/2:8,5,6".to_string()];
        assert!(retain_called_alternates(&mut spanning_deletion));
        assert_eq!(spanning_deletion.alt_allele, "T");
        assert_eq!(spanning_deletion.samples, ["0/1:8,5"]);
    }

    #[test]
    fn preprocessing_removes_uncalled_multi_alleles_before_primitive_split() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("input.vcf");
        let reference = directory.path().join("ref.fa");
        let output = directory.path().join("output.vcf.gz");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "##contig=<ID=chr1,length=20>\n",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "chr1\t2\t.\tA\tATC,ATCTC\t.\tPASS\t.\tGT:AD\t2|0:5,0,7\n",
                "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\tGT:AD\t0/0:9,0,0\n",
            ),
        )?;
        fs::write(&reference, ">chr1\nAAAAAAAAAAAAAAAAAAAA\n")?;

        let mut args = interval_args(&input, &output, &reference, None, None);
        args.decompose = true;
        args.window_size = 4096;
        run(args)?;

        let (_, records) = vcf::load_raw_vcf(&output)?;
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].pos, 2);
        assert_eq!(records[0].ref_allele, "A");
        assert_eq!(records[0].alt_allele, "ATCTC");
        assert_eq!(records[0].qual, "0");
        assert_eq!(records[0].samples, ["0/1:5,7:0:0"]);
        Ok(())
    }

    #[test]
    fn legacy_writer_unphases_and_canonicalizes_genotypes() {
        let mut record = make_record(".");
        record.format = Some("GT".to_string());
        record.samples = vec![
            "1|1".to_string(),
            "1|0".to_string(),
            "0|1".to_string(),
            "1|2".to_string(),
            "2/1".to_string(),
        ];
        canonicalize_legacy_genotypes(&mut record);
        assert_eq!(record.samples, ["1/1", "0/1", "0/1", "2/1", "2/1"]);
    }

    #[test]
    fn multi_allelic_order_remaps_gt_and_ad_like_legacy_aggregator() {
        let mut record = make_record(".");
        record.alt_allele = "T,G".to_string();
        record.format = Some("GT:AD".to_string());
        record.samples = vec!["2/1:10,11,9".to_string(), "0/1:14,14,0".to_string()];
        canonicalize_multi_allelic_order(&mut record);
        assert_eq!(record.alt_allele, "G,T");
        assert_eq!(record.samples[0], "2/1:10,9,11");
        assert_eq!(record.samples[1], "0/2:14,0,14");
    }

    #[test]
    fn secondary_samples_keep_core_fields_but_not_dynamic_annotations() {
        let mut record = make_record(".");
        record.format = Some("GT:AD:ADO:DP:GQ:TXT".to_string());
        record.samples = vec![
            "0/1:7,8:0:15:20:variant".to_string(),
            "0/0:16,0:0:16:25:normal".to_string(),
        ];
        let string_fields = BTreeSet::from(["GT".to_string(), "TXT".to_string()]);
        blank_secondary_sample_annotations(&mut record, false, &string_fields);
        assert_eq!(record.samples[0], "0/1:7,8:0:15:20:variant");
        assert_eq!(record.samples[1], "0/0:.,.:.:16:.:.");

        record.samples[1] = "0/0:16,0:0:16:25:normal".to_string();
        blank_secondary_sample_annotations(&mut record, true, &string_fields);
        assert_eq!(record.samples[1], "0/0:.,.:.:16:.:");
    }

    #[test]
    fn gvcf_conversion_trims_non_ref_and_unrequested_fields() {
        let mut record = make_record("END=10;AC=1;DP=40");
        record.alt_allele = "C,<NON_REF>".to_string();
        record.format = Some("GT:DP:GQ:AD:PL".to_string());
        record.samples = vec!["0/1:12:40:8,4,0:0,10,100,20,200,300".to_string()];

        assert!(convert_gvcf_record(&mut record));
        assert_eq!(record.alt_allele, "C");
        assert_eq!(record.info, ".");
        assert_eq!(record.format.as_deref(), Some("GT:DP:GQ"));
        assert_eq!(record.samples, ["0/1:12:40"]);
    }

    #[test]
    fn gvcf_conversion_strips_unrequested_info_and_format_headers() {
        let mut headers = vec![
            "##INFO=<ID=TAG,Number=1,Type=String,Description=\"Tag\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">".to_string(),
            "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"GQ\">".to_string(),
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">".to_string(),
            "##FORMAT=<ID=TXT,Number=1,Type=String,Description=\"TXT\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS".to_string(),
        ];
        filter_gvcf_headers(&mut headers);
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=GT,"))
        );
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=DP,"))
        );
        assert!(
            headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=GQ,"))
        );
        assert!(!headers.iter().any(|line| line.starts_with("##INFO=")));
        assert!(
            !headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=AD,"))
        );
        assert!(
            !headers
                .iter()
                .any(|line| line.starts_with("##FORMAT=<ID=TXT,"))
        );
    }

    #[test]
    fn gvcf_conversion_drops_single_alt_blocks_and_non_ref_calls() {
        let mut reference_block = make_record("END=10");
        reference_block.alt_allele = "<NON_REF>".to_string();
        reference_block.format = Some("GT:DP:GQ".to_string());
        reference_block.samples = vec!["0/0:12:40".to_string()];
        assert!(!convert_gvcf_record(&mut reference_block));

        let mut non_ref_call = make_record(".");
        non_ref_call.alt_allele = "C,<NON_REF>".to_string();
        non_ref_call.format = Some("GT:DP:GQ".to_string());
        non_ref_call.samples = vec!["0/2:12:40".to_string()];
        assert!(!convert_gvcf_record(&mut non_ref_call));
    }

    #[test]
    fn gvcf_conversion_retains_initial_multi_alt_homref_as_reference_record() {
        let mut record = make_record("END=10");
        record.alt_allele = "C,<NON_REF>".to_string();
        record.format = Some("GT:DP:GQ".to_string());
        record.samples = vec!["0/0:12:40".to_string()];

        assert!(convert_gvcf_record(&mut record));
        assert_eq!(record.alt_allele, ".");
        assert_eq!(record.info, ".");
        assert_eq!(record.samples, ["0/0:12:40"]);
    }

    #[test]
    fn parallel_blocksplit_omits_empty_partitions_but_falls_back_when_all_are_empty() -> Result<()>
    {
        let directory = tempdir()?;
        let reference = directory.path().join("ref.fa");
        let mixed_input = directory.path().join("mixed.vcf");
        let mixed_output = directory.path().join("mixed.out.vcf.gz");
        let empty_input = directory.path().join("empty.vcf");
        let empty_output = directory.path().join("empty.out.vcf.gz");
        fs::write(&reference, ">chr1\nAAAAA\n>chr2\nAAAAA\n")?;
        let header = concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=5>\n",
            "##contig=<ID=chr2,length=5>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n",
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n",
        );
        fs::write(
            &mixed_input,
            format!(
                "{header}{}{}",
                "chr1\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/1:12:40\n",
                "chr2\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/0:12:40\n",
            ),
        )?;
        fs::write(
            &empty_input,
            format!(
                "{header}{}{}",
                "chr1\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/0:12:40\n",
                "chr2\t2\t.\tA\tC,<NON_REF>\t30\tPASS\t.\tGT:DP:GQ\t0/0:12:40\n",
            ),
        )?;

        let mut mixed_args = interval_args(&mixed_input, &mixed_output, &reference, None, None);
        mixed_args.gender = PreprocessGender::None;
        mixed_args.threads = Some(2);
        mixed_args.convert_gvcf_to_vcf = true;
        run(mixed_args)?;
        let mixed_records = vcf::load_raw_vcf(&mixed_output)?.1;
        assert_eq!(
            mixed_records
                .iter()
                .map(|record| record.chrom.as_str())
                .collect::<Vec<_>>(),
            ["chr1"]
        );

        let mut empty_args = interval_args(&empty_input, &empty_output, &reference, None, None);
        empty_args.gender = PreprocessGender::None;
        empty_args.threads = Some(2);
        empty_args.convert_gvcf_to_vcf = true;
        run(empty_args)?;
        let empty_records = vcf::load_raw_vcf(&empty_output)?.1;
        assert_eq!(
            empty_records
                .iter()
                .map(|record| (record.chrom.as_str(), record.alt_allele.as_str()))
                .collect::<Vec<_>>(),
            [("chr1", "."), ("chr2", ".")]
        );
        Ok(())
    }
}
