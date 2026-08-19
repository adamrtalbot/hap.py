use flate2::read::MultiGzDecoder;
use serde_json::Value;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use tempfile::tempdir;

const SUBCOMMANDS: [&str; 6] = ["germline", "somatic", "pre", "ftx", "quantify", "validate"];

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn hap() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_hap"));
    command.current_dir(root());
    command
}

fn run(arguments: &[&str]) -> Output {
    hap().args(arguments).output().expect("run hap binary")
}

fn run_strings(arguments: &[String]) -> Output {
    hap().args(arguments).output().expect("run hap binary")
}

fn assert_success(output: &Output, context: &str) {
    assert!(
        output.status.success(),
        "{context} failed with {:?}\nstdout:\n{}\nstderr:\n{}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

fn assert_runtime_failure(output: &Output, context: &str) {
    assert_eq!(output.status.code(), Some(1), "{context}");
    assert!(output.stdout.is_empty(), "{context} stdout");
    assert!(
        String::from_utf8_lossy(&output.stderr).contains("Error:"),
        "{context} stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn path(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

fn read_gzip(path: &Path) -> Vec<u8> {
    let mut decoder = MultiGzDecoder::new(fs::File::open(path).expect("open gzip artifact"));
    let mut bytes = Vec::new();
    decoder
        .read_to_end(&mut bytes)
        .expect("decode gzip artifact");
    bytes
}

fn parse_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("read JSON artifact")).expect("valid JSON")
}

fn parse_gzip_json(path: &Path) -> Value {
    serde_json::from_slice(&read_gzip(path)).expect("valid gzipped JSON")
}

fn assert_metrics_document(document: &Value) {
    let root = document.as_object().expect("metrics root is an object");
    let run_info = root
        .get("runInfo")
        .and_then(Value::as_array)
        .expect("metrics runInfo is an array");
    assert!(run_info.iter().any(|entry| {
        entry.get("key") == Some(&Value::String("commandline".into()))
            && entry.get("value").is_some_and(Value::is_string)
    }));
    let metrics = root
        .get("metrics")
        .and_then(Value::as_array)
        .expect("metrics is an array");
    assert!(
        !metrics.is_empty(),
        "at least one metrics table is required"
    );
    for table in metrics {
        let table = table.as_object().expect("metrics table is an object");
        assert!(table.get("id").is_some_and(Value::is_string));
        assert!(table.get("data").is_some_and(Value::is_array));
    }
}

#[test]
fn help_version_and_failure_stream_contracts_are_pinned() {
    let help = run(&["--help"]);
    assert_eq!(help.status.code(), Some(0));
    assert!(help.stderr.is_empty());
    let help = String::from_utf8(help.stdout).expect("help is UTF-8");
    assert!(
        help.starts_with("Single-binary Rust haplotype comparison tool\n\nUsage: hap <COMMAND>")
    );
    for subcommand in SUBCOMMANDS {
        assert!(
            help.contains(&format!("  {subcommand}")),
            "missing {subcommand}"
        );
        let output = run(&[subcommand, "--help"]);
        assert_eq!(output.status.code(), Some(0), "{subcommand} help");
        assert!(output.stderr.is_empty(), "{subcommand} help stderr");
        assert!(
            String::from_utf8_lossy(&output.stdout).contains(&format!("Usage: hap {subcommand}")),
            "{subcommand} help usage"
        );
    }

    let version = run(&["--version"]);
    assert_eq!(version.status.code(), Some(0));
    assert_eq!(version.stdout, b"hap 0.1.0\n");
    assert!(version.stderr.is_empty());

    for subcommand in SUBCOMMANDS {
        let output = run(&[subcommand]);
        assert_eq!(output.status.code(), Some(2), "{subcommand}");
        assert!(output.stdout.is_empty(), "{subcommand} failure stdout");
        assert!(
            String::from_utf8_lossy(&output.stderr)
                .starts_with("error: the following required arguments were not provided:"),
            "{subcommand} failure stderr"
        );
    }
}

#[test]
fn legacy_alias_version_delimiter_and_unknown_option_contracts_are_pinned() {
    for alias in ["compare", "preprocess", "prepy", "qfy", "ftxpy", "vcfcheck"] {
        let output = run(&[alias, "--help"]);
        assert_eq!(output.status.code(), Some(0), "{alias} help");
        assert!(output.stderr.is_empty(), "{alias} help stderr");
    }

    for (arguments, expected) in [
        (&["germline", "-v"][..], "Hap.py \n"),
        (&["compare", "-v"][..], "Hap.py \n"),
        (&["pre", "-v", "in", "out"][..], "pre.py \n"),
        (&["preprocess", "-v", "in", "out"][..], "pre.py \n"),
        (&["prepy", "-v", "in", "out"][..], "pre.py \n"),
        (
            &["quantify", "-v", "in", "-o", "out", "-r", "ref"][..],
            "qfy.py \n",
        ),
        (
            &["qfy", "-v", "in", "-o", "out", "-r", "ref"][..],
            "qfy.py \n",
        ),
    ] {
        let output = run(arguments);
        assert_eq!(output.status.code(), Some(0), "{arguments:?}");
        assert_eq!(String::from_utf8_lossy(&output.stdout), expected);
        assert!(output.stderr.is_empty(), "{arguments:?}");
    }

    for subcommand in [
        "germline",
        "compare",
        "somatic",
        "pre",
        "preprocess",
        "prepy",
        "ftx",
        "ftxpy",
        "quantify",
        "qfy",
        "validate",
        "vcfcheck",
    ] {
        let output = run(&[subcommand, "--unknown-option"]);
        assert_eq!(output.status.code(), Some(2), "{subcommand}");
        assert!(String::from_utf8_lossy(&output.stderr).contains("unexpected argument"));
    }

    for alias in ["pre", "preprocess", "prepy"] {
        let delimiter = run(&[alias, "--", "--version", "out.vcf"]);
        assert_eq!(delimiter.status.code(), Some(1), "{alias}");
        assert!(delimiter.stdout.is_empty(), "{alias}");
        assert!(
            String::from_utf8_lossy(&delimiter.stderr).contains("Error:"),
            "{alias}"
        );
    }

    for alias in ["germline", "compare"] {
        let delimiter = run(&[alias, "--", "--version", "query.vcf"]);
        assert_eq!(delimiter.status.code(), Some(2), "{alias}");
        assert!(delimiter.stdout.is_empty(), "{alias}");
        assert!(
            String::from_utf8_lossy(&delimiter.stderr).contains("required arguments"),
            "{alias}"
        );
    }

    for alias in ["quantify", "qfy"] {
        let delimiter = run(&[alias, "--", "--version"]);
        assert_eq!(delimiter.status.code(), Some(2), "{alias}");
        assert!(delimiter.stdout.is_empty(), "{alias}");
        assert!(
            String::from_utf8_lossy(&delimiter.stderr).contains("required arguments"),
            "{alias}"
        );
    }
}

#[test]
fn every_public_subcommand_has_a_runtime_handler_failure_contract() {
    let directory = tempdir().expect("temporary output directory");
    let output = path(&directory.path().join("result"));
    let reference = "tests/fixtures/synth-snp-match/ref.fa";
    let missing = path(&directory.path().join("does-not-exist.vcf"));

    for (subcommand, arguments) in [
        (
            "germline",
            vec![
                "germline", &missing, &missing, "-r", reference, "-o", &output,
            ],
        ),
        (
            "somatic",
            vec![
                "somatic", &missing, &missing, "-r", reference, "-o", &output,
            ],
        ),
        ("pre", vec!["pre", &missing, &output, "-r", reference]),
        ("ftx", vec!["ftx", &missing, "-o", &output, "-r", reference]),
        (
            "quantify",
            vec!["quantify", &missing, "-o", &output, "-r", reference],
        ),
        ("validate", vec!["validate", &missing]),
    ] {
        assert_runtime_failure(&run(&arguments), subcommand);
    }
}

#[test]
fn every_public_subcommand_succeeds_as_a_black_box_and_json_is_structural() {
    let directory = tempdir().expect("temporary output directory");
    let output_root = directory.path();

    let germline_prefix = output_root.join("germline");
    let germline = run(&[
        "germline",
        "tests/fixtures/synth-snp-match/truth.vcf",
        "tests/fixtures/synth-snp-match/query.vcf",
        "-r",
        "tests/fixtures/synth-snp-match/ref.fa",
        "-o",
        &path(&germline_prefix),
    ]);
    assert_success(&germline, "germline");
    assert!(germline.stdout.is_empty());
    assert!(germline.stderr.is_empty());

    let runinfo = parse_json(&germline_prefix.with_extension("runinfo.json"));
    let runinfo = runinfo.as_object().expect("runinfo root is an object");
    assert_eq!(runinfo.get("name"), Some(&Value::String("hap.py".into())));
    assert!(runinfo.get("timestamp").is_some_and(Value::is_string));
    assert!(runinfo.get("environment").is_some_and(Value::is_object));
    assert!(runinfo.get("final_args").is_some_and(Value::is_object));
    assert!(runinfo.get("metadata").is_some_and(Value::is_object));
    assert_metrics_document(&parse_gzip_json(
        &germline_prefix.with_extension("metrics.json.gz"),
    ));

    let somatic_prefix = output_root.join("somatic");
    let somatic = run(&[
        "somatic",
        "tests/fixtures/synth-snp-match/truth.vcf",
        "tests/fixtures/synth-snp-match/query.vcf",
        "-r",
        "tests/fixtures/synth-snp-match/ref.fa",
        "-o",
        &path(&somatic_prefix),
        "--feature-table",
        "generic",
    ]);
    assert_success(&somatic, "somatic");
    assert!(somatic.stderr.is_empty());
    assert!(String::from_utf8_lossy(&somatic.stdout).contains("total.truth,total.query,tp,fp,fn"));
    assert_metrics_document(&parse_json(&somatic_prefix.with_extension("metrics.json")));

    let pre_output = output_root.join("pre.vcf.gz");
    let pre = run(&[
        "pre",
        "tests/fixtures/preprocess-pass-fixchr/input.vcf",
        &path(&pre_output),
        "-r",
        "tests/fixtures/preprocess-pass-fixchr/ref.fa",
        "-R",
        "tests/fixtures/preprocess-pass-fixchr/restrict.bed",
    ]);
    assert_success(&pre, "pre");
    assert!(pre.stdout.is_empty());
    assert!(pre.stderr.is_empty());
    assert!(PathBuf::from(format!("{}.tbi", pre_output.display())).is_file());

    let ftx_output = output_root.join("features.csv");
    let ftx = run(&[
        "ftx",
        "tests/fixtures/synth-snp-match/query.vcf",
        "-o",
        &path(&ftx_output),
        "-r",
        "tests/fixtures/synth-snp-match/ref.fa",
        "--feature-table",
        "generic",
    ]);
    assert_success(&ftx, "ftx");
    assert!(ftx.stdout.is_empty());
    assert!(ftx.stderr.is_empty());
    assert!(
        fs::read_to_string(ftx_output)
            .expect("FTX CSV")
            .starts_with(",CHROM,POS")
    );

    let quantify_prefix = output_root.join("quantify");
    let quantify_input = germline_prefix.with_extension("vcf.gz");
    let quantify = run(&[
        "quantify",
        &path(&quantify_input),
        "-o",
        &path(&quantify_prefix),
        "-r",
        "tests/fixtures/synth-snp-match/ref.fa",
    ]);
    assert_success(&quantify, "quantify");
    assert!(quantify.stdout.is_empty());
    assert!(quantify.stderr.is_empty());
    assert_metrics_document(&parse_gzip_json(
        &quantify_prefix.with_extension("metrics.json.gz"),
    ));

    let validate_output = output_root.join("validate.json");
    let validate = run(&[
        "validate",
        "tests/fixtures/validate-summary/input.vcf",
        "-o",
        &path(&validate_output),
    ]);
    assert_success(&validate, "validate");
    assert!(validate.stdout.is_empty());
    let validate_stderr = String::from_utf8(validate.stderr).expect("validate stderr is UTF-8");
    assert!(validate_stderr.contains("[I] Total VCF records:"));
    assert!(validate_stderr.contains("[I] Non-reference VCF records:"));
    let validation = parse_json(&validate_output);
    let validation = validation
        .as_object()
        .expect("validation summary is an object");
    for key in [
        "records",
        "nonref",
        "ref",
        "haploid",
        "diploid",
        "polyploid",
    ] {
        assert!(validation.get(key).is_some_and(Value::is_number), "{key}");
    }
    assert!(validation.get("male").is_some_and(Value::is_boolean));
}

#[test]
fn native_vcf_and_bcf_inputs_have_equivalent_preprocessing_behavior() {
    let directory = tempdir().expect("temporary output directory");
    let vcf_output = directory.path().join("from-vcf.vcf.gz");
    let bcf_output = directory.path().join("from-bcf.vcf.gz");
    for (input, output) in [
        (
            "verification/assets/fixtures/pre-matrix/input-bcf.source.vcf",
            &vcf_output,
        ),
        (
            "verification/assets/fixtures/pre-matrix/input.bcf",
            &bcf_output,
        ),
    ] {
        let result = run(&[
            "pre",
            input,
            &path(output),
            "-r",
            "verification/assets/fixtures/pre-matrix/reference.fa",
            "-R",
            "verification/assets/fixtures/pre-matrix/all.bed",
        ]);
        assert_success(&result, input);
    }

    let records = |artifact: &Path| {
        String::from_utf8(read_gzip(artifact))
            .expect("preprocessed VCF is UTF-8")
            .lines()
            .filter(|line| !line.starts_with("##bcftools_"))
            .map(str::to_owned)
            .collect::<Vec<_>>()
    };
    assert_eq!(records(&vcf_output), records(&bcf_output));
}

#[test]
fn complex_native_vcf_and_bcf_records_have_equivalent_behavior() {
    let directory = tempdir().expect("temporary fixture directory");
    let reference = directory.path().join("reference.fa");
    let index = directory.path().join("reference.fa.fai");
    let bed = directory.path().join("all.bed");
    let input = directory.path().join("complex.vcf");
    let native_bcf = directory.path().join("complex.bcf");
    fs::write(
        &reference,
        format!(">chr1\n{}\n>chr2\n{}\n", "A".repeat(80), "A".repeat(80)),
    )
    .unwrap();
    fs::write(&index, "chr1\t80\t6\t80\t81\nchr2\t80\t93\t80\t81\n").unwrap();
    fs::write(&bed, "chr1\t0\t80\nchr2\t0\t80\n").unwrap();
    fs::write(
        &input,
        concat!(
            "##fileformat=VCFv4.2\n",
            "##contig=<ID=chr1,length=80>\n##contig=<ID=chr2,length=80>\n",
            "##INFO=<ID=FLAG,Number=0,Type=Flag,Description=flag>\n",
            "##INFO=<ID=I,Number=1,Type=Integer,Description=int>\n",
            "##INFO=<ID=F,Number=1,Type=Float,Description=float>\n",
            "##INFO=<ID=V,Number=.,Type=Integer,Description=vector>\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=genotype>\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=depth>\n",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=depths>\n",
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=likelihoods>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPHASED\tHAPLOID\tPOLYPLOID\tMISSING\n",
            "chr1\t10\tmulti\tA\tC,G\t50\tPASS\tFLAG;I=7;F=1.5;V=1,2,3\tGT:DP:AD:PL\t1|2:18:3,7,8:60,50,40,30,20,10\t1:9:2,7\t0/1/2:21:4,8,9:90,80,70,60,50,40,30,20,10,0\t./.:.:.:.\n",
            "chr2\t20\tsnp\tA\tT\t.\tPASS\tI=-3;F=0.25;V=5,8\tGT:DP:AD\t0|1:12:6,6\t.:.\t1/1/1:30:0,30\t./.:.:.\n",
        ),
    ).unwrap();

    let common = [
        "-r".to_string(),
        path(&reference),
        "-R".to_string(),
        path(&bed),
        "--no-leftshift".to_string(),
        "--no-decompose".to_string(),
        "--gender".to_string(),
        "none".to_string(),
    ];
    let mut encode = vec![
        "pre".to_string(),
        path(&input),
        path(&native_bcf),
        "--bcf".to_string(),
    ];
    encode.extend(common.iter().cloned());
    assert_success(&run_strings(&encode), "encode complex BCF");

    let from_vcf = directory.path().join("from-vcf.vcf.gz");
    let from_bcf = directory.path().join("from-bcf.vcf.gz");
    for (source, destination) in [(&input, &from_vcf), (&native_bcf, &from_bcf)] {
        let mut arguments = vec!["pre".to_string(), path(source), path(destination)];
        arguments.extend(common.iter().cloned());
        assert_success(&run_strings(&arguments), &path(source));
    }
    let data_lines = |artifact: &Path| {
        String::from_utf8(read_gzip(artifact))
            .unwrap()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let mut columns = line.split('\t').map(str::to_owned).collect::<Vec<_>>();
                let width = columns[8].split(':').count();
                for sample in &mut columns[9..] {
                    let present = sample.split(':').count();
                    if present < width {
                        sample.push_str(&":.".repeat(width - present));
                    }
                }
                columns.join("\t")
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(data_lines(&from_vcf), data_lines(&from_bcf));
}
