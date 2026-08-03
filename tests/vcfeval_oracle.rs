#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/vcfeval-simple")
        .join(name)
}

#[test]
fn vcfeval_handoff_matches_saved_legacy_oracle() {
    let temp = tempfile::tempdir().expect("create temporary directory");
    let template = temp.path().join("template.sdf");
    fs::create_dir(&template).expect("create fake SDF template");
    let fake_rtg = temp.path().join("rtg");
    let handoff = fixture("rtg-output.vcf");
    let indexed_handoff = temp.path().join("rtg-output.vcf.gz");
    let (headers, records) = hap_rs::vcf::load_raw_vcf(&handoff).expect("load RTG handoff");
    let record_lines = records
        .iter()
        .map(hap_rs::vcf::RawVcfRecord::to_line)
        .collect::<Vec<_>>();
    hap_rs::vcf::write_indexed_vcf(
        &indexed_handoff,
        &headers,
        record_lines.iter().map(String::as_str),
    )
    .expect("index RTG handoff");
    fs::write(
        &fake_rtg,
        format!(
            "#!/bin/sh\nout=''\nwhile [ \"$#\" -gt 0 ]; do\n  if [ \"$1\" = '-o' ]; then out=\"$2\"; shift 2; else shift; fi\ndone\nmkdir -p \"$out\"\ncp '{}' \"$out/output.vcf.gz\"\ncp '{}.tbi' \"$out/output.vcf.gz.tbi\"\n",
            indexed_handoff.display(),
            indexed_handoff.display()
        ),
    )
    .expect("write fake RTG executable");
    fs::set_permissions(&fake_rtg, fs::Permissions::from_mode(0o755))
        .expect("make fake RTG executable runnable");

    let prefix = temp.path().join("result");
    let output = Command::new(env!("CARGO_BIN_EXE_hap"))
        .args([
            "germline",
            fixture("truth.vcf").to_str().unwrap(),
            fixture("query.vcf").to_str().unwrap(),
            "--reference",
            fixture("ref.fa").to_str().unwrap(),
            "--report-prefix",
            prefix.to_str().unwrap(),
            "--engine",
            "vcfeval",
            "--engine-vcfeval-path",
            fake_rtg.to_str().unwrap(),
            "--engine-vcfeval-template",
            template.to_str().unwrap(),
            "--threads",
            "1",
            "--no-json",
            "--no-roc",
        ])
        .output()
        .expect("run hap vcfeval contract");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    assert_eq!(
        fs::read(prefix.with_extension("summary.csv")).expect("read Rust summary"),
        fs::read(fixture("expected/summary.csv")).expect("read legacy summary")
    );
    let records = hap_rs::vcf::read_text(&prefix.with_extension("vcf.gz"))
        .expect("read Rust comparison VCF")
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| format!("{line}\n"))
        .collect::<String>();
    assert_eq!(
        records,
        fs::read_to_string(fixture("expected/vcf.records")).expect("read legacy VCF records")
    );
}
