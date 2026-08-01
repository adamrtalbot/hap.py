//! Pure-Rust parity harness used by the Nextflow evaluation pipeline.
//!
//! This module intentionally has no verifier-only dependencies.  It keeps the
//! comparison policy beside the implementation being evaluated. VCF/TBI
//! indexes are queried with `tabix`; BCF/CSI indexes are followed and decoded
//! internally so BCF verification has no external-tool dependency.

use anyhow::{Context, Result, anyhow, bail};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

const VOLATILE_VCF_PREFIXES: &[&str] = &[
    "##fileDate=",
    "##source=",
    "##commandline=",
    "##CL=",
    "##hap.py",
    "##som.py",
    "##reference=",
    "##bcftools_viewCommand=",
    "##bcftools_concatCommand=",
    "##bcftools_annotateCommand=",
    "##bcftools_normCommand=",
    "##bcftools_mergeCommand=",
    "##bcftools_viewVersion=",
    "##bcftools_concatVersion=",
    "##bcftools_annotateVersion=",
    "##bcftools_normVersion=",
    "##bcftools_mergeVersion=",
];

const VOLATILE_JSON_PATHS: &[&str] = &[
    "$/dist",
    "$/environment",
    "$/mac_ver",
    "$/python_implementation",
    "$/python_prefix",
    "$/python_version",
    "$/uname",
    "$/metadata/required/description",
];

const ROC_IDENTITY_FIELDS: &[&str] = &[
    "Type", "Subtype", "Subset", "Filter", "Genotype", "QQ.Field", "QQ",
];

#[derive(Clone, Debug, PartialEq)]
enum Json {
    Null,
    Bool(bool),
    Int(i64),
    Float(f64),
    String(String),
    Array(Vec<Json>),
    Object(BTreeMap<String, Json>),
}

impl Json {
    fn type_name(&self) -> &'static str {
        match self {
            Self::Null => "NoneType",
            Self::Bool(_) => "bool",
            Self::Int(_) => "int",
            Self::Float(_) => "float",
            Self::String(_) => "str",
            Self::Array(_) => "list",
            Self::Object(_) => "dict",
        }
    }

    fn as_object_mut(&mut self) -> Option<&mut BTreeMap<String, Json>> {
        match self {
            Self::Object(value) => Some(value),
            _ => None,
        }
    }

    fn as_array_mut(&mut self) -> Option<&mut Vec<Json>> {
        match self {
            Self::Array(value) => Some(value),
            _ => None,
        }
    }

    fn as_str(&self) -> Option<&str> {
        match self {
            Self::String(value) => Some(value),
            _ => None,
        }
    }

    fn as_i64(&self) -> Option<i64> {
        match self {
            Self::Int(value) => Some(*value),
            Self::Float(value) if value.is_finite() => Some(*value as i64),
            _ => None,
        }
    }
}

struct JsonParser<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> JsonParser<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            bytes: text.as_bytes(),
            offset: 0,
        }
    }

    fn parse(mut self) -> Result<Json> {
        let value = self.value()?;
        self.whitespace();
        if self.offset != self.bytes.len() {
            bail!("unexpected trailing JSON content at byte {}", self.offset);
        }
        Ok(value)
    }

    fn value(&mut self) -> Result<Json> {
        self.whitespace();
        match self.peek() {
            Some(b'n') => {
                self.literal(b"null")?;
                Ok(Json::Null)
            }
            Some(b't') => {
                self.literal(b"true")?;
                Ok(Json::Bool(true))
            }
            Some(b'f') => {
                self.literal(b"false")?;
                Ok(Json::Bool(false))
            }
            Some(b'"') => Ok(Json::String(self.string()?)),
            Some(b'[') => self.array(),
            Some(b'{') => self.object(),
            Some(b'-' | b'0'..=b'9') => self.number(),
            Some(byte) => bail!(
                "unexpected JSON byte {:?} at byte {}",
                char::from(byte),
                self.offset
            ),
            None => bail!("unexpected end of JSON"),
        }
    }

    fn array(&mut self) -> Result<Json> {
        self.expect(b'[')?;
        let mut values = Vec::new();
        self.whitespace();
        if self.consume(b']') {
            return Ok(Json::Array(values));
        }
        loop {
            values.push(self.value()?);
            self.whitespace();
            if self.consume(b']') {
                return Ok(Json::Array(values));
            }
            self.expect(b',')?;
        }
    }

    fn object(&mut self) -> Result<Json> {
        self.expect(b'{')?;
        let mut values = BTreeMap::new();
        self.whitespace();
        if self.consume(b'}') {
            return Ok(Json::Object(values));
        }
        loop {
            self.whitespace();
            let key = self.string()?;
            self.whitespace();
            self.expect(b':')?;
            values.insert(key, self.value()?);
            self.whitespace();
            if self.consume(b'}') {
                return Ok(Json::Object(values));
            }
            self.expect(b',')?;
        }
    }

    fn string(&mut self) -> Result<String> {
        self.expect(b'"')?;
        let mut out = String::new();
        while let Some(byte) = self.next() {
            match byte {
                b'"' => return Ok(out),
                b'\\' => match self
                    .next()
                    .ok_or_else(|| anyhow!("truncated JSON escape"))?
                {
                    b'"' => out.push('"'),
                    b'\\' => out.push('\\'),
                    b'/' => out.push('/'),
                    b'b' => out.push('\u{8}'),
                    b'f' => out.push('\u{c}'),
                    b'n' => out.push('\n'),
                    b'r' => out.push('\r'),
                    b't' => out.push('\t'),
                    b'u' => {
                        let first = self.hex_quad()?;
                        if (0xD800..=0xDBFF).contains(&first) {
                            self.expect(b'\\')?;
                            self.expect(b'u')?;
                            let second = self.hex_quad()?;
                            if !(0xDC00..=0xDFFF).contains(&second) {
                                bail!("invalid low surrogate in JSON string");
                            }
                            let scalar = 0x10000
                                + ((u32::from(first) - 0xD800) << 10)
                                + (u32::from(second) - 0xDC00);
                            out.push(
                                char::from_u32(scalar)
                                    .ok_or_else(|| anyhow!("invalid Unicode scalar"))?,
                            );
                        } else if (0xDC00..=0xDFFF).contains(&first) {
                            bail!("unpaired low surrogate in JSON string");
                        } else {
                            out.push(
                                char::from_u32(u32::from(first))
                                    .ok_or_else(|| anyhow!("invalid Unicode scalar"))?,
                            );
                        }
                    }
                    escape => bail!("invalid JSON escape {:?}", char::from(escape)),
                },
                0..=31 => bail!("unescaped control character in JSON string"),
                32..=127 => out.push(char::from(byte)),
                _ => {
                    let start = self.offset - 1;
                    let width =
                        utf8_width(byte).ok_or_else(|| anyhow!("invalid UTF-8 in JSON string"))?;
                    let end = start + width;
                    if end > self.bytes.len() {
                        bail!("truncated UTF-8 in JSON string");
                    }
                    let text = std::str::from_utf8(&self.bytes[start..end])?;
                    out.push_str(text);
                    self.offset = end;
                }
            }
        }
        bail!("unterminated JSON string")
    }

    fn hex_quad(&mut self) -> Result<u16> {
        let end = self.offset + 4;
        if end > self.bytes.len() {
            bail!("truncated Unicode escape");
        }
        let text = std::str::from_utf8(&self.bytes[self.offset..end])?;
        self.offset = end;
        u16::from_str_radix(text, 16).context("invalid Unicode escape")
    }

    fn number(&mut self) -> Result<Json> {
        let start = self.offset;
        self.consume(b'-');
        match self.peek() {
            Some(b'0') => self.offset += 1,
            Some(b'1'..=b'9') => {
                self.offset += 1;
                while matches!(self.peek(), Some(b'0'..=b'9')) {
                    self.offset += 1;
                }
            }
            _ => bail!("invalid JSON number at byte {start}"),
        }
        let mut float = false;
        if self.consume(b'.') {
            float = true;
            self.digits()?;
        }
        if matches!(self.peek(), Some(b'e' | b'E')) {
            float = true;
            self.offset += 1;
            if matches!(self.peek(), Some(b'+' | b'-')) {
                self.offset += 1;
            }
            self.digits()?;
        }
        let raw = std::str::from_utf8(&self.bytes[start..self.offset])?;
        if float {
            Ok(Json::Float(
                raw.parse::<f64>().context("invalid JSON float")?,
            ))
        } else {
            Ok(Json::Int(
                raw.parse::<i64>().context("invalid JSON integer")?,
            ))
        }
    }

    fn digits(&mut self) -> Result<()> {
        let start = self.offset;
        while matches!(self.peek(), Some(b'0'..=b'9')) {
            self.offset += 1;
        }
        if self.offset == start {
            bail!("JSON number requires a digit at byte {start}");
        }
        Ok(())
    }

    fn literal(&mut self, expected: &[u8]) -> Result<()> {
        if self.bytes.get(self.offset..self.offset + expected.len()) != Some(expected) {
            bail!("invalid JSON literal at byte {}", self.offset);
        }
        self.offset += expected.len();
        Ok(())
    }

    fn whitespace(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\n' | b'\r' | b'\t')) {
            self.offset += 1;
        }
    }

    fn expect(&mut self, expected: u8) -> Result<()> {
        if self.consume(expected) {
            Ok(())
        } else {
            bail!(
                "expected {:?} at byte {}",
                char::from(expected),
                self.offset
            )
        }
    }

    fn consume(&mut self, expected: u8) -> bool {
        if self.peek() == Some(expected) {
            self.offset += 1;
            true
        } else {
            false
        }
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.offset).copied()
    }

    fn next(&mut self) -> Option<u8> {
        let value = self.peek()?;
        self.offset += 1;
        Some(value)
    }
}

fn utf8_width(first: u8) -> Option<usize> {
    match first {
        0xC2..=0xDF => Some(2),
        0xE0..=0xEF => Some(3),
        0xF0..=0xF4 => Some(4),
        _ => None,
    }
}

#[derive(Clone, Debug)]
struct FileStatus {
    name: String,
    ok: bool,
    reason: String,
}

pub struct CompareOutputs<'a> {
    pub legacy_dir: &'a Path,
    pub rust_dir: &'a Path,
    pub prefix: &'a str,
    pub case: &'a str,
    pub sample: &'a str,
    pub report: &'a Path,
    pub status: &'a Path,
}

pub struct AggregateReport<'a> {
    pub image: &'a str,
    pub hap_bin: &'a str,
    pub markdown: &'a Path,
    pub csv: &'a Path,
    pub inputs: &'a [PathBuf],
}

pub fn compare_outputs(args: &CompareOutputs<'_>) -> Result<()> {
    let legacy = list_outputs(args.legacy_dir, args.prefix)?;
    let rust = list_outputs(args.rust_dir, args.prefix)?;
    let names = legacy
        .keys()
        .chain(rust.keys())
        .cloned()
        .collect::<BTreeSet<_>>();

    let mut files = Vec::new();
    let mut report = vec![
        format!("# parity diff: case={} sample={}", args.case, args.sample),
        format!("legacy_dir={}", absolute_display(args.legacy_dir).display()),
        format!("rust_dir={}", absolute_display(args.rust_dir).display()),
        String::new(),
    ];

    for name in names {
        let Some(left) = legacy.get(&name) else {
            files.push(FileStatus {
                name: name.clone(),
                ok: false,
                reason: "only-in-rust".to_string(),
            });
            report.push(format!("FAIL  {name:<40} (only in rust/)"));
            continue;
        };
        let Some(right) = rust.get(&name) else {
            files.push(FileStatus {
                name: name.clone(),
                ok: false,
                reason: "only-in-legacy".to_string(),
            });
            report.push(format!("FAIL  {name:<40} (only in legacy/)"));
            continue;
        };
        let (ok, detail) = compare_file(&name, left, right);
        let detail = detail.unwrap_or_else(|error| format!("comparison failed: {error:#}\n"));
        files.push(FileStatus {
            name: name.clone(),
            ok,
            reason: if ok {
                String::new()
            } else {
                detail
                    .trim()
                    .lines()
                    .next()
                    .unwrap_or("mismatch")
                    .to_string()
            },
        });
        report.push(format!("{}  {name}", if ok { "PASS" } else { "FAIL" }));
        if !detail.is_empty() {
            report.push(detail);
        }
    }

    if files.is_empty() {
        files.push(FileStatus {
            name: "<artifact-set>".to_string(),
            ok: false,
            reason: "no-matching-artifacts".to_string(),
        });
        report.push(format!(
            "FAIL  no artifacts begin with prefix {:?}",
            args.prefix
        ));
    }
    let failed = files.iter().filter(|entry| !entry.ok).count();
    create_parent(args.status)?;
    create_parent(args.report)?;
    fs::write(
        args.status,
        render_status(args.case, args.sample, failed, &files),
    )?;
    fs::write(args.report, format!("{}\n", report.join("\n")))?;
    Ok(())
}

fn absolute_display(path: &Path) -> PathBuf {
    path.canonicalize()
        .unwrap_or_else(|_| std::env::current_dir().unwrap_or_default().join(path))
}

fn create_parent(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

fn list_outputs(dir: &Path, prefix: &str) -> Result<BTreeMap<String, PathBuf>> {
    let mut out = BTreeMap::new();
    if !dir.is_dir() {
        return Ok(out);
    }
    for entry in fs::read_dir(dir)? {
        let path = entry?.path();
        let Some(name) = path.file_name().and_then(|value| value.to_str()) else {
            continue;
        };
        if path.is_file() && name.starts_with(prefix) {
            out.insert(name.to_string(), path);
        }
    }
    Ok(out)
}

fn compare_file(name: &str, left: &Path, right: &Path) -> (bool, Result<String>) {
    let result = if name.ends_with(".roc.all.csv.gz") {
        validate_roc(left)
            .map(|_| ())
            .with_context(|| "legacy ROC integrity check failed")
            .and_then(|()| compare_csv(left, right, "csv.gz"))
    } else if name.ends_with(".vcf.gz") || name.ends_with(".vcf") {
        compare_vcf(left, right)
    } else if name.ends_with(".bcf") {
        compare_bcf(left, right)
    } else if name.ends_with(".csv.gz") {
        compare_csv(left, right, "csv.gz")
    } else if name.ends_with(".csv") {
        compare_csv(left, right, "csv")
    } else if name.ends_with(".json.gz") || name.ends_with(".json") {
        compare_json(left, right)
    } else if name.ends_with(".tbi") || name.ends_with(".csi") {
        compare_index(left, right)
    } else {
        compare_bytes(left, right)
    };
    match result {
        Ok(detail) => (detail.is_empty(), Ok(detail)),
        Err(error) => (false, Err(error)),
    }
}

fn read_text(path: &Path) -> Result<String> {
    let bytes = fs::read(path).with_context(|| format!("failed to read {}", path.display()))?;
    let decoded = if path.extension().and_then(|value| value.to_str()) == Some("gz") {
        let mut decoder = flate2::read::MultiGzDecoder::new(bytes.as_slice());
        let mut decoded = Vec::new();
        decoder
            .read_to_end(&mut decoded)
            .with_context(|| format!("failed to decompress {}", path.display()))?;
        decoded
    } else {
        bytes
    };
    Ok(String::from_utf8_lossy(&decoded).into_owned())
}

fn compare_vcf(left: &Path, right: &Path) -> Result<String> {
    let (left_headers, left_records) = vcf_parts(left)?;
    let (right_headers, right_records) = vcf_parts(right)?;
    let mut detail = String::new();
    if left_headers != right_headers {
        detail.push_str(&diff_block(
            &left_headers,
            &right_headers,
            &format!("vcf headers {}", file_name(left)),
        ));
    }
    if left_records != right_records {
        detail.push_str(&diff_block(
            &left_records,
            &right_records,
            &format!("vcf records {}", file_name(left)),
        ));
    }
    Ok(detail)
}

fn compare_bcf(left: &Path, right: &Path) -> Result<String> {
    let (left_headers, left_records) = bcf_parts(left)?;
    let (right_headers, right_records) = bcf_parts(right)?;
    let mut detail = String::new();
    if left_headers != right_headers {
        detail.push_str(&diff_block(
            &left_headers,
            &right_headers,
            &format!("bcf headers {}", file_name(left)),
        ));
    }
    if left_records != right_records {
        detail.push_str(&diff_block(
            &left_records,
            &right_records,
            &format!("bcf records {}", file_name(left)),
        ));
    }
    Ok(detail)
}

fn bcf_parts(path: &Path) -> Result<(Vec<String>, Vec<String>)> {
    let (headers, records) = crate::vcf::load_raw_vcf(path)
        .with_context(|| format!("failed to decode BCF {}", path.display()))?;
    Ok((
        headers
            .into_iter()
            .filter(|line| {
                !VOLATILE_VCF_PREFIXES
                    .iter()
                    .any(|prefix| line.starts_with(prefix))
            })
            .collect(),
        records.into_iter().map(|record| record.to_line()).collect(),
    ))
}

fn vcf_parts(path: &Path) -> Result<(Vec<String>, Vec<String>)> {
    let mut headers = Vec::new();
    let mut records = Vec::new();
    for line in read_text(path)?.lines().filter(|line| !line.is_empty()) {
        if line.starts_with('#') {
            if !VOLATILE_VCF_PREFIXES
                .iter()
                .any(|prefix| line.starts_with(prefix))
            {
                headers.push(line.to_string());
            }
        } else {
            records.push(line.to_string());
        }
    }
    Ok((headers, records))
}

fn compare_csv(left: &Path, right: &Path, label: &str) -> Result<String> {
    let left_lines = read_text(left)?
        .lines()
        .map(str::to_string)
        .collect::<Vec<_>>();
    let right_lines = read_text(right)?
        .lines()
        .map(str::to_string)
        .collect::<Vec<_>>();
    Ok(diff_block(
        &left_lines,
        &right_lines,
        &format!("{label} {}", file_name(left)),
    ))
}

fn diff_block(left: &[String], right: &[String], label: &str) -> String {
    if left == right {
        return String::new();
    }
    let left_set = left.iter().collect::<BTreeSet<_>>();
    let right_set = right.iter().collect::<BTreeSet<_>>();
    let left_only = left_set.difference(&right_set).copied().collect::<Vec<_>>();
    let right_only = right_set.difference(&left_set).copied().collect::<Vec<_>>();
    let mut lines = vec![format!("---------- {label} ----------")];
    if !left_only.is_empty() {
        lines.push("only in legacy/:".to_string());
        lines.extend(left_only.iter().take(50).map(|line| format!("  -{line}")));
        if left_only.len() > 50 {
            lines.push(format!("  ... and {} more", left_only.len() - 50));
        }
    }
    if !right_only.is_empty() {
        lines.push("only in rust/:".to_string());
        lines.extend(right_only.iter().take(50).map(|line| format!("  +{line}")));
        if right_only.len() > 50 {
            lines.push(format!("  ... and {} more", right_only.len() - 50));
        }
    }
    if left_only.is_empty() && right_only.is_empty() {
        lines.push(format!(
            "same unique values but different ordering or duplicate counts: legacy={} rust={}",
            left.len(),
            right.len()
        ));
    }
    format!("{}\n", lines.join("\n"))
}

fn compare_json(left: &Path, right: &Path) -> Result<String> {
    let left_text = read_text(left)?;
    let right_text = read_text(right)?;
    let mut left_json = match JsonParser::new(&left_text).parse() {
        Ok(value) => value,
        Err(error) => {
            return Ok(format!(
                "legacy JSON parse failed for {}: {error}\n",
                file_name(left)
            ));
        }
    };
    let mut right_json = match JsonParser::new(&right_text).parse() {
        Ok(value) => value,
        Err(error) => {
            return Ok(format!(
                "rust JSON parse failed for {}: {error}\n",
                file_name(right)
            ));
        }
    };
    canonicalize_metric_json(&mut left_json);
    canonicalize_metric_json(&mut right_json);
    let mut differences = Vec::new();
    json_differences(&left_json, &right_json, "$", &mut differences);
    if differences.is_empty() {
        return Ok(String::new());
    }
    let remaining = differences.len().saturating_sub(100);
    differences.truncate(100);
    if remaining > 0 {
        differences.push(format!("... and {remaining} more JSON differences"));
    }
    Ok(format!(
        "---------- json {} ----------\n{}\n",
        file_name(left),
        differences.join("\n")
    ))
}

fn canonicalize_metric_json(value: &mut Json) {
    let Some(object) = value.as_object_mut() else {
        return;
    };
    if let Some(metrics) = object.get_mut("metrics").and_then(Json::as_array_mut) {
        let ids = metrics
            .iter()
            .map(|metric| {
                let Json::Object(metric) = metric else {
                    return None;
                };
                metric.get("id").and_then(Json::as_str).map(str::to_string)
            })
            .collect::<Option<Vec<_>>>();
        if let Some(ids) = ids
            && ids.iter().collect::<BTreeSet<_>>().len() == ids.len()
        {
            metrics.sort_by(|left, right| metric_id(left).cmp(metric_id(right)));
        }
        for metric in metrics {
            let Json::Object(metric) = metric else {
                continue;
            };
            let Some(Json::Array(data)) = metric.get_mut("data") else {
                continue;
            };
            for column in data {
                let Json::Object(column) = column else {
                    continue;
                };
                if column.get("id").and_then(Json::as_str) != Some("types") {
                    continue;
                }
                if let Some(Json::Array(values)) = column.get_mut("values") {
                    for (index, value) in values.iter_mut().enumerate() {
                        *value = Json::Int(index as i64);
                    }
                }
            }
        }
    }
    if let Some(Json::Array(run_info)) = object.get_mut("runInfo") {
        for entry in run_info {
            let Json::Object(entry) = entry else {
                continue;
            };
            if entry.get("key").and_then(Json::as_str) == Some("commandline") {
                entry.insert(
                    "value".to_string(),
                    Json::String("<normalized-commandline>".to_string()),
                );
            }
        }
    }
}

fn metric_id(value: &Json) -> &str {
    let Json::Object(value) = value else {
        return "";
    };
    value.get("id").and_then(Json::as_str).unwrap_or("")
}

fn json_differences(left: &Json, right: &Json, path: &str, out: &mut Vec<String>) {
    if left.type_name() != right.type_name() {
        out.push(format!(
            "{path}: type legacy={} rust={}",
            left.type_name(),
            right.type_name()
        ));
        return;
    }
    match (left, right) {
        (Json::Object(left), Json::Object(right)) => {
            for key in left.keys().chain(right.keys()).collect::<BTreeSet<_>>() {
                if json_path_is_volatile(path, key) {
                    continue;
                }
                let child = json_pointer(path, key);
                match (left.get(key), right.get(key)) {
                    (None, Some(_)) => out.push(format!("{child}: only in rust")),
                    (Some(_), None) => out.push(format!("{child}: only in legacy")),
                    (Some(left), Some(right)) => json_differences(left, right, &child, out),
                    (None, None) => {}
                }
            }
        }
        (Json::Array(left), Json::Array(right)) => {
            if left.len() != right.len() {
                out.push(format!(
                    "{path}: length legacy={} rust={}",
                    left.len(),
                    right.len()
                ));
            }
            for (index, (left, right)) in left.iter().zip(right).enumerate() {
                json_differences(left, right, &json_pointer(path, &index.to_string()), out);
            }
        }
        _ if left != right => out.push(format!(
            "{path}: legacy={} rust={}",
            json_repr(left),
            json_repr(right)
        )),
        _ => {}
    }
}

fn json_path_is_volatile(path: &str, key: &str) -> bool {
    key == "timestamp" || VOLATILE_JSON_PATHS.contains(&json_pointer(path, key).as_str())
}

fn json_pointer(path: &str, key: &str) -> String {
    format!("{path}/{}", key.replace('~', "~0").replace('/', "~1"))
}

fn json_repr(value: &Json) -> String {
    match value {
        Json::Null => "None".to_string(),
        Json::Bool(value) => if *value { "True" } else { "False" }.to_string(),
        Json::Int(value) => value.to_string(),
        Json::Float(value) => value.to_string(),
        Json::String(value) => format!("'{}'", value.replace('\\', "\\\\").replace('\'', "\\'")),
        Json::Array(_) => "<list>".to_string(),
        Json::Object(_) => "<dict>".to_string(),
    }
}

fn compare_index(left: &Path, right: &Path) -> Result<String> {
    if !left.is_file() || left.metadata()?.len() == 0 {
        return Ok(format!(
            "legacy index missing or empty: {}\n",
            left.display()
        ));
    }
    if !right.is_file() || right.metadata()?.len() == 0 {
        return Ok(format!(
            "rust index missing or empty: {}\n",
            right.display()
        ));
    }
    if let Err(error) = tabix_range_query(left) {
        return Ok(format!("legacy index query failed: {error}\n"));
    }
    if let Err(error) = tabix_range_query(right) {
        return Ok(format!("rust index query failed: {error}\n"));
    }
    Ok(String::new())
}

fn tabix_range_query(index: &Path) -> Result<()> {
    let extension = index.extension().and_then(|value| value.to_str());
    if !matches!(extension, Some("tbi" | "csi")) {
        bail!("unsupported query index: {}", file_name(index));
    }
    let mut data = index.to_path_buf();
    data.set_extension("");
    if !data.is_file() {
        bail!("indexed data file is absent: {}", data.display());
    }
    if data.extension().and_then(|value| value.to_str()) == Some("bcf") {
        return bcf_range_query(&data, index);
    }
    let listed = Command::new("tabix")
        .args(["-l"])
        .arg(&data)
        .output()
        .map_err(|error| {
            if error.kind() == std::io::ErrorKind::NotFound {
                anyhow!("tabix executable is required for complete index validation")
            } else {
                anyhow!(error)
            }
        })?;
    if !listed.status.success() {
        let stderr = String::from_utf8_lossy(&listed.stderr).trim().to_string();
        bail!(
            "{}",
            if stderr.is_empty() {
                "tabix could not list references"
            } else {
                &stderr
            }
        );
    }
    let contigs = String::from_utf8_lossy(&listed.stdout)
        .lines()
        .filter(|line| !line.is_empty())
        .map(str::to_string)
        .collect::<Vec<_>>();
    if contigs.is_empty() {
        bail!("index contains no queryable references");
    }
    let mut source: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for line in read_text(&data)?
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
    {
        let fields = line.splitn(3, '\t').collect::<Vec<_>>();
        if fields.len() < 2 || fields[1].parse::<usize>().is_err() {
            bail!("source VCF has no valid position: {}", truncate(line, 200));
        }
        source
            .entry(fields[0].to_string())
            .or_default()
            .push(line.to_string());
    }
    let indexed_contigs = contigs.iter().cloned().collect::<BTreeSet<_>>();
    let missing = source
        .keys()
        .filter(|name| !indexed_contigs.contains(*name))
        .cloned()
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        bail!("index omits source references: {}", missing.join(", "));
    }
    for contig in contigs {
        let queried = Command::new("tabix").arg(&data).arg(&contig).output()?;
        if !queried.status.success() {
            let stderr = String::from_utf8_lossy(&queried.stderr).trim().to_string();
            bail!(
                "{}",
                if stderr.is_empty() {
                    format!("tabix query failed for {contig}")
                } else {
                    stderr
                }
            );
        }
        let indexed = String::from_utf8_lossy(&queried.stdout)
            .lines()
            .filter(|line| !line.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        let expected = source.get(&contig).cloned().unwrap_or_default();
        if indexed != expected {
            bail!(
                "indexed contig {contig} returned {} of {} source records in source order",
                indexed.len(),
                expected.len()
            );
        }
    }
    Ok(())
}

fn bcf_range_query(data: &Path, index: &Path) -> Result<()> {
    let (_, records) = crate::vcf::load_raw_vcf(data)
        .with_context(|| format!("failed to decode indexed BCF {}", data.display()))?;
    let mut source: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for record in records {
        source
            .entry(record.chrom.clone())
            .or_default()
            .push(record.to_line());
    }
    let indexed = crate::bcf::read_indexed_records(data, index)?
        .into_iter()
        .map(|(contig, records)| {
            (
                contig,
                records
                    .iter()
                    .map(|record| record.to_line())
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    for (contig, expected) in source {
        let observed = indexed.get(&contig).cloned().unwrap_or_default();
        if observed != expected {
            bail!(
                "indexed BCF contig {contig} returned {} of {} source records in source order",
                observed.len(),
                expected.len()
            );
        }
    }
    Ok(())
}

fn truncate(value: &str, count: usize) -> String {
    value.chars().take(count).collect()
}

fn compare_bytes(left: &Path, right: &Path) -> Result<String> {
    let left_bytes = fs::read(left)?;
    let right_bytes = fs::read(right)?;
    if left_bytes == right_bytes {
        Ok(String::new())
    } else {
        Ok(format!(
            "byte diff {}: legacy={}B rust={}B\n",
            file_name(left),
            left_bytes.len(),
            right_bytes.len()
        ))
    }
}

fn file_name(path: &Path) -> &str {
    path.file_name()
        .and_then(|value| value.to_str())
        .unwrap_or("<non-utf8-path>")
}

fn render_status(case: &str, sample: &str, failed: usize, files: &[FileStatus]) -> String {
    let mut out = String::new();
    out.push_str("{\n");
    out.push_str(&format!("  \"case\": {},\n", json_string(case)));
    out.push_str(&format!("  \"sample\": {},\n", json_string(sample)));
    out.push_str(&format!("  \"total\": {},\n", files.len()));
    out.push_str(&format!("  \"failed\": {failed},\n"));
    out.push_str("  \"files\": [\n");
    for (index, file) in files.iter().enumerate() {
        out.push_str("    {\n");
        out.push_str(&format!("      \"name\": {},\n", json_string(&file.name)));
        out.push_str(&format!("      \"ok\": {},\n", file.ok));
        out.push_str(&format!(
            "      \"reason\": {}\n",
            json_string(&file.reason)
        ));
        out.push_str("    }");
        if index + 1 < files.len() {
            out.push(',');
        }
        out.push('\n');
    }
    out.push_str("  ]\n}\n");
    out
}

pub fn aggregate_report(args: &AggregateReport<'_>) -> Result<()> {
    #[derive(Debug)]
    struct Row {
        case: String,
        sample: String,
        total: i64,
        failed: i64,
    }
    let mut rows = Vec::new();
    for path in args.inputs {
        let parsed = fs::read_to_string(path)
            .ok()
            .and_then(|text| JsonParser::new(&text).parse().ok());
        let Some(Json::Object(payload)) = parsed else {
            rows.push(Row {
                case: "unknown".to_string(),
                sample: path.display().to_string(),
                total: 0,
                failed: 1,
            });
            continue;
        };
        rows.push(Row {
            case: payload
                .get("case")
                .and_then(Json::as_str)
                .unwrap_or("unknown")
                .to_string(),
            sample: payload
                .get("sample")
                .and_then(Json::as_str)
                .unwrap_or("unknown")
                .to_string(),
            total: payload.get("total").and_then(Json::as_i64).unwrap_or(0),
            failed: payload.get("failed").and_then(Json::as_i64).unwrap_or(0),
        });
    }
    rows.sort_by(|left, right| (&left.case, &left.sample).cmp(&(&right.case, &right.sample)));
    let mut markdown = vec![
        "# hap.py / hap-rs parity report".to_string(),
        String::new(),
        format!("- image: `{}`", args.image),
        format!("- hap:   `{}`", args.hap_bin),
        String::new(),
        "| case | sample | files | failed | status |".to_string(),
        "|------|--------|-------|--------|--------|".to_string(),
    ];
    let mut csv = vec!["case,sample,files,failed,status".to_string()];
    for row in rows {
        let status = if row.failed == 0 { "PASS" } else { "FAIL" };
        markdown.push(format!(
            "| {} | {} | {} | {} | {status} |",
            row.case, row.sample, row.total, row.failed
        ));
        csv.push(format!(
            "{},{},{},{},{status}",
            row.case, row.sample, row.total, row.failed
        ));
    }
    fs::write(args.markdown, format!("{}\n", markdown.join("\n")))?;
    fs::write(args.csv, format!("{}\n", csv.join("\n")))?;
    Ok(())
}

pub fn validate_roc(path: &Path) -> Result<usize> {
    let text = read_text(path)?;
    let mut records = csv_records(&text)?;
    let Some(header) = records.next() else {
        bail!("{}: missing CSV header", path.display());
    };
    let missing = ROC_IDENTITY_FIELDS
        .iter()
        .filter(|field| !header.iter().any(|column| column == **field))
        .copied()
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        bail!(
            "{}: missing required columns: {}",
            path.display(),
            missing.join(", ")
        );
    }
    let indexes = ROC_IDENTITY_FIELDS
        .iter()
        .map(|field| header.iter().position(|column| column == *field).unwrap())
        .collect::<Vec<_>>();
    let type_index = header.iter().position(|column| column == "Type").unwrap();
    let mut count = 0;
    for (offset, row) in records.enumerate() {
        let line = offset + 2;
        if row.len() > header.len() {
            bail!(
                "{}:{line}: row has more fields than the header",
                path.display()
            );
        }
        if indexes
            .iter()
            .any(|index| row.get(*index).is_none_or(String::is_empty))
        {
            bail!("{}:{line}: empty ROC identity field", path.display());
        }
        let variant_type = &row[type_index];
        if variant_type != "SNP" && variant_type != "INDEL" {
            bail!(
                "{}:{line}: invalid variant type {:?}",
                path.display(),
                variant_type
            );
        }
        count += 1;
    }
    if count == 0 {
        bail!("{}: no ROC data rows", path.display());
    }
    Ok(count)
}

fn csv_records(text: &str) -> Result<impl Iterator<Item = Vec<String>> + '_> {
    let mut rows = Vec::new();
    let mut row = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = text.chars().peekable();
    while let Some(character) = chars.next() {
        match character {
            '"' if quoted && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => row.push(std::mem::take(&mut field)),
            '\n' if !quoted => {
                if field.ends_with('\r') {
                    field.pop();
                }
                row.push(std::mem::take(&mut field));
                if row.len() != 1 || !row[0].is_empty() {
                    rows.push(std::mem::take(&mut row));
                } else {
                    row.clear();
                }
            }
            _ => field.push(character),
        }
    }
    if quoted {
        bail!("unterminated quoted CSV field");
    }
    if !field.is_empty() || !row.is_empty() {
        row.push(field);
        rows.push(row);
    }
    Ok(rows.into_iter())
}

fn json_string(value: &str) -> String {
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{8}' => out.push_str("\\b"),
            '\u{c}' => out.push_str("\\f"),
            character if character.is_control() => {
                out.push_str(&format!("\\u{:04x}", character as u32));
            }
            character => out.push(character),
        }
    }
    out.push('"');
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn compare_pair(extension: &str, left: &str, right: &str) -> (bool, String) {
        let dir = tempfile::tempdir().unwrap();
        let legacy = dir.path().join(format!("legacy.{extension}"));
        let rust = dir.path().join(format!("rust.{extension}"));
        fs::write(&legacy, left).unwrap();
        fs::write(&rust, right).unwrap();
        let (ok, detail) = compare_file(&format!("result.{extension}"), &legacy, &rust);
        (ok, detail.unwrap())
    }

    fn rewrite_test_csi(path: &Path, payload: &[u8]) {
        let mut writer = noodles_bgzf::io::Writer::new(fs::File::create(path).unwrap());
        writer.write_all(payload).unwrap();
        writer.finish().unwrap();
    }

    fn add_test_csi_count_metadata(path: &Path, mapped: u64, unmapped: u64) {
        let mut payload = crate::bcf::read_uncompressed(path).unwrap();
        assert_eq!(&payload[..4], b"CSI\x01");
        assert_eq!(i32::from_le_bytes(payload[12..16].try_into().unwrap()), 0);
        assert_eq!(i32::from_le_bytes(payload[16..20].try_into().unwrap()), 1);
        assert_eq!(i32::from_le_bytes(payload[20..24].try_into().unwrap()), 1);

        // The native test writer emits one data bin followed by n_no_coor.
        // Append the CSI metadata pseudo-bin whose second chunk stores
        // mapped/unmapped counts rather than virtual offsets.
        let data_start = u64::from_le_bytes(payload[40..48].try_into().unwrap());
        let data_end = u64::from_le_bytes(payload[48..56].try_into().unwrap());
        let trailing = payload.split_off(payload.len() - 8);
        payload[20..24].copy_from_slice(&2i32.to_le_bytes());
        payload.extend_from_slice(&37_450u32.to_le_bytes());
        payload.extend_from_slice(&0u64.to_le_bytes());
        payload.extend_from_slice(&2i32.to_le_bytes());
        payload.extend_from_slice(&data_start.to_le_bytes());
        payload.extend_from_slice(&data_end.to_le_bytes());
        payload.extend_from_slice(&mapped.to_le_bytes());
        payload.extend_from_slice(&unmapped.to_le_bytes());
        payload.extend_from_slice(&trailing);
        rewrite_test_csi(path, &payload);
    }

    fn reverse_test_csi_data_chunk(path: &Path) {
        let mut payload = crate::bcf::read_uncompressed(path).unwrap();
        let start = payload[40..48].to_vec();
        let end = payload[48..56].to_vec();
        payload[40..48].copy_from_slice(&end);
        payload[48..56].copy_from_slice(&start);
        rewrite_test_csi(path, &payload);
    }

    fn write_test_indexed_bcf(root: &Path) -> PathBuf {
        fs::create_dir_all(root).unwrap();
        let bcf = root.join("result.bcf");
        let headers = [
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=100>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let record =
            crate::vcf::RawVcfRecord::from_line("chr1\t7\t.\tA\tC\t.\tPASS\t.", &bcf).unwrap();
        crate::vcf::write_raw_vcf(&bcf, &headers, &[record]).unwrap();
        bcf
    }

    #[test]
    fn metric_order_types_index_commandline_and_volatility_are_canonicalized() {
        let left = r#"{"timestamp":"one","environment":{"host":"a"},"runInfo":[{"key":"commandline","value":"legacy"}],"metrics":[{"id":"summary.metrics","data":[{"id":"types","values":[3,1]},{"id":"TP","values":[4,5]}]},{"id":"roc.all","data":[{"id":"types","values":[9]},{"id":"FP","values":[2]}]}]}"#;
        let right = r#"{"timestamp":"two","environment":{"host":"b"},"runInfo":[{"key":"commandline","value":"rust"}],"metrics":[{"id":"roc.all","data":[{"id":"types","values":[0]},{"id":"FP","values":[2]}]},{"id":"summary.metrics","data":[{"id":"types","values":[0,1]},{"id":"TP","values":[4,5]}]}]}"#;
        let (ok, detail) = compare_pair("json", left, right);
        assert!(ok, "{detail}");
    }

    #[test]
    fn nested_metric_values_and_json_scalar_types_remain_semantic() {
        let (ok, detail) = compare_pair(
            "json",
            r#"{"metrics":[{"values":[1,2]}]}"#,
            r#"{"metrics":[{"values":[1,3.0]}]}"#,
        );
        assert!(!ok);
        assert!(detail.contains("$/metrics/0/values/1"));
        assert!(detail.contains("type legacy=int rust=float"));
    }

    #[test]
    fn semantic_runinfo_fields_and_nonvolatile_headers_are_not_hidden() {
        let (ok, detail) = compare_pair(
            "json",
            r#"{"final_args":{"pass_only":false}}"#,
            r#"{"final_args":{"pass_only":true}}"#,
        );
        assert!(!ok);
        assert!(detail.contains("$/final_args/pass_only"));

        let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let (ok, detail) = compare_pair(
            "vcf",
            &("##bcftools_viewCommandExtra=legacy\n".to_string() + header),
            &("##bcftools_viewCommandExtra=rust\n".to_string() + header),
        );
        assert!(!ok);
        assert!(detail.contains("bcftools_viewCommandExtra"));
    }

    #[test]
    fn exact_bcftools_provenance_headers_are_volatile() {
        let stable = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tC\t.\tPASS\t.\n";
        let legacy = format!(
            "##bcftools_normVersion=1.9\n##bcftools_normCommand=norm -f legacy.fa; Date=old\n{stable}"
        );
        let rust = format!(
            "##bcftools_normVersion=1.22\n##bcftools_normCommand=norm -f rust.fa; Date=new\n{stable}"
        );
        let (ok, detail) = compare_pair("vcf", &legacy, &rust);
        assert!(ok, "{detail}");

        let legacy_merge = format!(
            "##bcftools_mergeVersion=1.17\n##bcftools_mergeCommand=merge /tmp/legacy; Date=old\n{stable}"
        );
        let rust_merge = format!(
            "##bcftools_mergeVersion=1.22\n##bcftools_mergeCommand=merge native; Date=new\n{stable}"
        );
        let (ok, detail) = compare_pair("vcf", &legacy_merge, &rust_merge);
        assert!(ok, "{detail}");

        let (ok, detail) = compare_pair(
            "vcf",
            &("##bcftools_normCommandExtra=legacy\n".to_string() + stable),
            &("##bcftools_normCommandExtra=rust\n".to_string() + stable),
        );
        assert!(!ok);
        assert!(detail.contains("bcftools_normCommandExtra"));
    }

    #[test]
    fn vcf_header_and_record_order_contract_is_preserved() {
        let header = "##fileformat=VCFv4.2\n##commandline=VALUE\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let first = "chr1\t1\t.\tA\tC\t.\tPASS\t.\n";
        let second = "chr1\t2\t.\tG\tT\t.\tPASS\t.\n";
        let (ok, detail) = compare_pair(
            "vcf",
            &(header.to_string() + first + second),
            &(header.replace("VALUE", "OTHER") + second + first),
        );
        assert!(!ok);
        assert!(detail.contains("different ordering"));
    }

    #[test]
    fn symmetric_and_empty_artifact_sets_fail_in_status_schema() {
        let dir = tempfile::tempdir().unwrap();
        let legacy = dir.path().join("legacy");
        let rust = dir.path().join("rust");
        fs::create_dir_all(&legacy).unwrap();
        fs::create_dir_all(&rust).unwrap();
        fs::write(rust.join("result.extra"), "unexpected").unwrap();
        let report = dir.path().join("diff.log");
        let status = dir.path().join("status.json");
        compare_outputs(&CompareOutputs {
            legacy_dir: &legacy,
            rust_dir: &rust,
            prefix: "result",
            case: "test",
            sample: "sample",
            report: &report,
            status: &status,
        })
        .unwrap();
        let text = fs::read_to_string(&status).unwrap();
        assert!(text.contains("\"reason\": \"only-in-rust\""));

        fs::remove_file(rust.join("result.extra")).unwrap();
        compare_outputs(&CompareOutputs {
            legacy_dir: &legacy,
            rust_dir: &rust,
            prefix: "result",
            case: "test",
            sample: "sample",
            report: &report,
            status: &status,
        })
        .unwrap();
        assert!(
            fs::read_to_string(status)
                .unwrap()
                .contains("no-matching-artifacts")
        );
    }

    #[test]
    fn roc_integrity_accepts_complete_rows_and_rejects_corruption() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("result.roc.all.csv.gz");
        let valid = "Type,Subtype,Subset,Filter,Genotype,QQ.Field,QQ,TRUTH.TP\nSNP,*,*,ALL,*,QUAL,1.000000,1\n";
        let file = fs::File::create(&path).unwrap();
        let mut writer = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        use std::io::Write;
        writer.write_all(valid.as_bytes()).unwrap();
        writer.finish().unwrap();
        assert_eq!(validate_roc(&path).unwrap(), 1);

        let invalid =
            "Type,Subtype,Subset,Filter,Genotype,QQ.Field,QQ\nSNP,,*,ALL,*,QUAL,1.0,wide\n";
        let file = fs::File::create(&path).unwrap();
        let mut writer = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        writer.write_all(invalid.as_bytes()).unwrap();
        writer.finish().unwrap();
        assert!(
            validate_roc(&path)
                .unwrap_err()
                .to_string()
                .contains("more fields")
        );

        let invalid_type =
            "Type,Subtype,Subset,Filter,Genotype,QQ.Field,QQ\nSNP_CORRUPT,*,*,ALL,*,QUAL,1.0\n";
        let file = fs::File::create(&path).unwrap();
        let mut writer = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        writer.write_all(invalid_type.as_bytes()).unwrap();
        writer.finish().unwrap();
        assert!(
            validate_roc(&path)
                .unwrap_err()
                .to_string()
                .contains("invalid variant type")
        );

        let header_only = "Type,Subtype,Subset,Filter,Genotype,QQ.Field,QQ\n";
        let file = fs::File::create(&path).unwrap();
        let mut writer = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        writer.write_all(header_only.as_bytes()).unwrap();
        writer.finish().unwrap();
        assert!(
            validate_roc(&path)
                .unwrap_err()
                .to_string()
                .contains("no ROC data rows")
        );
    }

    #[test]
    fn aggregate_report_sorts_rows_and_marks_invalid_status_failed() {
        let dir = tempfile::tempdir().unwrap();
        let good = dir.path().join("good.json");
        let bad = dir.path().join("bad.json");
        fs::write(
            &good,
            r#"{"case":"happy","sample":"z","total":2,"failed":0}"#,
        )
        .unwrap();
        fs::write(&bad, "not json").unwrap();
        let markdown = dir.path().join("report.md");
        let csv = dir.path().join("report.csv");
        aggregate_report(&AggregateReport {
            image: "oracle",
            hap_bin: "hap",
            markdown: &markdown,
            csv: &csv,
            inputs: &[good, bad.clone()],
        })
        .unwrap();
        let report = fs::read_to_string(markdown).unwrap();
        assert!(report.contains("- image: `oracle`"));
        assert!(report.contains("| happy | z | 2 | 0 | PASS |"));
        assert!(report.contains(&format!("| unknown | {} | 0 | 1 | FAIL |", bad.display())));
    }

    #[test]
    fn bcf_comparison_is_semantic_not_bytewise() {
        let dir = tempfile::tempdir().unwrap();
        let compressed = dir.path().join("legacy.bcf");
        let uncompressed = dir.path().join("rust.bcf");
        let headers = [
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=100>".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE".to_string(),
        ];
        let record = crate::vcf::RawVcfRecord::from_line(
            "chr1\t7\t.\tA\tC\t30\tPASS\t.\tGT\t0/1",
            &compressed,
        )
        .unwrap();
        crate::vcf::write_raw_vcf(&compressed, &headers, &[record]).unwrap();
        fs::write(
            &uncompressed,
            crate::bcf::read_uncompressed(&compressed).unwrap(),
        )
        .unwrap();
        assert_ne!(
            fs::read(&compressed).unwrap(),
            fs::read(&uncompressed).unwrap()
        );
        let (ok, detail) = compare_file("result.bcf", &compressed, &uncompressed);
        assert!(ok, "{}", detail.unwrap());
    }

    #[test]
    fn functional_bcf_csi_query_detects_a_corrupt_index() {
        let dir = tempfile::tempdir().unwrap();
        for side in ["legacy", "rust"] {
            let root = dir.path().join(side);
            fs::create_dir_all(&root).unwrap();
            let bcf = root.join("result.bcf");
            let headers = [
                "##fileformat=VCFv4.2".to_string(),
                "##contig=<ID=chr1,length=100>".to_string(),
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
            ];
            let record =
                crate::vcf::RawVcfRecord::from_line("chr1\t7\t.\tA\tC\t.\tPASS\t.", &bcf).unwrap();
            crate::vcf::write_raw_vcf(&bcf, &headers, &[record]).unwrap();
        }
        let left = dir.path().join("legacy/result.bcf.csi");
        let right = dir.path().join("rust/result.bcf.csi");
        assert!(compare_index(&left, &right).unwrap().is_empty());
        fs::write(&right, "not an index").unwrap();
        assert!(
            compare_index(&left, &right)
                .unwrap()
                .contains("rust index query failed")
        );
    }

    #[test]
    fn functional_bcf_csi_query_accepts_count_metadata() {
        let dir = tempfile::tempdir().unwrap();
        let left_bcf = write_test_indexed_bcf(&dir.path().join("legacy"));
        write_test_indexed_bcf(&dir.path().join("rust"));
        let left = left_bcf.with_extension("bcf.csi");
        let right = dir.path().join("rust/result.bcf.csi");

        // In a CSI metadata pseudo-bin, the second chunk stores mapped and
        // unmapped record counts. A common (mapped > 0, unmapped == 0) pair
        // is not an inverted virtual-offset interval.
        add_test_csi_count_metadata(&left, 1, 0);
        assert!(compare_index(&left, &right).unwrap().is_empty());
    }

    #[test]
    fn functional_bcf_csi_query_rejects_reversed_data_chunk() {
        let dir = tempfile::tempdir().unwrap();
        let left_bcf = write_test_indexed_bcf(&dir.path().join("legacy"));
        let right_bcf = write_test_indexed_bcf(&dir.path().join("rust"));
        let left = left_bcf.with_extension("bcf.csi");
        let right = right_bcf.with_extension("bcf.csi");

        reverse_test_csi_data_chunk(&right);
        let detail = compare_index(&left, &right).unwrap();
        assert!(detail.contains("rust index query failed"), "{detail}");
        assert!(detail.contains("chunk start exceeds its end"), "{detail}");
    }

    #[test]
    fn functional_bcf_csi_query_detects_omitted_records() {
        let dir = tempfile::tempdir().unwrap();
        for side in ["legacy", "rust"] {
            let root = dir.path().join(side);
            fs::create_dir_all(&root).unwrap();
            let bcf = root.join("result.bcf");
            let headers = [
                "##fileformat=VCFv4.2".to_string(),
                "##contig=<ID=chr1,length=100>".to_string(),
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
            ];
            let records = [
                crate::vcf::RawVcfRecord::from_line("chr1\t7\tfirst\tA\tC\t.\tPASS\t.", &bcf)
                    .unwrap(),
                crate::vcf::RawVcfRecord::from_line("chr1\t11\tsecond\tG\tT\t.\tPASS\t.", &bcf)
                    .unwrap(),
            ];
            crate::vcf::write_raw_vcf(&bcf, &headers, &records).unwrap();
        }
        let left = dir.path().join("legacy/result.bcf.csi");
        let right = dir.path().join("rust/result.bcf.csi");
        assert!(compare_index(&left, &right).unwrap().is_empty());

        // This is a structurally valid CSI with the right reference count but
        // no chunks. Structural parsing alone would accept it; indexed record
        // retrieval must prove that both source records are missing.
        crate::bcf::write_csi(&right, &[None]).unwrap();
        let detail = compare_index(&left, &right).unwrap();
        assert!(detail.contains("rust index query failed"), "{detail}");
        assert!(
            detail.contains("returned 0 of 2 source records"),
            "{detail}"
        );
    }

    #[test]
    fn functional_tabix_query_detects_a_corrupt_index() {
        if Command::new("tabix").arg("--version").output().is_err() {
            return;
        }
        let dir = tempfile::tempdir().unwrap();
        for side in ["legacy", "rust"] {
            let root = dir.path().join(side);
            fs::create_dir_all(&root).unwrap();
            let vcf = root.join("result.vcf.gz");
            crate::vcf::write_indexed_vcf(
                &vcf,
                &[
                    "##fileformat=VCFv4.2".to_string(),
                    "##contig=<ID=chr1,length=100>".to_string(),
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
                ],
                ["chr1\t7\t.\tA\tC\t.\tPASS\t."],
            )
            .unwrap();
        }
        let left = dir.path().join("legacy/result.vcf.gz.tbi");
        let right = dir.path().join("rust/result.vcf.gz.tbi");
        assert!(compare_index(&left, &right).unwrap().is_empty());
        fs::write(&right, "not an index").unwrap();
        assert!(
            compare_index(&left, &right)
                .unwrap()
                .contains("rust index query failed")
        );
    }

    #[test]
    fn functional_tabix_query_rejects_an_index_omitting_a_source_contig() {
        if Command::new("tabix").arg("--version").output().is_err() {
            return;
        }
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("result.vcf.gz");
        let headers = [
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=100>".to_string(),
            "##contig=<ID=chr2,length=100>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        crate::vcf::write_indexed_vcf(&vcf, &headers, ["chr1\t7\t.\tA\tC\t.\tPASS\t."]).unwrap();
        let index = dir.path().join("result.vcf.gz.tbi");
        let stale = fs::read(&index).unwrap();
        crate::vcf::write_indexed_vcf(
            &vcf,
            &headers,
            [
                "chr1\t7\t.\tA\tC\t.\tPASS\t.",
                "chr2\t91\t.\tG\tT\t.\tPASS\t.",
            ],
        )
        .unwrap();
        fs::write(&index, stale).unwrap();
        assert!(
            tabix_range_query(&index)
                .unwrap_err()
                .to_string()
                .contains("index omits source references: chr2")
        );
    }
}
