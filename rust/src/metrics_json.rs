use anyhow::{Context, Result};
use flate2::Compression;
use flate2::write::GzEncoder;
use std::collections::BTreeMap;
use std::fs;
use std::io::Write;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

/// Inputs that differ between invocations of the legacy `hap.py` parser.
///
/// The remaining `final_args` entries are legacy defaults. Keeping those
/// defaults here makes their JSON scalar types explicit instead of relying on
/// stringly-typed maps assembled by callers.
pub struct CompareRunArgs<'a> {
    pub truth: &'a str,
    pub query: &'a str,
    pub reference: &'a str,
    pub reports_prefix: &'a str,
    pub annotation_type: Option<&'a str>,
    pub pass_only: bool,
    pub preprocessing_truth: bool,
    pub preprocessing_leftshift: bool,
    pub preprocessing_decompose: bool,
    pub regions_bedfile: Option<&'a str>,
    pub targets_bedfile: Option<&'a str>,
    pub fp_bedfile: Option<&'a str>,
    pub locations: Option<&'a str>,
    pub threads: usize,
    pub strat_tsv: Option<&'a str>,
    pub scratch_prefix: Option<&'a str>,
    pub keep_scratch: bool,
    pub bcf: bool,
    pub ci_alpha: f64,
    pub convert_gvcf_query: bool,
    pub convert_gvcf_to_vcf: bool,
    pub convert_gvcf_truth: bool,
    pub do_roc: bool,
    pub engine: &'a str,
    pub engine_scmp_distance: usize,
    pub engine_vcfeval: &'a str,
    pub engine_vcfeval_template: Option<&'a str>,
    pub filter_nonref: bool,
    pub filters_only: Option<&'a str>,
    pub fixchr: Option<bool>,
    pub fp_adjust_conf: bool,
    pub gender: &'a str,
    pub hb_expand: usize,
    pub logfile: Option<&'a str>,
    pub max_enum: usize,
    pub no_hc: bool,
    pub output_vtc: bool,
    pub preprocess_window: usize,
    pub preprocessing_norm: bool,
    pub preserve_info: bool,
    pub quiet: bool,
    pub roc: &'a str,
    pub roc_delta: f64,
    pub roc_filter: Option<&'a str>,
    pub roc_regions: &'a [String],
    pub somatic: bool,
    pub somatic_mode: Option<&'a str>,
    pub strat_fixchr: bool,
    pub strat_regions: &'a [String],
    pub usefiltered_truth: bool,
    pub verbose: bool,
    pub window: usize,
    pub write_counts: bool,
    pub write_json: bool,
    pub write_vcf: bool,
}

pub fn write_compare_runinfo(
    path: &Path,
    commandline: &str,
    args: &CompareRunArgs<'_>,
) -> Result<()> {
    create_parent(path)?;

    let mut body = String::new();
    body.push('{');
    body.push_str(&runinfo_platform_keys());
    push_field(&mut body, "name", &json_string("hap.py"));
    push_field(
        &mut body,
        "runInfo",
        &format!(
            "[{{\"value\":{},\"key\":\"commandline\"}}]",
            json_string(commandline)
        ),
    );
    push_field(&mut body, "version", &json_string(""));
    push_field(&mut body, "metadata", &metadata_json("hap.py", commandline));
    body.push_str("\"final_args\":");
    body.push_str(&final_args_json(args));
    body.push('}');

    fs::write(path, body).with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn final_args_json(args: &CompareRunArgs<'_>) -> String {
    let mut fields = Vec::with_capacity(57);
    let mut add = |key: &str, value: String| fields.push(format!("{}:{value}", json_string(key)));
    add(
        "_vcfs",
        format!("[{},{}]", json_string(args.truth), json_string(args.query)),
    );
    add("bcf", json_bool(args.bcf));
    add("ci_alpha", json_float(args.ci_alpha));
    add("convert_gvcf_query", json_bool(args.convert_gvcf_query));
    add("convert_gvcf_to_vcf", json_bool(args.convert_gvcf_to_vcf));
    add("convert_gvcf_truth", json_bool(args.convert_gvcf_truth));
    add("delete_scratch", json_bool(!args.keep_scratch));
    add("do_roc", json_bool(args.do_roc));
    add("engine", json_string(args.engine));
    add(
        "engine_scmp_distance",
        args.engine_scmp_distance.to_string(),
    );
    add("engine_vcfeval", json_string(args.engine_vcfeval));
    add(
        "engine_vcfeval_template",
        json_optional_string(args.engine_vcfeval_template),
    );
    add("filter_nonref", json_bool(args.filter_nonref));
    add("filters_only", json_string(args.filters_only.unwrap_or("")));
    add("fixchr", json_optional_bool(args.fixchr));
    add("force_interactive", json_bool(true));
    add("fp_bedfile", json_optional_string(args.fp_bedfile));
    add("gender", json_string(args.gender));
    add("hb_expand", args.hb_expand.to_string());
    add("locations", json_optional_string(args.locations));
    add("logfile", json_optional_string(args.logfile));
    add("max_enum", args.max_enum.to_string());
    add("no_hc", json_bool(args.no_hc));
    add("output_vtc", json_bool(args.output_vtc));
    add("pass_only", json_bool(args.pass_only));
    add("preprocess_window", args.preprocess_window.to_string());
    add(
        "preprocessing_decompose",
        json_bool(args.preprocessing_decompose),
    );
    add(
        "preprocessing_leftshift",
        json_bool(args.preprocessing_leftshift),
    );
    add("preprocessing_norm", json_bool(args.preprocessing_norm));
    add("preprocessing_truth", json_bool(args.preprocessing_truth));
    add(
        "preprocessing_truth_confregions",
        json_bool(args.fp_adjust_conf),
    );
    add("preserve_info", json_bool(args.preserve_info));
    add("quiet", json_bool(args.quiet));
    add("ref", json_string(args.reference));
    add(
        "regions_bedfile",
        json_optional_string(args.regions_bedfile),
    );
    add("reports_prefix", json_string(args.reports_prefix));
    add("roc", json_string(args.roc));
    add("roc_delta", json_float(args.roc_delta));
    add(
        "roc_filter",
        args.roc_filter
            .map(json_string)
            .unwrap_or_else(|| json_bool(false)),
    );
    add("roc_regions", json_string_array(args.roc_regions, "*"));
    add("scratch_prefix", json_optional_string(args.scratch_prefix));
    add(
        "somatic_allele_conversion",
        args.somatic_mode
            .map(json_string)
            .unwrap_or_else(|| json_bool(args.somatic)),
    );
    add(
        "strat_fixchr",
        if args.strat_fixchr {
            json_bool(true)
        } else {
            "null".to_string()
        },
    );
    add("strat_regions", json_string_array(args.strat_regions, ""));
    add("strat_tsv", json_optional_string(args.strat_tsv));
    add(
        "targets_bedfile",
        json_optional_string(args.targets_bedfile),
    );
    add("threads", args.threads.to_string());
    add("type", json_optional_string(args.annotation_type));
    add("usefiltered_truth", json_bool(args.usefiltered_truth));
    add("vcf1", json_string(args.truth));
    add("vcf2", json_string(args.query));
    add("verbose", json_bool(args.verbose));
    add("version", json_bool(false));
    add("window", args.window.to_string());
    add("write_counts", json_bool(args.write_counts));
    add("write_json", json_bool(args.write_json));
    add("write_vcf", json_bool(args.write_vcf));
    debug_assert_eq!(fields.len(), 57);
    format!("{{{}}}", fields.join(","))
}

fn runinfo_platform_keys() -> String {
    let mut out = String::new();
    push_field(&mut out, "dist", &json_string(""));
    push_field(&mut out, "python_version", &json_string(""));
    push_field(&mut out, "timestamp", &json_string(&timestamp_now()));
    push_field(&mut out, "python_implementation", &json_string(""));
    push_field(&mut out, "uname", &json_string(""));
    push_field(&mut out, "python_prefix", &json_string(""));
    push_field(&mut out, "environment", "{}");
    push_field(&mut out, "mac_ver", &json_string(""));
    out
}

fn metadata_json(module: &str, commandline: &str) -> String {
    let executable = commandline.split_whitespace().next().unwrap_or("hap.py");
    format!(
        "{{\"required\":{{\"version\":\"\",\"id\":\"haplotypes\",\"module\":{},\"description\":{}}}}}",
        json_string(module),
        json_string(&format!(
            "{executable} generated this JSON file via command line {commandline}"
        ))
    )
}

fn timestamp_now() -> String {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_secs().to_string())
        .unwrap_or_else(|_| "0".to_string())
}

pub fn write_metrics_gz(
    path: &Path,
    name: &str,
    commandline: &str,
    tables: &[(&str, &str, &Path)],
) -> Result<()> {
    write_metrics_gz_for_module(path, name, "hap.py", commandline, tables)
}

pub fn write_metrics_gz_for_module(
    path: &Path,
    name: &str,
    module: &str,
    commandline: &str,
    tables: &[(&str, &str, &Path)],
) -> Result<()> {
    write_metrics_gz_for_module_with_indices(path, name, module, commandline, tables, None)
}

pub fn write_metrics_gz_for_module_with_indices(
    path: &Path,
    name: &str,
    module: &str,
    commandline: &str,
    tables: &[(&str, &str, &Path)],
    indices: Option<&BTreeMap<String, Vec<usize>>>,
) -> Result<()> {
    create_parent(path)?;
    let file =
        fs::File::create(path).with_context(|| format!("failed to create {}", path.display()))?;
    let mut encoder = GzEncoder::new(file, Compression::default());
    encoder.write_all(
        metrics_json_for_module_with_indices(name, module, commandline, tables, indices)?
            .as_bytes(),
    )?;
    encoder.finish()?;
    Ok(())
}

pub fn write_metrics_json(
    path: &Path,
    name: &str,
    commandline: &str,
    tables: &[(&str, &str, &Path)],
) -> Result<()> {
    create_parent(path)?;
    fs::write(path, metrics_json(name, commandline, tables)?)
        .with_context(|| format!("failed to write {}", path.display()))?;
    Ok(())
}

fn metrics_json(name: &str, commandline: &str, tables: &[(&str, &str, &Path)]) -> Result<String> {
    metrics_json_for_module(name, "hap.py", commandline, tables)
}

fn metrics_json_for_module(
    name: &str,
    module: &str,
    commandline: &str,
    tables: &[(&str, &str, &Path)],
) -> Result<String> {
    metrics_json_for_module_with_indices(name, module, commandline, tables, None)
}

fn metrics_json_for_module_with_indices(
    name: &str,
    module: &str,
    commandline: &str,
    tables: &[(&str, &str, &Path)],
    indices: Option<&BTreeMap<String, Vec<usize>>>,
) -> Result<String> {
    let rendered_tables = tables
        .iter()
        .map(|(id, label, path)| {
            table_json_with_indices(
                id,
                label,
                path,
                indices
                    .and_then(|tables| tables.get(*id))
                    .map(Vec::as_slice),
            )
        })
        .collect::<Result<Vec<_>>>()?
        .join(",");
    Ok(format!(
        "{{\"runInfo\":[{{\"value\":{},\"key\":\"commandline\"}}],\"metrics\":[{rendered_tables}],\"version\":\"\",\"sampleInfo\":[],\"name\":{},\"parameters\":[],\"timestamp\":{},\"metadata\":{}}}",
        json_string(commandline),
        json_string(name),
        json_string(&timestamp_now()),
        metadata_json(module, commandline),
    ))
}

fn create_parent(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

#[cfg(test)]
fn table_json(id: &str, label: &str, path: &Path) -> Result<String> {
    table_json_with_indices(id, label, path, None)
}

fn table_json_with_indices(
    id: &str,
    label: &str,
    path: &Path,
    indices: Option<&[usize]>,
) -> Result<String> {
    let text = crate::vcf::read_text(path)
        .with_context(|| format!("failed to read {}", path.display()))?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("missing CSV header in {}", path.display()))?
        .split(',')
        .map(str::to_string)
        .collect::<Vec<_>>();
    let rows = lines
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.split(',').map(str::to_string).collect::<Vec<_>>())
        .collect::<Vec<_>>();

    let mut out = String::new();
    out.push_str("{\"data\":[");
    out.push_str(&index_column_json(rows.len(), indices));
    for (column, name) in header.iter().enumerate() {
        out.push(',');
        let values = rows
            .iter()
            .map(|row| row.get(column).map(String::as_str).unwrap_or(""))
            .collect::<Vec<_>>();
        let kind = legacy_column_type(id, name);
        out.push_str(&column_json(name, name, kind, &values));
    }
    out.push_str("],\"properties\":[],\"type\":\"Table\",\"id\":");
    out.push_str(&json_string(id));
    out.push_str(",\"label\":");
    out.push_str(&json_string(label));
    out.push('}');
    Ok(out)
}

fn index_column_json(count: usize, indices: Option<&[usize]>) -> String {
    let values = indices
        .filter(|values| values.len() == count)
        .map(|values| values.to_vec())
        .unwrap_or_else(|| (0..count).collect())
        .into_iter()
        .map(|index| index.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!("{{\"values\":[{values}],\"type\":\"string\",\"id\":\"types\",\"label\":\"types\"}}")
}

fn legacy_column_type(table_id: &str, column: &str) -> &'static str {
    if matches!(
        column,
        "TRUTH.TOTAL"
            | "TRUTH.TP"
            | "TRUTH.FN"
            | "QUERY.TOTAL"
            | "QUERY.TP"
            | "QUERY.FP"
            | "QUERY.UNK"
            | "FP.gt"
            | "FP.al"
    ) || (column == "Subset.Size" && table_id.starts_with("roc.Locations."))
    {
        "int64"
    } else if column.starts_with("METRIC.")
        || column.ends_with(".TiTv_ratio")
        || column.ends_with(".het_hom_ratio")
    {
        "double"
    } else {
        "string"
    }
}

fn column_json(id: &str, label: &str, kind: &str, values: &[&str]) -> String {
    let rendered = values
        .iter()
        .map(|value| match kind {
            "int64" => render_int(value),
            "double" => render_double(value),
            _ => json_string(value),
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"values\":[{rendered}],\"type\":{},\"id\":{},\"label\":{}}}",
        json_string(kind),
        json_string(id),
        json_string(label)
    )
}

fn render_int(value: &str) -> String {
    if value.is_empty() || value == "." {
        "null".to_string()
    } else {
        value
            .parse::<i64>()
            .map(|number| number.to_string())
            .unwrap_or_else(|_| "null".to_string())
    }
}

fn render_double(value: &str) -> String {
    if value.is_empty() || value == "." {
        return "null".to_string();
    }
    let Ok(number) = value.parse::<f64>() else {
        return "null".to_string();
    };
    if !number.is_finite() {
        return "null".to_string();
    }
    let mut rendered = number.to_string();
    if !rendered.contains(['.', 'e', 'E']) {
        rendered.push_str(".0");
    }
    rendered
}

fn push_field(out: &mut String, key: &str, value: &str) {
    out.push_str(&json_string(key));
    out.push(':');
    out.push_str(value);
    out.push(',');
}

fn json_optional_string(value: Option<&str>) -> String {
    value.map(json_string).unwrap_or_else(|| "null".to_string())
}

fn json_optional_bool(value: Option<bool>) -> String {
    value.map(json_bool).unwrap_or_else(|| "null".to_string())
}

fn json_float(value: f64) -> String {
    let mut rendered = value.to_string();
    if !rendered.contains(['.', 'e', 'E']) {
        rendered.push_str(".0");
    }
    rendered
}

fn json_string_array(values: &[String], default: &str) -> String {
    let rendered = if values.is_empty() && !default.is_empty() {
        vec![json_string(default)]
    } else {
        values.iter().map(|value| json_string(value)).collect()
    };
    format!("[{}]", rendered.join(","))
}

fn json_bool(value: bool) -> String {
    if value { "true" } else { "false" }.to_string()
}

fn json_string(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len() + 2);
    escaped.push('"');
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            '\u{08}' => escaped.push_str("\\b"),
            '\u{0c}' => escaped.push_str("\\f"),
            character if character.is_control() => {
                escaped.push_str(&format!("\\u{:04x}", character as u32));
            }
            character => escaped.push(character),
        }
    }
    escaped.push('"');
    escaped
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn metrics_table_uses_legacy_column_types_and_nulls() {
        let dir = tempfile::tempdir().unwrap();
        let csv = dir.path().join("table.csv");
        fs::write(
            &csv,
            "Type,Filter,TRUTH.TOTAL,METRIC.Recall,TRUTH.TOTAL.ti,TRUTH.TOTAL.tv,TRUTH.TOTAL.TiTv_ratio,Subset.Size,Subset.IS_CONF.Size\nSNP,ALL,3,0,3.000000,1.000000,.,42,21.000000\nINDEL,PASS,0,,,,1.5,42,\n",
        )
        .unwrap();

        let json = table_json("all.metrics", "all.metrics", &csv).unwrap();
        assert!(json.contains("\"values\":[0,1],\"type\":\"string\",\"id\":\"types\""));
        assert!(json.contains("\"values\":[3,0],\"type\":\"int64\",\"id\":\"TRUTH.TOTAL\""));
        assert!(
            json.contains("\"values\":[0.0,null],\"type\":\"double\",\"id\":\"METRIC.Recall\"")
        );
        assert!(json.contains(
            "\"values\":[\"3.000000\",\"\"],\"type\":\"string\",\"id\":\"TRUTH.TOTAL.ti\""
        ));
        assert!(json.contains(
            "\"values\":[\"1.000000\",\"\"],\"type\":\"string\",\"id\":\"TRUTH.TOTAL.tv\""
        ));
        assert!(json.contains(
            "\"values\":[null,1.5],\"type\":\"double\",\"id\":\"TRUTH.TOTAL.TiTv_ratio\""
        ));
        assert!(
            json.contains("\"values\":[\"42\",\"42\"],\"type\":\"string\",\"id\":\"Subset.Size\"")
        );
        assert!(json.contains(
            "\"values\":[\"21.000000\",\"\"],\"type\":\"string\",\"id\":\"Subset.IS_CONF.Size\""
        ));

        let locations = table_json("roc.Locations.SNP", "roc.Locations.SNP", &csv).unwrap();
        assert!(locations.contains("\"values\":[42,42],\"type\":\"int64\",\"id\":\"Subset.Size\""));
        assert!(locations.contains(
            "\"values\":[\"21.000000\",\"\"],\"type\":\"string\",\"id\":\"Subset.IS_CONF.Size\""
        ));

        let indexed =
            table_json_with_indices("all.metrics", "all.metrics", &csv, Some(&[9, 3])).unwrap();
        assert!(indexed.contains("\"values\":[9,3],\"type\":\"string\",\"id\":\"types\""));
    }

    #[test]
    fn runinfo_final_args_preserve_legacy_json_scalar_types() {
        let args = CompareRunArgs {
            truth: "truth.vcf.gz",
            query: "query.vcf.gz",
            reference: "ref.fa",
            reports_prefix: "result",
            annotation_type: Some("xcmp"),
            pass_only: true,
            preprocessing_truth: true,
            preprocessing_leftshift: false,
            preprocessing_decompose: false,
            regions_bedfile: None,
            targets_bedfile: Some("targets.bed"),
            fp_bedfile: None,
            locations: Some("chr1:1-2"),
            threads: 3,
            strat_tsv: None,
            scratch_prefix: None,
            keep_scratch: false,
            bcf: false,
            ci_alpha: 0.0,
            convert_gvcf_query: false,
            convert_gvcf_to_vcf: false,
            convert_gvcf_truth: false,
            do_roc: true,
            engine: "xcmp",
            engine_scmp_distance: 30,
            engine_vcfeval: "rtg",
            engine_vcfeval_template: None,
            filter_nonref: false,
            filters_only: None,
            fixchr: None,
            fp_adjust_conf: true,
            gender: "auto",
            hb_expand: 30,
            logfile: None,
            max_enum: 16_768,
            no_hc: false,
            output_vtc: false,
            preprocess_window: 10_000,
            preprocessing_norm: false,
            preserve_info: false,
            quiet: false,
            roc: "QUAL",
            roc_delta: 0.5,
            roc_filter: None,
            roc_regions: &[],
            somatic: false,
            somatic_mode: None,
            strat_fixchr: false,
            strat_regions: &[],
            usefiltered_truth: false,
            verbose: false,
            window: 50,
            write_counts: true,
            write_json: true,
            write_vcf: true,
        };
        let json = final_args_json(&args);

        assert!(json.contains("\"_vcfs\":[\"truth.vcf.gz\",\"query.vcf.gz\"]"));
        assert!(json.contains("\"pass_only\":true"));
        assert!(json.contains("\"preprocessing_truth\":true"));
        assert!(json.contains("\"preprocessing_leftshift\":false"));
        assert!(json.contains("\"preprocessing_decompose\":false"));
        assert!(json.contains("\"regions_bedfile\":null"));
        assert!(json.contains("\"targets_bedfile\":\"targets.bed\""));
        assert!(json.contains("\"threads\":3"));
        assert!(json.contains("\"ci_alpha\":0.0"));
        assert!(json.contains("\"roc_regions\":[\"*\"]"));
        assert_eq!(json.matches("\":").count(), 57);
    }

    #[test]
    fn metrics_document_matches_legacy_metadata_container_shapes() {
        let dir = tempfile::tempdir().unwrap();
        let csv = dir.path().join("summary.csv");
        fs::write(&csv, "Type,TRUTH.TOTAL\nSNP,1\n").unwrap();
        let json = metrics_json(
            "hap.py.comparison",
            "hap germline truth query",
            &[("summary.metrics", "summary.metrics", csv.as_path())],
        )
        .unwrap();

        assert!(json.contains("\"parameters\":[]"));
        assert!(json.contains("\"sampleInfo\":[]"));
        assert!(json.contains("\"name\":\"hap.py.comparison\""));
        assert!(json.contains(
            "\"metadata\":{\"required\":{\"version\":\"\",\"id\":\"haplotypes\",\"module\":\"hap.py\""
        ));
    }
}
