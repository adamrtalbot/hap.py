//! Cohesive preprocessing responsibility.

use super::*;

pub(super) fn validate_record_reference(
    record: &RawVcfRecord,
    reference_sequences: &std::collections::BTreeMap<String, String>,
) -> Result<()> {
    let reference = reference_sequences
        .get(&record.chrom)
        .ok_or_else(|| anyhow::anyhow!("reference contig {} not found", record.chrom))?;
    // DNA references are ASCII, so byte slicing is safe and avoids allocating a
    // 48 M-entry Vec<char> per record on chr-scale references.
    let bases = reference.as_bytes();
    let end = record.end_pos();
    if record.pos == 0 || end > bases.len() {
        bail!(
            "record {}:{} extends beyond the reference sequence",
            record.chrom,
            record.pos
        );
    }
    let observed = &bases[record.pos - 1..end];
    // Reference FASTA may contain soft-masked lowercase bases (repeat regions).
    // Legacy hap.py tolerates case differences between the VCF REF and the
    // reference sequence — any N in either side also matches anything.
    if !ref_bytes_equal(observed, record.ref_allele.as_bytes()) {
        bail!(
            "record {}:{} REF={} does not match reference {}",
            record.chrom,
            record.pos,
            record.ref_allele,
            String::from_utf8_lossy(observed)
        );
    }
    Ok(())
}

/// INFO keys whose semantics assume a specific allele layout. They go stale as
/// soon as preprocessing splits or decomposes variants, so the legacy C++
/// `preprocess` binary (via `VariantCallsOnly`) drops them on output. Keep this
/// list tight — any tag we strip permanently is one a consumer can't rely on.
pub(super) const STALE_INFO_KEYS: &[&str] = &["AC", "AN", "MLEAC", "MLEAF"];

/// Reorder `record.info` so its tags appear alphabetically, matching legacy
/// `std::map`-ordered serialisation. Flag-only entries (no `=`) keep their
/// position in the sort by bare tag name. Numeric values that round-trip
/// through htslib's float/int parsing are also collapsed to canonical form
/// (e.g. `-0` → `0`, matching legacy's `bcf_update_info_float` write path
/// which loses the sign on negative zero).
pub(super) fn sort_info_keys(record: &mut RawVcfRecord) {
    if record.info.is_empty() || record.info == "." {
        return;
    }
    let mut entries: Vec<String> = record
        .info
        .split(';')
        .map(canonicalise_info_entry)
        .collect();
    entries.sort_by_key(|entry| {
        let key = entry.split('=').next().unwrap_or(entry);
        key.to_string()
    });
    record.info = entries.join(";");
}

/// Collapse `-0`, `-0.0`, etc. in INFO values to `0`. Legacy hap.py reads
/// every INFO value through htslib's typed parsers (`bcf_update_info_float`
/// / `bcf_update_info_int32`) which lose the sign on negative zero before
/// re-emitting the record. Mirroring that here on the main code path keeps
/// our output byte-identical without a diff-only shim.
pub(super) fn canonicalise_info_entry(entry: &str) -> String {
    let Some((key, value)) = entry.split_once('=') else {
        return entry.to_string();
    };
    let normalised: Vec<String> = value
        .split(',')
        .map(|part| {
            // Only collapse a value that *parses* as a number whose float
            // representation is exactly zero. Leaves non-numeric strings
            // (e.g. `set=variant2`, `culprit=FS`) untouched.
            if let Ok(f) = part.parse::<f64>()
                && f == 0.0
                && (part.starts_with('-') || part.contains('-'))
            {
                // Preserve the original integer / float visual shape
                // (e.g. `-0` → `0`, `-0.0` → `0`) so downstream byte
                // comparison sees what htslib would emit.
                return if part.contains('.') || part.contains('e') || part.contains('E') {
                    // Keep float shape — bcftools emits "0" for any
                    // signed-zero float regardless of original
                    // precision; mirror that.
                    "0".to_string()
                } else {
                    "0".to_string()
                };
            }
            part.to_string()
        })
        .collect();
    format!("{key}={}", normalised.join(","))
}

/// Reduce every sample's PL cell to its last comma-separated value. Legacy
/// hap.py's `preprocess` stores PL as a scalar int per sample (see
/// `VariantWriter.cpp:563`) and htslib serialises that as the final value only.
pub(super) fn collapse_pl_to_last_value(record: &mut RawVcfRecord) {
    let Some(format) = &record.format else {
        return;
    };
    let Some(pl_index) = format.split(':').position(|field| field == "PL") else {
        return;
    };
    for sample in &mut record.samples {
        let mut fields: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        if let Some(cell) = fields.get_mut(pl_index)
            && let Some(last) = cell.rsplit(',').next()
        {
            *cell = last.to_string();
        }
        *sample = fields.join(":");
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum ScalarType {
    Integer = 0,
    Float = 1,
    String = 2,
}

/// Classify a per-sample value as integer, float, or string based on the
/// first non-missing scalar. Mirrors `VariantWriter.cpp`'s runtime dispatch
/// (`v.asInt()` / `v.isNumeric()` / `v.asString()`).
pub(super) fn classify_value(value: &str) -> ScalarType {
    if value.is_empty() || value == "." {
        return ScalarType::Integer;
    }
    // Multi-valued entries (`3279,0,716`, `32,126`) classify on the first
    // numeric element — mirrors legacy's int/float scalar decision.
    let first = value.split(',').next().unwrap_or(value);
    if first == "." || first.is_empty() {
        return ScalarType::Integer;
    }
    if first.parse::<i64>().is_ok() {
        ScalarType::Integer
    } else if first.parse::<f64>().is_ok() {
        ScalarType::Float
    } else {
        ScalarType::String
    }
}

/// Reorder FORMAT fields and sample columns so that legacy's type-bucketed
/// canonical order is used: GT, AD, ADO, DP come first (in that order);
/// everything else is partitioned into integer / float / string buckets by
/// runtime-inferred type and each bucket is sorted alphabetically.
pub(super) fn reorder_format_fields(record: &mut RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    let field_names: Vec<String> = format.split(':').map(|s| s.to_string()).collect();

    // Parse each sample into its fields and pad to FORMAT width with ".".
    let original_samples: Vec<Vec<String>> = record
        .samples
        .iter()
        .map(|s| {
            let mut cells: Vec<String> = s.split(':').map(|x| x.to_string()).collect();
            while cells.len() < field_names.len() {
                cells.push(".".to_string());
            }
            cells
        })
        .collect();

    const FIXED: &[&str] = &["GT", "AD", "ADO", "DP"];

    let (fixed_fields, rest_fields): (Vec<_>, Vec<_>) = field_names
        .iter()
        .enumerate()
        .partition(|(_, name)| FIXED.contains(&name.as_str()));

    // Preserve GT/AD/ADO/DP order exactly as listed in FIXED.
    let mut ordered: Vec<(usize, String)> = FIXED
        .iter()
        .filter_map(|target| {
            fixed_fields
                .iter()
                .find(|(_, name)| name.as_str() == *target)
                .map(|(idx, name)| (*idx, name.to_string()))
        })
        .collect();

    // Type-bucket the remaining fields. Use the first sample's value to pick.
    let probe_sample = original_samples.first();
    let mut tagged: Vec<(ScalarType, String, usize)> = rest_fields
        .iter()
        .map(|(idx, name)| {
            let value = probe_sample
                .and_then(|s| s.get(*idx))
                .map(String::as_str)
                .unwrap_or(".");
            (classify_value(value), name.to_string(), *idx)
        })
        .collect();
    tagged.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    ordered.extend(tagged.into_iter().map(|(_, name, idx)| (idx, name)));

    let new_format = ordered
        .iter()
        .map(|(_, name)| name.clone())
        .collect::<Vec<_>>()
        .join(":");
    let new_samples: Vec<String> = original_samples
        .iter()
        .map(|sample| {
            ordered
                .iter()
                .map(|(original_idx, _)| {
                    sample
                        .get(*original_idx)
                        .cloned()
                        .unwrap_or_else(|| ".".to_string())
                })
                .collect::<Vec<_>>()
                .join(":")
        })
        .collect();

    record.format = Some(new_format);
    record.samples = new_samples;
}

pub(super) fn canonicalize_multi_allelic_order(record: &mut RawVcfRecord) {
    let alts: Vec<String> = record.alt_allele.split(',').map(str::to_string).collect();
    if alts.len() < 2 || alts.iter().any(|alt| alt.starts_with('<') || alt == "*") {
        return;
    }
    let mut ordered: Vec<(usize, String)> = alts.iter().cloned().enumerate().collect();
    ordered.sort_by(|left, right| left.1.len().cmp(&right.1.len()).then(left.1.cmp(&right.1)));
    let mut mapping = vec![0usize; alts.len() + 1];
    for (new_offset, (old_offset, _)) in ordered.iter().enumerate() {
        mapping[old_offset + 1] = new_offset + 1;
    }
    let keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect())
        .unwrap_or_default();
    let gt_index = keys.iter().position(|key| *key == "GT");
    let ad_index = keys.iter().position(|key| *key == "AD");
    for sample in &mut record.samples {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        if let Some(index) = gt_index
            && let Some(gt) = cells.get_mut(index)
        {
            let remapped = remap_gt(gt, &mapping);
            if remapped.contains('/') {
                let mut alleles: Vec<&str> = remapped.split('/').collect();
                if alleles.len() == 2
                    && alleles[0] != alleles[1]
                    && alleles
                        .iter()
                        .all(|allele| *allele != "0" && *allele != ".")
                {
                    alleles.sort_by(|left, right| right.cmp(left));
                    *gt = alleles.join("/");
                } else {
                    *gt = remapped;
                }
            } else {
                *gt = remapped;
            }
        }
        if let Some(index) = ad_index
            && let Some(ad) = cells.get_mut(index)
        {
            let depths: Vec<&str> = ad.split(',').collect();
            if depths.len() == alts.len() + 1 {
                let mut reordered = vec![depths[0]];
                reordered.extend(ordered.iter().map(|(old_offset, _)| depths[old_offset + 1]));
                *ad = reordered.join(",");
            }
        }
        *sample = cells.join(":");
    }
    record.alt_allele = ordered
        .into_iter()
        .map(|(_, alt)| alt)
        .collect::<Vec<_>>()
        .join(",");
}

pub(super) fn string_format_fields(headers: &[String]) -> BTreeSet<String> {
    headers
        .iter()
        .filter(|header| header.starts_with("##FORMAT=<") && header.contains("Type=String"))
        .filter_map(|header| {
            header
                .strip_prefix("##FORMAT=<ID=")
                .and_then(|tail| tail.split_once(',').map(|(id, _)| id.to_string()))
        })
        .collect()
}

pub(super) fn blank_secondary_sample_annotations(
    record: &mut RawVcfRecord,
    bcf_output: bool,
    string_fields: &BTreeSet<String>,
) {
    let keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect())
        .unwrap_or_default();
    for sample in record.samples.iter_mut().skip(1) {
        let mut cells: Vec<String> = sample.split(':').map(str::to_string).collect();
        let homref = keys
            .iter()
            .position(|key| *key == "GT")
            .and_then(|index| cells.get(index))
            .is_some_and(|gt| gt == "0/0" || gt == "0|0");
        let missing_gt = keys
            .iter()
            .position(|key| *key == "GT")
            .and_then(|index| cells.get(index))
            .is_some_and(|gt| {
                gt.chars()
                    .all(|character| matches!(character, '.' | '/' | '|'))
            });
        for (index, key) in keys.iter().enumerate() {
            if missing_gt
                && *key != "GT"
                && let Some(value) = cells.get_mut(index)
            {
                *value = if *key == "AD" {
                    ".,.".to_string()
                } else {
                    ".".to_string()
                };
            } else if !matches!(*key, "GT" | "AD" | "ADO" | "DP" | "PL")
                && let Some(value) = cells.get_mut(index)
            {
                *value = if bcf_output && string_fields.contains(*key) {
                    String::new()
                } else {
                    ".".to_string()
                };
            } else if homref
                && *key == "AD"
                && let Some(value) = cells.get_mut(index)
            {
                let width = value.split(',').count().max(1);
                *value = vec!["."; width].join(",");
            } else if homref
                && *key == "ADO"
                && let Some(value) = cells.get_mut(index)
            {
                *value = ".".to_string();
            }
        }
        *sample = cells.join(":");
    }
}

/// Header records installed by the legacy C++ `VariantWriter` constructor.
/// Definitions here win over conflicting input definitions during merging.
pub(super) const LEGACY_BASE_HEADERS: &[&str] = &[
    "##fileformat=VCFv4.1",
    "##FILTER=<ID=PASS,Description=\"All filters passed\">",
    "##reference=hg19",
    "##contig=<ID=chr1,length=249250621>",
    "##contig=<ID=chr2,length=243199373>",
    "##contig=<ID=chr3,length=198022430>",
    "##contig=<ID=chr4,length=191154276>",
    "##contig=<ID=chr5,length=180915260>",
    "##contig=<ID=chr6,length=171115067>",
    "##contig=<ID=chr7,length=159138663>",
    "##contig=<ID=chr8,length=146364022>",
    "##contig=<ID=chr9,length=141213431>",
    "##contig=<ID=chr10,length=135534747>",
    "##contig=<ID=chr11,length=135006516>",
    "##contig=<ID=chr12,length=133851895>",
    "##contig=<ID=chr13,length=115169878>",
    "##contig=<ID=chr14,length=107349540>",
    "##contig=<ID=chr15,length=102531392>",
    "##contig=<ID=chr16,length=90354753>",
    "##contig=<ID=chr17,length=81195210>",
    "##contig=<ID=chr18,length=78077248>",
    "##contig=<ID=chr19,length=59128983>",
    "##contig=<ID=chr20,length=63025520>",
    "##contig=<ID=chr21,length=48129895>",
    "##contig=<ID=chr22,length=51304566>",
    "##contig=<ID=chrX,length=155270560>",
    "##INFO=<ID=END,Number=.,Type=Integer,Description=\"SV end position\">",
    "##INFO=<ID=IMPORT_FAIL,Number=.,Type=Flag,Description=\"Flag to identify variants that could not be imported.\">",
    "##FORMAT=<ID=AGT,Number=1,Type=String,Description=\"Genotypes at ambiguous locations\">",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">",
    "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allele Depths\">",
    "##FORMAT=<ID=ADO,Number=.,Type=Integer,Description=\"Summed depth of non-called alleles.\">",
];

/// Reproduce `VariantWriterImpl::writeHeader`: keep constructor-installed
/// records first, sort supplied meta-lines lexicographically, and retain only
/// the first structured declaration for each `(header key, ID)` pair.
pub(crate) fn canonicalize_legacy_headers(headers: &[String]) -> Vec<String> {
    let mut output: Vec<String> = LEGACY_BASE_HEADERS
        .iter()
        .map(|line| (*line).to_string())
        .collect();
    let mut identities: BTreeSet<String> = output
        .iter()
        .filter_map(|line| structured_header_identity(line))
        .collect();
    let mut exact: BTreeSet<String> = output.iter().cloned().collect();

    let mut supplied: Vec<&String> = headers
        .iter()
        .filter(|line| line.starts_with("##") && !line.starts_with("##fileformat="))
        .collect();
    supplied.sort_unstable();
    for line in supplied {
        if let Some(identity) = structured_header_identity(line) {
            if identities.insert(identity) {
                exact.insert(line.clone());
                output.push(line.clone());
            }
        } else if exact.insert(line.clone()) {
            output.push(line.clone());
        }
    }

    if let Some(chrom) = headers.iter().rev().find(|line| line.starts_with("#CHROM")) {
        output.push(chrom.clone());
    }
    output
}

/// Return htslib's merge identity for `<ID=...>` header records.
pub(crate) fn structured_header_identity(line: &str) -> Option<String> {
    let body = line.strip_prefix("##")?;
    let (key, value) = body.split_once('=')?;
    let id_value = value.strip_prefix("<ID=")?;
    let id = id_value.split([',', '>']).next()?;
    Some(format!("{key}:{id}"))
}

/// Insert `ADO` into a record's FORMAT immediately after `AD`. For each sample
/// compute ADO = sum of AD values for allele indices that are NOT present in
/// the sample's GT call (i.e. the depth of "other" alleles that were not the
/// genotype's choice). Matches `VariantWriter.cpp:603` semantics.
pub(super) fn insert_ado_format(record: &mut RawVcfRecord) {
    let Some(format) = record.format.clone() else {
        return;
    };
    let fields: Vec<&str> = format.split(':').collect();
    if fields.contains(&"ADO") {
        return;
    }
    let Some(ad_index) = fields.iter().position(|field| *field == "AD") else {
        return;
    };
    let gt_index = fields.iter().position(|field| *field == "GT");

    let mut new_format_fields: Vec<String> = fields.iter().map(|s| s.to_string()).collect();
    new_format_fields.insert(ad_index + 1, "ADO".to_string());
    record.format = Some(new_format_fields.join(":"));

    for sample in &mut record.samples {
        let mut sample_fields: Vec<String> = sample.split(':').map(|s| s.to_string()).collect();
        // Pad with missing cells if the sample is shorter than FORMAT (legal
        // VCF trailing-missing shorthand); we only need up to `ad_index + 1`.
        while sample_fields.len() <= ad_index {
            sample_fields.push(".".to_string());
        }
        let ado = compute_ado(
            sample_fields
                .get(ad_index)
                .map(String::as_str)
                .unwrap_or("."),
            gt_index.and_then(|i| sample_fields.get(i).map(String::as_str)),
        );
        sample_fields.insert(ad_index + 1, ado.to_string());
        *sample = sample_fields.join(":");
    }
}

/// Compute `ADO` = sum of allele depths for allele indices that don't appear
/// in the genotype call. Returns 0 when either AD or GT is missing, matching
/// legacy's `bcf_int32_missing` → 0 rendering for unknown.
pub(super) fn compute_ado(ad_cell: &str, gt_cell: Option<&str>) -> i64 {
    if ad_cell == "." || ad_cell.is_empty() {
        return 0;
    }
    let ad_values: Vec<i64> = ad_cell
        .split(',')
        .map(|v| v.parse::<i64>().unwrap_or(0))
        .collect();
    if ad_values.is_empty() {
        return 0;
    }
    let gt = match gt_cell {
        Some(gt) if !gt.is_empty() && gt != "." => gt,
        _ => return 0,
    };
    let called: std::collections::BTreeSet<usize> = gt
        .split(['/', '|'])
        .filter_map(|a| a.parse::<usize>().ok())
        .collect();
    if called.is_empty() {
        return 0;
    }
    ad_values
        .iter()
        .enumerate()
        .filter_map(|(index, depth)| (!called.contains(&index)).then_some(*depth))
        .sum()
}
