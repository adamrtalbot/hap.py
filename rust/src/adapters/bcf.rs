use crate::domain::RawVcfRecord;
use crate::output::{FailureOperation, OutputTransaction, fail_operation};
use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use std::collections::BTreeMap;
use std::fs;
use std::fs::File;
use std::io::{ErrorKind, Read, Write};
use std::path::{Path, PathBuf};

const BCF_MAGIC: &[u8; 5] = b"BCF\x02\x02";
/// Hard limits are deliberately below `u32::MAX`: corrupt length fields must
/// not be able to turn a tiny input into a multi-gigabyte allocation.
pub(crate) const MAX_BCF_HEADER_BYTES: usize = 64 * 1024 * 1024;
pub(crate) const MAX_BCF_RECORD_BYTES: usize = 64 * 1024 * 1024;
const MAX_BCF_SAMPLES: usize = 1_000_000;
const MAX_BCF_FIELDS: usize = 65_535;
const INT8_MISSING: i8 = i8::MIN;
const INT8_END: i8 = i8::MIN + 1;
const INT16_MISSING: i16 = i16::MIN;
const INT16_END: i16 = i16::MIN + 1;
const INT32_MISSING: i32 = i32::MIN;
const INT32_END: i32 = i32::MIN + 1;
const FLOAT_MISSING: u32 = 0x7f80_0001;
const FLOAT_END: u32 = 0x7f80_0002;

#[derive(Clone, Debug)]
struct HeaderDictionary {
    headers: Vec<String>,
    contigs: Vec<String>,
    keys: Vec<String>,
    info_types: BTreeMap<String, String>,
    format_types: BTreeMap<String, String>,
    sample_count: usize,
}

pub(crate) fn is_bcf_data(data: &[u8]) -> bool {
    data.starts_with(BCF_MAGIC) || data.starts_with(b"BCF\x02\x01")
}

#[cfg(test)]
pub(crate) fn decode(data: &[u8], path: &Path) -> Result<(Vec<String>, Vec<RawVcfRecord>)> {
    if !is_bcf_data(data) {
        bail!("{} is not BCF2", path.display());
    }
    let mut cursor = Cursor::new(data);
    cursor.take(5)?;
    let header_len = checked_header_len(cursor.u32()?, path)?;
    let header_bytes = cursor.take(header_len)?;
    let header_text = std::str::from_utf8(header_bytes)
        .with_context(|| format!("BCF header in {} is not UTF-8", path.display()))?
        .trim_end_matches('\0');
    let dictionary = parse_bcf_header(header_text)?;
    let mut records = Vec::new();
    while cursor.remaining() > 0 {
        if cursor.remaining() < 8 {
            bail!("truncated BCF record header in {}", path.display());
        }
        let record_number = records.len() + 1;
        let shared_len = cursor.u32()? as usize;
        let individual_len = cursor.u32()? as usize;
        checked_record_lengths(shared_len, individual_len, path, record_number)?;
        let shared = cursor.take(shared_len)?;
        let individual = cursor.take(individual_len)?;
        records.push(decode_record(shared, individual, &dictionary, path)?);
    }
    Ok((dictionary.headers, records))
}

pub(crate) struct RecordReader {
    reader: Box<dyn Read>,
    dictionary: HeaderDictionary,
    path: PathBuf,
    record_number: usize,
    payload: Vec<u8>,
}

impl RecordReader {
    pub(crate) fn new(mut reader: Box<dyn Read>, path: &Path) -> Result<Self> {
        let mut magic = [0; 5];
        reader
            .read_exact(&mut magic)
            .with_context(|| format!("truncated BCF magic in {}", path.display()))?;
        if magic != *BCF_MAGIC && magic != *b"BCF\x02\x01" {
            bail!("{} is not BCF2", path.display());
        }
        let header_len =
            checked_header_len(read_u32(&mut reader, "BCF header length", path)?, path)?;
        let mut header = vec![0; header_len];
        reader
            .read_exact(&mut header)
            .with_context(|| format!("truncated BCF header in {}", path.display()))?;
        let text = std::str::from_utf8(&header)
            .with_context(|| format!("BCF header in {} is not UTF-8", path.display()))?
            .trim_end_matches('\0');
        Ok(Self {
            reader,
            dictionary: parse_bcf_header(text)?,
            path: path.to_path_buf(),
            record_number: 0,
            payload: Vec::new(),
        })
    }

    pub(crate) fn headers(&self) -> &[String] {
        &self.dictionary.headers
    }

    pub(crate) fn next_record(&mut self) -> Result<Option<RawVcfRecord>> {
        let mut lengths = [0; 8];
        let mut read = 0;
        while read < lengths.len() {
            match self.reader.read(&mut lengths[read..]) {
                Ok(0) if read == 0 => return Ok(None),
                Ok(0) => bail!(
                    "truncated BCF record header at record {} in {}",
                    self.record_number + 1,
                    self.path.display()
                ),
                Ok(count) => read += count,
                Err(error) if error.kind() == ErrorKind::Interrupted => {}
                Err(error) => {
                    return Err(error)
                        .with_context(|| format!("failed to read {}", self.path.display()));
                }
            }
        }
        self.record_number += 1;
        let shared_len = u32::from_le_bytes(lengths[..4].try_into().unwrap()) as usize;
        let individual_len = u32::from_le_bytes(lengths[4..].try_into().unwrap()) as usize;
        checked_record_lengths(shared_len, individual_len, &self.path, self.record_number)?;
        self.payload.resize(shared_len + individual_len, 0);
        let (shared, individual) = self.payload.split_at_mut(shared_len);
        self.reader.read_exact(shared).with_context(|| {
            format!(
                "truncated shared BCF data at record {} in {}",
                self.record_number,
                self.path.display()
            )
        })?;
        self.reader.read_exact(individual).with_context(|| {
            format!(
                "truncated individual BCF data at record {} in {}",
                self.record_number,
                self.path.display()
            )
        })?;
        decode_record(shared, individual, &self.dictionary, &self.path)
            .with_context(|| {
                format!(
                    "failed to decode BCF record {} in {}",
                    self.record_number,
                    self.path.display()
                )
            })
            .map(Some)
    }
}

fn read_u32(reader: &mut dyn Read, label: &str, path: &Path) -> Result<u32> {
    let mut bytes = [0; 4];
    reader
        .read_exact(&mut bytes)
        .with_context(|| format!("truncated {label} in {}", path.display()))?;
    Ok(u32::from_le_bytes(bytes))
}

fn checked_header_len(raw: u32, path: &Path) -> Result<usize> {
    let length = raw as usize;
    if length == 0 || length > MAX_BCF_HEADER_BYTES {
        bail!(
            "BCF header length {length} exceeds limit {MAX_BCF_HEADER_BYTES} in {}",
            path.display()
        );
    }
    Ok(length)
}

fn checked_record_lengths(
    shared: usize,
    individual: usize,
    path: &Path,
    record: usize,
) -> Result<()> {
    let total = shared
        .checked_add(individual)
        .context("BCF record length overflow")?;
    if shared < 24 || total > MAX_BCF_RECORD_BYTES {
        bail!(
            "BCF record {record} in {} has invalid lengths shared={shared}, individual={individual}; maximum combined length is {MAX_BCF_RECORD_BYTES}",
            path.display()
        );
    }
    Ok(())
}

/// Retrieve every BCF record reachable through a companion CSI index.
///
/// The parity verifier uses this instead of an external `bcftools view -r`
/// process. CSI chunks contain BGZF virtual offsets, so following them proves
/// that the published index can actually seek into and decode the data file;
/// merely parsing the index structure would not detect omitted records or
/// invalid chunk boundaries.
#[cfg(test)]
pub(crate) fn read_indexed_records(
    path: &Path,
    index_path: &Path,
) -> Result<BTreeMap<String, Vec<RawVcfRecord>>> {
    let uncompressed = read_uncompressed(path)
        .with_context(|| format!("failed to read indexed BCF {}", path.display()))?;
    let dictionary = decode_header_dictionary(&uncompressed, path)?;
    let indexed_chunks = read_csi_chunks(index_path)?;
    if indexed_chunks.len() > dictionary.contigs.len() {
        bail!(
            "CSI contains {} references but BCF header contains {}",
            indexed_chunks.len(),
            dictionary.contigs.len()
        );
    }

    let mut indexed = dictionary
        .contigs
        .iter()
        .cloned()
        .map(|contig| (contig, Vec::new()))
        .collect::<BTreeMap<_, _>>();
    let mut reader = bgzf::io::Reader::new(File::open(path)?);
    for (rid, chunks) in indexed_chunks.into_iter().enumerate() {
        let contig = dictionary
            .contigs
            .get(rid)
            .context("CSI reference is absent from BCF header")?;
        for (start, end) in merge_chunks(chunks) {
            if start == end {
                continue;
            }
            reader
                .seek(bgzf::VirtualPosition::from(start))
                .with_context(|| format!("failed to seek CSI chunk for {contig}"))?;
            while u64::from(reader.virtual_position()) < end {
                let before = u64::from(reader.virtual_position());
                let record = read_record(&mut reader, &dictionary, path)?;
                let after = u64::from(reader.virtual_position());
                if after <= before {
                    bail!("BCF reader made no progress in CSI chunk for {contig}");
                }
                if record.chrom == *contig {
                    indexed
                        .get_mut(contig)
                        .expect("BCF contig map was initialized")
                        .push(record);
                }
            }
        }
    }
    Ok(indexed)
}

#[cfg(test)]
fn decode_header_dictionary(data: &[u8], path: &Path) -> Result<HeaderDictionary> {
    if !is_bcf_data(data) {
        bail!("{} is not BCF2", path.display());
    }
    let mut cursor = Cursor::new(data);
    cursor.take(5)?;
    let header_len = cursor.u32()? as usize;
    let header_bytes = cursor.take(header_len)?;
    let header_text = std::str::from_utf8(header_bytes)
        .with_context(|| format!("BCF header in {} is not UTF-8", path.display()))?
        .trim_end_matches('\0');
    parse_bcf_header(header_text)
}

#[cfg(test)]
fn read_csi_chunks(path: &Path) -> Result<Vec<Vec<(u64, u64)>>> {
    Ok(read_csi_bins(path)?
        .references
        .into_iter()
        .map(|bins| {
            bins.into_iter()
                .flat_map(|bin| bin.chunks)
                .collect::<Vec<_>>()
        })
        .collect())
}

#[cfg(test)]
struct CsiIndex {
    references: Vec<Vec<CsiBin>>,
}

#[cfg(test)]
struct CsiBin {
    chunks: Vec<(u64, u64)>,
}

#[cfg(test)]
fn read_csi_bins(path: &Path) -> Result<CsiIndex> {
    let payload = read_uncompressed(path)
        .with_context(|| format!("failed to read CSI index {}", path.display()))?;
    let mut cursor = Cursor::new(&payload);
    if cursor.take(4)? != b"CSI\x01" {
        bail!("{} is not a CSI index", path.display());
    }
    let min_shift = nonnegative_i32(&mut cursor, "minimum shift")?;
    let depth = nonnegative_i32(&mut cursor, "depth")?;
    let auxiliary_len = nonnegative_i32(&mut cursor, "auxiliary length")?;
    cursor.take(auxiliary_len)?;
    let reference_count = nonnegative_i32(&mut cursor, "reference count")?;
    let level_bits = depth
        .checked_add(1)
        .and_then(|levels| levels.checked_mul(3))
        .context("CSI depth overflow")?;
    if level_bits >= u64::BITS as usize || min_shift >= usize::BITS as usize {
        bail!("CSI indexing parameters are out of range");
    }
    let metadata_bin = (((1u64 << level_bits) - 1) / 7)
        .checked_add(1)
        .context("CSI metadata bin overflow")?;

    let mut references = Vec::with_capacity(reference_count);
    for _ in 0..reference_count {
        let bin_count = nonnegative_i32(&mut cursor, "bin count")?;
        let mut bins = Vec::new();
        for _ in 0..bin_count {
            let bin = u64::from(cursor.u32()?);
            let _loffset = cursor.u64()?;
            let chunk_count = nonnegative_i32(&mut cursor, "chunk count")?;
            if bin > metadata_bin {
                bail!("CSI bin {bin} exceeds metadata bin {metadata_bin}");
            }
            let mut chunks = Vec::new();
            for _ in 0..chunk_count {
                let start = cursor.u64()?;
                let end = cursor.u64()?;
                if bin != metadata_bin {
                    if start > end {
                        bail!("CSI chunk start exceeds its end");
                    }
                    chunks.push((start, end));
                }
            }
            if bin != metadata_bin {
                bins.push(CsiBin { chunks });
            }
        }
        references.push(bins);
    }
    match cursor.remaining() {
        0 => {}
        8 => {
            let _unplaced_records = cursor.u64()?;
        }
        trailing => bail!("CSI index has {trailing} trailing bytes"),
    }
    Ok(CsiIndex { references })
}

#[cfg(test)]
fn nonnegative_i32(cursor: &mut Cursor<'_>, field: &str) -> Result<usize> {
    usize::try_from(cursor.i32()?).with_context(|| format!("CSI {field} is negative"))
}

#[cfg(test)]
fn merge_chunks(mut chunks: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    chunks.sort_unstable_by_key(|chunk| chunk.0);
    let mut merged: Vec<(u64, u64)> = Vec::new();
    for (start, end) in chunks {
        if let Some(previous) = merged.last_mut()
            && start <= previous.1
        {
            previous.1 = previous.1.max(end);
        } else {
            merged.push((start, end));
        }
    }
    merged
}

#[cfg(test)]
fn read_record<R: Read>(
    reader: &mut R,
    dictionary: &HeaderDictionary,
    path: &Path,
) -> Result<RawVcfRecord> {
    let mut lengths = [0u8; 8];
    reader
        .read_exact(&mut lengths)
        .with_context(|| format!("truncated indexed BCF record in {}", path.display()))?;
    let shared_len = u32::from_le_bytes(lengths[..4].try_into().unwrap()) as usize;
    let individual_len = u32::from_le_bytes(lengths[4..].try_into().unwrap()) as usize;
    let total = shared_len
        .checked_add(individual_len)
        .context("indexed BCF record length overflow")?;
    if shared_len < 24 || total > MAX_BCF_RECORD_BYTES {
        bail!(
            "indexed BCF record in {} has invalid lengths shared={shared_len}, individual={individual_len}",
            path.display()
        );
    }
    let mut shared = vec![0; shared_len];
    let mut individual = vec![0; individual_len];
    reader.read_exact(&mut shared)?;
    reader.read_exact(&mut individual)?;
    decode_record(&shared, &individual, dictionary, path)
}

fn parse_bcf_header(text: &str) -> Result<HeaderDictionary> {
    let mut headers = Vec::new();
    let mut contig_slots = BTreeMap::new();
    let mut key_slots = BTreeMap::from([(0usize, "PASS".to_string())]);
    let mut info_types = BTreeMap::new();
    let mut format_types = BTreeMap::new();
    let mut sample_count = 0;

    for line in text.lines() {
        if line.starts_with("#CHROM\t") {
            sample_count = line.split('\t').count().saturating_sub(9);
            headers.push(strip_idx_attribute(line));
            continue;
        } else if let Some(body) = line.strip_prefix("##contig=<") {
            let id = header_attribute(body, "ID").unwrap_or_default();
            if let Some(index) = header_attribute(body, "IDX").and_then(|v| v.parse().ok()) {
                contig_slots.insert(index, id.to_string());
            }
        } else if let Some((kind, body)) =
            ["FILTER", "INFO", "FORMAT"].into_iter().find_map(|kind| {
                line.strip_prefix(&format!("##{kind}=<"))
                    .map(|body| (kind, body))
            })
        {
            let id = header_attribute(body, "ID").unwrap_or_default();
            if let Some(index) = header_attribute(body, "IDX").and_then(|v| v.parse().ok()) {
                key_slots.entry(index).or_insert_with(|| id.to_string());
            }
            if kind == "INFO" {
                info_types.insert(
                    id.to_string(),
                    header_attribute(body, "Type")
                        .unwrap_or("String")
                        .to_string(),
                );
            } else if kind == "FORMAT" {
                format_types.insert(
                    id.to_string(),
                    header_attribute(body, "Type")
                        .unwrap_or("String")
                        .to_string(),
                );
            }
        }
        headers.push(strip_idx_attribute(line));
    }

    let contigs = dense_dictionary(contig_slots, "contig")?;
    let keys = dense_dictionary(key_slots, "header key")?;
    Ok(HeaderDictionary {
        headers,
        contigs,
        keys,
        info_types,
        format_types,
        sample_count,
    })
}

fn dense_dictionary(slots: BTreeMap<usize, String>, label: &str) -> Result<Vec<String>> {
    let Some(maximum) = slots.keys().next_back().copied() else {
        return Ok(Vec::new());
    };
    let mut values = vec![String::new(); maximum + 1];
    for (index, value) in slots {
        values[index] = value;
    }
    if let Some(index) = values.iter().position(String::is_empty) {
        bail!("BCF {label} dictionary is missing index {index}");
    }
    Ok(values)
}

fn header_attribute<'a>(body: &'a str, name: &str) -> Option<&'a str> {
    let needle = format!("{name}=");
    let start = body.find(&needle)? + needle.len();
    let rest = &body[start..];
    if let Some(quoted) = rest.strip_prefix('"') {
        let end = quoted.find('"')?;
        Some(&quoted[..end])
    } else {
        let end = rest.find([',', '>']).unwrap_or(rest.len());
        Some(&rest[..end])
    }
}

fn strip_idx_attribute(line: &str) -> String {
    let Some(start) = line.find(",IDX=") else {
        return line.to_string();
    };
    let tail = &line[start + 5..];
    let end = tail.find([',', '>']).unwrap_or(tail.len());
    let mut cleaned = String::with_capacity(line.len());
    cleaned.push_str(&line[..start]);
    cleaned.push_str(&tail[end..]);
    cleaned
}

fn decode_record(
    shared: &[u8],
    individual: &[u8],
    dictionary: &HeaderDictionary,
    path: &Path,
) -> Result<RawVcfRecord> {
    let mut shared = Cursor::new(shared);
    let rid = shared.i32()?;
    let pos = shared.i32()?;
    let rlen = shared.i32()?;
    let qual_bits = shared.u32()?;
    let allele_info = shared.u32()?;
    let fmt_sample = shared.u32()?;
    if rid < 0 || pos < 0 || rlen < 0 {
        bail!("negative BCF RID/POS/RLEN in {}", path.display());
    }
    pos.checked_add(rlen.max(1))
        .context("BCF coordinate overflow")?;
    let chrom = dictionary
        .contigs
        .get(rid as usize)
        .with_context(|| format!("BCF RID {rid} is absent from header in {}", path.display()))?
        .clone();
    let n_info = (allele_info & 0xffff) as usize;
    let n_allele = (allele_info >> 16) as usize;
    let n_sample = (fmt_sample & 0x00ff_ffff) as usize;
    let n_fmt = (fmt_sample >> 24) as usize;
    if n_allele == 0
        || n_info > MAX_BCF_FIELDS
        || n_fmt > MAX_BCF_FIELDS
        || n_sample > MAX_BCF_SAMPLES
    {
        bail!(
            "BCF record counts exceed limits in {}: alleles={n_allele}, info={n_info}, format={n_fmt}, samples={n_sample}",
            path.display()
        );
    }
    if n_sample != dictionary.sample_count {
        bail!(
            "BCF sample count {n_sample} does not match header count {} in {}",
            dictionary.sample_count,
            path.display()
        );
    }
    if n_sample > individual.len() && n_fmt > 0 {
        bail!(
            "BCF sample count {n_sample} is inconsistent with {} individual bytes in {}",
            individual.len(),
            path.display()
        );
    }

    let id = decode_typed_string(&mut shared)?;
    let mut alleles = Vec::with_capacity(n_allele);
    for _ in 0..n_allele {
        alleles.push(decode_typed_string(&mut shared)?);
    }
    if alleles.is_empty() {
        bail!("BCF record has no REF allele in {}", path.display());
    }
    let filter_indexes = decode_typed_ints(&mut shared)?;
    let filter = if filter_indexes.is_empty() {
        ".".to_string()
    } else {
        filter_indexes
            .into_iter()
            .filter_map(|value| value.and_then(|value| dictionary.keys.get(value as usize)))
            .cloned()
            .collect::<Vec<_>>()
            .join(";")
    };

    let mut info = Vec::with_capacity(n_info);
    for _ in 0..n_info {
        let key = decode_typed_ints(&mut shared)?
            .first()
            .and_then(|value| *value)
            .and_then(|value| dictionary.keys.get(value as usize))
            .context("BCF INFO key is missing from header dictionary")?
            .clone();
        let declared = dictionary.info_types.get(&key).map(String::as_str);
        let value = decode_typed_value(&mut shared, declared, false)?;
        if declared == Some("Flag") {
            info.push(key);
        } else {
            info.push(format!("{key}={value}"));
        }
    }

    let mut indiv = Cursor::new(individual);
    let mut format = Vec::with_capacity(n_fmt);
    let mut samples = vec![Vec::with_capacity(n_fmt); n_sample];
    for _ in 0..n_fmt {
        let key = decode_typed_ints(&mut indiv)?
            .first()
            .and_then(|value| *value)
            .and_then(|value| dictionary.keys.get(value as usize))
            .context("BCF FORMAT key is missing from header dictionary")?
            .clone();
        let (width, atomic_type) = indiv.typed_header()?;
        format.push(key.clone());
        for sample in &mut samples {
            sample.push(decode_format_cell(
                &mut indiv,
                width,
                atomic_type,
                key == "GT",
            )?);
        }
    }

    Ok(RawVcfRecord {
        chrom,
        pos: pos as usize + 1,
        id: if id.is_empty() { ".".to_string() } else { id },
        ref_allele: alleles.remove(0),
        alt_allele: if alleles.is_empty() {
            ".".to_string()
        } else {
            alleles.join(",")
        },
        qual: if qual_bits == FLOAT_MISSING {
            ".".to_string()
        } else {
            format_float(f32::from_bits(qual_bits))
        },
        filter: if filter.is_empty() {
            ".".to_string()
        } else {
            filter
        },
        info: if info.is_empty() {
            ".".to_string()
        } else {
            info.join(";")
        },
        format: (!format.is_empty()).then(|| format.join(":")),
        samples: samples.into_iter().map(|cells| cells.join(":")).collect(),
    })
}

fn decode_typed_string(cursor: &mut Cursor<'_>) -> Result<String> {
    let (length, atomic_type) = cursor.typed_header()?;
    if atomic_type != 7 && length != 0 {
        bail!("BCF string uses unexpected atomic type {atomic_type}");
    }
    let bytes = cursor.take(length)?;
    Ok(String::from_utf8_lossy(bytes)
        .trim_end_matches('\0')
        .to_string())
}

fn decode_typed_ints(cursor: &mut Cursor<'_>) -> Result<Vec<Option<i32>>> {
    let (length, atomic_type) = cursor.typed_header()?;
    (0..length).map(|_| cursor.integer(atomic_type)).collect()
}

fn decode_typed_value(
    cursor: &mut Cursor<'_>,
    declared: Option<&str>,
    genotype: bool,
) -> Result<String> {
    let (length, atomic_type) = cursor.typed_header()?;
    if length == 0 {
        return Ok(String::new());
    }
    decode_format_cell(
        cursor,
        length,
        atomic_type,
        genotype || declared == Some("GT"),
    )
}

fn decode_format_cell(
    cursor: &mut Cursor<'_>,
    width: usize,
    atomic_type: u8,
    genotype: bool,
) -> Result<String> {
    if atomic_type == 7 {
        let bytes = cursor.take(width)?;
        return Ok(String::from_utf8_lossy(bytes)
            .trim_end_matches('\0')
            .to_string());
    }
    if atomic_type == 5 {
        let mut values = Vec::new();
        for _ in 0..width {
            let bits = cursor.u32()?;
            if bits == FLOAT_END {
                continue;
            }
            values.push(if bits == FLOAT_MISSING {
                ".".to_string()
            } else {
                format_float(f32::from_bits(bits))
            });
        }
        return Ok(if values.is_empty() {
            ".".to_string()
        } else {
            values.join(",")
        });
    }

    let mut values = Vec::new();
    for _ in 0..width {
        match cursor.integer_with_end(atomic_type)? {
            IntegerValue::End => {}
            IntegerValue::Missing => values.push(None),
            IntegerValue::Value(value) => values.push(Some(value)),
        }
    }
    if genotype {
        if values.is_empty() {
            return Ok(".".to_string());
        }
        let mut rendered = String::new();
        for (index, encoded) in values.into_iter().enumerate() {
            if index > 0 {
                rendered.push(if encoded.is_some_and(|value| value & 1 != 0) {
                    '|'
                } else {
                    '/'
                });
            }
            match encoded {
                None | Some(0) => rendered.push('.'),
                Some(value) => rendered.push_str(&((value >> 1) - 1).to_string()),
            }
        }
        Ok(rendered)
    } else {
        Ok(if values.is_empty() {
            ".".to_string()
        } else {
            values
                .into_iter()
                .map(|value| value.map_or_else(|| ".".to_string(), |value| value.to_string()))
                .collect::<Vec<_>>()
                .join(",")
        })
    }
}

fn format_float(value: f32) -> String {
    if value == 0.0 {
        "0".to_string()
    } else {
        value.to_string()
    }
}

#[derive(Copy, Clone)]
enum IntegerValue {
    Missing,
    End,
    Value(i32),
}

struct Cursor<'a> {
    data: &'a [u8],
    offset: usize,
}

impl<'a> Cursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, offset: 0 }
    }

    #[cfg(test)]
    fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.offset)
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8]> {
        let end = self
            .offset
            .checked_add(length)
            .context("BCF length overflow")?;
        if end > self.data.len() {
            bail!("truncated BCF value");
        }
        let value = &self.data[self.offset..end];
        self.offset = end;
        Ok(value)
    }

    fn u8(&mut self) -> Result<u8> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take(2)?.try_into().unwrap()))
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    #[cfg(test)]
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn typed_header(&mut self) -> Result<(usize, u8)> {
        let descriptor = self.u8()?;
        let atomic_type = descriptor & 0x0f;
        let mut length = usize::from(descriptor >> 4);
        if length == 15 {
            length = decode_length(self)?;
        }
        Ok((length, atomic_type))
    }

    fn integer(&mut self, atomic_type: u8) -> Result<Option<i32>> {
        Ok(match self.integer_with_end(atomic_type)? {
            IntegerValue::Value(value) => Some(value),
            IntegerValue::Missing | IntegerValue::End => None,
        })
    }

    fn integer_with_end(&mut self, atomic_type: u8) -> Result<IntegerValue> {
        Ok(match atomic_type {
            1 => match self.u8()? as i8 {
                INT8_MISSING => IntegerValue::Missing,
                INT8_END => IntegerValue::End,
                value => IntegerValue::Value(i32::from(value)),
            },
            2 => match self.u16()? as i16 {
                INT16_MISSING => IntegerValue::Missing,
                INT16_END => IntegerValue::End,
                value => IntegerValue::Value(i32::from(value)),
            },
            3 => match self.i32()? {
                INT32_MISSING => IntegerValue::Missing,
                INT32_END => IntegerValue::End,
                value => IntegerValue::Value(value),
            },
            other => bail!("unsupported BCF integer atomic type {other}"),
        })
    }
}

fn decode_length(cursor: &mut Cursor<'_>) -> Result<usize> {
    let (length, atomic_type) = cursor.typed_header()?;
    if length != 1 {
        bail!("BCF extended vector length is not scalar");
    }
    cursor
        .integer(atomic_type)?
        .and_then(|value| usize::try_from(value).ok())
        .context("invalid BCF extended vector length")
}

#[cfg(test)]
pub(crate) fn write(path: &Path, headers: &[String], records: &[RawVcfRecord]) -> Result<()> {
    write_iter(path, headers, records.iter().cloned().map(Ok))
}

pub(crate) fn write_iter<I>(path: &Path, headers: &[String], records: I) -> Result<()>
where
    I: IntoIterator<Item = Result<RawVcfRecord>>,
{
    let csi_path = path.with_extension("bcf.csi");
    let transaction = OutputTransaction::files(Vec::<PathBuf>::new(), [path, &csi_path])?;
    let staged_bcf = transaction.staged_file(path)?.to_path_buf();
    let staged_csi = transaction.staged_file(&csi_path)?.to_path_buf();
    write_inner(&staged_bcf, &staged_csi, path, &csi_path, headers, records).map_err(|error| {
        anyhow::anyhow!(
            "failed to write BCF destination {}: {error:#}",
            path.display()
        )
    })?;
    transaction.commit()
}

fn write_inner<I>(
    path: &Path,
    csi_path: &Path,
    logical_bcf: &Path,
    logical_csi: &Path,
    headers: &[String],
    records: I,
) -> Result<()>
where
    I: IntoIterator<Item = Result<RawVcfRecord>>,
{
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    let encoded_header = encode_header(headers)?;
    (|| {
        fail_operation(FailureOperation::Writer, logical_bcf)?;
        let mut writer = bgzf::io::Writer::new(
            File::create(path)
                .with_context(|| format!("failed to create {}", logical_bcf.display()))?,
        );
        writer.write_all(BCF_MAGIC)?;
        writer.write_all(&(encoded_header.text.len() as u32).to_le_bytes())?;
        writer.write_all(&encoded_header.text)?;

        let mut chunks = vec![None::<(u64, u64)>; encoded_header.dictionary.contigs.len()];
        let mut record_counts = vec![0_u64; encoded_header.dictionary.contigs.len()];
        for record in records {
            let record = record?;
            let rid = encoded_header
                .contig_indexes
                .get(&record.chrom)
                .copied()
                .with_context(|| format!("BCF output contig {} has no header", record.chrom))?;
            let start = u64::from(writer.virtual_position());
            let (shared, individual) = encode_record(&record, rid, &encoded_header)?;
            writer.write_all(&(shared.len() as u32).to_le_bytes())?;
            writer.write_all(&(individual.len() as u32).to_le_bytes())?;
            writer.write_all(&shared)?;
            writer.write_all(&individual)?;
            let end = u64::from(writer.virtual_position());
            chunks[rid as usize].get_or_insert((start, end)).1 = end;
            record_counts[rid as usize] += 1;
        }
        fail_operation(FailureOperation::Encoder, logical_bcf)?;
        writer
            .finish()
            .with_context(|| format!("failed to finish BCF encoder for {}", logical_bcf.display()))?
            .sync_all()
            .with_context(|| format!("failed to sync BCF artifact {}", logical_bcf.display()))?;
        fail_operation(FailureOperation::Index, logical_csi)?;
        write_csi(csi_path, &chunks, &record_counts)
            .with_context(|| format!("failed to write CSI artifact {}", logical_csi.display()))?;
        Ok::<(), anyhow::Error>(())
    })()
}

struct EncodedHeader {
    text: Vec<u8>,
    dictionary: HeaderDictionary,
    contig_indexes: BTreeMap<String, i32>,
    key_indexes: BTreeMap<String, i32>,
}

fn encode_header(headers: &[String]) -> Result<EncodedHeader> {
    let mut clean = headers.to_vec();
    if !clean
        .iter()
        .any(|line| line.starts_with("##FILTER=<ID=PASS,"))
    {
        clean.insert(
            1.min(clean.len()),
            "##FILTER=<ID=PASS,Description=\"All filters passed\">".into(),
        );
    }
    let mut contig_indexes = BTreeMap::new();
    let mut key_indexes = BTreeMap::from([("PASS".to_string(), 0i32)]);
    let mut next_key = 1i32;
    let mut encoded = Vec::with_capacity(clean.len());
    for line in &clean {
        let mut line = strip_idx_attribute(line);
        if let Some(body) = line.strip_prefix("##contig=<") {
            let id = header_attribute(body, "ID").context("BCF contig header lacks ID")?;
            let index = contig_indexes.len() as i32;
            contig_indexes.entry(id.to_string()).or_insert(index);
            line = append_idx(&line, index);
        } else if let Some(body) = ["FILTER", "INFO", "FORMAT"]
            .into_iter()
            .find_map(|kind| line.strip_prefix(&format!("##{kind}=<")))
        {
            let id = header_attribute(body, "ID").context("BCF dictionary header lacks ID")?;
            let index = if id == "PASS" {
                0
            } else if let Some(index) = key_indexes.get(id) {
                *index
            } else {
                let index = next_key;
                next_key += 1;
                key_indexes.insert(id.to_string(), index);
                index
            };
            line = append_idx(&line, index);
        }
        encoded.push(line);
    }
    let mut text = encoded.join("\n").into_bytes();
    text.push(b'\n');
    text.push(0);
    let dictionary = parse_bcf_header(std::str::from_utf8(&text)?.trim_end_matches('\0'))?;
    Ok(EncodedHeader {
        text,
        dictionary,
        contig_indexes,
        key_indexes,
    })
}

fn append_idx(line: &str, index: i32) -> String {
    line.strip_suffix('>')
        .map_or_else(|| line.to_string(), |body| format!("{body},IDX={index}>"))
}

fn encode_record(
    record: &RawVcfRecord,
    rid: i32,
    header: &EncodedHeader,
) -> Result<(Vec<u8>, Vec<u8>)> {
    let alts: Vec<&str> = record
        .alt_allele
        .split(',')
        .filter(|alt| *alt != "." && !alt.is_empty())
        .collect();
    let alleles = 1 + alts.len();
    let info_entries: Vec<&str> = if record.info == "." || record.info.is_empty() {
        Vec::new()
    } else {
        record.info.split(';').collect()
    };
    let format_keys: Vec<&str> = record
        .format
        .as_deref()
        .map(|format| format.split(':').collect())
        .unwrap_or_default();
    let mut shared = Vec::new();
    shared.extend_from_slice(&rid.to_le_bytes());
    shared.extend_from_slice(&(i32::try_from(record.pos)? - 1).to_le_bytes());
    shared.extend_from_slice(&(record.ref_allele.len() as i32).to_le_bytes());
    let qual = if record.qual == "." {
        FLOAT_MISSING
    } else {
        record.qual.parse::<f32>()?.to_bits()
    };
    shared.extend_from_slice(&qual.to_le_bytes());
    shared.extend_from_slice(&(((alleles as u32) << 16) | info_entries.len() as u32).to_le_bytes());
    shared.extend_from_slice(
        &(((format_keys.len() as u32) << 24) | record.samples.len() as u32).to_le_bytes(),
    );
    encode_string(&mut shared, if record.id == "." { "" } else { &record.id });
    encode_string(&mut shared, &record.ref_allele);
    for alt in alts {
        encode_string(&mut shared, alt);
    }
    let filters: Vec<i32> = if record.filter == "." || record.filter.is_empty() {
        Vec::new()
    } else {
        record
            .filter
            .split(';')
            .map(|name| {
                header
                    .key_indexes
                    .get(name)
                    .copied()
                    .with_context(|| format!("BCF FILTER {name} lacks header"))
            })
            .collect::<Result<_>>()?
    };
    encode_int_vector(
        &mut shared,
        &filters.iter().copied().map(Some).collect::<Vec<_>>(),
    );
    for entry in info_entries {
        let (key, value) = entry
            .split_once('=')
            .map_or((entry, None), |(key, value)| (key, Some(value)));
        let index = header
            .key_indexes
            .get(key)
            .copied()
            .with_context(|| format!("BCF INFO {key} lacks header"))?;
        encode_int_vector(&mut shared, &[Some(index)]);
        match header.dictionary.info_types.get(key).map(String::as_str) {
            Some("Flag") => shared.push(0),
            Some("Integer") => encode_integer_text(&mut shared, value.unwrap_or("."))?,
            Some("Float") => encode_float_text(&mut shared, value.unwrap_or("."))?,
            _ => encode_string(&mut shared, value.unwrap_or("")),
        }
    }

    let sample_cells: Vec<Vec<&str>> = record
        .samples
        .iter()
        .map(|sample| sample.split(':').collect())
        .collect();
    let mut individual = Vec::new();
    for (field_index, key) in format_keys.iter().enumerate() {
        let index = header
            .key_indexes
            .get(*key)
            .copied()
            .with_context(|| format!("BCF FORMAT {key} lacks header"))?;
        encode_int_vector(&mut individual, &[Some(index)]);
        let cells: Vec<&str> = sample_cells
            .iter()
            .map(|sample| sample.get(field_index).copied().unwrap_or("."))
            .collect();
        if *key == "GT" {
            encode_gt_cells(&mut individual, &cells)?;
        } else {
            match header.dictionary.format_types.get(*key).map(String::as_str) {
                Some("Integer") => encode_integer_cells(&mut individual, &cells)?,
                Some("Float") => encode_float_cells(&mut individual, &cells)?,
                _ => encode_string_cells(&mut individual, &cells),
            }
        }
    }
    Ok((shared, individual))
}

fn encode_type(output: &mut Vec<u8>, length: usize, atomic_type: u8) {
    if length < 15 {
        output.push(((length as u8) << 4) | atomic_type);
    } else {
        output.push(0xf0 | atomic_type);
        encode_int_vector(output, &[Some(length as i32)]);
    }
}

fn encode_string(output: &mut Vec<u8>, value: &str) {
    encode_type(output, value.len(), 7);
    output.extend_from_slice(value.as_bytes());
}

fn encode_int_vector(output: &mut Vec<u8>, values: &[Option<i32>]) {
    encode_type(output, values.len(), 3);
    for value in values {
        output.extend_from_slice(&value.unwrap_or(INT32_MISSING).to_le_bytes());
    }
}

fn encode_integer_text(output: &mut Vec<u8>, text: &str) -> Result<()> {
    let values = text
        .split(',')
        .map(|value| {
            if value == "." {
                Ok(None)
            } else {
                value.parse::<i32>().map(Some).map_err(Into::into)
            }
        })
        .collect::<Result<Vec<_>>>()?;
    encode_int_vector(output, &values);
    Ok(())
}

fn encode_float_text(output: &mut Vec<u8>, text: &str) -> Result<()> {
    let values = text
        .split(',')
        .map(|value| {
            if value == "." {
                Ok(None)
            } else {
                value.parse::<f32>().map(Some).map_err(Into::into)
            }
        })
        .collect::<Result<Vec<_>>>()?;
    encode_type(output, values.len(), 5);
    for value in values {
        output.extend_from_slice(&value.map_or(FLOAT_MISSING, f32::to_bits).to_le_bytes());
    }
    Ok(())
}

fn encode_gt_cells(output: &mut Vec<u8>, cells: &[&str]) -> Result<()> {
    let width = cells
        .iter()
        .map(|cell| cell.split(['/', '|']).count())
        .max()
        .unwrap_or(1);
    encode_type(output, width, 1);
    for cell in cells {
        let alleles: Vec<&str> = cell.split(['/', '|']).collect();
        let separators: Vec<char> = cell.chars().filter(|ch| matches!(ch, '/' | '|')).collect();
        for index in 0..width {
            let encoded = if let Some(allele) = alleles.get(index) {
                if *allele == "." {
                    0
                } else {
                    ((allele.parse::<i8>()? + 1) << 1)
                        | i8::from(index > 0 && separators.get(index - 1) == Some(&'|'))
                }
            } else {
                INT8_END
            };
            output.push(encoded as u8);
        }
    }
    Ok(())
}

fn encode_integer_cells(output: &mut Vec<u8>, cells: &[&str]) -> Result<()> {
    let parsed = cells
        .iter()
        .map(|cell| {
            cell.split(',')
                .map(|value| {
                    if value == "." {
                        Ok(None)
                    } else {
                        value.parse::<i32>().map(Some).map_err(Into::into)
                    }
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let width = parsed.iter().map(Vec::len).max().unwrap_or(1);
    encode_type(output, width, 3);
    for values in parsed {
        for index in 0..width {
            let value = values
                .get(index)
                .copied()
                .unwrap_or(Some(INT32_END))
                .unwrap_or(INT32_MISSING);
            output.extend_from_slice(&value.to_le_bytes());
        }
    }
    Ok(())
}

fn encode_float_cells(output: &mut Vec<u8>, cells: &[&str]) -> Result<()> {
    let parsed = cells
        .iter()
        .map(|cell| {
            cell.split(',')
                .map(|value| {
                    if value == "." {
                        Ok(None)
                    } else {
                        value.parse::<f32>().map(Some).map_err(Into::into)
                    }
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    let width = parsed.iter().map(Vec::len).max().unwrap_or(1);
    encode_type(output, width, 5);
    for values in parsed {
        for index in 0..width {
            let bits = values
                .get(index)
                .map_or(FLOAT_END, |value| value.map_or(FLOAT_MISSING, f32::to_bits));
            output.extend_from_slice(&bits.to_le_bytes());
        }
    }
    Ok(())
}

fn encode_string_cells(output: &mut Vec<u8>, cells: &[&str]) {
    let width = cells
        .iter()
        .map(|cell| cell.len())
        .max()
        .unwrap_or(1)
        .max(1);
    encode_type(output, width, 7);
    for cell in cells {
        output.extend_from_slice(cell.as_bytes());
        output.resize(output.len() + width.saturating_sub(cell.len()), 0);
    }
}

pub(crate) fn write_csi(
    path: &Path,
    chunks: &[Option<(u64, u64)>],
    record_counts: &[u64],
) -> Result<()> {
    if chunks.len() != record_counts.len() {
        bail!("CSI chunks and record counts have different reference counts");
    }
    let mut payload = Vec::new();
    payload.extend_from_slice(b"CSI\x01");
    payload.extend_from_slice(&14i32.to_le_bytes());
    payload.extend_from_slice(&5i32.to_le_bytes());
    payload.extend_from_slice(&0i32.to_le_bytes());
    payload.extend_from_slice(&(chunks.len() as i32).to_le_bytes());
    const METADATA_BIN: u32 = 37_450;
    for (chunk, record_count) in chunks.iter().zip(record_counts) {
        payload.extend_from_slice(&(if chunk.is_some() { 2_i32 } else { 0_i32 }).to_le_bytes());
        if let Some((start, end)) = chunk {
            payload.extend_from_slice(&0u32.to_le_bytes());
            payload.extend_from_slice(&start.to_le_bytes());
            payload.extend_from_slice(&1i32.to_le_bytes());
            payload.extend_from_slice(&start.to_le_bytes());
            payload.extend_from_slice(&end.to_le_bytes());

            // CSI inherits the BAI metadata pseudo-bin. Its first chunk is
            // the reference's virtual-offset range; its second stores mapped
            // and unmapped record counts. `bcftools index --stats` requires
            // this metadata even when ordinary region queries already work.
            payload.extend_from_slice(&METADATA_BIN.to_le_bytes());
            payload.extend_from_slice(&0u64.to_le_bytes());
            payload.extend_from_slice(&2i32.to_le_bytes());
            payload.extend_from_slice(&start.to_le_bytes());
            payload.extend_from_slice(&end.to_le_bytes());
            payload.extend_from_slice(&record_count.to_le_bytes());
            payload.extend_from_slice(&0u64.to_le_bytes());
        }
    }
    payload.extend_from_slice(&0u64.to_le_bytes());
    let mut writer = bgzf::io::Writer::new(File::create(path)?);
    writer.write_all(&payload)?;
    writer.finish()?.sync_all()?;
    Ok(())
}

#[cfg(test)]
pub(crate) fn read_uncompressed(path: &Path) -> Result<Vec<u8>> {
    let bytes = fs::read(path)?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut reader = bgzf::io::Reader::new(bytes.as_slice());
        let mut decoded = Vec::new();
        reader.read_to_end(&mut decoded)?;
        Ok(decoded)
    } else {
        Ok(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn strips_only_idx_header_attributes() {
        assert_eq!(
            strip_idx_attribute("##INFO=<ID=DP,Number=1,Type=Integer,IDX=4>"),
            "##INFO=<ID=DP,Number=1,Type=Integer>"
        );
        assert_eq!(strip_idx_attribute("##source=IDX=4"), "##source=IDX=4");
    }

    #[test]
    fn bcf_round_trip_preserves_raw_records_and_writes_csi() -> Result<()> {
        let directory = tempdir()?;
        let output = directory.path().join("roundtrip.bcf");
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "##INFO=<ID=TAG,Number=1,Type=String,Description=\"tag\">".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">".to_string(),
            "##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"GQ\">".to_string(),
            "##FORMAT=<ID=ST,Number=1,Type=String,Description=\"String\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB".to_string(),
        ];
        let records = vec![RawVcfRecord {
            chrom: "chr1".to_string(),
            pos: 2,
            id: "rs1".to_string(),
            ref_allele: "A".to_string(),
            alt_allele: "C,G".to_string(),
            qual: "42.5".to_string(),
            filter: "PASS".to_string(),
            info: "TAG=value".to_string(),
            format: Some("GT:AD:GQ:ST".to_string()),
            samples: vec![
                "1/2:3,4,5:20.5:.".to_string(),
                "0/1:8,9,0:.:value".to_string(),
            ],
        }];
        write(&output, &headers, &records)?;
        let data = read_uncompressed(&output)?;
        assert!(data.starts_with(BCF_MAGIC));
        let (decoded_headers, decoded_records) = decode(&data, &output)?;
        assert_eq!(decoded_headers, headers);
        assert_eq!(decoded_records.len(), 1);
        assert_eq!(decoded_records[0].to_line(), records[0].to_line());
        let csi = output.with_extension("bcf.csi");
        assert!(csi.metadata()?.len() > 0);
        let indexed = read_indexed_records(&output, &csi)?;
        assert_eq!(indexed["chr1"].len(), 1);
        assert_eq!(indexed["chr1"][0].to_line(), records[0].to_line());
        Ok(())
    }

    #[test]
    fn encoder_failure_preserves_existing_bcf_generation() -> Result<()> {
        let directory = tempdir()?;
        let output = directory.path().join("preserve.bcf");
        let index = output.with_extension("bcf.csi");
        fs::write(&output, "old-bcf")?;
        fs::write(&index, "old-csi")?;
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let record = RawVcfRecord::from_line("chr2\t1\t.\tA\tC\t.\tPASS\t.", Path::new("in.vcf"))?;
        let error = write(&output, &headers, &[record]).expect_err("unknown contig must fail");
        assert!(error.to_string().contains(&output.display().to_string()));
        assert_eq!(fs::read_to_string(output)?, "old-bcf");
        assert_eq!(fs::read_to_string(index)?, "old-csi");
        Ok(())
    }

    #[test]
    fn injected_bcf_operations_preserve_generation_and_cleanup() -> Result<()> {
        let directory = tempdir()?;
        let output = directory.path().join("preserve.bcf");
        let index = output.with_extension("bcf.csi");
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let record = RawVcfRecord::from_line("chr1\t1\t.\tA\tC\t.\tPASS\t.", Path::new("in.vcf"))?;
        for operation in [
            FailureOperation::Writer,
            FailureOperation::Encoder,
            FailureOperation::Index,
        ] {
            fs::write(&output, "old-bcf")?;
            fs::write(&index, "old-csi")?;
            crate::output::set_failure_operation(Some(operation));
            let error = write(&output, &headers, std::slice::from_ref(&record))
                .expect_err("injected BCF operation must fail");
            crate::output::set_failure_operation(None);
            let logical = if operation == FailureOperation::Index {
                &index
            } else {
                &output
            };
            assert!(error.to_string().contains(&logical.display().to_string()));
            assert_eq!(fs::read_to_string(&output)?, "old-bcf");
            assert_eq!(fs::read_to_string(&index)?, "old-csi");
            assert_eq!(fs::read_dir(directory.path())?.count(), 2);
        }
        Ok(())
    }

    fn minimal_header() -> Vec<u8> {
        b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\0".to_vec()
    }

    #[test]
    fn rejects_oversized_lengths_before_allocation() {
        let path = Path::new("adversarial.bcf");
        let mut oversized_header = BCF_MAGIC.to_vec();
        oversized_header.extend_from_slice(&((MAX_BCF_HEADER_BYTES as u32) + 1).to_le_bytes());
        let error = decode(&oversized_header, path).unwrap_err().to_string();
        assert!(error.contains("header length"), "{error}");

        let header = minimal_header();
        let mut oversized_record = BCF_MAGIC.to_vec();
        oversized_record.extend_from_slice(&(header.len() as u32).to_le_bytes());
        oversized_record.extend_from_slice(&header);
        oversized_record.extend_from_slice(&(MAX_BCF_RECORD_BYTES as u32).to_le_bytes());
        oversized_record.extend_from_slice(&1u32.to_le_bytes());
        let error = decode(&oversized_record, path).unwrap_err().to_string();
        assert!(error.contains("record 1"), "{error}");
        assert!(error.contains("maximum combined length"), "{error}");
    }

    #[test]
    fn rejects_truncated_records_with_context() {
        let header = minimal_header();
        let mut data = BCF_MAGIC.to_vec();
        data.extend_from_slice(&(header.len() as u32).to_le_bytes());
        data.extend_from_slice(&header);
        data.extend_from_slice(&24u32.to_le_bytes());
        data.extend_from_slice(&0u32.to_le_bytes());
        data.extend_from_slice(&[0; 5]);
        let error = decode(&data, Path::new("truncated.bcf"))
            .unwrap_err()
            .to_string();
        assert!(error.contains("truncated BCF value"), "{error}");
    }

    #[test]
    fn adversarial_truncation_corpus_never_panics() -> Result<()> {
        let directory = tempdir()?;
        let output = directory.path().join("seed.bcf");
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS".to_string(),
        ];
        let records = vec![RawVcfRecord {
            chrom: "chr1".into(),
            pos: 1,
            id: ".".into(),
            ref_allele: "A".into(),
            alt_allele: "C".into(),
            qual: "1".into(),
            filter: ".".into(),
            info: ".".into(),
            format: Some("GT".into()),
            samples: vec!["0/1".into()],
        }];
        write(&output, &headers, &records)?;
        let seed = read_uncompressed(&output)?;
        for end in 0..seed.len() {
            let result = std::panic::catch_unwind(|| decode(&seed[..end], Path::new("fuzz.bcf")));
            assert!(result.is_ok(), "decoder panicked for prefix length {end}");
            let streaming = std::panic::catch_unwind(|| -> Result<()> {
                let input = std::io::Cursor::new(seed[..end].to_vec());
                let mut reader = RecordReader::new(Box::new(input), Path::new("fuzz.bcf"))?;
                while reader.next_record()?.is_some() {}
                Ok(())
            });
            assert!(
                streaming.is_ok(),
                "streaming decoder panicked for prefix length {end}"
            );
        }
        Ok(())
    }

    #[test]
    fn adversarial_coordinate_and_count_fields_return_errors_without_panics() -> Result<()> {
        let directory = tempdir()?;
        let output = directory.path().join("seed.bcf");
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=10>".to_string(),
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS".to_string(),
        ];
        let records = vec![RawVcfRecord {
            chrom: "chr1".into(),
            pos: 1,
            id: ".".into(),
            ref_allele: "A".into(),
            alt_allele: "C".into(),
            qual: "1".into(),
            filter: ".".into(),
            info: ".".into(),
            format: Some("GT".into()),
            samples: vec!["0/1".into()],
        }];
        write(&output, &headers, &records)?;
        let seed = read_uncompressed(&output)?;
        let header_len = u32::from_le_bytes(seed[5..9].try_into().unwrap()) as usize;
        let shared = 9 + header_len + 8;
        let mutations = [
            (shared, (-1_i32).to_le_bytes()),
            (shared + 4, (-1_i32).to_le_bytes()),
            (shared + 8, (-1_i32).to_le_bytes()),
            (shared + 16, 0_u32.to_le_bytes()),
            (shared + 20, 0x00ff_ffff_u32.to_le_bytes()),
        ];
        for (offset, replacement) in mutations {
            let mut mutated = seed.clone();
            mutated[offset..offset + 4].copy_from_slice(&replacement);
            let decoded =
                std::panic::catch_unwind(|| decode(&mutated, Path::new("mutated-counts.bcf")));
            assert!(decoded.is_ok(), "slice decoder panicked at byte {offset}");
            assert!(
                decoded.unwrap().is_err(),
                "mutation at byte {offset} was accepted"
            );

            let streamed = std::panic::catch_unwind(|| -> Result<()> {
                let input = std::io::Cursor::new(mutated);
                let mut reader =
                    RecordReader::new(Box::new(input), Path::new("mutated-counts.bcf"))?;
                while reader.next_record()?.is_some() {}
                Ok(())
            });
            assert!(
                streamed.is_ok(),
                "streaming decoder panicked at byte {offset}"
            );
            assert!(
                streamed.unwrap().is_err(),
                "stream mutation at byte {offset} was accepted"
            );
        }
        Ok(())
    }
}
