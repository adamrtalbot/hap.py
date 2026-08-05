//! Minimal streaming BAM statistics for legacy FTX depth normalization.
//!
//! `Tools.bamstats.bamStats` computes `mapped * mean(read.rlen) / contig_len`,
//! sampling at most the first 10,001 mapped reads per contig. FTX averages
//! that coverage across repeated `--bam` inputs and multiplies it by three.

use anyhow::{Context, Result, bail};
use noodles_bgzf as bgzf;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, Read};
use std::path::{Path, PathBuf};

const BAM_MAGIC: &[u8; 4] = b"BAM\x01";
const BAI_MAGIC: &[u8; 4] = b"BAI\x01";
const BAI_METADATA_BIN: u32 = 37_450;
const CORE_SIZE: usize = 32;
const MAX_SAMPLED_READS: u64 = 10_001;

#[derive(Debug)]
struct IndexedRecord {
    bin: u32,
    virtual_position: u64,
}

#[derive(Debug)]
struct ReferenceStats {
    name: String,
    length: u64,
    mapped: u64,
    sampled_reads: u64,
    sampled_bases: u64,
    indexed_records: Vec<IndexedRecord>,
}

impl ReferenceStats {
    fn coverage(&self) -> f64 {
        if self.length == 0 || self.sampled_reads == 0 {
            0.0
        } else {
            let mean_read_len = self.sampled_bases as f64 / self.sampled_reads as f64;
            self.mapped as f64 * mean_read_len / self.length as f64
        }
    }
}

pub(super) fn normalization_depths(paths: &[String]) -> Result<BTreeMap<String, f64>> {
    let mut coverages: BTreeMap<String, (f64, usize)> = BTreeMap::new();
    for path in paths {
        // Legacy's indexed fetch accepts either common BAI filename and
        // leaves every reference at zero when the index is missing or cannot
        // retrieve the BAM records. htslib reports the index problem to
        // stderr, but bamstats catches the fetch failure and still publishes
        // a zero-depth feature table.
        let appended_index = PathBuf::from(format!("{path}.bai"));
        let stem_index = Path::new(path).with_extension("bai");
        let index = if appended_index.is_file() {
            Some(appended_index)
        } else if stem_index.is_file() {
            Some(stem_index)
        } else {
            None
        };
        let file = File::open(path).with_context(|| format!("failed to open BAM {path}"))?;
        let reader = bgzf::io::Reader::new(file);
        let references = if let Some(index) = index {
            let references =
                scan_indexed_bam(reader).with_context(|| format!("failed to read BAM {path}"))?;
            if validate_bai(Path::new(path), &index, &references).is_ok() {
                references
            } else {
                zero_depths(references)
            }
        } else {
            read_bam_references(reader).with_context(|| format!("failed to read BAM {path}"))?
        };
        for reference in references {
            let coverage = reference.coverage();
            let entry = coverages.entry(reference.name).or_default();
            entry.0 += coverage;
            entry.1 += 1;
        }
    }

    Ok(coverages
        .into_iter()
        .map(|(chrom, (sum, count))| (chrom, (sum / count as f64) * 3.0))
        .collect())
}

#[cfg(test)]
fn scan_bam<R: Read>(reader: R) -> Result<Vec<ReferenceStats>> {
    read_bam(reader, true, |_| None)
}

fn scan_indexed_bam(reader: bgzf::io::Reader<File>) -> Result<Vec<ReferenceStats>> {
    read_bam(reader, true, |reader| {
        Some(u64::from(reader.virtual_position()))
    })
}

fn read_bam_references<R: Read>(reader: R) -> Result<Vec<ReferenceStats>> {
    read_bam(reader, false, |_| None)
}

fn read_bam<R: Read, F: FnMut(&R) -> Option<u64>>(
    mut reader: R,
    scan_alignments: bool,
    mut virtual_position: F,
) -> Result<Vec<ReferenceStats>> {
    let mut magic = [0; 4];
    reader.read_exact(&mut magic)?;
    if &magic != BAM_MAGIC {
        bail!("invalid BAM magic");
    }

    let header_len = read_nonnegative_i32(&mut reader, "header length")?;
    discard(&mut reader, header_len)?;
    let reference_count = read_nonnegative_i32(&mut reader, "reference count")?;
    let mut references = Vec::with_capacity(reference_count);
    for _ in 0..reference_count {
        let name_len = read_nonnegative_i32(&mut reader, "reference name length")?;
        if name_len == 0 {
            bail!("BAM reference name is empty");
        }
        let mut encoded_name = vec![0; name_len];
        reader.read_exact(&mut encoded_name)?;
        if encoded_name.pop() != Some(0) {
            bail!("BAM reference name is not NUL-terminated");
        }
        let name = String::from_utf8(encoded_name).context("BAM reference name is not UTF-8")?;
        let length = read_nonnegative_i32(&mut reader, "reference length")? as u64;
        references.push(ReferenceStats {
            name,
            length,
            mapped: 0,
            sampled_reads: 0,
            sampled_bases: 0,
            indexed_records: Vec::new(),
        });
    }

    if !scan_alignments {
        return Ok(references);
    }

    loop {
        let record_position = virtual_position(&reader);
        let Some(block_size) = read_optional_i32(&mut reader)? else {
            break;
        };
        if block_size < CORE_SIZE as i32 {
            bail!("BAM alignment block is shorter than its core");
        }
        let mut block = vec![0; block_size as usize];
        reader.read_exact(&mut block)?;
        validate_record_layout(&block)?;

        let reference_id = i32::from_le_bytes(block[0..4].try_into().unwrap());
        let bin_mq_name = u32::from_le_bytes(block[8..12].try_into().unwrap());
        let bin = bin_mq_name >> 16;
        let flag_and_cigar = u32::from_le_bytes(block[12..16].try_into().unwrap());
        let flags = (flag_and_cigar >> 16) as u16;
        let query_len = i32::from_le_bytes(block[16..20].try_into().unwrap());
        if query_len < 0 {
            bail!("BAM alignment has a negative query length");
        }
        if flags & 0x4 != 0 || reference_id < 0 {
            continue;
        }
        let Some(reference) = references.get_mut(reference_id as usize) else {
            bail!("BAM alignment references unknown sequence id {reference_id}");
        };
        reference.mapped += 1;
        if let Some(virtual_position) = record_position {
            reference.indexed_records.push(IndexedRecord {
                bin,
                virtual_position,
            });
        }
        if reference.sampled_reads < MAX_SAMPLED_READS {
            reference.sampled_reads += 1;
            reference.sampled_bases += query_len as u64;
        }
    }

    Ok(references)
}

fn zero_depths(mut references: Vec<ReferenceStats>) -> Vec<ReferenceStats> {
    for reference in &mut references {
        reference.mapped = 0;
        reference.sampled_reads = 0;
        reference.sampled_bases = 0;
        reference.indexed_records.clear();
    }
    references
}

fn validate_bai(bam: &Path, index: &Path, references: &[ReferenceStats]) -> Result<()> {
    let bytes = std::fs::read(index)
        .with_context(|| format!("failed to read BAM index {}", index.display()))?;
    let mut cursor = BaiCursor::new(&bytes);
    if cursor.take(BAI_MAGIC.len())? != BAI_MAGIC {
        bail!("invalid BAI magic");
    }
    let reference_count = cursor.nonnegative_i32("reference count")?;
    if reference_count != references.len() {
        bail!(
            "BAI reference count {} does not match BAM reference count {}",
            reference_count,
            references.len()
        );
    }
    let bam_len = std::fs::metadata(bam)?.len();

    for reference in references {
        let bin_count = cursor.nonnegative_i32("bin count")?;
        let mut bins: BTreeMap<u32, Vec<(u64, u64)>> = BTreeMap::new();
        let mut mapped_count = None;
        for _ in 0..bin_count {
            let bin = cursor.u32("bin")?;
            let chunk_count = cursor.nonnegative_i32("chunk count")?;
            let mut chunks = Vec::with_capacity(chunk_count);
            for chunk_index in 0..chunk_count {
                let start = cursor.u64("chunk start")?;
                let end = cursor.u64("chunk end")?;
                if bin == BAI_METADATA_BIN && chunk_index == 1 {
                    mapped_count = Some(start);
                } else if bin != BAI_METADATA_BIN {
                    validate_virtual_offset(start, bam_len)?;
                    validate_virtual_offset(end, bam_len)?;
                    if start >= end {
                        bail!("BAI chunk start does not precede its end");
                    }
                    chunks.push((start, end));
                }
            }
            if bin != BAI_METADATA_BIN && bins.insert(bin, chunks).is_some() {
                bail!("BAI contains duplicate bin {bin}");
            }
        }

        let interval_count = cursor.nonnegative_i32("linear interval count")?;
        for _ in 0..interval_count {
            let offset = cursor.u64("linear interval offset")?;
            if offset != 0 {
                validate_virtual_offset(offset, bam_len)?;
            }
        }

        if mapped_count != Some(reference.mapped) {
            bail!(
                "BAI mapped count for {} is {:?}, BAM contains {}",
                reference.name,
                mapped_count,
                reference.mapped
            );
        }
        for record in &reference.indexed_records {
            let Some(chunks) = bins.get(&record.bin) else {
                bail!(
                    "BAI bin {} omits a mapped record on {}",
                    record.bin,
                    reference.name
                );
            };
            if !chunks.iter().any(|(start, end)| {
                *start <= record.virtual_position && record.virtual_position < *end
            }) {
                bail!(
                    "BAI chunks omit a mapped record on {} at virtual offset {}",
                    reference.name,
                    record.virtual_position
                );
            }
        }
    }

    match cursor.remaining() {
        0 => {}
        8 => {
            cursor.u64("unplaced-unmapped count")?;
        }
        remaining => bail!("BAI has {remaining} unexpected trailing bytes"),
    }
    Ok(())
}

fn validate_virtual_offset(offset: u64, bam_len: u64) -> Result<()> {
    let compressed = offset >> 16;
    if compressed >= bam_len {
        bail!("BAI virtual offset points beyond the BAM");
    }
    Ok(())
}

struct BaiCursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> BaiCursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.offset)
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8]> {
        let end = self
            .offset
            .checked_add(length)
            .context("BAI offset overflow")?;
        let value = self
            .bytes
            .get(self.offset..end)
            .context("BAI is truncated")?;
        self.offset = end;
        Ok(value)
    }

    fn u32(&mut self, field: &str) -> Result<u32> {
        Ok(u32::from_le_bytes(
            self.take(4)?
                .try_into()
                .with_context(|| format!("invalid BAI {field}"))?,
        ))
    }

    fn u64(&mut self, field: &str) -> Result<u64> {
        Ok(u64::from_le_bytes(
            self.take(8)?
                .try_into()
                .with_context(|| format!("invalid BAI {field}"))?,
        ))
    }

    fn nonnegative_i32(&mut self, field: &str) -> Result<usize> {
        let value = i32::from_le_bytes(
            self.take(4)?
                .try_into()
                .with_context(|| format!("invalid BAI {field}"))?,
        );
        if value < 0 {
            bail!("BAI {field} is negative");
        }
        Ok(value as usize)
    }
}

fn validate_record_layout(block: &[u8]) -> Result<()> {
    let bin_mq_name = u32::from_le_bytes(block[8..12].try_into().unwrap());
    let flag_and_cigar = u32::from_le_bytes(block[12..16].try_into().unwrap());
    let read_name_len = (bin_mq_name & 0xff) as usize;
    let cigar_count = (flag_and_cigar & 0xffff) as usize;
    let query_len = i32::from_le_bytes(block[16..20].try_into().unwrap());
    if query_len < 0 {
        bail!("BAM alignment has a negative query length");
    }
    let query_len = query_len as usize;
    let required = CORE_SIZE
        .checked_add(read_name_len)
        .and_then(|size| size.checked_add(cigar_count.checked_mul(4)?))
        .and_then(|size| size.checked_add(query_len.div_ceil(2)))
        .and_then(|size| size.checked_add(query_len))
        .context("BAM alignment block size overflow")?;
    if required > block.len() {
        bail!("BAM alignment block is truncated");
    }
    if read_name_len > 0 && block[CORE_SIZE + read_name_len - 1] != 0 {
        bail!("BAM read name is not NUL-terminated");
    }

    let cigar_start = CORE_SIZE + read_name_len;
    for index in 0..cigar_count {
        let offset = cigar_start + index * 4;
        let operation = u32::from_le_bytes(block[offset..offset + 4].try_into().unwrap());
        if operation & 0x0f > 9 {
            bail!("BAM alignment uses unknown CIGAR operation");
        }
    }
    Ok(())
}

fn read_optional_i32<R: Read>(reader: &mut R) -> io::Result<Option<i32>> {
    let mut bytes = [0; 4];
    let mut read = 0;
    while read < bytes.len() {
        match reader.read(&mut bytes[read..])? {
            0 if read == 0 => return Ok(None),
            0 => return Err(io::Error::from(io::ErrorKind::UnexpectedEof)),
            count => read += count,
        }
    }
    Ok(Some(i32::from_le_bytes(bytes)))
}

fn read_nonnegative_i32<R: Read>(reader: &mut R, field: &str) -> Result<usize> {
    let value = read_optional_i32(reader)?.with_context(|| format!("BAM ended before {field}"))?;
    if value < 0 {
        bail!("BAM {field} is negative");
    }
    Ok(value as usize)
}

fn discard<R: Read>(reader: &mut R, length: usize) -> io::Result<()> {
    io::copy(&mut reader.take(length as u64), &mut io::sink()).and_then(|read| {
        if read == length as u64 {
            Ok(())
        } else {
            Err(io::Error::from(io::ErrorKind::UnexpectedEof))
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Write};

    fn push_i32(bytes: &mut Vec<u8>, value: i32) {
        bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn push_record(bytes: &mut Vec<u8>, reference_id: i32, flags: u16, query_len: usize) {
        let name = b"r\0";
        let block_size = CORE_SIZE + name.len() + 4 + query_len.div_ceil(2) + query_len;
        push_i32(bytes, block_size as i32);
        push_i32(bytes, reference_id);
        push_i32(bytes, 0);
        bytes.extend_from_slice(&(name.len() as u32).to_le_bytes());
        bytes.extend_from_slice(&(((flags as u32) << 16) | 1).to_le_bytes());
        push_i32(bytes, query_len as i32);
        bytes.extend_from_slice(
            &[-1i32 as u32, -1i32 as u32, 0]
                .map(u32::to_le_bytes)
                .concat(),
        );
        bytes.extend_from_slice(name);
        bytes.extend_from_slice(&((query_len as u32) << 4).to_le_bytes());
        bytes.resize(bytes.len() + query_len.div_ceil(2) + query_len, 0);
    }

    fn bam_bytes(read_lengths: &[usize]) -> Vec<u8> {
        let mut bytes = BAM_MAGIC.to_vec();
        push_i32(&mut bytes, 0);
        push_i32(&mut bytes, 1);
        push_i32(&mut bytes, 5);
        bytes.extend_from_slice(b"chr1\0");
        push_i32(&mut bytes, 100);
        for length in read_lengths {
            push_record(&mut bytes, 0, 0, *length);
        }
        push_record(&mut bytes, 0, 0x4, 75);
        bytes
    }

    fn push_u32(bytes: &mut Vec<u8>, value: u32) {
        bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn push_u64(bytes: &mut Vec<u8>, value: u64) {
        bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn write_test_bai(bam: &Path, mapped_override: Option<u64>, bin_override: Option<u32>) {
        let reader = bgzf::io::Reader::new(File::open(bam).unwrap());
        let references = scan_indexed_bam(reader).unwrap();
        let bam_len = std::fs::metadata(bam).unwrap().len();
        let chunk_end = ((bam_len - 1) << 16) | u16::MAX as u64;
        let mut bytes = BAI_MAGIC.to_vec();
        push_i32(&mut bytes, references.len() as i32);
        for reference in references {
            let mut bins: BTreeMap<u32, u64> = BTreeMap::new();
            for record in &reference.indexed_records {
                bins.entry(bin_override.unwrap_or(record.bin))
                    .and_modify(|start| *start = (*start).min(record.virtual_position))
                    .or_insert(record.virtual_position);
            }
            push_i32(&mut bytes, (bins.len() + 1) as i32);
            for (bin, start) in &bins {
                push_u32(&mut bytes, *bin);
                push_i32(&mut bytes, 1);
                push_u64(&mut bytes, *start);
                push_u64(&mut bytes, chunk_end);
            }
            push_u32(&mut bytes, BAI_METADATA_BIN);
            push_i32(&mut bytes, 2);
            push_u64(&mut bytes, bins.values().next().copied().unwrap_or(0));
            push_u64(&mut bytes, chunk_end);
            push_u64(&mut bytes, mapped_override.unwrap_or(reference.mapped));
            push_u64(&mut bytes, 0);
            push_i32(&mut bytes, usize::from(!bins.is_empty()) as i32);
            if let Some(start) = bins.values().next() {
                push_u64(&mut bytes, *start);
            }
        }
        push_u64(&mut bytes, 0);
        std::fs::write(format!("{}.bai", bam.display()), bytes).unwrap();
    }

    #[test]
    fn scans_mapped_count_and_mean_query_length() {
        let stats = scan_bam(Cursor::new(bam_bytes(&[10, 20]))).unwrap();
        assert_eq!(stats[0].mapped, 2);
        assert_eq!(stats[0].sampled_reads, 2);
        assert_eq!(stats[0].sampled_bases, 30);
        assert_eq!(stats[0].coverage(), 0.3);
    }

    #[test]
    fn repeated_bgzf_bams_are_averaged_then_tripled() {
        let directory = tempfile::tempdir().unwrap();
        let mut paths = Vec::new();
        for (index, lengths) in [[10, 20].as_slice(), [20, 20].as_slice()]
            .iter()
            .enumerate()
        {
            let path = directory.path().join(format!("{index}.bam"));
            let file = File::create(&path).unwrap();
            let mut writer = bgzf::io::Writer::new(file);
            writer.write_all(&bam_bytes(lengths)).unwrap();
            writer.finish().unwrap();
            if index == 0 {
                write_test_bai(&path, None, None);
            } else {
                write_test_bai(&path, None, None);
                std::fs::rename(
                    format!("{}.bai", path.display()),
                    path.with_extension("bai"),
                )
                .unwrap();
            }
            paths.push(path.display().to_string());
        }
        let depths = normalization_depths(&paths).unwrap();
        assert!((depths["chr1"] - 1.05).abs() <= f64::EPSILON * 2.0);
    }

    #[test]
    fn missing_bai_yields_zero_coverage() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.bam");
        let file = File::create(&path).unwrap();
        let mut writer = bgzf::io::Writer::new(file);
        writer.write_all(&bam_bytes(&[20, 20])).unwrap();
        writer.finish().unwrap();

        let depths = normalization_depths(&[path.display().to_string()]).unwrap();

        assert_eq!(depths, BTreeMap::from([("chr1".to_string(), 0.0)]));
    }

    #[test]
    fn corrupt_existing_bai_yields_zero_coverage_like_legacy() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.bam");
        let file = File::create(&path).unwrap();
        let mut writer = bgzf::io::Writer::new(file);
        writer.write_all(&bam_bytes(&[20, 20])).unwrap();
        writer.finish().unwrap();
        std::fs::write(format!("{}.bai", path.display()), b"BAI\x01truncated").unwrap();

        let depths = normalization_depths(&[path.display().to_string()]).unwrap();

        assert_eq!(depths, BTreeMap::from([("chr1".to_string(), 0.0)]));
    }

    #[test]
    fn stale_bai_counts_yield_zero_coverage_like_legacy() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.bam");
        let file = File::create(&path).unwrap();
        let mut writer = bgzf::io::Writer::new(file);
        writer.write_all(&bam_bytes(&[20, 20])).unwrap();
        writer.finish().unwrap();
        write_test_bai(&path, Some(1), None);

        let depths = normalization_depths(&[path.display().to_string()]).unwrap();

        assert_eq!(depths, BTreeMap::from([("chr1".to_string(), 0.0)]));
    }

    #[test]
    fn wrong_bai_bin_yields_zero_coverage_like_legacy() {
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("reads.bam");
        let file = File::create(&path).unwrap();
        let mut writer = bgzf::io::Writer::new(file);
        writer.write_all(&bam_bytes(&[20, 20])).unwrap();
        writer.finish().unwrap();
        write_test_bai(&path, None, Some(1));

        let depths = normalization_depths(&[path.display().to_string()]).unwrap();

        assert_eq!(depths, BTreeMap::from([("chr1".to_string(), 0.0)]));
    }

    #[test]
    fn rejects_truncated_alignment_blocks() {
        let mut bytes = bam_bytes(&[10]);
        bytes.pop();
        assert!(scan_bam(Cursor::new(bytes)).is_err());
    }
}
