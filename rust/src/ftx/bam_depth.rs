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

const BAM_MAGIC: &[u8; 4] = b"BAM\x01";
const CORE_SIZE: usize = 32;
const MAX_SAMPLED_READS: u64 = 10_001;

#[derive(Debug)]
struct ReferenceStats {
    name: String,
    length: u64,
    mapped: u64,
    sampled_reads: u64,
    sampled_bases: u64,
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
        let file = File::open(path).with_context(|| format!("failed to open BAM {path}"))?;
        let reader = bgzf::io::Reader::new(file);
        for reference in scan_bam(reader).with_context(|| format!("failed to read BAM {path}"))? {
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

fn scan_bam<R: Read>(mut reader: R) -> Result<Vec<ReferenceStats>> {
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
        });
    }

    while let Some(block_size) = read_optional_i32(&mut reader)? {
        if block_size < CORE_SIZE as i32 {
            bail!("BAM alignment block is shorter than its core");
        }
        let mut block = vec![0; block_size as usize];
        reader.read_exact(&mut block)?;
        validate_record_layout(&block)?;

        let reference_id = i32::from_le_bytes(block[0..4].try_into().unwrap());
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
        if reference.sampled_reads < MAX_SAMPLED_READS {
            reference.sampled_reads += 1;
            reference.sampled_bases += query_len as u64;
        }
    }

    Ok(references)
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
            paths.push(path.display().to_string());
        }
        let depths = normalization_depths(&paths).unwrap();
        assert!((depths["chr1"] - 1.05).abs() <= f64::EPSILON * 2.0);
    }

    #[test]
    fn rejects_truncated_alignment_blocks() {
        let mut bytes = bam_bytes(&[10]);
        bytes.pop();
        assert!(scan_bam(Cursor::new(bytes)).is_err());
    }
}
