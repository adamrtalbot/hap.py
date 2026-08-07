use crate::output::{FailureOperation, OutputTransaction, fail_operation};
use anyhow::{Context, Result, bail};
use flate2::read::MultiGzDecoder;
use noodles_bgzf as bgzf;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File, OpenOptions};
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard};

const TBI_MAX_POSITION: usize = 1 << 29;
const TBI_LINEAR_SHIFT: usize = 14;
const TBI_METADATA_BIN: u32 = 37_450;
static TEMP_FILE_ID: AtomicU64 = AtomicU64::new(0);
static PUBLICATION_MUTEX: Mutex<()> = Mutex::new(());

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct VariantKey {
    pub chrom: String,
    pub pos: usize,
    pub ref_allele: String,
    pub alt_allele: String,
}

#[derive(Clone, Debug)]
pub struct Variant {
    pub key: VariantKey,
    pub qual: String,
    pub filter: String,
    pub gt: String,
}

impl Variant {
    pub fn is_pass(&self) -> bool {
        self.filter == "PASS" || self.filter == "."
    }

    pub fn is_single_base_snp(&self) -> bool {
        let active = self.active_alts();
        if active.is_empty() {
            // Fall back to the full-record shape when no allele is called.
            return self.key.ref_allele.len() == 1
                && self.key.alt_allele.len() == 1
                && self.key.alt_allele != "."
                && !self.key.alt_allele.contains(',');
        }
        self.key.ref_allele.len() == 1
            && active
                .iter()
                .all(|alt| alt.len() == 1 && *alt != "." && !alt.starts_with('<'))
    }

    pub fn primary_type(&self) -> &'static str {
        // Legacy hap.py classifies the VCF record's BVT by the allele that
        // the GT actually calls, not by scanning every ALT in a multi-allelic
        // record. A `T → C,TATC` truth with `GT=1|0` has only the `T→C` SNP
        // in play; the unused `TATC` alternate does NOT promote the record
        // to INDEL. Use `active_alts` so our classification matches that.
        let active = self.active_alts();
        if active.is_empty() {
            return if alt_type_matches_ref(&self.key.ref_allele, &self.key.alt_allele) {
                "SNP"
            } else {
                "INDEL"
            };
        }
        let ref_len = self.key.ref_allele.len();
        if active
            .iter()
            .all(|alt| *alt != "." && !alt.starts_with('<') && alt.len() == ref_len)
        {
            "SNP"
        } else {
            "INDEL"
        }
    }

    pub fn end_pos(&self) -> usize {
        self.key.pos + self.key.ref_allele.len().saturating_sub(1)
    }

    pub fn is_transition(&self) -> bool {
        let active = self.active_alts();
        if active.is_empty() {
            return matches!(
                (self.key.ref_allele.as_str(), self.key.alt_allele.as_str()),
                ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C")
            );
        }
        active.iter().all(|alt| {
            matches!(
                (self.key.ref_allele.as_str(), *alt),
                ("A", "G") | ("G", "A") | ("C", "T") | ("T", "C")
            )
        })
    }

    /// ALT sequences that the `GT` field actually names. Allele index 0 is
    /// the reference and is elided; unused indexes in a multi-allelic record
    /// are dropped entirely. Returns an empty vector for `GT=0/0`, half-
    /// calls (`./1`), or unparseable genotypes — callers should fall back
    /// to the full-record view in that case.
    pub fn active_alts(&self) -> Vec<&str> {
        let alts: Vec<&str> = self.key.alt_allele.split(',').collect();
        let mut indices: Vec<usize> = self
            .gt
            .split(['/', '|'])
            .filter_map(|value| value.parse::<usize>().ok())
            .filter(|index| *index > 0)
            .collect();
        indices.sort_unstable();
        indices.dedup();
        indices
            .into_iter()
            .filter_map(|index| alts.get(index.saturating_sub(1)).copied())
            .collect()
    }

    pub fn is_homalt(&self) -> bool {
        // Legacy's summary counts treat any diploid GT with two equal,
        // non-zero allele indexes as homalt — `2|2`, `3|3`, ... all count.
        // The truth-side multi-allelic records (e.g. `2|2` or `3|3` after
        // bcftools merge) require the wider rule; query-side records only
        // ever carry `1/1`, so the wider rule is a strict superset for
        // them and parity holds.
        let alleles = parse_gt_usize(&self.gt);
        alleles.len() == 2 && alleles[0] != 0 && alleles[0] == alleles[1]
    }

    pub fn is_het(&self) -> bool {
        // Legacy's summary counts treat a diploid GT as het when exactly
        // one allele is the reference (`0`). That includes `0|2`, `2|0`,
        // `0|3`, `3|0`, etc. — any heterozygous-with-ref configuration.
        // Hetalt (e.g. `1|2`) is intentionally excluded from het AND
        // homalt to mirror legacy's rendering of those rows.
        let alleles = parse_gt_usize(&self.gt);
        if alleles.len() != 2 {
            return false;
        }
        let zero_count = alleles.iter().filter(|a| **a == 0).count();
        zero_count == 1
    }
}

fn parse_gt_usize(gt: &str) -> Vec<usize> {
    gt.split(['/', '|'])
        .map(|part| part.parse::<usize>().unwrap_or(0))
        .collect()
}

fn alt_type_matches_ref(ref_allele: &str, alt_alleles: &str) -> bool {
    let ref_len = ref_allele.len();
    alt_alleles
        .split(',')
        .all(|alt| alt != "." && !alt.starts_with('<') && alt.len() == ref_len)
}

#[derive(Clone, Debug)]
pub struct RawVcfRecord {
    pub chrom: String,
    pub pos: usize,
    pub id: String,
    pub ref_allele: String,
    pub alt_allele: String,
    pub qual: String,
    pub filter: String,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
}

impl RawVcfRecord {
    pub fn from_line(line: &str, path: &Path) -> Result<Self> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 8 {
            bail!("VCF record has fewer than 8 fields in {}", path.display());
        }
        Ok(Self {
            chrom: fields[0].to_string(),
            pos: fields[1].parse::<usize>().with_context(|| {
                format!("invalid position '{}' in {}", fields[1], path.display())
            })?,
            id: fields[2].to_string(),
            ref_allele: fields[3].to_string(),
            alt_allele: fields[4].to_string(),
            qual: fields[5].to_string(),
            filter: fields[6].to_string(),
            info: fields[7].to_string(),
            format: fields.get(8).map(|value| value.to_string()),
            samples: if fields.len() > 9 {
                fields[9..]
                    .iter()
                    .map(|value| (*value).to_string())
                    .collect()
            } else {
                Vec::new()
            },
        })
    }

    pub fn to_line(&self) -> String {
        let mut fields = vec![
            self.chrom.clone(),
            self.pos.to_string(),
            self.id.clone(),
            self.ref_allele.clone(),
            self.alt_allele.clone(),
            self.qual.clone(),
            self.filter.clone(),
            self.info.clone(),
        ];
        if let Some(format) = &self.format {
            fields.push(format.clone());
            fields.extend(self.samples.iter().cloned());
        }
        fields.join("\t")
    }

    pub fn is_pass(&self) -> bool {
        self.filter == "PASS" || self.filter == "." || self.filter.is_empty()
    }

    pub fn end_pos(&self) -> usize {
        self.pos + self.ref_allele.len().saturating_sub(1)
    }

    /// Return the record's 1-based inclusive effective end using HTSlib's VCF
    /// interval rules. Numerically this is also the 0-based exclusive end used
    /// by BED overlap checks. Unlike [`Self::end_pos`], this honors INFO/END,
    /// INFO/SVLEN, and gVCF FORMAT/LEN spans.
    pub fn effective_end_pos(&self, path: &Path) -> Result<usize> {
        if self.pos == 0 {
            bail!("VCF positions must be >= 1 in {}", path.display());
        }
        effective_vcf_end(
            self.pos - 1,
            &self.ref_allele,
            &self.alt_allele,
            &self.info,
            self.format.as_deref(),
            &self.samples.iter().map(String::as_str).collect::<Vec<_>>(),
            path,
        )
    }

    pub fn format_keys(&self) -> Vec<&str> {
        self.format
            .as_deref()
            .map(|format| format.split(':').collect())
            .unwrap_or_default()
    }

    pub fn sample_map(&self, index: usize) -> BTreeMap<String, String> {
        let keys = self.format_keys();
        let Some(sample) = self.samples.get(index) else {
            return BTreeMap::new();
        };
        keys.into_iter()
            .zip(sample.split(':'))
            .map(|(key, value)| (key.to_string(), value.to_string()))
            .collect()
    }
}

pub fn read_text(path: &Path) -> Result<String> {
    let data = read_decoded_bytes(path)?;
    String::from_utf8(data).with_context(|| format!("{} is not valid UTF-8", path.display()))
}

fn read_decoded_bytes(path: &Path) -> Result<Vec<u8>> {
    let bytes = fs::read(path).with_context(|| format!("failed to read {}", path.display()))?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        if is_bgzf(&bytes) {
            let mut reader = bgzf::io::Reader::new(bytes.as_slice());
            let mut data = Vec::new();
            reader.read_to_end(&mut data)?;
            Ok(data)
        } else {
            let mut decoder = MultiGzDecoder::new(bytes.as_slice());
            let mut data = Vec::new();
            decoder.read_to_end(&mut data)?;
            Ok(data)
        }
    } else {
        Ok(bytes)
    }
}

fn is_bgzf(bytes: &[u8]) -> bool {
    // GZIP + FEXTRA flag + BGZF "BC" subfield marker.
    bytes.len() >= 16
        && bytes[0] == 0x1f
        && bytes[1] == 0x8b
        && bytes[2] == 0x08
        && (bytes[3] & 0x04) != 0
        && bytes[12] == b'B'
        && bytes[13] == b'C'
}

pub fn load_raw_vcf(path: &Path) -> Result<(Vec<String>, Vec<RawVcfRecord>)> {
    let data = read_decoded_bytes(path)?;
    if crate::bcf::is_bcf_data(&data) {
        return crate::bcf::decode(&data, path);
    }
    let text = String::from_utf8(data)
        .with_context(|| format!("{} is not valid UTF-8", path.display()))?;
    let mut headers = Vec::new();
    let mut records = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') {
            headers.push(line.to_string());
        } else if !line.trim().is_empty() {
            records.push(RawVcfRecord::from_line(line, path)?);
        }
    }
    Ok((headers, records))
}

pub fn write_raw_vcf(path: &Path, headers: &[String], records: &[RawVcfRecord]) -> Result<()> {
    if path.extension().and_then(|extension| extension.to_str()) == Some("bcf") {
        return crate::bcf::write(path, headers, records);
    }
    let lines: Vec<String> = records.iter().map(RawVcfRecord::to_line).collect();
    if path.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        write_indexed_vcf(path, headers, lines.iter().map(String::as_str))?;
    } else {
        let transaction = OutputTransaction::files(Vec::<PathBuf>::new(), [path])?;
        let temporary_path = transaction.staged_file(path)?.to_path_buf();
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary_path)
            .with_context(|| format!("failed to create VCF destination {}", path.display()))?;
        let mut writer = BufWriter::new(file);
        (|| {
            fail_operation(FailureOperation::Writer, path)?;
            write_vcf_lines(&mut writer, headers, lines.iter().map(String::as_str))?;
            writer.flush()?;
            writer
                .get_ref()
                .sync_all()
                .with_context(|| format!("failed to sync {}", temporary_path.display()))?;
            Ok::<(), anyhow::Error>(())
        })()
        .map_err(|error| {
            anyhow::anyhow!(
                "failed to write VCF destination {}: {error:#}",
                path.display()
            )
        })?;
        transaction.commit()?;
    }
    Ok(())
}

/// Writes a BGZF-compressed VCF and its Tabix v1 index without invoking
/// external executables.
///
/// Record lines are validated before any destination is replaced. Records
/// must be grouped by reference sequence and position-sorted within each
/// sequence, as required by the Tabix format.
pub fn write_indexed_vcf<'a, I>(path: &Path, headers: &[String], lines: I) -> Result<()>
where
    I: IntoIterator<Item = &'a str>,
{
    let sidecar_path = tabix_path(path);
    let transaction = OutputTransaction::files(Vec::<PathBuf>::new(), [path, &sidecar_path])?;
    let staged = transaction.staged_file(path)?.to_path_buf();
    write_indexed_vcf_inner(&staged, path, &sidecar_path, headers, lines).map_err(|error| {
        anyhow::anyhow!(
            "failed to write indexed VCF destination {} (index {}): {error:#}",
            path.display(),
            sidecar_path.display()
        )
    })?;
    transaction.commit()
}

fn write_indexed_vcf_inner<'a, I>(
    path: &Path,
    logical_vcf: &Path,
    logical_index: &Path,
    headers: &[String],
    lines: I,
) -> Result<()>
where
    I: IntoIterator<Item = &'a str>,
{
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create {}", parent.display()))?;
    }

    let records = parse_index_records(lines, path)?;
    let sidecar_path = tabix_path(path);
    let (temporary_vcf, vcf_file) = create_temporary_file(path)?;
    let (temporary_tbi, tbi_file) = match create_temporary_file(&sidecar_path) {
        Ok(value) => value,
        Err(error) => {
            let _ = fs::remove_file(&temporary_vcf);
            return Err(error);
        }
    };

    let result = (|| {
        fail_operation(FailureOperation::Writer, logical_vcf)?;
        let mut writer = bgzf::io::Writer::new(vcf_file);
        for header in headers {
            writeln!(writer, "{header}")?;
        }

        let mut index = TabixIndex::default();
        for record in &records {
            let chunk_start = u64::from(writer.virtual_position());
            writeln!(writer, "{}", record.line)?;
            let chunk_end = u64::from(writer.virtual_position());
            index.push(record, Chunk::new(chunk_start, chunk_end))?;
        }

        fail_operation(FailureOperation::Encoder, logical_vcf)?;
        let vcf_file = writer.finish()?;
        vcf_file
            .sync_all()
            .with_context(|| format!("failed to sync {}", temporary_vcf.display()))?;

        let payload = index.finish();
        fail_operation(FailureOperation::Index, logical_index)?;
        let mut index_writer = bgzf::io::Writer::new(tbi_file);
        index_writer.write_all(&payload)?;
        let tbi_file = index_writer.finish()?;
        tbi_file
            .sync_all()
            .with_context(|| format!("failed to sync {}", temporary_tbi.display()))?;

        publish_pair(&temporary_vcf, path, &temporary_tbi, &sidecar_path)
    })();

    if result.is_err() {
        let _ = fs::remove_file(&temporary_vcf);
        let _ = fs::remove_file(&temporary_tbi);
    }
    result
}

/// Publishes the VCF and index as one recoverable transaction.
///
/// POSIX filesystems do not offer a two-path atomic rename. Keeping any old
/// pair in adjacent backup files lets us roll both destinations back if either
/// publication rename fails, avoiding a VCF paired with an index from a
/// different generation.
fn publish_pair(
    temporary_vcf: &Path,
    destination_vcf: &Path,
    temporary_tbi: &Path,
    destination_tbi: &Path,
) -> Result<()> {
    // A stable lock keyed by the canonical destination coordinates separate
    // hap-rs processes without adding an artifact beside the requested pair.
    // The process-local mutex is also necessary because advisory lock
    // semantics for two descriptors owned by one process vary by platform.
    let _publication_lock = PairPublicationLock::acquire(destination_vcf)?;

    validate_regular_destination(destination_vcf)?;
    validate_regular_destination(destination_tbi)?;

    let backup_vcf = move_existing_to_backup(destination_vcf)?;
    let backup_tbi = match move_existing_to_backup(destination_tbi) {
        Ok(backup) => backup,
        Err(error) => {
            if let Some(backup) = &backup_vcf
                && let Err(restore_error) = fs::rename(backup, destination_vcf)
            {
                return Err(error.context(format!(
                    "also failed to restore {}: {restore_error}",
                    destination_vcf.display()
                )));
            }
            return Err(error);
        }
    };

    let mut published_tbi = false;
    let mut published_vcf = false;
    let publication = (|| {
        fs::rename(temporary_tbi, destination_tbi).with_context(|| {
            format!(
                "failed to publish {} as {}",
                temporary_tbi.display(),
                destination_tbi.display()
            )
        })?;
        published_tbi = true;

        fs::rename(temporary_vcf, destination_vcf).with_context(|| {
            format!(
                "failed to publish {} as {}",
                temporary_vcf.display(),
                destination_vcf.display()
            )
        })?;
        published_vcf = true;
        Ok(())
    })();

    if let Err(error) = publication {
        let mut rollback_errors = Vec::new();
        if published_vcf && let Err(rollback_error) = fs::remove_file(destination_vcf) {
            rollback_errors.push(format!(
                "failed to remove new {}: {rollback_error}",
                destination_vcf.display()
            ));
        }
        if published_tbi && let Err(rollback_error) = fs::remove_file(destination_tbi) {
            rollback_errors.push(format!(
                "failed to remove new {}: {rollback_error}",
                destination_tbi.display()
            ));
        }
        restore_backup(&backup_vcf, destination_vcf, &mut rollback_errors);
        restore_backup(&backup_tbi, destination_tbi, &mut rollback_errors);

        if rollback_errors.is_empty() {
            return Err(error);
        }
        return Err(error.context(format!("rollback errors: {}", rollback_errors.join("; "))));
    }

    remove_backup(&backup_vcf)?;
    remove_backup(&backup_tbi)?;
    Ok(())
}

struct PairPublicationLock {
    file: File,
    _process_guard: MutexGuard<'static, ()>,
}

impl PairPublicationLock {
    fn acquire(destination: &Path) -> Result<Self> {
        let process_guard = PUBLICATION_MUTEX
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let path = publication_lock_path(destination)?;
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).with_context(|| {
                format!(
                    "failed to create publication lock directory {}",
                    parent.display()
                )
            })?;
        }
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(&path)
            .with_context(|| format!("failed to open publication lock {}", path.display()))?;
        file.lock()
            .with_context(|| format!("failed to lock publication path {}", path.display()))?;
        Ok(Self {
            file,
            _process_guard: process_guard,
        })
    }
}

impl Drop for PairPublicationLock {
    fn drop(&mut self) {
        let _ = self.file.unlock();
    }
}

fn publication_lock_path(destination: &Path) -> Result<PathBuf> {
    // `Path::parent()` returns `Some("")` for a bare relative filename.
    // Canonicalize the working directory in that case, exactly as for an
    // explicit `./result.vcf.gz` destination.
    let parent = destination
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let canonical_parent = parent.canonicalize().with_context(|| {
        format!(
            "failed to resolve publication directory {}",
            parent.display()
        )
    })?;
    let canonical_destination = canonical_parent.join(
        destination
            .file_name()
            .unwrap_or_else(|| std::ffi::OsStr::new("output.vcf.gz")),
    );
    let mut hasher = std::hash::DefaultHasher::new();
    canonical_destination.hash(&mut hasher);
    Ok(std::env::temp_dir()
        .join("hap-rs-publication-locks")
        .join(format!("{:016x}.lock", hasher.finish())))
}

fn validate_regular_destination(path: &Path) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(metadata) if metadata.file_type().is_file() || metadata.file_type().is_symlink() => {
            Ok(())
        }
        Ok(_) => bail!(
            "refusing to replace non-file output destination {}",
            path.display()
        ),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).with_context(|| format!("failed to inspect {}", path.display())),
    }
}

fn move_existing_to_backup(destination: &Path) -> Result<Option<PathBuf>> {
    if !destination.try_exists().with_context(|| {
        format!(
            "failed to determine whether {} exists",
            destination.display()
        )
    })? {
        return Ok(None);
    }

    let (backup, placeholder) = create_temporary_file(destination)?;
    drop(placeholder);
    fs::remove_file(&backup)
        .with_context(|| format!("failed to prepare backup path {}", backup.display()))?;
    fs::rename(destination, &backup).with_context(|| {
        format!(
            "failed to preserve existing {} as {}",
            destination.display(),
            backup.display()
        )
    })?;
    Ok(Some(backup))
}

fn restore_backup(backup: &Option<PathBuf>, destination: &Path, errors: &mut Vec<String>) {
    let Some(backup) = backup else {
        return;
    };
    if let Err(error) = fs::rename(backup, destination) {
        errors.push(format!(
            "failed to restore {} as {}: {error}",
            backup.display(),
            destination.display()
        ));
    }
}

fn remove_backup(backup: &Option<PathBuf>) -> Result<()> {
    if let Some(path) = backup {
        fs::remove_file(path)
            .with_context(|| format!("failed to remove backup {}", path.display()))?;
    }
    Ok(())
}

fn write_vcf_lines<'a, W, I>(writer: &mut W, headers: &[String], lines: I) -> Result<()>
where
    W: Write,
    I: IntoIterator<Item = &'a str>,
{
    for header in headers {
        writeln!(writer, "{header}")?;
    }
    for line in lines {
        writeln!(writer, "{line}")?;
    }
    Ok(())
}

#[derive(Debug)]
struct IndexRecord<'a> {
    line: &'a str,
    chrom: &'a str,
    start: usize,
    end: usize,
}

fn parse_index_records<'a, I>(lines: I, path: &Path) -> Result<Vec<IndexRecord<'a>>>
where
    I: IntoIterator<Item = &'a str>,
{
    let mut records = Vec::new();
    let mut closed_chromosomes = BTreeSet::new();
    let mut current_chromosome: Option<&str> = None;
    let mut previous_start = 0;

    for line in lines {
        if line.contains(['\n', '\r']) {
            bail!(
                "VCF records must contain exactly one line in {}",
                path.display()
            );
        }
        let mut fields = line.split('\t');
        let chrom = fields
            .next()
            .filter(|value| !value.is_empty())
            .ok_or_else(|| anyhow::anyhow!("VCF record has no chromosome in {}", path.display()))?;
        let position = fields
            .next()
            .ok_or_else(|| anyhow::anyhow!("VCF record has no position in {}", path.display()))?
            .parse::<usize>()
            .with_context(|| format!("invalid VCF position in {}", path.display()))?;
        let _id = fields.next();
        let reference = fields.next().ok_or_else(|| {
            anyhow::anyhow!("VCF record has fewer than 4 fields in {}", path.display())
        })?;
        let alternate = fields.next().ok_or_else(|| {
            anyhow::anyhow!("VCF record has fewer than 5 fields in {}", path.display())
        })?;
        let _quality = fields.next().ok_or_else(|| {
            anyhow::anyhow!("VCF record has fewer than 6 fields in {}", path.display())
        })?;
        let _filter = fields.next().ok_or_else(|| {
            anyhow::anyhow!("VCF record has fewer than 7 fields in {}", path.display())
        })?;
        let info = fields.next().ok_or_else(|| {
            anyhow::anyhow!("VCF record has fewer than 8 fields in {}", path.display())
        })?;
        let format = fields.next();
        let samples: Vec<&str> = fields.collect();

        if position == 0 {
            bail!("VCF positions must be >= 1 in {}", path.display());
        }
        if reference.is_empty() {
            bail!("VCF REF alleles must not be empty in {}", path.display());
        }
        let start = position - 1;
        let end = effective_vcf_end(start, reference, alternate, info, format, &samples, path)?;
        if end > TBI_MAX_POSITION {
            bail!(
                "VCF record {chrom}:{position} exceeds the Tabix v1 coordinate limit (2^29); CSI is required"
            );
        }

        match current_chromosome {
            Some(current) if current == chrom => {
                if start < previous_start {
                    bail!(
                        "VCF records are not position-sorted at {chrom}:{position} in {}",
                        path.display()
                    );
                }
            }
            Some(current) => {
                closed_chromosomes.insert(current);
                if closed_chromosomes.contains(chrom) {
                    bail!(
                        "VCF records for chromosome {chrom} are not contiguous in {}",
                        path.display()
                    );
                }
                current_chromosome = Some(chrom);
            }
            None => current_chromosome = Some(chrom),
        }
        previous_start = start;
        records.push(IndexRecord {
            line,
            chrom,
            start,
            end,
        });
    }

    Ok(records)
}

/// Returns the record's 0-based exclusive end using HTSlib's Tabix VCF
/// interval rules. `INFO/END` applies to all records. `INFO/SVLEN` applies
/// to every ALT except the exact symbolic insertion `<INS>`, and
/// `FORMAT/LEN` applies only when ALT contains a gVCF `<*>` or `<NON_REF>`.
fn effective_vcf_end(
    start: usize,
    reference: &str,
    alternate: &str,
    info: &str,
    format: Option<&str>,
    samples: &[&str],
    path: &Path,
) -> Result<usize> {
    let mut end = checked_interval_end(start, reference.len().max(1), path)?;

    if let Some(value) = info_value(info, "END")
        && value != "."
        && let Some(value) = parse_vcf_integer(value, path)?
        && value > start as i128
    {
        let value = usize::try_from(value)
            .map_err(|_| anyhow::anyhow!("VCF coordinate overflow in {}", path.display()))?;
        end = end.max(value);
    }

    let alts: Vec<&str> = alternate.split(',').collect();
    if let Some(values) = info_value(info, "SVLEN") {
        let mut maximum = 0usize;
        for (alt, value) in alts.iter().zip(values.split(',')) {
            if !svlen_contributes_to_reference_span(alt) {
                continue;
            }
            let Some(value) = parse_vcf_integer(value, path)? else {
                continue;
            };
            let magnitude = value.unsigned_abs();
            let magnitude = usize::try_from(magnitude)
                .map_err(|_| anyhow::anyhow!("VCF coordinate overflow in {}", path.display()))?;
            maximum = maximum.max(magnitude);
        }
        if maximum > 0 {
            end = end.max(checked_interval_end(start, maximum, path)?);
        }
    }

    let is_gvcf = alts.iter().any(|alt| matches!(*alt, "<*>" | "<NON_REF>"));
    if is_gvcf
        && let Some(format) = format
        && let Some(length_index) = format.split(':').position(|key| key == "LEN")
    {
        let mut maximum = 0usize;
        for sample in samples {
            let Some(value) = sample.split(':').nth(length_index) else {
                continue;
            };
            let Some(value) = parse_vcf_integer(value, path)? else {
                continue;
            };
            if value <= 0 {
                continue;
            }
            let value = usize::try_from(value)
                .map_err(|_| anyhow::anyhow!("VCF coordinate overflow in {}", path.display()))?;
            maximum = maximum.max(value);
        }
        if maximum > 0 {
            end = end.max(checked_interval_end(start, maximum, path)?);
        }
    }

    Ok(end)
}

fn checked_interval_end(start: usize, length: usize, path: &Path) -> Result<usize> {
    start
        .checked_add(length)
        .ok_or_else(|| anyhow::anyhow!("VCF coordinate overflow in {}", path.display()))
}

fn info_value<'a>(info: &'a str, key: &str) -> Option<&'a str> {
    info.split(';')
        .find_map(|field| field.strip_prefix(key)?.strip_prefix('='))
}

fn parse_vcf_integer(value: &str, path: &Path) -> Result<Option<i128>> {
    let bytes = value.as_bytes();
    let digit_start = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let digit_end = bytes[digit_start..]
        .iter()
        .position(|byte| !byte.is_ascii_digit())
        .map_or(bytes.len(), |index| digit_start + index);
    if digit_end == digit_start {
        return Ok(None);
    }
    value[..digit_end]
        .parse()
        .map(Some)
        .map_err(|_| anyhow::anyhow!("VCF coordinate overflow in {}", path.display()))
}

fn svlen_contributes_to_reference_span(alt: &str) -> bool {
    alt != "<INS>"
}

#[derive(Clone, Copy, Debug)]
struct Chunk {
    start: u64,
    end: u64,
}

impl Chunk {
    fn new(start: u64, end: u64) -> Self {
        Self { start, end }
    }
}

#[derive(Debug, Default)]
struct ReferenceIndex {
    bins: BTreeMap<u32, Vec<Chunk>>,
    intervals: Vec<u64>,
    first_offset: Option<u64>,
    last_offset: u64,
    record_count: u64,
}

impl ReferenceIndex {
    fn push(&mut self, start: usize, end: usize, chunk: Chunk) {
        self.first_offset.get_or_insert(chunk.start);
        self.last_offset = chunk.end;
        self.record_count += 1;

        let chunks = self.bins.entry(reg2bin(start, end)).or_default();
        if let Some(last) = chunks.last_mut()
            && last.end == chunk.start
        {
            last.end = chunk.end;
        } else {
            chunks.push(chunk);
        }

        let first_window = start >> TBI_LINEAR_SHIFT;
        let last_window = (end - 1) >> TBI_LINEAR_SHIFT;
        if self.intervals.len() <= last_window {
            self.intervals.resize(last_window + 1, 0);
        }
        for offset in &mut self.intervals[first_window..=last_window] {
            if *offset == 0 || chunk.start < *offset {
                *offset = chunk.start;
            }
        }
    }

    fn fill_linear_gaps(&mut self) {
        let mut previous = 0;
        for offset in &mut self.intervals {
            if *offset == 0 {
                *offset = previous;
            } else {
                previous = *offset;
            }
        }
    }
}

#[derive(Debug, Default)]
struct TabixIndex<'a> {
    names: Vec<&'a str>,
    references: Vec<ReferenceIndex>,
}

impl<'a> TabixIndex<'a> {
    fn push(&mut self, record: &IndexRecord<'a>, chunk: Chunk) -> Result<()> {
        if self.names.last().copied() != Some(record.chrom) {
            self.names.push(record.chrom);
            self.references.push(ReferenceIndex::default());
        }
        let reference = self
            .references
            .last_mut()
            .ok_or_else(|| anyhow::anyhow!("missing Tabix reference index"))?;
        reference.push(record.start, record.end, chunk);
        Ok(())
    }

    fn finish(mut self) -> Vec<u8> {
        for reference in &mut self.references {
            reference.fill_linear_gaps();
        }

        let names_length: usize = self.names.iter().map(|name| name.len() + 1).sum();
        let mut payload = Vec::new();
        payload.extend_from_slice(b"TBI\x01");
        push_i32(&mut payload, self.references.len());
        push_i32(&mut payload, 2); // VCF preset.
        push_i32(&mut payload, 1); // CHROM column.
        push_i32(&mut payload, 2); // POS column.
        push_i32(&mut payload, 0); // End derived from REF.
        push_i32(&mut payload, b'#' as usize);
        push_i32(&mut payload, 0); // Header lines are meta-prefixed, not skipped by count.
        push_i32(&mut payload, names_length);
        for name in &self.names {
            payload.extend_from_slice(name.as_bytes());
            payload.push(0);
        }

        for reference in self.references {
            push_i32(&mut payload, reference.bins.len() + 1);
            for (bin, chunks) in reference.bins {
                payload.extend_from_slice(&bin.to_le_bytes());
                push_i32(&mut payload, chunks.len());
                for chunk in chunks {
                    payload.extend_from_slice(&chunk.start.to_le_bytes());
                    payload.extend_from_slice(&chunk.end.to_le_bytes());
                }
            }
            payload.extend_from_slice(&TBI_METADATA_BIN.to_le_bytes());
            push_i32(&mut payload, 2);
            payload.extend_from_slice(&reference.first_offset.unwrap_or_default().to_le_bytes());
            payload.extend_from_slice(&reference.last_offset.to_le_bytes());
            payload.extend_from_slice(&reference.record_count.to_le_bytes());
            payload.extend_from_slice(&0u64.to_le_bytes()); // No unplaced records.
            push_i32(&mut payload, reference.intervals.len());
            for offset in reference.intervals {
                payload.extend_from_slice(&offset.to_le_bytes());
            }
        }
        payload.extend_from_slice(&0u64.to_le_bytes()); // n_no_coor
        payload
    }
}

fn push_i32(payload: &mut Vec<u8>, value: usize) {
    payload.extend_from_slice(&(value as i32).to_le_bytes());
}

fn reg2bin(start: usize, end: usize) -> u32 {
    let end = end - 1;
    if start >> 14 == end >> 14 {
        4681 + (start >> 14) as u32
    } else if start >> 17 == end >> 17 {
        585 + (start >> 17) as u32
    } else if start >> 20 == end >> 20 {
        73 + (start >> 20) as u32
    } else if start >> 23 == end >> 23 {
        9 + (start >> 23) as u32
    } else if start >> 26 == end >> 26 {
        1 + (start >> 26) as u32
    } else {
        0
    }
}

fn tabix_path(vcf_path: &Path) -> PathBuf {
    let mut path = vcf_path.as_os_str().to_os_string();
    path.push(".tbi");
    PathBuf::from(path)
}

fn create_temporary_file(destination: &Path) -> Result<(PathBuf, File)> {
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let file_name = destination
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("output");
    for _ in 0..100 {
        let id = TEMP_FILE_ID.fetch_add(1, Ordering::Relaxed);
        let path = parent.join(format!(
            ".{file_name}.hap-rs.{}.{}.tmp",
            std::process::id(),
            id
        ));
        match OpenOptions::new().write(true).create_new(true).open(&path) {
            Ok(file) => return Ok((path, file)),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => {
                return Err(error).with_context(|| {
                    format!(
                        "failed to create temporary file near {}",
                        destination.display()
                    )
                });
            }
        }
    }
    bail!(
        "failed to allocate a temporary file near {}",
        destination.display()
    )
}

pub fn load_variants(
    path: &Path,
    reference_contigs: &BTreeSet<String>,
    pass_only: bool,
    regions: Option<&[BedInterval]>,
    targets: Option<&[BedInterval]>,
    locations: Option<&[LocationFilter]>,
) -> Result<Vec<Variant>> {
    let text = read_text(path)?;
    let mut variants = Vec::new();

    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let record = RawVcfRecord::from_line(line, path)?;
        if record.format.is_none() || record.samples.is_empty() {
            bail!("VCF record has fewer than 10 fields in {}", path.display());
        }
        let chrom = normalize_chrom(&record.chrom, reference_contigs);
        let format_keys: Vec<&str> = record
            .format
            .as_deref()
            .unwrap_or_default()
            .split(':')
            .collect();
        let sample_values: Vec<&str> = record.samples[0].split(':').collect();
        let gt = extract_gt(&format_keys, &sample_values)?.to_string();
        let effective_end = record.effective_end_pos(path)?;

        let variant = Variant {
            key: VariantKey {
                chrom,
                pos: record.pos,
                ref_allele: record.ref_allele.clone(),
                alt_allele: record.alt_allele.clone(),
            },
            qual: canonical_qual(&record.qual).to_string(),
            filter: record.filter,
            gt: canonical_gt(&gt),
        };

        if pass_only && !variant.is_pass() {
            continue;
        }
        if variant.key.alt_allele == "." {
            continue;
        }
        if !matches_interval_filters(
            &variant.key.chrom,
            variant.key.pos,
            effective_end,
            regions,
            targets,
            locations,
        ) {
            continue;
        }

        variants.push(variant);
    }

    Ok(variants)
}

fn extract_gt<'a>(format_keys: &[&str], sample_values: &'a [&str]) -> Result<&'a str> {
    for (index, key) in format_keys.iter().enumerate() {
        if *key == "GT" {
            return sample_values
                .get(index)
                .copied()
                .ok_or_else(|| anyhow::anyhow!("missing GT value"));
        }
    }
    Ok("./.")
}

pub fn canonical_gt(gt: &str) -> String {
    match gt {
        "1" => "1/1".to_string(),
        "0" => "0/0".to_string(),
        "." => "./.".to_string(),
        _ => gt.to_string(),
    }
}

fn canonical_qual(qual: &str) -> &str {
    if qual == "." { "0" } else { qual }
}

pub fn normalize_chrom(chrom: &str, reference_contigs: &BTreeSet<String>) -> String {
    if reference_contigs.contains(chrom) {
        return chrom.to_string();
    }
    if let Some(rest) = chrom.strip_prefix("chr")
        && reference_contigs.contains(rest)
    {
        return rest.to_string();
    }
    let prefixed = format!("chr{chrom}");
    if reference_contigs.contains(&prefixed) {
        return prefixed;
    }
    chrom.to_string()
}

#[derive(Clone, Debug)]
pub struct BedInterval {
    pub chrom: String,
    pub start: usize,
    pub end: usize,
}

impl BedInterval {
    pub fn matches(&self, chrom: &str, pos: usize) -> bool {
        self.chrom == chrom
            && pos.saturating_sub(1) >= self.start
            && pos.saturating_sub(1) < self.end
    }

    pub fn overlaps(&self, chrom: &str, start: usize, end: usize) -> bool {
        self.chrom == chrom && start.saturating_sub(1) < self.end && end > self.start
    }
}

/// Apply legacy bcftools selection semantics shared by the governed wrappers:
/// `-R` selects records whose effective reference span overlaps a BED interval,
/// while `-T` and `-l` select records by their start position. When multiple
/// selectors are supplied they are intersected, matching the legacy pipeline.
pub fn matches_interval_filters(
    chrom: &str,
    pos: usize,
    effective_end: usize,
    regions: Option<&[BedInterval]>,
    targets: Option<&[BedInterval]>,
    locations: Option<&[LocationFilter]>,
) -> bool {
    locations.is_none_or(|filters| filters.iter().any(|filter| filter.matches(chrom, pos)))
        && targets.is_none_or(|intervals| {
            intervals
                .iter()
                .any(|interval| interval.matches(chrom, pos))
        })
        && regions.is_none_or(|intervals| {
            intervals
                .iter()
                .any(|interval| interval.overlaps(chrom, pos, effective_end))
        })
}

pub fn load_bed(path: &Path, reference_contigs: &BTreeSet<String>) -> Result<Vec<BedInterval>> {
    let text = read_text(path).with_context(|| format!("failed to read BED {}", path.display()))?;
    let mut intervals = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            bail!("BED line has fewer than 3 columns in {}", path.display());
        }
        let chrom = normalize_chrom(fields[0], reference_contigs);
        let start = fields[1]
            .parse::<usize>()
            .with_context(|| format!("invalid BED start '{}' in {}", fields[1], path.display()))?;
        let end = fields[2]
            .parse::<usize>()
            .with_context(|| format!("invalid BED end '{}' in {}", fields[2], path.display()))?;
        intervals.push(BedInterval { chrom, start, end });
    }
    Ok(intervals)
}

#[derive(Clone, Debug)]
pub enum LocationFilter {
    Contig(String),
    Range {
        chrom: String,
        start: usize,
        end: usize,
    },
}

impl LocationFilter {
    pub fn matches(&self, chrom: &str, pos: usize) -> bool {
        match self {
            Self::Contig(expected) => expected == chrom,
            Self::Range {
                chrom: expected,
                start,
                end,
            } => expected == chrom && pos >= *start && pos <= *end,
        }
    }
}

pub fn parse_locations(
    text: &str,
    reference_contigs: &BTreeSet<String>,
) -> Result<Vec<LocationFilter>> {
    let mut filters = Vec::new();
    for token in text.split(',').filter(|token| !token.trim().is_empty()) {
        if let Some((chrom_part, position_part)) = token.split_once(':') {
            let (start, end) = position_part
                .split_once('-')
                .unwrap_or((position_part, position_part));
            filters.push(LocationFilter::Range {
                chrom: normalize_chrom(chrom_part, reference_contigs),
                start: start
                    .parse::<usize>()
                    .with_context(|| format!("invalid location start '{start}'"))?,
                end: end
                    .parse::<usize>()
                    .with_context(|| format!("invalid location end '{end}'"))?,
            });
        } else {
            filters.push(LocationFilter::Contig(normalize_chrom(
                token,
                reference_contigs,
            )));
        }
    }
    Ok(filters)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{BTreeMap, BTreeSet};
    use std::io::{BufRead, Cursor};
    use std::process::{Command, Stdio};
    use std::sync::{Arc, Barrier};
    use std::thread;
    use tempfile::tempdir;

    #[test]
    fn single_position_location_matches_only_that_position() -> Result<()> {
        let reference_contigs = BTreeSet::from(["1".to_string()]);
        let locations = parse_locations("1:7", &reference_contigs)?;

        assert_eq!(locations.len(), 1);
        assert!(locations[0].matches("1", 7));
        assert!(!locations[0].matches("1", 6));
        assert!(!locations[0].matches("1", 8));
        assert!(!locations[0].matches("chr1", 7));
        Ok(())
    }

    #[test]
    fn input_format_is_sniffed_independently_of_the_filename_suffix() -> Result<()> {
        let directory = tempdir()?;
        let headers = [
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=chr1,length=100>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let record = RawVcfRecord::from_line(
            "chr1\t7\tformat-sniff\tA\tC\t.\tPASS\t.",
            Path::new("source.vcf"),
        )?;

        let bcf = directory.path().join("source.bcf");
        write_raw_vcf(&bcf, &headers, std::slice::from_ref(&record))?;
        let misnamed_bcf = directory.path().join("misnamed.vcf");
        fs::copy(&bcf, &misnamed_bcf)?;
        assert_eq!(load_raw_vcf(&misnamed_bcf)?.1[0].id, "format-sniff");

        let bgzf = directory.path().join("source.vcf.gz");
        let record_line = record.to_line();
        write_indexed_vcf(&bgzf, &headers, [record_line.as_str()])?;
        let noncanonical_bgzf = directory.path().join("source.bgz");
        fs::copy(&bgzf, &noncanonical_bgzf)?;
        assert_eq!(load_raw_vcf(&noncanonical_bgzf)?.1[0].id, "format-sniff");
        Ok(())
    }

    #[derive(Debug)]
    struct ParsedReferenceIndex {
        bins: BTreeMap<u32, Vec<Chunk>>,
        intervals: Vec<u64>,
    }

    #[derive(Debug)]
    struct ParsedTabixIndex {
        names: Vec<String>,
        references: Vec<ParsedReferenceIndex>,
    }

    #[test]
    fn indexed_vcf_range_queries_match_linear_scan_across_bgzf_blocks() -> Result<()> {
        let directory = tempdir()?;
        let path = directory.path().join("records.vcf.gz");
        let headers = vec![
            "##fileformat=VCFv4.2".to_string(),
            "##contig=<ID=header_only,length=1000>".to_string(),
            "##contig=<ID=chr1,length=1000000>".to_string(),
            "##contig=<ID=chr2,length=1000000>".to_string(),
            "##contig=<ID=chr3,length=1000000>".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];

        let mut records = Vec::new();
        for i in 0usize..900 {
            let position = i * 100 + 1;
            records.push(format!(
                "chr1\t{position}\tchr1_{i}\tA\tC\t.\tPASS\tPAD={:016x}{}",
                i.wrapping_mul(0x9e37_79b9usize),
                "ACGT".repeat(40)
            ));
        }
        records.push(format!(
            "chr2\t16380\tspanning\t{}\tA\t.\tPASS\t.",
            "A".repeat(10)
        ));
        for i in 0..400 {
            let position = 20_001 + i * 100;
            records.push(format!(
                "chr2\t{position}\tchr2_{i}\tG\tT\t.\tPASS\tPAD={}",
                "TGCA".repeat(40)
            ));
        }
        records.push("chr3\t200000\tlast\tC\tG\t.\tPASS\t.".to_string());

        write_indexed_vcf(&path, &headers, records.iter().map(String::as_str))?;

        let index = read_tabix_index(&tabix_path(&path))?;
        assert_eq!(index.names, ["chr1", "chr2", "chr3"]);
        assert!(!index.names.iter().any(|name| name == "header_only"));
        assert!(
            index.references[0]
                .bins
                .iter()
                .filter(|(bin, _)| **bin != TBI_METADATA_BIN)
                .flat_map(|(_, chunks)| chunks)
                .any(|chunk| chunk.start >> 16 != chunk.end >> 16),
            "at least one indexed chunk should span BGZF blocks"
        );
        let metadata = &index.references[0].bins[&TBI_METADATA_BIN];
        assert_eq!(metadata.len(), 2);
        assert_eq!(metadata[1].start, 900);
        assert_eq!(metadata[1].end, 0);

        let queries = [
            ("chr1", 0, 1),
            ("chr1", 44_999, 45_001),
            ("chr1", 89_900, 90_000),
            ("chr1", 90_000, 90_100),
            ("chr2", 16_384, 16_385),
            ("chr2", 25_000, 25_001),
            ("chr3", 199_999, 200_000),
            ("chr3", 200_000, 200_001),
            ("chr3", 0, 1),
            ("missing", 0, 1),
        ];

        for (chrom, start, end) in queries {
            assert_eq!(
                query_with_tabix(&path, &index, chrom, start, end)?,
                query_by_linear_scan(&path, chrom, start, end)?,
                "query mismatch for {chrom}:{start}-{end}"
            );
        }
        Ok(())
    }

    #[test]
    fn invalid_sort_order_preserves_existing_pair_and_creates_no_temporaries() -> Result<()> {
        let directory = tempdir()?;
        let path = directory.path().join("output.vcf.gz");
        let index_path = tabix_path(&path);
        fs::write(&path, b"old-vcf")?;
        fs::write(&index_path, b"old-index")?;

        let records = [
            "chr1\t20\t.\tA\tC\t.\tPASS\t.",
            "chr1\t10\t.\tA\tG\t.\tPASS\t.",
        ];
        let error = write_indexed_vcf(&path, &[], records).unwrap_err();
        assert!(error.to_string().contains("not position-sorted"));
        assert_eq!(fs::read(&path)?, b"old-vcf");
        assert_eq!(fs::read(&index_path)?, b"old-index");
        assert_no_transaction_files(directory.path())?;
        Ok(())
    }

    #[test]
    fn injected_writer_and_index_failures_preserve_existing_pair() -> Result<()> {
        let directory = tempdir()?;
        let path = directory.path().join("output.vcf.gz");
        let index_path = tabix_path(&path);
        let headers = [
            "##fileformat=VCFv4.2".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let records = ["chr1\t1\t.\tA\tC\t.\tPASS\t."];
        for operation in [FailureOperation::Writer, FailureOperation::Index] {
            fs::write(&path, b"old-vcf")?;
            fs::write(&index_path, b"old-index")?;
            crate::output::set_failure_operation(Some(operation));
            let error = write_indexed_vcf(&path, &headers, records)
                .expect_err("injected output operation must fail");
            crate::output::set_failure_operation(None);
            assert!(
                error
                    .to_string()
                    .contains(&index_path.display().to_string())
                    || operation == FailureOperation::Writer
            );
            assert_eq!(fs::read(&path)?, b"old-vcf");
            assert_eq!(fs::read(&index_path)?, b"old-index");
            assert_no_transaction_files(directory.path())?;
        }
        Ok(())
    }

    #[test]
    fn indexed_ranges_honor_end_symbolic_svlen_and_gvcf_sample_len() -> Result<()> {
        let directory = tempdir()?;
        let path = directory.path().join("effective-ends.vcf.gz");
        let headers = [
            "##fileformat=VCFv4.5".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2".to_string(),
        ];
        let records = [
            "chr1\t16000\tfrom_end\tA\tC\t.\tPASS\tEND=16500\tGT\t0/1\t0/0",
            "chr1\t32700\tfrom_svlen\tN\t<DEL>\t.\tPASS\tSVLEN=-100\tGT\t0/1\t0/0",
            "chr1\t49100\tfrom_len\tN\t<NON_REF>\t.\tPASS\t.\tGT:LEN\t0/0:25\t0/0:75",
            // SVLEN for symbolic insertions does not describe reference span.
            "chr1\t65500\tinsertion\tN\t<INS>\t.\tPASS\tSVLEN=100\tGT\t0/1\t0/0",
            // HTSlib excludes only exact <INS>; subtyped insertions contribute.
            "chr1\t82000\tsubtyped_insertion\tN\t<INS:ME>\t.\tPASS\tSVLEN=50\tGT\t0/1\t0/0",
        ];
        write_indexed_vcf(&path, &headers, records)?;
        let index = read_tabix_index(&tabix_path(&path))?;

        let cases = [
            (16_499, 16_500, "from_end"),
            (32_798, 32_799, "from_svlen"),
            (49_173, 49_174, "from_len"),
            (82_048, 82_049, "subtyped_insertion"),
        ];
        for (start, end, expected_id) in cases {
            let matches = query_with_tabix(&path, &index, "chr1", start, end)?;
            assert_eq!(matches.len(), 1, "query chr1:{start}-{end}");
            assert!(matches.iter().any(|line| line.contains(expected_id)));
        }
        assert!(query_with_tabix(&path, &index, "chr1", 65_550, 65_551)?.is_empty());
        Ok(())
    }

    #[test]
    fn variant_loading_distinguishes_region_overlap_from_target_start() -> Result<()> {
        let directory = tempdir()?;
        let input = directory.path().join("spanning.vcf");
        fs::write(
            &input,
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
                "chr1\t3\t.\tC\t<NON_REF>\t.\tPASS\tEND=5\tGT\t0/0\n",
            ),
        )?;
        let contigs = BTreeSet::from(["chr1".to_string()]);
        let boundary = vec![BedInterval {
            chrom: "chr1".to_string(),
            start: 4,
            end: 5,
        }];
        let right_of_span = vec![BedInterval {
            chrom: "chr1".to_string(),
            start: 5,
            end: 6,
        }];

        let by_region = load_variants(&input, &contigs, false, Some(&boundary), None, None)?;
        let by_target = load_variants(&input, &contigs, false, None, Some(&boundary), None)?;
        let by_both = load_variants(
            &input,
            &contigs,
            false,
            Some(&boundary),
            Some(&boundary),
            None,
        )?;
        let after_region =
            load_variants(&input, &contigs, false, Some(&right_of_span), None, None)?;

        assert_eq!(by_region.len(), 1, "-R must use the INFO/END span");
        assert!(by_target.is_empty(), "-T must use the POS coordinate");
        assert!(
            by_both.is_empty(),
            "simultaneous -R/-T selectors must be intersected"
        );
        assert!(
            after_region.is_empty(),
            "the effective end is exclusive in BED coordinates"
        );
        Ok(())
    }

    #[test]
    fn concurrent_same_prefix_writers_publish_one_complete_generation() -> Result<()> {
        let directory = tempdir()?;
        let path = Arc::new(directory.path().join("shared.vcf.gz"));
        let writers = 12;
        let barrier = Arc::new(Barrier::new(writers));
        let mut handles = Vec::new();

        for generation in 0..writers {
            let path = Arc::clone(&path);
            let barrier = Arc::clone(&barrier);
            handles.push(thread::spawn(move || -> Result<()> {
                let headers = [
                    "##fileformat=VCFv4.2".to_string(),
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
                ];
                let chrom = format!("generation_{generation}");
                let records: Vec<String> = (0..200)
                    .map(|record| {
                        format!(
                            "{chrom}\t{}\t{chrom}_{record}\tA\tC\t.\tPASS\tPAD={}",
                            record + 1,
                            "ACGT".repeat(20)
                        )
                    })
                    .collect();
                barrier.wait();
                write_indexed_vcf(&path, &headers, records.iter().map(String::as_str))
            }));
        }

        for handle in handles {
            handle.join().expect("writer thread panicked")?;
        }

        let text = read_text(&path)?;
        let published_chrom = text
            .lines()
            .find(|line| !line.starts_with('#'))
            .and_then(|line| line.split('\t').next())
            .expect("published VCF has a record");
        let index = read_tabix_index(&tabix_path(&path))?;
        assert_eq!(index.names, [published_chrom]);
        assert_eq!(
            text.lines()
                .filter(|line| !line.starts_with('#'))
                .filter(|line| line.starts_with(published_chrom))
                .count(),
            200
        );
        assert_no_transaction_files(directory.path())?;
        Ok(())
    }

    #[test]
    fn concurrent_processes_publish_one_complete_generation() -> Result<()> {
        if let Ok(path) = std::env::var("HAP_RS_PUBLICATION_TEST_PATH") {
            return publication_process_helper(&path);
        }

        let directory = tempdir()?;
        let path = directory.path().join("shared-process.vcf.gz");
        let mut children = Vec::new();
        for generation in 0..8 {
            children.push(
                Command::new(std::env::current_exe()?)
                    .arg("--exact")
                    .arg("vcf::tests::concurrent_processes_publish_one_complete_generation")
                    .env("HAP_RS_PUBLICATION_TEST_PATH", &path)
                    .env("HAP_RS_PUBLICATION_TEST_GENERATION", generation.to_string())
                    .stdout(Stdio::null())
                    .stderr(Stdio::piped())
                    .spawn()?,
            );
        }
        for child in children {
            let output = child.wait_with_output()?;
            assert!(
                output.status.success(),
                "publication helper failed: {}",
                String::from_utf8_lossy(&output.stderr)
            );
        }

        let text = read_text(&path)?;
        let published_chrom = text
            .lines()
            .find(|line| !line.starts_with('#'))
            .and_then(|line| line.split('\t').next())
            .expect("published VCF has a record");
        assert_eq!(
            read_tabix_index(&tabix_path(&path))?.names,
            [published_chrom]
        );
        assert_eq!(
            text.lines()
                .filter(|line| !line.starts_with('#'))
                .filter(|line| line.starts_with(published_chrom))
                .count(),
            500
        );
        assert_no_transaction_files(directory.path())?;
        Ok(())
    }

    fn publication_process_helper(path: &str) -> Result<()> {
        let generation = std::env::var("HAP_RS_PUBLICATION_TEST_GENERATION")?;
        let chrom = format!("process_generation_{generation}");
        let headers = [
            "##fileformat=VCFv4.2".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let records: Vec<String> = (0..500)
            .map(|record| {
                format!(
                    "{chrom}\t{}\t{chrom}_{record}\tA\tC\t.\tPASS\tPAD={}",
                    record + 1,
                    "ACGT".repeat(50)
                )
            })
            .collect();
        write_indexed_vcf(
            Path::new(&path),
            &headers,
            records.iter().map(String::as_str),
        )
    }

    #[test]
    fn replacing_an_existing_pair_removes_transaction_files() -> Result<()> {
        let directory = tempdir()?;
        let path = directory.path().join("output.vcf.gz");
        let index_path = tabix_path(&path);
        fs::write(&path, b"old-vcf")?;
        fs::write(&index_path, b"old-index")?;

        let headers = [
            "##fileformat=VCFv4.2".to_string(),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO".to_string(),
        ];
        let records = ["chr1\t1\t.\tA\tC\t.\tPASS\t."];
        write_indexed_vcf(&path, &headers, records)?;

        assert!(read_text(&path)?.contains("chr1\t1"));
        assert_eq!(read_tabix_index(&index_path)?.names, ["chr1"]);
        assert_no_transaction_files(directory.path())?;
        Ok(())
    }

    #[test]
    fn bare_relative_destination_resolves_publication_lock_from_working_directory() -> Result<()> {
        assert_eq!(
            publication_lock_path(Path::new("result.vcf.gz"))?,
            publication_lock_path(Path::new("./result.vcf.gz"))?
        );
        Ok(())
    }

    #[test]
    fn publication_failure_rolls_back_both_existing_outputs() -> Result<()> {
        let directory = tempdir()?;
        let destination_vcf = directory.path().join("output.vcf.gz");
        let destination_tbi = tabix_path(&destination_vcf);
        let missing_temporary_vcf = directory.path().join("missing.vcf.gz.tmp");
        let temporary_tbi = directory.path().join("new.vcf.gz.tbi.tmp");
        fs::write(&destination_vcf, b"old-vcf")?;
        fs::write(&destination_tbi, b"old-index")?;
        fs::write(&temporary_tbi, b"new-index")?;

        let error = publish_pair(
            &missing_temporary_vcf,
            &destination_vcf,
            &temporary_tbi,
            &destination_tbi,
        )
        .unwrap_err();

        assert!(error.to_string().contains("failed to publish"));
        assert_eq!(fs::read(&destination_vcf)?, b"old-vcf");
        assert_eq!(fs::read(&destination_tbi)?, b"old-index");
        assert_no_transaction_files(directory.path())?;
        Ok(())
    }

    #[test]
    fn rejects_non_contiguous_chromosomes_and_tabix_coordinate_overflow() {
        let directory = tempdir().unwrap();
        let path = directory.path().join("output.vcf.gz");
        let non_contiguous = [
            "chr1\t1\t.\tA\tC\t.\tPASS\t.",
            "chr2\t1\t.\tA\tC\t.\tPASS\t.",
            "chr1\t2\t.\tA\tC\t.\tPASS\t.",
        ];
        assert!(
            write_indexed_vcf(&path, &[], non_contiguous)
                .unwrap_err()
                .to_string()
                .contains("not contiguous")
        );

        let overflow = format!("chr1\t{}\t.\tAA\tC\t.\tPASS\t.", TBI_MAX_POSITION);
        assert!(
            write_indexed_vcf(&path, &[], [overflow.as_str()])
                .unwrap_err()
                .to_string()
                .contains("CSI is required")
        );

        let end_overflow = format!(
            "chr1\t1\t.\tA\tC\t.\tPASS\tEND={}\tGT\t0/1",
            (usize::MAX as u128) + 1
        );
        assert!(
            write_indexed_vcf(&path, &[], [end_overflow.as_str()])
                .unwrap_err()
                .to_string()
                .contains("coordinate overflow")
        );

        let integer_overflow = format!(
            "chr1\t1\t.\tA\tC\t.\tPASS\tEND={}\tGT\t0/1",
            "9".repeat(100)
        );
        assert!(
            write_indexed_vcf(&path, &[], [integer_overflow.as_str()])
                .unwrap_err()
                .to_string()
                .contains("coordinate overflow")
        );

        let span_overflow = format!(
            "chr1\t2\t.\tN\t<DEL>\t.\tPASS\tSVLEN=-{}\tGT\t0/1",
            usize::MAX
        );
        assert!(
            write_indexed_vcf(&path, &[], [span_overflow.as_str()])
                .unwrap_err()
                .to_string()
                .contains("coordinate overflow")
        );

        let len_overflow = format!(
            "chr1\t2\t.\tN\t<NON_REF>\t.\tPASS\t.\tGT:LEN\t0/0:{}",
            usize::MAX
        );
        assert!(
            write_indexed_vcf(&path, &[], [len_overflow.as_str()])
                .unwrap_err()
                .to_string()
                .contains("coordinate overflow")
        );
    }

    fn assert_no_transaction_files(directory: &Path) -> Result<()> {
        let leftovers: Vec<_> = fs::read_dir(directory)?
            .filter_map(|entry| entry.ok())
            .filter(|entry| {
                let name = entry.file_name();
                let name = name.to_string_lossy();
                name.contains(".hap-rs.") && name.ends_with(".tmp")
            })
            .collect();
        assert!(
            leftovers.is_empty(),
            "leftover transaction files: {leftovers:?}"
        );
        Ok(())
    }

    fn read_tabix_index(path: &Path) -> Result<ParsedTabixIndex> {
        let mut reader = bgzf::io::Reader::new(File::open(path)?);
        let mut payload = Vec::new();
        reader.read_to_end(&mut payload)?;
        let mut cursor = Cursor::new(payload);
        let mut magic = [0; 4];
        cursor.read_exact(&mut magic)?;
        assert_eq!(&magic, b"TBI\x01");
        let reference_count = read_i32(&mut cursor)? as usize;
        for _ in 0..6 {
            let _ = read_i32(&mut cursor)?;
        }
        let names_length = read_i32(&mut cursor)? as usize;
        let mut raw_names = vec![0; names_length];
        cursor.read_exact(&mut raw_names)?;
        let names = raw_names
            .split(|byte| *byte == 0)
            .filter(|name| !name.is_empty())
            .map(|name| String::from_utf8(name.to_vec()))
            .collect::<std::result::Result<Vec<_>, _>>()?;

        let mut references = Vec::with_capacity(reference_count);
        for _ in 0..reference_count {
            let bin_count = read_i32(&mut cursor)? as usize;
            let mut bins = BTreeMap::new();
            for _ in 0..bin_count {
                let bin = read_u32(&mut cursor)?;
                let chunk_count = read_i32(&mut cursor)? as usize;
                let mut chunks = Vec::with_capacity(chunk_count);
                for _ in 0..chunk_count {
                    chunks.push(Chunk::new(read_u64(&mut cursor)?, read_u64(&mut cursor)?));
                }
                bins.insert(bin, chunks);
            }
            let interval_count = read_i32(&mut cursor)? as usize;
            let mut intervals = Vec::with_capacity(interval_count);
            for _ in 0..interval_count {
                intervals.push(read_u64(&mut cursor)?);
            }
            references.push(ParsedReferenceIndex { bins, intervals });
        }
        assert_eq!(read_u64(&mut cursor)?, 0);
        assert_eq!(names.len(), reference_count);
        Ok(ParsedTabixIndex { names, references })
    }

    fn query_with_tabix(
        path: &Path,
        index: &ParsedTabixIndex,
        chrom: &str,
        start: usize,
        end: usize,
    ) -> Result<BTreeSet<String>> {
        if start >= end {
            return Ok(BTreeSet::new());
        }
        let Some(reference_id) = index.names.iter().position(|name| name == chrom) else {
            return Ok(BTreeSet::new());
        };
        let reference = &index.references[reference_id];
        let bins = reg2bins(start, end);
        let minimum_offset = reference
            .intervals
            .get(start >> TBI_LINEAR_SHIFT)
            .copied()
            .unwrap_or(0);
        let mut chunks: Vec<Chunk> = bins
            .iter()
            .filter_map(|bin| reference.bins.get(bin))
            .flatten()
            .copied()
            .filter(|chunk| chunk.end > minimum_offset)
            .collect();
        chunks.sort_by_key(|chunk| chunk.start);
        let mut merged: Vec<Chunk> = Vec::new();
        for chunk in chunks {
            if let Some(previous) = merged.last_mut()
                && chunk.start <= previous.end
            {
                previous.end = previous.end.max(chunk.end);
            } else {
                merged.push(chunk);
            }
        }

        let mut matches = BTreeSet::new();
        let mut reader = bgzf::io::Reader::new(File::open(path)?);
        for chunk in merged {
            reader.seek(bgzf::VirtualPosition::from(chunk.start))?;
            while u64::from(reader.virtual_position()) < chunk.end {
                let mut line = String::new();
                if reader.read_line(&mut line)? == 0 {
                    break;
                }
                let line = line.trim_end_matches(['\n', '\r']);
                if record_overlaps(line, chrom, start, end)? {
                    matches.insert(line.to_string());
                }
            }
        }
        Ok(matches)
    }

    fn query_by_linear_scan(
        path: &Path,
        chrom: &str,
        start: usize,
        end: usize,
    ) -> Result<BTreeSet<String>> {
        let text = read_text(path)?;
        let mut matches = BTreeSet::new();
        for line in text.lines().filter(|line| !line.starts_with('#')) {
            if record_overlaps(line, chrom, start, end)? {
                matches.insert(line.to_string());
            }
        }
        Ok(matches)
    }

    fn record_overlaps(line: &str, chrom: &str, start: usize, end: usize) -> Result<bool> {
        let fields: Vec<&str> = line.split('\t').collect();
        let record_chrom = fields.first().copied().unwrap_or_default();
        let position = fields
            .get(1)
            .copied()
            .unwrap_or_default()
            .parse::<usize>()?;
        let reference = fields.get(3).copied().unwrap_or_default();
        let alternate = fields.get(4).copied().unwrap_or_default();
        let info = fields.get(7).copied().unwrap_or_default();
        let format = fields.get(8).copied();
        let samples = fields.get(9..).unwrap_or_default();
        let record_start = position.saturating_sub(1);
        let record_end = effective_vcf_end(
            record_start,
            reference,
            alternate,
            info,
            format,
            samples,
            Path::new("test.vcf"),
        )?;
        Ok(record_chrom == chrom && record_start < end && record_end > start)
    }

    fn reg2bins(start: usize, end: usize) -> BTreeSet<u32> {
        let end = end - 1;
        let mut bins = BTreeSet::from([0]);
        for bin in (1 + (start >> 26))..=(1 + (end >> 26)) {
            bins.insert(bin as u32);
        }
        for bin in (9 + (start >> 23))..=(9 + (end >> 23)) {
            bins.insert(bin as u32);
        }
        for bin in (73 + (start >> 20))..=(73 + (end >> 20)) {
            bins.insert(bin as u32);
        }
        for bin in (585 + (start >> 17))..=(585 + (end >> 17)) {
            bins.insert(bin as u32);
        }
        for bin in (4681 + (start >> 14))..=(4681 + (end >> 14)) {
            bins.insert(bin as u32);
        }
        bins
    }

    fn read_i32(reader: &mut impl Read) -> Result<i32> {
        let mut bytes = [0; 4];
        reader.read_exact(&mut bytes)?;
        Ok(i32::from_le_bytes(bytes))
    }

    fn read_u32(reader: &mut impl Read) -> Result<u32> {
        let mut bytes = [0; 4];
        reader.read_exact(&mut bytes)?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn read_u64(reader: &mut impl Read) -> Result<u64> {
        let mut bytes = [0; 8];
        reader.read_exact(&mut bytes)?;
        Ok(u64::from_le_bytes(bytes))
    }
}
