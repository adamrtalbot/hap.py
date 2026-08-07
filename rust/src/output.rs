//! Transactional publication for complete command output generations.

use anyhow::{Context, Result, bail};
use std::collections::{BTreeSet, HashMap, HashSet};
use std::ffi::{OsStr, OsString};
use std::fs::{self, File, OpenOptions};
#[cfg(unix)]
use std::io::{Seek, SeekFrom};
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, MutexGuard};
use std::time::{SystemTime, UNIX_EPOCH};

static TRANSACTION_ID: AtomicU64 = AtomicU64::new(0);
static PUBLICATION_MUTEX: Mutex<()> = Mutex::new(());
#[cfg(test)]
thread_local! {
    static FAIL_PUBLICATION_AFTER: std::cell::Cell<isize> = const { std::cell::Cell::new(-1) };
    static FAIL_BACKUP_CLEANUP: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static FAIL_OPERATION: std::cell::Cell<Option<FailureOperation>> = const { std::cell::Cell::new(None) };
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum FailureOperation {
    Writer,
    Encoder,
    Index,
}

pub(crate) fn fail_operation(operation: FailureOperation, destination: &Path) -> Result<()> {
    #[cfg(test)]
    if FAIL_OPERATION.get() == Some(operation) {
        bail!(
            "injected {operation:?} failure for destination {}",
            destination.display()
        );
    }
    #[cfg(not(test))]
    let _ = operation;
    #[cfg(not(test))]
    let _ = destination;
    Ok(())
}

#[cfg(test)]
pub(crate) fn set_failure_operation(operation: Option<FailureOperation>) {
    FAIL_OPERATION.set(operation);
}

#[derive(Clone, Debug)]
struct ArtifactId(String);

impl ArtifactId {
    fn new(value: &str) -> Result<Self> {
        validate_artifact_component("output artifact", value)?;
        Ok(Self(value.to_string()))
    }
}

impl AsRef<str> for ArtifactId {
    fn as_ref(&self) -> &str {
        &self.0
    }
}

#[derive(Debug)]
struct Target {
    destination: PathBuf,
    staged: PathBuf,
    required: bool,
}

/// A completely validated set of owned artifacts. Writers only see adjacent,
/// unique staging names; publication happens after every writer has finished.
pub(crate) struct OutputTransaction {
    targets: Vec<Target>,
    inputs: Vec<PathBuf>,
    staged_prefix: Option<PathBuf>,
    committed: bool,
}

impl OutputTransaction {
    pub(crate) fn files<I, O>(inputs: I, outputs: O) -> Result<Self>
    where
        I: IntoIterator,
        I::Item: AsRef<Path>,
        O: IntoIterator,
        O::Item: AsRef<Path>,
    {
        Self::new(
            inputs,
            outputs
                .into_iter()
                .map(|p| (p.as_ref().to_path_buf(), true)),
            None,
        )
    }

    /// Declares every artifact owned by a report prefix. Undeclared siblings
    /// are never inspected or removed. Family artifacts are optional because
    /// command flags and input content determine which ROC lanes are emitted.
    pub(crate) fn family<I, S>(inputs: I, prefix: &Path, suffixes: S) -> Result<Self>
    where
        I: IntoIterator,
        I::Item: AsRef<Path>,
        S: IntoIterator,
        S::Item: AsRef<str>,
    {
        let inputs = expand_input_paths(
            inputs
                .into_iter()
                .map(|p| p.as_ref().to_path_buf())
                .collect(),
        );
        let suffixes = suffixes
            .into_iter()
            .map(|suffix| ArtifactId::new(suffix.as_ref()))
            .collect::<Result<Vec<_>>>()?;
        let prospective = suffixes
            .iter()
            .map(|suffix| append_suffix(prefix, format!(".{}", suffix.as_ref()).as_ref()))
            .collect::<Vec<_>>();
        validate_plan(&inputs, prospective.iter())?;
        ensure_output_parents(&prospective)?;
        let prefix = anchored_path(prefix)?;
        let staged_prefix = staged_file_path(&prefix, &transaction_token());
        let targets = suffixes
            .into_iter()
            .map(|suffix| {
                let suffix = format!(".{}", suffix.as_ref());
                Target {
                    destination: append_suffix(&prefix, suffix.as_ref()),
                    staged: append_suffix(&staged_prefix, suffix.as_ref()),
                    required: false,
                }
            })
            .collect::<Vec<_>>();
        if targets.is_empty() {
            bail!("output plan is empty");
        }
        validate_plan(&inputs, targets.iter().map(|target| &target.destination))?;
        Ok(Self {
            targets,
            inputs,
            staged_prefix: Some(staged_prefix),
            committed: false,
        })
    }

    fn new<I, O>(inputs: I, outputs: O, staged_prefix: Option<PathBuf>) -> Result<Self>
    where
        I: IntoIterator,
        I::Item: AsRef<Path>,
        O: IntoIterator<Item = (PathBuf, bool)>,
    {
        let inputs = expand_input_paths(
            inputs
                .into_iter()
                .map(|p| p.as_ref().to_path_buf())
                .collect(),
        );
        let outputs = outputs.into_iter().collect::<Vec<_>>();
        if outputs.is_empty() {
            bail!("output plan is empty");
        }
        validate_plan(&inputs, outputs.iter().map(|(p, _)| p))?;
        let token = transaction_token();
        let targets = outputs
            .into_iter()
            .map(|(path, required)| {
                let destination = anchored_path(&path)?;
                let staged = staged_file_path(&destination, &token);
                Ok(Target {
                    destination,
                    staged,
                    required,
                })
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            targets,
            inputs,
            staged_prefix,
            committed: false,
        })
    }

    pub(crate) fn staged_file(&self, destination: &Path) -> Result<&Path> {
        let anchored = anchored_path(destination)?;
        self.targets
            .iter()
            .find(|t| t.destination == anchored)
            .map(|t| t.staged.as_path())
            .with_context(|| format!("{} is not in the output plan", destination.display()))
    }

    pub(crate) fn staged_prefix(&self) -> Result<&Path> {
        self.staged_prefix
            .as_deref()
            .context("output plan has no report family")
    }

    pub(crate) fn with_files<O>(mut self, outputs: O) -> Result<Self>
    where
        O: IntoIterator,
        O::Item: AsRef<Path>,
    {
        let outputs = outputs
            .into_iter()
            .map(|output| output.as_ref().to_path_buf())
            .collect::<Vec<_>>();
        let mut planned = self
            .targets
            .iter()
            .map(|target| target.destination.clone())
            .collect::<Vec<_>>();
        planned.extend(outputs.iter().cloned());
        validate_plan(&self.inputs, planned.iter())?;
        ensure_output_parents(&outputs)?;
        let token = transaction_token();
        for output in outputs {
            let destination = anchored_path(&output)?;
            self.targets.push(Target {
                staged: staged_file_path(&destination, &token),
                destination,
                required: true,
            });
        }
        Ok(self)
    }

    pub(crate) fn commit(mut self) -> Result<()> {
        let mut present = Vec::new();
        for target in &self.targets {
            if target.staged.is_file() {
                present.push((target.staged.clone(), target.destination.clone()));
            } else if target.required {
                bail!(
                    "writer did not finish destination {} (staging file {} is missing)",
                    target.destination.display(),
                    target.staged.display()
                );
            }
        }
        if present.is_empty() {
            bail!("output transaction produced no artifacts");
        }
        validate_plan(&self.inputs, self.targets.iter().map(|t| &t.destination))?;
        let destinations = self
            .targets
            .iter()
            .map(|t| t.destination.clone())
            .collect::<Vec<_>>();
        let _locks = PublicationLocks::acquire(&destinations)?;
        validate_plan(&self.inputs, self.targets.iter().map(|t| &t.destination))?;
        publish_generation(&present, &destinations)?;
        self.committed = true;
        Ok(())
    }
}

pub(crate) fn benchmark_artifacts(labels: impl IntoIterator<Item = String>) -> Vec<String> {
    let mut suffixes = vec![
        "summary.csv",
        "extended.csv",
        "vcf.gz",
        "vcf.gz.tbi",
        "bcf",
        "bcf.csi",
        "runinfo.json",
        "metrics.json.gz",
        "roc.all.csv.gz",
        "roc.tsv",
    ]
    .into_iter()
    .map(str::to_string)
    .collect::<Vec<_>>();
    let mut subsets = BTreeSet::from(["Locations".to_string()]);
    subsets.extend(labels);
    for subset in subsets {
        for variant in ["SNP", "INDEL"] {
            suffixes.push(format!("roc.{subset}.{variant}.csv.gz"));
            suffixes.push(format!("roc.{subset}.{variant}.PASS.csv.gz"));
            suffixes.push(format!("roc.{subset}.{variant}.SEL.csv.gz"));
        }
    }
    suffixes
}

fn validate_artifact_component(kind: &str, value: &str) -> Result<()> {
    let path = Path::new(value);
    if value.is_empty()
        || value == "."
        || value == ".."
        || value.contains('/')
        || value.contains('\\')
        || path.is_absolute()
        || path.components().count() != 1
    {
        bail!("invalid {kind} '{value}': expected a safe single path component");
    }
    Ok(())
}

pub(crate) fn stratification_inputs(
    tsv: Option<&str>,
    specs: &[String],
) -> Result<(Vec<PathBuf>, Vec<String>)> {
    let mut entries = Vec::<(String, PathBuf)>::new();
    let mut inputs = Vec::new();
    if let Some(tsv) = tsv {
        let tsv_path = PathBuf::from(tsv);
        let text = fs::read_to_string(&tsv_path)
            .with_context(|| format!("failed to read stratification TSV {}", tsv_path.display()))?;
        inputs.push(tsv_path.clone());
        for (line_index, line) in text.lines().enumerate() {
            if line.trim().is_empty() || line.starts_with('#') {
                continue;
            }
            let (name, raw_path) = line.split_once('\t').with_context(|| {
                format!(
                    "stratification TSV line {} has no region file",
                    line_index + 1
                )
            })?;
            let raw_path = PathBuf::from(raw_path.trim());
            let path = if raw_path.is_absolute() || raw_path.exists() {
                raw_path
            } else {
                tsv_path
                    .parent()
                    .unwrap_or_else(|| Path::new("."))
                    .join(raw_path)
            };
            entries.push((name.trim().to_string(), path));
        }
    }
    for spec in specs {
        let (name, path) = spec.split_once(':').with_context(|| {
            format!("invalid --stratification-region '{spec}'; expected NAME:BED")
        })?;
        entries.push((name.trim().to_string(), PathBuf::from(path.trim())));
    }
    let mut labels = BTreeSet::new();
    for (raw_name, path) in entries {
        if raw_name.is_empty() || path.as_os_str().is_empty() {
            bail!("invalid empty stratification region");
        }
        let fixed = raw_name.starts_with('=');
        let name = raw_name.trim_start_matches('=').to_string();
        validate_artifact_component("stratification name", &name)?;
        labels.insert(name.clone());
        let text = crate::adapters::vcf::read_text(&path).with_context(|| {
            format!("failed to preflight stratification BED {}", path.display())
        })?;
        if !fixed {
            for line in text
                .lines()
                .filter(|line| !line.is_empty() && !line.starts_with('#'))
            {
                if let Some(label) = line.split('\t').nth(3).filter(|label| !label.is_empty()) {
                    validate_artifact_component("stratification BED label", label)?;
                    labels.insert(format!("{name}_{label}"));
                }
            }
        }
        inputs.push(path);
    }
    Ok((inputs, labels.into_iter().collect()))
}

fn expand_input_paths(inputs: Vec<PathBuf>) -> Vec<PathBuf> {
    let mut expanded = Vec::new();
    for input in inputs {
        let name = input
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .to_ascii_lowercase();
        expanded.push(input.clone());
        let suffixes: &[&str] = if name.ends_with(".vcf.gz")
            || name.ends_with(".vcf.bgz")
            || name.ends_with(".vcf.bgzf")
        {
            &["tbi", "csi"]
        } else if name.ends_with(".bcf") {
            &["csi"]
        } else if name.ends_with(".fa") || name.ends_with(".fasta") || name.ends_with(".fna") {
            &["fai"]
        } else if name.ends_with(".bam") {
            &["bai"]
        } else {
            &[]
        };
        for suffix in suffixes {
            expanded.push(PathBuf::from(format!("{}.{suffix}", input.display())));
        }
    }
    expanded
}

impl Drop for OutputTransaction {
    fn drop(&mut self) {
        if !self.committed {
            for target in &self.targets {
                let _ = fs::remove_file(&target.staged);
            }
        }
    }
}

fn validate_plan<'a>(
    inputs: &[PathBuf],
    outputs: impl IntoIterator<Item = &'a PathBuf>,
) -> Result<()> {
    let mut input_keys = HashMap::<PathBuf, &Path>::new();
    let mut input_ids = HashMap::<ExistingFileId, &Path>::new();
    for input in inputs {
        input_keys.entry(resolved_path(input)?).or_insert(input);
        if let Some(identity) = existing_file_id(input) {
            input_ids.entry(identity).or_insert(input);
        }
    }
    let mut seen_keys = HashMap::<PathBuf, &Path>::new();
    let mut seen_ids = HashMap::<ExistingFileId, &Path>::new();
    for output in outputs {
        validate_destination_type(output)?;
        let key = resolved_path(output)?;
        let identity = existing_file_id(output);
        if let Some(first) = seen_keys
            .get(&key)
            .copied()
            .or_else(|| identity.and_then(|identity| seen_ids.get(&identity).copied()))
        {
            bail!(
                "output destinations must use distinct paths: {} and {} refer to the same file",
                first.display(),
                output.display()
            );
        }
        if let Some(input) = input_keys
            .get(&key)
            .copied()
            .or_else(|| identity.and_then(|identity| input_ids.get(&identity).copied()))
        {
            bail!(
                "output {} would overwrite input {}",
                output.display(),
                input.display()
            );
        }
        seen_keys.insert(key, output);
        if let Some(identity) = identity {
            seen_ids.insert(identity, output);
        }
    }
    Ok(())
}

fn validate_destination_type(path: &Path) -> Result<()> {
    match fs::symlink_metadata(path) {
        Ok(m) if m.file_type().is_file() || m.file_type().is_symlink() => Ok(()),
        Ok(_) => bail!(
            "refusing to replace non-file output destination {}",
            path.display()
        ),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e)
            .with_context(|| format!("failed to inspect output destination {}", path.display())),
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct ExistingFileId(u64, u64);

#[cfg(unix)]
fn existing_file_id(path: &Path) -> Option<ExistingFileId> {
    use std::os::unix::fs::MetadataExt;
    fs::metadata(path)
        .ok()
        .map(|metadata| ExistingFileId(metadata.dev(), metadata.ino()))
}
#[cfg(not(unix))]
fn existing_file_id(_: &Path) -> Option<ExistingFileId> {
    None
}

fn anchored_path(path: &Path) -> Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()?.join(path)
    };
    let normalized = lexical_normalize(&absolute);
    let parent = normalized
        .parent()
        .with_context(|| format!("path has no parent: {}", path.display()))?;
    let name = normalized
        .file_name()
        .with_context(|| format!("path has no file name: {}", path.display()))?;
    let mut ancestor = parent;
    let mut tail = Vec::new();
    while !ancestor.exists() {
        tail.push(
            ancestor
                .file_name()
                .context("output path has no existing ancestor")?
                .to_os_string(),
        );
        ancestor = ancestor
            .parent()
            .context("output path has no existing ancestor")?;
    }
    let mut parent = fs::canonicalize(ancestor).with_context(|| {
        format!(
            "failed to canonicalize output ancestor {}",
            ancestor.display()
        )
    })?;
    for component in tail.iter().rev() {
        parent.push(component);
    }
    Ok(parent.join(name))
}

fn ensure_output_parents(outputs: &[PathBuf]) -> Result<()> {
    let parents = outputs
        .iter()
        .filter_map(|path| path.parent())
        .collect::<BTreeSet<_>>();
    for parent in parents {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).with_context(|| {
                format!("failed to create output directory {}", parent.display())
            })?;
        }
    }
    Ok(())
}

fn resolved_path(path: &Path) -> Result<PathBuf> {
    let mut candidate = anchored_path(path)?;
    let mut visited = HashSet::new();
    for _ in 0..64 {
        if !visited.insert(candidate.clone()) {
            bail!("symlink cycle while resolving {}", path.display());
        }
        match fs::symlink_metadata(&candidate) {
            Ok(m) if m.file_type().is_symlink() => {
                let target = fs::read_link(&candidate).with_context(|| {
                    format!("failed to resolve symlink {}", candidate.display())
                })?;
                candidate = if target.is_absolute() {
                    lexical_normalize(&target)
                } else {
                    lexical_normalize(
                        &candidate
                            .parent()
                            .context("symlink has no parent")?
                            .join(target),
                    )
                };
            }
            Ok(_) => {
                return fs::canonicalize(&candidate)
                    .with_context(|| format!("failed to canonicalize {}", candidate.display()));
            }
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(candidate),
            Err(e) => {
                return Err(e)
                    .with_context(|| format!("failed to inspect {}", candidate.display()));
            }
        }
    }
    bail!("too many symlinks while resolving {}", path.display())
}

fn lexical_normalize(path: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for c in path.components() {
        match c {
            Component::CurDir => {}
            Component::ParentDir => {
                out.pop();
            }
            other => out.push(other.as_os_str()),
        }
    }
    out
}
fn transaction_token() -> String {
    let n = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    format!(
        "hap-rs.{}.{}.{}",
        std::process::id(),
        n,
        TRANSACTION_ID.fetch_add(1, Ordering::Relaxed)
    )
}
fn staged_file_path(destination: &Path, token: &str) -> PathBuf {
    let parent = destination.parent().unwrap_or_else(|| Path::new("."));
    let name = destination
        .file_name()
        .unwrap_or_else(|| OsStr::new("output"));
    let mut staged = OsString::from(".");
    staged.push(token);
    staged.push(".");
    staged.push(name);
    parent.join(staged)
}
fn append_suffix(prefix: &Path, suffix: &OsStr) -> PathBuf {
    let mut value = prefix.as_os_str().to_os_string();
    value.push(suffix);
    PathBuf::from(value)
}

struct PublicationLocks {
    file: File,
    _process_guard: MutexGuard<'static, ()>,
}
impl PublicationLocks {
    fn acquire(_destinations: &[PathBuf]) -> Result<Self> {
        let process_guard = PUBLICATION_MUTEX
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let directory = std::env::temp_dir().join("hap-rs-publication-locks");
        fs::create_dir_all(&directory).with_context(|| {
            format!(
                "failed to create publication lock directory {}",
                directory.display()
            )
        })?;
        let path = directory.join("global.lock");
        let file = OpenOptions::new()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(&path)
            .with_context(|| format!("failed to open publication lock {}", path.display()))?;
        file.lock().context("failed to lock output publication")?;
        Ok(Self {
            file,
            _process_guard: process_guard,
        })
    }
}
impl Drop for PublicationLocks {
    fn drop(&mut self) {
        let _ = self.file.unlock();
    }
}

fn publish_generation(publications: &[(PathBuf, PathBuf)], owned: &[PathBuf]) -> Result<()> {
    let token = transaction_token();
    let mut backups = Vec::new();
    for destination in owned {
        let backup = match fs::symlink_metadata(destination) {
            Ok(m) if m.file_type().is_file() || m.file_type().is_symlink() => {
                let backup = staged_file_path(destination, &format!("{token}.backup"));
                fs::rename(destination, &backup)
                    .with_context(|| {
                        format!(
                            "failed to preserve existing output {}",
                            destination.display()
                        )
                    })
                    .map(|()| Some((backup, destination.clone())))
            }
            Ok(_) => Err(anyhow::anyhow!(
                "refusing to replace non-file output destination {}",
                destination.display()
            )),
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(e).with_context(|| {
                format!(
                    "failed to inspect output destination {}",
                    destination.display()
                )
            }),
        };
        match backup {
            Ok(Some(backup)) => backups.push(backup),
            Ok(None) => {}
            Err(error) => {
                let mut rollback = Vec::new();
                restore_backups(&backups, &mut rollback);
                return if rollback.is_empty() {
                    Err(error)
                } else {
                    Err(error.context(format!("backup rollback errors: {}", rollback.join("; "))))
                };
            }
        }
    }
    let mut published: Vec<(PathBuf, PathBuf)> = Vec::new();
    for (index, (staged, destination)) in publications.iter().enumerate() {
        let result = if publication_failure_injected(index) {
            Err(anyhow::anyhow!(
                "injected rename failure for destination {}",
                destination.display()
            ))
        } else {
            fs::rename(staged, destination).with_context(|| {
                format!(
                    "failed to publish artifact {} as destination {}",
                    staged.display(),
                    destination.display()
                )
            })
        };
        if let Err(error) = result {
            let mut rollback = Vec::new();
            for (path, staged) in published.iter().rev() {
                if let Err(e) = fs::rename(path, staged) {
                    rollback.push(format!("failed to unpublish {}: {e}", path.display()));
                }
            }
            restore_backups(&backups, &mut rollback);
            return if rollback.is_empty() {
                Err(error)
            } else {
                Err(error.context(format!("rollback errors: {}", rollback.join("; "))))
            };
        }
        published.push((destination.clone(), staged.clone()));
    }
    cleanup_backups_or_rollback(&backups, &published)
}

#[cfg(test)]
fn publication_failure_injected(index: usize) -> bool {
    FAIL_PUBLICATION_AFTER.get() == index as isize
}
#[cfg(not(test))]
fn publication_failure_injected(_: usize) -> bool {
    false
}
fn cleanup_backups_or_rollback(
    backups: &[(PathBuf, PathBuf)],
    published: &[(PathBuf, PathBuf)],
) -> Result<()> {
    let snapshots = match backups
        .iter()
        .map(|(backup, _)| BackupSnapshot::capture(backup))
        .collect::<Result<Vec<_>>>()
    {
        Ok(snapshots) => snapshots,
        Err(error) => return rollback_after_cleanup_failure(backups, published, None, error),
    };
    for (index, (backup, _)) in backups.iter().enumerate() {
        if let Err(error) = fs::remove_file(backup) {
            return rollback_after_cleanup_failure(
                backups,
                published,
                Some(&snapshots),
                error.context(format!(
                    "failed to remove transaction backup {}",
                    backup.display()
                )),
            );
        }
        debug_assert!(!backup.exists(), "backup {index} was not removed");
        #[cfg(test)]
        if index == 0 && FAIL_BACKUP_CLEANUP.get() {
            return rollback_after_cleanup_failure(
                backups,
                published,
                Some(&snapshots),
                anyhow::anyhow!("injected persistent backup cleanup failure"),
            );
        }
    }
    Ok(())
}

fn rollback_after_cleanup_failure(
    backups: &[(PathBuf, PathBuf)],
    published: &[(PathBuf, PathBuf)],
    snapshots: Option<&[BackupSnapshot]>,
    error: anyhow::Error,
) -> Result<()> {
    let mut rollback_errors = Vec::new();
    for (destination, staged) in published.iter().rev() {
        if let Err(rollback) = fs::rename(destination, staged) {
            rollback_errors.push(format!(
                "failed to unpublish {}: {rollback}",
                destination.display()
            ));
        }
    }
    for (index, (backup, destination)) in backups.iter().enumerate().rev() {
        let restored = if backup.exists() {
            fs::rename(backup, destination).map_err(anyhow::Error::from)
        } else if let Some(snapshot) = snapshots.and_then(|values| values.get(index)) {
            snapshot.restore(destination)
        } else {
            Err(anyhow::anyhow!("transaction backup disappeared"))
        };
        if let Err(rollback) = restored {
            rollback_errors.push(format!(
                "failed to restore {}: {rollback:#}",
                destination.display()
            ));
        }
    }
    if rollback_errors.is_empty() {
        Err(error)
    } else {
        Err(error.context(format!(
            "cleanup rollback errors: {}",
            rollback_errors.join("; ")
        )))
    }
}

enum BackupSnapshot {
    File {
        #[cfg(unix)]
        source: File,
        #[cfg(not(unix))]
        contents: Vec<u8>,
        permissions: fs::Permissions,
    },
    Symlink(PathBuf),
}

impl BackupSnapshot {
    fn capture(path: &Path) -> Result<Self> {
        let metadata = fs::symlink_metadata(path)
            .with_context(|| format!("failed to inspect transaction backup {}", path.display()))?;
        if metadata.file_type().is_symlink() {
            Ok(Self::Symlink(fs::read_link(path).with_context(|| {
                format!(
                    "failed to read transaction backup symlink {}",
                    path.display()
                )
            })?))
        } else {
            Ok(Self::File {
                #[cfg(unix)]
                source: File::open(path).with_context(|| {
                    format!("failed to open transaction backup {}", path.display())
                })?,
                #[cfg(not(unix))]
                contents: fs::read(path).with_context(|| {
                    format!("failed to snapshot transaction backup {}", path.display())
                })?,
                permissions: metadata.permissions(),
            })
        }
    }

    fn restore(&self, destination: &Path) -> Result<()> {
        match self {
            Self::File {
                #[cfg(unix)]
                source,
                #[cfg(not(unix))]
                contents,
                permissions,
            } => {
                #[cfg(unix)]
                {
                    let mut source = source.try_clone().with_context(|| {
                        format!("failed to clone backup for {}", destination.display())
                    })?;
                    source.seek(SeekFrom::Start(0))?;
                    let mut output = File::create(destination).with_context(|| {
                        format!("failed to recreate output {}", destination.display())
                    })?;
                    std::io::copy(&mut source, &mut output).with_context(|| {
                        format!("failed to restore output {}", destination.display())
                    })?;
                    output.sync_all().with_context(|| {
                        format!("failed to sync restored output {}", destination.display())
                    })?;
                }
                #[cfg(not(unix))]
                fs::write(destination, contents).with_context(|| {
                    format!("failed to restore output {}", destination.display())
                })?;
                fs::set_permissions(destination, permissions.clone()).with_context(|| {
                    format!(
                        "failed to restore permissions for {}",
                        destination.display()
                    )
                })
            }
            Self::Symlink(target) => restore_symlink(target, destination),
        }
    }
}

#[cfg(unix)]
fn restore_symlink(target: &Path, destination: &Path) -> Result<()> {
    std::os::unix::fs::symlink(target, destination)
        .with_context(|| format!("failed to restore output symlink {}", destination.display()))
}

#[cfg(windows)]
fn restore_symlink(target: &Path, destination: &Path) -> Result<()> {
    std::os::windows::fs::symlink_file(target, destination)
        .with_context(|| format!("failed to restore output symlink {}", destination.display()))
}
fn restore_backups(backups: &[(PathBuf, PathBuf)], errors: &mut Vec<String>) {
    for (backup, destination) in backups.iter().rev() {
        if let Err(e) = fs::rename(backup, destination) {
            errors.push(format!("failed to restore {}: {e}", destination.display()));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::{Arc, Barrier};
    use std::thread;
    use tempfile::tempdir;

    #[test]
    fn artifact_components_reject_traversal_and_separators() {
        for invalid in ["", ".", "..", "../escape", "a/b", "a\\b", "/absolute"] {
            assert!(validate_artifact_component("test artifact", invalid).is_err());
        }
        for valid in ["summary.csv", "roc.Locations.SNP.csv.gz", "region_name"] {
            assert!(validate_artifact_component("test artifact", valid).is_ok());
        }
    }

    #[test]
    fn writer_and_index_injection_include_destination_context() {
        for operation in [FailureOperation::Writer, FailureOperation::Index] {
            set_failure_operation(Some(operation));
            let destination = Path::new("logical-output.vcf.gz.tbi");
            let error = fail_operation(operation, destination).expect_err("failure must inject");
            assert!(error.to_string().contains("logical-output.vcf.gz.tbi"));
        }
        set_failure_operation(None);
    }

    #[test]
    fn explicit_family_preserves_unowned_siblings() -> Result<()> {
        let dir = tempdir()?;
        let prefix = dir.path().join("run");
        let unrelated = dir.path().join("run.notes");
        fs::write(&unrelated, "mine")?;
        let tx = OutputTransaction::family(Vec::<PathBuf>::new(), &prefix, ["summary.csv"])?;
        fs::write(
            append_suffix(tx.staged_prefix()?, OsStr::new(".summary.csv")),
            "new",
        )?;
        tx.commit()?;
        assert_eq!(fs::read_to_string(unrelated)?, "mine");
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn rejects_input_symlink_and_hardlink_aliases_before_staging() -> Result<()> {
        let dir = tempdir()?;
        let input = dir.path().join("input.vcf");
        fs::write(&input, "input")?;
        let symlink_output = dir.path().join("symlink.vcf");
        std::os::unix::fs::symlink(&input, &symlink_output)?;
        let error = OutputTransaction::files([&input], [&symlink_output])
            .err()
            .context("expected symlink collision")?;
        assert!(error.to_string().contains("overwrite input"));

        fs::remove_file(&symlink_output)?;
        let hardlink_output = dir.path().join("hardlink.vcf");
        fs::hard_link(&input, &hardlink_output)?;
        let error = OutputTransaction::files([&input], [&hardlink_output])
            .err()
            .context("expected hardlink collision")?;
        assert!(error.to_string().contains("overwrite input"));
        assert_eq!(fs::read_to_string(input)?, "input");
        Ok(())
    }

    #[test]
    fn rollback_preserves_generation() -> Result<()> {
        let dir = tempdir()?;
        let a = dir.path().join("a");
        let b = dir.path().join("b");
        fs::write(&a, "old-a")?;
        fs::write(&b, "old-b")?;
        let tx = OutputTransaction::files(Vec::<PathBuf>::new(), [&a, &b])?;
        fs::write(tx.staged_file(&a)?, "new-a")?;
        fs::write(tx.staged_file(&b)?, "new-b")?;
        FAIL_PUBLICATION_AFTER.set(1);
        assert!(tx.commit().is_err());
        FAIL_PUBLICATION_AFTER.set(-1);
        assert_eq!(fs::read_to_string(a)?, "old-a");
        assert_eq!(fs::read_to_string(b)?, "old-b");
        Ok(())
    }

    #[test]
    fn persistent_backup_cleanup_failure_rolls_back_without_leaks() -> Result<()> {
        let dir = tempdir()?;
        let output = dir.path().join("report.csv");
        fs::write(&output, "old")?;
        let tx = OutputTransaction::files(Vec::<PathBuf>::new(), [&output])?;
        fs::write(tx.staged_file(&output)?, "new")?;
        FAIL_BACKUP_CLEANUP.set(true);
        let error = tx
            .commit()
            .expect_err("persistent cleanup failure must fail");
        FAIL_BACKUP_CLEANUP.set(false);
        assert!(
            error
                .to_string()
                .contains("persistent backup cleanup failure")
        );
        assert_eq!(fs::read_to_string(output)?, "old");
        assert!(fs::read_dir(dir.path())?.all(|entry| {
            !entry
                .expect("directory entry")
                .file_name()
                .to_string_lossy()
                .contains("backup")
        }));
        Ok(())
    }

    #[test]
    fn concurrent_generations_are_internally_consistent() -> Result<()> {
        let dir = tempdir()?;
        let prefix = dir.path().join("run");
        let barrier = Arc::new(Barrier::new(2));
        let mut workers = Vec::new();
        for generation in ["one", "two"] {
            let prefix = prefix.clone();
            let barrier = Arc::clone(&barrier);
            workers.push(thread::spawn(move || -> Result<()> {
                let tx = OutputTransaction::family(
                    Vec::<PathBuf>::new(),
                    &prefix,
                    ["summary.csv", "metrics.json"],
                )?;
                fs::write(
                    append_suffix(tx.staged_prefix()?, OsStr::new(".summary.csv")),
                    generation,
                )?;
                fs::write(
                    append_suffix(tx.staged_prefix()?, OsStr::new(".metrics.json")),
                    generation,
                )?;
                barrier.wait();
                tx.commit()
            }));
        }
        for worker in workers {
            worker.join().expect("worker panicked")?;
        }
        assert_eq!(
            fs::read_to_string(append_suffix(&prefix, OsStr::new(".summary.csv")))?,
            fs::read_to_string(append_suffix(&prefix, OsStr::new(".metrics.json")))?
        );
        assert!(fs::read_dir(dir.path())?.all(|entry| {
            !entry
                .expect("directory entry")
                .file_name()
                .to_string_lossy()
                .contains("hap-rs")
        }));
        Ok(())
    }

    #[test]
    fn overlapping_prefixes_do_not_share_a_wildcard_namespace() -> Result<()> {
        let dir = tempdir()?;
        let run = dir.path().join("run");
        let sample = dir.path().join("run.sample");
        let first = OutputTransaction::family(Vec::<PathBuf>::new(), &run, ["summary.csv"])?;
        fs::write(
            append_suffix(first.staged_prefix()?, OsStr::new(".summary.csv")),
            "run",
        )?;
        first.commit()?;
        let second = OutputTransaction::family(Vec::<PathBuf>::new(), &sample, ["summary.csv"])?;
        fs::write(
            append_suffix(second.staged_prefix()?, OsStr::new(".summary.csv")),
            "sample",
        )?;
        second.commit()?;
        assert_eq!(
            fs::read_to_string(append_suffix(&run, OsStr::new(".summary.csv")))?,
            "run"
        );
        assert_eq!(
            fs::read_to_string(append_suffix(&sample, OsStr::new(".summary.csv")))?,
            "sample"
        );
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn rejects_symlink_cycles() -> Result<()> {
        let dir = tempdir()?;
        let a = dir.path().join("a");
        let b = dir.path().join("b");
        std::os::unix::fs::symlink(&b, &a)?;
        std::os::unix::fs::symlink(&a, &b)?;
        let error = OutputTransaction::files(Vec::<PathBuf>::new(), [&a])
            .err()
            .context("expected cycle")?;
        assert!(error.to_string().contains("symlink cycle"));
        Ok(())
    }
}
