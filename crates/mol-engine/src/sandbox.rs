use std::fs::{File, OpenOptions};
use std::io::Read;
use std::path::{Component, Path, PathBuf};

use thiserror::Error;

/// Default maximum read size for input chemical files / archives (16 MiB).
pub const MAX_READ_BYTES: usize = 16 * 1024 * 1024;

#[derive(Debug, Error)]
pub enum SandboxError {
    #[error("read access denied: no --read-root configured")]
    NoReadRoot,
    #[error("write access denied: no --write-root configured")]
    NoWriteRoot,
    #[error("symbolic links are forbidden inside sandbox: {0}")]
    SymlinkForbidden(String),
    #[error("path escapes sandbox root: {0}")]
    PathEscape(String),
    #[error("refusing non-regular file (FIFO/device/socket): {0}")]
    NonRegularFile(String),
    #[error("file exceeds read limit ({1} bytes): {0}")]
    FileTooLarge(String, usize),
    #[error("I/O error at {0}: {1}")]
    Io(String, std::io::Error),
}

/// Hardened directory confinement (`symlink_metadata`, `O_NOFOLLOW | O_NONBLOCK`, no TOCTOU).
#[derive(Debug, Clone)]
pub struct PathSandbox {
    read_root: Option<PathBuf>,
    write_root: Option<PathBuf>,
}

impl PathSandbox {
    pub fn new(read_root: Option<&Path>, write_root: Option<&Path>) -> Result<Self, SandboxError> {
        let r = match read_root {
            Some(p) => Some(normalize_root(p, false)?),
            None => None,
        };
        let w = match write_root {
            Some(p) => Some(normalize_root(p, true)?),
            None => None,
        };
        Ok(Self {
            read_root: r,
            write_root: w,
        })
    }

    /// Resolve and read a bounded UTF-8 file inside `read_root` (no symlinks, no FIFOs).
    pub fn read_to_string_limited(
        &self,
        user_path: &Path,
        max_bytes: usize,
    ) -> Result<(PathBuf, String), SandboxError> {
        let bytes = self.read_bytes_limited(user_path, max_bytes)?;
        let text = String::from_utf8(bytes).map_err(|e| {
            SandboxError::Io(user_path.display().to_string(), std::io::Error::other(e))
        })?;
        let resolved = self.resolve_read_path(user_path)?;
        Ok((resolved, text))
    }

    /// Read bounded binary bytes (e.g. `.chemcraft` ZIP archives) safely.
    pub fn read_bytes_limited(
        &self,
        user_path: &Path,
        max_bytes: usize,
    ) -> Result<Vec<u8>, SandboxError> {
        let target = self.resolve_read_path(user_path)?;
        #[cfg(unix)]
        let mut file = {
            use std::os::unix::fs::OpenOptionsExt;
            OpenOptions::new()
                .read(true)
                .custom_flags(libc::O_NOFOLLOW | libc::O_NONBLOCK)
                .open(&target)
                .map_err(|e| SandboxError::Io(target.display().to_string(), e))?
        };
        #[cfg(not(unix))]
        let mut file = OpenOptions::new()
            .read(true)
            .open(&target)
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;

        let meta = file
            .metadata()
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;
        if !meta.is_file() {
            return Err(SandboxError::NonRegularFile(target.display().to_string()));
        }
        if meta.len() > max_bytes as u64 {
            return Err(SandboxError::FileTooLarge(
                target.display().to_string(),
                max_bytes,
            ));
        }

        let mut buf = Vec::with_capacity((meta.len() as usize).min(max_bytes));
        let mut limited = (&mut file).take((max_bytes + 1) as u64);
        limited
            .read_to_end(&mut buf)
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;
        if buf.len() > max_bytes {
            return Err(SandboxError::FileTooLarge(
                target.display().to_string(),
                max_bytes,
            ));
        }
        Ok(buf)
    }

    pub fn resolve_read_path(&self, user_path: &Path) -> Result<PathBuf, SandboxError> {
        let root = self.read_root.as_ref().ok_or(SandboxError::NoReadRoot)?;
        let candidate = resolve_relative_clean(root, user_path)?;
        check_no_symlinks_under(root, &candidate)?;
        let meta = std::fs::symlink_metadata(&candidate)
            .map_err(|e| SandboxError::Io(candidate.display().to_string(), e))?;
        if meta.file_type().is_symlink() {
            return Err(SandboxError::SymlinkForbidden(
                candidate.display().to_string(),
            ));
        }
        if !meta.is_file() {
            return Err(SandboxError::NonRegularFile(
                candidate.display().to_string(),
            ));
        }
        Ok(candidate)
    }

    /// Atomically validate and open an output file inside `write_root` with `O_NOFOLLOW`.
    pub fn open_write(&self, user_path: &Path) -> Result<(PathBuf, File), SandboxError> {
        let root = self.write_root.as_ref().ok_or(SandboxError::NoWriteRoot)?;
        let candidate = resolve_relative_clean(root, user_path)?;
        check_no_symlinks_under(root, &candidate)?;
        if let Some(parent) = candidate.parent()
            && !parent.exists()
        {
            std::fs::create_dir_all(parent)
                .map_err(|e| SandboxError::Io(parent.display().to_string(), e))?;
        }
        if let Ok(meta) = std::fs::symlink_metadata(&candidate) {
            if meta.file_type().is_symlink() {
                return Err(SandboxError::SymlinkForbidden(
                    candidate.display().to_string(),
                ));
            }
            if !meta.is_file() {
                return Err(SandboxError::NonRegularFile(
                    candidate.display().to_string(),
                ));
            }
        }

        #[cfg(unix)]
        let file = {
            use std::os::unix::fs::OpenOptionsExt;
            OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(true)
                .custom_flags(libc::O_NOFOLLOW | libc::O_NONBLOCK)
                .open(&candidate)
                .map_err(|e| SandboxError::Io(candidate.display().to_string(), e))?
        };
        #[cfg(not(unix))]
        let file = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&candidate)
            .map_err(|e| SandboxError::Io(candidate.display().to_string(), e))?;

        Ok((candidate, file))
    }

    /// Return resolved write path after verifying sandbox containment (for atomic rename workflows).
    pub fn validate_write_path(&self, user_path: &Path) -> Result<PathBuf, SandboxError> {
        let root = self.write_root.as_ref().ok_or(SandboxError::NoWriteRoot)?;
        let candidate = resolve_relative_clean(root, user_path)?;
        check_no_symlinks_under(root, &candidate)?;
        Ok(candidate)
    }
}

fn normalize_root(p: &Path, create_if_missing: bool) -> Result<PathBuf, SandboxError> {
    if create_if_missing && !p.exists() {
        std::fs::create_dir_all(p).map_err(|e| SandboxError::Io(p.display().to_string(), e))?;
    }
    std::fs::canonicalize(p).map_err(|e| SandboxError::Io(p.display().to_string(), e))
}

fn resolve_relative_clean(root: &Path, user_path: &Path) -> Result<PathBuf, SandboxError> {
    let mut cur = if user_path.is_absolute() {
        let stripped = user_path
            .strip_prefix(root)
            .map_err(|_| SandboxError::PathEscape(user_path.display().to_string()))?;
        let mut p = root.to_path_buf();
        for comp in stripped.components() {
            apply_comp(&mut p, root, comp, user_path)?;
        }
        p
    } else {
        let mut p = root.to_path_buf();
        for comp in user_path.components() {
            apply_comp(&mut p, root, comp, user_path)?;
        }
        p
    };
    if !cur.starts_with(root) {
        return Err(SandboxError::PathEscape(user_path.display().to_string()));
    }
    if cur == root {
        return Err(SandboxError::PathEscape(user_path.display().to_string()));
    }
    Ok(std::mem::take(&mut cur))
}

fn apply_comp(
    cur: &mut PathBuf,
    root: &Path,
    comp: Component<'_>,
    orig: &Path,
) -> Result<(), SandboxError> {
    match comp {
        Component::Normal(s) => cur.push(s),
        Component::CurDir => {}
        Component::ParentDir => {
            if cur == root || !cur.pop() || !cur.starts_with(root) {
                return Err(SandboxError::PathEscape(orig.display().to_string()));
            }
        }
        Component::RootDir | Component::Prefix(_) => {
            return Err(SandboxError::PathEscape(orig.display().to_string()));
        }
    }
    Ok(())
}

fn check_no_symlinks_under(root: &Path, target: &Path) -> Result<(), SandboxError> {
    let rel = target
        .strip_prefix(root)
        .map_err(|_| SandboxError::PathEscape(target.display().to_string()))?;
    let mut walk = root.to_path_buf();
    for comp in rel.components() {
        walk.push(comp);
        if let Ok(meta) = std::fs::symlink_metadata(&walk)
            && meta.file_type().is_symlink()
        {
            return Err(SandboxError::SymlinkForbidden(walk.display().to_string()));
        }
    }
    Ok(())
}
