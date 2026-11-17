use std::fs::File;
use std::io::Read;
use std::path::{Component, Path, PathBuf};
use std::sync::Arc;

#[cfg(unix)]
use std::os::fd::OwnedFd;

#[cfg(unix)]
use rustix::fs::{AtFlags, FileType, Mode, OFlags, ResolveFlags, mkdirat, openat, openat2, statat};
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
    #[error("refusing non-regular file (FIFO/device/socket/dir): {0}")]
    NonRegularFile(String),
    #[error("file exceeds read limit ({1} bytes): {0}")]
    FileTooLarge(String, usize),
    #[error("I/O error at {0}: {1}")]
    Io(String, std::io::Error),
}

/// Directory sandbox enforcing `--read-root` and `--write-root` isolation using pinned
/// root directory file descriptors and kernel `openat2(RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS)`
/// so intermediate-directory symlink swaps cannot escape the sandbox.
#[derive(Debug, Clone)]
pub struct PathSandbox {
    read_root: Option<PathBuf>,
    write_root: Option<PathBuf>,
    #[cfg(unix)]
    read_fd: Option<Arc<OwnedFd>>,
    #[cfg(unix)]
    write_fd: Option<Arc<OwnedFd>>,
}

impl PathSandbox {
    pub fn new(read_root: Option<&Path>, write_root: Option<&Path>) -> Result<Self, SandboxError> {
        let r = match read_root {
            Some(p) => Some(
                std::fs::canonicalize(p)
                    .map_err(|e| SandboxError::Io(p.display().to_string(), e))?,
            ),
            None => None,
        };
        let w = match write_root {
            Some(p) => Some(
                std::fs::canonicalize(p)
                    .map_err(|e| SandboxError::Io(p.display().to_string(), e))?,
            ),
            None => None,
        };
        #[cfg(unix)]
        let read_fd = match &r {
            Some(p) => Some(Arc::new(open_root_dir_fd(p)?)),
            None => None,
        };
        #[cfg(unix)]
        let write_fd = match &w {
            Some(p) => Some(Arc::new(open_root_dir_fd(p)?)),
            None => None,
        };
        Ok(Self {
            read_root: r,
            write_root: w,
            #[cfg(unix)]
            read_fd,
            #[cfg(unix)]
            write_fd,
        })
    }

    /// Atomically open a regular file inside `read_root` via `openat2` (`RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS`).
    pub fn open_read(&self, user_path: &Path) -> Result<(PathBuf, File), SandboxError> {
        let root = self.read_root.as_ref().ok_or(SandboxError::NoReadRoot)?;
        let rel = normalize_relative_under_root(root, user_path)?;
        let target = root.join(&rel);

        #[cfg(unix)]
        let file = {
            let root_fd = self.read_fd.as_ref().ok_or(SandboxError::NoReadRoot)?;
            let fd = open_beneath_no_symlinks(
                root_fd.as_ref(),
                &rel,
                OFlags::RDONLY | OFlags::NOFOLLOW | OFlags::NONBLOCK | OFlags::CLOEXEC,
                Mode::empty(),
            )
            .map_err(|e| map_unix_open_err(e, &target))?;
            File::from(fd)
        };
        #[cfg(not(unix))]
        let file =
            File::open(&target).map_err(|e| SandboxError::Io(target.display().to_string(), e))?;

        let meta = file
            .metadata()
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;
        if !meta.is_file() {
            return Err(SandboxError::NonRegularFile(target.display().to_string()));
        }
        Ok((target, file))
    }

    pub fn read_to_string_limited(
        &self,
        user_path: &Path,
        max_bytes: usize,
    ) -> Result<(PathBuf, String), SandboxError> {
        let (resolved, bytes) = self.read_bytes_limited_with_path(user_path, max_bytes)?;
        let text = String::from_utf8(bytes).map_err(|e| {
            SandboxError::Io(user_path.display().to_string(), std::io::Error::other(e))
        })?;
        Ok((resolved, text))
    }

    pub fn read_bytes_limited(
        &self,
        user_path: &Path,
        max_bytes: usize,
    ) -> Result<Vec<u8>, SandboxError> {
        let (_, bytes) = self.read_bytes_limited_with_path(user_path, max_bytes)?;
        Ok(bytes)
    }

    fn read_bytes_limited_with_path(
        &self,
        user_path: &Path,
        max_bytes: usize,
    ) -> Result<(PathBuf, Vec<u8>), SandboxError> {
        let (target, file) = self.open_read(user_path)?;
        let meta = file
            .metadata()
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;
        if meta.len() > max_bytes as u64 {
            return Err(SandboxError::FileTooLarge(
                target.display().to_string(),
                max_bytes,
            ));
        }
        let mut buf = Vec::with_capacity((meta.len() as usize).min(max_bytes));
        file.take((max_bytes + 1) as u64)
            .read_to_end(&mut buf)
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;
        if buf.len() > max_bytes {
            return Err(SandboxError::FileTooLarge(
                target.display().to_string(),
                max_bytes,
            ));
        }
        Ok((target, buf))
    }

    /// Atomically open/create a writable regular file inside `write_root` via `dirfd`-pinned
    /// parent traversal and `openat2(RESOLVE_BENEATH | RESOLVE_NO_SYMLINKS)`.
    pub fn open_write(&self, user_path: &Path) -> Result<(PathBuf, File), SandboxError> {
        let root = self.write_root.as_ref().ok_or(SandboxError::NoWriteRoot)?;
        let rel = normalize_relative_under_root(root, user_path)?;
        let target = root.join(&rel);

        #[cfg(unix)]
        let file = {
            let root_fd = self.write_fd.as_ref().ok_or(SandboxError::NoWriteRoot)?;
            let comps: Vec<&std::ffi::OsStr> = rel.iter().collect();
            let Some((file_name, parent_comps)) = comps.split_last() else {
                return Err(SandboxError::PathEscape(user_path.display().to_string()));
            };
            let mut parent_fd: Option<OwnedFd> = None;
            for comp in parent_comps {
                let base_fd = parent_fd.as_ref().map_or_else(|| root_fd.as_ref(), |f| f);
                let next_fd = open_or_create_subdir(base_fd, Path::new(comp), &target)?;
                parent_fd = Some(next_fd);
            }
            let dir_ref = parent_fd.as_ref().map_or_else(|| root_fd.as_ref(), |f| f);
            let fd = open_beneath_no_symlinks(
                dir_ref,
                Path::new(file_name),
                OFlags::WRONLY
                    | OFlags::CREATE
                    | OFlags::NOFOLLOW
                    | OFlags::NONBLOCK
                    | OFlags::CLOEXEC,
                Mode::from_bits_truncate(0o644),
            )
            .map_err(|e| map_unix_open_err(e, &target))?;
            File::from(fd)
        };
        #[cfg(not(unix))]
        let file = {
            if let Some(parent) = target.parent()
                && !parent.exists()
            {
                std::fs::create_dir_all(parent)
                    .map_err(|e| SandboxError::Io(parent.display().to_string(), e))?;
            }
            std::fs::OpenOptions::new()
                .write(true)
                .create(true)
                .open(&target)
                .map_err(|e| SandboxError::Io(target.display().to_string(), e))?
        };

        let meta = file
            .metadata()
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;
        if !meta.is_file() {
            return Err(SandboxError::NonRegularFile(target.display().to_string()));
        }
        file.set_len(0)
            .map_err(|e| SandboxError::Io(target.display().to_string(), e))?;
        Ok((target, file))
    }

    pub fn validate_write_path(&self, user_path: &Path) -> Result<PathBuf, SandboxError> {
        let root = self.write_root.as_ref().ok_or(SandboxError::NoWriteRoot)?;
        let rel = normalize_relative_under_root(root, user_path)?;
        Ok(root.join(rel))
    }
}

#[cfg(unix)]
fn open_root_dir_fd(root: &Path) -> Result<OwnedFd, SandboxError> {
    let f = std::os::unix::fs::OpenOptionsExt::custom_flags(
        std::fs::OpenOptions::new().read(true),
        libc::O_DIRECTORY | libc::O_NOFOLLOW | libc::O_CLOEXEC,
    )
    .open(root)
    .map_err(|e| SandboxError::Io(root.display().to_string(), e))?;
    Ok(OwnedFd::from(f))
}

#[cfg(unix)]
fn open_beneath_no_symlinks(
    dirfd: &OwnedFd,
    rel: &Path,
    flags: OFlags,
    mode: Mode,
) -> Result<OwnedFd, rustix::io::Errno> {
    match openat2(
        dirfd,
        rel,
        flags,
        mode,
        ResolveFlags::BENEATH | ResolveFlags::NO_SYMLINKS,
    ) {
        Ok(fd) => Ok(fd),
        Err(rustix::io::Errno::NOSYS) => {
            let comps: Vec<&std::ffi::OsStr> = rel.iter().collect();
            let Some((last, parents)) = comps.split_last() else {
                return Err(rustix::io::Errno::INVAL);
            };
            let mut cur: Option<OwnedFd> = None;
            for comp in parents {
                let base = cur.as_ref().unwrap_or(dirfd);
                let next = openat(
                    base,
                    Path::new(comp),
                    OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC,
                    Mode::empty(),
                )?;
                cur = Some(next);
            }
            let base = cur.as_ref().unwrap_or(dirfd);
            openat(base, Path::new(last), flags | OFlags::NOFOLLOW, mode)
        }
        Err(e) => Err(e),
    }
}

#[cfg(unix)]
fn open_or_create_subdir(
    dirfd: &OwnedFd,
    comp: &Path,
    full_target: &Path,
) -> Result<OwnedFd, SandboxError> {
    let dir_flags = OFlags::RDONLY | OFlags::DIRECTORY | OFlags::NOFOLLOW | OFlags::CLOEXEC;
    match open_beneath_no_symlinks(dirfd, comp, dir_flags, Mode::empty()) {
        Ok(fd) => Ok(fd),
        Err(rustix::io::Errno::NOENT) => {
            match mkdirat(dirfd, comp, Mode::from_bits_truncate(0o755)) {
                Ok(()) | Err(rustix::io::Errno::EXIST) => {}
                Err(e) => return Err(map_unix_open_err(e, full_target)),
            }
            open_beneath_no_symlinks(dirfd, comp, dir_flags, Mode::empty())
                .map_err(|e| map_unix_open_err(e, full_target))
        }
        Err(e) => {
            if e == rustix::io::Errno::NOTDIR
                && let Ok(st) = statat(dirfd, comp, AtFlags::SYMLINK_NOFOLLOW)
                && FileType::from_raw_mode(st.st_mode).is_symlink()
            {
                return Err(SandboxError::SymlinkForbidden(
                    full_target.display().to_string(),
                ));
            }
            Err(map_unix_open_err(e, full_target))
        }
    }
}

#[cfg(unix)]
fn map_unix_open_err(err: rustix::io::Errno, target: &Path) -> SandboxError {
    if err == rustix::io::Errno::LOOP {
        SandboxError::SymlinkForbidden(target.display().to_string())
    } else if err == rustix::io::Errno::XDEV {
        SandboxError::PathEscape(target.display().to_string())
    } else if matches!(
        err,
        rustix::io::Errno::NOTDIR | rustix::io::Errno::ISDIR | rustix::io::Errno::NXIO
    ) {
        SandboxError::NonRegularFile(target.display().to_string())
    } else {
        SandboxError::Io(
            target.display().to_string(),
            std::io::Error::from_raw_os_error(err.raw_os_error()),
        )
    }
}

fn normalize_relative_under_root(root: &Path, candidate: &Path) -> Result<PathBuf, SandboxError> {
    let joined = if candidate.is_absolute() {
        normalize_lexical(candidate)
    } else {
        normalize_lexical(&root.join(candidate))
    };
    let rel = joined
        .strip_prefix(root)
        .map_err(|_| SandboxError::PathEscape(candidate.display().to_string()))?
        .to_path_buf();
    if rel.as_os_str().is_empty() {
        return Err(SandboxError::PathEscape(candidate.display().to_string()));
    }
    Ok(rel)
}

fn normalize_lexical(path: &Path) -> PathBuf {
    let mut out = PathBuf::new();
    for comp in path.components() {
        match comp {
            Component::CurDir => {}
            Component::ParentDir => {
                out.pop();
            }
            other => out.push(other.as_os_str()),
        }
    }
    out
}
