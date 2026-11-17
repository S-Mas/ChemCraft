use std::fs::OpenOptions;
use std::io::{Read, Write};
use std::path::Path;

use crc32fast::Hasher;
use flate2::Compression;
use flate2::read::DeflateDecoder;
use flate2::write::DeflateEncoder;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use mol_graph::MolGraph;
use mol_sar::{SarFilter, SarTable};
use mol_vector::ReactionArrow;

/// Serialized project document inside `.chemcraft` ZIP+JSON container.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChemCraftDocument {
    pub format: String,
    pub version: String,
    pub canvas_mol: MolGraph,
    pub reaction_arrows: Vec<ReactionArrow>,
    pub sar_table: SarTable,
    pub sar_filter: SarFilter,
}

#[derive(Debug, Error)]
pub enum ProjectError {
    #[error("JSON serialization error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("invalid .chemcraft archive: {0}")]
    InvalidArchive(String),
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
}

/// Write a `.chemcraft` project atomically (`temp -> flush/sync_all -> rename`,
/// matching `craft-atomic` semantics) as a Deflate-compressed ZIP containing `manifest.json`.
pub fn save_chemcraft_atomic(path: &Path, doc: &ChemCraftDocument) -> Result<(), ProjectError> {
    let bytes = encode_chemcraft_zip(doc)?;
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    if !parent.exists() {
        std::fs::create_dir_all(parent)?;
    }
    let tmp_name = format!(
        ".{}.tmp.{}",
        path.file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("project"),
        std::process::id()
    );
    let tmp_path = parent.join(tmp_name);
    {
        #[cfg(unix)]
        let mut f = {
            use std::os::unix::fs::OpenOptionsExt;
            OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(true)
                .custom_flags(libc::O_NOFOLLOW)
                .open(&tmp_path)?
        };
        #[cfg(not(unix))]
        let mut f = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&tmp_path)?;
        f.write_all(&bytes)?;
        f.flush()?;
        f.sync_all()?;
    }
    std::fs::rename(&tmp_path, path)?;
    Ok(())
}

/// Encode `ChemCraftDocument` into pure-Rust ZIP bytes (`manifest.json` Deflate entry).
pub fn encode_chemcraft_zip(doc: &ChemCraftDocument) -> Result<Vec<u8>, ProjectError> {
    let json_bytes = serde_json::to_vec_pretty(doc)?;
    let mut hasher = Hasher::new();
    hasher.update(&json_bytes);
    let crc = hasher.finalize();

    let mut deflater = DeflateEncoder::new(Vec::new(), Compression::default());
    deflater.write_all(&json_bytes)?;
    let compressed = deflater.finish()?;

    let name = b"manifest.json";
    let mut zip = Vec::new();
    let local_offset = zip.len() as u32;

    // Local file header (30 bytes + filename)
    zip.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
    zip.extend_from_slice(&20u16.to_le_bytes()); // version needed
    zip.extend_from_slice(&0u16.to_le_bytes()); // flags
    zip.extend_from_slice(&8u16.to_le_bytes()); // method 8 = Deflate
    zip.extend_from_slice(&0u16.to_le_bytes()); // mod time
    zip.extend_from_slice(&0u16.to_le_bytes()); // mod date
    zip.extend_from_slice(&crc.to_le_bytes());
    zip.extend_from_slice(&(compressed.len() as u32).to_le_bytes());
    zip.extend_from_slice(&(json_bytes.len() as u32).to_le_bytes());
    zip.extend_from_slice(&(name.len() as u16).to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes()); // extra len
    zip.extend_from_slice(name);
    zip.extend_from_slice(&compressed);

    let cd_offset = zip.len() as u32;
    // Central directory header
    zip.extend_from_slice(&0x0201_4b50u32.to_le_bytes());
    zip.extend_from_slice(&20u16.to_le_bytes());
    zip.extend_from_slice(&20u16.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&8u16.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&crc.to_le_bytes());
    zip.extend_from_slice(&(compressed.len() as u32).to_le_bytes());
    zip.extend_from_slice(&(json_bytes.len() as u32).to_le_bytes());
    zip.extend_from_slice(&(name.len() as u16).to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&0u32.to_le_bytes());
    zip.extend_from_slice(&local_offset.to_le_bytes());
    zip.extend_from_slice(name);

    let cd_size = (zip.len() as u32).saturating_sub(cd_offset);
    // End of central directory record
    zip.extend_from_slice(&0x0605_4b50u32.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    zip.extend_from_slice(&1u16.to_le_bytes());
    zip.extend_from_slice(&1u16.to_le_bytes());
    zip.extend_from_slice(&cd_size.to_le_bytes());
    zip.extend_from_slice(&cd_offset.to_le_bytes());
    zip.extend_from_slice(&0u16.to_le_bytes());
    Ok(zip)
}

/// Decode a `.chemcraft` ZIP archive back into `ChemCraftDocument` (verifying CRC32).
pub fn decode_chemcraft_zip(bytes: &[u8]) -> Result<ChemCraftDocument, ProjectError> {
    if bytes.len() < 30 || &bytes[0..4] != b"PK\x03\x04" {
        return Err(ProjectError::InvalidArchive(
            "missing ZIP PK signature".to_string(),
        ));
    }
    let method = u16::from_le_bytes([bytes[8], bytes[9]]);
    let expected_crc = u32::from_le_bytes([bytes[14], bytes[15], bytes[16], bytes[17]]);
    let comp_len = u32::from_le_bytes([bytes[18], bytes[19], bytes[20], bytes[21]]) as usize;
    let name_len = u16::from_le_bytes([bytes[26], bytes[27]]) as usize;
    let extra_len = u16::from_le_bytes([bytes[28], bytes[29]]) as usize;
    let data_start = 30usize.saturating_add(name_len).saturating_add(extra_len);
    let data_end = data_start.saturating_add(comp_len);
    if data_end > bytes.len() {
        return Err(ProjectError::InvalidArchive(
            "truncated ZIP entry payload".to_string(),
        ));
    }
    /// Hard uncompressed payload limit (32 MiB) to neutralize ZIP decompression bombs.
    const MAX_ARCHIVE_UNCOMPRESSED_BYTES: usize = 32 * 1024 * 1024;
    let payload = &bytes[data_start..data_end];
    let raw = match method {
        0 => {
            if payload.len() > MAX_ARCHIVE_UNCOMPRESSED_BYTES {
                return Err(ProjectError::InvalidArchive(
                    "archive exceeds 32 MiB uncompressed safety limit".to_string(),
                ));
            }
            payload.to_vec()
        }
        8 => {
            let dec = DeflateDecoder::new(payload);
            let mut out = Vec::new();
            dec.take((MAX_ARCHIVE_UNCOMPRESSED_BYTES + 1) as u64)
                .read_to_end(&mut out)?;
            if out.len() > MAX_ARCHIVE_UNCOMPRESSED_BYTES {
                return Err(ProjectError::InvalidArchive(
                    "decompression bomb rejected: exceeds 32 MiB safety limit".to_string(),
                ));
            }
            out
        }
        other => {
            return Err(ProjectError::InvalidArchive(format!(
                "unsupported ZIP compression method {other}"
            )));
        }
    };
    let mut hasher = Hasher::new();
    hasher.update(&raw);
    if hasher.finalize() != expected_crc {
        return Err(ProjectError::InvalidArchive(
            "CRC32 checksum mismatch".to_string(),
        ));
    }
    Ok(serde_json::from_slice(&raw)?)
}
