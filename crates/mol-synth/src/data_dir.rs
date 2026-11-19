#![forbid(unsafe_code)]

use std::path::PathBuf;

/// Returns the centralized Craft bio-data directory (`CRAFT_DATA_DIR`, default
/// `/home/acs/Datos/CraftBio` when `/home/acs/Datos` exists, else `$XDG_DATA_HOME/CraftBio`).
#[must_use]
pub fn craft_data_dir() -> PathBuf {
    if let Ok(raw) = std::env::var("CRAFT_DATA_DIR") {
        let trimmed = raw.trim();
        if !trimmed.is_empty() {
            return PathBuf::from(trimmed);
        }
    }
    let datos = PathBuf::from("/home/acs/Datos");
    if datos.is_dir() {
        return datos.join("CraftBio");
    }
    if let Ok(xdg) = std::env::var("XDG_DATA_HOME") {
        let trimmed = xdg.trim();
        if !trimmed.is_empty() {
            return PathBuf::from(trimmed).join("CraftBio");
        }
    }
    if let Ok(home) = std::env::var("HOME") {
        return PathBuf::from(home).join(".local/share/CraftBio");
    }
    PathBuf::from("CraftBio")
}

/// Permanent building-block catalog directory (`CRAFT_DATA_DIR/building_blocks`).
#[must_use]
pub fn building_blocks_dir() -> PathBuf {
    craft_data_dir().join("building_blocks")
}

/// Default local commercial building-block TSV catalog (`CRAFT_DATA_DIR/building_blocks/medchem_stock.tsv`).
#[must_use]
pub fn default_stock_tsv_path() -> PathBuf {
    building_blocks_dir().join("medchem_stock.tsv")
}
