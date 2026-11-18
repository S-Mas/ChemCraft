//! DockCraft interop: export 3D SDF V2000 / PDB / SMILES bundles for single ligands
//! or multi-molecule libraries and optionally launch the standalone `dockcraft` studio.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use mol_graph::{MolGraph, write_pdb, write_sdf_v2000_3d};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum DockCraftError {
    #[error("no molecules available to send to DockCraft")]
    EmptyInput,
    #[error("DockCraft handoff I/O failed: {0}")]
    Io(String),
    #[error("could not launch `dockcraft` (3D library exported at {bundle_dir}): {msg}")]
    Launch { bundle_dir: String, msg: String },
}

#[derive(Debug, Clone)]
pub struct DockCraftBundle {
    pub dir: PathBuf,
    pub sdf_path: PathBuf,
    pub smi_path: PathBuf,
    pub pdb_path: Option<PathBuf>,
    pub molecule_count: usize,
}

/// Default handoff directory for ligands & screening libraries sent to `DockCraft`.
pub fn default_dockcraft_handoff_dir() -> PathBuf {
    if let Ok(xdg_data) = std::env::var("XDG_DATA_HOME")
        && !xdg_data.trim().is_empty()
    {
        return PathBuf::from(xdg_data)
            .join("DockCraft")
            .join("chemcraft_handoff");
    }
    if let Ok(home) = std::env::var("HOME")
        && !home.trim().is_empty()
    {
        return PathBuf::from(home)
            .join(".local")
            .join("share")
            .join("DockCraft")
            .join("chemcraft_handoff");
    }
    PathBuf::from("chemcraft_dockcraft")
}

/// Check if the `dockcraft` GUI binary exists on `$PATH`.
pub fn is_dockcraft_available() -> bool {
    binary_on_path("dockcraft")
}

/// Check if the `dockcraft-cli` binary exists on `$PATH`.
pub fn is_dockcraft_cli_available() -> bool {
    binary_on_path("dockcraft-cli")
}

fn binary_on_path(bin: &str) -> bool {
    if let Ok(path_var) = std::env::var("PATH") {
        for dir in std::env::split_paths(&path_var) {
            if dir.join(bin).is_file() {
                return true;
            }
        }
    }
    false
}

fn sanitize_stem(raw: &str) -> String {
    let s: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' || c == '-' {
                c
            } else {
                '_'
            }
        })
        .take(40)
        .collect();
    let trimmed = s.trim_matches('_').to_string();
    if trimmed.is_empty() {
        "chemcraft_docking".to_string()
    } else {
        trimmed
    }
}

/// Export one molecule or a multi-molecule screening library as 3D SDF V2000 + SMILES (+ PDB for single structure)
/// ready for `DockCraft` (`dockcraft` GUI or `dockcraft-cli screen --library <sdf>`).
pub fn write_dockcraft_bundle(
    mols: &[&MolGraph],
    dir: &Path,
    stem_override: Option<&str>,
) -> Result<DockCraftBundle, DockCraftError> {
    let non_empty: Vec<&MolGraph> = mols
        .iter()
        .copied()
        .filter(|m| !m.atoms.is_empty())
        .collect();
    if non_empty.is_empty() {
        return Err(DockCraftError::EmptyInput);
    }

    std::fs::create_dir_all(dir)
        .map_err(|e| DockCraftError::Io(format!("create_dir_all {}: {e}", dir.display())))?;

    let stem = sanitize_stem(
        stem_override
            .or_else(|| non_empty.first().and_then(|m| m.name.as_deref()))
            .unwrap_or("chemcraft_ligand"),
    );

    let sdf_path = dir.join(format!("{stem}_3d.sdf"));
    let smi_path = dir.join(format!("{stem}.smi"));

    let sdf_text = write_sdf_v2000_3d(&non_empty);
    std::fs::write(&sdf_path, sdf_text.as_bytes())
        .map_err(|e| DockCraftError::Io(format!("write {}: {e}", sdf_path.display())))?;

    let mut smi_text = String::new();
    for (idx, m) in non_empty.iter().enumerate() {
        let can = m.canonical_smiles();
        if can.is_empty() {
            continue;
        }
        let name = m
            .name
            .as_deref()
            .map(|n| n.replace(['\r', '\n'], " "))
            .unwrap_or_else(|| format!("Mol_{}", idx + 1));
        let _ = writeln!(smi_text, "{can} {name}");
    }
    std::fs::write(&smi_path, smi_text.as_bytes())
        .map_err(|e| DockCraftError::Io(format!("write {}: {e}", smi_path.display())))?;

    let pdb_path = if non_empty.len() == 1 {
        let p = dir.join(format!("{stem}.pdb"));
        std::fs::write(&p, write_pdb(non_empty[0]).as_bytes())
            .map_err(|e| DockCraftError::Io(format!("write {}: {e}", p.display())))?;
        Some(p)
    } else {
        None
    };

    Ok(DockCraftBundle {
        dir: dir.to_path_buf(),
        sdf_path,
        smi_path,
        pdb_path,
        molecule_count: non_empty.len(),
    })
}

/// Export a 3D ligand/library bundle and spawn the standalone `dockcraft` GUI desktop app.
pub fn launch_dockcraft_for_mols(
    mols: &[&MolGraph],
    dir: Option<&Path>,
    stem_override: Option<&str>,
) -> Result<DockCraftBundle, DockCraftError> {
    let target_dir = dir
        .map(Path::to_path_buf)
        .unwrap_or_else(default_dockcraft_handoff_dir);
    let bundle = write_dockcraft_bundle(mols, &target_dir, stem_override)?;
    Command::new("dockcraft")
        .current_dir(&bundle.dir)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| DockCraftError::Launch {
            bundle_dir: bundle.dir.display().to_string(),
            msg: e.to_string(),
        })?;
    Ok(bundle)
}
