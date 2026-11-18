//! DockCraft interop: export 3D SDF V2000 / PDB / SMILES bundles for single ligands
//! or multi-molecule libraries and optionally launch the standalone `dockcraft` studio.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use mol_graph::{MolGraph, write_pdb, write_sdf_v2000_3d};
use mol_sar::SarTable;
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

// ---------------------------------------------------------------------------------------------
// Return leg: folding DockCraft's results back into the SAR table.
// ---------------------------------------------------------------------------------------------

/// One compound's docking results, as DockCraft's `*.dockcraft.json` sidecar carries them.
///
/// Deserialised permissively on purpose: DockCraft owns this format, and a field it adds later
/// must not stop an older ChemCraft from reading the ones it already understands.
#[derive(Debug, Clone, serde::Deserialize)]
pub struct DockCraftResultRow {
    /// Canonical SMILES as *this* program spelled it when the library was exported. The join key.
    #[serde(default)]
    pub canonical_smiles: String,
    /// SAR row id from the project the library came from, when there was one.
    #[serde(default)]
    pub source_row_id: Option<usize>,
    #[serde(default)]
    pub name: String,
    #[serde(default)]
    pub rank: usize,
    #[serde(default)]
    pub affinity_kcal: f64,
    #[serde(default)]
    pub hyde_delta_g_kcal: f64,
    #[serde(default)]
    pub estimated_ki_nm: f64,
    #[serde(default)]
    pub ligand_efficiency: f64,
    #[serde(default)]
    pub mpo_score: f64,
    #[serde(default)]
    pub measured_pic50: Option<f64>,
    /// Path to the docked pose, relative to the sidecar file.
    #[serde(default)]
    pub pose_sdf: String,
}

/// How the docking run was set up.
#[derive(Debug, Clone, Default, serde::Deserialize)]
pub struct DockCraftRun {
    #[serde(default)]
    pub receptor_name: String,
    #[serde(default)]
    pub scoring: String,
    #[serde(default)]
    pub exhaustiveness: usize,
    #[serde(default)]
    pub box_center: [f64; 3],
    #[serde(default)]
    pub box_size: [f64; 3],
    #[serde(default)]
    pub dockcraft_version: String,
}

#[derive(Debug, Clone, serde::Deserialize)]
pub struct DockCraftResults {
    #[serde(default)]
    pub format: String,
    #[serde(default)]
    pub version: String,
    #[serde(default)]
    pub run: DockCraftRun,
    #[serde(default)]
    pub rows: Vec<DockCraftResultRow>,
}

/// What a merge did, so the caller can say so instead of guessing.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MergeOutcome {
    /// Rows of the SAR table that received docking results.
    pub merged: usize,
    /// Sidecar rows that matched nothing in the table.
    pub unmatched: Vec<String>,
}

/// Cap on a sidecar read, in the spirit of every other bounded reader here.
const MAX_SIDECAR_BYTES: u64 = 64 * 1024 * 1024;

/// Reads a `*.dockcraft.json` sidecar from disk.
pub fn read_dockcraft_results(path: &Path) -> Result<DockCraftResults, DockCraftError> {
    let meta = std::fs::metadata(path)
        .map_err(|e| DockCraftError::Io(format!("stat {}: {e}", path.display())))?;
    if meta.len() > MAX_SIDECAR_BYTES {
        return Err(DockCraftError::Io(format!(
            "{} is {} bytes, over the {MAX_SIDECAR_BYTES} cap",
            path.display(),
            meta.len()
        )));
    }
    let text = std::fs::read_to_string(path)
        .map_err(|e| DockCraftError::Io(format!("read {}: {e}", path.display())))?;
    let results: DockCraftResults = serde_json::from_str(&text)
        .map_err(|e| DockCraftError::Io(format!("parse {}: {e}", path.display())))?;
    if !results.format.is_empty() && results.format != "dockcraft-results" {
        return Err(DockCraftError::Io(format!(
            "{} declares format '{}', expected 'dockcraft-results'",
            path.display(),
            results.format
        )));
    }
    Ok(results)
}

/// Property keys the merge writes. Prefixed so a docking column is never confused with a measured
/// one, which is the whole point of keeping them side by side.
const KEY_AFFINITY: &str = "DockCraft_affinity_kcal";
const KEY_HYDE: &str = "DockCraft_hyde_dG_kcal";
const KEY_KI: &str = "DockCraft_estimated_Ki_nM";
const KEY_LE: &str = "DockCraft_ligand_efficiency";
const KEY_MPO: &str = "DockCraft_mpo";
const KEY_RANK: &str = "DockCraft_rank";
const KEY_POSE: &str = "DockCraft_pose_sdf";
const KEY_RECEPTOR: &str = "DockCraft_receptor";
const KEY_SCORING: &str = "DockCraft_scoring";

/// Folds docking results into a SAR table, matching on the SAR row id first and the canonical
/// SMILES second.
///
/// The results land in each row's `properties` map under `DockCraft_*` keys, so the SAR explorer
/// can filter, sort, and compute matched molecular pairs and activity cliffs on a docking score
/// exactly as it does on a measured one. Nothing measured is ever overwritten: a computed column
/// and an experimental column have to stay distinguishable.
///
/// `pose_dir` is the directory the sidecar lives in; pose paths are relative to it and are stored
/// resolved so the row points at a file that exists.
pub fn merge_dockcraft_results(
    table: &mut SarTable,
    results: &DockCraftResults,
    pose_dir: &Path,
) -> MergeOutcome {
    let mut outcome = MergeOutcome::default();
    for row in &results.rows {
        let by_id = table
            .rows
            .iter()
            .position(|r| row.source_row_id.is_some_and(|id| r.id == id));
        let idx = by_id.or_else(|| {
            if row.canonical_smiles.is_empty() {
                None
            } else {
                table
                    .rows
                    .iter()
                    .position(|r| r.canonical_smiles == row.canonical_smiles)
            }
        });
        let Some(idx) = idx else {
            let label = if row.name.is_empty() {
                row.canonical_smiles.clone()
            } else {
                row.name.clone()
            };
            outcome.unmatched.push(label);
            continue;
        };

        let props = &mut table.rows[idx].properties;
        props.insert(
            KEY_AFFINITY.to_string(),
            format!("{:.3}", row.affinity_kcal),
        );
        props.insert(
            KEY_HYDE.to_string(),
            format!("{:.3}", row.hyde_delta_g_kcal),
        );
        props.insert(KEY_KI.to_string(), format!("{:.4}", row.estimated_ki_nm));
        props.insert(KEY_LE.to_string(), format!("{:.4}", row.ligand_efficiency));
        props.insert(KEY_MPO.to_string(), format!("{:.4}", row.mpo_score));
        props.insert(KEY_RANK.to_string(), row.rank.to_string());
        if !row.pose_sdf.is_empty() {
            props.insert(
                KEY_POSE.to_string(),
                pose_dir.join(&row.pose_sdf).display().to_string(),
            );
        }
        if !results.run.receptor_name.is_empty() {
            props.insert(KEY_RECEPTOR.to_string(), results.run.receptor_name.clone());
        }
        if !results.run.scoring.is_empty() {
            props.insert(KEY_SCORING.to_string(), results.run.scoring.clone());
        }
        outcome.merged += 1;
    }
    outcome
}
