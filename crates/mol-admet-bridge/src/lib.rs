//! Layer 2 (`mol-admet-bridge`): Self-contained pure-Rust ADMET/druglikeness engine
//! behind a decoupled `AdmetProvider` trait, with optional `admetcraft-cli` IPC discovery
//! and built-in fallback so clean clones compile without sibling repositories.

#![forbid(unsafe_code)]

pub mod descriptors;

use std::fmt::Write as _;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

pub use descriptors::{
    AbsorptionLevel, BioavailabilityRadar, BoiledEggPoint, Druglikeness, Lipophilicity,
    MedicinalChemistry, MoleculeEvaluation, Pharmacokinetics, PhysicochemicalProperties,
    RuleResult, SolubilityClass, WaterSolubility, WaterSolubilityModel, evaluate_mol_builtin,
};
use flate2::Compression;
use flate2::write::ZlibEncoder;
use mol_graph::parse_smiles;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Supported ADMET export formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmetReportFormat {
    Csv,
    Json,
    Pdf,
    SvgBoiledEgg,
}

#[derive(Debug, Error)]
pub enum BridgeError {
    #[error("ADMET calculation failed: {0}")]
    Calc(String),
    #[error("ADMET report export failed: {0}")]
    Report(String),
}

/// Compact cached ADMET summary stored alongside each SAR row and displayed live.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AdmetSummary {
    pub formula: String,
    pub molecular_weight: f64,
    pub consensus_logp: f64,
    pub tpsa: f64,
    pub esol_log_s: f64,
    pub gi_absorption: String,
    pub bbb_permeant: bool,
    pub lipinski_violations: usize,
    pub pains_count: usize,
    pub brenk_count: usize,
    pub bioavailability_score: f64,
    pub radar_axes: [f64; 6],
}

impl From<&MoleculeEvaluation> for AdmetSummary {
    fn from(ev: &MoleculeEvaluation) -> Self {
        Self {
            formula: ev.physicochemical.formula.clone(),
            molecular_weight: ev.physicochemical.molecular_weight,
            consensus_logp: ev.lipophilicity.consensus_log_p,
            tpsa: ev.physicochemical.tpsa,
            esol_log_s: ev.water_solubility.esol.log_s,
            gi_absorption: ev.pharmacokinetics.gi_absorption.as_str().to_string(),
            bbb_permeant: ev.pharmacokinetics.bbb_permeant,
            lipinski_violations: ev.druglikeness.lipinski.violations,
            pains_count: ev.medicinal_chemistry.pains_alerts.len(),
            brenk_count: ev.medicinal_chemistry.brenk_alerts.len(),
            bioavailability_score: ev.druglikeness.bioavailability_score,
            radar_axes: ev.bioavailability_radar.normalized_axes,
        }
    }
}

/// Property delta produced when mutating a scaffold / applying a bioisostere swap.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AdmetDelta {
    pub old_smiles: String,
    pub new_smiles: String,
    pub delta_mw: f64,
    pub delta_logp: f64,
    pub delta_tpsa: f64,
    pub delta_esol_log_s: f64,
    pub old_gi: String,
    pub new_gi: String,
    pub old_bbb: bool,
    pub new_bbb: bool,
    pub old_lipinski_violations: usize,
    pub new_lipinski_violations: usize,
}

/// Decoupled trait for plugging in built-in descriptors or an external `admetcraft-cli` engine.
pub trait AdmetProvider: Send + Sync {
    fn evaluate_smiles(
        &self,
        smiles: &str,
        name: Option<&str>,
    ) -> Result<MoleculeEvaluation, BridgeError>;

    fn evaluate_batch(
        &self,
        items: &[(String, Option<String>)],
    ) -> Vec<Result<MoleculeEvaluation, BridgeError>> {
        items
            .par_iter()
            .map(|(smi, name)| self.evaluate_smiles(smi, name.as_deref()))
            .collect()
    }

    fn export_report_to_writer<W: Write>(
        &self,
        evals: &[MoleculeEvaluation],
        writer: &mut W,
        format: AdmetReportFormat,
    ) -> Result<(), BridgeError>;
}

/// Default self-contained ADMET provider (`ChemCraft` built-in engine, zero sibling-repo dependency).
/// Optionally delegates single evaluations to `admetcraft-cli eval` when `prefer_external_cli` is enabled.
#[derive(Debug, Default, Clone, Copy)]
pub struct AdmetCraftProvider {
    pub prefer_external_cli: bool,
}

impl AdmetCraftProvider {
    pub fn is_external_cli_available() -> bool {
        Command::new("admetcraft-cli")
            .arg("--version")
            .stdin(Stdio::null())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|s| s.success())
    }

    pub fn is_external_gui_available() -> bool {
        if let Ok(path_var) = std::env::var("PATH") {
            for dir in std::env::split_paths(&path_var) {
                if dir.join("admetcraft").is_file() {
                    return true;
                }
            }
        }
        false
    }

    pub fn evaluate_with_optional_cli(
        smiles: &str,
        name: Option<&str>,
    ) -> Option<MoleculeEvaluation> {
        let mut cmd = Command::new("admetcraft-cli");
        cmd.args(["eval", smiles]);
        if let Some(n) = name {
            cmd.args(["--name", n]);
        }
        let out = cmd
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .output()
            .ok()?;
        if !out.status.success() {
            return None;
        }
        let mut ev: MoleculeEvaluation = serde_json::from_slice(&out.stdout).ok()?;
        if let Some(n) = name {
            ev.name = n.to_string();
        }
        Some(ev)
    }

    /// Generate an official AdmetCraft report (`csv`/`json`/`pdf`/`svg` BOILED-Egg) using
    /// `admetcraft-cli batch` inside `work_dir` when installed, or fall back to built-in export.
    pub fn run_full_report_with_handoff(
        &self,
        items: &[(String, String)],
        work_dir: &Path,
        out_path: &Path,
        format: AdmetReportFormat,
    ) -> Result<String, BridgeError> {
        if items.is_empty() {
            return Err(BridgeError::Report("empty evaluation list".to_string()));
        }
        std::fs::create_dir_all(work_dir)
            .map_err(|e| BridgeError::Report(format!("mkdir {}: {e}", work_dir.display())))?;

        let smi_path = work_dir.join("chemcraft_handoff.smi");
        let mut smi_buf = String::new();
        for (smi, name) in items {
            let clean_name = name.replace(['\r', '\n'], " ");
            let _ = writeln!(smi_buf, "{smi} {clean_name}");
        }
        std::fs::write(&smi_path, smi_buf.as_bytes())
            .map_err(|e| BridgeError::Report(format!("write {}: {e}", smi_path.display())))?;

        let fmt_flag = match format {
            AdmetReportFormat::Csv => "csv",
            AdmetReportFormat::Json => "json",
            AdmetReportFormat::Pdf => "pdf",
            AdmetReportFormat::SvgBoiledEgg => "svg",
        };
        if let Ok(out) = Command::new("admetcraft-cli")
            .arg("batch")
            .arg(&smi_path)
            .arg("--out")
            .arg(out_path)
            .arg("--format")
            .arg(fmt_flag)
            .arg("--read-root")
            .arg(work_dir)
            .arg("--write-root")
            .arg(out_path.parent().unwrap_or(work_dir))
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()
            && out.status.success()
        {
            return Ok(format!(
                "Delegated to AdmetCraft CLI -> {}",
                out_path.display()
            ));
        }

        if format == AdmetReportFormat::SvgBoiledEgg {
            return Err(BridgeError::Report(
                "BOILED-Egg SVG export requires external `admetcraft-cli` on PATH".to_string(),
            ));
        }

        let batch_items: Vec<(String, Option<String>)> = items
            .iter()
            .map(|(s, n)| (s.clone(), Some(n.clone())))
            .collect();
        let evals = self.evaluate_batch(&batch_items);
        let ok_evals: Vec<_> = evals.into_iter().filter_map(Result::ok).collect();
        let mut file = std::fs::File::create(out_path)
            .map_err(|e| BridgeError::Report(format!("create {}: {e}", out_path.display())))?;
        self.export_report_to_writer(&ok_evals, &mut file, format)?;
        Ok(format!(
            "Exported built-in ADMET report -> {}",
            out_path.display()
        ))
    }
}

/// Write an optional `.smi` handoff file and spawn the standalone `admetcraft` GUI desktop app.
pub fn launch_admetcraft_gui(
    items: &[(String, String)],
    export_dir: &Path,
) -> Result<PathBuf, BridgeError> {
    std::fs::create_dir_all(export_dir)
        .map_err(|e| BridgeError::Report(format!("mkdir {}: {e}", export_dir.display())))?;
    let smi_path = export_dir.join("chemcraft_to_admetcraft.smi");
    let mut smi_buf = String::new();
    for (smi, name) in items {
        let _ = writeln!(smi_buf, "{smi} {}", name.replace(['\r', '\n'], " "));
    }
    std::fs::write(&smi_path, smi_buf.as_bytes())
        .map_err(|e| BridgeError::Report(format!("write {}: {e}", smi_path.display())))?;

    Command::new("admetcraft")
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| {
            BridgeError::Report(format!(
                "failed to launch `admetcraft` GUI (saved SMILES at {}): {e}",
                smi_path.display()
            ))
        })?;
    Ok(smi_path)
}

impl AdmetProvider for AdmetCraftProvider {
    fn evaluate_smiles(
        &self,
        smiles: &str,
        name: Option<&str>,
    ) -> Result<MoleculeEvaluation, BridgeError> {
        if self.prefer_external_cli
            && let Some(ev) = Self::evaluate_with_optional_cli(smiles, name)
        {
            return Ok(ev);
        }
        let mol = parse_smiles(smiles).map_err(|e| BridgeError::Calc(e.to_string()))?;
        Ok(evaluate_mol_builtin(&mol, name))
    }

    fn export_report_to_writer<W: Write>(
        &self,
        evals: &[MoleculeEvaluation],
        writer: &mut W,
        format: AdmetReportFormat,
    ) -> Result<(), BridgeError> {
        if evals.is_empty() {
            return Err(BridgeError::Report("empty evaluation list".to_string()));
        }
        match format {
            AdmetReportFormat::Csv => {
                let mut wtr = csv::Writer::from_writer(Vec::new());
                wtr.write_record([
                    "Molecule",
                    "Canonical SMILES",
                    "Formula",
                    "MW",
                    "WLOGP",
                    "XLOGP3",
                    "MLOGP",
                    "Consensus LogP",
                    "TPSA",
                    "ESOL LogS",
                    "GI Absorption",
                    "BBB Permeant",
                    "Lipinski #violations",
                    "PAINS #alerts",
                    "Brenk #alerts",
                    "Bioavailability Score",
                ])
                .map_err(|e| BridgeError::Report(e.to_string()))?;

                for e in evals {
                    let p = &e.physicochemical;
                    let l = &e.lipophilicity;
                    let w = &e.water_solubility;
                    let pk = &e.pharmacokinetics;
                    let d = &e.druglikeness;
                    let mc = &e.medicinal_chemistry;
                    wtr.write_record([
                        sanitize_csv(&e.name),
                        e.canonical_smiles.clone(),
                        p.formula.clone(),
                        format!("{:.2}", p.molecular_weight),
                        format!("{:.2}", l.wlogp),
                        format!("{:.2}", l.xlogp3),
                        format!("{:.2}", l.mlogp),
                        format!("{:.2}", l.consensus_log_p),
                        format!("{:.2}", p.tpsa),
                        format!("{:.2}", w.esol.log_s),
                        pk.gi_absorption.as_str().to_string(),
                        if pk.bbb_permeant { "Yes" } else { "No" }.to_string(),
                        d.lipinski.violations.to_string(),
                        mc.pains_alerts.len().to_string(),
                        mc.brenk_alerts.len().to_string(),
                        format!("{:.2}", d.bioavailability_score),
                    ])
                    .map_err(|err| BridgeError::Report(err.to_string()))?;
                }
                let buf = wtr
                    .into_inner()
                    .map_err(|e| BridgeError::Report(e.to_string()))?;
                writer
                    .write_all(&buf)
                    .map_err(|e| BridgeError::Report(e.to_string()))?;
            }
            AdmetReportFormat::Json => {
                let json = serde_json::to_string_pretty(evals)
                    .map_err(|e| BridgeError::Report(e.to_string()))?;
                writer
                    .write_all(json.as_bytes())
                    .map_err(|e| BridgeError::Report(e.to_string()))?;
            }
            AdmetReportFormat::Pdf => {
                let bytes = build_admet_summary_pdf(evals)?;
                writer
                    .write_all(&bytes)
                    .map_err(|e| BridgeError::Report(e.to_string()))?;
            }
            AdmetReportFormat::SvgBoiledEgg => {
                return Err(BridgeError::Report(
                    "BOILED-Egg SVG export requires external `admetcraft-cli` on PATH".to_string(),
                ));
            }
        }
        writer
            .flush()
            .map_err(|e| BridgeError::Report(e.to_string()))
    }
}

fn sanitize_csv(val: &str) -> String {
    if val
        .trim_start_matches(' ')
        .starts_with(['=', '+', '-', '@', '\t', '\r'])
    {
        format!("'{val}")
    } else {
        val.to_string()
    }
}

fn build_admet_summary_pdf(evals: &[MoleculeEvaluation]) -> Result<Vec<u8>, BridgeError> {
    let mut stream = String::new();
    let _ = writeln!(
        stream,
        "BT /F2 13 Tf 40 760 Td (ChemCraft ADMET Summary Report) Tj ET"
    );
    let _ = writeln!(
        stream,
        "BT /F1 9 Tf 40 742 Td (Compounds: {}) Tj ET",
        evals.len()
    );
    let mut y = 716.0_f64;
    for e in evals.iter().take(36) {
        let line = format!(
            "{} | {} | MW {:.1} | LogP {:.2} | TPSA {:.1} | GI {} | Lip {}",
            e.name,
            e.physicochemical.formula,
            e.physicochemical.molecular_weight,
            e.lipophilicity.consensus_log_p,
            e.physicochemical.tpsa,
            e.pharmacokinetics.gi_absorption.as_str(),
            e.druglikeness.lipinski.violations
        );
        let esc = line
            .replace('\\', "\\\\")
            .replace('(', "\\(")
            .replace(')', "\\)");
        let _ = writeln!(stream, "BT /F1 8.5 Tf 40 {y:.1} Td ({esc}) Tj ET");
        y -= 16.0;
        if y < 48.0 {
            break;
        }
    }

    let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
    enc.write_all(stream.as_bytes())
        .map_err(|e| BridgeError::Report(e.to_string()))?;
    let comp = enc
        .finish()
        .map_err(|e| BridgeError::Report(e.to_string()))?;

    let mut pdf = Vec::new();
    pdf.extend_from_slice(b"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");
    let mut off = vec![0usize];
    for obj in [
        b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n".as_slice(),
        b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n",
        b"3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 4 0 R /F2 5 0 R >> >> /Contents 6 0 R >>\nendobj\n",
        b"4 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>\nendobj\n",
        b"5 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>\nendobj\n",
    ] {
        off.push(pdf.len());
        pdf.extend_from_slice(obj);
    }
    off.push(pdf.len());
    pdf.extend_from_slice(
        format!(
            "6 0 obj\n<< /Length {} /Filter /FlateDecode >>\nstream\n",
            comp.len()
        )
        .as_bytes(),
    );
    pdf.extend_from_slice(&comp);
    pdf.extend_from_slice(b"\nendstream\nendobj\n");
    let xref = pdf.len();
    let _ = writeln!(&mut pdf, "xref\n0 {}\n0000000000 65535 f ", off.len());
    for &o in &off[1..] {
        let _ = writeln!(&mut pdf, "{o:010} 00000 n ");
    }
    let _ = writeln!(
        &mut pdf,
        "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF",
        off.len()
    );
    Ok(pdf)
}

/// Evaluate before/after ADMET impact of a scaffold mutation / bioisostere swap.
pub fn compute_admet_delta<P: AdmetProvider>(
    provider: &P,
    old_smiles: &str,
    new_smiles: &str,
    name: Option<&str>,
) -> Result<(AdmetDelta, MoleculeEvaluation), BridgeError> {
    let before = provider.evaluate_smiles(old_smiles, name)?;
    let after = provider.evaluate_smiles(new_smiles, name)?;
    let delta = AdmetDelta {
        old_smiles: before.canonical_smiles.clone(),
        new_smiles: after.canonical_smiles.clone(),
        delta_mw: after.physicochemical.molecular_weight - before.physicochemical.molecular_weight,
        delta_logp: after.lipophilicity.consensus_log_p - before.lipophilicity.consensus_log_p,
        delta_tpsa: after.physicochemical.tpsa - before.physicochemical.tpsa,
        delta_esol_log_s: after.water_solubility.esol.log_s - before.water_solubility.esol.log_s,
        old_gi: before.pharmacokinetics.gi_absorption.as_str().to_string(),
        new_gi: after.pharmacokinetics.gi_absorption.as_str().to_string(),
        old_bbb: before.pharmacokinetics.bbb_permeant,
        new_bbb: after.pharmacokinetics.bbb_permeant,
        old_lipinski_violations: before.druglikeness.lipinski.violations,
        new_lipinski_violations: after.druglikeness.lipinski.violations,
    };
    Ok((delta, after))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn admetcraft_provider_evaluates_aspirin_and_tetrazole_bioisostere_delta()
    -> Result<(), BridgeError> {
        let provider = AdmetCraftProvider::default();
        let (delta, after) = compute_admet_delta(
            &provider,
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(=O)Oc1ccccc1c1nnn[nH]1",
            Some("Tetrazole_Aspirin"),
        )?;
        assert!(delta.delta_tpsa > 15.0);
        assert!(after.physicochemical.tpsa > 70.0);
        Ok(())
    }
}
