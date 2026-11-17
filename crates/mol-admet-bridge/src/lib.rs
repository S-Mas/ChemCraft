//! Layer 2 (`mol-admet-bridge`): Decoupled `AdmetProvider` trait and direct in-process
//! connection to `AdmetCraft` (`admet-calc` and `admet-report`) for real-time SAR & canvas ADMET.

#![forbid(unsafe_code)]

use std::io::Write;

pub use admet_calc::{
    AbsorptionLevel, BioavailabilityRadar, BoiledEggPoint, Druglikeness, Lipophilicity,
    MedicinalChemistry, MoleculeEvaluation, Pharmacokinetics, PhysicochemicalProperties,
    SolubilityClass, WaterSolubility,
};
pub use admet_report::ExportFormat as AdmetReportFormat;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors produced by the `ChemCraft <-> AdmetCraft` bridge.
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

/// Decoupled trait for plugging in `AdmetCraft` (or mock/remote providers).
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

/// In-process provider linking directly to `../AdmetCraft/crates/admet-calc` & `admet-report`.
#[derive(Debug, Default, Clone, Copy)]
pub struct AdmetCraftProvider;

impl AdmetProvider for AdmetCraftProvider {
    fn evaluate_smiles(
        &self,
        smiles: &str,
        name: Option<&str>,
    ) -> Result<MoleculeEvaluation, BridgeError> {
        admet_calc::evaluate_smiles(smiles, name).map_err(|e| BridgeError::Calc(e.to_string()))
    }

    fn export_report_to_writer<W: Write>(
        &self,
        evals: &[MoleculeEvaluation],
        writer: &mut W,
        format: AdmetReportFormat,
    ) -> Result<(), BridgeError> {
        admet_report::write_report_to_writer(evals, writer, format)
            .map_err(|e| BridgeError::Report(e.to_string()))
    }
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
        let provider = AdmetCraftProvider;
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
