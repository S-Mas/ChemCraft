use std::collections::BTreeMap;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use mol_admet_bridge::{AdmetProvider, AdmetSummary, MoleculeEvaluation};
use mol_graph::{
    MolError, MolGraph, SmartsQuery, parse_csv_records, parse_dwar, parse_sdf, parse_smiles,
};

use crate::fp::MorganFp;

/// Single compound row in the in-memory SAR Explorer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SarRow {
    pub id: usize,
    pub name: String,
    pub canonical_smiles: String,
    pub formula: String,
    pub mw: f64,
    pub heavy_atoms: usize,
    pub mol: MolGraph,
    pub fp: MorganFp,
    pub pic50: Option<f64>,
    pub series: String,
    pub properties: BTreeMap<String, String>,
    pub admet: Option<AdmetSummary>,
    pub admet_full: Option<MoleculeEvaluation>,
}

impl SarRow {
    /// Ligand Efficiency `LE = 1.37 * pIC50 / heavy_atoms`.
    pub fn ligand_efficiency(&self) -> Option<f64> {
        let p = self.pic50?;
        if self.heavy_atoms == 0 {
            None
        } else {
            Some(1.37 * p / (self.heavy_atoms as f64))
        }
    }

    /// Lipophilic Ligand Efficiency `LLE = pIC50 - cLogP`.
    pub fn lipophilic_ligand_efficiency(&self) -> Option<f64> {
        let p = self.pic50?;
        let logp = self.admet.as_ref()?.consensus_logp;
        Some(p - logp)
    }
}

/// Interactive DataWarrior-style multi-criteria filter configuration.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SarFilter {
    pub smarts_substructure: Option<String>,
    pub reference_smiles: Option<String>,
    pub min_tanimoto: f64,
    pub min_pic50: Option<f64>,
    pub max_pic50: Option<f64>,
    pub max_mw: Option<f64>,
    pub max_logp: Option<f64>,
    pub lipinski_pass_only: bool,
    pub no_pains_only: bool,
}

/// In-memory SAR dataset with fast ECFP4 Tanimoto, SMARTS, and ADMET filtering.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SarTable {
    pub rows: Vec<SarRow>,
}

impl SarTable {
    pub fn new() -> Self {
        Self { rows: Vec::new() }
    }

    /// Build a row cleanly from a `MolGraph`, extracting `pIC50` or `IC50_nM` if present in `mol.properties`.
    pub fn add_molecule(
        &mut self,
        mol: MolGraph,
        pic50_override: Option<f64>,
        series: Option<String>,
    ) -> usize {
        let id = self.rows.len();
        let name = mol
            .name
            .clone()
            .unwrap_or_else(|| format!("Compound_{}", id + 1));
        let canonical_smiles = mol.canonical_smiles();
        let formula = mol.formula();
        let mw = mol.molecular_weight();
        let heavy_atoms = mol.heavy_atom_count();
        let fp = MorganFp::ecfp4(&mol);

        let parsed_pic50 = pic50_override.or_else(|| extract_pic50(&mol.properties));
        let ser = series
            .or_else(|| mol.properties.get("Series").cloned())
            .unwrap_or_else(|| "Primary".to_string());

        let properties = mol.properties.clone();
        self.rows.push(SarRow {
            id,
            name,
            canonical_smiles,
            formula,
            mw,
            heavy_atoms,
            mol,
            fp,
            pic50: parsed_pic50,
            series: ser,
            properties,
            admet: None,
            admet_full: None,
        });
        id
    }

    /// Import compounds from SDF (V2000/V3000) text.
    pub fn load_sdf_text(&mut self, text: &str) -> Result<usize, MolError> {
        let parsed = parse_sdf(text);
        let mut added = 0usize;
        for mol in parsed.into_iter().flatten() {
            self.add_molecule(mol, None, None);
            added += 1;
        }
        if added == 0 {
            return Err(MolError::EmptyInput);
        }
        Ok(added)
    }

    /// Import compounds from `.dwar` DataWarrior text.
    pub fn load_dwar_text(&mut self, text: &str) -> Result<usize, MolError> {
        let parsed = parse_dwar(text)?;
        let mut added = 0usize;
        for item in parsed.into_iter().flatten() {
            self.add_molecule(item.mol, None, None);
            added += 1;
        }
        if added == 0 {
            return Err(MolError::EmptyInput);
        }
        Ok(added)
    }

    /// Import compounds from CSV text.
    pub fn load_csv_text(&mut self, text: &str) -> Result<usize, MolError> {
        let parsed = parse_csv_records(text)?;
        let mut added = 0usize;
        for item in parsed.into_iter().flatten() {
            self.add_molecule(item.mol, None, None);
            added += 1;
        }
        if added == 0 {
            return Err(MolError::EmptyInput);
        }
        Ok(added)
    }

    /// Populate all rows in parallel using any `AdmetProvider` (`AdmetCraftProvider`).
    pub fn enrich_with_admet<P: AdmetProvider>(&mut self, provider: &P) {
        let inputs: Vec<(String, Option<String>)> = self
            .rows
            .iter()
            .map(|r| (r.canonical_smiles.clone(), Some(r.name.clone())))
            .collect();
        let results = provider.evaluate_batch(&inputs);
        for (row, res) in self.rows.iter_mut().zip(results) {
            if let Ok(ev) = res {
                row.admet = Some(AdmetSummary::from(&ev));
                row.admet_full = Some(ev);
            }
        }
    }

    /// Compute indices of rows matching `filter` in parallel.
    pub fn matching_indices(&self, filter: &SarFilter) -> Result<Vec<usize>, MolError> {
        let smarts = match filter.smarts_substructure.as_deref().map(str::trim) {
            Some(pat) if !pat.is_empty() => Some(SmartsQuery::parse(pat)?),
            _ => None,
        };
        let ref_fp = match filter.reference_smiles.as_deref().map(str::trim) {
            Some(smi) if !smi.is_empty() && filter.min_tanimoto > 0.0 => {
                let m = parse_smiles(smi)?;
                Some(MorganFp::ecfp4(&m))
            }
            _ => None,
        };

        let indices: Vec<usize> = self
            .rows
            .par_iter()
            .enumerate()
            .filter_map(|(idx, row)| {
                if let Some(ref q) = smarts
                    && !q.matches(&row.mol)
                {
                    return None;
                }
                if let Some(ref rfp) = ref_fp
                    && row.fp.tanimoto(rfp) < filter.min_tanimoto
                {
                    return None;
                }
                if let Some(min_p) = filter.min_pic50
                    && row.pic50.is_none_or(|p| p < min_p)
                {
                    return None;
                }
                if let Some(max_p) = filter.max_pic50
                    && row.pic50.is_none_or(|p| p > max_p)
                {
                    return None;
                }
                if let Some(max_mw) = filter.max_mw
                    && row.mw > max_mw
                {
                    return None;
                }
                if let Some(max_lp) = filter.max_logp
                    && row
                        .admet
                        .as_ref()
                        .is_some_and(|a| a.consensus_logp > max_lp)
                {
                    return None;
                }
                if filter.lipinski_pass_only
                    && row.admet.as_ref().is_none_or(|a| a.lipinski_violations > 0)
                {
                    return None;
                }
                if filter.no_pains_only && row.admet.as_ref().is_none_or(|a| a.pains_count > 0) {
                    return None;
                }
                Some(idx)
            })
            .collect();

        Ok(indices)
    }

    /// Generate a realistic 50-compound medicinal chemistry SAR series (COX / Kinase / NSAID
    /// carboxylic acid, tetrazole bioisostere, sulfonamide, and heterocyclic analogues).
    pub fn demo_50_compounds() -> Self {
        let entries: [(&str, &str, f64, &str); 50] = [
            (
                "Ibuprofen",
                "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
                6.45,
                "Arylpropionic",
            ),
            (
                "Ibuprofen-Tetrazole",
                "CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1",
                7.12,
                "Bioisostere",
            ),
            (
                "Naproxen",
                "COc1ccc2cc(C(C)C(=O)O)ccc2c1",
                7.05,
                "Arylpropionic",
            ),
            (
                "Flurbiprofen",
                "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1",
                7.82,
                "Arylpropionic",
            ),
            (
                "Ketoprofen",
                "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1",
                7.18,
                "Arylpropionic",
            ),
            (
                "Diclofenac",
                "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
                8.15,
                "Arylacetic",
            ),
            (
                "Indomethacin",
                "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1",
                8.32,
                "Indole-Acetic",
            ),
            ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O", 5.20, "Salicylate"),
            ("Salicylic_Acid", "O=C(O)c1ccccc1O", 4.95, "Salicylate"),
            (
                "Diflunisal",
                "O=C(O)c1cc(-c2ccc(F)cc2F)ccc1O",
                6.80,
                "Salicylate",
            ),
            (
                "Celecoxib",
                "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
                8.55,
                "Coxib",
            ),
            (
                "Rofecoxib",
                "CS(=O)(=O)c1ccc(-c2oc(=O)cc2-c2ccccc2)cc1",
                8.10,
                "Coxib",
            ),
            (
                "Valdecoxib",
                "Cc1onc(-c2ccccc2)c1-c1ccc(S(N)(=O)=O)cc1",
                8.40,
                "Coxib",
            ),
            (
                "Etoricoxib",
                "Cc1ccc(-c2nc(-c3ccc(S(C)(=O)=O)cc3)ccc2Cl)cn1",
                8.62,
                "Coxib",
            ),
            (
                "Meloxicam",
                "Cc1cnc(NC(=O)C2=C(O)c3ccccc3S(=O)(=O)N2C)s1",
                7.90,
                "Oxicam",
            ),
            (
                "Piroxicam",
                "CN1C(=C(O)Nc2ccccn2)C(=O)c2ccccc2S1(=O)=O",
                7.45,
                "Oxicam",
            ),
            (
                "Mefenamic_Acid",
                "Cc1cccc(Nc2ccccc2C(=O)O)c1C",
                6.92,
                "Fenamate",
            ),
            (
                "Meclofenamic_Acid",
                "Cc1ccc(Cl)c(Nc2ccccc2C(=O)O)c1Cl",
                7.65,
                "Fenamate",
            ),
            (
                "Tolfenamic_Acid",
                "Cc1cccc(Nc2ccccc2C(=O)O)c1Cl",
                7.38,
                "Fenamate",
            ),
            (
                "Flufenamic_Acid",
                "O=C(O)c1ccccc1Nc1cccc(C(F)(F)F)c1",
                7.15,
                "Fenamate",
            ),
            (
                "Lumiracoxib",
                "Cc1ccc(F)c(Cl)c1Nc1ccccc1CC(=O)O",
                8.05,
                "Arylacetic",
            ),
            (
                "Fenbufen",
                "O=C(O)CCC(=O)c1ccc(-c2ccccc2)cc1",
                6.30,
                "Biphenyl",
            ),
            (
                "Sulindac",
                "CC1=C(CC(=O)O)c2cc(F)ccc2C1=Cc1ccc(S(C)=O)cc1",
                7.22,
                "Indene",
            ),
            (
                "Etodolac",
                "CCc1cccc2c1[nH]c1c2CCOC1(CC)CC(=O)O",
                7.08,
                "Indole-Acetic",
            ),
            (
                "Tolmetin",
                "Cc1ccc(C(=O)c2ccc(CC(=O)O)n2C)cc1",
                6.74,
                "Pyrrole-Acetic",
            ),
            (
                "Ketorolac",
                "O=C(c1ccccc1)c1ccc2n1CCC2C(=O)O",
                7.95,
                "Pyrrole-Acetic",
            ),
            (
                "Aceclofenac",
                "O=C(O)COC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
                7.50,
                "Arylacetic",
            ),
            (
                "Bromfenac",
                "Nc1c(CC(=O)O)cccc1C(=O)c1ccc(Br)cc1",
                8.20,
                "Arylacetic",
            ),
            (
                "Nepafenac",
                "NC(=O)Cc1cccc(C(=O)c2ccccc2)c1N",
                6.60,
                "Arylacetic",
            ),
            (
                "Oxaprozin",
                "O=C(O)CCc1nc(-c2ccccc2)c(-c2ccccc2)o1",
                6.48,
                "Oxazole",
            ),
            (
                "Loxoprofen",
                "CC(C(=O)O)c1ccc(CC2CCCC2=O)cc1",
                6.88,
                "Arylpropionic",
            ),
            (
                "Zaltoprofen",
                "CC(C(=O)O)c1ccc2c(c1)Cc1ccccc1C2=O",
                7.30,
                "Arylpropionic",
            ),
            (
                "Carprofen",
                "CC(C(=O)O)c1ccc2[nH]c3ccc(Cl)cc3c2c1",
                7.62,
                "Carbazole",
            ),
            (
                "Vedaprofen",
                "CC(C(=O)O)c1ccc2c(C3CCCCC3)cccc2c1",
                7.40,
                "Naphthyl",
            ),
            (
                "Tiaprofenic_Acid",
                "CC(C(=O)O)c1ccc(C(=O)c2ccccc2)s1",
                7.10,
                "Thiophene",
            ),
            (
                "Suprofen",
                "CC(C(=O)O)c1ccc(C(=O)c2cccs2)cc1",
                6.95,
                "Thiophene",
            ),
            (
                "Tenoxicam",
                "CN1C(=C(O)Nc2ccccn2)C(=O)c2sccc2S1(=O)=O",
                7.25,
                "Oxicam",
            ),
            (
                "Lornoxicam",
                "CN1C(=C(O)Nc2ccccn2)C(=O)c2sc(Cl)cc2S1(=O)=O",
                7.92,
                "Oxicam",
            ),
            (
                "Nabumetone",
                "COc1ccc2cc(CCC(C)=O)ccc2c1",
                6.15,
                "Arylpropionic",
            ),
            (
                "Fenoprofen",
                "CC(C(=O)O)c1cccc(Oc2ccccc2)c1",
                6.65,
                "Arylpropionic",
            ),
            ("Phenacetin", "CCOc1ccc(NC(C)=O)cc1", 4.65, "Anilide"),
            ("Paracetamol", "CC(=O)Nc1ccc(O)cc1", 4.82, "Anilide"),
            (
                "Flunixin",
                "Cc1c(Nc2cccc(C(F)(F)F)c2)ncccc1C(=O)O",
                7.55,
                "Pyridine",
            ),
            (
                "Niflumic_Acid",
                "O=C(O)c1cccnc1Nc1cccc(C(F)(F)F)c1",
                7.12,
                "Pyridine",
            ),
            ("Clonixin", "Cc1cccc(Nc2cccnc2C(=O)O)c1Cl", 7.28, "Pyridine"),
            (
                "Valsartan_Core",
                "CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1)C(C(C)C)C(=O)O",
                8.70,
                "Bioisostere",
            ),
            (
                "Losartan_Core",
                "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1",
                8.45,
                "Bioisostere",
            ),
            (
                "Candesartan_Core",
                "CCOc1nc2cccc(C(=O)O)c2n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1",
                8.85,
                "Bioisostere",
            ),
            (
                "Nimesulide",
                "CS(=O)(=O)Nc1ccc([N+](=O)[O-])cc1Oc1ccccc1",
                7.35,
                "Sulfonanilide",
            ),
            (
                "NS398",
                "CS(=O)(=O)Nc1ccc([N+](=O)[O-])cc1OC1CCCCC1",
                7.75,
                "Sulfonanilide",
            ),
        ];

        let mut table = Self::new();
        for (name, smi, pic50, series) in entries {
            if let Ok(mut mol) = parse_smiles(smi) {
                mol.name = Some(name.to_string());
                mol.properties
                    .insert("pIC50".to_string(), format!("{pic50:.2}"));
                mol.properties
                    .insert("Series".to_string(), series.to_string());
                table.add_molecule(mol, Some(pic50), Some(series.to_string()));
            }
        }
        table
    }
}

fn extract_pic50(props: &BTreeMap<String, String>) -> Option<f64> {
    for (k, v) in props {
        let kl = k.to_ascii_lowercase();
        if kl == "pic50" || kl == "pki" || kl == "pic_50" {
            if let Ok(val) = v.trim().parse::<f64>()
                && val.is_finite()
            {
                return Some(val);
            }
        } else if (kl == "ic50" || kl == "ic50_nm" || kl == "ic50 (nm)")
            && let Ok(nm) = v.trim().parse::<f64>()
            && nm.is_finite()
            && nm > 0.0
        {
            return Some(9.0 - nm.log10());
        }
    }
    None
}
