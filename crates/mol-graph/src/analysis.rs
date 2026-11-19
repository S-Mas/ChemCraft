//! ChemDraw-equivalent Molecular Analysis Window: exact monoisotopic mass, IUPAC average MW,
//! charge-annotated Hill formula, elemental analysis (%C, %H, ...), Double Bond Equivalents (DBE),
//! and ESI/APCI Mass Spectrometry m/z adduct table.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::element::element_by_symbol;
use crate::mol::MolGraph;

const ELECTRON_MASS: f64 = 0.000_548_579_9;
const PROTON_MASS: f64 = 1.007_276_466_6;
const H_MONO: f64 = 1.007_825_032_2;
const NA_ION: f64 = 22.989_220_702;
const K_ION: f64 = 38.963_157_906;
const NH4_ION: f64 = 18.033_825_554;
const CL_ION: f64 = 34.969_401_262;
const FORMATE_ION: f64 = 44.998_202_852;

/// Monoisotopic exact mass (most abundant natural IUPAC isotope) in unified atomic mass units (Da).
pub fn monoisotopic_mass(symbol: &str) -> f64 {
    match symbol {
        "H" => H_MONO,
        "D" => 2.014_101_778,
        "T" => 3.016_049_281,
        "He" => 4.002_603_254,
        "Li" => 7.016_003_437,
        "Be" => 9.012_183_07,
        "B" => 11.009_305_36,
        "C" => 12.000_000_0,
        "N" => 14.003_074_004,
        "O" => 15.994_914_620,
        "F" => 18.998_403_163,
        "Na" => 22.989_769_282,
        "Mg" => 23.985_041_697,
        "Al" => 26.981_538_53,
        "Si" => 27.976_926_535,
        "P" => 30.973_761_998,
        "S" => 31.972_071_174,
        "Cl" => 34.968_852_682,
        "K" => 38.963_706_486,
        "Ca" => 39.962_590_863,
        "Ti" => 47.947_941_98,
        "Cr" => 51.940_506_23,
        "Mn" => 54.938_043_91,
        "Fe" => 55.934_936_33,
        "Co" => 58.933_194_29,
        "Ni" => 57.935_342_41,
        "Cu" => 62.929_597_72,
        "Zn" => 63.929_142_01,
        "Ga" => 68.925_573_5,
        "Ge" => 73.921_177_76,
        "As" => 74.921_594_57,
        "Se" => 79.916_521_8,
        "Br" => 78.918_337_6,
        "Ru" => 101.904_344_1,
        "Rh" => 102.905_498_0,
        "Pd" => 105.903_480_4,
        "Ag" => 106.905_091_6,
        "Sn" => 119.902_201_6,
        "Sb" => 120.903_812_0,
        "Te" => 129.906_222_7,
        "I" => 126.904_471_9,
        "Pt" => 194.964_791_7,
        "Au" => 196.966_568_7,
        "Hg" => 201.970_643_4,
        "Pb" => 207.976_652_5,
        "Bi" => 208.980_399_1,
        other => element_by_symbol(other).map_or(12.0, |e| e.atomic_weight),
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MassAdduct {
    pub label: String,
    pub mode: String,
    pub mz: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElementPercent {
    pub symbol: String,
    pub count: usize,
    pub percent: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MolecularAnalysis {
    pub formula: String,
    pub formula_with_charge: String,
    pub net_charge: i32,
    pub exact_mass: f64,
    pub molecular_weight: f64,
    pub dbe: f64,
    pub heavy_atoms: usize,
    pub total_atoms: usize,
    pub ring_count: usize,
    pub chiral_centers: usize,
    pub elemental_analysis: Vec<ElementPercent>,
    pub elemental_analysis_string: String,
    pub adducts: Vec<MassAdduct>,
}

/// Compute full ChemDraw Analysis Window metrics for `mol`.
pub fn analyze_molecule(mol: &MolGraph) -> MolecularAnalysis {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut net_charge: i32 = 0;
    let mut exact_mass_neutral_sum = 0.0_f64;
    let mut avg_mw = 0.0_f64;

    for a in &mol.atoms {
        *counts.entry(a.symbol.clone()).or_insert(0) += 1;
        net_charge += i32::from(a.formal_charge);
        let atom_exact = if let Some(iso) = a.isotope {
            f64::from(iso)
        } else {
            monoisotopic_mass(&a.symbol)
        };
        exact_mass_neutral_sum += atom_exact;
        let atom_avg = element_by_symbol(&a.symbol).map_or(12.0, |e| e.atomic_weight);
        avg_mw += atom_avg;

        let h_imp = usize::from(a.total_h());
        if h_imp > 0 {
            *counts.entry("H".to_string()).or_insert(0) += h_imp;
            exact_mass_neutral_sum += (h_imp as f64) * H_MONO;
            avg_mw += (h_imp as f64) * 1.008;
        }
    }

    // Electron mass correction for formally charged ions (matches RDKit / ChemDraw Exact Mass)
    let exact_mass = exact_mass_neutral_sum - (net_charge as f64) * ELECTRON_MASS;
    let formula = mol.formula();
    let formula_with_charge = if net_charge == 0 {
        formula.clone()
    } else if net_charge == 1 {
        format!("{formula}+")
    } else if net_charge == -1 {
        format!("{formula}-")
    } else if net_charge > 1 {
        format!("{formula}+{net_charge}")
    } else {
        format!("{formula}{net_charge}")
    };

    // Hill ordering: C, H, then alphabetical
    let mut ordered_syms: Vec<String> = Vec::with_capacity(counts.len());
    if counts.contains_key("C") {
        ordered_syms.push("C".to_string());
        if counts.contains_key("H") {
            ordered_syms.push("H".to_string());
        }
        for k in counts.keys() {
            if k != "C" && k != "H" {
                ordered_syms.push(k.clone());
            }
        }
    } else {
        ordered_syms.extend(counts.keys().cloned());
    }

    let mut elemental_analysis = Vec::with_capacity(ordered_syms.len());
    let mut ea_parts = Vec::with_capacity(ordered_syms.len());
    for sym in &ordered_syms {
        let cnt = counts.get(sym).copied().unwrap_or(0);
        let aw = if sym == "H" {
            1.008
        } else {
            element_by_symbol(sym).map_or(12.0, |e| e.atomic_weight)
        };
        let pct = if avg_mw > 1e-9 {
            100.0 * (cnt as f64) * aw / avg_mw
        } else {
            0.0
        };
        ea_parts.push(format!("{sym}, {pct:.2}"));
        elemental_analysis.push(ElementPercent {
            symbol: sym.clone(),
            count: cnt,
            percent: (pct * 100.0).round() / 100.0,
        });
    }
    let elemental_analysis_string = ea_parts.join("; ");

    // Double Bond Equivalents (Rings + Double Bonds): C + Si - (H + X)/2 + (N + P)/2 + 1
    let c_cnt = counts.get("C").copied().unwrap_or(0) as f64
        + counts.get("Si").copied().unwrap_or(0) as f64;
    let h_cnt = counts.get("H").copied().unwrap_or(0) as f64;
    let x_cnt = (counts.get("F").copied().unwrap_or(0)
        + counts.get("Cl").copied().unwrap_or(0)
        + counts.get("Br").copied().unwrap_or(0)
        + counts.get("I").copied().unwrap_or(0)) as f64;
    let n_cnt =
        (counts.get("N").copied().unwrap_or(0) + counts.get("P").copied().unwrap_or(0)) as f64;
    let dbe = if mol.atoms.is_empty() {
        0.0
    } else {
        c_cnt - 0.5 * (h_cnt + x_cnt) + 0.5 * n_cnt + 1.0
    };

    let m = exact_mass;
    let round4 = |v: f64| (v * 10_000.0).round() / 10_000.0;
    let adducts = if mol.atoms.is_empty() {
        Vec::new()
    } else {
        vec![
            MassAdduct {
                label: "[M+H]+".to_string(),
                mode: "ESI+".to_string(),
                mz: round4(m + PROTON_MASS),
            },
            MassAdduct {
                label: "[M+NH4]+".to_string(),
                mode: "ESI+".to_string(),
                mz: round4(m + NH4_ION),
            },
            MassAdduct {
                label: "[M+Na]+".to_string(),
                mode: "ESI+".to_string(),
                mz: round4(m + NA_ION),
            },
            MassAdduct {
                label: "[M+K]+".to_string(),
                mode: "ESI+".to_string(),
                mz: round4(m + K_ION),
            },
            MassAdduct {
                label: "[M+2H]2+".to_string(),
                mode: "ESI+".to_string(),
                mz: round4((m + 2.0 * PROTON_MASS) * 0.5),
            },
            MassAdduct {
                label: "[2M+H]+".to_string(),
                mode: "ESI+".to_string(),
                mz: round4(2.0 * m + PROTON_MASS),
            },
            MassAdduct {
                label: "[M-H]-".to_string(),
                mode: "ESI-".to_string(),
                mz: round4((m - PROTON_MASS).max(0.0)),
            },
            MassAdduct {
                label: "[M+Cl]-".to_string(),
                mode: "ESI-".to_string(),
                mz: round4(m + CL_ION),
            },
            MassAdduct {
                label: "[M+HCOO]-".to_string(),
                mode: "ESI-".to_string(),
                mz: round4(m + FORMATE_ION),
            },
        ]
    };

    let total_atoms: usize = counts.values().sum();
    let chiral_centers = mol.atoms.iter().filter(|a| a.cip_label.is_some()).count();

    MolecularAnalysis {
        formula,
        formula_with_charge,
        net_charge,
        exact_mass: round4(exact_mass),
        molecular_weight: (avg_mw * 100.0).round() / 100.0,
        dbe: (dbe * 10.0).round() / 10.0,
        heavy_atoms: mol.heavy_atom_count(),
        total_atoms,
        ring_count: mol.rings.len(),
        chiral_centers,
        elemental_analysis,
        elemental_analysis_string,
        adducts,
    }
}
