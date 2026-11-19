//! Clean-room ChemDraw-style `1H` and `13C` NMR Chemical Shift & Multiplicity Predictor.
//!
//! Uses additive substituent increments (Pretsch / Shoolery / aromatic ortho-meta-para shielding rules)
//! directly on `MolGraph` topology, grouping magnetically/topologically equivalent nuclei cleanly
//! and computing `1H` J-coupling multiplicities (`s`, `d`, `t`, `q`, `m`, `br s`).

use std::collections::BTreeMap;

use mol_graph::{BondOrder, MolGraph};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NmrPeak {
    pub shift_ppm: f64,
    pub multiplicity: String,
    pub integration: usize,
    pub environment: String,
    pub atom_indices: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NmrSpectrum {
    pub h1_peaks: Vec<NmrPeak>,
    pub c13_peaks: Vec<NmrPeak>,
}

/// Predict `1H` and `13C` NMR chemical shift tables for `mol`.
pub fn predict_nmr(mol: &MolGraph) -> NmrSpectrum {
    let mut c13_raw: Vec<(f64, usize, String)> = Vec::new();
    let mut h1_raw: Vec<(f64, String, usize, usize, String)> = Vec::new();

    for (idx, a) in mol.atoms.iter().enumerate() {
        if a.symbol == "C" {
            let (c_shift, c_env) = predict_c13_atom(mol, idx);
            c13_raw.push((round1(c_shift), idx, c_env));
        }

        let h_cnt = usize::from(a.total_h());
        if h_cnt > 0 && a.atomic_number != 1 {
            let (h_shift, mult, h_env) = predict_h1_group(mol, idx);
            h1_raw.push((round2(h_shift), mult, h_cnt, idx, h_env));
        }
    }

    // Group equivalent 13C shifts within 0.05 ppm & same label
    let mut c13_groups: BTreeMap<(i64, String), Vec<usize>> = BTreeMap::new();
    for (ppm, idx, env) in c13_raw {
        let key = ((ppm * 10.0).round() as i64, env);
        c13_groups.entry(key).or_default().push(idx);
    }
    let mut c13_peaks: Vec<NmrPeak> = c13_groups
        .into_iter()
        .map(|((k, env), atoms)| NmrPeak {
            shift_ppm: (k as f64) / 10.0,
            multiplicity: "s".to_string(),
            integration: atoms.len(),
            environment: env,
            atom_indices: atoms,
        })
        .collect();
    c13_peaks.sort_by(|a, b| b.shift_ppm.total_cmp(&a.shift_ppm));

    // Group equivalent 1H shifts with same shift (0.02 ppm) & multiplicity & environment
    let mut h1_groups: BTreeMap<(i64, String, String), (usize, Vec<usize>)> = BTreeMap::new();
    for (ppm, mult, h_cnt, idx, env) in h1_raw {
        let key = ((ppm * 50.0).round() as i64, mult, env);
        let entry = h1_groups.entry(key).or_insert((0, Vec::new()));
        entry.0 += h_cnt;
        entry.1.push(idx);
    }
    let mut h1_peaks: Vec<NmrPeak> = h1_groups
        .into_iter()
        .map(|((k, mult, env), (integ, atoms))| NmrPeak {
            shift_ppm: round2((k as f64) / 50.0),
            multiplicity: mult,
            integration: integ,
            environment: env,
            atom_indices: atoms,
        })
        .collect();
    h1_peaks.sort_by(|a, b| b.shift_ppm.total_cmp(&a.shift_ppm));

    NmrSpectrum {
        h1_peaks,
        c13_peaks,
    }
}

fn round1(v: f64) -> f64 {
    (v * 10.0).round() / 10.0
}

fn round2(v: f64) -> f64 {
    (v * 100.0).round() / 100.0
}

fn has_double_to(mol: &MolGraph, idx: usize, sym: &str) -> bool {
    mol.adj[idx]
        .iter()
        .any(|&(v, bid)| mol.bonds[bid].order == BondOrder::Double && mol.atoms[v].symbol == sym)
}

fn has_triple_to(mol: &MolGraph, idx: usize, sym: &str) -> bool {
    mol.adj[idx]
        .iter()
        .any(|&(v, bid)| mol.bonds[bid].order == BondOrder::Triple && mol.atoms[v].symbol == sym)
}

fn is_carbonyl_c(mol: &MolGraph, idx: usize) -> bool {
    mol.atoms[idx].symbol == "C" && has_double_to(mol, idx, "O")
}

fn predict_c13_atom(mol: &MolGraph, idx: usize) -> (f64, String) {
    let a = &mol.atoms[idx];
    let nbrs = &mol.adj[idx];

    // 1. Carbonyl / Thiocarbonyl C=O, C=S
    if has_double_to(mol, idx, "O") {
        let o_single = nbrs.iter().any(|&(v, bid)| {
            mol.atoms[v].symbol == "O" && mol.bonds[bid].order == BondOrder::Single
        });
        let n_single = nbrs.iter().any(|&(v, bid)| {
            mol.atoms[v].symbol == "N" && mol.bonds[bid].order == BondOrder::Single
        });
        let h_cnt = a.total_h();
        if a.aromatic {
            return (155.0, "Aromatic C=O (pyridone/xanthine)".to_string());
        }
        if h_cnt >= 1 {
            return (195.0, "Aldehyde C=O".to_string());
        }
        if o_single {
            let is_acid = nbrs.iter().any(|&(v, bid)| {
                mol.atoms[v].symbol == "O"
                    && mol.bonds[bid].order == BondOrder::Single
                    && mol.atoms[v].total_h() > 0
            });
            if is_acid {
                return (176.5, "Carboxylic acid C=O".to_string());
            }
            return (170.5, "Ester C=O".to_string());
        }
        if n_single {
            return (168.5, "Amide C=O".to_string());
        }
        return (205.0, "Ketone C=O".to_string());
    }

    // 2. Nitrile / Alkyne sp C
    if has_triple_to(mol, idx, "N") {
        return (117.5, "Nitrile C#N".to_string());
    }
    if has_triple_to(mol, idx, "C") {
        let shift = if a.total_h() > 0 { 71.0 } else { 81.5 };
        return (shift, "Alkyne sp C".to_string());
    }

    // 3. Aromatic C (benzene base 128.5 ppm + ipso/ortho/meta/para increments)
    if a.aromatic {
        let mut ppm = 128.5_f64;
        // Ipso substituent attached to `idx`
        for &(v, bid) in nbrs {
            if mol.bonds[bid].order == BondOrder::Aromatic {
                // Ring heteroatom (pyridine/pyrrole/furan/etc.)
                match mol.atoms[v].symbol.as_str() {
                    "N" => ppm += 11.5,
                    "O" => ppm += 14.0,
                    "S" => ppm -= 2.5,
                    _ => {}
                }
            } else {
                ppm += ipso_c13_increment(mol, v);
            }
        }
        // Ortho substituents on aromatic neighbors
        for &(nb, bid) in nbrs {
            if mol.bonds[bid].order != BondOrder::Aromatic || !mol.atoms[nb].aromatic {
                continue;
            }
            for &(sub, b2) in &mol.adj[nb] {
                if sub == idx || mol.bonds[b2].order == BondOrder::Aromatic {
                    continue;
                }
                ppm += ortho_c13_increment(mol, sub);
            }
        }
        return (ppm.clamp(95.0, 172.0), "Aromatic C".to_string());
    }

    // 4. Alkene / Imine sp2 C
    if nbrs
        .iter()
        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Double)
    {
        let mut ppm = 123.0_f64;
        for &(v, bid) in nbrs {
            if mol.bonds[bid].order == BondOrder::Double && mol.atoms[v].symbol == "N" {
                return (160.0, "Imine C=N".to_string());
            }
            if mol.bonds[bid].order == BondOrder::Single {
                match mol.atoms[v].symbol.as_str() {
                    "C" => ppm += 10.5,
                    "O" => ppm += 28.0,
                    "N" => ppm += 16.0,
                    "Cl" => ppm += 3.0,
                    "F" => ppm += 25.0,
                    _ => {}
                }
            }
        }
        return (ppm.clamp(90.0, 165.0), "Alkene sp2 C".to_string());
    }

    // 5. Aliphatic sp3 C (Grant-Paul additive rule: base -2.3 + alpha/beta increments)
    let mut ppm = -2.3_f64;
    for &(v, _) in nbrs {
        match mol.atoms[v].symbol.as_str() {
            "C" => {
                if mol.atoms[v].aromatic {
                    ppm += 23.0;
                } else if is_carbonyl_c(mol, v) {
                    ppm += 30.0;
                } else {
                    ppm += 9.1;
                }
            }
            "O" => ppm += 51.0,
            "N" => ppm += 31.0,
            "S" => ppm += 17.0,
            "F" => ppm += 68.0,
            "Cl" => ppm += 31.0,
            "Br" => ppm += 20.0,
            "I" => ppm -= 6.0,
            _ => ppm += 8.0,
        }
        // Beta neighbors
        for &(w, _) in &mol.adj[v] {
            if w == idx {
                continue;
            }
            match mol.atoms[w].symbol.as_str() {
                "C" => ppm += 7.5,
                "O" | "N" => ppm += 5.5,
                "F" | "Cl" | "Br" => ppm += 6.0,
                _ => {}
            }
        }
    }
    let kind = match a.total_h() {
        3 => "Aliphatic CH3",
        2 => "Aliphatic CH2",
        1 => "Aliphatic CH",
        _ => "Quaternary sp3 C",
    };
    (ppm.clamp(4.0, 128.0), kind.to_string())
}

fn ipso_c13_increment(mol: &MolGraph, v: usize) -> f64 {
    match mol.atoms[v].symbol.as_str() {
        "C" => {
            if is_carbonyl_c(mol, v) {
                2.5
            } else if has_triple_to(mol, v, "N") {
                -15.5
            } else {
                9.2
            }
        }
        "O" => 29.0,
        "N" => {
            if mol.atoms[v].formal_charge > 0 {
                19.6 // -NO2
            } else {
                18.5 // -NH2 / -NR2
            }
        }
        "F" => 35.0,
        "Cl" => 6.2,
        "Br" => -5.5,
        "I" => -34.0,
        "S" => 10.0,
        _ => 4.0,
    }
}

fn ortho_c13_increment(mol: &MolGraph, v: usize) -> f64 {
    match mol.atoms[v].symbol.as_str() {
        "O" => -12.8,
        "N" => {
            if mol.atoms[v].formal_charge > 0 {
                -5.3
            } else {
                -13.3
            }
        }
        "F" => -13.0,
        "Cl" => 0.4,
        "Br" => 3.4,
        "C" => {
            if is_carbonyl_c(mol, v) {
                1.5
            } else {
                0.7
            }
        }
        _ => 0.0,
    }
}

fn predict_h1_group(mol: &MolGraph, idx: usize) -> (f64, String, String) {
    let a = &mol.atoms[idx];
    let nbrs = &mol.adj[idx];

    // Heteroatom exchangeable protons (OH, NH, SH)
    match a.symbol.as_str() {
        "O" => {
            let on_acid = nbrs.iter().any(|&(v, _)| is_carbonyl_c(mol, v));
            let on_arom = nbrs.iter().any(|&(v, _)| mol.atoms[v].aromatic);
            if on_acid {
                return (11.50, "br s".to_string(), "Carboxylic -COOH".to_string());
            }
            if on_arom {
                return (5.80, "br s".to_string(), "Phenolic -OH".to_string());
            }
            return (2.60, "br s".to_string(), "Alcohol -OH".to_string());
        }
        "N" => {
            if a.aromatic {
                return (
                    10.50,
                    "br s".to_string(),
                    "Aromatic [nH] (pyrrole/indole)".to_string(),
                );
            }
            let on_amide = nbrs
                .iter()
                .any(|&(v, _)| is_carbonyl_c(mol, v) || mol.atoms[v].symbol == "S");
            let on_arom = nbrs.iter().any(|&(v, _)| mol.atoms[v].aromatic);
            if on_amide {
                return (
                    7.80,
                    "br s".to_string(),
                    "Amide/Sulfonamide -NH".to_string(),
                );
            }
            if on_arom {
                return (4.10, "br s".to_string(), "Aniline -NH".to_string());
            }
            return (1.60, "br s".to_string(), "Amine -NH".to_string());
        }
        "S" => return (1.50, "br s".to_string(), "Thiol -SH".to_string()),
        _ => {}
    }

    // Compute vicinal 3J coupling partner count (H atoms on adjacent C atoms)
    let mut vicinal_h_groups: Vec<usize> = Vec::new();
    for &(v, bid) in nbrs {
        if mol.atoms[v].symbol == "C" && mol.bonds[bid].order != BondOrder::Triple {
            let vh = usize::from(mol.atoms[v].total_h());
            if vh > 0 {
                vicinal_h_groups.push(vh);
            }
        }
    }
    let mult = format_multiplicity(&vicinal_h_groups);

    // Aldehyde C(=O)H
    if has_double_to(mol, idx, "O") {
        return (9.75, mult, "Aldehyde -CHO".to_string());
    }
    // Terminal alkyne C#CH
    if has_triple_to(mol, idx, "C") {
        return (2.35, "s".to_string(), "Alkyne C#CH".to_string());
    }
    // Aromatic C-H (benzene base 7.26 + ortho/meta increments)
    if a.aromatic {
        let mut ppm = 7.26_f64;
        for &(nb, bid) in nbrs {
            if mol.bonds[bid].order == BondOrder::Aromatic {
                match mol.atoms[nb].symbol.as_str() {
                    "N" => ppm += 1.20,
                    "O" => ppm += 0.15,
                    "S" => ppm -= 0.05,
                    _ => {}
                }
                for &(sub, b2) in &mol.adj[nb] {
                    if sub == idx || mol.bonds[b2].order == BondOrder::Aromatic {
                        continue;
                    }
                    ppm += ortho_h1_increment(mol, sub);
                }
            }
        }
        return (ppm.clamp(5.80, 9.60), mult, "Aromatic C-H".to_string());
    }
    // Alkene C=C-H
    if nbrs
        .iter()
        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Double)
    {
        let mut ppm = 5.25_f64;
        for &(v, bid) in nbrs {
            if mol.bonds[bid].order == BondOrder::Single {
                match mol.atoms[v].symbol.as_str() {
                    "C" => {
                        if mol.atoms[v].aromatic {
                            ppm += 1.38;
                        } else if is_carbonyl_c(mol, v) {
                            ppm += 1.05;
                        } else {
                            ppm += 0.45;
                        }
                    }
                    "O" => ppm += 1.22,
                    "N" => ppm += 0.80,
                    "Cl" => ppm += 1.08,
                    _ => {}
                }
            }
        }
        return (ppm.clamp(4.40, 8.20), mult, "Vinyl/Alkene =C-H".to_string());
    }

    // Aliphatic C-H (Shoolery's rule: base CH3=0.86, CH2=1.25, CH=1.55 + alpha substituent deltas)
    let mut ppm: f64 = match a.total_h() {
        3 => 0.86,
        2 => 1.25,
        _ => 1.55,
    };
    for &(v, _) in nbrs {
        match mol.atoms[v].symbol.as_str() {
            "C" => {
                if mol.atoms[v].aromatic {
                    ppm += 1.45;
                } else if is_carbonyl_c(mol, v) {
                    ppm += 1.15;
                } else if has_triple_to(mol, v, "N") {
                    ppm += 1.30;
                } else if mol.adj[v]
                    .iter()
                    .any(|&(_, b2)| mol.bonds[b2].order == BondOrder::Double)
                {
                    ppm += 0.80;
                }
            }
            "O" => {
                let is_ester_or_aryl = mol.adj[v]
                    .iter()
                    .any(|&(w, _)| w != idx && (is_carbonyl_c(mol, w) || mol.atoms[w].aromatic));
                ppm += if is_ester_or_aryl { 2.85 } else { 2.35 };
            }
            "N" => {
                let is_amide = mol.adj[v]
                    .iter()
                    .any(|&(w, _)| w != idx && (is_carbonyl_c(mol, w) || mol.atoms[w].aromatic));
                ppm += if is_amide { 2.15 } else { 1.55 };
            }
            "S" => ppm += 1.35,
            "F" => ppm += 3.15,
            "Cl" => ppm += 2.35,
            "Br" => ppm += 2.20,
            "I" => ppm += 2.00,
            _ => {}
        }
    }
    let env = match a.total_h() {
        3 => "Aliphatic -CH3",
        2 => "Aliphatic -CH2-",
        _ => "Aliphatic >CH-",
    };
    (ppm.clamp(0.50, 7.60), mult, env.to_string())
}

fn ortho_h1_increment(mol: &MolGraph, v: usize) -> f64 {
    match mol.atoms[v].symbol.as_str() {
        "O" => -0.48,
        "N" => {
            if mol.atoms[v].formal_charge > 0 {
                0.95 // -NO2 strongly deshielding
            } else {
                -0.65 // -NH2 shielding
            }
        }
        "F" => -0.26,
        "Cl" => 0.03,
        "Br" => 0.18,
        "C" => {
            if is_carbonyl_c(mol, v) {
                0.65
            } else if has_triple_to(mol, v, "N") {
                0.36
            } else {
                -0.18
            }
        }
        _ => 0.0,
    }
}

fn format_multiplicity(vicinal_groups: &[usize]) -> String {
    let total: usize = vicinal_groups.iter().sum();
    if total == 0 {
        return "s".to_string();
    }
    if vicinal_groups.len() == 1 {
        return match total {
            1 => "d".to_string(),
            2 => "t".to_string(),
            3 => "q".to_string(),
            4 => "quin".to_string(),
            5 => "sext".to_string(),
            6 => "sept".to_string(),
            _ => "m".to_string(),
        };
    }
    if vicinal_groups.len() == 2 && vicinal_groups[0] == 1 && vicinal_groups[1] == 1 {
        return "dd".to_string();
    }
    if total <= 3 {
        match total {
            2 => "t".to_string(),
            3 => "q".to_string(),
            _ => "m".to_string(),
        }
    } else {
        "m".to_string()
    }
}
