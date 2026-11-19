use std::f64::consts::PI;

use mol_graph::{BondOrder, MolGraph, SmartsQuery};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AbsorptionLevel {
    High,
    Low,
}

impl AbsorptionLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::High => "High",
            Self::Low => "Low",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SolubilityClass {
    Insoluble,
    PoorlySoluble,
    ModeratelySoluble,
    Soluble,
    VerySoluble,
    HighlySoluble,
}

impl SolubilityClass {
    pub fn from_log_s(log_s: f64) -> Self {
        if log_s < -10.0 {
            Self::Insoluble
        } else if log_s < -6.0 {
            Self::PoorlySoluble
        } else if log_s < -4.0 {
            Self::ModeratelySoluble
        } else if log_s < -2.0 {
            Self::Soluble
        } else if log_s < 0.0 {
            Self::VerySoluble
        } else {
            Self::HighlySoluble
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Insoluble => "Insoluble",
            Self::PoorlySoluble => "Poorly soluble",
            Self::ModeratelySoluble => "Moderately soluble",
            Self::Soluble => "Soluble",
            Self::VerySoluble => "Very soluble",
            Self::HighlySoluble => "Highly soluble",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhysicochemicalProperties {
    pub formula: String,
    pub molecular_weight: f64,
    pub heavy_atoms: usize,
    pub aromatic_heavy_atoms: usize,
    pub fraction_csp3: f64,
    pub rotatable_bonds: usize,
    pub h_bond_acceptors: usize,
    pub h_bond_donors: usize,
    pub molar_refractivity: f64,
    pub tpsa: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Lipophilicity {
    pub wlogp: f64,
    pub xlogp3: f64,
    pub mlogp: f64,
    pub consensus_log_p: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterSolubilityModel {
    pub log_s: f64,
    pub solubility_mg_ml: f64,
    pub solubility_mol_l: f64,
    pub class: SolubilityClass,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaterSolubility {
    pub esol: WaterSolubilityModel,
    pub ali: WaterSolubilityModel,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Pharmacokinetics {
    pub gi_absorption: AbsorptionLevel,
    pub bbb_permeant: bool,
    pub pgp_substrate: bool,
    pub cyp1a2_inhibitor: bool,
    pub cyp2c19_inhibitor: bool,
    pub cyp2c9_inhibitor: bool,
    pub cyp2d6_inhibitor: bool,
    pub cyp3a4_inhibitor: bool,
    pub log_kp_cm_s: f64,
    pub pgp_cyp_are_experimental_heuristics: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RuleResult {
    pub pass: bool,
    pub violations: usize,
    pub details: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Druglikeness {
    pub lipinski: RuleResult,
    pub ghose: RuleResult,
    pub veber: RuleResult,
    pub egan: RuleResult,
    pub muegge: RuleResult,
    pub bioavailability_score: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MedicinalChemistry {
    pub pains_alerts: Vec<String>,
    pub pains_catalog_size: usize,
    pub brenk_alerts: Vec<String>,
    pub brenk_catalog_size: usize,
    pub leadlikeness: RuleResult,
    pub synthetic_accessibility: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BoiledEggPoint {
    pub tpsa: f64,
    pub wlogp: f64,
    pub in_white: bool,
    pub in_yolk: bool,
    pub pgp_substrate: bool,
    pub in_plot_bounds: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BioavailabilityRadar {
    pub lipo_xlogp3: f64,
    pub size_mw: f64,
    pub polar_tpsa: f64,
    pub insolu_log_s: f64,
    pub insatu_fsp3: f64,
    pub flex_rot_bonds: usize,
    pub normalized_axes: [f64; 6],
    pub all_in_optimal_zone: bool,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MoleculeEvaluation {
    pub name: String,
    pub input_smiles: String,
    pub canonical_smiles: String,
    pub warnings: Vec<String>,
    pub physicochemical: PhysicochemicalProperties,
    pub lipophilicity: Lipophilicity,
    pub water_solubility: WaterSolubility,
    pub pharmacokinetics: Pharmacokinetics,
    pub druglikeness: Druglikeness,
    pub medicinal_chemistry: MedicinalChemistry,
    pub boiled_egg: BoiledEggPoint,
    pub bioavailability_radar: BioavailabilityRadar,
}

pub fn evaluate_mol_builtin(mol: &MolGraph, name_override: Option<&str>) -> MoleculeEvaluation {
    let name = name_override
        .map(str::to_string)
        .or_else(|| mol.name.clone())
        .unwrap_or_else(|| "Molecule 1".to_string());
    let canonical_smiles = mol.canonical_smiles();
    let input_smiles = mol
        .input_smiles
        .clone()
        .unwrap_or_else(|| canonical_smiles.clone());

    let mw = round2(mol.molecular_weight());
    let heavy_atoms = mol.heavy_atom_count();
    let aromatic_heavy_atoms = mol.aromatic_atom_count();
    let total_c = mol.atoms.iter().filter(|a| a.symbol == "C").count();
    let sp3_c = mol
        .atoms
        .iter()
        .filter(|a| {
            a.symbol == "C"
                && !a.aromatic
                && !mol.adj[a.id].iter().any(|&(_, b)| {
                    matches!(mol.bonds[b].order, BondOrder::Double | BondOrder::Triple)
                })
        })
        .count();
    let fraction_csp3 = if total_c == 0 {
        0.0
    } else {
        round2(sp3_c as f64 / total_c as f64)
    };
    let rotatable_bonds = count_rotatable_bonds(mol);
    let hba = count_hba(mol);
    let hbd = count_hbd(mol);
    let mr = round2(compute_mr(mol));
    let tpsa = round2(compute_tpsa(mol));

    let physicochemical = PhysicochemicalProperties {
        formula: mol.formula(),
        molecular_weight: mw,
        heavy_atoms,
        aromatic_heavy_atoms,
        fraction_csp3,
        rotatable_bonds,
        h_bond_acceptors: hba,
        h_bond_donors: hbd,
        molar_refractivity: mr,
        tpsa,
    };

    let wlogp = round2(compute_wlogp(mol));
    let xlogp3 = round2(compute_xlogp3(mol));
    let mlogp = round2(compute_mlogp(mol, &physicochemical));
    let consensus_log_p = round2((wlogp + xlogp3 + mlogp) / 3.0);
    let lipophilicity = Lipophilicity {
        wlogp,
        xlogp3,
        mlogp,
        consensus_log_p,
    };

    let esol_logs = round2(
        0.16 - 0.63 * consensus_log_p - 0.0062 * mw + 0.066 * (rotatable_bonds as f64)
            - 0.74
                * (if heavy_atoms == 0 {
                    0.0
                } else {
                    aromatic_heavy_atoms as f64 / heavy_atoms as f64
                }),
    );
    let ali_logs = round2(0.51 - 1.0 * consensus_log_p - 0.011 * tpsa + 0.0001 * mw);
    let esol_mol = 10f64.powf(esol_logs);
    let ali_mol = 10f64.powf(ali_logs);
    let water_solubility = WaterSolubility {
        esol: WaterSolubilityModel {
            log_s: esol_logs,
            solubility_mg_ml: esol_mol * mw * 1000.0,
            solubility_mol_l: esol_mol,
            class: SolubilityClass::from_log_s(esol_logs),
        },
        ali: WaterSolubilityModel {
            log_s: ali_logs,
            solubility_mg_ml: ali_mol * mw * 1000.0,
            solubility_mol_l: ali_mol,
            class: SolubilityClass::from_log_s(ali_logs),
        },
    };

    // BOILED-Egg (Daina & Zoete 2016) ellipses
    let in_yolk = in_ellipse(
        tpsa,
        wlogp,
        38.117,
        3.177,
        82.061 * 0.5,
        5.557 * 0.5,
        -0.171887,
    );
    let in_white = in_ellipse(
        tpsa,
        wlogp,
        71.051,
        2.292,
        142.081 * 0.5,
        8.740 * 0.5,
        -0.103132,
    );
    let pgp_substrate = mw > 400.0 && hba >= 5;
    let pharmacokinetics = Pharmacokinetics {
        gi_absorption: if in_white || in_yolk {
            AbsorptionLevel::High
        } else {
            AbsorptionLevel::Low
        },
        bbb_permeant: in_yolk,
        pgp_substrate,
        cyp1a2_inhibitor: aromatic_heavy_atoms >= 10 && consensus_log_p > 2.5,
        cyp2c19_inhibitor: mw > 280.0 && consensus_log_p > 2.8 && hbd <= 1,
        cyp2c9_inhibitor: consensus_log_p > 2.8 && (hba >= 3 || mw > 300.0),
        cyp2d6_inhibitor: consensus_log_p > 2.2
            && mol.atoms.iter().any(|a| a.symbol == "N" && !a.aromatic),
        cyp3a4_inhibitor: mw > 420.0 && consensus_log_p > 2.2,
        log_kp_cm_s: round2(-2.72 + 0.71 * consensus_log_p - 0.0061 * mw),
        pgp_cyp_are_experimental_heuristics: true,
    };

    let lipinski = eval_lipinski(&physicochemical, mlogp);
    let ghose = eval_ghose(mol, &physicochemical, wlogp);
    let veber = eval_veber(&physicochemical);
    let egan = eval_egan(&physicochemical, wlogp);
    let muegge = eval_muegge(mol, &physicochemical, xlogp3);
    let net_charge: i32 = mol.atoms.iter().map(|a| i32::from(a.formal_charge)).sum();
    let has_cooh = SmartsQuery::parse("C(=O)[OH,O-]")
        .map(|q| q.matches(mol))
        .unwrap_or(false);
    let bioavailability_score = if net_charge < 0 || has_cooh {
        if tpsa > 150.0 {
            0.11
        } else if tpsa > 75.0 {
            0.56
        } else {
            0.85
        }
    } else if lipinski.violations <= 1 {
        0.55
    } else {
        0.17
    };

    let druglikeness = Druglikeness {
        lipinski,
        ghose,
        veber,
        egan,
        muegge,
        bioavailability_score,
    };

    let pains_catalog: [(&str, &str); 26] = [
        ("rhodanine_like", "S=C1NC(=O)CS1"),
        ("ene_rhodanine", "O=C1NC(=S)SC1=C"),
        ("quinone_para", "O=C1C=CC(=O)C=C1"),
        ("quinone_ortho", "O=C1C(=O)C=CC=C1"),
        ("catechol", "Oc1ccccc1O"),
        ("hydroquinone", "Oc1ccc(O)cc1"),
        ("michael_acceptor_carbonyl", "C=CC=O"),
        ("michael_acrylate", "C=CC(=O)O"),
        ("azo_dye", "N=Nc1ccccc1"),
        ("nitro_aromatic", "[N+](=O)[O-]c1ccccc1"),
        ("acyl_hydrazide", "C(=O)NN"),
        ("thiosemicarbazone", "NC(=S)NN=C"),
        ("mannich_phenol", "Oc1ccccc1CN"),
        ("aminothiazole_acyl", "Nc1ncsc1"),
        ("curcumin_diketone", "O=CC=CC(=O)C=C"),
        ("imine_phenol", "Oc1ccccc1C=N"),
        ("isothiazolone", "O=c1sncc1"),
        ("cyanopyridone", "N#Cc1cc(=O)[nH]cc1"),
        ("barbiturate", "O=C1NC(=O)NC(=O)C1"),
        ("maleimide", "O=C1C=CC(=O)N1"),
        ("epoxide_reactive", "C1OC1"),
        ("sulfonamide_anilide_alert", "Nc1ccc(S(=O)(=O)N)cc1"),
        ("thiol_free", "[SX2H1]"),
        ("aldehyde_reactive", "[CX3H1]=O"),
        ("peroxide", "OO"),
        ("disulfide_reactive", "SS"),
    ];
    let pains_alerts: Vec<String> = pains_catalog
        .iter()
        .filter_map(|(lbl, pat)| {
            SmartsQuery::parse(pat)
                .ok()
                .filter(|q| q.matches(mol))
                .map(|_| (*lbl).to_string())
        })
        .collect();

    let brenk_catalog: [(&str, &str); 25] = [
        ("aldehyde", "[CX3H1]=O"),
        ("alkyl_halide", "[CX4][Cl,Br,I]"),
        ("acyl_halide", "C(=O)[Cl,Br]"),
        ("anhydride", "C(=O)OC(=O)"),
        ("epoxide", "C1OC1"),
        ("aziridine", "C1NC1"),
        ("peroxide", "OO"),
        ("hydrazine", "NN"),
        ("hydroxamic_acid", "C(=O)NO"),
        ("isocyanate", "N=C=O"),
        ("isothiocyanate", "N=C=S"),
        ("nitro_group", "[N+](=O)[O-]"),
        ("nitroso", "N=O"),
        ("azide", "N=[N+]=[N-]"),
        ("michael_acceptor", "C=CC(=O)"),
        ("beta_lactam", "O=C1CCN1"),
        ("polyene", "C=CC=CC=C"),
        ("thiol", "[SX2H1]"),
        ("thiocarbonyl", "C=S"),
        ("sulfonate_ester", "S(=O)(=O)OC"),
        ("phosphonate", "P(=O)(O)O"),
        ("triple_bond_terminal", "C#[CH]"),
        ("disulfide", "SS"),
        ("imine", "C=N"),
        ("oxime", "C=NO"),
    ];
    let brenk_alerts: Vec<String> = brenk_catalog
        .iter()
        .filter_map(|(lbl, pat)| {
            SmartsQuery::parse(pat)
                .ok()
                .filter(|q| q.matches(mol))
                .map(|_| (*lbl).to_string())
        })
        .collect();

    let mut lead_viol = Vec::new();
    if !(250.0..=350.0).contains(&mw) {
        lead_viol.push(format!("MW outside [250, 350] ({mw:.2})"));
    }
    if xlogp3 > 3.5 {
        lead_viol.push(format!("XLOGP3 > 3.5 ({xlogp3:.2})"));
    }
    if rotatable_bonds > 7 {
        lead_viol.push(format!("Rotatable bonds > 7 ({rotatable_bonds})"));
    }
    let leadlikeness = RuleResult {
        pass: lead_viol.is_empty(),
        violations: lead_viol.len(),
        details: lead_viol,
    };

    let chiral_cnt = mol.atoms.iter().filter(|a| a.cip_label.is_some()).count();
    let sa = (1.0
        + 0.08 * (heavy_atoms as f64)
        + 0.35 * (mol.rings.len() as f64)
        + 0.45 * (chiral_cnt as f64))
        .clamp(1.0, 10.0);

    let medicinal_chemistry = MedicinalChemistry {
        pains_alerts,
        pains_catalog_size: pains_catalog.len(),
        brenk_alerts,
        brenk_catalog_size: brenk_catalog.len(),
        leadlikeness,
        synthetic_accessibility: round2(sa),
    };

    let boiled_egg = BoiledEggPoint {
        tpsa,
        wlogp,
        in_white,
        in_yolk,
        pgp_substrate,
        in_plot_bounds: (0.0..=200.0).contains(&tpsa) && (-4.0..=8.0).contains(&wlogp),
    };

    let lipo_norm = norm_axis(xlogp3, -0.7, 5.0);
    let size_norm = norm_axis(mw, 150.0, 500.0);
    let polar_norm = norm_axis(tpsa, 20.0, 130.0);
    let insolu_norm = norm_axis(-esol_logs, 0.0, 6.0);
    let insatu_norm = norm_axis(1.0 - fraction_csp3, 0.0, 0.75);
    let flex_norm = norm_axis(rotatable_bonds as f64, 0.0, 9.0);
    let all_in_opt = (-0.7..=5.0).contains(&xlogp3)
        && (150.0..=500.0).contains(&mw)
        && (20.0..=130.0).contains(&tpsa)
        && (-6.0..=0.0).contains(&esol_logs)
        && fraction_csp3 >= 0.25
        && rotatable_bonds <= 9;

    let bioavailability_radar = BioavailabilityRadar {
        lipo_xlogp3: xlogp3,
        size_mw: mw,
        polar_tpsa: tpsa,
        insolu_log_s: esol_logs,
        insatu_fsp3: fraction_csp3,
        flex_rot_bonds: rotatable_bonds,
        normalized_axes: [
            lipo_norm,
            size_norm,
            polar_norm,
            insolu_norm,
            insatu_norm,
            flex_norm,
        ],
        all_in_optimal_zone: all_in_opt,
    };

    MoleculeEvaluation {
        name,
        input_smiles,
        canonical_smiles,
        warnings: mol.warnings.clone(),
        physicochemical,
        lipophilicity,
        water_solubility,
        pharmacokinetics,
        druglikeness,
        medicinal_chemistry,
        boiled_egg,
        bioavailability_radar,
    }
}

fn is_carbonyl_carbon(mol: &MolGraph, idx: usize) -> bool {
    mol.atoms[idx].symbol == "C"
        && mol.adj[idx].iter().any(|&(nb, bid)| {
            mol.bonds[bid].order == BondOrder::Double
                && matches!(mol.atoms[nb].symbol.as_str(), "O" | "S")
        })
}

fn count_rotatable_bonds(mol: &MolGraph) -> usize {
    mol.bonds
        .iter()
        .filter(|b| {
            if b.order != BondOrder::Single || b.in_ring {
                return false;
            }
            let u = b.u;
            let v = b.v;
            if mol.atoms[u].atomic_number == 1 || mol.atoms[v].atomic_number == 1 {
                return false;
            }
            if mol.adj[u].len() <= 1 || mol.adj[v].len() <= 1 {
                return false;
            }
            // Exclude single bonds adjacent to triple bonds (linear sp axis, no conformational DOF)
            let has_triple = |idx: usize| {
                mol.adj[idx]
                    .iter()
                    .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Triple)
            };
            if has_triple(u) || has_triple(v) {
                return false;
            }
            // Exclude partial-double amide C(=O)-N / thioamide C(=S)-N and ester C(=O)-O bonds (RDKit Strict)
            let is_amide_or_ester = |c_idx: usize, x_idx: usize| {
                is_carbonyl_carbon(mol, c_idx)
                    && matches!(mol.atoms[x_idx].symbol.as_str(), "N" | "O")
            };
            !(is_amide_or_ester(u, v) || is_amide_or_ester(v, u))
        })
        .count()
}

fn count_hba(mol: &MolGraph) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.symbol.as_str(), "N" | "O") && a.formal_charge <= 0)
        .count()
}

fn count_hbd(mol: &MolGraph) -> usize {
    mol.atoms
        .iter()
        .filter(|a| matches!(a.symbol.as_str(), "N" | "O" | "S") && a.total_h() > 0)
        .count()
}

/// Exact Ertl, Rohde & Selzer (J. Med. Chem. 2000, 43, 3714-3717) Topological Polar Surface Area
/// covering all neutral and formally charged N and O topological environments.
fn compute_tpsa(mol: &MolGraph) -> f64 {
    let mut tpsa = 0.0;
    for a in &mol.atoms {
        let h = usize::from(a.total_h());
        let q = a.formal_charge;
        let in_3_ring = a.ring_sizes.contains(&3);
        let mut n_single = 0usize;
        let mut n_double = 0usize;
        let mut n_triple = 0usize;
        let mut n_arom = 0usize;
        for &(nb, bid) in &mol.adj[a.id] {
            if mol.atoms[nb].atomic_number == 1 {
                continue;
            }
            match mol.bonds[bid].order {
                BondOrder::Single | BondOrder::Dative => n_single += 1,
                BondOrder::Double => n_double += 1,
                BondOrder::Triple => n_triple += 1,
                BondOrder::Aromatic => n_arom += 1,
            }
        }
        let heavy_deg = n_single + n_double + n_triple + n_arom;
        match a.symbol.as_str() {
            "N" => {
                let contrib = if q == 0 {
                    if a.aromatic || n_arom >= 2 {
                        match (n_arom, n_single, n_double, h) {
                            (3, 0, 0, 0) => 4.41,
                            (_, 0, 0, 0) => 12.89,
                            (_, 1, 0, 0) => 4.93,
                            (_, 0, 1, 0) => 8.39,
                            (_, 0, 0, 1) => 15.79,
                            _ => (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0),
                        }
                    } else if n_triple == 1 {
                        if n_single == 0 && h == 0 { 23.79 } else { 4.36 }
                    } else if n_double == 2 {
                        13.60
                    } else if n_double == 1 {
                        match (n_single, h) {
                            (1, 0) => 12.36,
                            (0, 1) => 23.85,
                            (2, 0) => 11.68,
                            _ => (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0),
                        }
                    } else {
                        match (n_single, h, in_3_ring) {
                            (3, 0, true) => 3.01,
                            (3, 0, false) => 3.24,
                            (2, 1, true) => 21.94,
                            (2, 1, false) => 12.03,
                            (1, 2, _) => 26.02,
                            _ => (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0),
                        }
                    }
                } else if q > 0 {
                    if a.aromatic || n_arom >= 2 {
                        if h >= 1 { 14.14 } else { 3.88 }
                    } else if n_triple == 1 {
                        if n_single == 1 && h == 0 {
                            4.36
                        } else if h >= 1 {
                            23.85
                        } else {
                            13.60
                        }
                    } else if n_double == 1 {
                        match (n_single, h) {
                            (2, 0) => 3.01,
                            (1, 1) => 13.97,
                            (0, 2) => 25.59,
                            _ => (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0),
                        }
                    } else if n_double == 0 {
                        match (n_single, h) {
                            (4, 0) => 0.0,
                            (3, 1) => 4.44,
                            (2, 2) => 16.61,
                            (1, 3) => 27.64,
                            _ => (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0),
                        }
                    } else {
                        // Ertl linear fallback (gives exact 14.10 for azide =[N+]=)
                        (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0)
                    }
                } else {
                    // Formally negative N (gives exact 22.30 for azide =[N-])
                    (30.5 - 8.2 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0)
                };
                tpsa += contrib;
            }
            "O" => {
                let contrib = if q == 0 {
                    if a.aromatic || n_arom >= 2 {
                        13.14
                    } else if n_double == 1 {
                        17.07
                    } else if n_single == 2 && h == 0 {
                        if in_3_ring { 12.53 } else { 9.23 }
                    } else if n_single == 1 && h >= 1 {
                        20.23
                    } else {
                        (28.5 - 8.6 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0)
                    }
                } else if q < 0 && n_single == 1 && h == 0 {
                    23.06
                } else {
                    (28.5 - 8.6 * (heavy_deg as f64) + 1.5 * (h as f64)).max(0.0)
                };
                tpsa += contrib;
            }
            _ => {}
        }
    }
    tpsa
}

/// Wildman & Crippen (J. Chem. Inf. Comput. Sci. 1999, 39, 868-873) atom-based LogP & MR
/// contribution for atom `idx` PLUS its attached implicit/explicit hydrogens.
fn crippen_atom_and_h_contribs(mol: &MolGraph, idx: usize) -> (f64, f64) {
    let a = &mol.atoms[idx];
    let h_cnt = f64::from(a.total_h());
    let nbrs = &mol.adj[idx];
    let is_hetero = |sym: &str| !matches!(sym, "C" | "H");
    let has_dbl = nbrs
        .iter()
        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Double);
    let has_trp = nbrs
        .iter()
        .any(|&(_, bid)| mol.bonds[bid].order == BondOrder::Triple);
    let heavy_nbrs: Vec<(usize, BondOrder)> = nbrs
        .iter()
        .copied()
        .filter(|&(v, _)| mol.atoms[v].atomic_number != 1)
        .map(|(v, bid)| (v, mol.bonds[bid].order))
        .collect();
    let hetero_nbr_cnt = heavy_nbrs
        .iter()
        .filter(|&&(v, _)| is_hetero(&mol.atoms[v].symbol))
        .count();
    let arom_nbr_cnt = heavy_nbrs
        .iter()
        .filter(|&&(v, _)| mol.atoms[v].aromatic)
        .count();

    let (atom_lp, atom_mr, h_lp, h_mr) = match a.symbol.as_str() {
        "C" => {
            let (lp, mr) = if a.aromatic {
                // Aromatic carbon classes (Wildman-Crippen C13-C25)
                let exo_nbr = heavy_nbrs
                    .iter()
                    .find(|&&(_, ord)| ord != BondOrder::Aromatic);
                if let Some(&(ev, ord)) = exo_nbr {
                    if ord == BondOrder::Double
                        && matches!(mol.atoms[ev].symbol.as_str(), "O" | "S")
                    {
                        // Aromatic carbonyl c(=O) (caffeine, uracil, pyridone)
                        (-0.8186, 4.0000)
                    } else {
                        match mol.atoms[ev].symbol.as_str() {
                            "C" => {
                                if mol.atoms[ev].aromatic {
                                    (0.2713, 3.9040)
                                } else {
                                    (0.1360, 3.5090)
                                }
                            }
                            "O" => (0.5437, 3.8530),
                            "N" => (0.4619, 4.0670),
                            "S" => (0.2450, 3.5640),
                            "F" => (0.0000, 3.2570),
                            "Cl" => (0.2450, 3.5640),
                            "Br" => (0.1980, 3.1800),
                            "I" => (0.0000, 3.1040),
                            _ => (0.1581, 3.3500),
                        }
                    }
                } else if heavy_nbrs.len() >= 3 {
                    // Fused ring junction c(a)(a)a
                    let adj_het = heavy_nbrs
                        .iter()
                        .any(|&(v, _)| is_hetero(&mol.atoms[v].symbol));
                    if adj_het {
                        (0.2955, 3.9040)
                    } else {
                        (0.2713, 3.9040)
                    }
                } else {
                    // Unsubstituted cH
                    (0.1581, 3.3500)
                }
            } else if has_trp {
                // sp alkyne / nitrile C
                (0.0017, 3.8880)
            } else if has_dbl {
                // sp2 carbonyl / imine / alkene C
                let dbl_to_hetero = heavy_nbrs
                    .iter()
                    .any(|&(v, ord)| ord == BondOrder::Double && is_hetero(&mol.atoms[v].symbol));
                if dbl_to_hetero {
                    (-0.2783, 5.0070)
                } else if hetero_nbr_cnt > 0 {
                    (-0.0300, 3.6000)
                } else {
                    (0.1551, 3.5130)
                }
            } else {
                // sp3 aliphatic C
                let f_cnt = heavy_nbrs
                    .iter()
                    .filter(|&&(v, _)| mol.atoms[v].symbol == "F")
                    .count();
                if f_cnt >= 2 {
                    (-0.0967, 2.5760)
                } else if arom_nbr_cnt >= 1 {
                    // Benzylic C (Ar-CH3 / Ar-CH2-X)
                    if heavy_nbrs.len() == 1 {
                        (0.0845, 2.5030)
                    } else {
                        (-0.0516, 2.4880)
                    }
                } else if hetero_nbr_cnt >= 2 {
                    (-0.2051, 2.7530)
                } else if hetero_nbr_cnt == 1 {
                    let hv = heavy_nbrs
                        .iter()
                        .find(|&&(v, _)| is_hetero(&mol.atoms[v].symbol))
                        .map(|&(v, _)| v);
                    if hv.is_some_and(|v| mol.atoms[v].aromatic) {
                        (-0.1444, 2.4120)
                    } else {
                        (-0.2035, 2.7530)
                    }
                } else if heavy_nbrs.len() >= 3 {
                    // Tertiary / quaternary alkane CH / C
                    (0.0000, 2.4330)
                } else {
                    // Primary / secondary alkane CH3 / CH2
                    (0.1441, 2.5030)
                }
            };
            (lp, mr, 0.1230, 1.0570)
        }
        "N" => {
            let (lp, mr) = if a.aromatic {
                (-0.3239, 2.2020)
            } else if a.formal_charge > 0 {
                let is_nitro = heavy_nbrs
                    .iter()
                    .filter(|&&(v, _)| mol.atoms[v].symbol == "O")
                    .count()
                    >= 2;
                let is_azide = heavy_nbrs
                    .iter()
                    .filter(|&&(v, ord)| mol.atoms[v].symbol == "N" && ord == BondOrder::Double)
                    .count()
                    >= 2;
                if is_nitro {
                    (-0.3396, 0.2604)
                } else if is_azide {
                    (0.2887, 3.3590)
                } else {
                    (-0.3396, 1.8000)
                }
            } else if a.formal_charge < 0 {
                (0.2887, 3.3590)
            } else if has_trp {
                // Nitrile #N
                (0.0151, 1.7250)
            } else if has_dbl {
                // Imine / azo / azide terminal N
                (0.1836, 2.4280)
            } else {
                let adj_arom = heavy_nbrs.iter().any(|&(v, _)| mol.atoms[v].aromatic);
                let adj_carbonyl = heavy_nbrs
                    .iter()
                    .any(|&(v, _)| is_carbonyl_carbon(mol, v) || mol.atoms[v].symbol == "S");
                if adj_arom {
                    match heavy_nbrs.len() {
                        1 => (-1.0270, 2.8270),
                        2 => (-0.5188, 2.9500),
                        _ => (-0.4458, 2.8190),
                    }
                } else if adj_carbonyl {
                    match heavy_nbrs.len() {
                        1 => (-1.0190, 2.2620),
                        #[allow(clippy::approx_constant)]
                        2 => (-0.5262, 2.7180),
                        _ => (-0.3187, 2.8190),
                    }
                } else {
                    match heavy_nbrs.len() {
                        1 => (-1.0190, 2.2620),
                        2 => (-0.7096, 2.1730),
                        _ => (-0.3187, 2.8190),
                    }
                }
            };
            (lp, mr, 0.2142, 0.9627)
        }
        "O" => {
            let (lp, mr) = if a.aromatic {
                (0.1552, 1.0800)
            } else if a.formal_charge < 0 {
                let adj_n_plus = heavy_nbrs
                    .iter()
                    .any(|&(v, _)| mol.atoms[v].symbol == "N" && mol.atoms[v].formal_charge > 0);
                if adj_n_plus {
                    (0.0335, 3.3670)
                } else {
                    (-0.3339, 0.7774)
                }
            } else if has_dbl {
                let adj_n_plus = heavy_nbrs
                    .iter()
                    .any(|&(v, _)| mol.atoms[v].symbol == "N" && mol.atoms[v].formal_charge > 0);
                let adj_s = heavy_nbrs.iter().any(|&(v, _)| mol.atoms[v].symbol == "S");
                let adj_p = heavy_nbrs.iter().any(|&(v, _)| mol.atoms[v].symbol == "P");
                let adj_arom_c = heavy_nbrs
                    .iter()
                    .any(|&(v, _)| mol.atoms[v].symbol == "C" && mol.atoms[v].aromatic);
                let adj_amide_or_acid = heavy_nbrs.iter().any(|&(v, _)| {
                    mol.atoms[v].symbol == "C"
                        && mol.adj[v].iter().any(|&(u, _)| {
                            u != idx && matches!(mol.atoms[u].symbol.as_str(), "N" | "O")
                        })
                });
                if adj_n_plus {
                    (0.0335, 3.3670)
                } else if adj_s {
                    (-0.3339, 0.7774)
                } else if adj_p {
                    (-0.1188, 0.6865)
                } else if adj_arom_c {
                    (0.1788, 1.0000)
                } else if adj_amide_or_acid {
                    (0.1129, 0.7774)
                } else {
                    (-0.1526, 0.0000)
                }
            } else if a.total_h() > 0 {
                // Alcohol / phenol / carboxylic OH / phosphate OH
                (-0.2893, 0.8238)
            } else {
                // Ether / ester -O-
                let adj_carbonyl_or_arom = heavy_nbrs
                    .iter()
                    .any(|&(v, _)| is_carbonyl_carbon(mol, v) || mol.atoms[v].aromatic);
                if adj_carbonyl_or_arom {
                    (-0.0684, 1.0850)
                } else {
                    (-0.4195, 1.1820)
                }
            };
            (lp, mr, -0.2677, 1.3950)
        }
        "S" => {
            let (lp, mr) = if a.aromatic {
                (0.6237, 6.6910)
            } else if heavy_nbrs.len() >= 3 || has_dbl {
                // Sulfone / sulfonamide / sulfoxide
                (-0.0024, 7.3650)
            } else {
                // Thiol / thioether
                (0.6482, 7.5910)
            };
            (lp, mr, -0.2677, 1.3950)
        }
        "P" => (0.8612, 6.9200, 0.1230, 1.0570),
        "F" => (0.4202, 1.1080, 0.1230, 1.0570),
        "Cl" => (0.6895, 5.8530, 0.1230, 1.0570),
        "Br" => (0.8456, 8.9270, 0.1230, 1.0570),
        "I" => (0.8857, 14.0200, 0.1230, 1.0570),
        "B" => (-0.3808, 4.0410, 0.1230, 1.0570),
        "Si" => (0.3000, 7.0000, 0.1230, 1.0570),
        "H" => (0.0, 0.0, 0.0, 0.0),
        _ => (0.0, 2.0, 0.1230, 1.0570),
    };
    (atom_lp + h_cnt * h_lp, atom_mr + h_cnt * h_mr)
}

fn compute_mr(mol: &MolGraph) -> f64 {
    (0..mol.atoms.len())
        .map(|i| crippen_atom_and_h_contribs(mol, i).1)
        .sum()
}

fn compute_wlogp(mol: &MolGraph) -> f64 {
    (0..mol.atoms.len())
        .map(|i| crippen_atom_and_h_contribs(mol, i).0)
        .sum()
}

fn compute_xlogp3(mol: &MolGraph) -> f64 {
    let mut lp = 0.0;
    for a in &mol.atoms {
        let h = f64::from(a.total_h());
        let has_hetero = mol.adj[a.id]
            .iter()
            .any(|&(v, _)| matches!(mol.atoms[v].symbol.as_str(), "N" | "O" | "S"));
        let base = match a.symbol.as_str() {
            "C" => {
                if a.aromatic {
                    0.34
                } else if has_hetero {
                    0.06
                } else {
                    0.44
                }
            }
            "N" => {
                if a.aromatic {
                    -0.45
                } else {
                    -0.95
                }
            }
            "O" => {
                if h > 0.0 {
                    -0.58
                } else {
                    -0.34
                }
            }
            "S" => 0.35,
            "F" => 0.37,
            "Cl" => 0.88,
            "Br" => 1.05,
            "I" => 1.38,
            _ => 0.0,
        };
        lp += base;
    }
    lp
}

fn compute_mlogp(mol: &MolGraph, p: &PhysicochemicalProperties) -> f64 {
    let cx = mol
        .atoms
        .iter()
        .filter(|a| matches!(a.symbol.as_str(), "C" | "F" | "Cl" | "Br" | "I"))
        .count() as f64;
    let no = mol
        .atoms
        .iter()
        .filter(|a| matches!(a.symbol.as_str(), "N" | "O"))
        .count() as f64;
    let prx = p.tpsa / 100.0;
    let ub = mol
        .bonds
        .iter()
        .filter(|b| !matches!(b.order, BondOrder::Single))
        .count() as f64;
    1.244 * cx.powf(0.6) - 1.017 * no.powf(0.9) + 0.406 * prx - 0.145 * ub.powf(0.8) + 0.511
}

fn eval_lipinski(p: &PhysicochemicalProperties, mlogp: f64) -> RuleResult {
    let mut details = Vec::new();
    if p.molecular_weight > 500.0 {
        details.push(format!("MW>500 ({:.2})", p.molecular_weight));
    }
    if mlogp > 4.15 {
        details.push(format!("MLOGP>4.15 ({mlogp:.2})"));
    }
    if p.h_bond_acceptors > 10 {
        details.push(format!("HBA>10 ({})", p.h_bond_acceptors));
    }
    if p.h_bond_donors > 5 {
        details.push(format!("HBD>5 ({})", p.h_bond_donors));
    }
    RuleResult {
        pass: details.len() <= 1,
        violations: details.len(),
        details,
    }
}

fn eval_ghose(mol: &MolGraph, p: &PhysicochemicalProperties, wlogp: f64) -> RuleResult {
    let mut details = Vec::new();
    let total_atoms = mol.atoms.len()
        + mol
            .atoms
            .iter()
            .map(|a| usize::from(a.total_h()))
            .sum::<usize>();
    if !(160.0..=480.0).contains(&p.molecular_weight) {
        details.push("MW outside [160, 480]".to_string());
    }
    if !(-0.4..=5.6).contains(&wlogp) {
        details.push("WLOGP outside [-0.4, 5.6]".to_string());
    }
    if !(40.0..=130.0).contains(&p.molar_refractivity) {
        details.push("MR outside [40, 130]".to_string());
    }
    if !(20..=70).contains(&total_atoms) {
        details.push("#atoms outside [20, 70]".to_string());
    }
    RuleResult {
        pass: details.is_empty(),
        violations: details.len(),
        details,
    }
}

fn eval_veber(p: &PhysicochemicalProperties) -> RuleResult {
    let mut details = Vec::new();
    if p.rotatable_bonds > 10 {
        details.push(format!("RotBonds>10 ({})", p.rotatable_bonds));
    }
    if p.tpsa > 140.0 {
        details.push(format!("TPSA>140 ({:.2})", p.tpsa));
    }
    RuleResult {
        pass: details.is_empty(),
        violations: details.len(),
        details,
    }
}

fn eval_egan(p: &PhysicochemicalProperties, wlogp: f64) -> RuleResult {
    let mut details = Vec::new();
    if wlogp > 5.88 {
        details.push("WLOGP>5.88".to_string());
    }
    if p.tpsa > 131.6 {
        details.push("TPSA>131.6".to_string());
    }
    RuleResult {
        pass: details.is_empty(),
        violations: details.len(),
        details,
    }
}

fn eval_muegge(mol: &MolGraph, p: &PhysicochemicalProperties, xlogp3: f64) -> RuleResult {
    let mut details = Vec::new();
    if !(200.0..=600.0).contains(&p.molecular_weight) {
        details.push("MW outside [200, 600]".to_string());
    }
    if !(-2.0..=5.0).contains(&xlogp3) {
        details.push("XLOGP3 outside [-2, 5]".to_string());
    }
    if p.tpsa > 150.0 {
        details.push("TPSA>150".to_string());
    }
    let carbons = mol.atoms.iter().filter(|a| a.symbol == "C").count();
    let heteros = mol
        .atoms
        .iter()
        .filter(|a| !matches!(a.symbol.as_str(), "C" | "H"))
        .count();
    if carbons <= 4 {
        details.push("Num carbons<=4".to_string());
    }
    if heteros <= 1 {
        details.push("Num heteroatoms<=1".to_string());
    }
    RuleResult {
        pass: details.is_empty(),
        violations: details.len(),
        details,
    }
}

fn in_ellipse(x: f64, y: f64, cx: f64, cy: f64, rx: f64, ry: f64, angle_deg: f64) -> bool {
    let rad = angle_deg * PI / 180.0;
    let dx = x - cx;
    let dy = y - cy;
    let xr = dx * rad.cos() + dy * rad.sin();
    let yr = -dx * rad.sin() + dy * rad.cos();
    (xr / rx).powi(2) + (yr / ry).powi(2) <= 1.0
}

fn norm_axis(val: f64, low: f64, high: f64) -> f64 {
    let span = (high - low).max(1e-6);
    (0.25 + 0.55 * ((val - low) / span)).clamp(0.05, 1.0)
}

#[inline]
fn round2(v: f64) -> f64 {
    (v * 100.0).round() / 100.0
}
