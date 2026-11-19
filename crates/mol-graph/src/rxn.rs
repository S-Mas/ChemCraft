//! Clean-room MDL `.rxn` V2000 chemical reaction reader/writer, reaction SMILES parser,
//! and ChemDraw-equivalent Reaction Stoichiometry / Atom Economy / Yield Calculator.

use std::fmt::Write;

use serde::{Deserialize, Serialize};

use crate::error::MolError;
use crate::mol::MolGraph;
use crate::sdf::{parse_mol_block, write_mol_v2000_block};
use crate::smiles::parse_smiles;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReactionRole {
    Reactant,
    Agent,
    Product,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionScheme {
    pub name: String,
    pub reactants: Vec<MolGraph>,
    pub agents: Vec<MolGraph>,
    pub products: Vec<MolGraph>,
    pub conditions: String,
}

impl ReactionScheme {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            reactants: Vec::new(),
            agents: Vec::new(),
            products: Vec::new(),
            conditions: String::new(),
        }
    }

    /// Format as standard Daylight Reaction SMILES (`R1.R2>A1>P1.P2`).
    pub fn reaction_smiles(&self) -> String {
        let r = self
            .reactants
            .iter()
            .map(MolGraph::canonical_smiles)
            .collect::<Vec<_>>()
            .join(".");
        let a = self
            .agents
            .iter()
            .map(MolGraph::canonical_smiles)
            .collect::<Vec<_>>()
            .join(".");
        let p = self
            .products
            .iter()
            .map(MolGraph::canonical_smiles)
            .collect::<Vec<_>>()
            .join(".");
        format!("{r}>{a}>{p}")
    }
}

/// Parse Daylight Reaction SMILES (`reactants>agents>products`, e.g. `CC(=O)O.CCO>>CC(=O)OCC.O`).
pub fn parse_reaction_smiles(input: &str) -> Result<ReactionScheme, MolError> {
    let trimmed = input.trim();
    let parts: Vec<&str> = trimmed.split('>').collect();
    if parts.len() != 3 {
        return Err(MolError::InvalidFormat(
            "Reaction SMILES must have exactly two '>' delimiters (reactants>agents>products)"
                .to_string(),
        ));
    }
    let parse_side = |s: &str| -> Result<Vec<MolGraph>, MolError> {
        let mut out = Vec::new();
        for frag in s.split('.').map(str::trim).filter(|f| !f.is_empty()) {
            out.push(parse_smiles(frag)?);
        }
        Ok(out)
    };
    let reactants = parse_side(parts[0])?;
    let agents = parse_side(parts[1])?;
    let products = parse_side(parts[2])?;
    if reactants.is_empty() && products.is_empty() {
        return Err(MolError::EmptyInput);
    }
    Ok(ReactionScheme {
        name: "Reaction".to_string(),
        reactants,
        agents,
        products,
        conditions: String::new(),
    })
}

/// Serialize `rxn` to standard MDL `$RXN` V2000 format (readable by ChemDraw, RDKit, and Indigo).
pub fn write_rxn_v2000(rxn: &ReactionScheme) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "$RXN");
    let _ = writeln!(out, "{}", rxn.name);
    let _ = writeln!(out, "  ChemCraft 1119260900");
    let _ = writeln!(out, "{}", rxn.conditions);
    if rxn.agents.is_empty() {
        let _ = writeln!(out, "{:>3}{:>3}", rxn.reactants.len(), rxn.products.len());
    } else {
        let _ = writeln!(
            out,
            "{:>3}{:>3}{:>3}",
            rxn.reactants.len(),
            rxn.products.len(),
            rxn.agents.len()
        );
    }
    for m in &rxn.reactants {
        let _ = writeln!(out, "$MOL");
        out.push_str(&write_mol_v2000_block(m));
    }
    for m in &rxn.products {
        let _ = writeln!(out, "$MOL");
        out.push_str(&write_mol_v2000_block(m));
    }
    for m in &rxn.agents {
        let _ = writeln!(out, "$MOL");
        out.push_str(&write_mol_v2000_block(m));
    }
    out
}

/// Parse an MDL `$RXN` V2000 reaction file into a `ReactionScheme`.
pub fn parse_rxn_v2000(text: &str) -> Result<ReactionScheme, MolError> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 5 || !lines[0].trim_start().starts_with("$RXN") {
        return Err(MolError::InvalidFormat(
            "Input is not an MDL $RXN file".to_string(),
        ));
    }
    let name = lines
        .get(1)
        .copied()
        .unwrap_or("Reaction")
        .trim()
        .to_string();
    let conditions = lines.get(3).copied().unwrap_or("").trim().to_string();
    let counts_line = lines.get(4).copied().unwrap_or("");
    let n_reactants = counts_line
        .get(0..3)
        .unwrap_or("")
        .trim()
        .parse::<usize>()
        .or_else(|_| {
            counts_line
                .split_whitespace()
                .next()
                .unwrap_or("0")
                .parse::<usize>()
        })
        .map_err(|_| MolError::InvalidFormat("Invalid $RXN reactant count".to_string()))?;
    let n_products = counts_line
        .get(3..6)
        .unwrap_or("")
        .trim()
        .parse::<usize>()
        .or_else(|_| {
            counts_line
                .split_whitespace()
                .nth(1)
                .unwrap_or("0")
                .parse::<usize>()
        })
        .map_err(|_| MolError::InvalidFormat("Invalid $RXN product count".to_string()))?;
    let n_agents = counts_line
        .get(6..9)
        .unwrap_or("")
        .trim()
        .parse::<usize>()
        .unwrap_or(0);

    // Split remaining text by `$MOL` delimiters
    let rest = lines[5..].join("\n");
    let mut blocks: Vec<&str> = Vec::new();
    for chunk in rest.split("$MOL") {
        let c = chunk.trim_start_matches(['\r', '\n']);
        if !c.trim().is_empty() {
            blocks.push(c);
        }
    }
    if blocks.len() < n_reactants + n_products {
        return Err(MolError::InvalidFormat(format!(
            "Expected {} $MOL blocks ({n_reactants} reactants + {n_products} products), found {}",
            n_reactants + n_products,
            blocks.len()
        )));
    }

    let mut reactants = Vec::with_capacity(n_reactants);
    for blk in blocks.iter().take(n_reactants) {
        reactants.push(parse_mol_block(blk)?);
    }
    let mut products = Vec::with_capacity(n_products);
    for blk in blocks.iter().skip(n_reactants).take(n_products) {
        products.push(parse_mol_block(blk)?);
    }
    let mut agents = Vec::with_capacity(n_agents);
    for blk in blocks.iter().skip(n_reactants + n_products).take(n_agents) {
        agents.push(parse_mol_block(blk)?);
    }

    Ok(ReactionScheme {
        name: if name.is_empty() {
            "Reaction".to_string()
        } else {
            name
        },
        reactants,
        agents,
        products,
        conditions,
    })
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StoichRow {
    pub role: ReactionRole,
    pub name: String,
    pub formula: String,
    pub mw: f64,
    pub equivalents: f64,
    pub mmol: f64,
    pub mass_mg: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReactionStoichiometry {
    pub rows: Vec<StoichRow>,
    pub atom_economy_pct: f64,
    pub theoretical_yield_mg: f64,
    pub isolated_yield_pct: Option<f64>,
}

/// Compute ChemDraw-style reaction stoichiometry table, atom economy, theoretical yield, and % yield.
pub fn compute_stoichiometry(
    rxn: &ReactionScheme,
    limiting_mmol: f64,
    actual_yield_mg: Option<f64>,
) -> ReactionStoichiometry {
    let mmol = limiting_mmol.max(0.0);
    let mut rows = Vec::new();
    let mut total_reactant_mw = 0.0_f64;

    for (i, m) in rxn.reactants.iter().enumerate() {
        let mw = m.molecular_weight();
        total_reactant_mw += mw;
        rows.push(StoichRow {
            role: ReactionRole::Reactant,
            name: m
                .name
                .clone()
                .unwrap_or_else(|| format!("Reactant {}", i + 1)),
            formula: m.formula(),
            mw: (mw * 100.0).round() / 100.0,
            equivalents: 1.0,
            mmol: (mmol * 1000.0).round() / 1000.0,
            mass_mg: (mmol * mw * 100.0).round() / 100.0,
        });
    }
    for (i, m) in rxn.agents.iter().enumerate() {
        let mw = m.molecular_weight();
        rows.push(StoichRow {
            role: ReactionRole::Agent,
            name: m.name.clone().unwrap_or_else(|| format!("Agent {}", i + 1)),
            formula: m.formula(),
            mw: (mw * 100.0).round() / 100.0,
            equivalents: 0.1,
            mmol: (mmol * 0.1 * 1000.0).round() / 1000.0,
            mass_mg: (mmol * 0.1 * mw * 100.0).round() / 100.0,
        });
    }
    let mut main_product_mw = 0.0_f64;
    for (i, m) in rxn.products.iter().enumerate() {
        let mw = m.molecular_weight();
        if i == 0 {
            main_product_mw = mw;
        }
        rows.push(StoichRow {
            role: ReactionRole::Product,
            name: m
                .name
                .clone()
                .unwrap_or_else(|| format!("Product {}", i + 1)),
            formula: m.formula(),
            mw: (mw * 100.0).round() / 100.0,
            equivalents: 1.0,
            mmol: (mmol * 1000.0).round() / 1000.0,
            mass_mg: (mmol * mw * 100.0).round() / 100.0,
        });
    }

    let atom_economy_pct = if total_reactant_mw > 1e-6 {
        ((100.0 * main_product_mw / total_reactant_mw) * 10.0).round() / 10.0
    } else {
        0.0
    };
    let theoretical_yield_mg = ((mmol * main_product_mw) * 100.0).round() / 100.0;
    let isolated_yield_pct = actual_yield_mg.and_then(|actual| {
        if theoretical_yield_mg > 1e-6 {
            Some(((100.0 * actual / theoretical_yield_mg) * 10.0).round() / 10.0)
        } else {
            None
        }
    });

    ReactionStoichiometry {
        rows,
        atom_economy_pct,
        theoretical_yield_mg,
        isolated_yield_pct,
    }
}
