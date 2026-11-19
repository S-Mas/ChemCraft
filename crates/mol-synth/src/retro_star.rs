#![forbid(unsafe_code)]

use std::collections::{HashMap, HashSet};
use std::time::Instant;

use mol_graph::{
    ChiralTag, MolGraph, ReactionScheme, compute_stoichiometry, parse_smiles, write_rxn_v2000,
};
use serde::{Deserialize, Serialize};

use crate::stock::{StockAnalog, StockCatalog, StockMatch, StockTier};
use crate::templates::{Disconnection, enumerate_disconnections};

/// Configuration for Retro* multi-step retrosynthesis planning.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RetroConfig {
    /// Maximum synthetic depth (number of steps in longest linear sequence, default `4`).
    pub max_steps: usize,
    /// Maximum number of top ranked routes to return (default `5`).
    pub max_routes: usize,
    /// Maximum reaction branches explored per intermediate node (default `8`).
    pub branch_limit: usize,
    /// Maximum total molecular node expansions to bound runtime (default `160`).
    pub max_expansions: usize,
    /// Target synthesis scale in milligrams for stoichiometry calculation (default `100.0` mg).
    pub target_scale_mg: f64,
}

impl Default for RetroConfig {
    fn default() -> Self {
        Self {
            max_steps: 4,
            max_routes: 5,
            branch_limit: 8,
            max_expansions: 160,
            target_scale_mg: 100.0,
        }
    }
}

/// One starting material on the Buy-Before-Make shopping list.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StartingMaterial {
    pub canonical_smiles: String,
    pub name: String,
    pub formula: String,
    pub mw: f64,
    pub tier: StockTier,
    pub catalog_id: Option<String>,
    pub supplier: Option<String>,
    pub est_usd_per_g: Option<f64>,
    pub required_mg: f64,
}

/// One forward reaction step inside a planned synthetic route.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SynthesisStep {
    pub step_number: usize,
    pub template_id: String,
    pub reaction_name: String,
    pub reaction_class: String,
    pub reactant_smiles: Vec<String>,
    pub product_smiles: String,
    pub reaction_smiles: String,
    pub reagents: String,
    pub solvent: String,
    pub temperature: String,
    pub est_yield_pct: f64,
    pub atom_economy_pct: f64,
    pub product_scale_mg: f64,
    pub reactant_masses_mg: Vec<(String, f64, f64)>, // (smiles, mmol, mass_mg)
    pub workup_notes: String,
    pub warnings: Vec<String>,
}

impl SynthesisStep {
    /// Builds a `mol_graph::ReactionScheme` for this step and exports MDL `.rxn` V2000 text.
    #[must_use]
    pub fn to_rxn_v2000(&self) -> Option<String> {
        let mut rxn = ReactionScheme::new(&self.reaction_name);
        rxn.conditions = format!(
            "{} | {} | {}",
            self.reagents, self.solvent, self.temperature
        );
        for rsmi in &self.reactant_smiles {
            let m = parse_smiles(rsmi).ok()?;
            rxn.reactants.push(m);
        }
        rxn.products.push(parse_smiles(&self.product_smiles).ok()?);
        Some(write_rxn_v2000(&rxn))
    }
}

/// Complete ranked multi-step synthetic route from commercial/simple building blocks to target.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SynthesisRoute {
    pub rank: usize,
    pub total_steps: usize,
    pub longest_linear_seq: usize,
    pub fully_resolved_in_stock: bool,
    pub stock_coverage_pct: f64,
    pub overall_est_yield_pct: f64,
    pub mean_atom_economy_pct: f64,
    pub composite_score: f64,
    pub warning_count: usize,
    pub steps: Vec<SynthesisStep>,
    pub starting_materials: Vec<StartingMaterial>,
}

/// Full output of the `mol-synth` Buy-Before-Make + Retro* retrosynthesis planner.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SynthesisPlan {
    pub target_name: String,
    pub target_smiles: String,
    pub target_formula: String,
    pub target_mw: f64,
    pub target_scale_mg: f64,
    pub target_stock_match: StockMatch,
    pub nearest_commercial_analogs: Vec<StockAnalog>,
    pub catalog_size: usize,
    pub nodes_expanded: usize,
    pub elapsed_ms: u128,
    pub routes: Vec<SynthesisRoute>,
}

impl SynthesisPlan {
    /// Renders a clean, ELN-ready Markdown retrosynthesis & lab protocol report.
    #[must_use]
    pub fn render_markdown(&self) -> String {
        let mut md = String::new();
        md.push_str(&format!("# Retrosynthesis Plan: {}\n\n", self.target_name));
        md.push_str(&format!(
            "- **Target SMILES:** `{}`\n- **Formula / MW:** `{}` ({:.2} Da) | **Target Scale:** `{:.1} mg`\n- **Buy-Before-Make Check:** **{}**",
            self.target_smiles,
            self.target_formula,
            self.target_mw,
            self.target_scale_mg,
            self.target_stock_match.tier.as_str()
        ));
        if let Some(ref e) = self.target_stock_match.entry {
            md.push_str(&format!(
                " — `{}` *{}* (~${:.2}/g, {})\n",
                e.catalog_id, e.name, e.est_usd_per_g, e.supplier
            ));
        } else {
            md.push('\n');
        }
        md.push_str(&format!(
            "- **Catalog Indexed:** `{}` commercial building blocks | **Search:** `{} nodes` in `{} ms`\n\n",
            self.catalog_size, self.nodes_expanded, self.elapsed_ms
        ));

        if !self.nearest_commercial_analogs.is_empty() {
            md.push_str("## Nearest Commercial Building Blocks (Buy-Before-Make)\n\n");
            md.push_str(
                "| Catalog ID | Name | Tanimoto | Est. $/g | SMILES |\n|---|---|---:|---:|---|\n",
            );
            for a in &self.nearest_commercial_analogs {
                md.push_str(&format!(
                    "| `{}` | {} | {:.3} | ${:.2} | `{}` |\n",
                    a.catalog_id, a.name, a.tanimoto, a.est_usd_per_g, a.canonical_smiles
                ));
            }
            md.push('\n');
        }

        if self.routes.is_empty() {
            md.push_str(
                "_No retrosynthetic disconnections matched within the configured step budget._\n",
            );
            return md;
        }

        for r in &self.routes {
            md.push_str(&format!(
                "## Route #{} — {} step(s) (LLS {}), {:.0}% stock coverage, est. {:.1}% overall yield, {:.1}% atom economy\n\n",
                r.rank,
                r.total_steps,
                r.longest_linear_seq,
                r.stock_coverage_pct,
                r.overall_est_yield_pct,
                r.mean_atom_economy_pct
            ));

            md.push_str("### Forward Synthesis Protocol\n\n");
            for st in &r.steps {
                md.push_str(&format!(
                    "#### Step {}: {} (`{}`)\n",
                    st.step_number, st.reaction_name, st.template_id
                ));
                md.push_str(&format!(
                    "- **Reaction SMILES:** `{}`\n",
                    st.reaction_smiles
                ));
                md.push_str(&format!(
                    "- **Conditions:** {} | **Solvent:** {} | **Temp/Time:** {}\n",
                    st.reagents, st.solvent, st.temperature
                ));
                md.push_str(&format!(
                    "- **Est. Yield:** `{:.0}%` | **Atom Economy:** `{:.1}%` | **Product Scale:** `{:.1} mg`\n",
                    st.est_yield_pct, st.atom_economy_pct, st.product_scale_mg
                ));
                for (i, (rsmi, mmol, mg)) in st.reactant_masses_mg.iter().enumerate() {
                    md.push_str(&format!(
                        "  - Reactant {}: `{}` — **{:.1} mg** ({:.3} mmol)\n",
                        i + 1,
                        rsmi,
                        mg,
                        mmol
                    ));
                }
                md.push_str(&format!("- **Workup:** {}\n", st.workup_notes));
                for w in &st.warnings {
                    md.push_str(&format!("- **Chemoselectivity / PG Alert:** {}\n", w));
                }
                md.push('\n');
            }

            md.push_str("### Starting Materials Shopping List\n\n");
            md.push_str("| Status | Name / ID | Formula | MW | Req. Mass | Est. $/g | SMILES |\n|---|---|---|---:|---:|---:|---|\n");
            for sm in &r.starting_materials {
                let id_str = sm
                    .catalog_id
                    .as_deref()
                    .map_or_else(|| sm.name.clone(), |cid| format!("{} (`{cid}`)", sm.name));
                let price_str = sm
                    .est_usd_per_g
                    .map_or_else(|| "—".to_string(), |p| format!("${p:.2}"));
                md.push_str(&format!(
                    "| {} | {} | `{}` | {:.1} | {:.1} mg | {} | `{}` |\n",
                    sm.tier.as_str(),
                    id_str,
                    sm.formula,
                    sm.mw,
                    sm.required_mg,
                    price_str,
                    sm.canonical_smiles
                ));
            }
            md.push('\n');
        }
        md
    }
}

/// Internal AND-OR solved subtree used during Retro* search.
#[derive(Debug, Clone)]
enum AndOrTree {
    Leaf {
        smiles: String,
        mol: MolGraph,
        stock: StockMatch,
        cost: f64,
    },
    Step {
        product_smiles: String,
        product_mol: MolGraph,
        disc: Disconnection,
        children: Vec<AndOrTree>,
        total_cost: f64,
        total_steps: usize,
        lls: usize,
        unresolved_leaves: usize,
        total_leaves: usize,
        overall_yield: f64,
        warning_count: usize,
    },
}

impl AndOrTree {
    fn total_cost(&self) -> f64 {
        match self {
            Self::Leaf { cost, .. } => *cost,
            Self::Step { total_cost, .. } => *total_cost,
        }
    }

    fn total_steps(&self) -> usize {
        match self {
            Self::Leaf { .. } => 0,
            Self::Step { total_steps, .. } => *total_steps,
        }
    }

    fn lls(&self) -> usize {
        match self {
            Self::Leaf { .. } => 0,
            Self::Step { lls, .. } => *lls,
        }
    }

    fn unresolved_leaves(&self) -> usize {
        match self {
            Self::Leaf { stock, .. } => usize::from(!stock.tier.is_buyable()),
            Self::Step {
                unresolved_leaves, ..
            } => *unresolved_leaves,
        }
    }

    fn total_leaves(&self) -> usize {
        match self {
            Self::Leaf { .. } => 1,
            Self::Step { total_leaves, .. } => *total_leaves,
        }
    }

    fn overall_yield(&self) -> f64 {
        match self {
            Self::Leaf { .. } => 1.0,
            Self::Step { overall_yield, .. } => *overall_yield,
        }
    }

    fn warning_count(&self) -> usize {
        match self {
            Self::Leaf { .. } => 0,
            Self::Step { warning_count, .. } => *warning_count,
        }
    }

    fn signature(&self) -> String {
        match self {
            Self::Leaf { smiles, .. } => format!("L:{smiles}"),
            Self::Step { disc, children, .. } => {
                let mut cs: Vec<String> = children.iter().map(Self::signature).collect();
                cs.sort();
                format!("S:{}({})", disc.template_id, cs.join("+"))
            }
        }
    }
}

fn heuristic_cost(mol: &MolGraph, stock: &StockCatalog) -> f64 {
    let sm = stock.lookup(mol);
    match sm.tier {
        StockTier::InCatalog => {
            let price = sm.entry.map_or(5.0, |e| e.est_usd_per_g);
            0.05 * (1.0 + price).ln()
        }
        StockTier::SimpleFragmentEst => 0.35,
        StockTier::NotInStock => {
            let ha = mol.heavy_atom_count() as f64;
            let rings = mol.rings.len() as f64;
            let chiral = mol
                .atoms
                .iter()
                .filter(|a| a.chiral != ChiralTag::None)
                .count() as f64;
            1.20 + 0.45 * rings + 0.08 * ha + 0.30 * chiral
        }
    }
}

struct SearchCtx<'a> {
    stock: &'a StockCatalog,
    cfg: &'a RetroConfig,
    nodes_expanded: usize,
    memo: HashMap<(String, usize), Vec<AndOrTree>>,
}

impl<'a> SearchCtx<'a> {
    fn solve(
        &mut self,
        mol: &MolGraph,
        depth_left: usize,
        ancestors: &mut HashSet<String>,
        is_root: bool,
    ) -> Vec<AndOrTree> {
        let key = mol.structural_key();
        let smi = mol.canonical_smiles();
        let stock_match = self.stock.lookup(mol);

        // Non-root buyable stock molecules terminate immediately (Buy-Before-Make!)
        if !is_root && stock_match.tier.is_buyable() {
            let cost = heuristic_cost(mol, self.stock);
            return vec![AndOrTree::Leaf {
                smiles: smi,
                mol: mol.clone(),
                stock: stock_match,
                cost,
            }];
        }

        if depth_left == 0 || self.nodes_expanded >= self.cfg.max_expansions {
            let cost = heuristic_cost(mol, self.stock);
            return vec![AndOrTree::Leaf {
                smiles: smi,
                mol: mol.clone(),
                stock: stock_match,
                cost,
            }];
        }

        let memo_key = (key.clone(), depth_left);
        if !is_root && let Some(cached) = self.memo.get(&memo_key) {
            return cached.clone();
        }

        self.nodes_expanded += 1;
        ancestors.insert(key.clone());

        let mut disconnections = enumerate_disconnections(mol);
        // Retro* value-priority ordering: V(r) = step_cost + sum h(precursor)
        let mut ranked_discs: Vec<(f64, Disconnection)> = disconnections
            .drain(..)
            .filter(|d| {
                // Reject cycles back to an ancestor
                !d.precursors
                    .iter()
                    .any(|p| ancestors.contains(&p.structural_key()))
            })
            .map(|d| {
                let child_h: f64 = d
                    .precursors
                    .iter()
                    .map(|p| heuristic_cost(p, self.stock))
                    .sum();
                (d.step_cost + child_h, d)
            })
            .collect();
        ranked_discs.sort_by(|a, b| a.0.total_cmp(&b.0));
        let branch_cap = if is_root {
            self.cfg.branch_limit + 4
        } else {
            self.cfg.branch_limit
        };
        ranked_discs.truncate(branch_cap);

        let mut candidates: Vec<AndOrTree> = Vec::new();
        for (_, disc) in ranked_discs {
            if disc.precursors.is_empty() {
                continue;
            }
            // Solve each precursor recursively (keeping top 2 distinct sub-trees per child)
            let mut child_options: Vec<Vec<AndOrTree>> = Vec::with_capacity(disc.precursors.len());
            for p in &disc.precursors {
                let mut sub = self.solve(p, depth_left - 1, ancestors, false);
                sub.truncate(2);
                if sub.is_empty() {
                    let sm = self.stock.lookup(p);
                    sub.push(AndOrTree::Leaf {
                        smiles: p.canonical_smiles(),
                        mol: p.clone(),
                        cost: heuristic_cost(p, self.stock),
                        stock: sm,
                    });
                }
                child_options.push(sub);
            }

            // Combine child sub-trees (1 precursor or 2 precursors)
            let combos: Vec<Vec<AndOrTree>> = if child_options.len() == 1 {
                child_options[0].iter().cloned().map(|c| vec![c]).collect()
            } else {
                let mut v = Vec::new();
                for c0 in &child_options[0] {
                    for c1 in &child_options[1] {
                        v.push(vec![c0.clone(), c1.clone()]);
                    }
                }
                v
            };

            for ch in combos {
                let total_steps = 1 + ch.iter().map(AndOrTree::total_steps).sum::<usize>();
                let lls = 1 + ch.iter().map(AndOrTree::lls).max().unwrap_or(0);
                let unresolved_leaves = ch.iter().map(AndOrTree::unresolved_leaves).sum::<usize>();
                let total_leaves = ch.iter().map(AndOrTree::total_leaves).sum::<usize>();
                let overall_yield =
                    disc.est_yield * ch.iter().map(AndOrTree::overall_yield).product::<f64>();
                let warning_count =
                    disc.warnings.len() + ch.iter().map(AndOrTree::warning_count).sum::<usize>();
                let total_cost = disc.step_cost
                    + ch.iter().map(AndOrTree::total_cost).sum::<f64>()
                    + 0.10 * (lls as f64);

                candidates.push(AndOrTree::Step {
                    product_smiles: smi.clone(),
                    product_mol: mol.clone(),
                    disc: disc.clone(),
                    children: ch,
                    total_cost,
                    total_steps,
                    lls,
                    unresolved_leaves,
                    total_leaves,
                    overall_yield,
                    warning_count,
                });
            }
        }

        ancestors.remove(&key);

        // Sort best routes first: fewest unresolved leaves, then lowest Retro* total_cost
        candidates.sort_by(|a, b| {
            a.unresolved_leaves()
                .cmp(&b.unresolved_leaves())
                .then_with(|| a.total_cost().total_cmp(&b.total_cost()))
                .then_with(|| a.total_steps().cmp(&b.total_steps()))
        });

        // Deduplicate by tree signature
        let mut unique = Vec::new();
        let mut sigs = HashSet::new();
        for c in candidates {
            if sigs.insert(c.signature()) {
                unique.push(c);
                if unique.len() >= if is_root { self.cfg.max_routes } else { 3 } {
                    break;
                }
            }
        }

        // If no disconnection succeeded on a non-root node, keep it as an unsolved leaf
        if !is_root && unique.is_empty() {
            unique.push(AndOrTree::Leaf {
                smiles: smi,
                mol: mol.clone(),
                cost: heuristic_cost(mol, self.stock),
                stock: stock_match,
            });
        }

        if !is_root {
            self.memo.insert(memo_key, unique.clone());
        }
        unique
    }
}

/// Linearizes an AND-OR tree into ordered forward synthesis steps + exact backward-propagated
/// stoichiometry (mg / mmol) + deduplicated Buy-Before-Make shopping list.
fn materialize_route(
    rank: usize,
    tree: &AndOrTree,
    target_scale_mg: f64,
) -> Option<SynthesisRoute> {
    let mut steps = Vec::new();
    let mut shopping: HashMap<String, StartingMaterial> = HashMap::new();

    fn walk(
        node: &AndOrTree,
        req_product_mg: f64,
        steps: &mut Vec<SynthesisStep>,
        shopping: &mut HashMap<String, StartingMaterial>,
    ) {
        match node {
            AndOrTree::Leaf {
                smiles, mol, stock, ..
            } => {
                let entry = shopping
                    .entry(smiles.clone())
                    .or_insert_with(|| StartingMaterial {
                        canonical_smiles: smiles.clone(),
                        name: stock.entry.as_ref().map_or_else(
                            || format!("Fragment ({})", mol.formula()),
                            |e| e.name.clone(),
                        ),
                        formula: mol.formula(),
                        mw: round2(mol.molecular_weight()),
                        tier: stock.tier,
                        catalog_id: stock.entry.as_ref().map(|e| e.catalog_id.clone()),
                        supplier: stock.entry.as_ref().map(|e| e.supplier.clone()),
                        est_usd_per_g: stock.entry.as_ref().map(|e| e.est_usd_per_g),
                        required_mg: 0.0,
                    });
                entry.required_mg = round1(entry.required_mg + req_product_mg);
            }
            AndOrTree::Step {
                product_smiles,
                product_mol,
                disc,
                children,
                ..
            } => {
                let prod_mw = product_mol.molecular_weight().max(1.0);
                let prod_mmol = req_product_mg / prod_mw;
                let req_limit_mmol = prod_mmol / disc.est_yield.clamp(0.15, 0.98);

                // Compute reactant masses (1.0 eq limiting + 1.15 eq partner) & walk children first (forward order!)
                let mut rxn = ReactionScheme::new(&disc.reaction_name);
                rxn.conditions = format!(
                    "{} | {} | {}",
                    disc.reagents, disc.solvent, disc.temperature
                );
                rxn.products.push(product_mol.clone());

                let mut reactant_masses = Vec::with_capacity(children.len());
                let mut reactant_smis = Vec::with_capacity(children.len());
                for (ci, child) in children.iter().enumerate() {
                    let (csmi, cmol) = match child {
                        AndOrTree::Leaf { smiles, mol, .. } => (smiles.clone(), mol.clone()),
                        AndOrTree::Step {
                            product_smiles,
                            product_mol,
                            ..
                        } => (product_smiles.clone(), product_mol.clone()),
                    };
                    let eq = if ci == 0 { 1.0 } else { 1.15 };
                    let r_mmol = req_limit_mmol * eq;
                    let r_mg = r_mmol * cmol.molecular_weight().max(1.0);
                    walk(child, r_mg, steps, shopping);
                    rxn.reactants.push(cmol);
                    reactant_smis.push(csmi.clone());
                    reactant_masses.push((csmi, round3(r_mmol), round1(r_mg)));
                }

                let stoich = compute_stoichiometry(&rxn, req_limit_mmol, Some(req_product_mg));
                // Include implicit stoichiometric reagent mass (CO, H2, H2O, HN3, CDI, [O])
                // for FGI / cycloaddition steps so atom economy stays physically bounded in [0, 100]%.
                let implicit_reagent_mw = match disc.template_id.as_str() {
                    "BENZYLIC_CARBONYLATION" => 28.010,      // CO
                    "KETONE_REDUCTION" => 2.016,             // H2 / hydride equiv.
                    "NITRILE_HYDROLYSIS" => 36.030,          // 2 H2O -> COOH + NH3
                    "TETRAZOLE_CLICK" => 43.028,             // HN3 (from NaN3)
                    "SULFONE_OX" => 31.998,                  // 2 [O]
                    "UREA_CDI" | "CARBAMATE_CDI" => 162.150, // CDI carbonyl donor
                    _ => 0.0,
                };
                let total_r_mw = rxn
                    .reactants
                    .iter()
                    .map(MolGraph::molecular_weight)
                    .sum::<f64>()
                    + implicit_reagent_mw;
                let ae_pct = if total_r_mw > 0.0 {
                    (100.0 * prod_mw / total_r_mw).clamp(0.0, 100.0)
                } else {
                    stoich.atom_economy_pct.clamp(0.0, 100.0)
                };
                let step_no = steps.len() + 1;
                steps.push(SynthesisStep {
                    step_number: step_no,
                    template_id: disc.template_id.clone(),
                    reaction_name: disc.reaction_name.clone(),
                    reaction_class: disc.reaction_class.clone(),
                    reactant_smiles: reactant_smis,
                    product_smiles: product_smiles.clone(),
                    reaction_smiles: rxn.reaction_smiles(),
                    reagents: disc.reagents.clone(),
                    solvent: disc.solvent.clone(),
                    temperature: disc.temperature.clone(),
                    est_yield_pct: round1(disc.est_yield * 100.0),
                    atom_economy_pct: round1(ae_pct),
                    product_scale_mg: round1(req_product_mg),
                    reactant_masses_mg: reactant_masses,
                    workup_notes: disc.workup_notes.clone(),
                    warnings: disc.warnings.clone(),
                });
            }
        }
    }

    walk(tree, target_scale_mg.max(1.0), &mut steps, &mut shopping);
    if steps.is_empty() {
        return None;
    }

    let mut starting_materials: Vec<StartingMaterial> = shopping.into_values().collect();
    starting_materials.sort_by(|a, b| a.canonical_smiles.cmp(&b.canonical_smiles));

    let total_leaves = tree.total_leaves().max(1);
    let resolved_leaves = total_leaves.saturating_sub(tree.unresolved_leaves());
    let stock_coverage_pct = round1(100.0 * (resolved_leaves as f64) / (total_leaves as f64));
    let mean_ae =
        round1(steps.iter().map(|s| s.atom_economy_pct).sum::<f64>() / (steps.len() as f64));
    let overall_yield_pct = round1(tree.overall_yield() * 100.0);

    // Composite score in [0, 100]: rewards stock coverage, yield, short LLS, atom economy, few PG alerts
    let score = (0.40 * stock_coverage_pct
        + 0.30 * overall_yield_pct
        + 0.15 * mean_ae
        + 0.15 * (100.0 / (tree.lls().max(1) as f64)).min(100.0)
        - 4.0 * (tree.warning_count() as f64))
        .clamp(0.0, 100.0);

    Some(SynthesisRoute {
        rank,
        total_steps: tree.total_steps(),
        longest_linear_seq: tree.lls(),
        fully_resolved_in_stock: tree.unresolved_leaves() == 0,
        stock_coverage_pct,
        overall_est_yield_pct: overall_yield_pct,
        mean_atom_economy_pct: mean_ae,
        composite_score: round1(score),
        warning_count: tree.warning_count(),
        steps,
        starting_materials,
    })
}

/// Plans multi-step retrosynthetic routes for `target` against `stock` using Retro* AND-OR search.
#[must_use]
pub fn plan_retrosynthesis(
    target: &MolGraph,
    stock: &StockCatalog,
    cfg: &RetroConfig,
) -> SynthesisPlan {
    let t0 = Instant::now();
    let target_name = target
        .name
        .clone()
        .unwrap_or_else(|| "Target Molecule".to_string());
    let target_smiles = target.canonical_smiles();
    let target_formula = target.formula();
    let target_mw = round2(target.molecular_weight());
    let target_stock_match = stock.lookup(target);
    let nearest_commercial_analogs = stock.nearest_analogs(target, 5);

    let mut ctx = SearchCtx {
        stock,
        cfg,
        nodes_expanded: 0,
        memo: HashMap::new(),
    };
    let mut ancestors = HashSet::new();
    let trees = ctx.solve(target, cfg.max_steps.clamp(1, 6), &mut ancestors, true);

    let mut routes: Vec<SynthesisRoute> = trees
        .iter()
        .enumerate()
        .filter_map(|(i, t)| materialize_route(i + 1, t, cfg.target_scale_mg))
        .collect();
    routes.sort_by(|a, b| {
        b.fully_resolved_in_stock
            .cmp(&a.fully_resolved_in_stock)
            .then_with(|| b.composite_score.total_cmp(&a.composite_score))
            .then_with(|| a.total_steps.cmp(&b.total_steps))
    });
    for (i, r) in routes.iter_mut().enumerate() {
        r.rank = i + 1;
    }

    SynthesisPlan {
        target_name,
        target_smiles,
        target_formula,
        target_mw,
        target_scale_mg: cfg.target_scale_mg,
        target_stock_match,
        nearest_commercial_analogs,
        catalog_size: stock.len(),
        nodes_expanded: ctx.nodes_expanded,
        elapsed_ms: t0.elapsed().as_millis(),
        routes,
    }
}

#[inline]
fn round1(v: f64) -> f64 {
    (v * 10.0).round() / 10.0
}

#[inline]
fn round2(v: f64) -> f64 {
    (v * 100.0).round() / 100.0
}

#[inline]
fn round3(v: f64) -> f64 {
    (v * 1000.0).round() / 1000.0
}
