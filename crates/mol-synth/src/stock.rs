#![forbid(unsafe_code)]

use std::collections::HashMap;
use std::path::Path;
use std::sync::OnceLock;

use mol_graph::{MolGraph, parse_smiles};
use mol_sar::MorganFp;
use serde::{Deserialize, Serialize};

use crate::data_dir::building_blocks_dir;

/// Commercial availability tier for a starting material or target compound.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StockTier {
    /// Verified exact match in the local commercial building-block catalog.
    InCatalog,
    /// Small, simple organic fragment (<= 10 heavy atoms, <= 1 standard ring) estimated commercially available.
    SimpleFragmentEst,
    /// Requires multi-step synthesis (not found in stock).
    NotInStock,
}

impl StockTier {
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::InCatalog => "In Stock (Catalog)",
            Self::SimpleFragmentEst => "Est. Buyable (Simple)",
            Self::NotInStock => "Synthesize",
        }
    }

    #[must_use]
    pub fn is_buyable(self) -> bool {
        matches!(self, Self::InCatalog | Self::SimpleFragmentEst)
    }
}

/// Commercial catalog entry for a building block.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StockEntry {
    pub catalog_id: String,
    pub name: String,
    pub canonical_smiles: String,
    pub structural_key: String,
    pub supplier: String,
    pub est_usd_per_g: f64,
    pub mw: f64,
    pub heavy_atoms: usize,
    #[serde(skip)]
    pub fp: Option<MorganFp>,
}

/// Nearest buyable commercial analog found via ECFP4 Tanimoto similarity (Buy-Before-Make).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StockAnalog {
    pub catalog_id: String,
    pub name: String,
    pub canonical_smiles: String,
    pub supplier: String,
    pub est_usd_per_g: f64,
    pub tanimoto: f64,
}

/// Result of querying a molecule against the commercial building-block catalog.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StockMatch {
    pub tier: StockTier,
    pub entry: Option<StockEntry>,
}

/// In-memory indexed commercial building-block catalog (built-in MedChem core + `CRAFT_DATA_DIR/building_blocks`).
#[derive(Debug, Clone)]
pub struct StockCatalog {
    by_key: HashMap<String, usize>,
    by_smiles: HashMap<String, usize>,
    entries: Vec<StockEntry>,
}

impl StockCatalog {
    /// Builds the curated built-in MedChem building block catalog and merges any `.tsv`/`.smi` files
    /// found in `CRAFT_DATA_DIR/building_blocks`.
    #[must_use]
    pub fn load_default() -> Self {
        let mut cat = Self::builtin_medchem();
        let dir = building_blocks_dir();
        if let Ok(rd) = std::fs::read_dir(&dir) {
            let mut paths: Vec<_> = rd
                .filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|p| {
                    p.is_file()
                        && p.extension().and_then(|s| s.to_str()).is_some_and(|e| {
                            matches!(
                                e.to_ascii_lowercase().as_str(),
                                "tsv" | "smi" | "ism" | "csv"
                            )
                        })
                })
                .collect();
            paths.sort();
            for p in paths {
                let _ = cat.merge_file(&p);
            }
        }
        cat
    }

    /// Number of indexed commercial building blocks.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Lookup exact commercial stock match or simple-fragment commerciality.
    #[must_use]
    pub fn lookup(&self, mol: &MolGraph) -> StockMatch {
        let key = mol.structural_key();
        if let Some(&idx) = self.by_key.get(&key) {
            return StockMatch {
                tier: StockTier::InCatalog,
                entry: self.entries.get(idx).cloned(),
            };
        }
        let smi = mol.canonical_smiles();
        if let Some(&idx) = self.by_smiles.get(&smi) {
            return StockMatch {
                tier: StockTier::InCatalog,
                entry: self.entries.get(idx).cloned(),
            };
        }
        if is_simple_buyable_fragment(mol) {
            return StockMatch {
                tier: StockTier::SimpleFragmentEst,
                entry: None,
            };
        }
        StockMatch {
            tier: StockTier::NotInStock,
            entry: None,
        }
    }

    /// Finds the top-K structurally nearest commercial building blocks via ECFP4 Tanimoto.
    #[must_use]
    pub fn nearest_analogs(&self, mol: &MolGraph, top_k: usize) -> Vec<StockAnalog> {
        let qfp = MorganFp::ecfp4(mol);
        let mut scored: Vec<(f64, &StockEntry)> = self
            .entries
            .iter()
            .filter_map(|e| {
                let fp = e.fp.as_ref()?;
                let sim = qfp.tanimoto(fp);
                if sim >= 0.25 { Some((sim, e)) } else { None }
            })
            .collect();
        scored.sort_by(|a, b| b.0.total_cmp(&a.0));
        scored
            .into_iter()
            .take(top_k)
            .map(|(sim, e)| StockAnalog {
                catalog_id: e.catalog_id.clone(),
                name: e.name.clone(),
                canonical_smiles: e.canonical_smiles.clone(),
                supplier: e.supplier.clone(),
                est_usd_per_g: e.est_usd_per_g,
                tanimoto: (sim * 1000.0).round() / 1000.0,
            })
            .collect()
    }

    /// Adds one building block parsed from SMILES (deduplicating by `structural_key`).
    pub fn insert_smiles(
        &mut self,
        catalog_id: &str,
        name: &str,
        smiles: &str,
        supplier: &str,
        est_usd_per_g: f64,
    ) {
        let Ok(mol) = parse_smiles(smiles) else {
            return;
        };
        if mol.atoms.is_empty() {
            return;
        }
        let key = mol.structural_key();
        if self.by_key.contains_key(&key) {
            return;
        }
        let canon = mol.canonical_smiles();
        let entry = StockEntry {
            catalog_id: catalog_id.to_string(),
            name: name.to_string(),
            canonical_smiles: canon.clone(),
            structural_key: key.clone(),
            supplier: supplier.to_string(),
            est_usd_per_g,
            mw: (mol.molecular_weight() * 100.0).round() / 100.0,
            heavy_atoms: mol.heavy_atom_count(),
            fp: Some(MorganFp::ecfp4(&mol)),
        };
        let idx = self.entries.len();
        self.by_key.insert(key, idx);
        self.by_smiles.insert(canon, idx);
        self.entries.push(entry);
    }

    /// Merges a `.tsv` (`catalog_id\tname\tsmiles\tsupplier\tusd_per_g`) or `.smi` (`SMILES NAME`)
    /// catalog file into this in-memory stock index.
    pub fn merge_file(&mut self, path: &Path) -> std::io::Result<usize> {
        let raw = std::fs::read_to_string(path)?;
        let before = self.entries.len();
        let is_tsv = path
            .extension()
            .and_then(|s| s.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("tsv"));
        for (line_no, line) in raw.lines().enumerate() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if is_tsv {
                let cols: Vec<&str> = trimmed.split('\t').collect();
                if line_no == 0
                    && cols.first().is_some_and(|c| {
                        c.eq_ignore_ascii_case("catalog_id") || c.eq_ignore_ascii_case("id")
                    })
                {
                    continue;
                }
                if cols.len() >= 3 {
                    let cid = cols[0].trim();
                    let nm = cols[1].trim();
                    let smi = cols[2].trim();
                    let sup = cols.get(3).map_or("LocalStock", |s| s.trim());
                    let price = cols
                        .get(4)
                        .and_then(|s| s.trim().parse::<f64>().ok())
                        .unwrap_or(8.0);
                    self.insert_smiles(cid, nm, smi, sup, price);
                }
            } else {
                let mut parts = trimmed.split_whitespace();
                let Some(smi) = parts.next() else {
                    continue;
                };
                let nm = parts.next().unwrap_or("LocalBB");
                let cid = format!("LOC-{}", line_no + 1);
                self.insert_smiles(&cid, nm, smi, "LocalSMI", 10.0);
            }
        }
        Ok(self.entries.len().saturating_sub(before))
    }

    /// Writes the current catalog to a clean editable TSV file (for `CRAFT_DATA_DIR/building_blocks`).
    pub fn write_tsv(&self, path: &Path) -> std::io::Result<()> {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let mut buf = String::from("catalog_id\tname\tsmiles\tsupplier\tusd_per_g\n");
        for e in &self.entries {
            buf.push_str(&format!(
                "{}\t{}\t{}\t{}\t{:.2}\n",
                e.catalog_id, e.name, e.canonical_smiles, e.supplier, e.est_usd_per_g
            ));
        }
        std::fs::write(path, buf)
    }

    /// Curated core medicinal chemistry commercial building block library (~170 high-frequency
    /// commodity & MedChem reagents spanning acids, amines, aryl/heteroaryl halides, boronic acids,
    /// heterocycles, PEG linkers, sulfonyl chlorides, aldehydes, ketones, phenols & alcohols).
    #[must_use]
    pub fn builtin_medchem() -> Self {
        let mut cat = Self {
            by_key: HashMap::new(),
            by_smiles: HashMap::new(),
            entries: Vec::new(),
        };
        for &(cid, name, smi, price) in BUILTIN_BUILDING_BLOCKS {
            cat.insert_smiles(cid, name, smi, "MedChemCore", price);
        }
        cat
    }
}

static DEFAULT_STOCK_CACHE: OnceLock<StockCatalog> = OnceLock::new();

/// Returns a lazily cached reference to the default building-block catalog.
#[must_use]
pub fn default_stock_catalog() -> &'static StockCatalog {
    DEFAULT_STOCK_CACHE.get_or_init(StockCatalog::load_default)
}

/// Conservative simple-fragment heuristic for small commodity organic building blocks
/// (`<= 10` heavy atoms, `<= 1` standard 5/6-membered ring, standard organic elements only).
fn is_simple_buyable_fragment(mol: &MolGraph) -> bool {
    let ha = mol.heavy_atom_count();
    if ha == 0 || ha > 10 {
        return false;
    }
    if mol.rings.len() > 1 {
        return false;
    }
    if mol.rings.iter().any(|r| r.len() < 5 || r.len() > 6) {
        return false;
    }
    let chiral_cnt = mol
        .atoms
        .iter()
        .filter(|a| a.chiral != mol_graph::ChiralTag::None)
        .count();
    if chiral_cnt > 1 {
        return false;
    }
    mol.atoms.iter().all(|a| {
        matches!(
            a.atomic_number,
            1 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
        )
    })
}

/// Built-in curated MedChem commercial building blocks (`(catalog_id, name, smiles, est_usd_per_g)`).
const BUILTIN_BUILDING_BLOCKS: &[(&str, &str, &str, f64)] = &[
    // --- Carboxylic acids & acyl chlorides ---
    ("BB-001", "Acetic acid", "CC(=O)O", 0.2),
    ("BB-002", "Propionic acid", "CCC(=O)O", 0.3),
    ("BB-003", "Isobutyric acid", "CC(C)C(=O)O", 0.5),
    ("BB-004", "Pivalic acid", "CC(C)(C)C(=O)O", 0.6),
    ("BB-005", "Cyclopropanecarboxylic acid", "O=C(O)C1CC1", 1.2),
    (
        "BB-006",
        "Cyclohexanecarboxylic acid",
        "O=C(O)C1CCCCC1",
        0.8,
    ),
    ("BB-007", "Benzoic acid", "O=C(O)c1ccccc1", 0.3),
    ("BB-008", "4-Methylbenzoic acid", "Cc1ccc(C(=O)O)cc1", 0.6),
    ("BB-009", "4-Fluorobenzoic acid", "O=C(O)c1ccc(F)cc1", 0.9),
    ("BB-010", "4-Chlorobenzoic acid", "O=C(O)c1ccc(Cl)cc1", 0.7),
    ("BB-011", "4-Bromobenzoic acid", "O=C(O)c1ccc(Br)cc1", 1.0),
    ("BB-012", "4-Methoxybenzoic acid", "COc1ccc(C(=O)O)cc1", 0.7),
    (
        "BB-013",
        "4-(Trifluoromethyl)benzoic acid",
        "O=C(O)c1ccc(C(F)(F)F)cc1",
        1.8,
    ),
    ("BB-014", "Salicylic acid", "O=C(O)c1ccccc1O", 0.4),
    ("BB-015", "Phenylacetic acid", "O=C(O)Cc1ccccc1", 0.5),
    (
        "BB-016",
        "2-(4-Isobutylphenyl)propanoic acid (Ibuprofen)",
        "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        1.5,
    ),
    ("BB-017", "Nicotinic acid", "O=C(O)c1cccnc1", 0.6),
    ("BB-018", "Isonicotinic acid", "O=C(O)c1ccncc1", 0.7),
    ("BB-019", "Picolinic acid", "O=C(O)c1ccccn1", 0.8),
    (
        "BB-020",
        "Pyrazine-2-carboxylic acid",
        "O=C(O)c1cnccn1",
        1.5,
    ),
    (
        "BB-021",
        "Indole-3-carboxylic acid",
        "O=C(O)c1c[nH]c2ccccc12",
        2.5,
    ),
    ("BB-022", "1-Naphthoic acid", "O=C(O)c1cccc2ccccc12", 1.4),
    ("BB-023", "2-Naphthoic acid", "O=C(O)c1ccc2ccccc2c1", 1.4),
    ("BB-024", "Acetyl chloride", "CC(=O)Cl", 0.4),
    ("BB-025", "Benzoyl chloride", "O=C(Cl)c1ccccc1", 0.5),
    ("BB-026", "4-Aminobenzoic acid", "Nc1ccc(C(=O)O)cc1", 0.7),
    ("BB-027", "3-Aminobenzoic acid", "Nc1cccc(C(=O)O)cc1", 0.9),
    (
        "BB-028",
        "4-Nitrobenzoic acid",
        "O=C(O)c1ccc([N+](=O)[O-])cc1",
        0.8,
    ),
    ("BB-029", "4-Cyanobenzoic acid", "N#Cc1ccc(C(=O)O)cc1", 1.4),
    (
        "BB-030",
        "3,4-Dimethoxybenzoic acid",
        "COc1ccc(C(=O)O)cc1OC",
        1.1,
    ),
    // --- Primary & secondary amines, anilines, cyclic amines ---
    ("BB-031", "Ammonia", "N", 0.1),
    ("BB-032", "Methylamine", "CN", 0.3),
    ("BB-033", "Dimethylamine", "CNC", 0.3),
    ("BB-034", "Ethylamine", "CCN", 0.3),
    ("BB-035", "Diethylamine", "CCNCC", 0.3),
    ("BB-036", "Isopropylamine", "CC(C)N", 0.4),
    ("BB-037", "Cyclopropylamine", "NC1CC1", 0.9),
    ("BB-038", "Cyclohexylamine", "NC1CCCCC1", 0.5),
    ("BB-039", "Benzylamine", "NCc1ccccc1", 0.5),
    ("BB-040", "2-Phenylethylamine", "NCCc1ccccc1", 0.6),
    ("BB-041", "Aniline", "Nc1ccccc1", 0.3),
    ("BB-042", "p-Toluidine", "Cc1ccc(N)cc1", 0.5),
    ("BB-043", "4-Fluoroaniline", "Nc1ccc(F)cc1", 0.7),
    ("BB-044", "4-Chloroaniline", "Nc1ccc(Cl)cc1", 0.6),
    ("BB-045", "4-Bromoaniline", "Nc1ccc(Br)cc1", 0.9),
    ("BB-046", "p-Anisidine", "COc1ccc(N)cc1", 0.6),
    ("BB-047", "4-Aminophenol", "Nc1ccc(O)cc1", 0.5),
    (
        "BB-048",
        "4-(Trifluoromethyl)aniline",
        "Nc1ccc(C(F)(F)F)cc1",
        1.5,
    ),
    (
        "BB-049",
        "4-Aminobenzenesulfonamide (Sulfanilamide)",
        "NS(=O)(=O)c1ccc(N)cc1",
        0.9,
    ),
    ("BB-050", "3-Aminopyridine", "Nc1cccnc1", 1.2),
    ("BB-051", "4-Aminopyridine", "Nc1ccncc1", 1.0),
    ("BB-052", "2-Aminopyridine", "Nc1ccccn1", 0.8),
    ("BB-053", "2-Aminopyrimidine", "Nc1ncccn1", 1.2),
    ("BB-054", "Pyrrolidine", "C1CCNC1", 0.4),
    ("BB-055", "Piperidine", "C1CCNCC1", 0.4),
    ("BB-056", "Morpholine", "C1COCCN1", 0.4),
    ("BB-057", "Piperazine", "C1CNCCN1", 0.5),
    ("BB-058", "1-Methylpiperazine", "CN1CCNCC1", 0.6),
    ("BB-059", "1-Boc-piperazine", "CC(C)(C)OC(=O)N1CCNCC1", 1.8),
    ("BB-060", "4-Hydroxypiperidine", "OC1CCNCC1", 1.0),
    ("BB-061", "Azetidine", "C1CNC1", 2.5),
    ("BB-062", "2-Aminoethanol", "NCCO", 0.3),
    ("BB-063", "2-Methoxyethylamine", "COCCN", 0.6),
    ("BB-064", "2-(Dimethylamino)ethylamine", "CN(C)CCN", 0.7),
    ("BB-065", "3-(Dimethylamino)propylamine", "CN(C)CCCN", 0.6),
    ("BB-066", "4-Aminocyclohexanol", "NC1CCC(O)CC1", 2.2),
    ("BB-067", "2-Amino-4-methylpyridine", "Cc1ccnc(N)c1", 1.4),
    ("BB-068", "N-Methylaniline", "CNc1ccccc1", 0.5),
    // --- Aryl / heteroaryl halides ---
    ("BB-069", "Bromobenzene", "Brc1ccccc1", 0.4),
    ("BB-070", "Iodobenzene", "Ic1ccccc1", 0.8),
    ("BB-071", "4-Bromotoluene", "Cc1ccc(Br)cc1", 0.6),
    ("BB-072", "4-Bromoanisole", "COc1ccc(Br)cc1", 0.7),
    ("BB-073", "1-Bromo-4-fluorobenzene", "Fc1ccc(Br)cc1", 0.7),
    ("BB-074", "1-Bromo-4-chlorobenzene", "Clc1ccc(Br)cc1", 0.7),
    ("BB-075", "4-Bromobenzonitrile", "N#Cc1ccc(Br)cc1", 1.2),
    (
        "BB-076",
        "1-Bromo-4-(trifluoromethyl)benzene",
        "FC(F)(F)c1ccc(Br)cc1",
        1.4,
    ),
    ("BB-077", "4-Bromophenol", "Oc1ccc(Br)cc1", 0.7),
    ("BB-078", "2-Bromopyridine", "Brc1ccccn1", 0.9),
    ("BB-079", "3-Bromopyridine", "Brc1cccnc1", 1.0),
    ("BB-080", "4-Bromopyridine", "Brc1ccncc1", 1.5),
    ("BB-081", "2-Chloropyrimidine", "Clc1ncccn1", 1.4),
    ("BB-082", "2,4-Dichloropyrimidine", "Clc1ccnc(Cl)n1", 1.8),
    ("BB-083", "4,6-Dichloropyrimidine", "Clc1cc(Cl)ncn1", 1.9),
    ("BB-084", "6-Chloropurine", "Clc1ncnc2[nH]cnc12", 2.5),
    ("BB-085", "4-Chloroquinazoline", "Clc1ncnc2ccccc12", 3.0),
    ("BB-086", "2-Chloro-4-aminopyrimidine", "Nc1ccnc(Cl)n1", 2.8),
    ("BB-087", "2-Chloropyridine", "Clc1ccccn1", 0.7),
    ("BB-088", "4-Chloropyridine", "Clc1ccncc1", 1.2),
    ("BB-089", "2-Bromo-5-methylpyridine", "Cc1ccc(Br)nc1", 1.6),
    ("BB-090", "3-Bromoquinoline", "Brc1cnc2ccccc2c1", 2.4),
    ("BB-091", "1-Bromonaphthalene", "Brc1cccc2ccccc12", 0.9),
    ("BB-092", "2-Bromothiophene", "Brc1cccs1", 1.1),
    // --- Boronic acids (Suzuki-Miyaura partners) ---
    ("BB-093", "Phenylboronic acid", "OB(O)c1ccccc1", 0.9),
    ("BB-094", "4-Tolylboronic acid", "Cc1ccc(B(O)O)cc1", 1.1),
    (
        "BB-095",
        "4-Fluorophenylboronic acid",
        "OB(O)c1ccc(F)cc1",
        1.2,
    ),
    (
        "BB-096",
        "4-Chlorophenylboronic acid",
        "OB(O)c1ccc(Cl)cc1",
        1.1,
    ),
    (
        "BB-097",
        "4-Methoxyphenylboronic acid",
        "COc1ccc(B(O)O)cc1",
        1.1,
    ),
    (
        "BB-098",
        "4-Hydroxyphenylboronic acid",
        "OB(O)c1ccc(O)cc1",
        1.6,
    ),
    (
        "BB-099",
        "4-Cyanophenylboronic acid",
        "N#Cc1ccc(B(O)O)cc1",
        1.8,
    ),
    (
        "BB-100",
        "4-(Trifluoromethyl)phenylboronic acid",
        "OB(O)c1ccc(C(F)(F)F)cc1",
        1.9,
    ),
    ("BB-101", "3-Pyridinylboronic acid", "OB(O)c1cccnc1", 2.2),
    ("BB-102", "4-Pyridinylboronic acid", "OB(O)c1ccncc1", 2.4),
    ("BB-103", "Pyrimidin-5-ylboronic acid", "OB(O)c1cncnc1", 3.5),
    ("BB-104", "2-Thienylboronic acid", "OB(O)c1cccs1", 2.0),
    (
        "BB-105",
        "1H-Indol-5-ylboronic acid",
        "OB(O)c1ccc2[nH]ccc2c1",
        3.8,
    ),
    (
        "BB-106",
        "1H-Pyrazol-4-ylboronic acid",
        "OB(O)c1cn[nH]c1",
        3.5,
    ),
    // --- Phenols, alcohols, alkyl halides & PEG linkers ---
    ("BB-107", "Methanol", "CO", 0.1),
    ("BB-108", "Ethanol", "CCO", 0.1),
    ("BB-109", "Isopropanol", "CC(C)O", 0.1),
    ("BB-110", "Phenol", "Oc1ccccc1", 0.3),
    ("BB-111", "4-Methoxyphenol", "COc1ccc(O)cc1", 0.6),
    ("BB-112", "4-Fluorophenol", "Oc1ccc(F)cc1", 0.8),
    ("BB-113", "4-Chlorophenol", "Oc1ccc(Cl)cc1", 0.6),
    ("BB-114", "1-Naphthol", "Oc1cccc2ccccc12", 0.9),
    ("BB-115", "2-Naphthol", "Oc1ccc2ccccc2c1", 0.8),
    ("BB-116", "Benzyl alcohol", "OCc1ccccc1", 0.4),
    ("BB-117", "2-Phenylethanol", "OCCc1ccccc1", 0.5),
    ("BB-118", "Iodomethane", "CI", 0.5),
    ("BB-119", "Bromoethane", "CCBr", 0.4),
    ("BB-120", "1-Bromopropane", "CCCBr", 0.4),
    ("BB-121", "Benzyl bromide", "BrCc1ccccc1", 0.6),
    ("BB-122", "2-Bromoethyl methyl ether", "COCCBr", 1.1),
    ("BB-123", "1,2-Dibromoethane", "BrCCBr", 0.5),
    ("BB-124", "1,3-Dibromopropane", "BrCCCBr", 0.6),
    ("BB-125", "Bis(2-bromoethyl) ether", "BrCCOCCBr", 1.5),
    (
        "BB-126",
        "1,2-Bis(2-bromoethoxy)ethane",
        "BrCCOCCOCCBr",
        2.5,
    ),
    (
        "BB-127",
        "Tetraethylene glycol dibromide",
        "BrCCOCCOCCOCCBr",
        3.5,
    ),
    ("BB-128", "Triethylene glycol", "OCCOCCOCCO", 0.6),
    ("BB-129", "Tetraethylene glycol", "OCCOCCOCCOCCO", 0.8),
    (
        "BB-130",
        "2-Chloro-N,N-dimethylethylamine",
        "CN(C)CCCl",
        1.2,
    ),
    ("BB-131", "4-(2-Chloroethyl)morpholine", "ClCCN1CCOCC1", 1.5),
    // --- Aldehydes & ketones ---
    ("BB-132", "Formaldehyde", "C=O", 0.1),
    ("BB-133", "Acetaldehyde", "CC=O", 0.2),
    ("BB-134", "Acetone", "CC(C)=O", 0.1),
    ("BB-135", "Benzaldehyde", "O=Cc1ccccc1", 0.4),
    ("BB-136", "4-Fluorobenzaldehyde", "O=Cc1ccc(F)cc1", 0.8),
    ("BB-137", "4-Chlorobenzaldehyde", "O=Cc1ccc(Cl)cc1", 0.6),
    ("BB-138", "4-Methoxybenzaldehyde", "COc1ccc(C=O)cc1", 0.6),
    ("BB-139", "4-Hydroxybenzaldehyde", "O=Cc1ccc(O)cc1", 0.7),
    (
        "BB-140",
        "4-(Trifluoromethyl)benzaldehyde",
        "O=Cc1ccc(C(F)(F)F)cc1",
        1.5,
    ),
    ("BB-141", "Nicotinaldehyde", "O=Cc1cccnc1", 1.1),
    ("BB-142", "Isonicotinaldehyde", "O=Cc1ccncc1", 1.2),
    ("BB-143", "Acetophenone", "CC(=O)c1ccccc1", 0.4),
    ("BB-144", "Cyclohexanone", "O=C1CCCCC1", 0.4),
    (
        "BB-145",
        "N-Boc-4-piperidone",
        "CC(C)(C)OC(=O)N1CCC(=O)CC1",
        1.8,
    ),
    ("BB-146", "Tetrahydro-4H-pyran-4-one", "O=C1CCOCC1", 2.0),
    // --- Sulfonyl chlorides, heterocycles, indoles, alkynes & azides ---
    ("BB-147", "Methanesulfonyl chloride", "CS(=O)(=O)Cl", 0.5),
    (
        "BB-148",
        "Benzenesulfonyl chloride",
        "O=S(=O)(Cl)c1ccccc1",
        0.6,
    ),
    (
        "BB-149",
        "p-Toluenesulfonyl chloride",
        "Cc1ccc(S(=O)(=O)Cl)cc1",
        0.5,
    ),
    (
        "BB-150",
        "4-Fluorobenzenesulfonyl chloride",
        "O=S(=O)(Cl)c1ccc(F)cc1",
        1.2,
    ),
    (
        "BB-151",
        "4-Methoxybenzenesulfonyl chloride",
        "COc1ccc(S(=O)(=O)Cl)cc1",
        1.1,
    ),
    ("BB-152", "Indole", "c1ccc2[nH]ccc2c1", 0.9),
    ("BB-153", "1H-Indazole", "c1ccc2[nH]ncc2c1", 2.2),
    ("BB-154", "Benzimidazole", "c1ccc2[nH]cnc2c1", 1.0),
    ("BB-155", "Imidazole", "c1c[nH]cn1", 0.4),
    ("BB-156", "Pyrazole", "c1cn[nH]c1", 0.9),
    (
        "BB-157",
        "3,4-Bis(1H-indol-3-yl)-1H-pyrrole-2,5-dione (Arcyriarubin A)",
        "O=C1NC(=O)C(c2c[nH]c3ccccc23)=C1c1c[nH]c2ccccc12",
        18.0,
    ),
    (
        "BB-157B",
        "3,4-Bis(1H-indol-3-yl)-1H-pyrrole-2,5-dione (Arcyriarubin A)",
        "c1ccc2c(c1)c(c[nH]2)c1c(c2c[nH]c3ccccc23)c([nH]c1=O)=O",
        18.0,
    ),
    (
        "BB-158",
        "3,4-Dichloromaleimide",
        "O=C1NC(=O)C(Cl)=C1Cl",
        3.5,
    ),
    ("BB-159", "Maleimide", "O=C1C=CC(=O)N1", 0.9),
    ("BB-160", "Phenylacetylene", "C#Cc1ccccc1", 1.0),
    ("BB-161", "Propargyl bromide", "C#CCBr", 0.9),
    ("BB-162", "Benzyl azide", "[N-]=[N+]=NCc1ccccc1", 3.0),
    ("BB-163", "Benzenethiol", "Sc1ccccc1", 0.7),
    ("BB-164", "4-Fluorobenzenethiol", "Fc1ccc(S)cc1", 1.4),
    ("BB-165", "Nitrobenzene", "O=[N+]([O-])c1ccccc1", 0.3),
    ("BB-166", "Benzonitrile", "N#Cc1ccccc1", 0.4),
    ("BB-167", "Styrene", "C=Cc1ccccc1", 0.3),
    ("BB-168", "4-Nitroaniline", "Nc1ccc([N+](=O)[O-])cc1", 0.6),
    ("BB-169", "2-Aminobenzimidazole", "Nc1nc2ccccc2[nH]1", 1.4),
    ("BB-170", "4-Aminobenzonitrile", "N#Cc1ccc(N)cc1", 1.2),
    ("BB-171", "Isobutylbenzene", "CC(C)Cc1ccccc1", 0.3),
    ("BB-172", "Acetyl chloride", "CC(=O)Cl", 0.3),
];
