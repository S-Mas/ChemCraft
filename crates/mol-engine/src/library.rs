//! Chemical library management, bounded streaming ingestion (`.sdf`, `.smi`, `.gz`, `.csv`, `.dwar`),
//! parallel *in silico* ADMET/SMARTS/ECFP4 screening funnel, physiological pH 7.4 preparation,
//! diversity picking (MaxMin), and 3D SDF export for DockCraft virtual screening.

use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use mol_admet_bridge::{AdmetCraftProvider, AdmetProvider, MoleculeEvaluation};
use mol_graph::{
    BondOrder, MolGraph, SmartsQuery, ensure_3d_coords, parse_csv_records, parse_dwar,
    parse_mol_block, parse_smiles, write_sdf_v2000_3d,
};
use mol_sar::MorganFp;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Safety cap for streaming decompressed library files (`128 MiB`).
pub const MAX_LIBRARY_UNCOMPRESSED_BYTES: u64 = 128 * 1024 * 1024;
/// Default interactive safety cap on parsed molecules per load.
pub const DEFAULT_MAX_LIBRARY_RECORDS: usize = 50_000;

#[derive(Debug, Error)]
pub enum LibraryError {
    #[error("library I/O error: {0}")]
    Io(String),
    #[error("library parse error: {0}")]
    Parse(String),
    #[error("invalid filter query: {0}")]
    InvalidQuery(String),
    #[error("empty library or no molecules passed the filter")]
    EmptyResult,
}

/// Built-in curated chemical screening libraries embedded cleanly offline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum BuiltinLibraryPreset {
    FdaApprovedMedicines,
    FragmentRo3,
    KinaseAndGpcrScaffolds,
}

impl BuiltinLibraryPreset {
    pub fn label(self) -> &'static str {
        match self {
            Self::FdaApprovedMedicines => "FDA Approved Medicines (Curated)",
            Self::FragmentRo3 => "Fragment Ro3 Library (<300 Da)",
            Self::KinaseAndGpcrScaffolds => "Kinase & GPCR Scaffolds",
        }
    }

    pub fn filename_stem(self) -> &'static str {
        match self {
            Self::FdaApprovedMedicines => "fda_approved_curated",
            Self::FragmentRo3 => "fragment_ro3_curated",
            Self::KinaseAndGpcrScaffolds => "kinase_gpcr_scaffolds",
        }
    }

    pub fn from_str_lossy(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "fda" | "fda_approved" | "drugs" | "approved" => Some(Self::FdaApprovedMedicines),
            "frag" | "fragment" | "fragments" | "ro3" => Some(Self::FragmentRo3),
            "kinase" | "gpcr" | "scaffolds" => Some(Self::KinaseAndGpcrScaffolds),
            _ => None,
        }
    }
}

/// Medicinal chemistry pre-filter funnel before 3D docking.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum ScreeningFunnel {
    All,
    #[default]
    RuleOf5,
    LeadLike,
    FragmentRo3,
    BbbCns,
}

impl ScreeningFunnel {
    pub fn label(self) -> &'static str {
        match self {
            Self::All => "All Valid (No Property Filter)",
            Self::RuleOf5 => "Lipinski Ro5 + No PAINS",
            Self::LeadLike => "Lead-Like (MW 200-350, LogP -1..3.5, Clean MedChem)",
            Self::FragmentRo3 => "Fragment Ro3 (MW <= 300, LogP <= 3, TPSA <= 60)",
            Self::BbbCns => "CNS / BBB Permeant (MW <= 400, TPSA <= 90, BBB+)",
        }
    }

    pub fn from_str_lossy(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "all" | "none" => Some(Self::All),
            "ro5" | "lipinski" => Some(Self::RuleOf5),
            "lead" | "lead-like" | "lead_like" => Some(Self::LeadLike),
            "ro3" | "fragment" | "frag" => Some(Self::FragmentRo3),
            "bbb" | "cns" => Some(Self::BbbCns),
            _ => None,
        }
    }
}

/// Complete *in silico* library filtering and preparation configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibraryFilterConfig {
    pub funnel: ScreeningFunnel,
    pub smarts_substructure: Option<String>,
    pub reference_smiles: Option<String>,
    pub min_tanimoto: f64,
    pub diversity_pick_n: Option<usize>,
    pub max_hits: usize,
    pub protonate_ph74: bool,
    pub generate_3d: bool,
}

impl Default for LibraryFilterConfig {
    fn default() -> Self {
        Self {
            funnel: ScreeningFunnel::RuleOf5,
            smarts_substructure: None,
            reference_smiles: None,
            min_tanimoto: 0.0,
            diversity_pick_n: None,
            max_hits: 500,
            protonate_ph74: true,
            generate_3d: true,
        }
    }
}

/// Screened library entry ready for inspection, SAR import, 3D viewing, or DockCraft screening.
#[derive(Debug, Clone)]
pub struct ScreenedLibraryHit {
    pub name: String,
    pub canonical_smiles: String,
    pub mol: MolGraph,
    pub eval: MoleculeEvaluation,
    pub tanimoto_to_ref: Option<f64>,
    pub fp: MorganFp,
}

/// Report returned by a parallel library filtering + 3D prep run.
#[derive(Debug, Clone)]
pub struct LibraryScreeningResult {
    pub source_label: String,
    pub total_loaded: usize,
    pub passed_funnel: usize,
    pub hits: Vec<ScreenedLibraryHit>,
    pub elapsed_ms: u128,
}

/// Default directory for user chemical screening libraries (`~/.local/share/chemcraft/libraries`).
pub fn default_library_dir() -> PathBuf {
    if let Ok(xdg_data) = std::env::var("XDG_DATA_HOME")
        && !xdg_data.trim().is_empty()
    {
        return PathBuf::from(xdg_data).join("chemcraft").join("libraries");
    }
    if let Ok(home) = std::env::var("HOME")
        && !home.trim().is_empty()
    {
        return PathBuf::from(home)
            .join(".local")
            .join("share")
            .join("chemcraft")
            .join("libraries");
    }
    PathBuf::from("chemcraft_libraries")
}

/// Discover `.sdf`, `.sdf.gz`, `.smi`, `.smi.gz`, `.csv`, `.dwar` files inside `dir`.
pub fn discover_local_libraries(dir: &Path) -> Vec<PathBuf> {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for ent in entries.flatten() {
        let p = ent.path();
        if !p.is_file() {
            continue;
        }
        let name = ent.file_name().to_string_lossy().to_ascii_lowercase();
        if name.ends_with(".sdf")
            || name.ends_with(".sd")
            || name.ends_with(".sdf.gz")
            || name.ends_with(".smi")
            || name.ends_with(".smiles")
            || name.ends_with(".smi.gz")
            || name.ends_with(".csv")
            || name.ends_with(".dwar")
        {
            out.push(p);
        }
    }
    out.sort();
    out
}

/// Load a built-in curated screening preset as parsed `MolGraph` records.
pub fn load_builtin_library(preset: BuiltinLibraryPreset) -> Vec<MolGraph> {
    let pairs: &[(&str, &str)] = match preset {
        BuiltinLibraryPreset::FdaApprovedMedicines => &[
            ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
            ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
            ("Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
            ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
            (
                "Celecoxib",
                "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1",
            ),
            (
                "Atorvastatin",
                "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
            ),
            (
                "Simvastatin",
                "CCC(C)(C)C(=O)OC1CC(C)C=C2C=CC(C)C(CCC3CC(O)CC(=O)O3)C12",
            ),
            (
                "Losartan",
                "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nnn[nH]2)cc1",
            ),
            ("Metformin", "CN(C)C(=N)NC(=N)N"),
            ("Omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"),
            (
                "Sildenafil",
                "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12",
            ),
            (
                "Imatinib",
                "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
            ),
            ("Erlotinib", "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"),
            (
                "Gefitinib",
                "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
            ),
            (
                "Sorafenib",
                "CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)cc2)ccn1",
            ),
            (
                "Ciprofloxacin",
                "O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O",
            ),
            (
                "Amoxicillin",
                "CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O",
            ),
            ("Fluconazole", "Oc(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"),
            (
                "Voriconazole",
                "C[C@@H](c1ncncc1F)[C@](O)(Cn1cncn1)c1ccc(F)cc1F",
            ),
            ("Oseltamivir", "CCOC(=O)C1=CC(OC(CC)CC)C(NC(=O)C)C(N)C1"),
            ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
            ("Fluoxetine", "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1"),
            ("Sertraline", "CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc21"),
            (
                "Haloperidol",
                "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1",
            ),
            ("Clozapine", "CN1CCN(C2=Nc3cc(Cl)ccc3Nc3ccccc32)CC1"),
            ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
            (
                "Rivaroxaban",
                "O=C1OCC(CNC(=O)c2ccc(Cl)s2)N1c1ccc(N2CCOCC2=O)cc1",
            ),
            ("Metoprolol", "COCCc1ccc(OCC(O)CNC(C)C)cc1"),
            (
                "Amlodipine",
                "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl",
            ),
            (
                "Dexamethasone",
                "CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO",
            ),
            ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
            ("Nicotine", "CN1CCC[C@H]1c1cccnc1"),
            ("Tamoxifen", "CCC(=C(c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1"),
            ("Lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
            ("Propranolol", "CC(C)NCC(O)COc1cccc2ccccc12"),
        ],
        BuiltinLibraryPreset::FragmentRo3 => &[
            ("Frag_Indole", "c1ccc2[nH]ccc2c1"),
            ("Frag_7Azaindole", "c1cnc2[nH]ccc2c1"),
            ("Frag_Benzimidazole", "c1ccc2[nH]cnc2c1"),
            ("Frag_Benzoxazole", "c1ccc2ocnc2c1"),
            ("Frag_Benzothiazole", "c1ccc2scnc2c1"),
            ("Frag_Quinoline", "c1ccc2ncccc2c1"),
            ("Frag_Quinazoline", "c1ccc2ncncc2c1"),
            ("Frag_2Aminopyridine", "Nc1ccccn1"),
            ("Frag_2Aminopyrimidine", "Nc1ncccn1"),
            ("Frag_Pyrazole", "c1cn[nH]c1"),
            ("Frag_3Methylpyrazole", "Cc1cc[nH]n1"),
            ("Frag_Isoxazole", "c1conc1"),
            ("Frag_124Oxadiazole", "c1ncno1"),
            ("Frag_Morpholine", "C1COCCN1"),
            ("Frag_Piperazine", "C1CNCCN1"),
            ("Frag_Piperidine", "C1CCNCC1"),
            ("Frag_Pyrrolidine", "C1CCNC1"),
            ("Frag_Tetrahydropyran", "C1CCOCC1"),
            ("Frag_CyclopropylAmine", "NC1CC1"),
            ("Frag_BenzoicAcid", "O=C(O)c1ccccc1"),
            ("Frag_4Fluorophenol", "Oc1ccc(F)cc1"),
            ("Frag_4Chloroaniline", "Nc1ccc(Cl)cc1"),
            ("Frag_BenzylAmine", "NCc1ccccc1"),
            ("Frag_Thiophene2Carboxamide", "NC(=O)c1cccs1"),
            ("Frag_Furan2CarboxylicAcid", "O=C(O)c1ccco1"),
            ("Frag_SulfonamideBenzene", "NS(=O)(=O)c1ccccc1"),
            ("Frag_TetrazolePhenyl", "c1ccc(-c2nnn[nH]2)cc1"),
            ("Frag_Nicotinamide", "NC(=O)c1cccnc1"),
            ("Frag_Indane1One", "O=C1CCc2ccccc21"),
            ("Frag_Chromane", "c1ccc2c(c1)CCCO2"),
        ],
        BuiltinLibraryPreset::KinaseAndGpcrScaffolds => &[
            ("Kin_4Aminoquinazoline", "Nc1ncnc2ccccc12"),
            ("Kin_Anilinoquinazoline", "c1ccc(Nc2ncnc3ccccc23)cc1"),
            ("Kin_4Phenylaminopyrimidine", "Nc1ccnc(Nc2ccccc2)n1"),
            ("Kin_Pyrrolo23dPyrimidine", "c1cnc2[nH]ccc2n1"),
            ("Kin_Pyrazolo34dPyrimidine", "c1nc2c([nH]n1)ncnc2"),
            ("Kin_3Aminopyrazole", "Nc1cc[nH]n1"),
            ("Kin_IndazoleCarboxamide", "NC(=O)c1nn[nH]c2ccccc12"),
            ("Kin_OxindoleCore", "O=C1Cc2ccccc2N1"),
            ("Kin_2Aminothiazole", "Nc1ncsc1"),
            ("Kin_Imidazo12aPyridine", "c1cccn2ccnc12"),
            ("Kin_DiphenylUrea", "O=C(Nc1ccccc1)Nc1ccccc1"),
            ("Kin_DiarylAmide", "O=C(Nc1ccc(F)cc1)c1ccc(C)cc1"),
            ("Gpcr_Tryptamine", "NCCc1c[nH]c2ccccc12"),
            ("Gpcr_Phenethylamine", "NCCc1ccccc1"),
            ("Gpcr_Dopamine", "NCCc1ccc(O)c(O)c1"),
            ("Gpcr_Histamine", "NCCc1c[nH]cn1"),
            ("Gpcr_Aryloxypropanolamine", "CC(C)NCC(O)COc1ccccc1"),
            ("Gpcr_Arylpiperazine", "c1ccc(N2CCNCC2)cc1"),
            ("Gpcr_Phenylpiperidine", "c1ccc(C2CCNCC2)cc1"),
            (
                "Gpcr_14Dihydropyridine",
                "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1",
            ),
            ("Gpcr_BenzodiazepineCore", "O=C1CN=C(c2ccccc2)c2ccccc2N1"),
            ("Gpcr_SulfonylPiperazine", "CS(=O)(=O)N1CCN(c2ccccc2)CC1"),
            ("Macro_Captopril", "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O"),
            (
                "Macro_Enalaprilat",
                "C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O",
            ),
        ],
    };

    pairs
        .iter()
        .filter_map(|&(name, smi)| {
            // Skip disconnected mixture fragments if any slipped in
            if smi.contains('.') {
                return None;
            }
            let mut m = parse_smiles(smi).ok()?;
            m.name = Some(name.to_string());
            Some(m)
        })
        .collect()
}

/// Write any built-in library preset to `~/.local/share/chemcraft/libraries/<stem>.smi`
/// (or compressed `<stem>.sdf.gz`) so users can customize/inspect catalogs on disk.
pub fn materialize_builtin_library(
    preset: BuiltinLibraryPreset,
    dir: &Path,
) -> Result<PathBuf, LibraryError> {
    use std::fmt::Write as _;
    std::fs::create_dir_all(dir).map_err(|e| LibraryError::Io(e.to_string()))?;
    let out = dir.join(format!("{}.smi", preset.filename_stem()));
    let mols = load_builtin_library(preset);
    let mut text = String::new();
    for m in &mols {
        let _ = writeln!(
            text,
            "{} {}",
            m.canonical_smiles(),
            m.name.as_deref().unwrap_or("Mol")
        );
    }
    std::fs::write(&out, text.as_bytes()).map_err(|e| LibraryError::Io(e.to_string()))?;
    Ok(out)
}

/// Stream-parse a chemical library file (`.sdf`, `.sdf.gz`, `.smi`, `.smi.gz`, `.csv`, `.dwar`)
/// up to `max_records` without reading an unbounded decompressed blob into memory.
pub fn stream_load_library(path: &Path, max_records: usize) -> Result<Vec<MolGraph>, LibraryError> {
    let fname = path
        .file_name()
        .map(|s| s.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();
    let file = std::fs::File::open(path)
        .map_err(|e| LibraryError::Io(format!("open {}: {e}", path.display())))?;
    stream_load_library_from_file(file, &fname, max_records)
}

/// Stream-parse an already-opened regular `File` (e.g. from `PathSandbox::open_read`) safely.
pub fn stream_load_library_from_file(
    file: std::fs::File,
    filename_hint: &str,
    max_records: usize,
) -> Result<Vec<MolGraph>, LibraryError> {
    let limit = max_records.max(1);
    let fname = filename_hint.to_ascii_lowercase();

    if fname.ends_with(".sdf.gz") || fname.ends_with(".sd.gz") {
        let dec = GzDecoder::new(file).take(MAX_LIBRARY_UNCOMPRESSED_BYTES);
        return stream_parse_sdf_reader(BufReader::new(dec), limit);
    }
    if fname.ends_with(".smi.gz") || fname.ends_with(".smiles.gz") {
        let dec = GzDecoder::new(file).take(MAX_LIBRARY_UNCOMPRESSED_BYTES);
        return stream_parse_smi_reader(BufReader::new(dec), limit);
    }
    if fname.ends_with(".sdf") || fname.ends_with(".sd") {
        let bounded = file.take(MAX_LIBRARY_UNCOMPRESSED_BYTES);
        return stream_parse_sdf_reader(BufReader::new(bounded), limit);
    }
    if fname.ends_with(".smi") || fname.ends_with(".smiles") {
        let bounded = file.take(MAX_LIBRARY_UNCOMPRESSED_BYTES);
        return stream_parse_smi_reader(BufReader::new(bounded), limit);
    }

    // CSV or DWAR bounded read
    let mut txt = String::new();
    file.take(MAX_LIBRARY_UNCOMPRESSED_BYTES)
        .read_to_string(&mut txt)
        .map_err(|e| LibraryError::Io(format!("read {filename_hint}: {e}")))?;

    if fname.ends_with(".dwar") {
        let recs = parse_dwar(&txt).map_err(|e| LibraryError::Parse(e.to_string()))?;
        return Ok(recs
            .into_iter()
            .filter_map(Result::ok)
            .take(limit)
            .map(|r| {
                let mut m = r.mol;
                if m.name.is_none() && !r.name.is_empty() {
                    m.name = Some(r.name);
                }
                m
            })
            .collect());
    }
    if fname.ends_with(".csv") {
        let recs = parse_csv_records(&txt).map_err(|e| LibraryError::Parse(e.to_string()))?;
        return Ok(recs
            .into_iter()
            .filter_map(Result::ok)
            .take(limit)
            .map(|r| {
                let mut m = r.mol;
                if m.name.is_none() && !r.name.is_empty() {
                    m.name = Some(r.name);
                }
                m
            })
            .collect());
    }

    Err(LibraryError::Parse(format!(
        "unsupported library extension on {filename_hint} (expected .sdf/.sdf.gz/.smi/.smi.gz/.csv/.dwar)"
    )))
}

fn stream_parse_sdf_reader<R: BufRead>(
    mut reader: R,
    max_records: usize,
) -> Result<Vec<MolGraph>, LibraryError> {
    let mut out = Vec::new();
    let mut block = String::new();
    let mut line = String::new();
    while out.len() < max_records {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| LibraryError::Io(e.to_string()))?;
        if n == 0 {
            if !block.trim().is_empty()
                && let Ok(m) = parse_mol_block(&block)
            {
                out.push(m);
            }
            break;
        }
        if line.trim() == "$$$$" {
            if !block.trim().is_empty()
                && let Ok(mut m) = parse_mol_block(&block)
            {
                if m.name.is_none() {
                    m.name = Some(format!("Lib_{}", out.len() + 1));
                }
                out.push(m);
            }
            block.clear();
        } else {
            block.push_str(&line);
        }
    }
    Ok(out)
}

fn stream_parse_smi_reader<R: BufRead>(
    mut reader: R,
    max_records: usize,
) -> Result<Vec<MolGraph>, LibraryError> {
    let mut out = Vec::new();
    let mut line = String::new();
    while out.len() < max_records {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| LibraryError::Io(e.to_string()))?;
        if n == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut parts = trimmed.splitn(2, char::is_whitespace);
        let smi = parts.next().unwrap_or("");
        if smi.eq_ignore_ascii_case("SMILES") {
            continue;
        }
        let name = parts
            .next()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string);
        if let Ok(mut m) = parse_smiles(smi) {
            m.name = Some(name.unwrap_or_else(|| format!("Lib_{}", out.len() + 1)));
            out.push(m);
        }
    }
    Ok(out)
}

/// Apply approximate physiological pH 7.4 ionization rules to a `MolGraph`:
/// - Carboxylic / sulfonic acids (`-C(=O)OH`, `-S(=O)2OH`) -> `-O-` (`formal_charge = -1`)
/// - Tetrazole NH -> `[n-]`
/// - Aliphatic amines (non-aromatic `N` not conjugated to `C=O`, `S=O`, `P=O`, aryl, or `C=N`) -> `[N+]`
pub fn protonate_at_ph74(mol: &mut MolGraph) {
    let n = mol.atoms.len();
    let mut deprot_o = Vec::new();
    let mut deprot_n = Vec::new();
    let mut prot_n = Vec::new();

    for i in 0..n {
        let a = &mol.atoms[i];
        if a.formal_charge != 0 || a.valence_exceeded {
            continue;
        }
        // Carboxylic/sulfonic acid -OH
        if a.symbol == "O" && mol.adj[i].len() == 1 && a.total_h() == 1 {
            let parent = mol.adj[i][0].0;
            let has_carbonyl_or_sulfonyl =
                matches!(mol.atoms[parent].symbol.as_str(), "C" | "S" | "P")
                    && mol.adj[parent].iter().any(|&(nb, bid)| {
                        nb != i
                            && mol.atoms[nb].symbol == "O"
                            && mol.bonds[bid].order == BondOrder::Double
                    });
            if has_carbonyl_or_sulfonyl {
                deprot_o.push(i);
            }
        }
        // Acidic tetrazole NH in 5-membered ring with >=3 N atoms
        if a.symbol == "N" && a.aromatic && a.total_h() == 1 {
            let in_tetrazole = mol.rings.iter().any(|r| {
                r.len() == 5
                    && r.contains(&i)
                    && r.iter().filter(|&&x| mol.atoms[x].symbol == "N").count() >= 4
            });
            if in_tetrazole {
                deprot_n.push(i);
            }
        }
        // Basic aliphatic amine
        if a.symbol == "N"
            && !a.aromatic
            && mol.adj[i]
                .iter()
                .all(|&(_, bid)| mol.bonds[bid].order == BondOrder::Single)
        {
            let conjugated_or_ewg = mol.adj[i].iter().any(|&(nb, _)| {
                mol.atoms[nb].aromatic
                    || mol.adj[nb].iter().any(|&(nb2, bid2)| {
                        nb2 != i
                            && matches!(mol.atoms[nb2].symbol.as_str(), "O" | "S" | "N")
                            && matches!(
                                mol.bonds[bid2].order,
                                BondOrder::Double | BondOrder::Triple | BondOrder::Aromatic
                            )
                    })
            });
            if !conjugated_or_ewg {
                prot_n.push(i);
            }
        }
    }

    let changed = !deprot_o.is_empty() || !deprot_n.is_empty() || !prot_n.is_empty();
    for i in deprot_o {
        mol.atoms[i].formal_charge = -1;
        mol.atoms[i].explicit_h = 0;
    }
    for i in deprot_n {
        mol.atoms[i].formal_charge = -1;
        mol.atoms[i].explicit_h = 0;
    }
    for i in prot_n {
        mol.atoms[i].formal_charge = 1;
    }
    if changed {
        let _ = mol.perceive();
    }
}

fn passes_funnel(eval: &MoleculeEvaluation, funnel: ScreeningFunnel) -> bool {
    let mw = eval.physicochemical.molecular_weight;
    let logp = eval.lipophilicity.consensus_log_p;
    let tpsa = eval.physicochemical.tpsa;
    let hbd = eval.physicochemical.h_bond_donors;
    let hba = eval.physicochemical.h_bond_acceptors;
    let pains = eval.medicinal_chemistry.pains_alerts.len();
    let brenk = eval.medicinal_chemistry.brenk_alerts.len();

    match funnel {
        ScreeningFunnel::All => true,
        ScreeningFunnel::RuleOf5 => eval.druglikeness.lipinski.violations <= 1 && pains == 0,
        ScreeningFunnel::LeadLike => {
            (200.0..=350.0).contains(&mw)
                && (-1.0..=3.5).contains(&logp)
                && eval.physicochemical.rotatable_bonds <= 7
                && pains == 0
                && brenk == 0
        }
        ScreeningFunnel::FragmentRo3 => {
            mw <= 300.0 && logp <= 3.0 && hbd <= 3 && hba <= 3 && tpsa <= 60.0 && pains == 0
        }
        ScreeningFunnel::BbbCns => {
            mw <= 400.0 && tpsa <= 90.0 && eval.pharmacokinetics.bbb_permeant && pains == 0
        }
    }
}

/// Greedy MaxMin diverse subset selection in ECFP4 Tanimoto space.
pub fn maxmin_diversity_pick(
    hits: Vec<ScreenedLibraryHit>,
    pick_n: usize,
) -> Vec<ScreenedLibraryHit> {
    if pick_n == 0 || hits.len() <= pick_n {
        return hits;
    }
    let total = hits.len();
    let mut selected_idx = Vec::with_capacity(pick_n);
    let mut min_sim = vec![0.0_f64; total];
    let mut chosen = vec![false; total];

    // Seed with first hit (highest similarity to ref if present, else catalog order)
    selected_idx.push(0);
    chosen[0] = true;
    for i in 0..total {
        if i != 0 {
            min_sim[i] = hits[i].fp.tanimoto(&hits[0].fp);
        }
    }

    while selected_idx.len() < pick_n {
        // Pick compound with SMALLEST max-similarity (= largest MaxMin distance 1 - max_sim)
        let mut best_i = None;
        let mut best_sim = f64::INFINITY;
        for i in 0..total {
            if !chosen[i] && min_sim[i] < best_sim {
                best_sim = min_sim[i];
                best_i = Some(i);
            }
        }
        let Some(pick) = best_i else { break };
        chosen[pick] = true;
        selected_idx.push(pick);
        for i in 0..total {
            if !chosen[i] {
                let s = hits[i].fp.tanimoto(&hits[pick].fp);
                if s > min_sim[i] {
                    min_sim[i] = s;
                }
            }
        }
    }

    // Move out cleanly preserving selection order
    let mut opt: Vec<Option<ScreenedLibraryHit>> = hits.into_iter().map(Some).collect();
    selected_idx
        .into_iter()
        .filter_map(|i| opt[i].take())
        .collect()
}

/// Run the parallel *in silico* library filtering + pH 7.4 ionization + 3D preparation funnel.
pub fn filter_and_prepare_library(
    source_label: &str,
    mols: Vec<MolGraph>,
    cfg: &LibraryFilterConfig,
) -> Result<LibraryScreeningResult, LibraryError> {
    let t0 = std::time::Instant::now();
    let total_loaded = mols.len();
    if total_loaded == 0 {
        return Err(LibraryError::EmptyResult);
    }

    let smarts_q = match cfg
        .smarts_substructure
        .as_deref()
        .map(str::trim)
        .filter(|s| !s.is_empty())
    {
        Some(q) => {
            Some(SmartsQuery::parse(q).map_err(|e| LibraryError::InvalidQuery(e.to_string()))?)
        }
        None => None,
    };

    let ref_fp = match cfg
        .reference_smiles
        .as_deref()
        .map(str::trim)
        .filter(|s| !s.is_empty())
    {
        Some(smi) => {
            let m = parse_smiles(smi).map_err(|e| LibraryError::InvalidQuery(e.to_string()))?;
            Some(MorganFp::ecfp4(&m))
        }
        None => None,
    };

    let provider = AdmetCraftProvider::default();

    // Parallel ADMET + SMARTS + Tanimoto evaluation
    let mut hits: Vec<ScreenedLibraryHit> = mols
        .into_par_iter()
        .enumerate()
        .filter_map(|(idx, mut mol)| {
            if mol.atoms.is_empty() {
                return None;
            }
            if let Some(ref sq) = smarts_q
                && !sq.matches(&mol)
            {
                return None;
            }
            let fp = MorganFp::ecfp4(&mol);
            let tanimoto_to_ref = ref_fp.as_ref().map(|rf| fp.tanimoto(rf));
            if let Some(t) = tanimoto_to_ref
                && t < cfg.min_tanimoto
            {
                return None;
            }
            let can = mol.canonical_smiles();
            if can.is_empty() {
                return None;
            }
            let name = mol
                .name
                .clone()
                .filter(|s| !s.is_empty())
                .unwrap_or_else(|| format!("Compound_{}", idx + 1));
            mol.name = Some(name.clone());

            let eval = provider.evaluate_smiles(&can, Some(&name)).ok()?;
            if !passes_funnel(&eval, cfg.funnel) {
                return None;
            }
            if cfg.protonate_ph74 {
                protonate_at_ph74(&mut mol);
            }
            if cfg.generate_3d {
                ensure_3d_coords(&mut mol);
            }
            Some(ScreenedLibraryHit {
                name,
                canonical_smiles: can,
                mol,
                eval,
                tanimoto_to_ref,
                fp,
            })
        })
        .collect();

    let passed_funnel = hits.len();
    if hits.is_empty() {
        return Err(LibraryError::EmptyResult);
    }

    if ref_fp.is_some() {
        hits.sort_by(|a, b| {
            b.tanimoto_to_ref
                .unwrap_or(0.0)
                .total_cmp(&a.tanimoto_to_ref.unwrap_or(0.0))
        });
    }

    if let Some(k) = cfg.diversity_pick_n
        && k > 0
        && hits.len() > k
    {
        hits = maxmin_diversity_pick(hits, k);
    }

    let cap = cfg.max_hits.max(1);
    if hits.len() > cap {
        hits.truncate(cap);
    }

    Ok(LibraryScreeningResult {
        source_label: source_label.to_string(),
        total_loaded,
        passed_funnel,
        hits,
        elapsed_ms: t0.elapsed().as_millis(),
    })
}

/// Write a `LibraryScreeningResult`'s hits as a multi-record 3D SDF V2000 library file.
pub fn export_library_hits_sdf3d(
    res: &LibraryScreeningResult,
    out: &Path,
) -> Result<usize, LibraryError> {
    if res.hits.is_empty() {
        return Err(LibraryError::EmptyResult);
    }
    if let Some(parent) = out.parent()
        && !parent.as_os_str().is_empty()
    {
        std::fs::create_dir_all(parent).map_err(|e| LibraryError::Io(e.to_string()))?;
    }
    let refs: Vec<&MolGraph> = res.hits.iter().map(|h| &h.mol).collect();
    let sdf = write_sdf_v2000_3d(&refs);
    std::fs::write(out, sdf.as_bytes()).map_err(|e| LibraryError::Io(e.to_string()))?;
    Ok(refs.len())
}
