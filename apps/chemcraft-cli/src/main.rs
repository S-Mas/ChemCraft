#![forbid(unsafe_code)]

mod mcp;

use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand, ValueEnum};
use mol_admet_bridge::{AdmetCraftProvider, AdmetProvider, AdmetReportFormat};
use mol_engine::{
    BuiltinLibraryPreset, DEFAULT_MAX_LIBRARY_RECORDS, LibraryFilterConfig, PathSandbox,
    ScreeningFunnel, default_chimerax_export_dir, default_dockcraft_handoff_dir,
    export_library_hits_sdf3d, filter_and_prepare_library, launch_chimerax_for_mol,
    launch_dockcraft_for_mols, load_builtin_library, stream_load_library, write_chimerax_bundle,
    write_dockcraft_bundle,
};
use mol_graph::{
    MolGraph, SecondaryStructure, build_peptide, ensure_3d_coords, generate_3d_coords,
    parse_smiles, write_pdb, write_sdf_v2000, write_sdf_v2000_3d,
};
use mol_sar::{
    SarFilter, SarTable, find_activity_cliffs, find_duplicate_compounds,
    find_matched_molecular_pairs,
};
use mol_vector::{
    SceneStyle, build_mol_scene, generate_2d_coords, scene_to_eps, scene_to_pdf_bytes,
    scene_to_png_bytes, scene_to_svg,
};

#[derive(Debug, Parser)]
#[command(
    name = "chemcraft-cli",
    version,
    about = "ChemCraft CLI & MCP server — Clean-room 2D vector chemistry, DataWarrior-style SAR, and live AdmetCraft integration"
)]
struct Cli {
    #[command(subcommand)]
    cmd: SubCmd,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum RenderFormat {
    Svg,
    Eps,
    Pdf,
    Png,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum EnrichFormat {
    Csv,
    Json,
    Pdf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum PeptideSsCli {
    Beta,
    Alpha,
    Extended,
}

impl From<PeptideSsCli> for SecondaryStructure {
    fn from(value: PeptideSsCli) -> Self {
        match value {
            PeptideSsCli::Beta => SecondaryStructure::BetaStrand,
            PeptideSsCli::Alpha => SecondaryStructure::AlphaHelix,
            PeptideSsCli::Extended => SecondaryStructure::Extended,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum Export3dFormat {
    Pdb,
    Sdf3d,
}

#[derive(Debug, Subcommand)]
enum SubCmd {
    /// Render a SMILES string to vector SVG, EPS, PDF, or raster PNG.
    Render {
        smiles: String,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, value_enum, default_value_t = RenderFormat::Svg)]
        format: RenderFormat,
        #[arg(long, default_value_t = 520.0)]
        width: f64,
        #[arg(long, default_value_t = 380.0)]
        height: f64,
        #[arg(long)]
        dark: bool,
    },
    /// Filter a chemical library (.sdf, .dwar, .csv) by SMARTS substructure, Tanimoto, and Lipinski/PAINS.
    Filter {
        input: PathBuf,
        #[arg(long)]
        substructure: Option<String>,
        #[arg(long)]
        reference: Option<String>,
        #[arg(long, default_value_t = 0.0)]
        min_tanimoto: f64,
        #[arg(long)]
        lipinski_only: bool,
        #[arg(long)]
        no_pains: bool,
        #[arg(long)]
        out: PathBuf,
    },
    /// Enrich an SDF/DWAR/CSV library with live AdmetCraft predictions and export CSV/JSON/PDF.
    Enrich {
        input: PathBuf,
        #[arg(long)]
        with_admet: bool,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, value_enum, default_value_t = EnrichFormat::Csv)]
        format: EnrichFormat,
    },
    /// Compute Matched Molecular Pairs (MMP) and Activity Cliffs (SALI).
    Mmp {
        input: PathBuf,
        #[arg(long, default_value_t = 0.40)]
        min_tanimoto: f64,
        #[arg(long)]
        out: Option<PathBuf>,
    },
    /// Build a peptide in 2D + 3D from 1-letter or 3-letter sequence (`Ac-Ala-Gly-Ser-NH2`, `AGYS`) and export SVG/PDB/3D-SDF.
    Peptide {
        sequence: String,
        #[arg(long, value_enum, default_value_t = PeptideSsCli::Beta)]
        ss: PeptideSsCli,
        #[arg(long)]
        svg: Option<PathBuf>,
        #[arg(long)]
        pdb: Option<PathBuf>,
        #[arg(long)]
        sdf3d: Option<PathBuf>,
    },
    /// Generate a 3D conformer from SMILES and export PDB or 3D SDF V2000.
    Export3d {
        smiles: String,
        #[arg(long)]
        out: PathBuf,
        #[arg(long, value_enum, default_value_t = Export3dFormat::Pdb)]
        format: Export3dFormat,
    },
    /// Export a self-contained ChimeraX 3D bundle (.pdb/.sdf + .cxc) and optionally launch ChimeraX.
    Chimerax {
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        peptide: Option<String>,
        #[arg(long, value_enum, default_value_t = PeptideSsCli::Alpha)]
        ss: PeptideSsCli,
        #[arg(long)]
        out_dir: Option<PathBuf>,
        #[arg(long)]
        launch: bool,
    },
    /// Export a 3D ligand or chemical library (.sdf/.smi/.pdb) for DockCraft virtual screening and optionally launch DockCraft.
    Dockcraft {
        #[arg(long)]
        smiles: Option<String>,
        #[arg(long)]
        peptide: Option<String>,
        #[arg(long)]
        library: Option<PathBuf>,
        #[arg(long, value_enum, default_value_t = PeptideSsCli::Beta)]
        ss: PeptideSsCli,
        #[arg(long)]
        out_dir: Option<PathBuf>,
        #[arg(long)]
        launch: bool,
    },
    /// Stream, pre-filter (Ro5/Lead-Like/Ro3/BBB/SMARTS/MaxMin), protonate (pH 7.4), and prepare a 3D screening library for DockCraft.
    Library {
        /// Built-in preset ('fda', 'fragments', 'kinase')
        #[arg(long)]
        preset: Option<String>,
        /// Local library file (.sdf, .sdf.gz, .smi, .smi.gz, .csv, .dwar)
        #[arg(long)]
        input: Option<PathBuf>,
        /// Funnel ('all', 'ro5', 'lead', 'ro3', 'bbb')
        #[arg(long, default_value = "ro5")]
        funnel: String,
        #[arg(long)]
        smarts: Option<String>,
        #[arg(long)]
        reference: Option<String>,
        #[arg(long, default_value_t = 0.0)]
        min_tanimoto: f64,
        /// Optional MaxMin ECFP4 diversity subset size
        #[arg(long)]
        diverse_n: Option<usize>,
        #[arg(long, default_value_t = 500)]
        max_hits: usize,
        #[arg(long, default_value_t = true)]
        ph74: bool,
        #[arg(long)]
        out: PathBuf,
        #[arg(long)]
        send_dockcraft: bool,
    },
    /// Inspect SMILES structure, key, coordinates, and AdmetCraft profile as JSON.
    Inspect { smiles: String },
    /// Start the hardened MCP server (stdio or TCP with token + sandbox).
    Mcp {
        #[arg(long)]
        read_root: Option<PathBuf>,
        #[arg(long)]
        write_root: Option<PathBuf>,
        #[arg(long)]
        tcp: Option<String>,
        #[arg(long)]
        allow_remote: bool,
        #[arg(long)]
        mcp_token: Option<String>,
    },
}

fn load_table_from_path(path: &Path) -> Result<SarTable, Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(path)?;
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let mut table = SarTable::new();
    if ext == "dwar" {
        table.load_dwar_text(&text)?;
    } else if ext == "csv" {
        table.load_csv_text(&text)?;
    } else {
        table.load_sdf_text(&text)?;
    }
    Ok(table)
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();
    match cli.cmd {
        SubCmd::Render {
            smiles,
            out,
            format,
            width,
            height,
            dark,
        } => {
            let mut mol = parse_smiles(&smiles)?;
            generate_2d_coords(&mut mol);
            let style = SceneStyle {
                width,
                height,
                dark_mode: dark,
                ..SceneStyle::default()
            };
            let scene = build_mol_scene(&mol, &style, &[]);
            let bytes = match format {
                RenderFormat::Svg => scene_to_svg(&scene).into_bytes(),
                RenderFormat::Eps => scene_to_eps(&scene).into_bytes(),
                RenderFormat::Pdf => scene_to_pdf_bytes(&scene)?,
                RenderFormat::Png => scene_to_png_bytes(&scene, 2.0)?,
            };
            std::fs::write(&out, bytes)?;
            eprintln!("Rendered {} -> {}", mol.canonical_smiles(), out.display());
        }
        SubCmd::Filter {
            input,
            substructure,
            reference,
            min_tanimoto,
            lipinski_only,
            no_pains,
            out,
        } => {
            let mut table = load_table_from_path(&input)?;
            if lipinski_only || no_pains {
                table.enrich_with_admet(&AdmetCraftProvider::default());
            }
            let filter = SarFilter {
                smarts_substructure: substructure,
                reference_smiles: reference,
                min_tanimoto,
                lipinski_pass_only: lipinski_only,
                no_pains_only: no_pains,
                ..SarFilter::default()
            };
            let hit_indices = table.matching_indices(&filter)?;
            let refs: Vec<&mol_graph::MolGraph> = hit_indices
                .iter()
                .map(|&idx| &table.rows[idx].mol)
                .collect();
            let sdf_out = write_sdf_v2000(&refs);
            std::fs::write(&out, sdf_out)?;
            eprintln!(
                "Filtered {}/{} compounds -> {}",
                hit_indices.len(),
                table.rows.len(),
                out.display()
            );
        }
        SubCmd::Enrich {
            input,
            with_admet,
            out,
            format,
        } => {
            let mut table = load_table_from_path(&input)?;
            run_enrich_command(&mut table, with_admet, &out, format)?;
        }
        SubCmd::Mmp {
            input,
            min_tanimoto,
            out,
        } => {
            let mut table = load_table_from_path(&input)?;
            table.enrich_with_admet(&AdmetCraftProvider::default());
            let duplicates = find_duplicate_compounds(&table.rows);
            let cliffs = find_activity_cliffs(&table.rows, min_tanimoto);
            let mmps = find_matched_molecular_pairs(&table.rows);
            let payload = serde_json::to_string_pretty(&serde_json::json!({
                "compounds": table.rows.len(),
                "duplicates": duplicates,
                "activity_cliffs": cliffs,
                "matched_molecular_pairs": mmps,
            }))?;
            if let Some(p) = out {
                std::fs::write(p, payload)?;
            } else {
                println!("{payload}");
            }
        }
        SubCmd::Inspect { smiles } => {
            let mut mol = parse_smiles(&smiles)?;
            generate_2d_coords(&mut mol);
            let eval = AdmetCraftProvider::default().evaluate_smiles(&smiles, None)?;
            println!(
                "{}",
                serde_json::to_string_pretty(&serde_json::json!({
                    "canonical_smiles": mol.canonical_smiles(),
                    "formula": mol.formula(),
                    "structural_key": mol.structural_key(),
                    "warnings": mol.warnings,
                    "admet": eval,
                }))?
            );
        }
        SubCmd::Mcp {
            read_root,
            write_root,
            tcp,
            allow_remote,
            mcp_token,
        } => {
            let sandbox = PathSandbox::new(read_root.as_deref(), write_root.as_deref())?;
            if let Some(addr) = tcp {
                let tok = mcp_token.or_else(|| std::env::var("CHEMCRAFT_MCP_TOKEN").ok());
                mcp::run_mcp_tcp(&addr, &sandbox, allow_remote, tok)?;
            } else {
                mcp::run_mcp_stdio(&sandbox)?;
            }
        }
        SubCmd::Peptide {
            sequence,
            ss,
            svg,
            pdb,
            sdf3d,
        } => {
            let mol = build_peptide(&sequence, ss.into())?;
            if let Some(svg_out) = svg {
                let scene = build_mol_scene(&mol, &SceneStyle::default(), &[]);
                std::fs::write(&svg_out, scene_to_svg(&scene).as_bytes())?;
            }
            if let Some(pdb_out) = pdb {
                std::fs::write(&pdb_out, write_pdb(&mol).as_bytes())?;
            }
            if let Some(sdf_out) = sdf3d {
                std::fs::write(&sdf_out, write_sdf_v2000_3d(&[&mol]).as_bytes())?;
            }
            let summary = serde_json::json!({
                "name": mol.name,
                "sequence": mol.properties.get("PEPTIDE_SEQUENCE"),
                "one_letter": mol.properties.get("PEPTIDE_1LETTER"),
                "secondary_structure": mol.properties.get("SECONDARY_STRUCTURE"),
                "residues": mol.properties.get("PEPTIDE_RESIDUES"),
                "est_charge_ph74": mol.properties.get("EST_CHARGE_PH74"),
                "formula": mol.formula(),
                "mw": mol.molecular_weight(),
                "canonical_smiles": mol.canonical_smiles(),
            });
            println!("{}", serde_json::to_string_pretty(&summary)?);
        }
        SubCmd::Export3d {
            smiles,
            out,
            format,
        } => {
            let mut mol = parse_smiles(&smiles)?;
            generate_3d_coords(&mut mol);
            let txt = match format {
                Export3dFormat::Pdb => write_pdb(&mol),
                Export3dFormat::Sdf3d => write_sdf_v2000_3d(&[&mol]),
            };
            std::fs::write(&out, txt.as_bytes())?;
            eprintln!("Exported 3D -> {}", out.display());
        }
        SubCmd::Chimerax {
            smiles,
            peptide,
            ss,
            out_dir,
            launch,
        } => {
            let mol = if let Some(seq) = peptide {
                build_peptide(&seq, ss.into())?
            } else if let Some(smi) = smiles {
                let mut m = parse_smiles(&smi)?;
                generate_3d_coords(&mut m);
                m
            } else {
                return Err("provide --smiles or --peptide".into());
            };
            let target_dir = out_dir.unwrap_or_else(default_chimerax_export_dir);
            let bundle = if launch {
                launch_chimerax_for_mol(&mol, Some(&target_dir))?
            } else {
                write_chimerax_bundle(&mol, &target_dir, None)?
            };
            let info = serde_json::json!({
                "dir": bundle.dir.display().to_string(),
                "script": bundle.script_path.display().to_string(),
                "structure": bundle.structure_path.display().to_string(),
                "is_peptide": bundle.is_peptide,
                "launched": launch,
            });
            println!("{}", serde_json::to_string_pretty(&info)?);
        }
        SubCmd::Dockcraft {
            smiles,
            peptide,
            library,
            ss,
            out_dir,
            launch,
        } => {
            let mut mols: Vec<MolGraph> = if let Some(lib_path) = library {
                let table = load_table_from_path(&lib_path)?;
                table
                    .rows
                    .into_iter()
                    .map(|r| {
                        let mut m = r.mol;
                        ensure_3d_coords(&mut m);
                        m
                    })
                    .collect()
            } else if let Some(seq) = peptide {
                vec![build_peptide(&seq, ss.into())?]
            } else if let Some(smi) = smiles {
                let mut m = parse_smiles(&smi)?;
                generate_3d_coords(&mut m);
                vec![m]
            } else {
                return Err("provide --smiles, --peptide, or --library".into());
            };
            for m in &mut mols {
                ensure_3d_coords(m);
            }
            let refs: Vec<&MolGraph> = mols.iter().collect();
            let target_dir = out_dir.unwrap_or_else(default_dockcraft_handoff_dir);
            let bundle = if launch {
                launch_dockcraft_for_mols(&refs, Some(&target_dir), None)?
            } else {
                write_dockcraft_bundle(&refs, &target_dir, None)?
            };
            let info = serde_json::json!({
                "dir": bundle.dir.display().to_string(),
                "sdf_3d": bundle.sdf_path.display().to_string(),
                "smi": bundle.smi_path.display().to_string(),
                "pdb": bundle.pdb_path.as_ref().map(|p| p.display().to_string()),
                "molecule_count": bundle.molecule_count,
                "launched": launch,
            });
            println!("{}", serde_json::to_string_pretty(&info)?);
        }
        SubCmd::Library {
            preset,
            input,
            funnel,
            smarts,
            reference,
            min_tanimoto,
            diverse_n,
            max_hits,
            ph74,
            out,
            send_dockcraft,
        } => {
            let (source_label, mols) = if let Some(p) = input {
                (
                    p.display().to_string(),
                    stream_load_library(&p, DEFAULT_MAX_LIBRARY_RECORDS)?,
                )
            } else if let Some(ref pre) = preset {
                let p = BuiltinLibraryPreset::from_str_lossy(pre).ok_or_else(|| {
                    format!("unknown library preset '{pre}' (expected fda|fragments|kinase)")
                })?;
                (p.label().to_string(), load_builtin_library(p))
            } else {
                return Err("provide --preset or --input".into());
            };
            let funnel_enum = ScreeningFunnel::from_str_lossy(&funnel).ok_or_else(|| {
                format!("unknown funnel '{funnel}' (expected all|ro5|lead|ro3|bbb)")
            })?;
            let cfg = LibraryFilterConfig {
                funnel: funnel_enum,
                smarts_substructure: smarts,
                reference_smiles: reference,
                min_tanimoto,
                diversity_pick_n: diverse_n,
                max_hits,
                protonate_ph74: ph74,
                generate_3d: true,
            };
            let res = filter_and_prepare_library(&source_label, mols, &cfg)?;
            let count = export_library_hits_sdf3d(&res, &out)?;
            if send_dockcraft {
                let refs: Vec<&MolGraph> = res.hits.iter().map(|h| &h.mol).collect();
                let _ = write_dockcraft_bundle(
                    &refs,
                    &default_dockcraft_handoff_dir(),
                    Some("chemcraft_screen_lib"),
                )?;
            }
            let summary = serde_json::json!({
                "source": res.source_label,
                "total_loaded": res.total_loaded,
                "passed_funnel": res.passed_funnel,
                "selected_count": count,
                "elapsed_ms": res.elapsed_ms,
                "out_sdf_3d": out.display().to_string(),
            });
            println!("{}", serde_json::to_string_pretty(&summary)?);
        }
    }
    Ok(())
}

fn escape_csv_cell(raw: &str) -> String {
    let trimmed = raw.trim_start();
    let safe = if trimmed.starts_with(['=', '+', '-', '@']) {
        format!("'{raw}")
    } else {
        raw.to_string()
    };
    if safe.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", safe.replace('"', "\"\""))
    } else {
        safe
    }
}

fn run_enrich_command(
    table: &mut SarTable,
    with_admet: bool,
    out: &Path,
    format: EnrichFormat,
) -> Result<(), Box<dyn std::error::Error>> {
    if with_admet {
        let provider = AdmetCraftProvider::default();
        table.enrich_with_admet(&provider);
        let evals: Vec<_> = table
            .rows
            .iter()
            .filter_map(|r| r.admet_full.clone())
            .collect();
        let fmt = match format {
            EnrichFormat::Csv => AdmetReportFormat::Csv,
            EnrichFormat::Json => AdmetReportFormat::Json,
            EnrichFormat::Pdf => AdmetReportFormat::Pdf,
        };
        let mut f = std::fs::File::create(out)?;
        provider.export_report_to_writer(&evals, &mut f, fmt)?;
        eprintln!(
            "Wrote ADMET report ({} compounds) -> {}",
            evals.len(),
            out.display()
        );
        return Ok(());
    }

    match format {
        EnrichFormat::Csv => {
            use std::fmt::Write as _;
            let mut csv_out = String::from(
                "Molecule,Canonical SMILES,Structural Key,Formula,MW,HeavyAtoms,pIC50,Series\n",
            );
            for r in &table.rows {
                let pic50 = r.pic50.map_or_else(String::new, |v| format!("{v:.3}"));
                let _ = writeln!(
                    csv_out,
                    "{},{},{},{},{:.2},{},{},{}",
                    escape_csv_cell(&r.name),
                    escape_csv_cell(&r.canonical_smiles),
                    escape_csv_cell(&r.structural_key),
                    escape_csv_cell(&r.formula),
                    r.mw,
                    r.heavy_atoms,
                    pic50,
                    escape_csv_cell(&r.series)
                );
            }
            std::fs::write(out, csv_out)?;
        }
        EnrichFormat::Json => {
            let rows: Vec<_> = table
                .rows
                .iter()
                .map(|r| {
                    serde_json::json!({
                        "name": r.name,
                        "canonical_smiles": r.canonical_smiles,
                        "structural_key": r.structural_key,
                        "formula": r.formula,
                        "mw": r.mw,
                        "heavy_atoms": r.heavy_atoms,
                        "pic50": r.pic50,
                        "series": r.series,
                    })
                })
                .collect();
            std::fs::write(out, serde_json::to_string_pretty(&rows)?)?;
        }
        EnrichFormat::Pdf => {
            return Err("PDF ADMET report requires --with-admet".into());
        }
    }
    eprintln!(
        "Wrote table metadata ({} compounds) -> {}",
        table.rows.len(),
        out.display()
    );
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("chemcraft-cli error: {e}");
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn e2e_50_compounds_sdf_substructure_filter_tetrazole_bioisostere_svg_and_admet()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_e2e_integration_test");
        let _ = std::fs::remove_dir_all(&base);
        let rdir = base.join("read");
        let wdir = base.join("write");
        std::fs::create_dir_all(&rdir)?;
        std::fs::create_dir_all(&wdir)?;

        // 1. Export realistic 50-compound library as SDF into read sandbox
        let demo = SarTable::demo_50_compounds();
        assert_eq!(demo.rows.len(), 50);
        let refs: Vec<&mol_graph::MolGraph> = demo.rows.iter().map(|r| &r.mol).collect();
        std::fs::write(rdir.join("library_50.sdf"), write_sdf_v2000(&refs))?;

        // 2. Drive full workflow via hardened MCP server: initialize -> search_substructure -> mutate_scaffold (tetrazole) -> export_vector_figure -> mol.analyze_admet
        let sb = PathSandbox::new(Some(&rdir), Some(&wdir))?;
        let requests = [
            serde_json::json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{}}).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"search_substructure",
                    "arguments":{"file":"library_50.sdf","smarts":"c1ccccc1","output_sdf":"aromatic_hits.sdf"}
                }
            }).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":3,"method":"tools/call",
                "params":{
                    "name":"mutate_scaffold",
                    "arguments":{"smiles":"CC(C)Cc1ccc(C(C)C(=O)O)cc1","group":"tetrazole","name":"Ibuprofen-Tetrazole"}
                }
            }).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":4,"method":"tools/call",
                "params":{
                    "name":"export_vector_figure",
                    "arguments":{"smiles":["CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1 Ibuprofen-Tetrazole"],"output_path":"bioisostere.svg","format":"svg"}
                }
            }).to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":5,"method":"tools/call",
                "params":{
                    "name":"mol.analyze_admet",
                    "arguments":{"smiles":["CC(C)Cc1ccc(C(C)C(=O)O)cc1 Ibuprofen","CC(C)Cc1ccc(C(C)c2nnn[nH]2)cc1 Ibuprofen-Tetrazole"],"output_report":"admet_delta.csv","format":"csv"}
                }
            }).to_string(),
            // Strict JSON-RPC 2.0 check: missing jsonrpc or invalid params -> -32602 / -32600
            serde_json::json!({"jsonrpc":"2.0","id":6,"method":"tools/call","params":{"name":"draw_from_smiles","arguments":{}}}).to_string(),
        ].join("\n") + "\n";

        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(requests.as_bytes()), &mut out_buf, &sb, None)?;
        let out_str = String::from_utf8(out_buf)?;
        let lines: Vec<serde_json::Value> = out_str
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 6);

        // Aromatic hits from 50-compound SDF
        let total_loaded = lines[1]["result"]["structuredContent"]["total_loaded"]
            .as_u64()
            .ok_or("missing total_loaded")?;
        let matched = lines[1]["result"]["structuredContent"]["matched_count"]
            .as_u64()
            .ok_or("missing matched_count")?;
        assert!(
            matched >= 45,
            "expected >=45 aromatic hits, got {matched}/{total_loaded}; full resp: {}",
            lines[1]
        );
        assert!(wdir.join("aromatic_hits.sdf").exists());

        // Tetrazole bioisostere mutation changes logP, TPSA (+17.16 A^2), and preserves GI absorption
        let delta_tpsa = lines[2]["result"]["structuredContent"]["delta_tpsa"]
            .as_f64()
            .ok_or("missing delta_tpsa")?;
        assert!(
            delta_tpsa > 15.0,
            "expected tetrazole to raise TPSA, got {delta_tpsa}"
        );
        let new_gi = lines[2]["result"]["structuredContent"]["new_gi"]
            .as_str()
            .ok_or("missing new_gi")?;
        assert_eq!(new_gi, "High");

        // SVG valid + zero overlaps + AdmetCraft CSV generated
        let svg_text = std::fs::read_to_string(wdir.join("bioisostere.svg"))?;
        assert!(svg_text.contains("<svg") && svg_text.contains("</svg>"));
        assert_eq!(lines[3]["result"]["structuredContent"]["overlap_count"], 0);
        assert!(wdir.join("admet_delta.csv").exists());
        assert_eq!(lines[5]["error"]["code"], -32602);

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn enrich_with_admet_differs_from_without_admet() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_enrich_flag_test");
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(&base)?;

        let mut table_plain = SarTable::new();
        table_plain.load_csv_text(
            "name,smiles,pIC50\nAspirin,CC(=O)Oc1ccccc1C(=O)O,5.2\nCaffeine,Cn1cnc2c1c(=O)n(C)c(=O)n2C,4.8\n",
        )?;
        let mut table_admet = table_plain.clone();

        let plain_csv = base.join("plain.csv");
        let admet_csv = base.join("admet.csv");
        run_enrich_command(&mut table_plain, false, &plain_csv, EnrichFormat::Csv)?;
        run_enrich_command(&mut table_admet, true, &admet_csv, EnrichFormat::Csv)?;

        let plain_text = std::fs::read_to_string(&plain_csv)?;
        let admet_text = std::fs::read_to_string(&admet_csv)?;
        assert_ne!(plain_text, admet_text);
        assert!(admet_text.contains("Consensus LogP"));
        assert!(!plain_text.contains("Consensus LogP"));

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_export_sar_grid_eps_and_png_include_multiple_compounds()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_grid_eps_png_test");
        let _ = std::fs::remove_dir_all(&base);
        let wdir = base.join("write");
        std::fs::create_dir_all(&wdir)?;

        let sb = PathSandbox::new(None, Some(&wdir))?;
        let smiles = [
            "CC(=O)Oc1ccccc1C(=O)O Aspirin",
            "Cn1cnc2c1c(=O)n(C)c(=O)n2C Caffeine",
            "CC(C)Cc1ccc(C(C)C(=O)O)cc1 Ibuprofen",
            "CC(=O)Nc1ccc(O)cc1 Paracetamol",
            "c1ccccc1 Benzene",
            "CCO Ethanol",
        ];
        let requests = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{
                    "name":"export_vector_figure",
                    "arguments":{"smiles":smiles,"output_path":"grid.eps","format":"eps"}
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"export_vector_figure",
                    "arguments":{"smiles":smiles,"output_path":"grid.png","format":"png"}
                }
            })
            .to_string(),
        ]
        .join("\n")
            + "\n";

        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(requests.as_bytes()), &mut out_buf, &sb, None)?;
        let eps = std::fs::read_to_string(wdir.join("grid.eps"))?;
        for label in [
            "Aspirin",
            "Caffeine",
            "Ibuprofen",
            "Paracetamol",
            "Benzene",
            "Ethanol",
        ] {
            assert!(
                eps.contains(label),
                "EPS grid missing compound label {label}"
            );
        }
        let png = std::fs::read(wdir.join("grid.png"))?;
        assert!(png.starts_with(b"\x89PNG\r\n\x1a\n"));
        assert!(png.len() > 2000);

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_search_substructure_rejects_oversized_csv() -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_oversized_csv_test");
        let _ = std::fs::remove_dir_all(&base);
        let rdir = base.join("read");
        std::fs::create_dir_all(&rdir)?;

        let mut big = String::from("name,smiles\n");
        for i in 0..80_000 {
            big.push_str(&format!("Mol_{i},CC(=O)Oc1ccccc1C(=O)O\n"));
        }
        assert!(big.len() > 2 * 1024 * 1024);
        std::fs::write(rdir.join("huge.csv"), big)?;

        let sb = PathSandbox::new(Some(&rdir), None)?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{"name":"search_substructure","arguments":{"file":"huge.csv","smarts":"c1ccccc1"}}
            })
            .to_string(),
        ]
        .join("\n")
            + "\n";
        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(req.as_bytes()), &mut out_buf, &sb, None)?;
        let lines: Vec<serde_json::Value> = String::from_utf8(out_buf)?
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1]["result"]["isError"], true);

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_draw_from_smiles_1000_atoms_avoids_quadratic_exponential_dos()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_1000_atom_dos_test");
        let _ = std::fs::remove_dir_all(&base);
        std::fs::create_dir_all(&base)?;

        // 1,000-atom poly-chiral/diene branched cyclic SMILES (must finish in < 2s)
        let mut big_smi = String::from("N[C@@H](C(=O)O)");
        for _ in 0..80 {
            big_smi.push_str("/C=C/C[C@@H](c1ccccc1)");
        }
        let sb = PathSandbox::new(None, None)?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{"name":"draw_from_smiles","arguments":{"smiles":big_smi}}
            })
            .to_string(),
        ]
        .join("\n")
            + "\n";

        let t0 = std::time::Instant::now();
        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(req.as_bytes()), &mut out_buf, &sb, None)?;
        let elapsed = t0.elapsed();
        let lines: Vec<serde_json::Value> = String::from_utf8(out_buf)?
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[1]["result"]["isError"], false);
        assert!(
            elapsed.as_secs_f64() < 2.5,
            "1000-atom MCP draw_from_smiles took too long: {elapsed:?}"
        );

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }

    #[test]
    fn mcp_build_peptide_export_3d_and_chimerax_bundle_work_inside_sandbox()
    -> Result<(), Box<dyn std::error::Error>> {
        let base = std::env::temp_dir().join("chemcraft_mcp_peptide_chimerax_test");
        let _ = std::fs::remove_dir_all(&base);
        let wdir = base.join("write");
        std::fs::create_dir_all(&wdir)?;

        let sb = PathSandbox::new(None, Some(&wdir))?;
        let req = [
            serde_json::json!({"jsonrpc":"2.0","id":0,"method":"initialize","params":{}})
                .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":1,"method":"tools/call",
                "params":{
                    "name":"build_peptide",
                    "arguments":{
                        "sequence":"Ac-Ala-Gly-Ser-NH2",
                        "secondary_structure":"alpha-helix",
                        "output_pdb":"pep.pdb",
                        "output_svg":"pep.svg"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":2,"method":"tools/call",
                "params":{
                    "name":"export_3d",
                    "arguments":{
                        "smiles":"CC(=O)Oc1ccccc1C(=O)O",
                        "output_path":"aspirin_3d.sdf",
                        "format":"sdf3d"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":3,"method":"tools/call",
                "params":{
                    "name":"export_chimerax",
                    "arguments":{
                        "peptide":"Tyr-Gly-Gly-Phe-Met",
                        "secondary_structure":"beta-strand",
                        "stem":"enkephalin"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":4,"method":"tools/call",
                "params":{
                    "name":"export_3d",
                    "arguments":{
                        "smiles":"CCO",
                        "output_path":"../outside.pdb",
                        "format":"pdb"
                    }
                }
            })
            .to_string(),
            serde_json::json!({
                "jsonrpc":"2.0","id":5,"method":"tools/call",
                "params":{
                    "name":"screen_library",
                    "arguments":{
                        "preset":"fragments",
                        "funnel":"ro3",
                        "diverse_n":8,
                        "output_sdf":"frags_3d.sdf"
                    }
                }
            })
            .to_string(),
        ]
        .join("\n")
            + "\n";

        let mut out_buf = Vec::new();
        mcp::serve_stream_for_test(BufReader::new(req.as_bytes()), &mut out_buf, &sb, None)?;
        let lines: Vec<serde_json::Value> = String::from_utf8(out_buf)?
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(serde_json::from_str)
            .collect::<Result<_, _>>()?;
        assert_eq!(lines.len(), 6);
        assert_eq!(lines[1]["result"]["isError"], false);
        assert_eq!(lines[2]["result"]["isError"], false);
        assert_eq!(lines[3]["result"]["isError"], false);
        // Path traversal rejected by sandbox
        assert_eq!(lines[4]["result"]["isError"], true);
        assert_eq!(lines[5]["result"]["isError"], false);
        assert_eq!(lines[5]["result"]["structuredContent"]["selected_count"], 8);

        assert!(wdir.join("pep.pdb").exists());
        assert!(wdir.join("pep.svg").exists());
        assert!(wdir.join("aspirin_3d.sdf").exists());
        assert!(wdir.join("frags_3d.sdf").exists());
        assert!(wdir.join("enkephalin.cxc").exists());
        assert!(wdir.join("enkephalin.pdb").exists());
        let cxc = std::fs::read_to_string(wdir.join("enkephalin.cxc"))?;
        assert!(cxc.contains("open enkephalin.pdb"));
        assert!(cxc.contains("cartoon"));

        let _ = std::fs::remove_dir_all(&base);
        Ok(())
    }
}
