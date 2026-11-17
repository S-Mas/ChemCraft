#![forbid(unsafe_code)]

mod mcp;

use std::path::{Path, PathBuf};

use clap::{Parser, Subcommand, ValueEnum};
use mol_admet_bridge::{AdmetCraftProvider, AdmetProvider, AdmetReportFormat};
use mol_engine::PathSandbox;
use mol_graph::{parse_smiles, write_sdf_v2000};
use mol_sar::{SarFilter, SarTable, find_activity_cliffs, find_matched_molecular_pairs};
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
                table.enrich_with_admet(&AdmetCraftProvider);
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
            let _ = with_admet;
            table.enrich_with_admet(&AdmetCraftProvider);
            let evals: Vec<_> = table
                .rows
                .into_iter()
                .filter_map(|r| r.admet_full)
                .collect();
            let fmt = match format {
                EnrichFormat::Csv => AdmetReportFormat::Csv,
                EnrichFormat::Json => AdmetReportFormat::Json,
                EnrichFormat::Pdf => AdmetReportFormat::Pdf,
            };
            let mut f = std::fs::File::create(&out)?;
            AdmetCraftProvider.export_report_to_writer(&evals, &mut f, fmt)?;
            eprintln!("Enriched {} compounds -> {}", evals.len(), out.display());
        }
        SubCmd::Mmp {
            input,
            min_tanimoto,
            out,
        } => {
            let mut table = load_table_from_path(&input)?;
            table.enrich_with_admet(&AdmetCraftProvider);
            let cliffs = find_activity_cliffs(&table.rows, min_tanimoto);
            let mmps = find_matched_molecular_pairs(&table.rows);
            let payload = serde_json::to_string_pretty(&serde_json::json!({
                "compounds": table.rows.len(),
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
            let eval = AdmetCraftProvider.evaluate_smiles(&smiles, None)?;
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
    }
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
}
