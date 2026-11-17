# ChemCraft

**Local-First, Clean-Room 2D Vector Molecular Design & SAR Studio in Pure Rust**  
*(Integrated in-process with [`AdmetCraft`](../AdmetCraft))*

`ChemCraft` combines ChemDraw/VectorCraft-style 2D chemical drawing, DataWarrior-style in-memory SAR exploration (ECFP4 Tanimoto, SMARTS, Matched Molecular Pairs & Activity Cliffs), and live `< 5 ms` ADMET/Bioavailability predictions from `AdmetCraft` — all written in **100% safe, offline pure Rust** (`#![forbid(unsafe_code)]`, zero `.unwrap()` / `.expect()` / `panic!()`).

---

## Architecture (7-Crate Cargo Workspace)

| Layer | Crate / App | Role |
|---|---|---|
| L0 | `crates/mol-graph` | Clean-room editable `MolGraph`, SSSR ring basis, Hückel aromaticity, valence diagnostics, E/Z & R/S hints, canonical SMILES, SMARTS subgraph matching, SDF **V2000 & V3000**, DataWarrior `.dwar`, CSV, bioisostere & substituent mutation |
| L1 | `crates/mol-vector` | 2D layout generator (fused/spiro/bridged rings + `Clean Structure 2D`), ring-snapping templates, ring-aligned double/triple/aromatic bonds, stereo wedges/hashes/wavy/dative arrows, label collision resolution, and **SVG / EPS / FlateDecode PDF / PNG** export |
| L2 | `crates/mol-admet-bridge` | Decoupled `AdmetProvider` trait + direct `AdmetCraftProvider` linking `../AdmetCraft/crates/{admet-calc,admet-report}` |
| L2 | `crates/mol-sar` | In-memory SAR compound table, 2048-bit ECFP4 (`MorganFp`) fingerprints, `rayon`-parallel Tanimoto/SMARTS/numeric filtering, Matched Molecular Pairs (MMP), Activity Cliffs (`SALI`), and Ligand Efficiency (`LE`, `LLE`) |
| L3 | `crates/mol-engine` | Session state, Copy-on-Write structural Undo/Redo (`Arc` snapshots), `CommandRegistry`, atomic `.chemcraft` ZIP+JSON project archives, and hardened `PathSandbox` (`symlink_metadata`, `O_NOFOLLOW | O_NONBLOCK`, no TOCTOU) |
| L4 | `apps/chemcraft-cli` | Headless CLI (`render`, `filter`, `enrich`, `mmp`, `inspect`) + strict JSON-RPC 2.0 **MCP server** (`stdio` & `--tcp`) |
| L5 | `apps/chemcraft` | Native `egui` / `eframe` desktop studio with drawing palette, interactive snapping canvas, SAR Explorer tabs (Table, Scatter Plot, MMP/Cliffs, Filters), Command Palette (`Ctrl+K`), and Live ADMET Inspector |

---

## Building & Running

```bash
cargo test --workspace --offline
cargo build --workspace --release --offline
./packaging/install.sh
```

### Desktop Studio
```bash
chemcraft
```

### CLI & MCP Examples
```bash
# Render vector SVG, PDF, EPS, or PNG
chemcraft-cli render "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C" --out testosterone.svg --format svg
chemcraft-cli render "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12" --out sildenafil.pdf --format pdf

# Substructure + Lipinski filter
chemcraft-cli filter library.sdf --substructure "c1ccccc1" --lipinski-only --out hits.sdf

# Enrich an SDF/DWAR/CSV dataset with AdmetCraft and export PDF/CSV
chemcraft-cli enrich library.sdf --with-admet --out admet_report.pdf --format pdf

# Hardened MCP Server (stdio or loopback TCP with token)
chemcraft-cli mcp --read-root ./data --write-root ./exports
chemcraft-cli mcp --read-root ./data --write-root ./exports --tcp 127.0.0.1:7878 --mcp-token "$TOKEN"
```

### MCP Tools Exposed
- `draw_from_smiles`: Computes clean 2D coordinates, valence warnings, deterministic open connectivity key (`CCKEY-...`), and SVG output.
- `mutate_scaffold`: Replaces carboxylic acids (`-COOH` -> `tetrazole`) or attaches substituents and returns the **AdmetCraft** delta (`ΔMW`, `ΔLogP`, `ΔTPSA`, `ΔLogS`, `GI`, `BBB`).
- `search_substructure`: Filters SDF/DWAR/CSV libraries inside `--read-root` by SMARTS and ECFP4 Tanimoto.
- `export_vector_figure`: Exports single molecules or multi-molecule SAR grids (`svg`, `eps`, `pdf`, `png`) inside `--write-root`.
- `mol.analyze_admet`: Runs the in-process `AdmetCraft` engine and optionally writes CSV/JSON/PDF ADMET reports.

## Clean-Room & Scientific Scope Notes
- Zero third-party proprietary code (`ChemDraw`, `MarvinSketch`, `DataWarrior`).
- `structural_key()` produces a deterministic topological `CCKEY-...` digest from canonical SMILES (does not invoke external C InChI binaries).
- Live ADMET descriptors inherit `AdmetCraft v0.2.0` (`admet-calc` & `admet-report`).
