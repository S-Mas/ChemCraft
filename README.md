# ChemCraft

**Local-First, Clean-Room 2D Vector + 3D Molecular Design, Peptide Builder & SAR Studio in Pure Rust**  
*(Works side-by-side with [`AdmetCraft`](../AdmetCraft) for deep ADMET/PK & [`ChimeraX`](https://www.cgl.ucsf.edu/chimerax/) for macromolecular 3D)*

`ChemCraft` combines ChemDraw/VectorCraft-style 2D chemical drawing, interactive pure-Rust 3D molecular visualization, a 2D+3D canonical L-peptide builder, DataWarrior-style in-memory SAR exploration (ECFP4 Tanimoto, SMARTS, Matched Molecular Pairs & Activity Cliffs), one-click **UCSF ChimeraX** 3D handoff, and live `< 5 ms` ADMET previews with optional delegation to external **`AdmetCraft`** — all written in **100% safe, offline pure Rust** (`#![forbid(unsafe_code)]`, zero `.unwrap()` / `.expect()` / `panic!()`).

---

## ChemCraft vs AdmetCraft vs ChimeraX

| Tool | Primary Focus | Interoperability |
|---|---|---|
| **ChemCraft** | 2D vector drawing, ring snapping, 3D interactive preview, simple 2D+3D L-peptides, MMP/SAR exploration, SVG/EPS/PDF/PNG/PDB export | Sends SMILES/SAR libraries to **AdmetCraft** & exports 3D `.pdb`/`.sdf` + `.cxc` bundles to **ChimeraX** |
| **AdmetCraft** | Dedicated ADMET, pharmacokinetics (GI/BBB/P-gp/CYP), ESOL/Ali/Silicos-IT solubility, Lipinski/Ghose/Veber/Egan/Muegge, PAINS/Brenk/SA, Bioavailability Radar & BOILED-Egg SVGs | Invoked optionally from ChemCraft via `admetcraft-cli eval/batch` or launched standalone (`admetcraft`) |
| **ChimeraX** | Full macromolecular 3D rendering, publication cartoons, molecular surfaces, electrostatic potentials & H-bond networks | Opened directly from ChemCraft (`Visualizar en ChimeraX` / `chemcraft-cli chimerax`) using sandbox-safe relative `.cxc` bundles |

---

## Architecture (7-Crate Cargo Workspace)

| Layer | Crate / App | Role |
|---|---|---|
| L0 | `crates/mol-graph` | Clean-room editable `MolGraph`, SSSR ring basis, iterative Hückel aromaticity, valence diagnostics, Cahn-Ingold-Prelog `R/S` & `E/Z` stereo, canonical SMILES, SMARTS subgraph matching, SDF **V2000 & V3000**, DataWarrior `.dwar`, CSV, bioisostere mutation, pure-Rust **3D conformer generator & PDB/3D-SDF writer**, and **20-AA 2D+3D Peptide Builder** |
| L1 | `crates/mol-vector` | 2D layout generator (fused/spiro/bridged rings + `Clean Structure 2D`), ring-snapping templates, ring-aligned double/triple/aromatic bonds, stereo wedges/hashes/wavy/dative arrows, label collision resolution, and **SVG / EPS / FlateDecode PDF / PNG** export (single or multi-molecule grid) |
| L2 | `crates/mol-admet-bridge` | Decoupled `AdmetProvider` trait + standalone pure-Rust ADMET preview engine, external `admetcraft-cli` delegation (`eval`/`batch` + BOILED-Egg SVG handoff), and `admetcraft` GUI launcher — zero sibling path dependencies |
| L2 | `crates/mol-sar` | In-memory SAR compound table, 2048-bit ECFP4 (`MorganFp`) fingerprints, `rayon`-parallel Tanimoto/SMARTS/numeric filtering, deterministic Matched Molecular Pairs (MMP), duplicate detection, Activity Cliffs (`SALI`), and Ligand Efficiency (`LE`, `LLE`) |
| L3 | `crates/mol-engine` | Session state, Copy-on-Write structural Undo/Redo (`Arc` snapshots), `CommandRegistry` (including `BuildPeptide`), bounded `.chemcraft` ZIP+JSON archives, **ChimeraX `.pdb`/`.sdf`+`.cxc` bundle launcher**, **DockCraft 3D SDF/SMILES handoff**, and kernel-confined `PathSandbox` (`openat2` with `RESOLVE_BENEATH \| RESOLVE_NO_SYMLINKS`) |
| L4 | `apps/chemcraft-cli` | Headless CLI (`render`, `peptide`, `export3d`, `chimerax`, `dockcraft`, `filter`, `enrich`, `mmp`, `inspect`) + strict JSON-RPC 2.0 **MCP server** (`stdio` & `--tcp`) |
| L5 | `apps/chemcraft` | Native `egui` / `eframe` desktop studio focused on small-molecule & peptide 2D/3D design (`2D Vector` / `3D Interactive` / `Split 2D + 3D`), optional toggleable **SAR Explorer**, **"Visualizar en ChimeraX"**, **"Enviar a DockCraft"**, Command Palette (`Ctrl+K`), and AdmetCraft Inspector + handoff |

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

# Build a canonical L-peptide in 2D + 3D and export SVG + PDB
chemcraft-cli peptide "Ac-Ala-Gly-Tyr-Ser-NH2" --ss alpha --svg peptide.svg --pdb peptide.pdb

# Export a 3D conformer from SMILES to PDB or 3D SDF V2000
chemcraft-cli export3d "CC(=O)Oc1ccccc1C(=O)O" --out aspirin_3d.pdb --format pdb

# Export a ChimeraX bundle (.pdb/.sdf + .cxc) and launch ChimeraX
chemcraft-cli chimerax --peptide "Tyr-Gly-Gly-Phe-Met" --ss beta --launch

# Export 3D ligand or SDF library to DockCraft and optionally launch DockCraft
chemcraft-cli dockcraft --smiles "CC(C)Cc1ccc(C(C)C(=O)O)cc1" --launch
chemcraft-cli dockcraft --library hits.sdf --out-dir ./dock_handoff

# Stream, pre-filter (Ro5/Lead-Like/Ro3/BBB/MaxMin), protonate pH 7.4, and prepare 3D SDF libraries
chemcraft-cli library --preset fragments --funnel ro3 --diverse-n 15 --out frags_3d.sdf --send-dockcraft
chemcraft-cli library --input ~/.local/share/chemcraft/libraries/zinc_subset.sdf.gz --funnel lead --out lead_hits_3d.sdf

# Hardened MCP Server (stdio or loopback TCP with token)
chemcraft-cli mcp --read-root ./data --write-root ./exports
chemcraft-cli mcp --read-root ./data --write-root ./exports --tcp 127.0.0.1:7878 --mcp-token "$TOKEN"
```

### MCP Tools Exposed
- `draw_from_smiles`: Computes clean 2D coordinates, valence warnings, deterministic open connectivity key (`CCKEY-...`), and SVG output.
- `mutate_scaffold`: Replaces carboxylic acids (`-COOH` -> `tetrazole`) or attaches substituents and returns the **AdmetCraft** delta (`ΔMW`, `ΔLogP`, `ΔTPSA`, `ΔLogS`, `GI`, `BBB`).
- `search_substructure`: Filters SDF/DWAR/CSV libraries inside `--read-root` by SMARTS and ECFP4 Tanimoto.
- `export_vector_figure`: Exports single molecules or multi-molecule SAR grids (`svg`, `eps`, `pdf`, `png`) inside `--write-root`.
- `mol.analyze_admet`: Runs the built-in ADMET engine (or optional `admetcraft-cli`) and optionally writes CSV/JSON/PDF ADMET reports.
- `build_peptide`: Builds 2D+3D peptides (`Beta-Strand`, `Alpha-Helix`, `Extended`) from 1-letter or 3-letter sequences with optional `.pdb`/`.svg` export inside `--write-root`.
- `export_3d`: Generates 3D conformers from SMILES or peptides and writes PDB or 3D SDF inside `--write-root`.
- `export_chimerax`: Exports self-contained ChimeraX `.pdb`/`.sdf` + `.cxc` scripts inside `--write-root`.
- `screen_library`: Streams built-in presets or sandboxed `.sdf`/`.sdf.gz`/`.smi`/`.smi.gz`/`.csv`/`.dwar` chemical libraries through ADMET/SMARTS/MaxMin funnels, ionizes at pH 7.4, generates 3D coordinates, and writes DockCraft-ready 3D SDF files inside `--write-root`.

## Clean-Room & Scientific Scope Notes
- Zero third-party proprietary code (`ChemDraw`, `MarvinSketch`, `DataWarrior`).
- `structural_key()` produces a deterministic topological `CCKEY-...` digest from canonical SMILES (does not invoke external C InChI binaries).
- Live ADMET descriptors run cleanly standalone inside `mol-admet-bridge` (with optional `admetcraft-cli` IPC discovery when enabled).
