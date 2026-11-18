//! Optional In-Silico Chemical Library Manager & Virtual Screening Funnel window for `ChemCraft`.

use std::path::PathBuf;
use std::process::{Command, Stdio};

use eframe::egui;
use egui::RichText;
use mol_engine::{
    BuiltinLibraryPreset, DEFAULT_MAX_LIBRARY_RECORDS, LibraryFilterConfig, LibraryScreeningResult,
    ScreeningFunnel, default_dockcraft_handoff_dir, default_library_dir, discover_local_libraries,
    export_library_hits_sdf3d, filter_and_prepare_library, is_dockcraft_cli_available,
    launch_dockcraft_for_mols, load_builtin_library, materialize_builtin_library,
    stream_load_library,
};
use mol_graph::MolGraph;
use mol_sar::SarTable;

pub struct LibraryWindowState {
    pub open: bool,
    pub selected_preset: BuiltinLibraryPreset,
    pub custom_path: String,
    pub use_custom_path: bool,
    pub filter: LibraryFilterConfig,
    pub diversity_count: usize,
    pub local_files: Vec<PathBuf>,
    pub last_result: Option<LibraryScreeningResult>,
    pub receptor_path: String,
    pub box_center: [f64; 3],
    pub box_size: [f64; 3],
    pub dockcraft_cli_available: bool,
}

impl Default for LibraryWindowState {
    fn default() -> Self {
        let lib_dir = default_library_dir();
        let local_files = discover_local_libraries(&lib_dir);
        Self {
            open: false,
            selected_preset: BuiltinLibraryPreset::FdaApprovedMedicines,
            custom_path: String::new(),
            use_custom_path: false,
            filter: LibraryFilterConfig::default(),
            diversity_count: 0,
            local_files,
            last_result: None,
            receptor_path: String::new(),
            box_center: [0.0, 0.0, 0.0],
            box_size: [22.0, 22.0, 22.0],
            dockcraft_cli_available: is_dockcraft_cli_available(),
        }
    }
}

pub enum LibraryAction {
    None,
    LoadHitOnCanvas(MolGraph),
    SendResultToSar(SarTable),
    Status(String),
}

pub fn render_library_window(
    ctx: &egui::Context,
    state: &mut LibraryWindowState,
    export_dir: &str,
) -> LibraryAction {
    if !state.open {
        return LibraryAction::None;
    }

    let mut action = LibraryAction::None;
    let mut open = state.open;

    egui::Window::new("Screening Libraries (In Silico Pipeline)")
        .open(&mut open)
        .resizable(true)
        .default_width(680.0)
        .default_height(520.0)
        .show(ctx, |ui| {
            ui.label(
                RichText::new(
                    "Stream .sdf/.sdf.gz/.smi/.smi.gz/.csv/.dwar -> ADMET/SMARTS funnel -> pH 7.4 + 3D -> DockCraft",
                )
                .small()
                .weak(),
            );
            ui.separator();

            // 1. Library Source Selection
            ui.horizontal_wrapped(|ui| {
                ui.radio_value(&mut state.use_custom_path, false, "Built-in Preset:");
                for p in [
                    BuiltinLibraryPreset::FdaApprovedMedicines,
                    BuiltinLibraryPreset::FragmentRo3,
                    BuiltinLibraryPreset::KinaseAndGpcrScaffolds,
                ] {
                    if ui
                        .selectable_label(
                            !state.use_custom_path && state.selected_preset == p,
                            p.label(),
                        )
                        .clicked()
                    {
                        state.selected_preset = p;
                        state.use_custom_path = false;
                    }
                }
            });

            ui.horizontal(|ui| {
                ui.radio_value(&mut state.use_custom_path, true, "Local File:");
                ui.add(
                    egui::TextEdit::singleline(&mut state.custom_path)
                        .desired_width(320.0)
                        .hint_text("~/.local/share/chemcraft/libraries/zinc_subset.sdf.gz"),
                );
                if ui.button("Save Presets to Library Dir").clicked() {
                    let dir = default_library_dir();
                    let mut saved = 0usize;
                    for p in [
                        BuiltinLibraryPreset::FdaApprovedMedicines,
                        BuiltinLibraryPreset::FragmentRo3,
                        BuiltinLibraryPreset::KinaseAndGpcrScaffolds,
                    ] {
                        if materialize_builtin_library(p, &dir).is_ok() {
                            saved += 1;
                        }
                    }
                    state.local_files = discover_local_libraries(&dir);
                    action = LibraryAction::Status(format!(
                        "Saved {saved} curated .smi libraries into {}",
                        dir.display()
                    ));
                }
            });

            if !state.local_files.is_empty() {
                ui.horizontal_wrapped(|ui| {
                    ui.label(RichText::new("Discovered in library dir:").small().weak());
                    for lf in state.local_files.clone() {
                        let fname = lf
                            .file_name()
                            .map(|s| s.to_string_lossy().to_string())
                            .unwrap_or_default();
                        if ui.small_button(fname).clicked() {
                            state.custom_path = lf.display().to_string();
                            state.use_custom_path = true;
                        }
                    }
                });
            }

            ui.separator();

            // 2. Funnel & Preparation Controls
            ui.horizontal_wrapped(|ui| {
                ui.strong("Pre-Filter Funnel:");
                for f in [
                    ScreeningFunnel::All,
                    ScreeningFunnel::RuleOf5,
                    ScreeningFunnel::LeadLike,
                    ScreeningFunnel::FragmentRo3,
                    ScreeningFunnel::BbbCns,
                ] {
                    ui.selectable_value(&mut state.filter.funnel, f, f.label());
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.label("SMARTS Substructure:");
                let mut sm = state
                    .filter
                    .smarts_substructure
                    .clone()
                    .unwrap_or_default();
                if ui
                    .add(
                        egui::TextEdit::singleline(&mut sm)
                            .desired_width(140.0)
                            .hint_text("e.g. c1ccccc1"),
                    )
                    .changed()
                {
                    state.filter.smarts_substructure =
                        if sm.trim().is_empty() { None } else { Some(sm) };
                }
                ui.label("MaxMin Diverse K:");
                ui.add(egui::DragValue::new(&mut state.diversity_count).range(0..=500));
                ui.label("(0=off)");
                ui.checkbox(&mut state.filter.protonate_ph74, "Ionize pH 7.4");
                ui.checkbox(&mut state.filter.generate_3d, "Build 3D");
            });

            ui.horizontal_wrapped(|ui| {
                if ui.button("Filter & Prepare 3D Library").clicked() {
                    state.filter.diversity_pick_n = if state.diversity_count > 0 {
                        Some(state.diversity_count)
                    } else {
                        None
                    };
                    let load_res = if state.use_custom_path {
                        let p = PathBuf::from(state.custom_path.trim());
                        stream_load_library(&p, DEFAULT_MAX_LIBRARY_RECORDS)
                            .map(|mols| (p.display().to_string(), mols))
                    } else {
                        Ok((
                            state.selected_preset.label().to_string(),
                            load_builtin_library(state.selected_preset),
                        ))
                    };

                    match load_res {
                        Ok((label, mols)) => {
                            match filter_and_prepare_library(&label, mols, &state.filter) {
                                Ok(res) => {
                                    action = LibraryAction::Status(format!(
                                        "Library '{}': {}/{} passed funnel -> {} ready in {} ms",
                                        res.source_label,
                                        res.passed_funnel,
                                        res.total_loaded,
                                        res.hits.len(),
                                        res.elapsed_ms
                                    ));
                                    state.last_result = Some(res);
                                }
                                Err(e) => {
                                    action = LibraryAction::Status(format!("Library filter: {e}"));
                                }
                            }
                        }
                        Err(e) => {
                            action = LibraryAction::Status(format!("Library load error: {e}"));
                        }
                    }
                }

                if let Some(ref res) = state.last_result {
                    if ui.button("Send Hits to DockCraft GUI").clicked() {
                        let refs: Vec<&MolGraph> = res.hits.iter().map(|h| &h.mol).collect();
                        match launch_dockcraft_for_mols(&refs, None, Some("chemcraft_screen_lib"))
                        {
                            Ok(b) => {
                                action = LibraryAction::Status(format!(
                                    "Sent {} prepared 3D ligands to DockCraft -> {}",
                                    b.molecule_count,
                                    b.sdf_path.display()
                                ));
                            }
                            Err(e) => action = LibraryAction::Status(format!("DockCraft: {e}")),
                        }
                    }
                    if ui.button("Export 3D SDF").clicked() {
                        let out = PathBuf::from(export_dir).join("chemcraft_filtered_3d.sdf");
                        match export_library_hits_sdf3d(res, &out) {
                            Ok(n) => {
                                action = LibraryAction::Status(format!(
                                    "Exported {n} prepared 3D compounds -> {}",
                                    out.display()
                                ));
                            }
                            Err(e) => action = LibraryAction::Status(format!("Export error: {e}")),
                        }
                    }
                    if ui.button("Load into SAR Explorer").clicked() {
                        let mut table = SarTable::new();
                        for h in &res.hits {
                            table.add_molecule(h.mol.clone(), None, Some(res.source_label.clone()));
                        }
                        action = LibraryAction::SendResultToSar(table);
                    }
                }
            });

            // Optional direct headless DockCraft screening runner when `dockcraft-cli` is installed
            if state.dockcraft_cli_available && state.last_result.is_some() {
                ui.collapsing("Run Headless DockCraft Virtual Screening (`dockcraft-cli`)", |ui| {
                    ui.horizontal(|ui| {
                        ui.label("Receptor PDB:");
                        ui.add(
                            egui::TextEdit::singleline(&mut state.receptor_path)
                                .desired_width(250.0)
                                .hint_text("/path/to/receptor.pdb"),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.label("Center:");
                        for c in &mut state.box_center {
                            ui.add(egui::DragValue::new(c).speed(0.2));
                        }
                        ui.label("Box:");
                        for s in &mut state.box_size {
                            ui.add(egui::DragValue::new(s).range(6.0..=60.0).speed(0.2));
                        }
                        if ui.button("Run Screening").clicked()
                            && !state.receptor_path.trim().is_empty()
                            && let Some(ref res) = state.last_result
                        {
                            let handoff_dir = default_dockcraft_handoff_dir();
                                let lib_sdf = handoff_dir.join("chemcraft_screen_input_3d.sdf");
                                let out_csv = PathBuf::from(export_dir).join("dockcraft_screening_hits.csv");
                                match export_library_hits_sdf3d(res, &lib_sdf) {
                                    Ok(_) => {
                                        let center_arg = format!(
                                            "{},{},{}",
                                            state.box_center[0], state.box_center[1], state.box_center[2]
                                        );
                                        let size_arg = format!(
                                            "{},{},{}",
                                            state.box_size[0], state.box_size[1], state.box_size[2]
                                        );
                                        match Command::new("dockcraft-cli")
                                            .args([
                                                "screen",
                                                "--receptor",
                                                state.receptor_path.trim(),
                                                "--library",
                                            ])
                                            .arg(&lib_sdf)
                                            .arg("--out")
                                            .arg(&out_csv)
                                            .args(["--center", &center_arg, "--size", &size_arg])
                                            .stdin(Stdio::null())
                                            .stdout(Stdio::piped())
                                            .stderr(Stdio::piped())
                                            .output()
                                        {
                                            Ok(o) if o.status.success() => {
                                                action = LibraryAction::Status(format!(
                                                    "DockCraft screening complete -> {}",
                                                    out_csv.display()
                                                ));
                                            }
                                            Ok(o) => {
                                                let err = String::from_utf8_lossy(&o.stderr);
                                                action = LibraryAction::Status(format!(
                                                    "dockcraft-cli screen failed: {}",
                                                    err.trim()
                                                ));
                                            }
                                            Err(e) => {
                                                action = LibraryAction::Status(format!(
                                                    "dockcraft-cli spawn failed: {e}"
                                                ));
                                            }
                                        }
                                    }
                                    Err(e) => {
                                        action = LibraryAction::Status(format!("Library write error: {e}"));
                                    }
                                }
                            }
                        });
                });
            }

            ui.separator();

            // 3. Screened Hits Table
            if let Some(ref res) = state.last_result {
                ui.strong(format!(
                    "Hits: {} selected ({} passed filter / {} loaded) in {} ms",
                    res.hits.len(),
                    res.passed_funnel,
                    res.total_loaded,
                    res.elapsed_ms
                ));
                egui::ScrollArea::both().show(ui, |ui| {
                    egui::Grid::new("library_hits_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            for h in ["Action", "Compound", "Formula", "MW", "cLogP", "TPSA", "GI", "BBB", "Lipinski"] {
                                ui.strong(h);
                            }
                            ui.end_row();
                            for hit in res.hits.iter().take(200) {
                                if ui.small_button("View 2D+3D").clicked() {
                                    action = LibraryAction::LoadHitOnCanvas(hit.mol.clone());
                                }
                                ui.label(&hit.name);
                                ui.monospace(&hit.eval.physicochemical.formula);
                                ui.label(format!("{:.1}", hit.eval.physicochemical.molecular_weight));
                                ui.label(format!("{:.2}", hit.eval.lipophilicity.consensus_log_p));
                                ui.label(format!("{:.1}", hit.eval.physicochemical.tpsa));
                                ui.label(hit.eval.pharmacokinetics.gi_absorption.as_str());
                                ui.label(if hit.eval.pharmacokinetics.bbb_permeant { "Yes" } else { "No" });
                                ui.label(format!("{}v", hit.eval.druglikeness.lipinski.violations));
                                ui.end_row();
                            }
                        });
                });
            } else {
                ui.label("Choose a library preset or local file and click 'Filter & Prepare 3D Library'.");
            }
        });

    state.open = open;
    action
}
