#![forbid(unsafe_code)]

mod viewer3d;
mod widgets;

use std::path::PathBuf;

use eframe::egui;
use egui::{Color32, Pos2, RichText, Sense, Stroke, Vec2};
use mol_admet_bridge::{
    AdmetCraftProvider, AdmetProvider, AdmetReportFormat, launch_admetcraft_gui,
};
use mol_engine::{
    ChemCraftDocument, ChemSession, Command, CommandRegistry, MAX_READ_BYTES, PathSandbox,
    decode_chemcraft_zip, is_chimerax_available, launch_chimerax_for_mol, save_chemcraft_atomic,
};
use mol_graph::{
    BondOrder, BondStereo, SecondaryStructure, canonical_amino_acids, generate_3d_coords,
    parse_csv_records, parse_dwar, parse_sdf, write_pdb, write_sdf_v2000_3d,
};
use mol_vector::{
    ReactionArrow, RingTemplate, SceneStyle, VectorScene, build_mol_scene, scene_to_eps,
    scene_to_pdf_bytes, scene_to_png_bytes, scene_to_svg,
};
use viewer3d::{Viewer3dState, render_3d_viewport};
use widgets::{
    ScatterXAxis, draw_admet_radar, draw_sar_scatter, paint_mol_thumbnail, paint_vector_scene,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ViewMode {
    View2D,
    View3D,
    Split2D3D,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DrawTool {
    SelectMove,
    BondOrder(BondOrder),
    BondStereo(BondStereo),
    Ring(RingTemplate),
    ReactionArrow,
    Eraser,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SarTab {
    Table,
    Scatter,
    MmpCliffs,
    Filters,
}

struct ChemCraftApp {
    session: ChemSession,
    tool: DrawTool,
    view_mode: ViewMode,
    viewer_3d: Viewer3dState,
    peptide_seq: String,
    peptide_ss: SecondaryStructure,
    peptide_drawer_open: bool,
    chimerax_available: bool,
    admet_cli_available: bool,
    admet_gui_available: bool,
    active_element: &'static str,
    smiles_bar: String,
    file_path_bar: String,
    export_dir: String,
    status_msg: String,
    dark_mode: bool,
    sar_tab: SarTab,
    scatter_x: ScatterXAxis,
    palette_open: bool,
    palette_query: String,
    hover_atom: Option<usize>,
    drag_start_atom: Option<usize>,
    last_canvas_scene: Option<VectorScene>,
    arrow_start: Option<[f64; 2]>,
}

impl ChemCraftApp {
    fn new() -> Self {
        let mut session = ChemSession::new();
        let _ = CommandRegistry::execute(&mut session, Command::LoadDemo50Sar);
        let export_dir = std::env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .display()
            .to_string();
        let smiles_bar = session.canvas_mol.canonical_smiles();
        let chimerax_available = is_chimerax_available();
        let admet_cli_available = AdmetCraftProvider::is_external_cli_available();
        let admet_gui_available = AdmetCraftProvider::is_external_gui_available();
        Self {
            session,
            tool: DrawTool::BondOrder(BondOrder::Single),
            view_mode: ViewMode::Split2D3D,
            viewer_3d: Viewer3dState::default(),
            peptide_seq: "Ac-Ala-Gly-Tyr-Ser-NH2".to_string(),
            peptide_ss: SecondaryStructure::AlphaHelix,
            peptide_drawer_open: false,
            chimerax_available,
            admet_cli_available,
            admet_gui_available,
            active_element: "C",
            smiles_bar,
            file_path_bar: String::new(),
            export_dir,
            status_msg:
                "Ready — 2D Vector + 3D Conformer + Peptide Builder + ChimeraX + AdmetCraft"
                    .to_string(),
            dark_mode: false,
            sar_tab: SarTab::Table,
            scatter_x: ScatterXAxis::ConsensusLogP,
            palette_open: false,
            palette_query: String::new(),
            hover_atom: None,
            drag_start_atom: None,
            last_canvas_scene: None,
            arrow_start: None,
        }
    }

    fn run_cmd(&mut self, cmd: Command) {
        match CommandRegistry::execute(&mut self.session, cmd) {
            Ok(out) => {
                self.status_msg = out.message;
                self.smiles_bar = self.session.canvas_mol.canonical_smiles();
            }
            Err(e) => {
                self.status_msg = format!("Error: {e}");
            }
        }
    }

    fn export_canvas(&mut self, fmt: &str) {
        if self.session.canvas_mol.atoms.is_empty() {
            self.status_msg = "Canvas is empty".to_string();
            return;
        }
        let root = PathBuf::from(&self.export_dir);
        let sb = match PathSandbox::new(None, Some(&root)) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Export sandbox error: {e}");
                return;
            }
        };
        let scene = build_mol_scene(
            &self.session.canvas_mol,
            &SceneStyle {
                width: 640.0,
                height: 460.0,
                dark_mode: false,
                ..SceneStyle::default()
            },
            &self.session.reaction_arrows,
        );
        let filename = format!("chemcraft_structure.{fmt}");
        let bytes: Vec<u8> = match fmt {
            "svg" => scene_to_svg(&scene).into_bytes(),
            "eps" => scene_to_eps(&scene).into_bytes(),
            "pdf" => match scene_to_pdf_bytes(&scene) {
                Ok(b) => b,
                Err(e) => {
                    self.status_msg = format!("PDF error: {e}");
                    return;
                }
            },
            "png" => match scene_to_png_bytes(&scene, 2.0) {
                Ok(b) => b,
                Err(e) => {
                    self.status_msg = format!("PNG error: {e}");
                    return;
                }
            },
            _ => return,
        };
        match sb.open_write(std::path::Path::new(&filename)) {
            Ok((resolved, mut f)) => {
                use std::io::Write;
                if let Err(e) = f.write_all(&bytes) {
                    self.status_msg = format!("Write error: {e}");
                } else {
                    self.status_msg = format!("Exported {}", resolved.display());
                }
            }
            Err(e) => self.status_msg = format!("Sandbox write error: {e}"),
        }
    }

    fn export_structure_3d(&mut self, fmt: &str) {
        if self.session.canvas_mol.atoms.is_empty() {
            self.status_msg = "Canvas is empty".to_string();
            return;
        }
        let root = PathBuf::from(&self.export_dir);
        let sb = match PathSandbox::new(None, Some(&root)) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Export sandbox error: {e}");
                return;
            }
        };
        let (filename, content) = match fmt {
            "pdb" => (
                "chemcraft_structure.pdb".to_string(),
                write_pdb(&self.session.canvas_mol),
            ),
            "sdf3d" => (
                "chemcraft_structure_3d.sdf".to_string(),
                write_sdf_v2000_3d(&[&self.session.canvas_mol]),
            ),
            _ => return,
        };
        match sb.open_write(std::path::Path::new(&filename)) {
            Ok((resolved, mut f)) => {
                use std::io::Write;
                if let Err(e) = f.write_all(content.as_bytes()) {
                    self.status_msg = format!("Write error: {e}");
                } else {
                    self.status_msg = format!("Exported 3D -> {}", resolved.display());
                }
            }
            Err(e) => self.status_msg = format!("Sandbox write error: {e}"),
        }
    }

    fn open_in_chimerax(&mut self) {
        match launch_chimerax_for_mol(&self.session.canvas_mol, None) {
            Ok(bundle) => {
                self.status_msg = format!("Launched ChimeraX ({})", bundle.script_path.display());
            }
            Err(e) => {
                self.status_msg = format!("ChimeraX: {e}");
            }
        }
    }

    fn collect_handoff_items(&self) -> Vec<(String, String)> {
        let mut items = Vec::new();
        if !self.session.canvas_mol.atoms.is_empty() {
            let smi = self.session.canvas_mol.canonical_smiles();
            if !smi.is_empty() {
                items.push((
                    smi,
                    self.session
                        .canvas_mol
                        .name
                        .clone()
                        .unwrap_or_else(|| "Canvas_Mol".to_string()),
                ));
            }
        }
        for row in &self.session.sar_table.rows {
            if !row.canonical_smiles.is_empty() {
                items.push((row.canonical_smiles.clone(), row.name.clone()));
            }
        }
        items
    }

    fn open_in_admetcraft(&mut self) {
        let items = self.collect_handoff_items();
        if items.is_empty() {
            self.status_msg = "Nothing to send to AdmetCraft".to_string();
            return;
        }
        let dir = PathBuf::from(&self.export_dir);
        match launch_admetcraft_gui(&items, &dir) {
            Ok(smi_path) => {
                self.status_msg = format!(
                    "Launched AdmetCraft GUI (handoff SMILES saved at {})",
                    smi_path.display()
                );
            }
            Err(e) => {
                self.status_msg = format!("AdmetCraft launch: {e}");
            }
        }
    }

    fn export_admet_boiled_egg(&mut self) {
        let items = self.collect_handoff_items();
        if items.is_empty() {
            self.status_msg = "No compounds available for BOILED-Egg plot".to_string();
            return;
        }
        let work_dir = PathBuf::from(&self.export_dir);
        let out_path = work_dir.join("admetcraft_boiled_egg.svg");
        match self.session.admet_provider().run_full_report_with_handoff(
            &items,
            &work_dir,
            &out_path,
            AdmetReportFormat::SvgBoiledEgg,
        ) {
            Ok(msg) => self.status_msg = msg,
            Err(e) => self.status_msg = format!("AdmetCraft BOILED-Egg error: {e}"),
        }
    }

    fn export_admet_report(&mut self, fmt: AdmetReportFormat, filename: &str) {
        let evals: Vec<_> = self
            .session
            .sar_table
            .rows
            .iter()
            .filter_map(|r| r.admet_full.clone())
            .collect();
        if evals.is_empty() {
            self.status_msg = "No ADMET profiles to export".to_string();
            return;
        }
        let root = PathBuf::from(&self.export_dir);
        let sb = match PathSandbox::new(None, Some(&root)) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Export root error: {e}");
                return;
            }
        };
        match sb.open_write(std::path::Path::new(filename)) {
            Ok((resolved, mut f)) => {
                match self
                    .session
                    .admet_provider()
                    .export_report_to_writer(&evals, &mut f, fmt)
                {
                    Ok(()) => {
                        self.status_msg = format!("Exported ADMET report -> {}", resolved.display())
                    }
                    Err(e) => self.status_msg = format!("ADMET export error: {e}"),
                }
            }
            Err(e) => self.status_msg = format!("Sandbox error: {e}"),
        }
    }

    fn render_peptide_drawer(&mut self, ctx: &egui::Context) {
        if !self.peptide_drawer_open {
            return;
        }
        let mut open = self.peptide_drawer_open;
        let mut trigger_build = false;
        egui::Window::new("Peptide Builder (2D + 3D + PDB)")
            .open(&mut open)
            .resizable(true)
            .default_width(410.0)
            .show(ctx, |ui| {
                ui.label(
                    RichText::new("Build canonical L-peptides in 2D + 3D from 1-letter or 3-letter sequences:")
                        .small(),
                );
                ui.horizontal(|ui| {
                    ui.label("Sequence:");
                    let edit = ui.add(
                        egui::TextEdit::singleline(&mut self.peptide_seq)
                            .desired_width(240.0)
                            .font(egui::TextStyle::Monospace),
                    );
                    if edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                        trigger_build = true;
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Secondary structure:");
                    ui.selectable_value(
                        &mut self.peptide_ss,
                        SecondaryStructure::BetaStrand,
                        "Beta-Strand",
                    );
                    ui.selectable_value(
                        &mut self.peptide_ss,
                        SecondaryStructure::AlphaHelix,
                        "Alpha-Helix",
                    );
                    ui.selectable_value(
                        &mut self.peptide_ss,
                        SecondaryStructure::Extended,
                        "Extended",
                    );
                });
                ui.separator();
                ui.label(RichText::new("Append residue:").small().weak());
                ui.horizontal_wrapped(|ui| {
                    for (c1, c3, full) in canonical_amino_acids() {
                        if ui
                            .small_button(format!("{c1} ({c3})"))
                            .on_hover_text(full)
                            .clicked()
                        {
                            let s = self.peptide_seq.trim().to_string();
                            if s.is_empty() {
                                self.peptide_seq = c3.to_string();
                            } else if let Some(prefix) = s.strip_suffix("-NH2") {
                                if prefix.contains('-') {
                                    self.peptide_seq = format!("{prefix}-{c3}-NH2");
                                } else {
                                    self.peptide_seq = format!("{prefix}{c1}-NH2");
                                }
                            } else if s.contains('-') {
                                self.peptide_seq = format!("{s}-{c3}");
                            } else {
                                self.peptide_seq = format!("{s}{c1}");
                            }
                        }
                    }
                });
                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Preset: Ac-AGYS-NH2").clicked() {
                        self.peptide_seq = "Ac-Ala-Gly-Tyr-Ser-NH2".to_string();
                        trigger_build = true;
                    }
                    if ui.button("Preset: RGD (cyclic/strand)").clicked() {
                        self.peptide_seq = "Arg-Gly-Asp".to_string();
                        self.peptide_ss = SecondaryStructure::BetaStrand;
                        trigger_build = true;
                    }
                    if ui.button("Preset: Poly-Ala Helix (A8)").clicked() {
                        self.peptide_seq = "AAAAAAAA".to_string();
                        self.peptide_ss = SecondaryStructure::AlphaHelix;
                        trigger_build = true;
                    }
                });
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    if ui.button("Build Peptide (2D + 3D)").clicked() {
                        trigger_build = true;
                    }
                    if ui.button("Clear Seq").clicked() {
                        self.peptide_seq.clear();
                    }
                });
            });
        self.peptide_drawer_open = open;
        if trigger_build {
            let seq = self.peptide_seq.clone();
            let ss = self.peptide_ss;
            self.run_cmd(Command::BuildPeptide {
                sequence: seq,
                secondary_structure: ss,
            });
            self.viewer_3d.reset();
            if self.view_mode == ViewMode::View2D {
                self.view_mode = ViewMode::Split2D3D;
            }
        }
    }

    fn draw_2d_canvas_at_rect(&mut self, ui: &mut egui::Ui, rect: egui::Rect) {
        let resp = ui.interact(
            rect,
            ui.id().with("chemcraft_2d_canvas"),
            Sense::click_and_drag(),
        );
        let painter = ui.painter_at(rect);
        let bg = if self.dark_mode {
            Color32::from_rgb(22, 24, 30)
        } else {
            Color32::from_rgb(255, 255, 255)
        };
        painter.rect_filled(rect, 4.0, bg);
        painter.rect_stroke(
            rect,
            4.0,
            Stroke::new(1.0_f32, Color32::from_gray(165)),
            egui::StrokeKind::Inside,
        );

        let style = SceneStyle {
            width: f64::from(rect.width()),
            height: f64::from(rect.height()),
            padding: 30.0,
            stroke_width: 2.0,
            font_size: 14.0,
            dark_mode: self.dark_mode,
            show_stereo_labels: true,
        };
        let mut scene = build_mol_scene(
            &self.session.canvas_mol,
            &style,
            &self.session.reaction_arrows,
        );
        scene.background = None;
        paint_vector_scene(&painter, &scene, rect.min);

        self.hover_atom = None;
        let mut hover_bond = None;
        if let Some(hp) = resp.hover_pos() {
            let lx = f64::from(hp.x - rect.min.x);
            let ly = f64::from(hp.y - rect.min.y);
            let mut best_d = 16.0_f64;
            for (aid, sp) in scene.atom_screen_pos.iter().enumerate() {
                let d = (sp[0] - lx).hypot(sp[1] - ly);
                if d < best_d {
                    best_d = d;
                    self.hover_atom = Some(aid);
                }
            }
            if self.hover_atom.is_none() {
                let mut best_bd = 11.0_f64;
                for (bid, b) in self.session.canvas_mol.bonds.iter().enumerate() {
                    if let (Some(p1), Some(p2)) = (
                        scene.atom_screen_pos.get(b.u),
                        scene.atom_screen_pos.get(b.v),
                    ) {
                        let d = point_seg_dist([lx, ly], *p1, *p2);
                        if d < best_bd {
                            best_bd = d;
                            hover_bond = Some(bid);
                        }
                    }
                }
            }
        }

        if let Some(aid) = self.hover_atom
            && let Some(sp) = scene.atom_screen_pos.get(aid)
        {
            painter.circle_stroke(
                Pos2::new(rect.min.x + sp[0] as f32, rect.min.y + sp[1] as f32),
                11.0,
                Stroke::new(1.8_f32, Color32::from_rgb(34, 139, 230)),
            );
        } else if let Some(bid) = hover_bond
            && let (Some(p1), Some(p2)) = (
                scene
                    .atom_screen_pos
                    .get(self.session.canvas_mol.bonds[bid].u),
                scene
                    .atom_screen_pos
                    .get(self.session.canvas_mol.bonds[bid].v),
            )
        {
            painter.line_segment(
                [
                    Pos2::new(rect.min.x + p1[0] as f32, rect.min.y + p1[1] as f32),
                    Pos2::new(rect.min.x + p2[0] as f32, rect.min.y + p2[1] as f32),
                ],
                Stroke::new(4.5_f32, Color32::from_rgba_unmultiplied(34, 139, 230, 110)),
            );
        }

        if self.tool == DrawTool::SelectMove {
            if resp.drag_started() {
                self.drag_start_atom = self.hover_atom;
            }
            if resp.dragged()
                && let Some(aid) = self.drag_start_atom
            {
                let delta = resp.drag_delta();
                let scale = scene.world_scale.max(1e-6);
                let mut mol = (*self.session.canvas_mol).clone();
                if let Some(Some(coords)) = mol.atoms.get_mut(aid).map(|a| &mut a.coords_2d) {
                    coords[0] += f64::from(delta.x) / scale;
                    coords[1] -= f64::from(delta.y) / scale;
                }
                self.session.canvas_mol = std::sync::Arc::new(mol);
            }
            if resp.drag_stopped() {
                self.drag_start_atom = None;
            }
        }

        if resp.clicked()
            && let Some(cp) = resp.interact_pointer_pos()
        {
            let [wx, wy] =
                scene.screen_to_world(f64::from(cp.x - rect.min.x), f64::from(cp.y - rect.min.y));
            match self.tool {
                DrawTool::Ring(r) => {
                    self.run_cmd(Command::StampRing {
                        template_name: r.label().to_string(),
                        x: wx,
                        y: wy,
                        snap_atom: self.hover_atom,
                        snap_bond: hover_bond,
                    });
                }
                DrawTool::BondOrder(order) => {
                    let new_pos = [wx, wy];
                    self.run_cmd(Command::AddAtomAndBond {
                        from_atom: self.hover_atom.or_else(|| {
                            if self.session.canvas_mol.atoms.is_empty() {
                                None
                            } else {
                                Some(self.session.canvas_mol.atoms.len() - 1)
                            }
                        }),
                        symbol: self.active_element.to_string(),
                        x: new_pos[0],
                        y: new_pos[1],
                        order,
                        stereo: BondStereo::None,
                    });
                }
                DrawTool::BondStereo(stereo) => {
                    let mut mol = (*self.session.canvas_mol).clone();
                    if let Some(bid) = hover_bond {
                        mol.bonds[bid].stereo = stereo;
                        let _ = mol.perceive();
                        self.session.canvas_mol = std::sync::Arc::new(mol);
                        let _ = self.session.refresh_live_admet();
                    } else if let Some(aid) = self.hover_atom {
                        self.run_cmd(Command::AddAtomAndBond {
                            from_atom: Some(aid),
                            symbol: self.active_element.to_string(),
                            x: wx,
                            y: wy,
                            order: BondOrder::Single,
                            stereo,
                        });
                    }
                }
                DrawTool::Eraser => {
                    let mut mol = (*self.session.canvas_mol).clone();
                    if let Some(aid) = self.hover_atom {
                        let _ = mol.remove_atom(aid);
                        self.session.canvas_mol = std::sync::Arc::new(mol);
                        let _ = self.session.refresh_live_admet();
                        self.smiles_bar = self.session.canvas_mol.canonical_smiles();
                    } else if let Some(bid) = hover_bond {
                        let _ = mol.remove_bond(bid);
                        self.session.canvas_mol = std::sync::Arc::new(mol);
                        let _ = self.session.refresh_live_admet();
                        self.smiles_bar = self.session.canvas_mol.canonical_smiles();
                    }
                }
                DrawTool::ReactionArrow => {
                    let sx = f64::from(cp.x - rect.min.x);
                    let sy = f64::from(cp.y - rect.min.y);
                    if let Some(start) = self.arrow_start.take() {
                        let mut arrows = (*self.session.reaction_arrows).clone();
                        arrows.push(ReactionArrow {
                            start,
                            end: [sx, sy],
                            label_above: "Reagent".to_string(),
                            label_below: String::new(),
                        });
                        self.session.reaction_arrows = std::sync::Arc::new(arrows);
                    } else {
                        self.arrow_start = Some([sx, sy]);
                    }
                }
                DrawTool::SelectMove => {}
            }
        }
        self.last_canvas_scene = Some(scene);
    }

    fn save_project(&mut self) {
        let root = PathBuf::from(&self.export_dir);
        let sb = match PathSandbox::new(None, Some(&root)) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Save root error: {e}");
                return;
            }
        };
        let target = match sb.validate_write_path(std::path::Path::new("workspace.chemcraft")) {
            Ok(p) => p,
            Err(e) => {
                self.status_msg = format!("Save sandbox error: {e}");
                return;
            }
        };
        let doc = ChemCraftDocument {
            format: "chemcraft".to_string(),
            version: env!("CARGO_PKG_VERSION").to_string(),
            canvas_mol: (*self.session.canvas_mol).clone(),
            reaction_arrows: (*self.session.reaction_arrows).clone(),
            sar_table: (*self.session.sar_table).clone(),
            sar_filter: self.session.sar_filter.clone(),
        };
        match save_chemcraft_atomic(&target, &doc) {
            Ok(()) => self.status_msg = format!("Saved project -> {}", target.display()),
            Err(e) => self.status_msg = format!("Save error: {e}"),
        }
    }

    fn import_file(&mut self, path_str: &str) {
        let p = PathBuf::from(path_str.trim());
        let root = p.parent().unwrap_or_else(|| std::path::Path::new("."));
        let sb = match PathSandbox::new(Some(root), None) {
            Ok(s) => s,
            Err(e) => {
                self.status_msg = format!("Import root error: {e}");
                return;
            }
        };
        let rel = p
            .file_name()
            .map(std::path::Path::new)
            .unwrap_or(p.as_path());
        let ext = p
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();

        if ext == "chemcraft" {
            match sb.read_bytes_limited(rel, MAX_READ_BYTES).and_then(|b| {
                decode_chemcraft_zip(&b).map_err(|e| {
                    mol_engine::SandboxError::Io(p.display().to_string(), std::io::Error::other(e))
                })
            }) {
                Ok(doc) => {
                    self.session.canvas_mol = std::sync::Arc::new(doc.canvas_mol);
                    self.session.reaction_arrows = std::sync::Arc::new(doc.reaction_arrows);
                    self.session.sar_table = std::sync::Arc::new(doc.sar_table);
                    self.session.sar_filter = doc.sar_filter;
                    self.smiles_bar = self.session.canvas_mol.canonical_smiles();
                    let _ = self.session.refresh_live_admet();
                    self.status_msg = format!("Loaded {}", p.display());
                }
                Err(e) => self.status_msg = format!("Load .chemcraft error: {e}"),
            }
            return;
        }

        let (_, text) = match sb.read_to_string_limited(rel, MAX_READ_BYTES) {
            Ok(r) => r,
            Err(e) => {
                self.status_msg = format!("Import sandbox error: {e}");
                return;
            }
        };
        let mut table = mol_sar::SarTable::new();
        let res = if ext == "dwar" {
            parse_dwar(&text).map(|recs| {
                for r in recs.into_iter().flatten() {
                    table.add_molecule(r.mol, None, None);
                }
            })
        } else if ext == "csv" {
            parse_csv_records(&text).map(|recs| {
                for r in recs.into_iter().flatten() {
                    table.add_molecule(r.mol, None, None);
                }
            })
        } else {
            for r in parse_sdf(&text).into_iter().flatten() {
                table.add_molecule(r, None, None);
            }
            Ok(())
        };
        match res {
            Ok(()) if !table.rows.is_empty() => {
                self.session.sar_table = std::sync::Arc::new(table);
                self.session.recompute_sar_analytics();
                self.status_msg = format!(
                    "Imported {} compounds from {}",
                    self.session.sar_table.rows.len(),
                    p.display()
                );
            }
            Ok(()) => self.status_msg = "No valid molecules found in file".to_string(),
            Err(e) => self.status_msg = format!("Parse error: {e}"),
        }
    }

    fn handle_shortcuts(&mut self, ctx: &egui::Context) {
        let wants_kb = ctx.wants_keyboard_input();
        let mut toggle_palette = false;
        let mut do_undo = false;
        let mut do_redo = false;
        let mut do_clean = false;
        let mut do_demo = false;
        let mut atom_hotkey: Option<&'static str> = None;

        ctx.input(|i| {
            if i.modifiers.ctrl && i.key_pressed(egui::Key::K) {
                toggle_palette = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::Z) {
                do_undo = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::Y) {
                do_redo = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::L) {
                do_clean = true;
            }
            if !wants_kb && i.modifiers.ctrl && i.key_pressed(egui::Key::E) {
                do_demo = true;
            }
            if !wants_kb && !i.modifiers.ctrl {
                let keys: [(egui::Key, &'static str); 8] = [
                    (egui::Key::C, "C"),
                    (egui::Key::N, "N"),
                    (egui::Key::O, "O"),
                    (egui::Key::S, "S"),
                    (egui::Key::P, "P"),
                    (egui::Key::F, "F"),
                    (egui::Key::B, "Br"),
                    (egui::Key::I, "I"),
                ];
                for (k, sym) in keys {
                    if i.key_pressed(k) {
                        atom_hotkey = Some(sym);
                    }
                }
            }
        });

        if toggle_palette {
            self.palette_open = !self.palette_open;
            self.palette_query.clear();
        }
        if do_undo {
            self.run_cmd(Command::Undo);
        }
        if do_redo {
            self.run_cmd(Command::Redo);
        }
        if do_clean {
            self.run_cmd(Command::CleanStructure2D);
        }
        if do_demo {
            self.run_cmd(Command::LoadDemo50Sar);
        }
        if let Some(sym) = atom_hotkey {
            self.active_element = sym;
            if let Some(aid) = self.hover_atom {
                self.run_cmd(Command::SetAtomElement {
                    atom_idx: aid,
                    symbol: sym.to_string(),
                });
            }
        }
    }

    fn render_palette(&mut self, ctx: &egui::Context) {
        if !self.palette_open {
            return;
        }
        let mut close = false;
        let mut chosen = None;
        egui::Window::new("ChemCraft Command Palette (Ctrl+K)")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_TOP, [0.0, 64.0])
            .fixed_size([440.0, 290.0])
            .show(ctx, |ui| {
                let resp = ui.add(
                    egui::TextEdit::singleline(&mut self.palette_query)
                        .hint_text("Type a command..."),
                );
                resp.request_focus();
                let q = self.palette_query.trim().to_ascii_lowercase();
                ui.separator();
                egui::ScrollArea::vertical().show(ui, |ui| {
                    for info in CommandRegistry::catalog() {
                        if !q.is_empty() && !info.label.to_ascii_lowercase().contains(&q) {
                            continue;
                        }
                        let lbl = match info.shortcut {
                            Some(sc) => format!("{}   [{sc}]", info.label),
                            None => info.label.to_string(),
                        };
                        if ui.button(lbl).clicked() {
                            chosen = Some(info.id);
                            close = true;
                        }
                    }
                });
                if ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                    close = true;
                }
            });
        if let Some(id) = chosen {
            match id {
                "clean_2d" => self.run_cmd(Command::CleanStructure2D),
                "load_demo_50" => self.run_cmd(Command::LoadDemo50Sar),
                "swap_tetrazole" => self.run_cmd(Command::MutateBioisostere {
                    atom_index: None,
                    group: "tetrazole".to_string(),
                }),
                "add_to_sar" => self.run_cmd(Command::AddCanvasToSar { pic50: None }),
                "filter_canvas" => self.run_cmd(Command::FilterSarByCanvasSubstructure),
                "clear_filter" => self.run_cmd(Command::ClearSarFilter),
                "clear_canvas" => self.run_cmd(Command::ClearCanvas),
                "undo" => self.run_cmd(Command::Undo),
                "redo" => self.run_cmd(Command::Redo),
                _ => {}
            }
        }
        if close {
            self.palette_open = false;
        }
    }
}

impl eframe::App for ChemCraftApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.dark_mode {
            ctx.set_visuals(egui::Visuals::dark());
        } else {
            ctx.set_visuals(egui::Visuals::light());
        }
        self.handle_shortcuts(ctx);
        self.render_palette(ctx);
        self.render_peptide_drawer(ctx);

        // 1. Top Toolbar
        egui::TopBottomPanel::top("top_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading(RichText::new("ChemCraft").strong());
                ui.label(RichText::new("2D/3D + Peptides + SAR + ChimeraX + AdmetCraft").weak());
                ui.separator();
                if ui.button("Demo 50 SAR").clicked() {
                    self.run_cmd(Command::LoadDemo50Sar);
                }
                if ui.button("Clean 2D").clicked() {
                    self.run_cmd(Command::CleanStructure2D);
                }
                if ui.button("Peptide Builder").clicked() {
                    self.peptide_drawer_open = !self.peptide_drawer_open;
                }
                if ui.button("-COOH -> Tetrazole").clicked() {
                    self.run_cmd(Command::MutateBioisostere {
                        atom_index: None,
                        group: "tetrazole".to_string(),
                    });
                }
                if ui
                    .add_enabled(self.session.can_undo(), egui::Button::new("Undo"))
                    .clicked()
                {
                    self.run_cmd(Command::Undo);
                }
                if ui
                    .add_enabled(self.session.can_redo(), egui::Button::new("Redo"))
                    .clicked()
                {
                    self.run_cmd(Command::Redo);
                }
                ui.separator();
                if ui.button("Save .chemcraft").clicked() {
                    self.save_project();
                }
                if ui.button("SVG").clicked() {
                    self.export_canvas("svg");
                }
                if ui.button("PDF").clicked() {
                    self.export_canvas("pdf");
                }
                if ui.button("EPS").clicked() {
                    self.export_canvas("eps");
                }
                if ui.button("PNG").clicked() {
                    self.export_canvas("png");
                }
                if ui.button("PDB").clicked() {
                    self.export_structure_3d("pdb");
                }
                if ui.button("3D SDF").clicked() {
                    self.export_structure_3d("sdf3d");
                }
                ui.separator();
                let cx_btn = ui.button("Visualizar en ChimeraX").on_hover_text(if self.chimerax_available {
                    "Exporta .pdb/.sdf + .cxc en la carpeta sandbox de ChimeraX y lo abre en 3D avanzado"
                } else {
                    "Exporta bundle .pdb/.sdf + .cxc para abrir en UCSF ChimeraX"
                });
                if cx_btn.clicked() {
                    self.open_in_chimerax();
                }
                if ui
                    .button(if self.dark_mode { "Light" } else { "Dark" })
                    .clicked()
                {
                    self.dark_mode = !self.dark_mode;
                }
                if ui.button("Palette (Ctrl+K)").clicked() {
                    self.palette_open = true;
                }
            });
        });

        // Bottom Status Bar
        egui::TopBottomPanel::bottom("status_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(&self.status_msg);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(format!("Export dir: {}", self.export_dir));
                });
            });
        });

        // 2. Left Panel: Vector Drawing Tools & Ring Templates
        egui::SidePanel::left("tools_panel")
            .default_width(185.0)
            .show(ctx, |ui| {
                ui.strong("Drawing Tools");
                ui.separator();
                let bond_tools = [
                    ("Select / Move", DrawTool::SelectMove),
                    ("Single (-)", DrawTool::BondOrder(BondOrder::Single)),
                    ("Double (=)", DrawTool::BondOrder(BondOrder::Double)),
                    ("Triple (#)", DrawTool::BondOrder(BondOrder::Triple)),
                    ("Aromatic (:)", DrawTool::BondOrder(BondOrder::Aromatic)),
                    ("Dative (->)", DrawTool::BondOrder(BondOrder::Dative)),
                    ("Up Wedge", DrawTool::BondStereo(BondStereo::Up)),
                    ("Down Hash", DrawTool::BondStereo(BondStereo::Down)),
                    ("Wavy Stereo", DrawTool::BondStereo(BondStereo::Wavy)),
                    ("Reaction Arrow", DrawTool::ReactionArrow),
                    ("Eraser", DrawTool::Eraser),
                ];
                for (lbl, t) in bond_tools {
                    if ui.selectable_label(self.tool == t, lbl).clicked() {
                        self.tool = t;
                    }
                }

                ui.add_space(6.0);
                ui.strong("Ring Templates (Snap)");
                ui.separator();
                let rings = [
                    RingTemplate::Benzene,
                    RingTemplate::Cyclohexane,
                    RingTemplate::Cyclopentane,
                    RingTemplate::Pyridine,
                    RingTemplate::Indole,
                    RingTemplate::Naphthalene,
                    RingTemplate::Tetrazole,
                    RingTemplate::Cyclobutane,
                    RingTemplate::Cyclopropane,
                ];
                for r in rings {
                    if ui
                        .selectable_label(self.tool == DrawTool::Ring(r), r.label())
                        .clicked()
                    {
                        self.tool = DrawTool::Ring(r);
                    }
                }

                ui.add_space(6.0);
                ui.strong("Element Hotkeys");
                ui.separator();
                ui.horizontal_wrapped(|ui| {
                    for sym in ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I"] {
                        if ui
                            .selectable_label(self.active_element == sym, sym)
                            .clicked()
                        {
                            self.active_element = sym;
                        }
                    }
                });
                ui.label(
                    RichText::new("Hover atom + key changes element")
                        .small()
                        .weak(),
                );

                ui.add_space(6.0);
                ui.strong("Quick Peptide (2D + 3D)");
                ui.separator();
                ui.add(
                    egui::TextEdit::singleline(&mut self.peptide_seq)
                        .hint_text("Ala-Gly-Tyr-Ser / AGYS"),
                );
                ui.horizontal(|ui| {
                    if ui.button("Build Peptide").clicked() {
                        let seq = self.peptide_seq.clone();
                        let ss = self.peptide_ss;
                        self.run_cmd(Command::BuildPeptide {
                            sequence: seq,
                            secondary_structure: ss,
                        });
                        self.viewer_3d.reset();
                        if self.view_mode == ViewMode::View2D {
                            self.view_mode = ViewMode::Split2D3D;
                        }
                    }
                    if ui.button("20-AA Pad...").clicked() {
                        self.peptide_drawer_open = true;
                    }
                });

                ui.add_space(8.0);
                ui.strong("Import (.sdf/.dwar/.csv/.chemcraft)");
                ui.add(
                    egui::TextEdit::singleline(&mut self.file_path_bar)
                        .hint_text("/path/to/file.sdf"),
                );
                if ui.button("Import File").clicked() && !self.file_path_bar.trim().is_empty() {
                    let p = self.file_path_bar.clone();
                    self.import_file(&p);
                }
            });

        // 3. Right Dockable Panel: Live ADMET Inspector (ChemCraft <-> AdmetCraft Bridge)
        egui::SidePanel::right("admet_inspector")
            .default_width(305.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Live ADMET Inspector");
                    ui.label(
                        RichText::new(
                            if self.session.prefer_external_admet && self.admet_cli_available {
                                "Engine: External AdmetCraft CLI (`admetcraft-cli`)"
                            } else {
                                "Engine: ChemCraft Quick ADMET (Built-in)"
                            },
                        )
                        .small()
                        .weak(),
                    );
                    if self.admet_cli_available {
                        let mut ext = self.session.prefer_external_admet;
                        if ui
                            .checkbox(&mut ext, "Delegate evals to external `admetcraft-cli`")
                            .changed()
                        {
                            self.session.prefer_external_admet = ext;
                            let _ = self.session.refresh_live_admet();
                        }
                    }
                    ui.separator();

                    if let Some(ref eval) = self.session.live_admet {
                        draw_admet_radar(ui, eval, self.dark_mode);
                        ui.separator();
                        let p = &eval.physicochemical;
                        let l = &eval.lipophilicity;
                        let pk = &eval.pharmacokinetics;
                        let d = &eval.druglikeness;
                        let mc = &eval.medicinal_chemistry;

                        egui::Grid::new("admet_kv").striped(true).show(ui, |ui| {
                            ui.label("Formula");
                            ui.monospace(&p.formula);
                            ui.end_row();
                            ui.label("MW");
                            ui.label(format!("{:.2} g/mol", p.molecular_weight));
                            ui.end_row();
                            ui.label("Consensus LogP");
                            ui.label(format!("{:.2}", l.consensus_log_p));
                            ui.end_row();
                            ui.label("TPSA");
                            ui.label(format!("{:.2} A^2", p.tpsa));
                            ui.end_row();
                            ui.label("HBD / HBA");
                            ui.label(format!("{} / {}", p.h_bond_donors, p.h_bond_acceptors));
                            ui.end_row();
                            ui.label("ESOL LogS");
                            ui.label(format!(
                                "{:.2} ({})",
                                eval.water_solubility.esol.log_s,
                                eval.water_solubility.esol.class.as_str()
                            ));
                            ui.end_row();
                            ui.label("GI Absorption");
                            ui.strong(pk.gi_absorption.as_str());
                            ui.end_row();
                            ui.label("BBB Permeant");
                            ui.label(if pk.bbb_permeant { "Yes" } else { "No" });
                            ui.end_row();
                            ui.label("Lipinski");
                            ui.label(format!("{} violation(s)", d.lipinski.violations));
                            ui.end_row();
                            ui.label("PAINS / Brenk");
                            ui.label(format!(
                                "{} / {} alert(s)",
                                mc.pains_alerts.len(),
                                mc.brenk_alerts.len()
                            ));
                            ui.end_row();
                        });

                        if let Some(ref delta) = self.session.last_mutation_delta {
                            ui.add_space(6.0);
                            ui.group(|ui| {
                                ui.strong("Bioisostere ADMET Delta");
                                ui.label(format!(
                                    "ΔMW: {:+.2} | ΔLogP: {:+.2} | ΔTPSA: {:+.2}",
                                    delta.delta_mw, delta.delta_logp, delta.delta_tpsa
                                ));
                                ui.label(format!("GI: {} -> {}", delta.old_gi, delta.new_gi));
                            });
                        }
                    } else {
                        ui.label("Draw a structure or build a peptide to inspect properties.");
                    }

                    ui.add_space(10.0);
                    ui.separator();
                    ui.strong("Handoff to AdmetCraft (Dedicated PK Studio)");
                    ui.horizontal_wrapped(|ui| {
                        if ui.button("Export ADMET CSV").clicked() {
                            self.export_admet_report(AdmetReportFormat::Csv, "chemcraft_admet.csv");
                        }
                        if ui.button("Export ADMET PDF").clicked() {
                            self.export_admet_report(AdmetReportFormat::Pdf, "chemcraft_admet.pdf");
                        }
                        if ui
                            .add_enabled(
                                self.admet_cli_available,
                                egui::Button::new("BOILED-Egg SVG (AdmetCraft)"),
                            )
                            .clicked()
                        {
                            self.export_admet_boiled_egg();
                        }
                        if ui
                            .add_enabled(
                                self.admet_gui_available,
                                egui::Button::new("Open in AdmetCraft GUI"),
                            )
                            .clicked()
                        {
                            self.open_in_admetcraft();
                        }
                    });
                });
            });

        // 4. Central Panel: Top 2D Drawing Canvas + Bottom SAR Explorer
        egui::CentralPanel::default().show(ctx, |ui| {
            // Top controls row for canvas
            ui.horizontal(|ui| {
                ui.label("SMILES:");
                let smi_edit = ui.add(
                    egui::TextEdit::singleline(&mut self.smiles_bar)
                        .desired_width(340.0)
                        .font(egui::TextStyle::Monospace),
                );
                if (smi_edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)))
                    || ui.button("Load").clicked()
                {
                    let smi = self.smiles_bar.clone();
                    self.run_cmd(Command::LoadCanvasSmiles {
                        smiles: smi,
                        name: None,
                    });
                }
                if ui.button("Add to SAR Table").clicked() {
                    self.run_cmd(Command::AddCanvasToSar { pic50: None });
                }
                if ui.button("Filter SAR by Fragment").clicked() {
                    self.run_cmd(Command::FilterSarByCanvasSubstructure);
                }
                if ui.button("Clear Canvas").clicked() {
                    self.run_cmd(Command::ClearCanvas);
                }
            });

            ui.horizontal_wrapped(|ui| {
                ui.strong("View:");
                ui.selectable_value(&mut self.view_mode, ViewMode::View2D, "2D Vector");
                ui.selectable_value(&mut self.view_mode, ViewMode::View3D, "3D Interactive");
                ui.selectable_value(&mut self.view_mode, ViewMode::Split2D3D, "Split 2D + 3D");
                if self.view_mode != ViewMode::View2D {
                    ui.separator();
                    ui.checkbox(&mut self.viewer_3d.show_ribbon, "Peptide Ribbon");
                    ui.checkbox(&mut self.viewer_3d.show_labels, "Atom Labels");
                    ui.checkbox(&mut self.viewer_3d.auto_spin, "Spin");
                    if ui.button("Rebuild 3D").clicked() {
                        let mut mol = (*self.session.canvas_mol).clone();
                        generate_3d_coords(&mut mol);
                        self.session.canvas_mol = std::sync::Arc::new(mol);
                        self.status_msg = "Regenerated 3D conformer coordinates".to_string();
                    }
                    if ui.button("Reset Cam").clicked() {
                        self.viewer_3d.reset();
                    }
                    if ui.button("Visualizar en ChimeraX").clicked() {
                        self.open_in_chimerax();
                    }
                }
            });

            if !self.session.canvas_mol.warnings.is_empty() {
                ui.colored_label(
                    Color32::from_rgb(220, 48, 48),
                    format!("Warnings: {}", self.session.canvas_mol.warnings.join(" | ")),
                );
            }

            // Center-top interactive viewport (2D Vector / 3D Conformer / Split 2D + 3D)
            let canvas_h = (ui.available_height() * 0.44).max(220.0);
            let canvas_w = ui.available_width().max(320.0);
            let (full_rect, _) =
                ui.allocate_exact_size(Vec2::new(canvas_w, canvas_h), Sense::hover());

            match self.view_mode {
                ViewMode::View2D => {
                    self.draw_2d_canvas_at_rect(ui, full_rect);
                }
                ViewMode::View3D => {
                    let mol_ref = std::sync::Arc::clone(&self.session.canvas_mol);
                    render_3d_viewport(
                        ui,
                        &mol_ref,
                        &mut self.viewer_3d,
                        full_rect,
                        self.dark_mode,
                    );
                }
                ViewMode::Split2D3D => {
                    let gap = 6.0_f32;
                    let half_w = ((full_rect.width() - gap) * 0.5).max(120.0);
                    let left_rect = egui::Rect::from_min_size(
                        full_rect.min,
                        Vec2::new(half_w, full_rect.height()),
                    );
                    let right_rect = egui::Rect::from_min_size(
                        Pos2::new(full_rect.min.x + half_w + gap, full_rect.min.y),
                        Vec2::new(half_w, full_rect.height()),
                    );
                    self.draw_2d_canvas_at_rect(ui, left_rect);
                    let mol_ref = std::sync::Arc::clone(&self.session.canvas_mol);
                    render_3d_viewport(
                        ui,
                        &mol_ref,
                        &mut self.viewer_3d,
                        right_rect,
                        self.dark_mode,
                    );
                }
            }

            ui.separator();
            // 5. Center-Bottom Panel: Interactive DataWarrior-Style SAR Explorer
            let matching = self
                .session
                .sar_table
                .matching_indices(&self.session.sar_filter)
                .unwrap_or_else(|_| (0..self.session.sar_table.rows.len()).collect());

            ui.horizontal(|ui| {
                ui.strong(format!(
                    "SAR Explorer ({}/{} shown)",
                    matching.len(),
                    self.session.sar_table.rows.len()
                ));
                ui.separator();
                ui.selectable_value(&mut self.sar_tab, SarTab::Table, "SAR Table");
                ui.selectable_value(&mut self.sar_tab, SarTab::Scatter, "Scatter Plot");
                ui.selectable_value(
                    &mut self.sar_tab,
                    SarTab::MmpCliffs,
                    "MMP & Activity Cliffs",
                );
                ui.selectable_value(&mut self.sar_tab, SarTab::Filters, "Dynamic Filters");
            });
            ui.separator();

            match self.sar_tab {
                SarTab::Table => {
                    egui::ScrollArea::both().show(ui, |ui| {
                        let mut select_row = None;
                        egui::Grid::new("sar_grid").striped(true).show(ui, |ui| {
                            for h in [
                                "2D Structure",
                                "Compound",
                                "Series",
                                "pIC50",
                                "LE",
                                "MW",
                                "cLogP",
                                "TPSA",
                                "GI",
                                "PAINS",
                            ] {
                                ui.strong(h);
                            }
                            ui.end_row();

                            for &idx in matching.iter().take(100) {
                                let row = &self.session.sar_table.rows[idx];
                                let _ = paint_mol_thumbnail(
                                    ui,
                                    &row.mol,
                                    &[],
                                    Vec2::new(84.0, 52.0),
                                    self.dark_mode,
                                );
                                let is_sel = self.session.selected_sar_idx == Some(idx);
                                if ui.selectable_label(is_sel, &row.name).clicked() {
                                    select_row = Some(idx);
                                }
                                ui.label(&row.series);
                                ui.label(row.pic50.map_or("-".to_string(), |v| format!("{v:.2}")));
                                ui.label(
                                    row.ligand_efficiency()
                                        .map_or("-".to_string(), |v| format!("{v:.2}")),
                                );
                                ui.label(format!("{:.1}", row.mw));
                                if let Some(ref a) = row.admet {
                                    ui.label(format!("{:.2}", a.consensus_logp));
                                    ui.label(format!("{:.1}", a.tpsa));
                                    ui.label(&a.gi_absorption);
                                    ui.label(a.pains_count.to_string());
                                } else {
                                    ui.label("-");
                                    ui.label("-");
                                    ui.label("-");
                                    ui.label("-");
                                }
                                ui.end_row();
                            }
                        });
                        if let Some(idx) = select_row {
                            self.run_cmd(Command::SelectSarRow { idx });
                        }
                    });
                }
                SarTab::Scatter => {
                    ui.horizontal(|ui| {
                        ui.label("X Axis:");
                        for ax in [
                            ScatterXAxis::ConsensusLogP,
                            ScatterXAxis::Tpsa,
                            ScatterXAxis::MolecularWeight,
                            ScatterXAxis::LigandEfficiency,
                        ] {
                            ui.selectable_value(&mut self.scatter_x, ax, ax.label());
                        }
                    });
                    if let Some(clicked) = draw_sar_scatter(
                        ui,
                        &self.session.sar_table.rows,
                        &matching,
                        self.scatter_x,
                        self.session.selected_sar_idx,
                        self.dark_mode,
                    ) {
                        self.run_cmd(Command::SelectSarRow { idx: clicked });
                    }
                }
                SarTab::MmpCliffs => {
                    egui::ScrollArea::vertical().show(ui, |ui| {
                        ui.strong(format!(
                            "Activity Cliffs (SALI) — {} pairs",
                            self.session.cliffs.len()
                        ));
                        egui::Grid::new("cliffs_grid").striped(true).show(ui, |ui| {
                            ui.strong("Mol A");
                            ui.strong("Mol B");
                            ui.strong("Tanimoto");
                            ui.strong("|ΔpIC50|");
                            ui.strong("SALI");
                            ui.end_row();
                            for c in self.session.cliffs.iter().take(25) {
                                ui.label(&c.mol_a_name);
                                ui.label(&c.mol_b_name);
                                ui.label(format!("{:.2}", c.tanimoto));
                                ui.label(format!("{:.2}", c.delta_pic50));
                                ui.strong(format!("{:.1}", c.sali));
                                ui.end_row();
                            }
                        });
                        ui.add_space(8.0);
                        ui.strong(format!(
                            "Matched Molecular Pairs — {} transforms",
                            self.session.mmps.len()
                        ));
                        egui::Grid::new("mmp_grid").striped(true).show(ui, |ui| {
                            ui.strong("Mol A -> Mol B");
                            ui.strong("R1 -> R2");
                            ui.strong("ΔpIC50");
                            ui.strong("ΔcLogP");
                            ui.strong("ΔTPSA");
                            ui.end_row();
                            for m in self.session.mmps.iter().take(25) {
                                ui.label(format!("{} -> {}", m.mol_a_name, m.mol_b_name));
                                ui.monospace(format!("{} -> {}", m.r_group_a, m.r_group_b));
                                ui.label(
                                    m.delta_pic50
                                        .map_or("-".to_string(), |v| format!("{v:+.2}")),
                                );
                                ui.label(
                                    m.delta_clogp
                                        .map_or("-".to_string(), |v| format!("{v:+.2}")),
                                );
                                ui.label(
                                    m.delta_tpsa.map_or("-".to_string(), |v| format!("{v:+.1}")),
                                );
                                ui.end_row();
                            }
                        });
                    });
                }
                SarTab::Filters => {
                    ui.horizontal(|ui| {
                        ui.label("SMARTS Substructure:");
                        let mut sm = self
                            .session
                            .sar_filter
                            .smarts_substructure
                            .clone()
                            .unwrap_or_default();
                        if ui.text_edit_singleline(&mut sm).changed() {
                            self.session.sar_filter.smarts_substructure =
                                if sm.trim().is_empty() { None } else { Some(sm) };
                        }
                        if ui.button("Use Canvas").clicked() {
                            self.run_cmd(Command::FilterSarByCanvasSubstructure);
                        }
                        if ui.button("Clear All").clicked() {
                            self.run_cmd(Command::ClearSarFilter);
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("Tanimoto Reference SMILES:");
                        let mut rf = self
                            .session
                            .sar_filter
                            .reference_smiles
                            .clone()
                            .unwrap_or_else(|| self.session.canvas_mol.canonical_smiles());
                        if ui.text_edit_singleline(&mut rf).changed() {
                            self.session.sar_filter.reference_smiles = Some(rf);
                        }
                        ui.add(
                            egui::Slider::new(&mut self.session.sar_filter.min_tanimoto, 0.0..=1.0)
                                .text("Min Tanimoto"),
                        );
                    });
                    ui.horizontal(|ui| {
                        ui.checkbox(
                            &mut self.session.sar_filter.lipinski_pass_only,
                            "Lipinski 0 violations",
                        );
                        ui.checkbox(
                            &mut self.session.sar_filter.no_pains_only,
                            "No PAINS alerts",
                        );
                    });
                }
            }
        });
    }
}

fn point_seg_dist(p: [f64; 2], a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-6 {
        return (p[0] - a[0]).hypot(p[1] - a[1]);
    }
    let t = (((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2).clamp(0.0, 1.0);
    let px = a[0] + t * dx;
    let py = a[1] + t * dy;
    (p[0] - px).hypot(p[1] - py)
}

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1340.0, 860.0])
            .with_title("ChemCraft — Vector Molecular Design & SAR Studio"),
        ..Default::default()
    };
    eframe::run_native(
        "ChemCraft",
        options,
        Box::new(|_cc| Ok(Box::new(ChemCraftApp::new()))),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn app_initializes_with_50_sar_and_runs_bioisostere_and_undo() {
        let mut app = ChemCraftApp::new();
        assert_eq!(app.session.sar_table.rows.len(), 50);
        assert!(app.session.live_admet.is_some());
        app.run_cmd(Command::MutateBioisostere {
            atom_index: None,
            group: "tetrazole".to_string(),
        });
        assert!(app.session.last_mutation_delta.is_some());
        assert!(app.session.can_undo());
        app.run_cmd(Command::Undo);
        assert!(app.session.can_redo());
    }
}
