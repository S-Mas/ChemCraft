//! Pure-Rust interactive 3D molecular & peptide ball-and-stick + backbone ribbon viewer
//! rendered directly on `egui::Painter` (zero extra GPU/shader dependencies).

use eframe::egui;
use egui::{Color32, Pos2, Rect, Response, Sense, Stroke, Vec2};
use mol_graph::{BondOrder, MolGraph, compute_3d_coords};

#[derive(Debug, Clone)]
pub struct Viewer3dState {
    pub yaw: f32,
    pub pitch: f32,
    pub zoom: f32,
    pub show_ribbon: bool,
    pub show_labels: bool,
    pub auto_spin: bool,
}

impl Default for Viewer3dState {
    fn default() -> Self {
        Self {
            yaw: 0.45,
            pitch: 0.25,
            zoom: 1.0,
            show_ribbon: true,
            show_labels: true,
            auto_spin: false,
        }
    }
}

impl Viewer3dState {
    pub fn reset(&mut self) {
        self.yaw = 0.45;
        self.pitch = 0.25;
        self.zoom = 1.0;
    }
}

fn elem_color(symbol: &str, dark_mode: bool) -> Color32 {
    match symbol {
        "C" => {
            if dark_mode {
                Color32::from_rgb(155, 165, 178)
            } else {
                Color32::from_rgb(76, 82, 92)
            }
        }
        "N" => Color32::from_rgb(56, 128, 246),
        "O" => Color32::from_rgb(230, 62, 56),
        "S" => Color32::from_rgb(232, 188, 42),
        "P" => Color32::from_rgb(235, 132, 34),
        "F" => Color32::from_rgb(104, 214, 114),
        "Cl" => Color32::from_rgb(44, 196, 86),
        "Br" => Color32::from_rgb(168, 48, 42),
        "I" => Color32::from_rgb(148, 48, 178),
        "H" => Color32::from_rgb(225, 228, 232),
        "Fe" => Color32::from_rgb(208, 108, 52),
        _ => Color32::from_rgb(190, 120, 215),
    }
}

fn elem_radius(symbol: &str) -> f32 {
    match symbol {
        "H" => 0.26,
        "F" | "O" | "N" => 0.37,
        "C" => 0.40,
        "S" | "P" | "Cl" => 0.48,
        "Br" | "I" | "Fe" => 0.54,
        _ => 0.42,
    }
}

fn scale_rgb(c: Color32, f: f32) -> Color32 {
    let f = f.clamp(0.0, 2.2);
    Color32::from_rgb(
        ((c.r() as f32 * f).clamp(0.0, 255.0)) as u8,
        ((c.g() as f32 * f).clamp(0.0, 255.0)) as u8,
        ((c.b() as f32 * f).clamp(0.0, 255.0)) as u8,
    )
}

fn mix_white(c: Color32, t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    Color32::from_rgb(
        (c.r() as f32 * (1.0 - t) + 255.0 * t) as u8,
        (c.g() as f32 * (1.0 - t) + 255.0 * t) as u8,
        (c.b() as f32 * (1.0 - t) + 255.0 * t) as u8,
    )
}

fn rainbow_color(t: f32) -> Color32 {
    let t = t.clamp(0.0, 1.0);
    // Blue -> Cyan -> Green -> Yellow -> Red N-to-C peptide ramp
    let (r, g, b) = if t < 0.25 {
        let u = t / 0.25;
        (45.0, 110.0 + 125.0 * u, 240.0)
    } else if t < 0.50 {
        let u = (t - 0.25) / 0.25;
        (45.0, 235.0, 240.0 - 140.0 * u)
    } else if t < 0.75 {
        let u = (t - 0.50) / 0.25;
        (45.0 + 200.0 * u, 235.0 - 25.0 * u, 100.0 - 65.0 * u)
    } else {
        let u = (t - 0.75) / 0.25;
        (245.0, 210.0 - 155.0 * u, 35.0)
    };
    Color32::from_rgb(r as u8, g as u8, b as u8)
}

#[derive(Clone)]
enum DrawItem {
    RibbonSeg {
        z: f32,
        p1: Pos2,
        p2: Pos2,
        width: f32,
        color: Color32,
    },
    Bond {
        z: f32,
        p1: Pos2,
        p2: Pos2,
        order: BondOrder,
        width: f32,
        c1: Color32,
        c2: Color32,
    },
    Atom {
        z: f32,
        pos: Pos2,
        radius: f32,
        color: Color32,
        label: Option<String>,
        hovered: bool,
    },
}

impl DrawItem {
    fn depth(&self) -> f32 {
        match self {
            Self::RibbonSeg { z, .. } | Self::Bond { z, .. } | Self::Atom { z, .. } => *z,
        }
    }
}

/// Draw an interactive 3D viewport into `rect`, returning the `egui::Response` for drag/zoom handling.
pub fn render_3d_viewport(
    ui: &mut egui::Ui,
    mol: &MolGraph,
    state: &mut Viewer3dState,
    rect: Rect,
    dark_mode: bool,
) -> Response {
    let resp = ui.interact(
        rect,
        ui.id().with("chemcraft_3d_view"),
        Sense::click_and_drag(),
    );

    if resp.dragged() {
        let d = resp.drag_delta();
        state.yaw += d.x * 0.011;
        state.pitch = (state.pitch + d.y * 0.011).clamp(-1.52, 1.52);
    }
    if resp.double_clicked() {
        state.reset();
    }
    if resp.hovered() {
        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll.abs() > 0.1 {
            state.zoom = (state.zoom * (1.0 + scroll * 0.0022)).clamp(0.35, 4.5);
        }
    }
    if state.auto_spin {
        state.yaw += 0.012;
        ui.ctx().request_repaint();
    }

    let painter = ui.painter_at(rect);
    let bg = if dark_mode {
        Color32::from_rgb(16, 19, 26)
    } else {
        Color32::from_rgb(248, 250, 253)
    };
    painter.rect_filled(rect, 4.0, bg);
    painter.rect_stroke(
        rect,
        4.0,
        Stroke::new(
            1.0_f32,
            Color32::from_gray(if dark_mode { 75 } else { 170 }),
        ),
        egui::StrokeKind::Inside,
    );

    if mol.atoms.is_empty() {
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            "3D Viewer — Load SMILES, Draw, or Build a Peptide",
            egui::FontId::proportional(14.0),
            Color32::GRAY,
        );
        return resp;
    }

    let coords = compute_3d_coords(mol);
    // Compute centered radius for automatic framing
    let mut cent = [0.0_f64; 3];
    for p in &coords {
        cent[0] += p[0];
        cent[1] += p[1];
        cent[2] += p[2];
    }
    let inv_n = 1.0 / coords.len() as f64;
    cent[0] *= inv_n;
    cent[1] *= inv_n;
    cent[2] *= inv_n;

    let mut max_r = 2.0_f64;
    for p in &coords {
        let dx = p[0] - cent[0];
        let dy = p[1] - cent[1];
        let dz = p[2] - cent[2];
        max_r = max_r.max((dx * dx + dy * dy + dz * dz).sqrt());
    }

    let (sy, cy) = state.yaw.sin_cos();
    let (sp, cp) = state.pitch.sin_cos();
    let cam_dist = (max_r as f32 * 2.6 + 8.0) / state.zoom;
    let base_scale = rect.width().min(rect.height()) * 0.39 * state.zoom / (max_r as f32);
    let focal = cam_dist;

    let mut proj: Vec<(Pos2, f32, f32)> = Vec::with_capacity(coords.len());
    for p in &coords {
        let x0 = (p[0] - cent[0]) as f32;
        let y0 = (p[1] - cent[1]) as f32;
        let z0 = (p[2] - cent[2]) as f32;

        // Rotate around Y (yaw) then X (pitch)
        let x1 = x0 * cy + z0 * sy;
        let z1 = -x0 * sy + z0 * cy;
        let y2 = y0 * cp - z1 * sp;
        let z2 = y0 * sp + z1 * cp;

        let depth = (cam_dist - z2).max(1.5);
        let persp = focal / depth;
        let sx = rect.center().x + x1 * base_scale * persp;
        let sy_scr = rect.center().y - y2 * base_scale * persp;
        proj.push((Pos2::new(sx, sy_scr), z2, persp));
    }

    // Hover detection (closest frontmost projected atom)
    let mut hovered_atom: Option<(usize, f32)> = None;
    if let Some(hp) = resp.hover_pos() {
        for (i, (pos, z, persp)) in proj.iter().enumerate() {
            let r = (elem_radius(&mol.atoms[i].symbol) * base_scale * 0.42 * persp).max(7.0);
            if pos.distance(hp) <= r + 3.0 && hovered_atom.is_none_or(|(_, best_z)| *z > best_z) {
                hovered_atom = Some((i, *z));
            }
        }
    }

    let mut items: Vec<DrawItem> = Vec::with_capacity(mol.atoms.len() + mol.bonds.len() + 32);

    // Optional peptide C-alpha backbone spline ribbon
    if state.show_ribbon {
        let ca_indices: Vec<usize> = mol
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(idx, a)| {
                a.pdb_meta
                    .as_ref()
                    .filter(|m| !m.is_hetatm && m.atom_name == "CA")
                    .map(|_| idx)
            })
            .collect();
        if ca_indices.len() >= 2 {
            for w in 0..(ca_indices.len() - 1) {
                let a = ca_indices[w];
                let b = ca_indices[w + 1];
                let (p1, z1, per1) = proj[a];
                let (p2, z2, per2) = proj[b];
                let t = w as f32 / (ca_indices.len() - 1).max(1) as f32;
                let width = (5.2 * (per1 + per2) * 0.5).clamp(2.5, 11.0);
                items.push(DrawItem::RibbonSeg {
                    z: (z1 + z2) * 0.5 - 0.15,
                    p1,
                    p2,
                    width,
                    color: rainbow_color(t),
                });
            }
        }
    }

    for b in &mol.bonds {
        let (p1, z1, per1) = proj[b.u];
        let (p2, z2, per2) = proj[b.v];
        let avg_z = (z1 + z2) * 0.5;
        let depth_shade =
            (0.72 + (avg_z / (max_r as f32 * 2.0)).clamp(-0.28, 0.28)).clamp(0.42, 1.0);
        let c1 = scale_rgb(elem_color(&mol.atoms[b.u].symbol, dark_mode), depth_shade);
        let c2 = scale_rgb(elem_color(&mol.atoms[b.v].symbol, dark_mode), depth_shade);
        let w = (0.14 * base_scale * (per1 + per2) * 0.5).clamp(1.8, 6.5);
        items.push(DrawItem::Bond {
            z: avg_z,
            p1,
            p2,
            order: b.order,
            width: w,
            c1,
            c2,
        });
    }

    for (i, a) in mol.atoms.iter().enumerate() {
        let (pos, z, persp) = proj[i];
        let r = (elem_radius(&a.symbol) * base_scale * 0.42 * persp).clamp(4.0, 24.0);
        let depth_shade = (0.85 + (z / (max_r as f32 * 2.2)).clamp(-0.30, 0.25)).clamp(0.45, 1.05);
        let col = scale_rgb(elem_color(&a.symbol, dark_mode), depth_shade);
        let label = if state.show_labels && a.symbol != "C" {
            Some(a.symbol.clone())
        } else {
            None
        };
        items.push(DrawItem::Atom {
            z,
            pos,
            radius: r,
            color: col,
            label,
            hovered: hovered_atom.is_some_and(|(hi, _)| hi == i),
        });
    }

    // Back-to-front painter's sort (smaller z = further away; larger z = closer to camera)
    items.sort_by(|a, b| a.depth().total_cmp(&b.depth()));

    let bond_outline = if dark_mode {
        Color32::from_rgba_unmultiplied(245, 245, 255, 45)
    } else {
        Color32::from_rgba_unmultiplied(30, 35, 45, 55)
    };

    for item in items {
        match item {
            DrawItem::RibbonSeg {
                p1,
                p2,
                width,
                color,
                ..
            } => {
                painter.line_segment(
                    [p1, p2],
                    Stroke::new(width + 2.0, color.linear_multiply(0.35)),
                );
                painter.line_segment([p1, p2], Stroke::new(width, color));
            }
            DrawItem::Bond {
                p1,
                p2,
                order,
                width,
                c1,
                c2,
                ..
            } => {
                let d = p2 - p1;
                let len = d.length().max(1e-3);
                let perp = Vec2::new(-d.y / len, d.x / len);
                let mid = p1 + d * 0.5;
                match order {
                    BondOrder::Double | BondOrder::Aromatic => {
                        let off = perp * (width * 0.72).max(2.1);
                        painter.line_segment([p1 + off, mid + off], Stroke::new(width * 0.68, c1));
                        painter.line_segment([mid + off, p2 + off], Stroke::new(width * 0.68, c2));
                        painter.line_segment([p1 - off, mid - off], Stroke::new(width * 0.68, c1));
                        painter.line_segment([mid - off, p2 - off], Stroke::new(width * 0.68, c2));
                    }
                    BondOrder::Triple => {
                        let off = perp * (width * 0.95).max(2.6);
                        painter.line_segment([p1, mid], Stroke::new(width * 0.6, c1));
                        painter.line_segment([mid, p2], Stroke::new(width * 0.6, c2));
                        painter.line_segment([p1 + off, mid + off], Stroke::new(width * 0.45, c1));
                        painter.line_segment([mid + off, p2 + off], Stroke::new(width * 0.45, c2));
                        painter.line_segment([p1 - off, mid - off], Stroke::new(width * 0.45, c1));
                        painter.line_segment([mid - off, p2 - off], Stroke::new(width * 0.45, c2));
                    }
                    _ => {
                        painter.line_segment([p1, p2], Stroke::new(width + 1.2, bond_outline));
                        painter.line_segment([p1, mid], Stroke::new(width, c1));
                        painter.line_segment([mid, p2], Stroke::new(width, c2));
                    }
                }
            }
            DrawItem::Atom {
                pos,
                radius,
                color,
                label,
                hovered,
                ..
            } => {
                // Shaded sphere (ambient base + diffuse offset + specular glint)
                painter.circle_filled(pos, radius, scale_rgb(color, 0.72));
                painter.circle_filled(
                    pos + Vec2::new(-radius * 0.16, -radius * 0.16),
                    radius * 0.78,
                    color,
                );
                painter.circle_filled(
                    pos + Vec2::new(-radius * 0.28, -radius * 0.28),
                    radius * 0.40,
                    mix_white(color, 0.45),
                );
                painter.circle_filled(
                    pos + Vec2::new(-radius * 0.34, -radius * 0.34),
                    radius * 0.17,
                    mix_white(color, 0.85),
                );
                painter.circle_stroke(
                    pos,
                    radius,
                    Stroke::new(
                        if hovered { 2.2_f32 } else { 0.8_f32 },
                        if hovered {
                            Color32::from_rgb(34, 139, 230)
                        } else {
                            Color32::from_rgba_unmultiplied(15, 18, 25, 110)
                        },
                    ),
                );
                if let Some(lbl) = label
                    && radius >= 6.5
                {
                    painter.text(
                        pos,
                        egui::Align2::CENTER_CENTER,
                        lbl,
                        egui::FontId::proportional((radius * 0.95).clamp(8.0, 12.5)),
                        Color32::WHITE,
                    );
                }
            }
        }
    }

    // Corner HUD hints + hover atom inspector badge
    let hud_text = if let Some((hi, _)) = hovered_atom {
        let a = &mol.atoms[hi];
        if let Some(meta) = &a.pdb_meta {
            format!(
                "Atom #{}: {} ({} {} {})",
                hi + 1,
                a.symbol,
                meta.atom_name,
                meta.res_name,
                meta.res_seq
            )
        } else {
            format!(
                "Atom #{}: {}{}",
                hi + 1,
                a.symbol,
                a.cip_label.map_or(String::new(), |c| format!(" ({c})"))
            )
        }
    } else {
        "Drag: Orbit | Scroll: Zoom | Double-Click: Reset | 3D Conformer".to_string()
    };
    painter.text(
        rect.left_bottom() + Vec2::new(8.0, -6.0),
        egui::Align2::LEFT_BOTTOM,
        hud_text,
        egui::FontId::proportional(11.0),
        Color32::from_gray(if dark_mode { 170 } else { 100 }),
    );

    resp
}
