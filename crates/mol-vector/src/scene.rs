use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use mol_graph::{BondOrder, BondStereo, MolAtom, MolGraph};

use crate::layout::ensure_2d_coords;

/// Palette configuration for vector chemical rendering (light/dark publication themes).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SceneStyle {
    pub width: f64,
    pub height: f64,
    pub padding: f64,
    pub stroke_width: f64,
    pub font_size: f64,
    pub dark_mode: bool,
    pub show_stereo_labels: bool,
}

impl Default for SceneStyle {
    fn default() -> Self {
        Self {
            width: 480.0,
            height: 360.0,
            padding: 28.0,
            stroke_width: 1.85,
            font_size: 13.5,
            dark_mode: false,
            show_stereo_labels: true,
        }
    }
}

/// Optional chemical reaction arrow primitive (VectorCraft smooth Bézier / shaft).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionArrow {
    pub start: [f64; 2],
    pub end: [f64; 2],
    pub label_above: String,
    pub label_below: String,
}

/// RGBA color tuple `0..=255`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Rgba(pub u8, pub u8, pub u8, pub u8);

/// 2D bounding box in canvas points.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Rect2 {
    pub x: f64,
    pub y: f64,
    pub w: f64,
    pub h: f64,
}

impl Rect2 {
    pub fn overlaps(&self, other: &Self, margin: f64) -> bool {
        self.x - margin < other.x + other.w
            && self.x + self.w + margin > other.x
            && self.y - margin < other.y + other.h
            && self.y + self.h + margin > other.y
    }
}

/// Vector primitive element in `VectorScene`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum VectorElement {
    Line {
        x1: f64,
        y1: f64,
        x2: f64,
        y2: f64,
        width: f64,
        color: Rgba,
        dashed: bool,
    },
    Polyline {
        points: Vec<[f64; 2]>,
        width: f64,
        color: Rgba,
        closed: bool,
    },
    Polygon {
        points: Vec<[f64; 2]>,
        fill: Rgba,
    },
    Circle {
        cx: f64,
        cy: f64,
        r: f64,
        fill: Option<Rgba>,
        stroke: Option<(Rgba, f64)>,
    },
    CubicCurve {
        p0: [f64; 2],
        p1: [f64; 2],
        p2: [f64; 2],
        p3: [f64; 2],
        width: f64,
        color: Rgba,
    },
    Label {
        text: String,
        cx: f64,
        cy: f64,
        font_size: f64,
        color: Rgba,
        bbox: Rect2,
    },
}

/// Clean-room 2D vector IR produced from a `MolGraph` (inspired by VectorCraft).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VectorScene {
    pub width: f64,
    pub height: f64,
    pub background: Option<Rgba>,
    pub elements: Vec<VectorElement>,
    pub atom_screen_pos: Vec<[f64; 2]>,
    pub label_boxes: Vec<Rect2>,
}

impl VectorScene {
    /// Count any overlapping atom/stereo label bounding boxes (publication QA metric: ideally 0).
    pub fn label_overlap_count(&self) -> usize {
        let mut count = 0usize;
        for i in 0..self.label_boxes.len() {
            for j in (i + 1)..self.label_boxes.len() {
                if self.label_boxes[i].overlaps(&self.label_boxes[j], 0.5) {
                    count += 1;
                }
            }
        }
        count
    }
}

/// Build a publication-grade `VectorScene` from a `MolGraph`.
pub fn build_mol_scene(
    mol_in: &MolGraph,
    style: &SceneStyle,
    arrows: &[ReactionArrow],
) -> VectorScene {
    let mut mol = mol_in.clone();
    ensure_2d_coords(&mut mol);

    let bg = if style.dark_mode {
        Some(Rgba(24, 26, 32, 255))
    } else {
        Some(Rgba(255, 255, 255, 255))
    };
    let fg = if style.dark_mode {
        Rgba(232, 236, 242, 255)
    } else {
        Rgba(28, 30, 36, 255)
    };

    if mol.atoms.is_empty() && arrows.is_empty() {
        return VectorScene {
            width: style.width,
            height: style.height,
            background: bg,
            elements: Vec::new(),
            atom_screen_pos: Vec::new(),
            label_boxes: Vec::new(),
        };
    }

    // Compute world bounds
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for a in &mol.atoms {
        let [x, y] = a.coords_2d.unwrap_or([0.0, 0.0]);
        min_x = min_x.min(x);
        max_x = max_x.max(x);
        min_y = min_y.min(y);
        max_y = max_y.max(y);
    }
    if !min_x.is_finite() {
        min_x = 0.0;
        max_x = 1.0;
        min_y = 0.0;
        max_y = 1.0;
    }

    let span_x = (max_x - min_x).max(2.0);
    let span_y = (max_y - min_y).max(2.0);
    let inner_w = (style.width - 2.0 * style.padding).max(20.0);
    let inner_h = (style.height - 2.0 * style.padding).max(20.0);
    let scale = (inner_w / span_x).min(inner_h / span_y).min(38.0);
    let cx_world = (min_x + max_x) * 0.5;
    let cy_world = (min_y + max_y) * 0.5;
    let cx_scr = style.width * 0.5;
    let cy_scr = style.height * 0.5;

    let to_screen = |wx: f64, wy: f64| -> [f64; 2] {
        [
            cx_scr + (wx - cx_world) * scale,
            cy_scr - (wy - cy_world) * scale,
        ]
    };

    let atom_pos: Vec<[f64; 2]> = mol
        .atoms
        .iter()
        .map(|a| {
            let [wx, wy] = a.coords_2d.unwrap_or([0.0, 0.0]);
            to_screen(wx, wy)
        })
        .collect();

    // Scale font size slightly for compact/dense molecules so labels don't crowd!
    let eff_font = (style.font_size * (scale / 28.0).clamp(0.72, 1.15)).clamp(8.5, 16.0);

    // 1. Prepare visible atom label boxes + collision resolution offsets
    type PendingLabel = (String, Rgba, [f64; 2], Rect2);
    let mut visible_labels: Vec<Option<PendingLabel>> = vec![None; mol.atoms.len()];
    for (idx, a) in mol.atoms.iter().enumerate() {
        let Some(text) = format_atom_label(a, &mol, &atom_pos) else {
            continue;
        };
        let p = atom_pos[idx];
        let est_w = (text.len() as f64) * eff_font * 0.62 + 3.0;
        let est_h = eff_font + 3.0;
        let color = element_color(&a.symbol, style.dark_mode);
        let bbox = Rect2 {
            x: p[0] - est_w * 0.5,
            y: p[1] - est_h * 0.5,
            w: est_w,
            h: est_h,
        };
        visible_labels[idx] = Some((text, color, p, bbox));
    }

    // Check and gently separate colliding label boxes (only slightly so bonds stay accurate)
    for _pass in 0..3 {
        for i in 0..visible_labels.len() {
            for j in (i + 1)..visible_labels.len() {
                if let (Some((_, _, p1, b1)), Some((_, _, p2, b2))) =
                    (&visible_labels[i], &visible_labels[j])
                    && b1.overlaps(b2, 1.0)
                {
                    let dx = p2[0] - p1[0];
                    let dy = p2[1] - p1[1];
                    let d = dx.hypot(dy).max(0.5);
                    let shift = 2.0;
                    let ux = dx / d * shift;
                    let uy = dy / d * shift;
                    if let Some((_, _, p_i, b_i)) = &mut visible_labels[i] {
                        p_i[0] -= ux;
                        p_i[1] -= uy;
                        b_i.x -= ux;
                        b_i.y -= uy;
                    }
                    if let Some((_, _, p_j, b_j)) = &mut visible_labels[j] {
                        p_j[0] += ux;
                        p_j[1] += uy;
                        b_j.x += ux;
                        b_j.y += uy;
                    }
                }
            }
        }
    }

    let mut elements = Vec::new();
    let mut label_boxes = Vec::new();

    // Valence error highlight halos underneath bonds
    for (idx, a) in mol.atoms.iter().enumerate() {
        if a.valence_exceeded {
            let [cx, cy] = atom_pos[idx];
            elements.push(VectorElement::Circle {
                cx,
                cy,
                r: (eff_font * 0.85).max(10.0),
                fill: Some(Rgba(225, 45, 45, 55)),
                stroke: Some((Rgba(220, 38, 38, 220), 1.5)),
            });
        }
    }

    // Ring centroids on screen for ring-inward double/aromatic bond offsets
    let ring_info: Vec<(BTreeSet<usize>, [f64; 2])> = mol
        .rings
        .iter()
        .filter(|r| r.len() >= 3)
        .map(|r| {
            let set: BTreeSet<usize> = r.iter().copied().collect();
            let mut sx = 0.0;
            let mut sy = 0.0;
            for &vid in r {
                sx += atom_pos[vid][0];
                sy += atom_pos[vid][1];
            }
            let k = r.len() as f64;
            (set, [sx / k, sy / k])
        })
        .collect();

    // 2. Render bonds trimmed cleanly around atom labels
    let stroke = (style.stroke_width * (scale / 28.0).clamp(0.75, 1.2)).clamp(1.1, 2.5);
    for b in &mol.bonds {
        let p1_raw = atom_pos[b.u];
        let p2_raw = atom_pos[b.v];
        let (p1, p2) = trim_segment_to_labels(
            p1_raw,
            p2_raw,
            visible_labels[b.u].as_ref().map(|x| &x.3),
            visible_labels[b.v].as_ref().map(|x| &x.3),
        );
        let dx = p2[0] - p1[0];
        let dy = p2[1] - p1[1];
        let len = dx.hypot(dy);
        if len < 1e-3 {
            continue;
        }
        let ux = dx / len;
        let uy = dy / len;
        let mut nx = -uy;
        let mut ny = ux;

        // If edge belongs to a ring, point normal `(nx, ny)` toward the ring centroid!
        let mut in_ring_centroid = None;
        for (rset, rc) in &ring_info {
            if rset.contains(&b.u) && rset.contains(&b.v) {
                in_ring_centroid = Some(*rc);
                break;
            }
        }
        if let Some(rc) = in_ring_centroid {
            let mid = [(p1[0] + p2[0]) * 0.5, (p1[1] + p2[1]) * 0.5];
            if (rc[0] - mid[0]) * nx + (rc[1] - mid[1]) * ny < 0.0 {
                nx = -nx;
                ny = -ny;
            }
        }

        // Stereo wedge / hash / wavy takes precedence on Single bonds
        match b.stereo {
            BondStereo::Up => {
                let w = (scale * 0.14).clamp(2.5, 5.2);
                elements.push(VectorElement::Polygon {
                    points: vec![
                        p1,
                        [p2[0] + nx * w, p2[1] + ny * w],
                        [p2[0] - nx * w, p2[1] - ny * w],
                    ],
                    fill: fg,
                });
                continue;
            }
            BondStereo::Down => {
                let w_max = (scale * 0.15).clamp(2.6, 5.4);
                let steps = 7usize;
                for s in 1..=steps {
                    let t = s as f64 / steps as f64;
                    let w = w_max * t;
                    let cx = p1[0] + dx * t;
                    let cy = p1[1] + dy * t;
                    elements.push(VectorElement::Line {
                        x1: cx + nx * w,
                        y1: cy + ny * w,
                        x2: cx - nx * w,
                        y2: cy - ny * w,
                        width: (stroke * 0.85).max(1.0),
                        color: fg,
                        dashed: false,
                    });
                }
                continue;
            }
            BondStereo::Wavy => {
                let amp = (scale * 0.09).clamp(1.8, 3.4);
                let segs = 12usize;
                let mut pts = Vec::with_capacity(segs + 1);
                for s in 0..=segs {
                    let t = s as f64 / segs as f64;
                    let wave = (t * std::f64::consts::PI * 3.0).sin() * amp;
                    pts.push([p1[0] + dx * t + nx * wave, p1[1] + dy * t + ny * wave]);
                }
                elements.push(VectorElement::Polyline {
                    points: pts,
                    width: stroke,
                    color: fg,
                    closed: false,
                });
                continue;
            }
            _ => {}
        }

        match b.order {
            BondOrder::Single => {
                elements.push(VectorElement::Line {
                    x1: p1[0],
                    y1: p1[1],
                    x2: p2[0],
                    y2: p2[1],
                    width: stroke,
                    color: fg,
                    dashed: false,
                });
            }
            BondOrder::Dative => {
                elements.push(VectorElement::Line {
                    x1: p1[0],
                    y1: p1[1],
                    x2: p2[0],
                    y2: p2[1],
                    width: stroke,
                    color: fg,
                    dashed: false,
                });
                let head = (scale * 0.22).clamp(4.0, 7.5);
                let tip = p2;
                let base = [p2[0] - ux * head, p2[1] - uy * head];
                elements.push(VectorElement::Polygon {
                    points: vec![
                        tip,
                        [base[0] + nx * (head * 0.45), base[1] + ny * (head * 0.45)],
                        [base[0] - nx * (head * 0.45), base[1] - ny * (head * 0.45)],
                    ],
                    fill: fg,
                });
            }
            BondOrder::Double => {
                let off = (scale * 0.12).clamp(2.2, 4.2);
                if in_ring_centroid.is_some() {
                    // Ring-aligned double bond: outer main stroke + shorter inner parallel stroke
                    elements.push(VectorElement::Line {
                        x1: p1[0],
                        y1: p1[1],
                        x2: p2[0],
                        y2: p2[1],
                        width: stroke,
                        color: fg,
                        dashed: false,
                    });
                    let inset = (len * 0.14).min(5.0);
                    elements.push(VectorElement::Line {
                        x1: p1[0] + ux * inset + nx * (off * 1.45),
                        y1: p1[1] + uy * inset + ny * (off * 1.45),
                        x2: p2[0] - ux * inset + nx * (off * 1.45),
                        y2: p2[1] - uy * inset + ny * (off * 1.45),
                        width: stroke * 0.9,
                        color: fg,
                        dashed: false,
                    });
                } else {
                    // Symmetric double bond
                    elements.push(VectorElement::Line {
                        x1: p1[0] + nx * (off * 0.75),
                        y1: p1[1] + ny * (off * 0.75),
                        x2: p2[0] + nx * (off * 0.75),
                        y2: p2[1] + ny * (off * 0.75),
                        width: stroke * 0.9,
                        color: fg,
                        dashed: false,
                    });
                    elements.push(VectorElement::Line {
                        x1: p1[0] - nx * (off * 0.75),
                        y1: p1[1] - ny * (off * 0.75),
                        x2: p2[0] - nx * (off * 0.75),
                        y2: p2[1] - ny * (off * 0.75),
                        width: stroke * 0.9,
                        color: fg,
                        dashed: false,
                    });
                }
            }
            BondOrder::Triple => {
                let off = (scale * 0.13).clamp(2.3, 4.5);
                elements.push(VectorElement::Line {
                    x1: p1[0],
                    y1: p1[1],
                    x2: p2[0],
                    y2: p2[1],
                    width: stroke * 0.9,
                    color: fg,
                    dashed: false,
                });
                elements.push(VectorElement::Line {
                    x1: p1[0] + nx * off,
                    y1: p1[1] + ny * off,
                    x2: p2[0] + nx * off,
                    y2: p2[1] + ny * off,
                    width: stroke * 0.85,
                    color: fg,
                    dashed: false,
                });
                elements.push(VectorElement::Line {
                    x1: p1[0] - nx * off,
                    y1: p1[1] - ny * off,
                    x2: p2[0] - nx * off,
                    y2: p2[1] - ny * off,
                    width: stroke * 0.85,
                    color: fg,
                    dashed: false,
                });
            }
            BondOrder::Aromatic => {
                let off = (scale * 0.13).clamp(2.2, 4.3);
                elements.push(VectorElement::Line {
                    x1: p1[0],
                    y1: p1[1],
                    x2: p2[0],
                    y2: p2[1],
                    width: stroke,
                    color: fg,
                    dashed: false,
                });
                let inset = (len * 0.14).min(5.0);
                elements.push(VectorElement::Line {
                    x1: p1[0] + ux * inset + nx * (off * 1.4),
                    y1: p1[1] + uy * inset + ny * (off * 1.4),
                    x2: p2[0] - ux * inset + nx * (off * 1.4),
                    y2: p2[1] - uy * inset + ny * (off * 1.4),
                    width: stroke * 0.85,
                    color: fg,
                    dashed: true,
                });
            }
        }
    }

    // 3. Emit atom labels on top
    for item in visible_labels.into_iter().flatten() {
        let (text, color, p, bbox) = item;
        label_boxes.push(bbox);
        elements.push(VectorElement::Label {
            text,
            cx: p[0],
            cy: p[1],
            font_size: eff_font,
            color,
            bbox,
        });
    }

    // 4. Optional Reaction Arrows (with smooth VectorCraft cubic arc / shaft)
    for arr in arrows {
        let p0 = arr.start;
        let p3 = arr.end;
        let dx = p3[0] - p0[0];
        let dy = p3[1] - p0[1];
        let d = dx.hypot(dy).max(1.0);
        let ux = dx / d;
        let uy = dy / d;
        let nx = -uy;
        let ny = ux;
        elements.push(VectorElement::CubicCurve {
            p0,
            p1: [p0[0] + dx * 0.33, p0[1] + dy * 0.33],
            p2: [p0[0] + dx * 0.66, p0[1] + dy * 0.66],
            p3,
            width: stroke,
            color: fg,
        });
        let head = 8.0;
        elements.push(VectorElement::Polygon {
            points: vec![
                p3,
                [p3[0] - ux * head + nx * 4.0, p3[1] - uy * head + ny * 4.0],
                [p3[0] - ux * head - nx * 4.0, p3[1] - uy * head - ny * 4.0],
            ],
            fill: fg,
        });
        if !arr.label_above.is_empty() {
            let mid = [(p0[0] + p3[0]) * 0.5, (p0[1] + p3[1]) * 0.5 - 12.0];
            let bw = arr.label_above.len() as f64 * 6.5;
            let bbox = Rect2 {
                x: mid[0] - bw * 0.5,
                y: mid[1] - 6.0,
                w: bw,
                h: 12.0,
            };
            label_boxes.push(bbox);
            elements.push(VectorElement::Label {
                text: arr.label_above.clone(),
                cx: mid[0],
                cy: mid[1],
                font_size: 10.5,
                color: fg,
                bbox,
            });
        }
    }

    VectorScene {
        width: style.width,
        height: style.height,
        background: bg,
        elements,
        atom_screen_pos: atom_pos,
        label_boxes,
    }
}

fn format_atom_label(a: &MolAtom, mol: &MolGraph, pos: &[[f64; 2]]) -> Option<String> {
    let show_carbon = mol.adj[a.id].is_empty() || a.formal_charge != 0 || a.isotope.is_some();
    if a.symbol == "C" && !show_carbon {
        return None;
    }

    // Decide whether implicit/explicit H sits on left (HO, H2N) vs right (OH, NH2)
    let h_count = a.total_h();
    let mut sum_dx = 0.0;
    for &(nbr, _) in &mol.adj[a.id] {
        sum_dx += pos[nbr][0] - pos[a.id][0];
    }
    let h_on_left = mol.adj[a.id].len() == 1 && sum_dx > 2.0;

    let h_str = match h_count {
        0 => String::new(),
        1 => "H".to_string(),
        n => format!("H{n}"),
    };
    let mut out = if h_on_left && !h_str.is_empty() {
        format!("{h_str}{}", a.symbol)
    } else {
        format!("{}{h_str}", a.symbol)
    };

    if a.formal_charge > 0 {
        if a.formal_charge == 1 {
            out.push('+');
        } else {
            out.push_str(&format!("+{}", a.formal_charge));
        }
    } else if a.formal_charge < 0 {
        let mag = i16::from(a.formal_charge).abs();
        if mag == 1 {
            out.push('-');
        } else {
            out.push_str(&format!("-{mag}"));
        }
    }
    Some(out)
}

fn trim_segment_to_labels(
    p1: [f64; 2],
    p2: [f64; 2],
    b1: Option<&Rect2>,
    b2: Option<&Rect2>,
) -> ([f64; 2], [f64; 2]) {
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let d = dx.hypot(dy).max(1e-6);
    let ux = dx / d;
    let uy = dy / d;

    let trim1 = b1.map_or(0.0, |r| {
        let rx = (r.w * 0.55) / ux.abs().max(0.25);
        let ry = (r.h * 0.55) / uy.abs().max(0.25);
        rx.min(ry).min(d * 0.35)
    });
    let trim2 = b2.map_or(0.0, |r| {
        let rx = (r.w * 0.55) / ux.abs().max(0.25);
        let ry = (r.h * 0.55) / uy.abs().max(0.25);
        rx.min(ry).min(d * 0.35)
    });

    (
        [p1[0] + ux * trim1, p1[1] + uy * trim1],
        [p2[0] - ux * trim2, p2[1] - uy * trim2],
    )
}

fn element_color(sym: &str, dark: bool) -> Rgba {
    match sym {
        "N" => {
            if dark {
                Rgba(116, 178, 255, 255)
            } else {
                Rgba(24, 92, 204, 255)
            }
        }
        "O" => {
            if dark {
                Rgba(255, 118, 118, 255)
            } else {
                Rgba(204, 34, 34, 255)
            }
        }
        "S" => {
            if dark {
                Rgba(235, 205, 84, 255)
            } else {
                Rgba(184, 134, 11, 255)
            }
        }
        "P" => {
            if dark {
                Rgba(255, 168, 76, 255)
            } else {
                Rgba(204, 102, 0, 255)
            }
        }
        "F" | "Cl" => {
            if dark {
                Rgba(102, 217, 134, 255)
            } else {
                Rgba(22, 148, 56, 255)
            }
        }
        "Br" => {
            if dark {
                Rgba(224, 128, 96, 255)
            } else {
                Rgba(158, 42, 43, 255)
            }
        }
        "I" => {
            if dark {
                Rgba(198, 132, 240, 255)
            } else {
                Rgba(130, 44, 184, 255)
            }
        }
        _ => {
            if dark {
                Rgba(232, 236, 242, 255)
            } else {
                Rgba(28, 30, 36, 255)
            }
        }
    }
}
