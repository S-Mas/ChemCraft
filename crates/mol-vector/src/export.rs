use std::fmt::Write as _;
use std::io::Write as _;

use flate2::Compression;
use flate2::write::ZlibEncoder;
use thiserror::Error;
use tiny_skia::{
    Color, FillRule, LineCap, LineJoin, Paint, PathBuilder, Pixmap, Stroke, StrokeDash, Transform,
};

use mol_graph::MolGraph;

use crate::scene::{Rgba, SceneStyle, VectorElement, VectorScene, build_mol_scene};

/// Errors produced by `mol-vector` exports.
#[derive(Debug, Error)]
pub enum VectorError {
    #[error("compression error: {0}")]
    Compression(String),
    #[error("PNG encoding error: {0}")]
    Png(String),
    #[error("invalid raster dimensions")]
    Dimensions,
}

/// Export a single `VectorScene` as clean Inkscape/VectorCraft-compatible SVG.
pub fn scene_to_svg(scene: &VectorScene) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.1}\" height=\"{:.1}\" viewBox=\"0 0 {:.1} {:.1}\">",
        scene.width, scene.height, scene.width, scene.height
    );
    if let Some(bg) = scene.background {
        let _ = writeln!(
            out,
            "<rect width=\"100%\" height=\"100%\" fill=\"{}\"/>",
            hex_color(bg)
        );
    }
    append_svg_elements(&mut out, scene, 0.0, 0.0);
    let _ = writeln!(out, "</svg>");
    out
}

/// Export a multi-molecule SAR grid figure as clean SVG.
pub fn export_sar_grid_svg(
    items: &[(&str, &MolGraph)],
    cols: usize,
    cell_w: f64,
    cell_h: f64,
    dark_mode: bool,
) -> String {
    let cols = cols.max(1);
    let rows = items.len().div_ceil(cols).max(1);
    let total_w = cols as f64 * cell_w;
    let total_h = rows as f64 * cell_h;

    let mut out = String::new();
    let _ = writeln!(
        out,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{total_w:.1}\" height=\"{total_h:.1}\" viewBox=\"0 0 {total_w:.1} {total_h:.1}\">"
    );
    let bg = if dark_mode { "#181a20" } else { "#ffffff" };
    let border = if dark_mode { "#343a40" } else { "#dee2e6" };
    let title_col = if dark_mode { "#e8ecf2" } else { "#1c1e24" };
    let _ = writeln!(out, "<rect width=\"100%\" height=\"100%\" fill=\"{bg}\"/>");

    let style = SceneStyle {
        width: cell_w,
        height: cell_h - 24.0,
        padding: 18.0,
        stroke_width: 1.5,
        font_size: 11.5,
        dark_mode,
        show_stereo_labels: true,
        ..SceneStyle::default()
    };

    for (idx, (title, mol)) in items.iter().enumerate() {
        let r = idx / cols;
        let c = idx % cols;
        let ox = c as f64 * cell_w;
        let oy = r as f64 * cell_h;
        let _ = writeln!(
            out,
            "<rect x=\"{ox:.1}\" y=\"{oy:.1}\" width=\"{cell_w:.1}\" height=\"{cell_h:.1}\" fill=\"none\" stroke=\"{border}\" stroke-width=\"1\"/>"
        );
        let scene = build_mol_scene(mol, &style, &[]);
        append_svg_elements(&mut out, &scene, ox, oy);
        let safe_title = escape_xml(&truncate_str(title, 28));
        let _ = writeln!(
            out,
            "<text x=\"{:.1}\" y=\"{:.1}\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"10.5\" font-weight=\"600\" text-anchor=\"middle\" fill=\"{title_col}\">{safe_title}</text>",
            ox + cell_w * 0.5,
            oy + cell_h - 8.0
        );
    }

    let _ = writeln!(out, "</svg>");
    out
}

fn append_svg_elements(out: &mut String, scene: &VectorScene, ox: f64, oy: f64) {
    for el in &scene.elements {
        match el {
            VectorElement::Line {
                x1,
                y1,
                x2,
                y2,
                width,
                color,
                dashed,
            } => {
                let dash_attr = if *dashed {
                    " stroke-dasharray=\"3.5,3\""
                } else {
                    ""
                };
                let _ = writeln!(
                    out,
                    "<line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\" stroke=\"{}\" stroke-width=\"{width:.2}\" stroke-linecap=\"round\"{dash_attr}/>",
                    x1 + ox,
                    y1 + oy,
                    x2 + ox,
                    y2 + oy,
                    hex_color(*color)
                );
            }
            VectorElement::Polyline {
                points,
                width,
                color,
                closed,
            } => {
                let pts: String = points
                    .iter()
                    .map(|p| format!("{:.2},{:.2}", p[0] + ox, p[1] + oy))
                    .collect::<Vec<_>>()
                    .join(" ");
                let tag = if *closed { "polygon" } else { "polyline" };
                let _ = writeln!(
                    out,
                    "<{tag} points=\"{pts}\" fill=\"none\" stroke=\"{}\" stroke-width=\"{width:.2}\" stroke-linejoin=\"round\"/>",
                    hex_color(*color)
                );
            }
            VectorElement::Polygon { points, fill } => {
                let pts: String = points
                    .iter()
                    .map(|p| format!("{:.2},{:.2}", p[0] + ox, p[1] + oy))
                    .collect::<Vec<_>>()
                    .join(" ");
                let _ = writeln!(
                    out,
                    "<polygon points=\"{pts}\" fill=\"{}\"/>",
                    hex_color(*fill)
                );
            }
            VectorElement::Circle {
                cx,
                cy,
                r,
                fill,
                stroke,
            } => {
                let fill_s = fill.map_or_else(|| "none".to_string(), hex_color);
                let stroke_attr = stroke.map_or_else(String::new, |(c, w)| {
                    format!(" stroke=\"{}\" stroke-width=\"{w:.2}\"", hex_color(c))
                });
                let _ = writeln!(
                    out,
                    "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{r:.2}\" fill=\"{fill_s}\"{stroke_attr}/>",
                    cx + ox,
                    cy + oy
                );
            }
            VectorElement::CubicCurve {
                p0,
                p1,
                p2,
                p3,
                width,
                color,
            } => {
                let _ = writeln!(
                    out,
                    "<path d=\"M {:.2} {:.2} C {:.2} {:.2}, {:.2} {:.2}, {:.2} {:.2}\" fill=\"none\" stroke=\"{}\" stroke-width=\"{width:.2}\" stroke-linecap=\"round\"/>",
                    p0[0] + ox,
                    p0[1] + oy,
                    p1[0] + ox,
                    p1[1] + oy,
                    p2[0] + ox,
                    p2[1] + oy,
                    p3[0] + ox,
                    p3[1] + oy,
                    hex_color(*color)
                );
            }
            VectorElement::Label {
                text,
                cx,
                cy,
                font_size,
                color,
                ..
            } => {
                let safe = escape_xml(text);
                let _ = writeln!(
                    out,
                    "<text x=\"{:.2}\" y=\"{:.2}\" font-family=\"Helvetica,Arial,sans-serif\" font-size=\"{font_size:.1}\" font-weight=\"700\" text-anchor=\"middle\" dominant-baseline=\"central\" fill=\"{}\">{safe}</text>",
                    cx + ox,
                    cy + oy,
                    hex_color(*color)
                );
            }
        }
    }
}

/// Export a `VectorScene` as Encapsulated PostScript (EPS 3.0) vector graphics.
pub fn scene_to_eps(scene: &VectorScene) -> String {
    let w = scene.width;
    let h = scene.height;
    let mut out = String::new();
    let _ = writeln!(out, "%!PS-Adobe-3.0 EPSF-3.0");
    let _ = writeln!(
        out,
        "%%BoundingBox: 0 0 {} {}",
        w.ceil() as i64,
        h.ceil() as i64
    );
    let _ = writeln!(out, "%%Creator: ChemCraft mol-vector");
    let _ = writeln!(out, "1 setlinecap 1 setlinejoin");
    if let Some(bg) = scene.background {
        let _ = writeln!(
            out,
            "{:.3} {:.3} {:.3} setrgbcolor 0 0 {w:.2} {h:.2} rectfill",
            f64::from(bg.0) / 255.0,
            f64::from(bg.1) / 255.0,
            f64::from(bg.2) / 255.0
        );
    }

    append_eps_scene(&mut out, scene, 0.0, 0.0, h);
    let _ = writeln!(out, "showpage\n%%EOF");
    out
}

/// Export a multi-molecule SAR grid as Encapsulated PostScript (EPS 3.0).
pub fn export_sar_grid_eps(
    items: &[(&str, &MolGraph)],
    cols: usize,
    cell_w: f64,
    cell_h: f64,
    dark_mode: bool,
) -> String {
    let cols = cols.max(1);
    let rows = items.len().div_ceil(cols).max(1);
    let total_w = cols as f64 * cell_w;
    let total_h = rows as f64 * cell_h;

    let mut out = String::new();
    let _ = writeln!(out, "%!PS-Adobe-3.0 EPSF-3.0");
    let _ = writeln!(
        out,
        "%%BoundingBox: 0 0 {} {}",
        total_w.ceil() as i64,
        total_h.ceil() as i64
    );
    let _ = writeln!(out, "%%Creator: ChemCraft mol-vector SAR Grid");
    let _ = writeln!(out, "1 setlinecap 1 setlinejoin");
    if dark_mode {
        let _ = writeln!(
            out,
            "0.094 0.102 0.125 setrgbcolor 0 0 {total_w:.2} {total_h:.2} rectfill"
        );
    } else {
        let _ = writeln!(
            out,
            "1 1 1 setrgbcolor 0 0 {total_w:.2} {total_h:.2} rectfill"
        );
    }

    let style = SceneStyle {
        width: cell_w,
        height: cell_h - 24.0,
        padding: 18.0,
        stroke_width: 1.4,
        font_size: 11.0,
        dark_mode,
        show_stereo_labels: true,
        ..SceneStyle::default()
    };

    for (idx, (title, mol)) in items.iter().enumerate() {
        let r = idx / cols;
        let c = idx % cols;
        let ox = c as f64 * cell_w;
        let oy = r as f64 * cell_h;
        let py_bot = total_h - (oy + cell_h);
        let _ = writeln!(
            out,
            "0.82 0.85 0.88 setrgbcolor 0.8 setlinewidth newpath {ox:.2} {py_bot:.2} {cell_w:.2} {cell_h:.2} rectstroke"
        );
        let scene = build_mol_scene(mol, &style, &[]);
        append_eps_scene(&mut out, &scene, ox, oy, total_h);
        let safe = escape_pdf_ps_text(&truncate_str(title, 26));
        let title_col = if dark_mode {
            "0.91 0.92 0.95"
        } else {
            "0.11 0.12 0.14"
        };
        let tx = ox + cell_w * 0.5;
        let ty = total_h - (oy + cell_h - 9.0);
        let _ = writeln!(
            out,
            "/Helvetica-Bold findfont 10 scalefont setfont {title_col} setrgbcolor {tx:.2} {ty:.2} moveto ({safe}) dup stringwidth pop -0.5 mul 0 rmoveto show"
        );
    }

    let _ = writeln!(out, "showpage\n%%EOF");
    out
}

fn append_eps_scene(out: &mut String, scene: &VectorScene, ox: f64, oy: f64, h: f64) {
    for el in &scene.elements {
        match el {
            VectorElement::Line {
                x1,
                y1,
                x2,
                y2,
                width,
                color,
                dashed,
            } => {
                if *dashed {
                    let _ = writeln!(out, "[3 3] 0 setdash");
                }
                let _ = writeln!(
                    out,
                    "{:.3} {:.3} {:.3} setrgbcolor {width:.2} setlinewidth newpath {:.2} {:.2} moveto {:.2} {:.2} lineto stroke",
                    f64::from(color.0) / 255.0,
                    f64::from(color.1) / 255.0,
                    f64::from(color.2) / 255.0,
                    x1 + ox,
                    h - (y1 + oy),
                    x2 + ox,
                    h - (y2 + oy)
                );
                if *dashed {
                    let _ = writeln!(out, "[] 0 setdash");
                }
            }
            VectorElement::Polyline {
                points,
                width,
                color,
                ..
            } => {
                if let Some(first) = points.first() {
                    let _ = writeln!(
                        out,
                        "{:.3} {:.3} {:.3} setrgbcolor {width:.2} setlinewidth newpath {:.2} {:.2} moveto",
                        f64::from(color.0) / 255.0,
                        f64::from(color.1) / 255.0,
                        f64::from(color.2) / 255.0,
                        first[0] + ox,
                        h - (first[1] + oy)
                    );
                    for p in &points[1..] {
                        let _ = writeln!(out, "{:.2} {:.2} lineto", p[0] + ox, h - (p[1] + oy));
                    }
                    let _ = writeln!(out, "stroke");
                }
            }
            VectorElement::Polygon { points, fill } => {
                if let Some(first) = points.first() {
                    let _ = writeln!(
                        out,
                        "{:.3} {:.3} {:.3} setrgbcolor newpath {:.2} {:.2} moveto",
                        f64::from(fill.0) / 255.0,
                        f64::from(fill.1) / 255.0,
                        f64::from(fill.2) / 255.0,
                        first[0] + ox,
                        h - (first[1] + oy)
                    );
                    for p in &points[1..] {
                        let _ = writeln!(out, "{:.2} {:.2} lineto", p[0] + ox, h - (p[1] + oy));
                    }
                    let _ = writeln!(out, "closepath fill");
                }
            }
            VectorElement::CubicCurve {
                p0,
                p1,
                p2,
                p3,
                width,
                color,
            } => {
                let _ = writeln!(
                    out,
                    "{:.3} {:.3} {:.3} setrgbcolor {width:.2} setlinewidth newpath {:.2} {:.2} moveto {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} curveto stroke",
                    f64::from(color.0) / 255.0,
                    f64::from(color.1) / 255.0,
                    f64::from(color.2) / 255.0,
                    p0[0] + ox,
                    h - (p0[1] + oy),
                    p1[0] + ox,
                    h - (p1[1] + oy),
                    p2[0] + ox,
                    h - (p2[1] + oy),
                    p3[0] + ox,
                    h - (p3[1] + oy)
                );
            }
            VectorElement::Circle {
                cx,
                cy,
                r,
                fill,
                stroke,
            } => {
                if let Some(f) = fill {
                    let _ = writeln!(
                        out,
                        "{:.3} {:.3} {:.3} setrgbcolor newpath {:.2} {:.2} {r:.2} 0 360 arc fill",
                        f64::from(f.0) / 255.0,
                        f64::from(f.1) / 255.0,
                        f64::from(f.2) / 255.0,
                        cx + ox,
                        h - (cy + oy)
                    );
                }
                if let Some((c, w)) = stroke {
                    let _ = writeln!(
                        out,
                        "{:.3} {:.3} {:.3} setrgbcolor {w:.2} setlinewidth newpath {:.2} {:.2} {r:.2} 0 360 arc stroke",
                        f64::from(c.0) / 255.0,
                        f64::from(c.1) / 255.0,
                        f64::from(c.2) / 255.0,
                        cx + ox,
                        h - (cy + oy)
                    );
                }
            }
            VectorElement::Label {
                text,
                cx,
                cy,
                font_size,
                color,
                ..
            } => {
                let esc = escape_pdf_ps_text(text);
                let _ = writeln!(
                    out,
                    "/Helvetica-Bold findfont {font_size:.1} scalefont setfont {:.3} {:.3} {:.3} setrgbcolor {:.2} {:.2} moveto ({esc}) dup stringwidth pop -0.5 mul {:.2} rmoveto show",
                    f64::from(color.0) / 255.0,
                    f64::from(color.1) / 255.0,
                    f64::from(color.2) / 255.0,
                    cx + ox,
                    h - (cy + oy),
                    -font_size * 0.32
                );
            }
        }
    }
}

/// Export a `VectorScene` as a compressed pure-Rust `FlateDecode` vector PDF 1.4.
pub fn scene_to_pdf_bytes(scene: &VectorScene) -> Result<Vec<u8>, VectorError> {
    let mut stream = String::new();
    let h = scene.height;
    let _ = writeln!(stream, "1 J 1 j");
    if let Some(bg) = scene.background {
        let _ = writeln!(
            stream,
            "{:.3} {:.3} {:.3} rg 0 0 {:.2} {:.2} re f",
            f64::from(bg.0) / 255.0,
            f64::from(bg.1) / 255.0,
            f64::from(bg.2) / 255.0,
            scene.width,
            scene.height
        );
    }
    append_pdf_scene(&mut stream, scene, 0.0, 0.0, h);
    build_single_page_pdf(scene.width, scene.height, stream.as_bytes())
}

/// Export a multi-molecule SAR grid as compressed vector PDF.
pub fn export_sar_grid_pdf_bytes(
    items: &[(&str, &MolGraph)],
    cols: usize,
    cell_w: f64,
    cell_h: f64,
) -> Result<Vec<u8>, VectorError> {
    let cols = cols.max(1);
    let rows = items.len().div_ceil(cols).max(1);
    let total_w = cols as f64 * cell_w;
    let total_h = rows as f64 * cell_h;

    let style = SceneStyle {
        width: cell_w,
        height: cell_h - 24.0,
        padding: 18.0,
        stroke_width: 1.4,
        font_size: 11.0,
        dark_mode: false,
        show_stereo_labels: true,
        ..SceneStyle::default()
    };

    let mut stream = String::new();
    let _ = writeln!(
        stream,
        "1 J 1 j 1 1 1 rg 0 0 {total_w:.2} {total_h:.2} re f"
    );
    for (idx, (title, mol)) in items.iter().enumerate() {
        let r = idx / cols;
        let c = idx % cols;
        let ox = c as f64 * cell_w;
        let oy = r as f64 * cell_h;
        let scene = build_mol_scene(mol, &style, &[]);
        let _ = writeln!(
            stream,
            "0.85 0.87 0.90 RG 0.8 w {ox:.2} {:.2} {cell_w:.2} {cell_h:.2} re S",
            total_h - (oy + cell_h)
        );
        append_pdf_scene(&mut stream, &scene, ox, oy, total_h);
        let safe = escape_pdf_ps_text(&truncate_str(title, 26));
        let tx = ox + cell_w * 0.5 - (safe.len() as f64 * 2.7);
        let ty = total_h - (oy + cell_h - 9.0);
        let _ = writeln!(
            stream,
            "0.12 0.14 0.18 rg BT /F2 9.5 Tf {tx:.2} {ty:.2} Td ({safe}) Tj ET"
        );
    }

    build_single_page_pdf(total_w, total_h, stream.as_bytes())
}

fn append_pdf_scene(stream: &mut String, scene: &VectorScene, ox: f64, oy: f64, page_h: f64) {
    for el in &scene.elements {
        match el {
            VectorElement::Line {
                x1,
                y1,
                x2,
                y2,
                width,
                color,
                dashed,
            } => {
                if *dashed {
                    let _ = writeln!(stream, "[3 2.5] 0 d");
                }
                let _ = writeln!(
                    stream,
                    "{:.3} {:.3} {:.3} RG {width:.2} w {:.2} {:.2} m {:.2} {:.2} l S",
                    f64::from(color.0) / 255.0,
                    f64::from(color.1) / 255.0,
                    f64::from(color.2) / 255.0,
                    x1 + ox,
                    page_h - (y1 + oy),
                    x2 + ox,
                    page_h - (y2 + oy)
                );
                if *dashed {
                    let _ = writeln!(stream, "[] 0 d");
                }
            }
            VectorElement::Polyline {
                points,
                width,
                color,
                ..
            } => {
                if let Some(first) = points.first() {
                    let _ = writeln!(
                        stream,
                        "{:.3} {:.3} {:.3} RG {width:.2} w {:.2} {:.2} m",
                        f64::from(color.0) / 255.0,
                        f64::from(color.1) / 255.0,
                        f64::from(color.2) / 255.0,
                        first[0] + ox,
                        page_h - (first[1] + oy)
                    );
                    for p in &points[1..] {
                        let _ = writeln!(stream, "{:.2} {:.2} l", p[0] + ox, page_h - (p[1] + oy));
                    }
                    let _ = writeln!(stream, "S");
                }
            }
            VectorElement::Polygon { points, fill } => {
                if let Some(first) = points.first() {
                    let _ = writeln!(
                        stream,
                        "{:.3} {:.3} {:.3} rg {:.2} {:.2} m",
                        f64::from(fill.0) / 255.0,
                        f64::from(fill.1) / 255.0,
                        f64::from(fill.2) / 255.0,
                        first[0] + ox,
                        page_h - (first[1] + oy)
                    );
                    for p in &points[1..] {
                        let _ = writeln!(stream, "{:.2} {:.2} l", p[0] + ox, page_h - (p[1] + oy));
                    }
                    let _ = writeln!(stream, "f");
                }
            }
            VectorElement::CubicCurve {
                p0,
                p1,
                p2,
                p3,
                width,
                color,
            } => {
                let _ = writeln!(
                    stream,
                    "{:.3} {:.3} {:.3} RG {width:.2} w {:.2} {:.2} m {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c S",
                    f64::from(color.0) / 255.0,
                    f64::from(color.1) / 255.0,
                    f64::from(color.2) / 255.0,
                    p0[0] + ox,
                    page_h - (p0[1] + oy),
                    p1[0] + ox,
                    page_h - (p1[1] + oy),
                    p2[0] + ox,
                    page_h - (p2[1] + oy),
                    p3[0] + ox,
                    page_h - (p3[1] + oy)
                );
            }
            VectorElement::Circle {
                cx,
                cy,
                r,
                fill,
                stroke,
            } => {
                let x = cx + ox;
                let y = page_h - (cy + oy);
                let k = r * 0.552_284_749_8;
                if let Some(f) = fill {
                    let _ = writeln!(
                        stream,
                        "{:.3} {:.3} {:.3} rg {:.2} {:.2} m {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c f",
                        f64::from(f.0) / 255.0,
                        f64::from(f.1) / 255.0,
                        f64::from(f.2) / 255.0,
                        x + r,
                        y,
                        x + r,
                        y + k,
                        x + k,
                        y + r,
                        x,
                        y + r,
                        x - k,
                        y + r,
                        x - r,
                        y + k,
                        x - r,
                        y,
                        x - r,
                        y - k,
                        x - k,
                        y - r,
                        x,
                        y - r,
                        x + k,
                        y - r,
                        x + r,
                        y - k,
                        x + r,
                        y
                    );
                }
                if let Some((c, w)) = stroke {
                    let _ = writeln!(
                        stream,
                        "{:.3} {:.3} {:.3} RG {w:.2} w {:.2} {:.2} m {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {:.2} c S",
                        f64::from(c.0) / 255.0,
                        f64::from(c.1) / 255.0,
                        f64::from(c.2) / 255.0,
                        x + r,
                        y,
                        x + r,
                        y + k,
                        x + k,
                        y + r,
                        x,
                        y + r,
                        x - k,
                        y + r,
                        x - r,
                        y + k,
                        x - r,
                        y,
                        x - r,
                        y - k,
                        x - k,
                        y - r,
                        x,
                        y - r,
                        x + k,
                        y - r,
                        x + r,
                        y - k,
                        x + r,
                        y
                    );
                }
            }
            VectorElement::Label {
                text,
                cx,
                cy,
                font_size,
                color,
                ..
            } => {
                let esc = escape_pdf_ps_text(text);
                let tx = cx + ox - (text.len() as f64 * font_size * 0.30);
                let ty = page_h - (cy + oy) - font_size * 0.32;
                let _ = writeln!(
                    stream,
                    "{:.3} {:.3} {:.3} rg BT /F2 {font_size:.1} Tf {tx:.2} {ty:.2} Td ({esc}) Tj ET",
                    f64::from(color.0) / 255.0,
                    f64::from(color.1) / 255.0,
                    f64::from(color.2) / 255.0
                );
            }
        }
    }
}

fn build_single_page_pdf(w: f64, h: f64, content_raw: &[u8]) -> Result<Vec<u8>, VectorError> {
    let mut enc = ZlibEncoder::new(Vec::new(), Compression::default());
    enc.write_all(content_raw)
        .map_err(|e| VectorError::Compression(e.to_string()))?;
    let compressed = enc
        .finish()
        .map_err(|e| VectorError::Compression(e.to_string()))?;

    let mut pdf = Vec::new();
    pdf.extend_from_slice(b"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n");
    let mut offsets = vec![0usize];

    let obj1 = b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n";
    offsets.push(pdf.len());
    pdf.extend_from_slice(obj1);

    let obj2 = b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n";
    offsets.push(pdf.len());
    pdf.extend_from_slice(obj2);

    let obj3 = format!(
        "3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {w:.2} {h:.2}] /Resources << /Font << /F1 4 0 R /F2 5 0 R >> >> /Contents 6 0 R >>\nendobj\n"
    );
    offsets.push(pdf.len());
    pdf.extend_from_slice(obj3.as_bytes());

    let obj4 = b"4 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>\nendobj\n";
    offsets.push(pdf.len());
    pdf.extend_from_slice(obj4);

    let obj5 = b"5 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold /Encoding /WinAnsiEncoding >>\nendobj\n";
    offsets.push(pdf.len());
    pdf.extend_from_slice(obj5);

    offsets.push(pdf.len());
    let head6 = format!(
        "6 0 obj\n<< /Length {} /Filter /FlateDecode >>\nstream\n",
        compressed.len()
    );
    pdf.extend_from_slice(head6.as_bytes());
    pdf.extend_from_slice(&compressed);
    pdf.extend_from_slice(b"\nendstream\nendobj\n");

    let xref_pos = pdf.len();
    let _ = write_pdf_xref(&mut pdf, &offsets, xref_pos);
    Ok(pdf)
}

fn write_pdf_xref(
    pdf: &mut Vec<u8>,
    offsets: &[usize],
    xref_pos: usize,
) -> Result<(), VectorError> {
    let count = offsets.len();
    let mut s = format!("xref\n0 {count}\n0000000000 65535 f \n");
    for &off in &offsets[1..] {
        let _ = writeln!(s, "{off:010} 00000 n ");
    }
    let _ = writeln!(
        s,
        "trailer\n<< /Size {count} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF"
    );
    pdf.extend_from_slice(s.as_bytes());
    Ok(())
}

/// Rasterize a `VectorScene` to an anti-aliased RGBA PNG buffer using pure-Rust `tiny-skia` + `png`.
pub fn scene_to_png_bytes(scene: &VectorScene, scale: f64) -> Result<Vec<u8>, VectorError> {
    let sc = scale.clamp(0.5, 4.0);
    let px_w = ((scene.width * sc).round() as u32).max(16);
    let px_h = ((scene.height * sc).round() as u32).max(16);
    let mut pixmap = Pixmap::new(px_w, px_h).ok_or(VectorError::Dimensions)?;
    if let Some(bg) = scene.background {
        pixmap.fill(Color::from_rgba8(bg.0, bg.1, bg.2, bg.3));
    }
    let xform = Transform::from_scale(sc as f32, sc as f32);
    render_scene_onto_pixmap(&mut pixmap, scene, xform);
    encode_pixmap_png(&pixmap)
}

/// Rasterize a multi-molecule SAR grid to an anti-aliased RGBA PNG buffer.
pub fn export_sar_grid_png_bytes(
    items: &[(&str, &MolGraph)],
    cols: usize,
    cell_w: f64,
    cell_h: f64,
    dark_mode: bool,
    scale: f64,
) -> Result<Vec<u8>, VectorError> {
    let cols = cols.max(1);
    let rows = items.len().div_ceil(cols).max(1);
    let total_w = cols as f64 * cell_w;
    let total_h = rows as f64 * cell_h;
    let sc = scale.clamp(0.5, 4.0);
    let px_w = ((total_w * sc).round() as u32).max(16);
    let px_h = ((total_h * sc).round() as u32).max(16);
    let mut pixmap = Pixmap::new(px_w, px_h).ok_or(VectorError::Dimensions)?;
    let bg = if dark_mode {
        Color::from_rgba8(24, 26, 32, 255)
    } else {
        Color::from_rgba8(255, 255, 255, 255)
    };
    pixmap.fill(bg);

    let style = SceneStyle {
        width: cell_w,
        height: cell_h - 24.0,
        padding: 18.0,
        stroke_width: 1.4,
        font_size: 11.0,
        dark_mode,
        show_stereo_labels: true,
        ..SceneStyle::default()
    };

    let border_col = Rgba(210, 215, 222, 255);
    let title_col = if dark_mode {
        Rgba(232, 236, 242, 255)
    } else {
        Rgba(28, 30, 36, 255)
    };

    for (idx, (title, mol)) in items.iter().enumerate() {
        let r = idx / cols;
        let c = idx % cols;
        let ox = c as f64 * cell_w;
        let oy = r as f64 * cell_h;
        let xform = Transform::from_scale(sc as f32, sc as f32).pre_translate(ox as f32, oy as f32);

        let mut pb = PathBuilder::new();
        pb.push_rect(
            tiny_skia::Rect::from_xywh(0.0, 0.0, cell_w as f32, cell_h as f32)
                .ok_or(VectorError::Dimensions)?,
        );
        let mut paint = Paint::default();
        paint.set_color_rgba8(border_col.0, border_col.1, border_col.2, border_col.3);
        if let Some(p) = pb.finish() {
            pixmap.stroke_path(
                &p,
                &paint,
                &Stroke {
                    width: 1.0,
                    ..Stroke::default()
                },
                xform,
                None,
            );
        }
        let scene = build_mol_scene(mol, &style, &[]);
        render_scene_onto_pixmap(&mut pixmap, &scene, xform);
        draw_vector_glyph_label(
            &mut pixmap,
            &truncate_str(title, 22),
            cell_w * 0.5,
            cell_h - 11.0,
            9.5,
            title_col,
            xform,
        );
    }

    encode_pixmap_png(&pixmap)
}

fn encode_pixmap_png(pixmap: &Pixmap) -> Result<Vec<u8>, VectorError> {
    let mut png_bytes = Vec::new();
    {
        let mut encoder = png::Encoder::new(&mut png_bytes, pixmap.width(), pixmap.height());
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        let mut writer = encoder
            .write_header()
            .map_err(|e| VectorError::Png(e.to_string()))?;
        writer
            .write_image_data(pixmap.data())
            .map_err(|e| VectorError::Png(e.to_string()))?;
    }
    Ok(png_bytes)
}

fn render_scene_onto_pixmap(pixmap: &mut Pixmap, scene: &VectorScene, xform: Transform) {
    for el in &scene.elements {
        match el {
            VectorElement::Line {
                x1,
                y1,
                x2,
                y2,
                width,
                color,
                dashed,
            } => {
                let mut pb = PathBuilder::new();
                pb.move_to(*x1 as f32, *y1 as f32);
                pb.line_to(*x2 as f32, *y2 as f32);
                if let Some(path) = pb.finish() {
                    let mut paint = Paint::default();
                    paint.set_color_rgba8(color.0, color.1, color.2, color.3);
                    paint.anti_alias = true;
                    let stroke = Stroke {
                        width: *width as f32,
                        line_cap: LineCap::Round,
                        line_join: LineJoin::Round,
                        dash: if *dashed {
                            StrokeDash::new(vec![3.5, 3.0], 0.0)
                        } else {
                            None
                        },
                        ..Stroke::default()
                    };
                    pixmap.stroke_path(&path, &paint, &stroke, xform, None);
                }
            }
            VectorElement::Polyline {
                points,
                width,
                color,
                closed,
            } => {
                if let Some(first) = points.first() {
                    let mut pb = PathBuilder::new();
                    pb.move_to(first[0] as f32, first[1] as f32);
                    for p in &points[1..] {
                        pb.line_to(p[0] as f32, p[1] as f32);
                    }
                    if *closed {
                        pb.close();
                    }
                    if let Some(path) = pb.finish() {
                        let mut paint = Paint::default();
                        paint.set_color_rgba8(color.0, color.1, color.2, color.3);
                        paint.anti_alias = true;
                        let stroke = Stroke {
                            width: *width as f32,
                            line_cap: LineCap::Round,
                            line_join: LineJoin::Round,
                            ..Stroke::default()
                        };
                        pixmap.stroke_path(&path, &paint, &stroke, xform, None);
                    }
                }
            }
            VectorElement::Polygon { points, fill } => {
                if let Some(first) = points.first() {
                    let mut pb = PathBuilder::new();
                    pb.move_to(first[0] as f32, first[1] as f32);
                    for p in &points[1..] {
                        pb.line_to(p[0] as f32, p[1] as f32);
                    }
                    pb.close();
                    if let Some(path) = pb.finish() {
                        let mut paint = Paint::default();
                        paint.set_color_rgba8(fill.0, fill.1, fill.2, fill.3);
                        paint.anti_alias = true;
                        pixmap.fill_path(&path, &paint, FillRule::Winding, xform, None);
                    }
                }
            }
            VectorElement::CubicCurve {
                p0,
                p1,
                p2,
                p3,
                width,
                color,
            } => {
                let mut pb = PathBuilder::new();
                pb.move_to(p0[0] as f32, p0[1] as f32);
                pb.cubic_to(
                    p1[0] as f32,
                    p1[1] as f32,
                    p2[0] as f32,
                    p2[1] as f32,
                    p3[0] as f32,
                    p3[1] as f32,
                );
                if let Some(path) = pb.finish() {
                    let mut paint = Paint::default();
                    paint.set_color_rgba8(color.0, color.1, color.2, color.3);
                    paint.anti_alias = true;
                    let stroke = Stroke {
                        width: *width as f32,
                        line_cap: LineCap::Round,
                        ..Stroke::default()
                    };
                    pixmap.stroke_path(&path, &paint, &stroke, xform, None);
                }
            }
            VectorElement::Circle {
                cx,
                cy,
                r,
                fill,
                stroke,
            } => {
                if let Some(path) = PathBuilder::from_circle(*cx as f32, *cy as f32, *r as f32) {
                    if let Some(f) = fill {
                        let mut paint = Paint::default();
                        paint.set_color_rgba8(f.0, f.1, f.2, f.3);
                        paint.anti_alias = true;
                        pixmap.fill_path(&path, &paint, FillRule::Winding, xform, None);
                    }
                    if let Some((c, w)) = stroke {
                        let mut paint = Paint::default();
                        paint.set_color_rgba8(c.0, c.1, c.2, c.3);
                        paint.anti_alias = true;
                        let st = Stroke {
                            width: *w as f32,
                            ..Stroke::default()
                        };
                        pixmap.stroke_path(&path, &paint, &st, xform, None);
                    }
                }
            }
            VectorElement::Label {
                text,
                cx,
                cy,
                font_size,
                color,
                ..
            } => {
                draw_vector_glyph_label(pixmap, text, *cx, *cy, *font_size, *color, xform);
            }
        }
    }
}

/// Normalize common Latin-1 accented letters to `(ascii_base, diacritic_kind)` so
/// stroke titles remain legible in PNG export with zero system font dependencies.
fn decompose_latin_glyph(ch: char) -> (char, u8) {
    match ch {
        'á' | 'à' | 'â' | 'ã' | 'ä' => ('a', 1),
        'Á' | 'À' | 'Â' | 'Ã' | 'Ä' => ('A', 1),
        'é' | 'è' | 'ê' | 'ë' => ('e', 1),
        'É' | 'È' | 'Ê' | 'Ë' => ('E', 1),
        'í' | 'ì' | 'î' | 'ï' => ('i', 1),
        'Í' | 'Ì' | 'Î' | 'Ï' => ('I', 1),
        'ó' | 'ò' | 'ô' | 'õ' | 'ö' => ('o', 1),
        'Ó' | 'Ò' | 'Ô' | 'Õ' | 'Ö' => ('O', 1),
        'ú' | 'ù' | 'û' | 'ü' => ('u', 1),
        'Ú' | 'Ù' | 'Û' | 'Ü' => ('U', 1),
        'ñ' => ('n', 2),
        'Ñ' => ('N', 2),
        'ç' => ('c', 0),
        'Ç' => ('C', 0),
        other => (other, 0),
    }
}

/// Normalized `[0..1] x [0..1]` vector stroke polylines covering ASCII `A-Z`, `a-z`, `0-9`, and symbols.
fn glyph_polylines(ch: char) -> &'static [&'static [(f32, f32)]] {
    match ch {
        'A' => &[
            &[(0.0, 1.0), (0.5, 0.0), (1.0, 1.0)],
            &[(0.2, 0.65), (0.8, 0.65)],
        ],
        'B' => &[&[
            (0.0, 0.0),
            (0.0, 1.0),
            (0.8, 1.0),
            (1.0, 0.8),
            (1.0, 0.62),
            (0.75, 0.48),
            (0.0, 0.48),
            (0.75, 0.48),
            (0.95, 0.34),
            (0.95, 0.14),
            (0.75, 0.0),
            (0.0, 0.0),
        ]],
        'C' => &[&[
            (1.0, 0.18),
            (0.65, 0.0),
            (0.25, 0.0),
            (0.0, 0.25),
            (0.0, 0.75),
            (0.25, 1.0),
            (0.65, 1.0),
            (1.0, 0.82),
        ]],
        'D' => &[&[
            (0.0, 0.0),
            (0.0, 1.0),
            (0.65, 1.0),
            (1.0, 0.72),
            (1.0, 0.28),
            (0.65, 0.0),
            (0.0, 0.0),
        ]],
        'E' => &[
            &[(1.0, 0.0), (0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
            &[(0.0, 0.5), (0.8, 0.5)],
        ],
        'F' => &[
            &[(1.0, 0.0), (0.0, 0.0), (0.0, 1.0)],
            &[(0.0, 0.5), (0.8, 0.5)],
        ],
        'G' => &[&[
            (0.95, 0.18),
            (0.65, 0.0),
            (0.25, 0.0),
            (0.0, 0.25),
            (0.0, 0.75),
            (0.25, 1.0),
            (0.75, 1.0),
            (1.0, 0.75),
            (1.0, 0.52),
            (0.55, 0.52),
        ]],
        'H' => &[
            &[(0.0, 0.0), (0.0, 1.0)],
            &[(1.0, 0.0), (1.0, 1.0)],
            &[(0.0, 0.5), (1.0, 0.5)],
        ],
        'I' => &[
            &[(0.5, 0.0), (0.5, 1.0)],
            &[(0.15, 0.0), (0.85, 0.0)],
            &[(0.15, 1.0), (0.85, 1.0)],
        ],
        'J' => &[&[
            (0.85, 0.0),
            (0.85, 0.78),
            (0.6, 1.0),
            (0.25, 1.0),
            (0.05, 0.8),
        ]],
        'K' => &[
            &[(0.0, 0.0), (0.0, 1.0)],
            &[(0.95, 0.0), (0.0, 0.55), (1.0, 1.0)],
        ],
        'L' => &[&[(0.0, 0.0), (0.0, 1.0), (0.95, 1.0)]],
        'M' => &[&[(0.0, 1.0), (0.0, 0.0), (0.5, 0.6), (1.0, 0.0), (1.0, 1.0)]],
        'N' => &[&[(0.0, 1.0), (0.0, 0.0), (1.0, 1.0), (1.0, 0.0)]],
        'O' | '0' => &[&[
            (0.5, 0.0),
            (0.95, 0.25),
            (0.95, 0.75),
            (0.5, 1.0),
            (0.05, 0.75),
            (0.05, 0.25),
            (0.5, 0.0),
        ]],
        'P' => &[&[
            (0.0, 1.0),
            (0.0, 0.0),
            (0.8, 0.0),
            (1.0, 0.16),
            (1.0, 0.38),
            (0.8, 0.52),
            (0.0, 0.52),
        ]],
        'Q' => &[
            &[
                (0.5, 0.0),
                (0.95, 0.25),
                (0.95, 0.72),
                (0.5, 0.96),
                (0.05, 0.72),
                (0.05, 0.25),
                (0.5, 0.0),
            ],
            &[(0.6, 0.68), (1.0, 1.04)],
        ],
        'R' => &[
            &[
                (0.0, 1.0),
                (0.0, 0.0),
                (0.8, 0.0),
                (1.0, 0.16),
                (1.0, 0.36),
                (0.8, 0.5),
                (0.0, 0.5),
            ],
            &[(0.55, 0.5), (1.0, 1.0)],
        ],
        'S' => &[&[
            (0.95, 0.16),
            (0.7, 0.0),
            (0.25, 0.0),
            (0.05, 0.18),
            (0.05, 0.38),
            (0.25, 0.48),
            (0.75, 0.52),
            (0.95, 0.64),
            (0.95, 0.84),
            (0.7, 1.0),
            (0.25, 1.0),
            (0.05, 0.84),
        ]],
        'T' => &[&[(0.0, 0.0), (1.0, 0.0)], &[(0.5, 0.0), (0.5, 1.0)]],
        'U' => &[&[
            (0.0, 0.0),
            (0.0, 0.78),
            (0.25, 1.0),
            (0.75, 1.0),
            (1.0, 0.78),
            (1.0, 0.0),
        ]],
        'V' => &[&[(0.0, 0.0), (0.5, 1.0), (1.0, 0.0)]],
        'W' => &[&[
            (0.0, 0.0),
            (0.25, 1.0),
            (0.5, 0.45),
            (0.75, 1.0),
            (1.0, 0.0),
        ]],
        'X' => &[&[(0.0, 0.0), (1.0, 1.0)], &[(1.0, 0.0), (0.0, 1.0)]],
        'Y' => &[
            &[(0.0, 0.0), (0.5, 0.5), (1.0, 0.0)],
            &[(0.5, 0.5), (0.5, 1.0)],
        ],
        'Z' => &[&[(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]],
        'a' => &[
            &[(0.15, 0.42), (0.55, 0.34), (0.9, 0.48), (0.9, 1.0)],
            &[
                (0.9, 0.65),
                (0.25, 0.65),
                (0.05, 0.82),
                (0.25, 1.0),
                (0.9, 0.92),
            ],
        ],
        'b' => &[
            &[(0.1, 0.0), (0.1, 1.0)],
            &[
                (0.1, 0.48),
                (0.5, 0.34),
                (0.92, 0.5),
                (0.92, 0.84),
                (0.5, 1.0),
                (0.1, 0.88),
            ],
        ],
        'c' => &[&[
            (0.95, 0.46),
            (0.65, 0.34),
            (0.25, 0.34),
            (0.05, 0.56),
            (0.05, 0.8),
            (0.25, 1.0),
            (0.65, 1.0),
            (0.95, 0.88),
        ]],
        'd' => &[
            &[(0.9, 0.0), (0.9, 1.0)],
            &[
                (0.9, 0.48),
                (0.5, 0.34),
                (0.08, 0.5),
                (0.08, 0.84),
                (0.5, 1.0),
                (0.9, 0.88),
            ],
        ],
        'e' => &[&[
            (0.08, 0.66),
            (0.95, 0.66),
            (0.95, 0.5),
            (0.65, 0.34),
            (0.28, 0.34),
            (0.05, 0.54),
            (0.05, 0.82),
            (0.3, 1.0),
            (0.85, 0.94),
        ]],
        'f' => &[
            &[(0.85, 0.08), (0.6, 0.0), (0.4, 0.14), (0.4, 1.0)],
            &[(0.12, 0.42), (0.78, 0.42)],
        ],
        'g' => &[
            &[(0.9, 0.38), (0.9, 1.05), (0.65, 1.2), (0.2, 1.15)],
            &[
                (0.9, 0.48),
                (0.5, 0.34),
                (0.1, 0.5),
                (0.1, 0.82),
                (0.5, 0.96),
                (0.9, 0.84),
            ],
        ],
        'h' => &[
            &[(0.1, 0.0), (0.1, 1.0)],
            &[(0.1, 0.54), (0.45, 0.34), (0.85, 0.48), (0.85, 1.0)],
        ],
        'i' => &[&[(0.5, 0.36), (0.5, 1.0)], &[(0.5, 0.12), (0.5, 0.2)]],
        'j' => &[
            &[(0.65, 0.36), (0.65, 1.06), (0.42, 1.2), (0.15, 1.14)],
            &[(0.65, 0.12), (0.65, 0.2)],
        ],
        'k' => &[
            &[(0.15, 0.0), (0.15, 1.0)],
            &[(0.9, 0.36), (0.15, 0.72), (0.92, 1.0)],
        ],
        'l' => &[&[(0.45, 0.0), (0.45, 0.92), (0.72, 1.0)]],
        'm' => &[
            &[
                (0.02, 1.0),
                (0.02, 0.36),
                (0.32, 0.36),
                (0.5, 0.55),
                (0.5, 1.0),
            ],
            &[(0.5, 0.48), (0.78, 0.36), (0.98, 0.52), (0.98, 1.0)],
        ],
        'n' => &[
            &[(0.12, 1.0), (0.12, 0.34)],
            &[(0.12, 0.54), (0.45, 0.34), (0.88, 0.48), (0.88, 1.0)],
        ],
        'o' => &[&[
            (0.5, 0.34),
            (0.94, 0.52),
            (0.94, 0.82),
            (0.5, 1.0),
            (0.06, 0.82),
            (0.06, 0.52),
            (0.5, 0.34),
        ]],
        'p' => &[
            &[(0.12, 0.34), (0.12, 1.2)],
            &[
                (0.12, 0.46),
                (0.52, 0.34),
                (0.92, 0.52),
                (0.92, 0.82),
                (0.52, 0.96),
                (0.12, 0.84),
            ],
        ],
        'q' => &[
            &[(0.9, 0.34), (0.9, 1.2)],
            &[
                (0.9, 0.46),
                (0.48, 0.34),
                (0.08, 0.52),
                (0.08, 0.82),
                (0.48, 0.96),
                (0.9, 0.84),
            ],
        ],
        'r' => &[
            &[(0.2, 0.34), (0.2, 1.0)],
            &[(0.2, 0.54), (0.52, 0.36), (0.9, 0.44)],
        ],
        's' => &[&[
            (0.9, 0.42),
            (0.6, 0.34),
            (0.18, 0.38),
            (0.1, 0.54),
            (0.3, 0.64),
            (0.78, 0.7),
            (0.92, 0.84),
            (0.7, 1.0),
            (0.25, 1.0),
            (0.08, 0.9),
        ]],
        't' => &[
            &[(0.45, 0.08), (0.45, 0.9), (0.68, 1.0), (0.9, 0.94)],
            &[(0.12, 0.38), (0.84, 0.38)],
        ],
        'u' => &[
            &[(0.12, 0.36), (0.12, 0.84), (0.4, 1.0), (0.88, 0.84)],
            &[(0.88, 0.36), (0.88, 1.0)],
        ],
        'v' => &[&[(0.06, 0.36), (0.5, 1.0), (0.94, 0.36)]],
        'w' => &[&[
            (0.02, 0.36),
            (0.26, 1.0),
            (0.5, 0.58),
            (0.74, 1.0),
            (0.98, 0.36),
        ]],
        'x' => &[&[(0.1, 0.36), (0.9, 1.0)], &[(0.9, 0.36), (0.1, 1.0)]],
        'y' => &[
            &[(0.08, 0.36), (0.5, 0.94)],
            &[(0.92, 0.36), (0.42, 1.12), (0.16, 1.2)],
        ],
        'z' => &[&[(0.1, 0.36), (0.9, 0.36), (0.1, 1.0), (0.9, 1.0)]],
        '1' => &[
            &[(0.25, 0.2), (0.52, 0.0), (0.52, 1.0)],
            &[(0.2, 1.0), (0.85, 1.0)],
        ],
        '2' => &[&[
            (0.08, 0.22),
            (0.35, 0.0),
            (0.78, 0.0),
            (0.96, 0.2),
            (0.92, 0.44),
            (0.06, 1.0),
            (0.98, 1.0),
        ]],
        '3' => &[&[
            (0.1, 0.14),
            (0.4, 0.0),
            (0.8, 0.0),
            (0.96, 0.2),
            (0.85, 0.44),
            (0.4, 0.48),
            (0.85, 0.54),
            (0.98, 0.78),
            (0.75, 1.0),
            (0.3, 1.0),
            (0.08, 0.86),
        ]],
        '4' => &[&[(0.82, 1.0), (0.82, 0.0), (0.08, 0.68), (0.98, 0.68)]],
        '5' => &[&[
            (0.94, 0.0),
            (0.14, 0.0),
            (0.1, 0.44),
            (0.72, 0.44),
            (0.96, 0.64),
            (0.94, 0.86),
            (0.65, 1.0),
            (0.2, 1.0),
            (0.08, 0.86),
        ]],
        '6' => &[&[
            (0.86, 0.12),
            (0.58, 0.0),
            (0.22, 0.14),
            (0.06, 0.46),
            (0.06, 0.84),
            (0.34, 1.0),
            (0.78, 1.0),
            (0.96, 0.78),
            (0.88, 0.52),
            (0.48, 0.44),
            (0.08, 0.58),
        ]],
        '7' => &[&[(0.06, 0.0), (0.98, 0.0), (0.36, 1.0)]],
        '8' => &[&[
            (0.5, 0.48),
            (0.15, 0.34),
            (0.15, 0.12),
            (0.5, 0.0),
            (0.85, 0.12),
            (0.85, 0.34),
            (0.5, 0.48),
            (0.1, 0.64),
            (0.1, 0.88),
            (0.5, 1.0),
            (0.9, 0.88),
            (0.9, 0.64),
            (0.5, 0.48),
        ]],
        '9' => &[&[
            (0.14, 0.88),
            (0.44, 1.0),
            (0.82, 0.86),
            (0.94, 0.5),
            (0.94, 0.18),
            (0.66, 0.0),
            (0.22, 0.0),
            (0.06, 0.24),
            (0.14, 0.48),
            (0.56, 0.56),
            (0.94, 0.44),
        ]],
        '+' => &[&[(0.5, 0.2), (0.5, 0.8)], &[(0.15, 0.5), (0.85, 0.5)]],
        '-' => &[&[(0.15, 0.5), (0.85, 0.5)]],
        '_' => &[&[(0.05, 1.04), (0.95, 1.04)]],
        '.' => &[&[(0.48, 0.92), (0.54, 1.0)]],
        ',' => &[&[(0.55, 0.92), (0.36, 1.14)]],
        ':' => &[&[(0.5, 0.38), (0.52, 0.44)], &[(0.5, 0.88), (0.52, 0.96)]],
        '/' => &[&[(0.12, 1.0), (0.88, 0.0)]],
        '(' => &[&[(0.72, 0.0), (0.36, 0.26), (0.36, 0.74), (0.72, 1.0)]],
        ')' => &[&[(0.28, 0.0), (0.64, 0.26), (0.64, 0.74), (0.28, 1.0)]],
        '[' => &[&[(0.74, 0.0), (0.32, 0.0), (0.32, 1.0), (0.74, 1.0)]],
        ']' => &[&[(0.26, 0.0), (0.68, 0.0), (0.68, 1.0), (0.26, 1.0)]],
        '=' => &[&[(0.14, 0.36), (0.86, 0.36)], &[(0.14, 0.66), (0.86, 0.66)]],
        ' ' => &[],
        _ => &[
            &[(0.2, 0.15), (0.8, 0.15), (0.5, 0.58)],
            &[(0.5, 0.86), (0.5, 0.96)],
        ],
    }
}

/// Crisp geometric stroke-font renderer for PNG export (zero system font dependencies).
fn draw_vector_glyph_label(
    pixmap: &mut Pixmap,
    text: &str,
    cx: f64,
    cy: f64,
    font_size: f64,
    color: Rgba,
    xform: Transform,
) {
    let gw = font_size * 0.56;
    let gh = font_size * 0.78;
    let step_x = gw * 1.38;
    let char_count = text.chars().count();
    let total_w = if char_count == 0 {
        0.0
    } else {
        (char_count - 1) as f64 * step_x + gw
    };
    let mut x0 = cx - total_w * 0.5;
    let y0 = cy - gh * 0.5;

    let mut paint = Paint::default();
    paint.set_color_rgba8(color.0, color.1, color.2, color.3);
    paint.anti_alias = true;
    let stroke = Stroke {
        width: (font_size * 0.11).clamp(1.0, 2.0) as f32,
        line_cap: LineCap::Round,
        line_join: LineJoin::Round,
        ..Stroke::default()
    };

    for raw_ch in text.chars() {
        let (base_ch, diacritic) = decompose_latin_glyph(raw_ch);
        let mut pb = PathBuilder::new();
        for poly in glyph_polylines(base_ch) {
            if let Some(&(px0, py0)) = poly.first() {
                pb.move_to(
                    (x0 + gw * f64::from(px0)) as f32,
                    (y0 + gh * f64::from(py0)) as f32,
                );
                for &(px, py) in &poly[1..] {
                    pb.line_to(
                        (x0 + gw * f64::from(px)) as f32,
                        (y0 + gh * f64::from(py)) as f32,
                    );
                }
            }
        }
        if diacritic == 1 {
            pb.move_to((x0 + gw * 0.35) as f32, (y0 - gh * 0.06) as f32);
            pb.line_to((x0 + gw * 0.7) as f32, (y0 - gh * 0.22) as f32);
        } else if diacritic == 2 {
            pb.move_to((x0 + gw * 0.18) as f32, (y0 - gh * 0.08) as f32);
            pb.line_to((x0 + gw * 0.45) as f32, (y0 - gh * 0.18) as f32);
            pb.line_to((x0 + gw * 0.65) as f32, (y0 - gh * 0.08) as f32);
            pb.line_to((x0 + gw * 0.85) as f32, (y0 - gh * 0.18) as f32);
        }
        if let Some(path) = pb.finish() {
            pixmap.stroke_path(&path, &paint, &stroke, xform, None);
        }
        x0 += step_x;
    }
}

fn hex_color(c: Rgba) -> String {
    format!("#{:02x}{:02x}{:02x}", c.0, c.1, c.2)
}

fn escape_xml(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        if c.is_control() || matches!(c, '\u{FFFE}' | '\u{FFFF}') {
            out.push(' ');
            continue;
        }
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

fn escape_pdf_ps_text(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '(' => out.push_str("\\("),
            ')' => out.push_str("\\)"),
            '\\' => out.push_str("\\\\"),
            c if c.is_ascii() && !c.is_control() => out.push(c),
            '\u{00A0}'..='\u{00FF}' => {
                let code = c as u32;
                let _ = write!(out, "\\{code:03o}");
            }
            _ => out.push('?'),
        }
    }
    out
}

fn truncate_str(s: &str, max_chars: usize) -> String {
    let mut chars = s.chars();
    let prefix: String = chars.by_ref().take(max_chars).collect();
    if chars.next().is_some() {
        format!("{prefix}...")
    } else {
        prefix
    }
}
