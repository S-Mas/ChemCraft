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
                    x1,
                    h - y1,
                    x2,
                    h - y2
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
                        first[0],
                        h - first[1]
                    );
                    for p in &points[1..] {
                        let _ = writeln!(out, "{:.2} {:.2} lineto", p[0], h - p[1]);
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
                        first[0],
                        h - first[1]
                    );
                    for p in &points[1..] {
                        let _ = writeln!(out, "{:.2} {:.2} lineto", p[0], h - p[1]);
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
                    p0[0],
                    h - p0[1],
                    p1[0],
                    h - p1[1],
                    p2[0],
                    h - p2[1],
                    p3[0],
                    h - p3[1]
                );
            }
            VectorElement::Circle {
                cx, cy, r, stroke, ..
            } => {
                if let Some((c, w)) = stroke {
                    let _ = writeln!(
                        out,
                        "{:.3} {:.3} {:.3} setrgbcolor {w:.2} setlinewidth newpath {:.2} {:.2} {r:.2} 0 360 arc stroke",
                        f64::from(c.0) / 255.0,
                        f64::from(c.1) / 255.0,
                        f64::from(c.2) / 255.0,
                        cx,
                        h - cy
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
                    cx,
                    h - cy,
                    -font_size * 0.32
                );
            }
        }
    }
    let _ = writeln!(out, "showpage\n%%EOF");
    out
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
            VectorElement::Circle { .. } => {}
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
                draw_vector_glyph_label(&mut pixmap, text, *cx, *cy, *font_size, *color, xform);
            }
        }
    }

    let mut png_bytes = Vec::new();
    {
        let mut encoder = png::Encoder::new(&mut png_bytes, px_w, px_h);
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
    let total_w = text.len() as f64 * (gw * 1.15);
    let mut x0 = cx - total_w * 0.5;
    let y0 = cy - gh * 0.5;

    let mut paint = Paint::default();
    paint.set_color_rgba8(color.0, color.1, color.2, color.3);
    paint.anti_alias = true;
    let stroke = Stroke {
        width: (font_size * 0.13).clamp(1.2, 2.4) as f32,
        line_cap: LineCap::Round,
        line_join: LineJoin::Round,
        ..Stroke::default()
    };

    for ch in text.chars() {
        let mut pb = PathBuilder::new();
        let l = x0 as f32;
        let r = (x0 + gw) as f32;
        let m = (x0 + gw * 0.5) as f32;
        let t = y0 as f32;
        let b = (y0 + gh) as f32;
        let mid_y = (y0 + gh * 0.5) as f32;

        match ch {
            'N' => {
                pb.move_to(l, b);
                pb.line_to(l, t);
                pb.line_to(r, b);
                pb.line_to(r, t);
            }
            'O' | '0' => {
                pb.move_to(m, t);
                pb.line_to(r, t + gh as f32 * 0.25);
                pb.line_to(r, b - gh as f32 * 0.25);
                pb.line_to(m, b);
                pb.line_to(l, b - gh as f32 * 0.25);
                pb.line_to(l, t + gh as f32 * 0.25);
                pb.close();
            }
            'H' => {
                pb.move_to(l, t);
                pb.line_to(l, b);
                pb.move_to(r, t);
                pb.line_to(r, b);
                pb.move_to(l, mid_y);
                pb.line_to(r, mid_y);
            }
            'C' => {
                pb.move_to(r, t + 1.5);
                pb.line_to(m, t);
                pb.line_to(l, t + 2.5);
                pb.line_to(l, b - 2.5);
                pb.line_to(m, b);
                pb.line_to(r, b - 1.5);
            }
            'S' => {
                pb.move_to(r, t + 1.5);
                pb.line_to(l, t + 2.0);
                pb.line_to(l, mid_y);
                pb.line_to(r, mid_y);
                pb.line_to(r, b - 2.0);
                pb.line_to(l, b - 1.5);
            }
            'P' => {
                pb.move_to(l, b);
                pb.line_to(l, t);
                pb.line_to(r, t);
                pb.line_to(r, mid_y);
                pb.line_to(l, mid_y);
            }
            'F' => {
                pb.move_to(r, t);
                pb.line_to(l, t);
                pb.line_to(l, b);
                pb.move_to(l, mid_y);
                pb.line_to(m + 1.0, mid_y);
            }
            'I' => {
                pb.move_to(m, t);
                pb.line_to(m, b);
                pb.move_to(l, t);
                pb.line_to(r, t);
                pb.move_to(l, b);
                pb.line_to(r, b);
            }
            'B' => {
                pb.move_to(l, t);
                pb.line_to(l, b);
                pb.line_to(r, b);
                pb.line_to(r, mid_y);
                pb.line_to(l, mid_y);
                pb.line_to(r, mid_y);
                pb.line_to(r, t);
                pb.close();
            }
            'l' => {
                pb.move_to(m, t);
                pb.line_to(m, b);
            }
            'r' => {
                pb.move_to(l, mid_y);
                pb.line_to(l, b);
                pb.move_to(l, mid_y + 1.0);
                pb.line_to(r, mid_y);
            }
            '+' => {
                pb.move_to(m, t + 2.0);
                pb.line_to(m, b - 2.0);
                pb.move_to(l, mid_y);
                pb.line_to(r, mid_y);
            }
            '-' => {
                pb.move_to(l, mid_y);
                pb.line_to(r, mid_y);
            }
            '2' => {
                pb.move_to(l, t + 1.5);
                pb.line_to(r, t + 1.5);
                pb.line_to(r, mid_y);
                pb.line_to(l, mid_y);
                pb.line_to(l, b);
                pb.line_to(r, b);
            }
            '3' => {
                pb.move_to(l, t);
                pb.line_to(r, t);
                pb.line_to(r, b);
                pb.line_to(l, b);
                pb.move_to(l, mid_y);
                pb.line_to(r, mid_y);
            }
            _ => {
                pb.move_to(l, mid_y);
                pb.line_to(r, mid_y);
            }
        }
        if let Some(path) = pb.finish() {
            pixmap.stroke_path(&path, &paint, &stroke, xform, None);
        }
        x0 += gw * 1.15;
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
