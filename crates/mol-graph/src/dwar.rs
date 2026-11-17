use std::collections::BTreeMap;
use std::fmt::Write;

use crate::error::MolError;
use crate::mol::MolGraph;
use crate::smiles::parse_smiles;

/// Parsed compound record from `.dwar` or `.csv`.
#[derive(Debug, Clone)]
pub struct TabularMolRecord {
    pub name: String,
    pub mol: MolGraph,
    pub columns: BTreeMap<String, String>,
}

/// Clean-room DataWarrior `.dwar` parser (skips optional `<datawarrior-fileinfo>` /
/// `<column properties>` metadata headers, locates the SMILES & Name columns, and parses
/// every row into a `MolGraph` plus numeric/categorical SAR metadata).
pub fn parse_dwar(input: &str) -> Result<Vec<Result<TabularMolRecord, MolError>>, MolError> {
    if input.trim().is_empty() {
        return Err(MolError::EmptyInput);
    }
    let normalized = input.replace("\r\n", "\n");
    let lines = normalized.lines().enumerate();
    let mut in_block: Option<&'static str> = None;
    let mut header_cols: Option<Vec<String>> = None;
    let mut records = Vec::new();

    for (ln, raw) in lines {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == "<datawarrior-fileinfo>" {
            in_block = Some("</datawarrior-fileinfo>");
            continue;
        }
        if trimmed == "<column properties>" {
            in_block = Some("</column properties>");
            continue;
        }
        if let Some(end_tag) = in_block {
            if trimmed == end_tag {
                in_block = None;
            }
            continue;
        }
        if trimmed.starts_with('<') && trimmed.ends_with('>') {
            // Stop before embedded hitlists/views trailers if present
            break;
        }

        let fields: Vec<String> = raw.split('\t').map(|s| s.trim().to_string()).collect();
        if header_cols.is_none() {
            header_cols = Some(fields);
            continue;
        }

        let Some(ref headers) = header_cols else {
            return Err(MolError::InvalidDwar {
                line: ln + 1,
                msg: "missing header row".to_string(),
            });
        };

        let smi_idx = headers
            .iter()
            .position(|h| {
                let norm = h.to_ascii_lowercase();
                norm == "smiles"
                    || norm == "canonical smiles"
                    || norm == "canonical_smiles"
                    || norm.starts_with("structure of ")
                    || norm == "structure"
            })
            .unwrap_or(0);
        let name_idx = headers.iter().position(|h| {
            let norm = h.to_ascii_lowercase();
            matches!(
                norm.as_str(),
                "id" | "name" | "molecule" | "compound" | "compound_id" | "molecule name"
            )
        });

        let smi = fields.get(smi_idx).map(String::as_str).unwrap_or("").trim();
        if smi.is_empty() {
            continue;
        }
        let default_name = format!("Row_{}", records.len() + 1);
        let name = name_idx
            .and_then(|idx| fields.get(idx))
            .map(String::as_str)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .unwrap_or(&default_name)
            .to_string();

        let mut cols = BTreeMap::new();
        for (idx, h) in headers.iter().enumerate() {
            if idx == smi_idx {
                continue;
            }
            if let Some(val) = fields.get(idx)
                && !val.is_empty()
            {
                cols.insert(h.clone(), val.clone());
            }
        }

        let res = parse_smiles(smi).map(|mut mol| {
            mol.name = Some(name.clone());
            mol.properties = cols.clone();
            TabularMolRecord {
                name,
                mol,
                columns: cols,
            }
        });
        records.push(res);
    }

    if header_cols.is_none() {
        return Err(MolError::InvalidDwar {
            line: 1,
            msg: "no tabular header found in .dwar input".to_string(),
        });
    }
    Ok(records)
}

/// Export molecules and SAR columns in clean-room `.dwar` format.
pub fn write_dwar(rows: &[(&MolGraph, &BTreeMap<String, String>)]) -> String {
    let mut all_cols: Vec<String> = Vec::new();
    for (_, props) in rows {
        for k in props.keys() {
            if k != "SMILES" && k != "Name" && !all_cols.contains(k) {
                all_cols.push(k.clone());
            }
        }
    }

    let mut out = String::new();
    let _ = writeln!(out, "<datawarrior-fileinfo>");
    let _ = writeln!(out, "<version=\"3.3\">");
    let _ = writeln!(out, "<rowcount=\"{}\">", rows.len());
    let _ = writeln!(out, "</datawarrior-fileinfo>");
    let _ = write!(out, "Name\tSMILES");
    for c in &all_cols {
        let _ = write!(out, "\t{}", c.replace(['\t', '\n'], " "));
    }
    let _ = writeln!(out);

    for (idx, (mol, props)) in rows.iter().enumerate() {
        let name = mol
            .name
            .clone()
            .unwrap_or_else(|| format!("Compound_{}", idx + 1));
        let smi = mol.canonical_smiles();
        let _ = write!(
            out,
            "{}\t{}",
            name.replace(['\t', '\n'], " "),
            smi.replace(['\t', '\n'], " ")
        );
        for c in &all_cols {
            let v = props
                .get(c)
                .map(|s| s.replace(['\t', '\n'], " "))
                .unwrap_or_default();
            let _ = write!(out, "\t{v}");
        }
        let _ = writeln!(out);
    }
    out
}

/// Parse RFC-4180 CSV into `TabularMolRecord`s using `smiles`/`Canonical SMILES` and SAR columns.
pub fn parse_csv_records(input: &str) -> Result<Vec<Result<TabularMolRecord, MolError>>, MolError> {
    if input.trim().is_empty() {
        return Err(MolError::EmptyInput);
    }
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .flexible(true)
        .from_reader(input.as_bytes());
    let headers = rdr
        .headers()
        .map_err(|e| MolError::InvalidCsv(e.to_string()))?
        .clone();

    let smi_col = headers
        .iter()
        .position(|c| {
            let norm = c.trim().to_ascii_lowercase().replace(' ', "_");
            matches!(norm.as_str(), "smiles" | "canonical_smiles" | "smi")
        })
        .unwrap_or(0);
    let name_col = headers.iter().position(|c| {
        let norm = c.trim().to_ascii_lowercase().replace(' ', "_");
        matches!(
            norm.as_str(),
            "name" | "molecule" | "molecule_name" | "id" | "compound"
        )
    });

    let mut out = Vec::new();
    for (idx, rec_res) in rdr.records().enumerate() {
        let rec = match rec_res {
            Ok(r) => r,
            Err(e) => {
                out.push(Err(MolError::InvalidCsv(format!("row {}: {e}", idx + 2))));
                continue;
            }
        };
        let Some(smi) = rec.get(smi_col).map(str::trim).filter(|s| !s.is_empty()) else {
            continue;
        };
        let default_name = format!("Molecule_{}", idx + 1);
        let name = name_col
            .and_then(|c| rec.get(c))
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .unwrap_or(&default_name)
            .to_string();
        let mut cols = BTreeMap::new();
        for (ci, h) in headers.iter().enumerate() {
            if ci == smi_col {
                continue;
            }
            if let Some(val) = rec.get(ci).map(str::trim).filter(|s| !s.is_empty()) {
                cols.insert(h.trim().to_string(), val.to_string());
            }
        }
        let res = parse_smiles(smi).map(|mut mol| {
            mol.name = Some(name.clone());
            mol.properties = cols.clone();
            TabularMolRecord {
                name,
                mol,
                columns: cols,
            }
        });
        out.push(res);
    }
    Ok(out)
}
