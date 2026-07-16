//! Basis-set + molecule representation and an NWChem-format parser.
//!
//! Pure Rust (shared). The libcint array builder (atm/bas/env) lives in
//! crates/cint-sys; this module just holds the data and parses basis text so
//! ANY molecule/basis can be fed to the engine (needed for the 27 emitters at
//! def2-SVP).

use std::collections::HashMap;

pub const ANG_TO_BOHR: f64 = 1.8897259886;

/// One contracted shell of a single angular momentum.
#[derive(Clone, Debug)]
pub struct Shell {
    pub l: i32,
    pub exps: Vec<f64>,
    pub coeffs: Vec<f64>,
}

#[derive(Clone, Debug, Default)]
pub struct ElementBasis {
    pub shells: Vec<Shell>,
}

/// Basis set keyed by atomic number Z.
#[derive(Clone, Debug, Default)]
pub struct BasisSet {
    pub elements: HashMap<i32, ElementBasis>,
}

#[derive(Clone, Debug)]
pub struct Atom {
    pub z: i32,
    pub coord: [f64; 3], // bohr
}

#[derive(Clone, Debug, Default)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
}

impl Molecule {
    pub fn n_electrons(&self) -> i32 {
        self.atoms.iter().map(|a| a.z).sum()
    }

    /// Parse a standard `.xyz` file (Ångström) into a Molecule (bohr).
    /// Format: line1 = atom count, line2 = comment, then `El x y z` rows.
    pub fn from_xyz_str(text: &str) -> Result<Molecule, String> {
        let mut lines = text.lines();
        let n: usize = lines
            .next()
            .ok_or("empty xyz")?
            .trim()
            .parse()
            .map_err(|_| "xyz: bad atom count on line 1")?;
        let _comment = lines.next();
        let mut atoms = Vec::with_capacity(n);
        for (i, l) in lines.enumerate().take(n) {
            let t: Vec<&str> = l.split_whitespace().collect();
            if t.len() < 4 {
                return Err(format!("xyz line {}: expected 'El x y z'", i + 3));
            }
            let z = symbol_to_z(t[0]).ok_or_else(|| format!("unknown element '{}'", t[0]))?;
            let parse = |s: &str| s.parse::<f64>().map_err(|_| format!("xyz line {}: bad number", i + 3));
            atoms.push(Atom {
                z,
                coord: [
                    parse(t[1])? * ANG_TO_BOHR,
                    parse(t[2])? * ANG_TO_BOHR,
                    parse(t[3])? * ANG_TO_BOHR,
                ],
            });
        }
        if atoms.len() != n {
            return Err(format!("xyz: header says {n} atoms, found {}", atoms.len()));
        }
        Ok(Molecule { atoms })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_xyz_water() {
        let xyz = "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.0 0.9584\nH 0.9276 0.0 -0.2400\n";
        let m = Molecule::from_xyz_str(xyz).unwrap();
        assert_eq!(m.atoms.len(), 3);
        assert_eq!(m.atoms[0].z, 8);
        assert_eq!(m.atoms[1].z, 1);
        assert_eq!(m.n_electrons(), 10);
        // 0.9584 Å -> bohr
        assert!((m.atoms[1].coord[2] - 0.9584 * ANG_TO_BOHR).abs() < 1e-9);
    }

    #[test]
    fn parse_sto3g_basis() {
        let txt = "O S\n130.7 0.15\n23.8 0.53\nO SP\n5.0 -0.1 0.15\n1.2 0.4 0.6\n";
        let bs = BasisSet::parse_nwchem(txt).unwrap();
        let o = bs.elements.get(&8).unwrap();
        // S + SP -> 3 shells (one S, one S from SP, one P from SP)
        assert_eq!(o.shells.len(), 3);
        assert_eq!(o.shells[0].l, 0);
        assert_eq!(o.shells[1].l, 0);
        assert_eq!(o.shells[2].l, 1);
    }
}

pub fn symbol_to_z(sym: &str) -> Option<i32> {
    const T: &[&str] = &[
        "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P",
        "S", "Cl", "Ar", "K", "Ca",
    ];
    let s = sym.trim();
    T.iter().position(|&e| e.eq_ignore_ascii_case(s)).map(|i| i as i32)
}

fn ang_letters(label: &str) -> Vec<i32> {
    // "S"->[0], "P"->[1], "SP"->[0,1], "D"->[2], "F"->[3] ...
    label
        .chars()
        .filter_map(|c| match c.to_ascii_uppercase() {
            'S' => Some(0),
            'P' => Some(1),
            'D' => Some(2),
            'F' => Some(3),
            'G' => Some(4),
            _ => None,
        })
        .collect()
}

impl BasisSet {
    /// Parse NWChem/BSE-style basis text. Lines are either a header
    /// (`<Element> <AngLabel>`) or primitive rows (`exp coeff [coeff2 ...]`).
    /// Comments (`#`), `BASIS ...`, `END`, blanks are ignored.
    pub fn parse_nwchem(text: &str) -> Result<BasisSet, String> {
        let mut bs = BasisSet::default();
        let mut cur_z: Option<i32> = None;
        let mut cur_ls: Vec<i32> = Vec::new();
        let mut rows: Vec<Vec<f64>> = Vec::new();

        let flush = |bs: &mut BasisSet, z: Option<i32>, ls: &[i32], rows: &[Vec<f64>]| {
            if let Some(z) = z {
                if ls.is_empty() || rows.is_empty() {
                    return;
                }
                let eb = bs.elements.entry(z).or_default();
                for (ci, &l) in ls.iter().enumerate() {
                    let exps: Vec<f64> = rows.iter().map(|r| r[0]).collect();
                    let coeffs: Vec<f64> = rows.iter().map(|r| r[1 + ci]).collect();
                    eb.shells.push(Shell { l, exps, coeffs });
                }
            }
        };

        for raw in text.lines() {
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let up = line.to_ascii_uppercase();
            if up.starts_with("BASIS") || up == "END" || up.starts_with("****") {
                continue;
            }
            let toks: Vec<&str> = line.split_whitespace().collect();
            // Header line: <symbol> <anglabel>, first token not a number.
            if toks.len() == 2 && toks[0].parse::<f64>().is_err() {
                if let Some(z) = symbol_to_z(toks[0]) {
                    // flush previous block
                    flush(&mut bs, cur_z, &cur_ls, &rows);
                    rows.clear();
                    cur_z = Some(z);
                    cur_ls = ang_letters(toks[1]);
                    continue;
                }
            }
            // Data row of numbers.
            let nums: Result<Vec<f64>, _> = toks.iter().map(|t| t.replace('D', "E").parse::<f64>()).collect();
            match nums {
                Ok(v) if v.len() >= 2 => rows.push(v),
                _ => {} // ignore unrecognized lines
            }
        }
        flush(&mut bs, cur_z, &cur_ls, &rows);
        if bs.elements.is_empty() {
            return Err("no basis shells parsed".into());
        }
        Ok(bs)
    }
}
