//! Minimal, dependency-free CSV reader.
//!
//! Handles double-quoted fields containing commas and escaped quotes ("").
//! Assumes records do not span multiple physical lines (true for our
//! reference set). Deterministic, no allocation surprises.

pub struct Table {
    pub header: Vec<String>,
    pub rows: Vec<Vec<String>>,
}

impl Table {
    /// Index of a column by header name, if present.
    pub fn col(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|h| h == name)
    }
}

/// Parse one CSV line into fields, honoring quotes.
fn parse_line(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        match c {
            '"' => {
                if in_quotes {
                    // Escaped quote inside a quoted field?
                    if chars.peek() == Some(&'"') {
                        cur.push('"');
                        chars.next();
                    } else {
                        in_quotes = false;
                    }
                } else {
                    in_quotes = true;
                }
            }
            ',' if !in_quotes => {
                fields.push(std::mem::take(&mut cur));
            }
            _ => cur.push(c),
        }
    }
    fields.push(cur);
    fields
}

/// Parse a CSV string into a `Table`. Empty/whitespace-only lines are skipped.
pub fn parse(text: &str) -> Result<Table, String> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let header = match lines.next() {
        Some(h) => parse_line(h),
        None => return Err("empty CSV".to_string()),
    };
    let ncol = header.len();
    let mut rows = Vec::new();
    for (i, line) in lines.enumerate() {
        let rec = parse_line(line);
        if rec.len() != ncol {
            return Err(format!(
                "row {}: expected {} fields, got {}",
                i + 2,
                ncol,
                rec.len()
            ));
        }
        rows.push(rec);
    }
    Ok(Table { header, rows })
}
