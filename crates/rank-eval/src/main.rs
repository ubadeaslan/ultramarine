//! M0 ranking-validation harness for the deep-blue TADF screening engine.
//!
//! Primary question this answers: *does a predictor rank known emitters'
//! ΔEST like experiment does?* Ranking accuracy (Spearman / Kendall tau-b)
//! is the headline metric, evaluated SEPARATELY for the two physics buckets
//! (donor-acceptor vs MR-TADF), because mixing them is physically meaningless.
//!
//! There is no quantum-chemistry engine yet. So this binary:
//!   * loads the experimental reference set,
//!   * proves the metric machinery with a deterministic self-test, and
//!   * if given `--pred <file>` (lines: `id,value`), scores those predictions
//!     against the experimental ΔEST.
//!
//! Usage:
//!   cargo run -p rank-eval -- [data/reference_emitters.csv] [--pred preds.csv]

mod csvmini;
mod stats;

use std::collections::HashMap;
use std::fs;
use std::process::ExitCode;

/// eV per nm conversion constant (E[eV] = HC_EV_NM / lambda[nm]).
const HC_EV_NM: f64 = 1239.841984;

/// Experimental emission-energy reference uncertainty (eV), from the
/// Bonn/Herbert STGABS27 work. Pairs closer than this are unrankable.
const EMISSION_TIE_EV: f64 = 0.2;

#[derive(Clone)]
struct Emitter {
    id: String,
    name: String,
    bucket: Bucket,
    delta_est: Option<f64>,
    lambda_nm: Option<f64>,
    verification: String,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Bucket {
    DonorAcceptor,
    MultiResonance,
}

impl Bucket {
    fn label(self) -> &'static str {
        match self {
            Bucket::DonorAcceptor => "D-A",
            Bucket::MultiResonance => "MR-TADF",
        }
    }
}

fn parse_opt_f64(s: &str) -> Option<f64> {
    let t = s.trim();
    if t.is_empty() {
        None
    } else {
        t.parse::<f64>().ok()
    }
}

fn load_emitters(path: &str) -> Result<Vec<Emitter>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
    let table = csvmini::parse(&text)?;
    let c_id = table.col("id").ok_or("missing column: id")?;
    let c_name = table.col("name").ok_or("missing column: name")?;
    let c_class = table.col("class").ok_or("missing column: class")?;
    let c_dest = table.col("deltaEST_eV").ok_or("missing column: deltaEST_eV")?;
    let c_lam = table.col("lambda_em_nm").ok_or("missing column: lambda_em_nm")?;
    let c_ver = table.col("verification").ok_or("missing column: verification")?;

    let mut out = Vec::new();
    for row in &table.rows {
        let class = &row[c_class];
        let bucket = if class.contains("MR") {
            Bucket::MultiResonance
        } else {
            Bucket::DonorAcceptor
        };
        out.push(Emitter {
            id: row[c_id].clone(),
            name: row[c_name].clone(),
            bucket,
            delta_est: parse_opt_f64(&row[c_dest]),
            lambda_nm: parse_opt_f64(&row[c_lam]),
            verification: row[c_ver].clone(),
        });
    }
    Ok(out)
}

/// Deterministic shuffle (LCG, fixed seed) for the self-test only.
fn deterministic_shuffle<T>(v: &mut [T]) {
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    for i in (1..v.len()).rev() {
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        let j = (state >> 33) as usize % (i + 1);
        v.swap(i, j);
    }
}

fn print_bucket_report(bucket: Bucket, emitters: &[Emitter]) {
    let members: Vec<&Emitter> = emitters.iter().filter(|e| e.bucket == bucket).collect();
    let with_dest: Vec<&Emitter> = members.iter().copied().filter(|e| e.delta_est.is_some()).collect();
    let with_lam: Vec<&Emitter> = members.iter().copied().filter(|e| e.lambda_nm.is_some()).collect();

    println!("\n=== Bucket: {} ===", bucket.label());
    println!("  members: {}", members.len());
    println!("  with experimental ΔEST: {}", with_dest.len());
    println!("  with experimental λ_em: {}", with_lam.len());
    for e in &with_dest {
        println!(
            "    - {:<8} {:<10} ΔEST={:.3} eV  [{}]",
            e.id,
            e.name,
            e.delta_est.unwrap(),
            e.verification
        );
    }
    if with_dest.len() < 4 {
        println!(
            "  ⚠  n={} too small for a meaningful rank correlation (need >=4). \
             Awaiting STGABS27 import.",
            with_dest.len()
        );
    }
}

/// Self-test sized to the bucket's usable ΔEST count. Uses a clean synthetic
/// monotone reference (1..=n) so the metric unambiguously returns +1 for a
/// perfect predictor, -1 for a reversed one, and a deterministic value for a
/// scramble. Proves the plumbing independent of real-data ties.
fn self_test(bucket: Bucket, emitters: &[Emitter]) {
    let n = emitters
        .iter()
        .filter(|e| e.bucket == bucket)
        .filter(|e| e.delta_est.is_some())
        .count();
    if n < 2 {
        return;
    }
    let exp: Vec<f64> = (1..=n).map(|v| v as f64).collect();
    let perfect = exp.clone();
    let reversed: Vec<f64> = exp.iter().rev().cloned().collect();
    let mut scrambled = exp.clone();
    deterministic_shuffle(&mut scrambled);

    println!("\n  self-test (synthetic, {n} pts):");
    for (tag, pred) in [("perfect", &perfect), ("reversed", &reversed), ("scrambled", &scrambled)] {
        let rho = stats::spearman(&exp, pred);
        let tau = stats::kendall_tau_b(&exp, pred, 0.0);
        println!(
            "    {:<10} spearman={:>7}  kendall_tau_b={:>7}",
            tag,
            fmt_opt(rho),
            fmt_opt(tau)
        );
    }
}

fn fmt_opt(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x:+.3}"),
        None => "  n/a".to_string(),
    }
}

/// Score an external predictions file (`id,value`) against experimental ΔEST,
/// per bucket. This is the real evaluation path once an engine exists.
fn score_predictions(emitters: &[Emitter], pred_path: &str) -> Result<(), String> {
    let text = fs::read_to_string(pred_path).map_err(|e| format!("read {pred_path}: {e}"))?;
    let mut preds: HashMap<String, f64> = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("id,") {
            continue;
        }
        let (id, val) = line
            .split_once(',')
            .ok_or_else(|| format!("{pred_path}:{}: expected `id,value`", i + 1))?;
        let v = val
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("{pred_path}:{}: bad number `{val}`", i + 1))?;
        preds.insert(id.trim().to_string(), v);
    }

    println!("\n########## SCORING PREDICTIONS: {pred_path} ##########");
    for bucket in [Bucket::DonorAcceptor, Bucket::MultiResonance] {
        let mut exp = Vec::new();
        let mut prd = Vec::new();
        for e in emitters.iter().filter(|e| e.bucket == bucket) {
            if let (Some(d), Some(p)) = (e.delta_est, preds.get(&e.id)) {
                exp.push(d);
                prd.push(*p);
            }
        }
        println!("\n  [{}] matched pairs: {}", bucket.label(), exp.len());
        if exp.len() < 2 {
            println!("    not enough matched data to score.");
            continue;
        }
        println!("    ΔEST spearman      = {}", fmt_opt(stats::spearman(&exp, &prd)));
        println!("    ΔEST kendall_tau_b = {}", fmt_opt(stats::kendall_tau_b(&exp, &prd, 0.0)));
    }
    Ok(())
}

// --- STGABS27 ranking benchmark -------------------------------------------

struct StgRow {
    exp_rank: f64,
    name: String,
    key: String,
    bucket: Bucket,
    outlier: bool,
}

fn load_stgabs27(path: &str) -> Result<Vec<StgRow>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
    let table = csvmini::parse(&text)?;
    let c_rank = table.col("exp_rank").ok_or("missing column: exp_rank")?;
    let c_name = table.col("name").ok_or("missing column: name")?;
    let c_key = table.col("key").ok_or("missing column: key")?;
    let c_bucket = table.col("bucket").ok_or("missing column: bucket")?;
    let c_out = table.col("tddft_outlier").ok_or("missing column: tddft_outlier")?;
    let mut rows = Vec::new();
    for r in &table.rows {
        rows.push(StgRow {
            exp_rank: r[c_rank].trim().parse::<f64>().map_err(|_| "bad exp_rank")?,
            name: r[c_name].clone(),
            key: r[c_key].clone(),
            bucket: if r[c_bucket].contains("MR") {
                Bucket::MultiResonance
            } else {
                Bucket::DonorAcceptor
            },
            outlier: r[c_out].trim().eq_ignore_ascii_case("Y"),
        });
    }
    Ok(rows)
}

/// Read a name-keyed predictions file (`name,value`).
fn load_named_preds(path: &str) -> Result<HashMap<String, f64>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("read {path}: {e}"))?;
    let mut m = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("name,") {
            continue;
        }
        let (k, v) = line
            .split_once(',')
            .ok_or_else(|| format!("{path}:{}: expected `name,value`", i + 1))?;
        let val = v
            .trim()
            .parse::<f64>()
            .map_err(|_| format!("{path}:{}: bad number `{v}`", i + 1))?;
        m.insert(k.trim().to_string(), val);
    }
    Ok(m)
}

fn corr_line(label: &str, refs: &[f64], pred: &[f64]) {
    println!(
        "    {:<28} n={:<3} spearman={}  kendall_tau_b={}",
        label,
        refs.len(),
        fmt_opt(stats::spearman(refs, pred)),
        fmt_opt(stats::kendall_tau_b(refs, pred, 0.0))
    );
}

/// Summarize STGABS27 and (optionally) score predictions against its
/// experimental rank order. Because Spearman/Kendall operate on ranks,
/// correlating predicted ΔEST against the experimental *rank* is exactly
/// equivalent to correlating against the (here unavailable) absolute eV
/// values.
fn handle_stgabs27(path: &str, pred_path: Option<&str>) -> Result<(), String> {
    let rows = load_stgabs27(path)?;
    let da = rows.iter().filter(|r| r.bucket == Bucket::DonorAcceptor).count();
    let mr = rows.iter().filter(|r| r.bucket == Bucket::MultiResonance).count();
    let out = rows.iter().filter(|r| r.outlier).count();
    println!("\n########## STGABS27 ranking benchmark ({path}) ##########");
    println!("  molecules: {}  (D-A: {da}, MR-TADF: {mr})", rows.len());
    println!("  experimental reference = ascending ΔEST rank order (1..n)");
    println!("  known TDDFT outliers (method failures): {out}");

    match pred_path {
        None => {
            // self-test on the rank order
            let refs: Vec<f64> = rows.iter().map(|r| r.exp_rank).collect();
            let mut scrambled = refs.clone();
            deterministic_shuffle(&mut scrambled);
            println!("  self-test:");
            corr_line("perfect (pred=rank)", &refs, &refs.clone());
            corr_line("scrambled", &refs, &scrambled);
            println!(
                "  -> supply predictions with `--stgabs-pred <name,value.csv>` to score an engine."
            );
        }
        Some(pp) => {
            let preds = load_named_preds(pp)?;
            println!("\n  scoring predictions: {pp}");
            let cases: [(&str, fn(&StgRow) -> bool); 4] = [
                ("all", |_r| true),
                ("D-A only", |r| r.bucket == Bucket::DonorAcceptor),
                ("MR-TADF only", |r| r.bucket == Bucket::MultiResonance),
                ("excl. TDDFT outliers", |r| !r.outlier),
            ];
            for (label, filter) in cases {
                let mut refs = Vec::new();
                let mut prd = Vec::new();
                let mut missing = 0;
                for r in rows.iter().filter(|r| filter(r)) {
                    // predictions may be keyed by folder id (`key`) or canonical name
                    match preds.get(&r.key).or_else(|| preds.get(&r.name)) {
                        Some(v) => {
                            refs.push(r.exp_rank);
                            prd.push(*v);
                        }
                        None => missing += 1,
                    }
                }
                if refs.len() >= 2 {
                    corr_line(label, &refs, &prd);
                    if missing > 0 {
                        println!("        ({missing} molecules had no prediction)");
                    }
                } else {
                    println!("    {label:<28} insufficient matched predictions");
                }
            }
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mut data_path = "data/reference_emitters.csv".to_string();
    let mut stgabs_path = "data/stgabs27.csv".to_string();
    let mut pred_path: Option<String> = None;
    let mut stgabs_pred: Option<String> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--pred" => {
                i += 1;
                match args.get(i) {
                    Some(p) => pred_path = Some(p.clone()),
                    None => {
                        eprintln!("error: --pred needs a file path");
                        return ExitCode::from(2);
                    }
                }
            }
            "--stgabs" => {
                i += 1;
                match args.get(i) {
                    Some(p) => stgabs_path = p.clone(),
                    None => {
                        eprintln!("error: --stgabs needs a file path");
                        return ExitCode::from(2);
                    }
                }
            }
            "--stgabs-pred" => {
                i += 1;
                match args.get(i) {
                    Some(p) => stgabs_pred = Some(p.clone()),
                    None => {
                        eprintln!("error: --stgabs-pred needs a file path");
                        return ExitCode::from(2);
                    }
                }
            }
            other => data_path = other.to_string(),
        }
        i += 1;
    }

    let emitters = match load_emitters(&data_path) {
        Ok(e) => e,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::from(1);
        }
    };

    println!("Loaded {} emitters from {}", emitters.len(), data_path);
    println!(
        "Emission tie-floor for λ ranking: ±{:.2} eV (experimental reference uncertainty).",
        EMISSION_TIE_EV
    );
    let _ = HC_EV_NM; // reserved for λ(nm)->eV ranking once enough λ data lands

    for bucket in [Bucket::DonorAcceptor, Bucket::MultiResonance] {
        print_bucket_report(bucket, &emitters);
        self_test(bucket, &emitters);
    }

    // STGABS27 ranking benchmark (the real D-A/MR ordering reference).
    if std::path::Path::new(&stgabs_path).exists() {
        if let Err(e) = handle_stgabs27(&stgabs_path, stgabs_pred.as_deref()) {
            eprintln!("error: {e}");
            return ExitCode::from(1);
        }
    }

    if let Some(p) = pred_path {
        if let Err(e) = score_predictions(&emitters, &p) {
            eprintln!("error: {e}");
            return ExitCode::from(1);
        }
    } else if stgabs_pred.is_none() {
        println!(
            "\nNo --pred file given. Metric machinery verified via self-test above.\n\
             Once the QC engine emits ΔEST per id, run:\n\
             \tcargo run -p rank-eval -- --pred <engine_output.csv>"
        );
    }

    ExitCode::SUCCESS
}
