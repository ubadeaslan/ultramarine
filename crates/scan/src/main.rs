//! Rust orchestration layer (Seçenek A).
//!
//! Discovers `.xyz` geometries, runs the PySCF physics kernel
//! (`engine/pyscf_kernel.py`) on each in parallel, and writes:
//!   * `<out>`        — predictions CSV `name,deltaEST_eV` for `rank-eval`.
//!   * `<out>.jsonl`  — full per-molecule JSON records (provenance kept).
//!
//! The kernel command is configurable so the same binary works whether PySCF
//! lives in native Python or in WSL:
//!   --python "python"        (native; fails on Windows — PySCF needs Linux)
//!   --python "wsl python3"   (WSL; the working path on this machine)
//!
//! Usage:
//!   cargo run -p scan -- --geom-dir data/raw/geom --out data/stg_pred.csv \
//!       --python "wsl python3" --functional pbe0 --basis def2-svp [--jobs 4]
//!   cargo run -p rank-eval -- --stgabs-pred data/stg_pred.csv

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::sync::{Arc, Mutex};
use std::thread;

mod json {
    /// Extract a top-level string value for `key` from a flat JSON line.
    pub fn str_field(s: &str, key: &str) -> Option<String> {
        let pat = format!("\"{key}\"");
        let i = s.find(&pat)? + pat.len();
        let rest = &s[i..];
        let c = rest.find(':')? + 1;
        let rest = rest[c..].trim_start();
        let rest = rest.strip_prefix('"')?;
        let end = rest.find('"')?;
        Some(rest[..end].to_string())
    }

    /// Extract a top-level numeric value for `key` from a flat JSON line.
    pub fn num_field(s: &str, key: &str) -> Option<f64> {
        let pat = format!("\"{key}\"");
        let i = s.find(&pat)? + pat.len();
        let rest = &s[i..];
        let c = rest.find(':')? + 1;
        let tok: String = rest[c..]
            .trim_start()
            .chars()
            .take_while(|ch| ch.is_ascii_digit() || matches!(ch, '.' | '-' | '+' | 'e' | 'E'))
            .collect();
        tok.parse::<f64>().ok()
    }
}

struct Args {
    geom_dir: String,
    out: String,
    python: Vec<String>,
    functional: String,
    basis: String,
    solvent_eps: Option<String>,
    jobs: usize,
    kernel: String,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        geom_dir: "data/raw/geom".into(),
        out: "data/stg_pred.csv".into(),
        python: vec!["python".into()],
        functional: "pbe0".into(),
        basis: "def2-svp".into(),
        solvent_eps: None,
        jobs: 2,
        kernel: "engine/pyscf_kernel.py".into(),
    };
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let need = |i: &mut usize| -> Result<String, String> {
            *i += 1;
            argv.get(*i).cloned().ok_or_else(|| format!("{} needs a value", argv[*i - 1]))
        };
        match argv[i].as_str() {
            "--geom-dir" => a.geom_dir = need(&mut i)?,
            "--out" => a.out = need(&mut i)?,
            "--python" => a.python = need(&mut i)?.split_whitespace().map(String::from).collect(),
            "--functional" => a.functional = need(&mut i)?,
            "--basis" => a.basis = need(&mut i)?,
            "--solvent-eps" => a.solvent_eps = Some(need(&mut i)?),
            "--jobs" => a.jobs = need(&mut i)?.parse().map_err(|_| "bad --jobs")?,
            "--kernel" => a.kernel = need(&mut i)?,
            other => return Err(format!("unknown arg: {other}")),
        }
        i += 1;
    }
    if a.python.is_empty() {
        return Err("--python cannot be empty".into());
    }
    Ok(a)
}

struct Job {
    name: String,
    s1: PathBuf,
    t1: Option<PathBuf>, // present -> adiabatic ΔEST
}

/// Discover jobs. Two layouts (deterministic order):
///   * subdir per molecule containing `s1.xyz` (+ optional `t1.xyz`) -> adiabatic
///   * flat `*.xyz` files -> vertical (single geometry)
fn collect_jobs(dir: &str) -> Result<Vec<Job>, String> {
    let mut jobs = Vec::new();
    let mut entries: Vec<PathBuf> = fs::read_dir(dir)
        .map_err(|e| format!("read_dir {dir}: {e}"))?
        .map(|e| e.map(|e| e.path()).map_err(|e| e.to_string()))
        .collect::<Result<_, _>>()?;
    entries.sort();
    for p in entries {
        if p.is_dir() {
            let s1 = p.join("s1.xyz");
            if s1.exists() {
                let t1 = p.join("t1.xyz");
                jobs.push(Job {
                    name: p.file_name().and_then(|s| s.to_str()).unwrap_or("unnamed").to_string(),
                    s1,
                    t1: t1.exists().then_some(t1),
                });
            }
        } else if p.extension().and_then(|x| x.to_str()) == Some("xyz") {
            jobs.push(Job {
                name: p.file_stem().and_then(|s| s.to_str()).unwrap_or("unnamed").to_string(),
                s1: p,
                t1: None,
            });
        }
    }
    Ok(jobs)
}

#[derive(Debug)]
struct Outcome {
    name: String,
    json: String,
    delta_est: Option<f64>,
    error: Option<String>,
}

fn run_one(a: &Args, job: &Job) -> Outcome {
    let name = job.name.clone();
    let mut cmd = Command::new(&a.python[0]);
    cmd.args(&a.python[1..]);
    cmd.arg(&a.kernel)
        .arg("--xyz").arg(&job.s1)
        .arg("--name").arg(&name)
        .arg("--functional").arg(&a.functional)
        .arg("--basis").arg(&a.basis);
    if let Some(t1) = &job.t1 {
        cmd.arg("--xyz-t1").arg(t1);
    }
    if let Some(eps) = &a.solvent_eps {
        cmd.arg("--solvent-eps").arg(eps);
    }
    match cmd.output() {
        Ok(o) => {
            let line = String::from_utf8_lossy(&o.stdout).trim().to_string();
            let delta_est = json::num_field(&line, "deltaEST_eV");
            let error = json::str_field(&line, "error");
            Outcome { name, json: line, delta_est, error }
        }
        Err(e) => Outcome {
            name: name.clone(),
            json: format!("{{\"name\": \"{name}\", \"error\": \"spawn failed: {e}\"}}"),
            delta_est: None,
            error: Some(format!("spawn failed: {e}")),
        },
    }
}

fn main() -> std::process::ExitCode {
    let a = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("error: {e}");
            return std::process::ExitCode::from(2);
        }
    };
    let jobs = match collect_jobs(&a.geom_dir) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("error: {e}");
            return std::process::ExitCode::from(1);
        }
    };
    if jobs.is_empty() {
        eprintln!(
            "no geometries in {} — expected subdirs with s1.xyz or flat *.xyz files.",
            a.geom_dir
        );
        return std::process::ExitCode::from(1);
    }
    let n_adia = jobs.iter().filter(|j| j.t1.is_some()).count();
    println!(
        "scan: {} molecules ({} adiabatic), kernel='{} {}', {}/{} basis, {} jobs",
        jobs.len(),
        n_adia,
        a.python.join(" "),
        a.kernel,
        a.functional,
        a.basis,
        a.jobs
    );

    let a = Arc::new(a);
    let queue = Arc::new(Mutex::new(jobs.into_iter().enumerate().collect::<Vec<_>>()));
    let results: Arc<Mutex<Vec<(usize, Outcome)>>> = Arc::new(Mutex::new(Vec::new()));

    let mut handles = Vec::new();
    for _ in 0..a.jobs.max(1) {
        let a = Arc::clone(&a);
        let queue = Arc::clone(&queue);
        let results = Arc::clone(&results);
        handles.push(thread::spawn(move || loop {
            let job = queue.lock().unwrap().pop();
            let Some((idx, j)) = job else { break };
            let out = run_one(&a, &j);
            let tag = match (&out.error, out.delta_est) {
                (Some(e), _) => format!("ERROR: {e}"),
                (None, Some(d)) => format!("ΔEST={d:.3} eV"),
                (None, None) => "no ΔEST parsed".to_string(),
            };
            println!("  [{}] {}", out.name, tag);
            results.lock().unwrap().push((idx, out));
        }));
    }
    for h in handles {
        let _ = h.join();
    }

    let mut res = Arc::try_unwrap(results).unwrap().into_inner().unwrap();
    res.sort_by_key(|(i, _)| *i); // restore deterministic input order

    let mut pred = String::from("name,deltaEST_eV\n");
    let mut jsonl = String::new();
    let mut ok = 0;
    for (_, o) in &res {
        jsonl.push_str(&o.json);
        jsonl.push('\n');
        if let Some(d) = o.delta_est {
            pred.push_str(&format!("{},{:.4}\n", o.name, d));
            ok += 1;
        }
    }
    if let Err(e) = fs::write(&a.out, pred) {
        eprintln!("error writing {}: {e}", a.out);
        return std::process::ExitCode::from(1);
    }
    let jsonl_path = format!("{}.jsonl", a.out);
    if let Err(e) = fs::File::create(&jsonl_path).and_then(|mut f| f.write_all(jsonl.as_bytes())) {
        eprintln!("error writing {jsonl_path}: {e}");
        return std::process::ExitCode::from(1);
    }
    println!("\nwrote {} predictions -> {}  (full records -> {})", ok, a.out, jsonl_path);
    println!("next: cargo run -p rank-eval -- --stgabs-pred {}", a.out);
    std::process::ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::json;

    #[test]
    fn parse_success_line() {
        let line = r#"{"name": "MCz-XT", "natom": 42, "deltaEST_eV": 0.083, "E_S1_eV": 2.71, "f_osc_S1": 0.012, "provenance": {"functional": "pbe0"}}"#;
        assert_eq!(json::str_field(line, "name").as_deref(), Some("MCz-XT"));
        assert!((json::num_field(line, "deltaEST_eV").unwrap() - 0.083).abs() < 1e-9);
        assert!(json::str_field(line, "error").is_none());
    }

    #[test]
    fn parse_error_line() {
        let line = r#"{"name": "H2TEST", "error": "pyscf not available: No module named 'pyscf'"}"#;
        assert_eq!(json::str_field(line, "error").as_deref(), Some("pyscf not available: No module named 'pyscf'"));
        assert!(json::num_field(line, "deltaEST_eV").is_none());
    }

    #[test]
    fn parse_negative_number() {
        let line = r#"{"deltaEST_eV": -0.0123}"#;
        assert!((json::num_field(line, "deltaEST_eV").unwrap() + 0.0123).abs() < 1e-9);
    }
}
