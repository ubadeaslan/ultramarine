//! Owned quantum-chemistry core (Seçenek B).
//!
//! The SCF driver (`scf`) is integral-source-agnostic: it consumes precomputed
//! S, Hcore and ERI matrices. Integrals come either from the pure-Rust Stage-0
//! `integrals` module (Windows-native) or from libcint (crates/cint-sys, WSL).

pub mod aoeval;
pub mod basis;
pub mod grid;
pub mod integrals;
pub mod linalg;
pub mod scf;
pub mod tda;
