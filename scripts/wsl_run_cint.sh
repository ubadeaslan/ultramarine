#!/bin/bash
# Build+run the cint-sys FFI test in WSL. Uses a WSL-native target dir
# (building under /mnt/c is slow and flaky for cargo).
set -e
. "$HOME/.cargo/env"
export CARGO_TARGET_DIR="$HOME/cargo-target/cint-sys"
cd /mnt/c/Users/ubade/oled-mavi/crates/cint-sys
cargo run --release 2>&1 | tail -30
