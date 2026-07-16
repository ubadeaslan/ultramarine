#!/bin/bash
# Install Rust toolchain in WSL (no sudo; installs to ~/.cargo).
set -e
if [ ! -x "$HOME/.cargo/bin/rustc" ] && ! command -v rustc >/dev/null 2>&1; then
  echo ">> installing rustup (minimal profile)"
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal
fi
. "$HOME/.cargo/env" 2>/dev/null || true
echo "rustc: $(rustc --version)"
echo "cargo: $(cargo --version)"
