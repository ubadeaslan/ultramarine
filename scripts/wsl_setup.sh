#!/bin/bash
# One-time WSL build toolchain for the owned QC engine (Seçenek B).
# Run WITH sudo in your Ubuntu terminal (asks your password once):
#   sudo bash /mnt/c/Users/ubade/oled-mavi/scripts/wsl_setup.sh
set -e
echo ">> apt update"
apt-get update
echo ">> installing build toolchain + BLAS/LAPACK + python3-pip"
apt-get install -y \
  build-essential \
  cmake \
  gfortran \
  pkg-config \
  libopenblas-dev \
  liblapack-dev \
  python3-pip \
  python3-venv \
  git curl
echo
echo "DONE. Toolchain installed:"
for t in gcc g++ cmake make gfortran; do printf '  %-8s %s\n' "$t" "$($t --version 2>/dev/null | head -1)"; done
echo
echo "Now tell the assistant: 'toolchain kuruldu'."
