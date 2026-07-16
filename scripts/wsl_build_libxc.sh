#!/bin/bash
# Build libxc (exchange-correlation functionals) into ~/.local (no sudo).
# Needed for Stage 2 (KS-DFT). Built in WSL-native home.
set -e
PREFIX="$HOME/.local"
SRC="$HOME/qc-deps"
mkdir -p "$SRC" "$PREFIX"

if [ ! -e "$PREFIX/lib/libxc.so" ] && [ ! -e "$PREFIX/lib64/libxc.so" ] && \
   [ ! -e "$PREFIX/lib/libxc.a" ]; then
  cd "$SRC"
  if [ ! -d libxc ]; then
    echo ">> cloning libxc"
    git clone --depth 1 https://gitlab.com/libxc/libxc.git
  fi
  cd libxc
  rm -rf objdir
  echo ">> cmake configure"
  cmake -H. -Bobjdir \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_FORTRAN=OFF \
        -DENABLE_PYTHON=OFF \
        -DBUILD_TESTING=OFF \
        -DBUILD_SHARED_LIBS=ON
  echo ">> build (-j$(nproc))"
  cmake --build objdir -j"$(nproc)"
  cmake --install objdir
fi

echo "=== libxc artifacts ==="
ls -la "$PREFIX"/lib*/libxc* 2>/dev/null || echo NONE
ls "$PREFIX"/include/xc.h 2>/dev/null || echo "no xc.h"
