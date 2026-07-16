#!/bin/bash
# Build libcint (Gaussian integrals) from source into ~/.local (no sudo).
# Built in WSL-native home (fast); /mnt/c is slow under WSL.
set -e
PREFIX="$HOME/.local"
SRC="$HOME/qc-deps"
mkdir -p "$SRC" "$PREFIX"

if [ ! -e "$PREFIX/lib/libcint.so" ] && [ ! -e "$PREFIX/lib64/libcint.so" ]; then
  cd "$SRC"
  if [ ! -d libcint ]; then
    echo ">> cloning libcint"
    git clone --depth 1 https://github.com/sunqm/libcint.git
  fi
  cd libcint
  rm -rf build && mkdir build && cd build
  echo ">> cmake"
  cmake -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DCMAKE_BUILD_TYPE=Release \
        -DWITH_FORTRAN=0 \
        -DWITH_CINT2_INTERFACE=1 \
        -DBUILD_SHARED_LIBS=1 ..
  echo ">> make (-j$(nproc))"
  make -j"$(nproc)"
  make install
fi

echo "=== libcint artifacts ==="
ls -la "$PREFIX"/lib*/libcint* 2>/dev/null || echo "NONE"
echo "=== headers ==="
ls "$PREFIX"/include/cint*.h 2>/dev/null || echo "NONE"
