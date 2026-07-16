#!/bin/bash
# WSL toolchain probe for the owned QC engine build.
for t in gcc g++ cmake make gfortran git python3 pip3 curl; do
  if command -v "$t" >/dev/null 2>&1; then
    printf '%-9s %s\n' "$t" "$("$t" --version 2>/dev/null | head -1)"
  else
    printf '%-9s MISSING\n' "$t"
  fi
done
echo "--- os ---"
. /etc/os-release && echo "$PRETTY_NAME"
echo "--- blas/lapack (apt) ---"
for lib in libopenblas-dev liblapack-dev; do
  dpkg -s "$lib" >/dev/null 2>&1 && echo "$lib installed" || echo "$lib MISSING"
done
