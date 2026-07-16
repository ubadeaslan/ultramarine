# Motor Derleme (Seçenek B, WSL)

Sahip olunan QC çekirdeği libcint (integraller) + libxc (fonksiyoneller) C
kütüphanelerine FFI ile bağlanır → **WSL/Linux'ta derlenir.** Saf-Rust Stage-0
spike (`qc-core`) ayrıca Windows'ta native çalışır.

## Ortam (bir kez)

```bash
# 1) derleme araçları (sudo, kullanıcı şifresi):
sudo bash /mnt/c/Users/ubade/oled-mavi/scripts/wsl_setup.sh
#    -> build-essential, cmake, gfortran, openblas, lapack, pip

# 2) Rust (sudo'suz):
bash /mnt/c/Users/ubade/oled-mavi/scripts/wsl_rust.sh        # rustup -> ~/.cargo

# 3) libcint (sudo'suz, ~/.local'a):
bash /mnt/c/Users/ubade/oled-mavi/scripts/wsl_build_libcint.sh
#    -> ~/.local/lib/libcint.so (+ headers).  libxc Stage 2'de eklenecek.
```

Doğrulanmış sürümler: Ubuntu 26.04, gcc 15.2, cmake 4.2, rustc 1.96, libcint 6.1.3.

## Derleme/çalıştırma (WSL)

`cint-sys` crate'i libcint'e linklenir; Windows workspace'inden `exclude` edilir.
WSL-native target dizini kullan (/mnt/c yavaş):

```bash
bash /mnt/c/Users/ubade/oled-mavi/scripts/wsl_run_cint.sh
# = cargo run --release  (CARGO_TARGET_DIR=~/cargo-target)
```

## Windows'tan WSL çağırırken

Git Bash `/mnt/...` yollarını bozar → komutun başına `MSYS_NO_PATHCONV=1` koy:

```sh
MSYS_NO_PATHCONV=1 wsl -- bash /mnt/c/Users/ubade/oled-mavi/scripts/<script>.sh
```

## Durum

- ✅ libcint FFI: H2/STO-3G örtüşme matrisi = ders kitabı (S12=0.6593).
- Sıradaki: T, V, ERI (cint2e_sph) → qc-core'un el-yazımı integrallerini libcint ile
  değiştir → H2 RHF'i libcint'le -1.1167 Ha doğrula → gerçek moleküller (def2-SVP).
