# Ultramarine

**Saf-fizik (ML YOK), deterministik bir derin-mavi TADF emitör tarama motoru — kendi
kuantum kimyası çekirdeğiyle.**

> *Ultramarine, tarihin en pahalı pigmentiydi — lapis lazuli'den öğütülür, bir dönem
> altından değerli sayılırdı. Kararlı derin-mavi emitör de OLED'in en değerli,
> 20 yıldır çözülememiş problemi.*

---

## Ne yapar?

İridyumsuz (metalsiz), derin-mavi (~450–470 nm) **TADF** OLED emitör adaylarını,
sentezden önce, **fizikle hesaplayıp** eler. Her aday için:

| Büyüklük | Ne demek | Hedef |
|---|---|---|
| **ΔEST** | S₁–T₁ enerji aralığı | küçük (TADF'in çalışma şartı) |
| **f** | osçilatör gücü (ışıma) | yüksek |
| **λ** | emisyon dalga boyu (renk) | 450–470 nm |
| *(planlanan)* kararlılık | bozunma vekil-göstergeleri | yüksek |

Çıktı: küçük-ΔEST ↔ yüksek-ışıma **Pareto cephesi** + renk penceresi + (planlanan)
kararlılık filtresi — hepsi **tam provenance** ile (hangi fonksiyonel, hangi baz set,
hangi yöntem).

**Felsefe: ML yok. Tahmin değil, türetme.** Her sayı bir denklemden gelir ve izlenebilir.

## Neden kendi motorumuz?

Hazır bir QC motorunu (PySCF vb.) çağırmak yerine çekirdeği kendimiz yazdık. Aynı
temel kayayı (libcint = integraller, libxc = fonksiyoneller — bunları herkes kullanır,
yeniden yazmak akılsızlık) kullanıp **SCF / KS-DFT / TDDFT-response mantığını** kendimiz
yazıyoruz. Sebep: başkasının motoruyla o alanda en iyi olunamaz; özelleşerek olunur.

## Durum: çekirdek kanıtlı doğru ✅

Her katman **PySCF ile birebir doğrulandı** (altın referans / oracle):

| Yöntem | Ultramarine | PySCF | Δ |
|---|---|---|---|
| H₂ PBE0 | -1.154311 | -1.154311 | **0** |
| H₂O PBE0 | -75.245849 | -75.245851 | 2e-6 |
| H₂O **CAM-B3LYP** (menzil-ayrılmış) | -75.279842 | -75.279857 | 1.5e-5 |
| H₂O PBE0/def2-SVP | -76.276179 | -76.276274 | 9.5e-5 |
| H₂O **CIS ΔEST / f / λ** | 2.0991 eV / 0.0035 / 94.2 nm | aynı | **0** |
| H₂O **TDA-DFT ΔEST** | 1.9077 eV | 1.9077 | **0** |

En büyük sapma 9.5e-5 Ha = **0.06 kcal/mol** (kimyasal doğruluğun ~16× altında).

**Yetenekler:** HF · LDA · PBE · PBE0 · CAM-B3LYP(RSH) · CIS · TDA-DFT → ΔEST + f + λ
**Altyapı:** libcint (+erf menzil-ayrımı) · libxc · Becke/Lebedev ızgara · direct-SCF ·
8-kat ERI simetrisi · Schwarz eleme · AO-grid eleme · DIIS · rayon (16 çekirdek)

## Yapı

```
crates/
  qc-core/     kuantum çekirdeği: basis · grid · aoeval · integrals · scf · tda · linalg
  cint-sys/    libcint + libxc FFI, sürücü (WSL/Linux)
  rank-eval/   doğrulama harness'ı: Spearman/Kendall, iki fizik kovası ayrı
  scan/        orkestrasyon: adaylar üzerinde paralel koşturma
engine/        pyscf_kernel.py — PySCF (artık doğrulama oracle'ı)
data/          referans set · STGABS27 sıralaması · 27 emitör geometrisi · bazlar
docs/          YOLCULUK.md (baştan sona anlatı) · ENGINE_BUILD.md · REFERENCE_SET.md
scripts/       WSL kurulum · libcint/libxc derleme · PySCF karşılaştırma · baz dökümü
```

## Kurulum & çalıştırma

libcint/libxc C kütüphaneleri gerektiği için motor **WSL/Linux**'ta derlenir
(saf-Rust katman Windows'ta da çalışır). Adımlar: [docs/ENGINE_BUILD.md](docs/ENGINE_BUILD.md)

```sh
cargo test                                   # saf-Rust katman testleri (13/13)
bash scripts/wsl_setup.sh                    # (sudo) derleme araçları
bash scripts/wsl_rust.sh                     # Rust
bash scripts/wsl_build_libcint.sh            # libcint
bash scripts/wsl_build_libxc.sh              # libxc
bash scripts/wsl_run_cint.sh                 # motor + doğrulama bataryası
python3 scripts/pyscf_ref.py                 # PySCF referansı (karşılaştırma)
```

## Yol haritası

- [ ] **SAD başlangıç tahmini** — büyük-molekül SCF yakınsaması (tek blokaj)
- [ ] **RSH-TDA + GGA TDDFT kernel** — CT-doğru ΔEST
- [ ] **27 STGABS27 emitörünü sırala** → deneysel sırayla karşılaştır (nihai test)
- [ ] **Kararlılık modülü** — BDE(T₁/katyon) vekil-göstergeleri + dark-data kalibrasyonu

## Veri & atıf

Referans set **STGABS27** (Kunze, Hansen, Grimme, Mewes) ve geometriler yayınlanmış
Supporting Information'lardan türetilmiştir; kaynaklar ve DOI'ler:
[data/SOURCES.md](data/SOURCES.md), [docs/REFERENCE_SET.md](docs/REFERENCE_SET.md).
Ham yayıncı PDF'leri/arşivleri repoya dahil değildir (telif + boyut).

## Hikâye

Sıfırdan bu motora nasıl gelindiği, hiçbir şey bilmeyen birine anlatır gibi:
**[docs/YOLCULUK.md](docs/YOLCULUK.md)**
