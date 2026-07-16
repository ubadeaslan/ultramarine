# Referans Set — Metodoloji ve Yol Haritası

Bu doküman, derin-mavi TADF tarama motorunun **doğrulama harness'ının** temelini
oluşturan deneysel referans setini tanımlar. Motorun değeri buna bağlı: motor
"bilinen emitörleri doğru sıralıyor mu?" sorusu bu setle yanıtlanır.

## Felsefe

- **Sıralama doğruluğu > mutlak doğruluk.** Birincil metrik: motorun ΔEST
  sıralaması ile deneysel ΔEST sıralaması arasındaki Spearman/Kendall korelasyonu.
- **Her satır izlenebilir.** Değer + ölçüm ortamı (çözücü/host) + DOI + doğrulama
  durumu. Ortam kritik: literatür ΔEST'i 77 K başlangıç vs. kinetik fit vs. host'a
  göre 50+ meV oynayabilir. Provenance olmadan "doğru sıralama" testi kendi
  gürültüsünde boğulur.
- **Açık erişim öncelikli.** Tüm değerler açık kaynaklardan doğrulanır; paywall
  SI'lar yalnızca erişim varsa kalite artırmak için kullanılır.

## İki kova (KARIŞTIRMA)

D-A ve MR-TADF **farklı fiziktir** ve tek sıralama metriğinde karıştırılırsa
sonuç anlamsızlaşır:

1. **D-A tipi** — S1 yük-transfer (CT) karakterli. ΔEST, uzaysal olarak ayrılmış
   HOMO/LUMO arasındaki değişim integraliyle yönetilir. Yöntem: menzil-ayrılmış
   hibrit (CAM-B3LYP / wB97X-D, ideali IP-tuned RSH) + TDA.
2. **MR-TADF tipi** (DABNA ailesi) — kısa-menzilli korelasyon hakim, CT değil.
   **Standart TDDFT (RSH dahil) burada başarısız.** Yöntem: **ΔDFT (ΔSCF tipi)** —
   dalga-fonksiyonu doğruluğu, DFT maliyetinde.

## İki kritik bilimsel kısıt (mimariyi belirler)

### 1. "Ya doğru ΔEST, ya doğru geometri — ikisi birden değil"
*JCTC* 2022, 10.1021/acs.jctc.2c00905 (STGABS27 üzerinde gösterildi): TD-DFT bir
fonksiyonelle ΔEST'i tutturursa relaxed-S1 geometrisini (→ emisyon λ) kaçırır,
tersi de geçerli.

**Sonuç:** Pipeline ΔEST ve λ için **farklı yöntem/fonksiyonel** kullanmalı.
Provenance katmanı property-başına method-tag taşımalı (tek "doğru fonksiyonel"
yok). Emisyon λ relaxed-S1 geometrisinden gelir — S0 geometrisi absorpsiyon verir,
emisyon değil.

### 2. Kararlılık halkası dürüstçe ayrı
Operasyonel ömür (LT50/LT95) cihaz-bağımlı bir bozunma süreci ve tek-molekül
temel-hal DFT'sinin erişemediği yerde. Fizikle yapılabilecek en fazlası vekil
göstergeler (T1/katyon halinde bağ ayrışma enerjileri). Ömre bağlama = dark
data'ya fit. Bu, motorun ana çıktısı değil, dürüstçe işaretlenmiş yan filtredir.
**"ML yok" sınırı:** bu fit'in az-parametreli/fiziksel mi yoksa istatistiksel mi
sayılacağı netleştirilmeli (açık karar bekliyor).

## Mevcut durum (v0.1)

`data/reference_emitters.csv` — 10 tohum emitör:
- D-A: 4CzIPN (yeşil çapa), 2CzPN, 4CzPN, DMAC-DPS, DMAC-TRZ, PXZ-TRZ (yeşil
  renk-kontrolü), ACRSA
- MR-TADF: DABNA-1, DABNA-2, nu-DABNA

OPEN_VERIFIED satırlar test edilebilir; NEEDS_* satırlar değer/DOI bekliyor.

## STGABS27 — içe aktarıldı (sıralama referansı)

`data/stgabs27.csv` — 27 molekül, **deneysel ΔEST'e göre artan rank sırası** ile.
Kaynak: Froitzheim/Grimme/Mewes JCTC 2022 SI figürleri (chemrxiv-2022-rdtxc) +
özgün STGABS27 (Kunze 2021, 10.1021/acs.jpclett.1c02299).

- **Mutlak eV değerleri henüz yok:** SI'da yalnızca figür olarak var, metin tablosu
  değil → uydurmuyoruz. Ama **rank sırası** sıralama-metriği (Spearman/Kendall) için
  yeterli ve matematiksel olarak tam doğru referans.
- Kovalar: 25 D-A + 2 MR-TADF (DABNA-1/2). Birkaç boron-içeren molekül (TMCz-BO,
  p/m'-AC-DBNA) MR olabilir → `notes`'ta "verify" işaretli.
- **TDDFT aykırıları:** 8, 9, 21, 23, 26, 27 (DABNA'lar multireferans → standart
  TDDFT başarısız). Harness "excl. TDDFT outliers" alt-skorunu ayrı verir.
- Skorlama: `cargo run -p rank-eval -- --stgabs-pred <name,value.csv>`

Mutlak eV değerlerini tamamlamak için: özgün 2021 PCM-ROKS SI tablosu (kullanıcı
indirip `data/raw/`'a koyacak) veya figür vektör-noktalarından kalibrasyon.

## İçe-aktarma kuyruğu (kalan)

1. STGABS27 mutlak eV değerleri (rank yerine) — opsiyonel hassasiyet artışı.
2. **747-set mavi alt-örneği (+15-20)** — repo yayınlanınca λ < 480 nm filtrele.
3. **MR-TADF genişletme** — 747-set MR + DABNA türevleri (n=2 çok az).

## Doğrulanmış metodoloji (Bonn/Herbert SI'dan — açık erişim)

Froitzheim, Kunze, Grimme, Herbert, Mewes (*J. Phys. Chem. A* 2024, 128, 6324;
SI açık) ve önceki STGABS27 çalışmalarından çıkarılan, motor reçetesini doğrudan
belirleyen pratik sonuçlar:

- **STGABS27'nin sahibi:** Bonn (Mulliken Center) + Herbert grubu — Kunze, Hansen,
  Grimme, Mewes (JPCL 2021, 10.1021/acs.jpclett.1c02299). (Önceki "Marian" notu
  düzeltildi.)
- **En iyi fonksiyonel (D-A, ΔEST):** optimally-tuned **OT-ωB97M-V**. ω parametresi
  sistem-başına ayarlanır (sabit RSH değil). Bu, "IP/OT-tuned RSH şart" tezini
  doğruluyor.
- **Emisyon için ΔDFT > TDDFT:** ΔUKS / ΔROKS, TDA-DFT'yi emisyon enerjilerinde
  yener. MR-TADF *ve* D-A için ΔDFT dalı ciddi avantaj.
- **Baz set:** **def2-SVP yeterli** — bu CT halleri için triple-ζ'ya (def2-TZVPP)
  göre fark < 0.1 eV. Büyük maliyet tasarrufu; tarama motoru için kritik.
- **Çözücü:** IEF-PCM, **ölçüme-özgü ε ve n** kullan (CBP host için ε ≈ 3.5).
  CT ve LE halleri dielektrik ortama farklı tepki verir; state-specific (ptSS-PCM)
  düzeltme uygulanmalı.
- **Gürültü tabanı:** deneysel **emisyon** referansının belirsizliği **±0.2 eV**.
  Yani λ/renk sıralamasında ±0.2 eV içindeki farklar fiziksel olarak ayrışmaz —
  bu motorun kusuru değil, doğanın sınırı. (ΔEST referansları STGABS27'de daha sıkı,
  "kimyasal doğruluk" iddiasıyla.) M0 metriği bu tabanı hesaba katmalı: ±eşik
  içindeki çiftleri "berabere" say.

## İnsan-kapılı tek blokaj (kullanıcı aksiyonu)

STGABS27'nin sayısal tablosu + optimize geometrileri yayıncı SI .zip arşivinde;
ChemRxiv/ACS Cloudflare+paywall arkasında, otomatik (bot) erişim engelleniyor.
Tarayıcıdan (EKUAL + Cloudflare'i sen geçersin) indirilecek 2 dosya işi bitirir:

1. **STGABS27 SI .zip** — JPCL 2021 (10.1021/acs.jpclett.1c02299) veya JCTC 2022
   (10.1021/acs.jctc.2c00905) "Supporting Information". İçinde: optimize XYZ
   geometriler (→ geometri-opt adımını atlatır) + deneysel ΔEST tablosu.
2. İndirip `data/raw/` içine koy; ben parse edip `reference_emitters.csv`'ye
   27 satırı işlerim.

## Sıradaki adım

- [ ] OPEN_VERIFIED satırlarla M0 sıralama testinin iskeletini kur (henüz hesap yok,
      sadece harness + metrik).
- [ ] STGABS27 üye listesini açık ChemRxiv SI'sından çek.
- [ ] Kullanıcı ACS/Wiley erişimini doğrularsa NEEDS_* satırları SI'dan tamamla.
