# Kaynaklar — Referans Emitör Seti (Açık Erişim)

Tüm kayıtlar açık kaynaklardan doğrulandı veya doğrulama bekliyor. `verification`
sütunu kuralları:

- `OPEN_VERIFIED` — değer bu oturumda açık bir kaynaktan (yayıncı HTML, arXiv,
  açık PMC, üretici veri sayfası + birincil makaleyle çapraz kontrol) okundu.
- `SECONDARY` — değer ikincil kaynaktan (derleme, üretici sayfası); birincil
  SI ile doğrulanması gerekiyor.
- `NEEDS_SI` — değer henüz çekilmedi; birincil makale SI'sından alınacak.
- `NEEDS_DOI` — değer var ama birincil DOI doğrulanmadı.

## Birincil makaleler (DOI)

| Emitör(ler) | Makale | DOI | Erişim |
|---|---|---|---|
| 4CzIPN, 2CzPN, 4CzPN | Uoyama et al., *Nature* 492, 234 (2012) | 10.1038/nature11687 | Paywall (Nature) — değerler açık ikincil kaynaklarla doğrulandı |
| DMAC-DPS | Zhang et al., *Nat. Photonics* (2014) | 10.1038/nphoton.2014.12 | DOI doğrulanacak |
| DMAC-TRZ | Tsai et al., *Chem. Commun.* 51, 13662 (2015) | 10.1039/c5cc05022g | RSC ✓ |
| PXZ-TRZ | Tanaka et al., *Chem. Commun.* 48, 11392 (2012) | 10.1039/c2cc36237f | RSC ✓ |
| ACRSA | Nasu et al., *Chem. Commun.* 49, 10385 (2013) | 10.1039/c3cc44179b | RSC ✓ |
| DABNA-1, DABNA-2 | Hatakeyama et al., *Adv. Mater.* 28, 2777 (2016) | 10.1002/adma.201505491 | Wiley — EKUAL üzerinden erişilebilir olabilir |
| nu-DABNA | Kondo et al., *Nat. Photonics* 13, 678 (2019) | 10.1038/s41566-019-0476-5 | Paywall |

## Benchmark setleri (içe-aktarma kuyruğu)

| Set | İçerik | Kaynak | Erişim |
|---|---|---|---|
| **STGABS27** | 27 TADF emitör, deneysel adiyabatik ΔEST | PCM-ROKS, *JPCL* 2021, 10.1021/acs.jpclett.1c02299; ayrıca JCTC 2022 10.1021/acs.jctc.2c00905 | **İçe aktarıldı** (rank sırası) → `data/stgabs27.csv`. Kaynak: chemrxiv-2022-rdtxc SI figürleri (açık). Mutlak eV pending. |
| **747-set** | 747 emitör, 312 deneysel ΔEST + 213 deneysel λ | arXiv 2511.00922 (2025) | Açık; veri reposu "yayında gelecek" |
| **xTB-8** | 8 temsili emitör (DMAC-TRZ, DMAC-DPS, PSPCz, 4CzIPN, Px2BP, CzS2, 2TCz-DPS, TDBA-DI) | arXiv 2502.20410 (2025) | Açık HTML |

## Geometriler (motor girdisi)

`data/geom/stgabs27/<NN_ad>/{s1.xyz,t1.xyz}` — 27 molekülün PBE0-D4 ile optimize
edilmiş **S1 (ROKS)** ve **T1 (UKS)** geometrileri. Kaynak: JPCA 2024 SI arşivi
(`jp4c03273_si_002.zip`, 10.1021/acs.jpca.4c03273). Motor adiyabatik ΔEST'i
E(S1@S1geom) − E(T1@T1geom) olarak bu iki geometriden hesaplar. (S0 geometrisi
arşivde yok; adiyabatik yaklaşım zaten STGABS27'nin doğru tanımı.) Zip ayrıca her
fonksiyonel için Q-Chem job.out'larını içerir → ileride bizim ΔEST'imizi onların
hesabıyla çapraz-doğrulamak için kullanılabilir.

## Kararlılık / ömür verisi (cihaz-bağımlı, ayrı kova)

Ömür (LT50/LT95) tek-molekül DFT'nin dışında — dürüstçe ayrı tutuluyor
(bkz. REFERENCE_SET.md). Açık kaynaklardan toplanan örnek değerler:

- 4CzBN: LT50 = 62 saat, 4TCzBN: LT50 = 167 saat @ L0 = 500 cd/m^2
  (Nat. Commun. 2023, 10.1038/s41467-023-39697-7)
- t-DABNA (sadece TE): LT95 = 19.8 saat, EQE1000 = %23.7
- Hedef referans: mavi TADF için T95 ~ 500 saat @ 1000 cd/m^2 (ticari eşik)
