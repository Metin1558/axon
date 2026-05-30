# Axon

**MEA tabanlı organoid elektrofizyolojisi analiz pipeline'ı.**

NWB formatındaki beyin organoidi MEA kayıtlarını DANDI Archive'dan indirip uçtan uca analiz eder. Per-unit ateşleme hızı, CV, STTC senkronizasyon matrisleri, adaptif ağ burst tespiti, graph-teorik topoloji ölçütleri ve SpikeInterface ile spike sorting — tek CLI komutuyla.

DANDI:001603 (van der Molen et al., Nature Neuroscience 2026) ve DANDI:001132 (Andrews et al.) üzerinde doğrulandı.

---

## Dosya Yapısı

```
axon/
├── KURULUM.py                  ← İlk kurulum (bir kez çalıştır)
├── organoid_cli.py             ← Ana CLI: indir + analiz
├── dandi_ara.py                ← DANDI dataset arama
├── graph_analiz.py             ← Graph theory analizi (otomatik subject tespiti)
├── gruplu_analiz.py            ← Grup karşılaştırması (otomatik gruplama)
├── yas_metadata_cek.py         ← NWB metadata okuyucu
├── si_sorting.py               ← SpikeInterface spike sorting modülü
│
└── v6_3/                       ← Analiz modülleri
    ├── organoid_units_analiz.py ← Per-unit analiz + STTC + burst + plot
    ├── organoid_io.py           ← NWB okuma, ham sinyal, RAM koruması
    ├── organoid_signal.py       ← Bandpass + MAD spike detection
    ├── organoid_metrics.py      ← ISI, CV, STTC hesapları
    ├── organoid_qc.py           ← Refractory violation, Wilson CI, SNR
    ├── organoid_burst_ref.py    ← Bakkum 2013 referans burst karşılaştırması
    ├── organoid_replikasyon.py  ← Multi-session stabilite analizi
    ├── organoid_compare.py      ← İki kayıt karşılaştırması
    ├── organoid_plot.py         ← Görselleştirme modülü
    └── organoid_types.py        ← Veri tipleri
```

---

## Kurulum

```bash
# 1. Repoyu klonla
git clone https://github.com/metin/axon
cd axon

# 2. Kurulum (Python bağımlılıkları + SpikeInterface + algoritmalar)
python KURULUM.py
```

**Gereksinimler:** Python 3.9+, internet bağlantısı

**Otomatik kurulan paketler:**
`dandi`, `pynwb`, `h5py`, `numpy`, `scipy`, `matplotlib`, `networkx`, `remfile`, `pandas`, `spikeinterface[full]`, `mountainsort5`, `tridesclous2`

---

## Hızlı Başlangıç

```bash
# DANDI dataset ara
python dandi_ara.py organoid

# Dataset içindeki dosyaları listele
python dandi_ara.py --dataset 001603

# Subject indir ve analiz et (sorted units)
python organoid_cli.py 001603 sub-HO2

# Ham sinyal + SpikeInterface sorting
python organoid_cli.py 001132 sub-5C-1 --max-mb 750
# → Seçenek 2 → Algoritma seç (1: mountainsort5 önerilir)

# Graph theory analizi (sonuclar/ içindekilerden otomatik)
python graph_analiz.py

# Grup karşılaştırması
python gruplu_analiz.py

# Metadata özeti
python yas_metadata_cek.py
```

---

## Çıktılar

Her analiz sonrası `sonuclar/<subject>/` içine:

| Dosya | İçerik |
|-------|--------|
| `*_per_unit.csv` | Per-nöron metrikler: Hz, CV, STTC, refractory violation |
| `*_ozet.csv` | Organoid özeti: n_unit, burst rate, medyan STTC |
| `*_quicklook.png` | 4 panel: raster, population rate, Hz histogramı, STTC matrisi |

---

## Desteklenen Algoritmalar (SpikeInterface)

| # | Algoritma | Notlar |
|---|-----------|--------|
| 1 | mountainsort5 | Önerilen — hızlı, stabil |
| 2 | spykingcircus2 | Büyük kanallar için |
| 3 | tridesclous2 | Alternatif |
| 4 | simple | En hafif, düşük RAM |

---

## Doğrulama

**Sorted units pipeline:** 15/15 sentetik test geçti

**SpikeInterface entegrasyonu:** DANDI:001132 sub-5C-1 üzerinde doğrulandı
- 942 kanal HD-MEA, 50.1 sn kayıt
- mountainsort5 ile 9 unit sorted
- Medyan Hz: 1.15, STTC: -0.013

**DANDI:001603 tam dataset analizi:** 18 subject, ~1,610 nöron

---

## DANDI:001603 Özet Bulgular

| Kategori | n | Medyan Hz | Medyan CV | STTC | Burst/dk |
|----------|---|-----------|-----------|------|----------|
| İnsan organoidi (HO1–HO8) | 8 | 0.42 | 1.92 | 0.159 | 9.84 |
| Fare organoidi (MO10–MO14) | 5 | 4.42 | 1.15 | 0.014 | 5.40 |
| Ex vivo fare korteksi (M1S1–M3S2) | 5 | 0.39 | 2.83 | 0.144 | 10.10 |

**Yaş grubu farkı (insan organoidi):** Mann-Whitney p=0.029 — STTC ve CV

**Graph theory:** HO6 ve HO7 (her ikisi P100D) güçlü small-world topoloji (σ ≈ 4.3)

---

## Companion

**organoid-oi v3** — Kapalı döngü organoid zekası simülasyon sistemi  
→ [github.com/metin/organoid-oi](https://github.com/metin/organoid-oi)

---

*Axon v6.3 · MEA tabanlı organoid elektrofizyolojisi · Not peer reviewed · May 2026*
