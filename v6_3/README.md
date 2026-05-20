# organoid v6.3

Organoid elektrofizyoloji analiz aracı.

## v6.3 Değişiklikleri (v6.2 üzerinden)

1. **Byte tabanlı RAM güvenliği (organoid_io.py)** — Eski `MAX_CHUNK_SN = 300`
   tek kanal varsayımıyla konulmuştu. 256 kanallı MEA ile 300 sn = 9 GB RAM.
   Yeni: `MAX_CHUNK_BYTES = 500 MB`. `chunk_sn × n_kanal × sr × 4` formülüyle
   dinamik kontrol.

2. **STTC sirali=True parametre (organoid_metrics.py)** — Spike trenler
   zaten sıralı geliyor, gereksiz `np.sort` çağrıları kaldırıldı.

3. **Welch's t-test default (organoid_metrics.py)** — `equal_var=False`
   varsayılan. Eşit olmayan varyans varsayar — daha sağlam.

4. **Spike default çift yönlü (organoid_signal.py)** — `sadece_negatif=False`
   default. `--sadece-negatif` flag ile sadece negatif eşik.

5. **`--acf-esik` parametresi (organoid_compare.py)** — Hardcoded 0.3 yerine
   kullanıcı seçimi.

6. **SQLite BLOB depolama (organoid_db.py)** — `spike_arrays` tablosu,
   numpy serialize. 250 spike = 2000 byte tek satır.

7. **Dataclass veri taşıyıcıları (organoid_types.py)** — `KayitBilgisi`,
   `SpikeVerisi`, `SinyalKalitesi`, `AnalizParametre`. SADECE veri,
   davranış yok. Modülerlik korunur.

8. **`--plotly` opsiyonel grafik (organoid_plot.py)** — Default matplotlib
   PNG. `--plotly` ile interaktif HTML.

## Reddedilen Tavsiyeler

- **numpy.memmap** — NWB HDF5 formatında, h5py zaten chunk-based lazy
  loading yapıyor. Çakışır. `es.data[bas:bit, kanal]` zaten yeterli.

## Felsefe

1. **Modüller yorum yapmaz** — sadece sayısal ölçüm verir.
2. **Kontrol grubu hardcoded yok** — `organoid_compare.py` ile.
3. **SQLite ölçüm için, JSON metadata için** — karıştırma.
4. **Tüm istatistik organoid_metrics.py'de** — modüler bütünlük.
5. **Ham p değerleri** — manipülasyon yok, alt sınır yok.

## Modüller

```
organoid_io.py        — NWB lazy okuma, byte safety, multi-channel
organoid_signal.py    — band-pass, spike (default çift yönlü), MAD
organoid_qc.py        — kalite metrikleri + Wilson CI
organoid_metrics.py   — ISI/CV/FFT/Welch/permütasyon/Shapiro/STTC
organoid_sorting.py   — Tridesclous2 sarmalı (opsiyonel)
organoid_db.py        — SQLite WAL + BLOB depolama
organoid_plot.py      — quick-look (matplotlib, opsiyonel plotly)
organoid_output.py    — terminal + CSV
organoid_types.py     — dataclass veri taşıyıcılar
organoid.py           — ana komut
organoid_compare.py   — karşılaştırma + ACF + circular perm
```

## Kurulum

```
pip install pynwb scipy matplotlib numpy
pip install spikeinterface  # spike sorting için (opsiyonel)
pip install plotly          # interaktif grafik için (opsiyonel)
```

## Kullanım

```
python organoid.py <dosya.nwb> <organoid> <kayit> [seçenekler]
```

### Seçenekler

```
--kanal N         Okunacak kanal (varsayılan 0)
--snr SIGMA       Spike detection sigma (varsayılan 5)
--bandpass LOW HI Band-pass aralığı Hz (varsayılan 300 3000)
--chunk SN        Chunk boyutu sn (byte-safety var, 500 MB max)
--mad-esik        MAD/0.6745 ile gürültü std (Quiroga 2004)
--sadece-negatif  Sadece negatif spike eşiği (default: çift yönlü)
--sort            Tridesclous2 sorting (başarısızsa exit 3)
--sort-or-mua     Sort başarısızsa MUA fallback (etiketli)
--csv DOSYA       CSV çıktısı
--grafik DIR      Grafik klasörü
--no-grafik       Grafik oluşturma
--plotly          Plotly HTML çıktı (default: matplotlib PNG)
--metadata J      Metadata JSON
--refractory MS   Refractory eşiği ms (varsayılan 1.0)
--pencere SN      CV/aktiflik penceresi (varsayılan 30)
```

### Karşılaştırma

```
python organoid_compare.py <organoid1> <kayit1> <organoid2> <kayit2> [...]
                            [--acf-esik 0.3]
```

## Çıkış Kodları

- 0 : Başarılı
- 1 : Genel hata
- 2 : SpikeInterface yok (--sort)
- 3 : Spike sorting başarısız (--sort, fallback yok)

## İstatistik Yorumu

| Test | Ne zaman? | Yorum |
|------|-----------|-------|
| Welch t-test | Veri normal dağılıyorsa | Shapiro normal=True ise |
| Mann-Whitney | Normal dağılım yoksa | Genelde tercih edilir |
| Circular perm | Pencere bazlı autocorrelated veride | t-test'ten sağlam |
| ISI shuffle perm | İçsel ritim testi (tek kayıt) | Otokorelasyonu korur |

**Genel kural**: Birden fazla test "anlamlı" diyorsa güvenli.
Sadece biri diyorsa kuşkulan.

## Mimari Tarih

- v5: Tek dosya, kalite skoru, hardcoded kontrol, PDF
- v6: 10 modül, yorum yok, ham p, modüler ayrım
- v6.1: Negatif eşik, Wilson CI, Shapiro, sort fallback yok
- v6.2: NWB tek açılış, STTC, MAD, sort modları, circular perm
- v6.3: Byte RAM, BLOB, dataclass, Welch default, ACF param, plotly
