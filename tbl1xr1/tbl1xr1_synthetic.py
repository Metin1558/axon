"""
tbl1xr1_synthetic.py — Mastrototaro 2021'in GERÇEK Sayılarına Kalibre Sentetik Veri
=======================================================================================

axon-ndd'deki sentetik_organoid.py'den farkı: orada genel/keyfi
parametrelerle "bir yön üretmeye çalışıyorduk". Burada, doğrudan
Mastrototaro'nun makalesinde YAYINLANMIŞ sayıları (WT 4.0±0.3 Hz,
KO 6.3±0.4 Hz, n=19 hücre) kullanıyoruz — sentetik veri, gerçek bir
deneyin gerçek sonucunun bir simülasyonu, keyfi bir hedef değil.

Bu, "bu pipeline gerçek yayınlanmış bir bulguyu doğru tespit
edebiliyor mu" sorusuna cevap veren bir SANITY CHECK'tir — organoid
verisini TAHMİN ETMEZ, sadece pipeline'ın matematiğini doğrular.
"""

import numpy as np


def poisson_spike_treni(rate_hz, sure_sn, rng):
    n_beklenen = int(rate_hz * sure_sn * 1.5) + 10
    isi = rng.exponential(1.0 / rate_hz, size=n_beklenen)
    zaman = np.cumsum(isi)
    return zaman[zaman < sure_sn]


# Mastrototaro 2021, Figure 2C — GERÇEK yayınlanmış değerler
MASTROTOTARO_WT_ORTALAMA_HZ = 4.0
MASTROTOTARO_WT_SEM = 0.3
MASTROTOTARO_KO_ORTALAMA_HZ = 6.3
MASTROTOTARO_KO_SEM = 0.4
MASTROTOTARO_N_HUCRE = 19  # her grupta


def tbl1xr1_sentetik_organoid_uret(n_unit, sure_sn, genotip_artis_mi,
                                    hucre_arasi_varyasyon_sd=None, seed=42):
    """
    genotip_artis_mi:
        True  -> KO/mutant kalibrasyonu (ortalama 6.3 Hz, Mastrototaro SEM'i
                 hücreler-arası SD'ye kabaca çevrilmiş: SEM*sqrt(n) ≈ SD)
        False -> WT/kontrol kalibrasyonu (ortalama 4.0 Hz)
        None  -> BELİRSİZ genotipler (G70N/L282P/Y446C) için — literatürde
                 sayı YOK, bu fonksiyon bilerek NotImplementedError fırlatır,
                 sahte bir sayı uydurmuyoruz.
    """
    if genotip_artis_mi is None:
        raise NotImplementedError(
            'Bu genotip için Mastrototaro\'da (ya da başka bir kaynakta) '
            'yayınlanmış bir sayı YOK. Sentetik kalibrasyon uydurmak yerine '
            'bilerek burada duruyoruz — gerçek veri gelmeden bu genotip için '
            'pipeline\'ı "doğrulayamayız", sadece ÇALIŞTIĞINI test edebiliriz '
            '(bkz. test_tbl1xr1_pipeline.py\'deki kesif-modu testi).'
        )

    rng = np.random.default_rng(seed)
    if genotip_artis_mi:
        ortalama = MASTROTOTARO_KO_ORTALAMA_HZ
        sd = hucre_arasi_varyasyon_sd or (MASTROTOTARO_KO_SEM * np.sqrt(MASTROTOTARO_N_HUCRE))
    else:
        ortalama = MASTROTOTARO_WT_ORTALAMA_HZ
        sd = hucre_arasi_varyasyon_sd or (MASTROTOTARO_WT_SEM * np.sqrt(MASTROTOTARO_N_HUCRE))

    spike_listesi = []
    for _ in range(n_unit):
        hucre_hizi = max(0.1, rng.normal(ortalama, sd))
        spike_listesi.append(poisson_spike_treni(hucre_hizi, sure_sn, rng))
    return spike_listesi
