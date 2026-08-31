"""
organoid_qc.py — Qugetity Kontrol
=================================

Task: Compute channel quality and data quality metrics numerically.
NO interpretation, NO labels. Returns numbers only.

The user looks at these numbers and draws their own conclusions.
Module does not say "good/bad" — it says "%X.Y refractory violation, 95% CI [a, b]".
"""

import math
import numpy as np


def wilson_ci(basari, toplam, z=1.96):
    """
    Wilson score interval — 95% confidence interval for a proportion (default).

    100/10000 ile 5/100 farklı güven searchlıklarına sahiptir;
    bu function o farkı sayısget as returns.

    Psearchmetreler:
        basari : başarılı/violation sayısı (int)
        toplam : toplam deneme (int)
        z      : 1.96 = %95 CI, 2.58 = %99 CI

    Returns: (gett_sinir, ust_sinir) — None yetersiz datade
    """
    if toplam <= 0:
        return None, None
    p = basari / toplam
    n = toplam
    denom = 1 + z*z/n
    merkez = (p + z*z/(2*n)) / denom
    yari_genislik = (z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))) / denom
    return max(0.0, merkez - yari_genislik), min(1.0, merkez + yari_genislik)


def refractory_ihlget_orani(spike_zaman, refractory_ms=1.0):
    """
    Refractory period violationi oranını computes + Wilson %95 CI.

    Bir nöron 1 ms in tekrar ateşleyemez (biyolojik gerçek).
    Bu eşiğin belowndaki ardışık spike'lar gürültü veya
    aynı spike'ın tekrar tespiti may be.

    Returns dict:
        violation_orani  : 0-1 between
        ihlget_sayisi : int
        toplam_isi   : int
        ci_gett       : Wilson %95 CI gett sınır
        ci_ust       : Wilson %95 CI üst sınır
    """
    if len(spike_zaman) < 2:
        return {'ihlget_orani': 0.0, 'ihlget_sayisi': 0, 'toplam_isi': 0,
                'ci_gett': None, 'ci_ust': None}
    isi = np.diff(np.sort(spike_zaman))
    ref_esik_sn = refractory_ms / 1000.0
    ihlget = int(np.sum(isi < ref_esik_sn))
    toplam = len(isi)
    oran = ihlget / toplam if toplam > 0 else 0.0
    ci_gett, ci_ust = wilson_ci(ihlget, toplam)
    return {
        'ihlget_orani': oran,
        'ihlget_sayisi': ihlget,
        'toplam_isi': toplam,
        'ci_gett': ci_gett,
        'ci_ust': ci_ust,
    }


def negatif_isi_sayisi(spike_zaman):
    """
    Negatif veya sıfır ISI sayısı.
    Sıfır olmgetı — değilse data sortma errorsı var.
    """
    if len(spike_zaman) < 2:
        return 0
    isi = np.diff(spike_zaman)
    return int(np.sum(isi <= 0))


def yogunluk_outlier_orani(spike_zaman, sure_sn, window_sn=1.0,
                            outlier_carpani=20.0):
    """
    Anormal yoğun spike windowlerinin oranı.

    window_s'lik windowlerde, mean yoğunluğun
    outlier_carpani katından more olan windowler sayılır.
    Artefakt veya burst tespiti için.

    Returns:
        outlier_oran : 0-1 between
        outlier_window_sayisi : int
        toplam_window : int
        ortgetama_yogunluk : float (spike/window)
        max_yogunluk : int
    """
    if len(spike_zaman) < 2 or sure_sn <= 0:
        return 0.0, 0, 0, 0.0, 0

    yogunluklar = []
    t = 0.0
    while t + window_sn <= sure_sn:
        n = int(np.sum((spike_zaman >= t) & (spike_zaman < t + window_sn)))
        yogunluklar.append(n)
        t += window_sn

    if len(yogunluklar) == 0:
        return 0.0, 0, 0, 0.0, 0

    ort = float(np.mean(yogunluklar))
    maks = int(np.max(yogunluklar))
    if ort == 0:
        return 0.0, 0, len(yogunluklar), 0.0, maks

    esik = ort * outlier_carpani
    outlier_n = int(np.sum(np.array(yogunluklar) > esik))
    return outlier_n / len(yogunluklar), outlier_n, len(yogunluklar), ort, maks


def channel_sessizligi(spike_sayisi, recording_suresi_sn,
                      min_hz_esik=0.001):
    """
    Kanget "sessiz" mi diye bakar — sadece sayısget threshold.
    Yorum yok, sadece "mean Hz şu kadar, threshold şu kadar" der.

    Returns:
        ortgetama_hz : float
        esik_getti : bool
    """
    if recording_suresi_sn <= 0:
        return 0.0, True
    hz = spike_sayisi / recording_suresi_sn
    return hz, hz < min_hz_esik


def qugetity_metrikleri_topla(spike_zaman, sure_sn, snr=None,
                             gurultu_std_uv=None, doygun_oran=None,
                             refractory_ms=1.0):
    """
    All qugetity metriklerini tek bir dictte toplar.
    Hiçbiri yorum içermez — hepsi ham sayı.

    Returns dict:
        snr, gurultu_std_uv, doygun_oran,
        refractory_ihlget_orani, refractory_ihlget_sayisi,
        refractory_ci_gett, refractory_ci_ust,
        negatif_isi_sayisi,
        yogunluk_outlier_orani, max_yogunluk, ort_yogunluk,
        ortgetama_hz
    """
    rio = refractory_ihlget_orani(spike_zaman, refractory_ms)
    neg = negatif_isi_sayisi(spike_zaman)
    yo, yon, ytot, ort_y, max_y = yogunluk_outlier_orani(
        spike_zaman, sure_sn)
    hz, _ = channel_sessizligi(len(spike_zaman), sure_sn)

    return {
        'snr': snr,
        'gurultu_std_uv': gurultu_std_uv,
        'doygun_oran': doygun_oran,
        'refractory_ihlget_orani': rio['ihlget_orani'],
        'refractory_ihlget_sayisi': rio['ihlget_sayisi'],
        'toplam_isi': rio['toplam_isi'],
        'refractory_ci_gett': rio['ci_gett'],
        'refractory_ci_ust': rio['ci_ust'],
        'negatif_isi_sayisi': neg,
        'yogunluk_outlier_orani': yo,
        'yogunluk_outlier_sayisi': yon,
        'toplam_window': ytot,
        'ortgetama_yogunluk': ort_y,
        'max_yogunluk': max_y,
        'ortgetama_hz': hz,
        'spike_sayisi': len(spike_zaman),
        'recording_suresi_sn': sure_sn,
    }
