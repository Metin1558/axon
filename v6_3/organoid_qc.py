"""
organoid_qc.py — Kalite Kontrol
=================================

Görev: Kanal sağlığı ve veri kalitesi metriklerini sayısal olarak hesaplar.
Yorum YOK, etiket YOK. Sadece sayı döndürür.

Kullanıcı bu sayılara bakıp kendi kararını verir.
Modül "iyi/kötü" demez — "%X.Y refractory ihlali, %95 CI [a, b]" der.
"""

import math
import numpy as np


def wilson_ci(basari, toplam, z=1.96):
    """
    Wilson score interval — bir oran için %95 güven aralığı (varsayılan).

    100/10000 ile 5/100 farklı güven aralıklarına sahiptir;
    bu fonksiyon o farkı sayısal olarak verir.

    Parametreler:
        basari : başarılı/ihlal sayısı (int)
        toplam : toplam deneme (int)
        z      : 1.96 = %95 CI, 2.58 = %99 CI

    Dönen: (alt_sinir, ust_sinir) — None yetersiz veride
    """
    if toplam <= 0:
        return None, None
    p = basari / toplam
    n = toplam
    denom = 1 + z*z/n
    merkez = (p + z*z/(2*n)) / denom
    yari_genislik = (z * math.sqrt(p*(1-p)/n + z*z/(4*n*n))) / denom
    return max(0.0, merkez - yari_genislik), min(1.0, merkez + yari_genislik)


def refractory_ihlal_orani(spike_zaman, refractory_ms=1.0):
    """
    Refractory period ihlali oranını hesaplar + Wilson %95 CI.

    Bir nöron 1 ms içinde tekrar ateşleyemez (biyolojik gerçek).
    Bu eşiğin altındaki ardışık spike'lar gürültü veya
    aynı spike'ın tekrar tespiti olabilir.

    Dönen sözlük:
        ihlal_orani  : 0-1 arası
        ihlal_sayisi : int
        toplam_isi   : int
        ci_alt       : Wilson %95 CI alt sınır
        ci_ust       : Wilson %95 CI üst sınır
    """
    if len(spike_zaman) < 2:
        return {'ihlal_orani': 0.0, 'ihlal_sayisi': 0, 'toplam_isi': 0,
                'ci_alt': None, 'ci_ust': None}
    isi = np.diff(np.sort(spike_zaman))
    ref_esik_sn = refractory_ms / 1000.0
    ihlal = int(np.sum(isi < ref_esik_sn))
    toplam = len(isi)
    oran = ihlal / toplam if toplam > 0 else 0.0
    ci_alt, ci_ust = wilson_ci(ihlal, toplam)
    return {
        'ihlal_orani': oran,
        'ihlal_sayisi': ihlal,
        'toplam_isi': toplam,
        'ci_alt': ci_alt,
        'ci_ust': ci_ust,
    }


def negatif_isi_sayisi(spike_zaman):
    """
    Negatif veya sıfır ISI sayısı.
    Sıfır olmalı — değilse veri sıralama hatası var.
    """
    if len(spike_zaman) < 2:
        return 0
    isi = np.diff(spike_zaman)
    return int(np.sum(isi <= 0))


def yogunluk_outlier_orani(spike_zaman, sure_sn, pencere_sn=1.0,
                            outlier_carpani=20.0):
    """
    Anormal yoğun spike pencerelerinin oranı.

    pencere_sn'lik pencerelerde, ortalama yoğunluğun
    outlier_carpani katından fazla olan pencereler sayılır.
    Artefakt veya patlama tespiti için.

    Dönen değer:
        outlier_oran : 0-1 arası
        outlier_pencere_sayisi : int
        toplam_pencere : int
        ortalama_yogunluk : float (spike/pencere)
        max_yogunluk : int
    """
    if len(spike_zaman) < 2 or sure_sn <= 0:
        return 0.0, 0, 0, 0.0, 0

    yogunluklar = []
    t = 0.0
    while t + pencere_sn <= sure_sn:
        n = int(np.sum((spike_zaman >= t) & (spike_zaman < t + pencere_sn)))
        yogunluklar.append(n)
        t += pencere_sn

    if len(yogunluklar) == 0:
        return 0.0, 0, 0, 0.0, 0

    ort = float(np.mean(yogunluklar))
    maks = int(np.max(yogunluklar))
    if ort == 0:
        return 0.0, 0, len(yogunluklar), 0.0, maks

    esik = ort * outlier_carpani
    outlier_n = int(np.sum(np.array(yogunluklar) > esik))
    return outlier_n / len(yogunluklar), outlier_n, len(yogunluklar), ort, maks


def kanal_sessizligi(spike_sayisi, kayit_suresi_sn,
                      min_hz_esik=0.001):
    """
    Kanal "sessiz" mi diye bakar — sadece sayısal eşik.
    Yorum yok, sadece "ortalama Hz şu kadar, eşik şu kadar" der.

    Dönen değer:
        ortalama_hz : float
        esik_alti : bool
    """
    if kayit_suresi_sn <= 0:
        return 0.0, True
    hz = spike_sayisi / kayit_suresi_sn
    return hz, hz < min_hz_esik


def kalite_metrikleri_topla(spike_zaman, sure_sn, snr=None,
                             gurultu_std_uv=None, doygun_oran=None,
                             refractory_ms=1.0):
    """
    Tüm kalite metriklerini tek bir sözlükte toplar.
    Hiçbiri yorum içermez — hepsi ham sayı.

    Dönen sözlük:
        snr, gurultu_std_uv, doygun_oran,
        refractory_ihlal_orani, refractory_ihlal_sayisi,
        refractory_ci_alt, refractory_ci_ust,
        negatif_isi_sayisi,
        yogunluk_outlier_orani, max_yogunluk, ort_yogunluk,
        ortalama_hz
    """
    rio = refractory_ihlal_orani(spike_zaman, refractory_ms)
    neg = negatif_isi_sayisi(spike_zaman)
    yo, yon, ytot, ort_y, max_y = yogunluk_outlier_orani(
        spike_zaman, sure_sn)
    hz, _ = kanal_sessizligi(len(spike_zaman), sure_sn)

    return {
        'snr': snr,
        'gurultu_std_uv': gurultu_std_uv,
        'doygun_oran': doygun_oran,
        'refractory_ihlal_orani': rio['ihlal_orani'],
        'refractory_ihlal_sayisi': rio['ihlal_sayisi'],
        'toplam_isi': rio['toplam_isi'],
        'refractory_ci_alt': rio['ci_alt'],
        'refractory_ci_ust': rio['ci_ust'],
        'negatif_isi_sayisi': neg,
        'yogunluk_outlier_orani': yo,
        'yogunluk_outlier_sayisi': yon,
        'toplam_pencere': ytot,
        'ortalama_yogunluk': ort_y,
        'max_yogunluk': max_y,
        'ortalama_hz': hz,
        'spike_sayisi': len(spike_zaman),
        'kayit_suresi_sn': sure_sn,
    }
