"""
tbl1xr1_metrics.py — TBL1XR1'e Özel İstatistik Katmanı
============================================================

TASARIM KARARI (axon-ndd'den farkı): Genel çerçevede burst_rate
"birincil metrik" gibi ele alınıp sonra bunun yanlış olduğu
sentetik testte ortaya çıkmıştı. Burada, Mastrototaro 2021'in
GERÇEK bulgusundan (sEPSC frekansı, hücre-seviyesi) yola çıkarak
BAŞTAN doğru metriği seçiyoruz: birim-başı ortalama ateşleme hızı
(ortalama_unit_ateşleme_hizi), popülasyon burst sayısı değil.

Network burst sayısı İKİNCİL/keşifsel bir metrik olarak hâlâ
hesaplanır (organoid_units_analiz.network_burst_tespit ile,
dokunulmadan) ama TBL1XR1 hipotezini TEST ETMEK için kullanılmaz.
"""

import numpy as np
from scipy import stats
from tbl1xr1_types import Tbl1xr1Genotip, profil_getir


def ortalama_unit_ateşleme_hizi(spike_listesi, sure_sn):
    """
    BİRİNCİL METRİK. Her unit'in kendi ateşleme hızını (Hz) hesaplar,
    sonra unit'ler arası ORTALAMAYI döner — Mastrototaro'nun sEPSC
    frekansı ölçümünün (hücre-başına, WT 4.0 Hz vs KO 6.3 Hz) en
    yakın MEA karşılığı budur.

    NOT: sEPSC bir SİNAPTİK AKIM ölçümü (postsinaptik potansiyel
    sıklığı), MEA ise EKSTRASELÜLER SPIKE ölçümü — bu ikisi aynı
    şey değil, sadece "birim-başı olay sıklığı" kavramsal olarak
    en yakın eşleşme. Bu çeviri varsayımının kendisi
    doğrulanmamıştır (bkz. tbl1xr1_types.PROFIL_LOF.notlar).

    Dönen: dict {ortalama_hz, medyan_hz, unit_bazli_hizlar, n_unit}
    """
    n_unit = len(spike_listesi)
    if n_unit == 0 or sure_sn <= 0:
        return None
    unit_hizlari = np.array([len(sp) / sure_sn for sp in spike_listesi])
    return {
        'ortalama_hz': float(np.mean(unit_hizlari)),
        'medyan_hz': float(np.median(unit_hizlari)),
        'unit_bazli_hizlar': unit_hizlari.tolist(),
        'n_unit': n_unit,
    }


def wilcoxon_izogenik_karsilastirma(mutant_degerler, kontrol_degerler):
    """
    İzogenik çift karşılaştırması — Wilcoxon signed-rank.
    (axon-ndd'deki organoid_metrics_ext.wilcoxon_paired_iki_grup ile
    aynı algoritma, burada TBL1XR1 bağlamına özel isimlendirmeyle
    yeniden yazıldı — bağımsız çalışabilsin diye.)
    """
    m = np.asarray(mutant_degerler, dtype=np.float64)
    k = np.asarray(kontrol_degerler, dtype=np.float64)

    if len(m) != len(k):
        raise ValueError(
            f'Eşleştirilmiş test için n eşit olmalı: mutant n={len(m)}, kontrol n={len(k)}'
        )
    if len(m) < 5:
        return None

    fark = m - k
    if np.all(fark == 0):
        return {'test': 'Wilcoxon signed-rank', 'stat': 0.0, 'p_deger': 1.0,
                'n_cift': len(m), 'medyan_fark': 0.0, 'yon': 'fark yok'}

    stat, p = stats.wilcoxon(m, k)
    return {
        'test': 'Wilcoxon signed-rank', 'stat': float(stat), 'p_deger': float(p),
        'n_cift': len(m), 'medyan_fark': float(np.median(fark)),
        'yon': 'mutant > kontrol' if np.median(fark) > 0
               else 'mutant < kontrol' if np.median(fark) < 0 else 'fark yok',
    }


def sure_confound_kontrolu(grup1_sureler, grup2_sureler, alpha=0.05):
    """axon-ndd'deki ile aynı mantık — TBL1XR1 pipeline'ında da zorunlu adım."""
    g1 = np.asarray(grup1_sureler, dtype=np.float64)
    g2 = np.asarray(grup2_sureler, dtype=np.float64)
    if len(g1) < 2 or len(g2) < 2:
        return {'confound_riski_var': None, 'uyari': 'Yetersiz veri (n<2).'}
    stat, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
    confound_var = bool(p < alpha)
    return {
        'confound_riski_var': confound_var,
        'p_deger': float(p),
        'uyari': (f'UYARI: kayıt-süresi farkı da anlamlı (p={p:.4f}) — genotip '
                  f'etkisi süre etkisinden ayrıştırılamaz.') if confound_var else None,
    }


def genotip_hipotezini_test_et(genotip: Tbl1xr1Genotip, mutant_hizlari, kontrol_hizlari):
    """
    ANA GİRİŞ NOKTASI. tbl1xr1_types.py'deki literatür-temelli profili
    ile gerçek ölçülen yönü karşılaştırır.

    BILINMEYEN_MISSENSE ya da KONTROL için hipotez testi YAPILMAZ —
    profil_getir zaten bunu engelliyor (BILINMEYEN_MISSENSE hata
    fırlatır); burada da açıkça yakalayıp "keşif modu" sonucu
    döndürüyoruz, sessizce geçmiyoruz.
    """
    if genotip in (Tbl1xr1Genotip.BILINMEYEN_MISSENSE, Tbl1xr1Genotip.KONTROL):
        return {
            'mod': 'kesif',
            'not': ('Bu genotip için literatürden hipotez yok. Sadece '
                    'tanımlayıcı istatistik + gerçek yön raporlanıyor, '
                    'onay/red iddiası YAPILMIYOR.'),
            'olculen_yon': ('artis' if np.median(mutant_hizlari) > np.median(kontrol_hizlari)
                             else 'azalis' if np.median(mutant_hizlari) < np.median(kontrol_hizlari)
                             else 'fark_yok'),
        }

    profil = profil_getir(genotip)
    wilcoxon_sonuc = wilcoxon_izogenik_karsilastirma(mutant_hizlari, kontrol_hizlari)

    olculen_yon = ('artis' if np.median(mutant_hizlari) > np.median(kontrol_hizlari)
                   else 'azalis' if np.median(mutant_hizlari) < np.median(kontrol_hizlari)
                   else 'fark_yok')

    beklenen = profil.beklenen_hucre_seviyesi_ateşleme_yonu

    if beklenen is None or beklenen == 'belirsiz':
        uyum = 'imza_belirsiz_karsilastirilamiyor'
    elif olculen_yon == beklenen:
        uyum = 'imza_ile_uyumlu'
    else:
        uyum = 'imza_ile_uyumsuz'

    return {
        'mod': 'hipotez_testi',
        'genotip': genotip.value,
        'beklenen_yon': beklenen,
        'olculen_yon': olculen_yon,
        'uyum': uyum,
        'kanit_seviyesi': profil.kanit_seviyesi.value,
        'wilcoxon': wilcoxon_sonuc,
        'kaynak': profil.kaynak,
        'ceviri_uyarisi': profil.notlar,
    }
