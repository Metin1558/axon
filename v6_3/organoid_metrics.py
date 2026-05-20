"""
organoid_metrics.py — Metrik ve İstatistik Hesaplama
======================================================

Görev: ISI, CV, FFT, Welch, ISI shuffle permütasyon,
       t-test, Mann-Whitney, ANOVA gibi tüm istatistiksel testler.

Yorum YOK. Sadece sayı. Anlamlı/anlamsız demez.

Tüm istatistiksel testler bu modülde toplanır — başka modül kendisi
istatistik yapmaz, buradan çağırır.
"""

import numpy as np
from scipy import stats
from scipy.signal import welch


# ─────────────────────────────────────────────────────────────
# TEMEL SPIKE METRİKLERİ
# ─────────────────────────────────────────────────────────────

def isi_metrikleri(spike_zaman):
    """
    Inter-spike interval temel istatistikleri.

    Dönen sözlük (yetersiz veride None'lar):
        ort_isi, std_isi, medyan_isi, cv,
        isi_q25, isi_q75, isi_min, isi_max
    """
    if len(spike_zaman) < 2:
        return {
            'ort_isi': None, 'std_isi': None, 'medyan_isi': None,
            'cv': None, 'isi_q25': None, 'isi_q75': None,
            'isi_min': None, 'isi_max': None,
        }
    isi = np.diff(np.sort(spike_zaman))
    ort = float(np.mean(isi))
    std = float(np.std(isi))
    return {
        'ort_isi': ort,
        'std_isi': std,
        'medyan_isi': float(np.median(isi)),
        'cv': std / ort if ort > 0 else None,
        'isi_q25': float(np.percentile(isi, 25)),
        'isi_q75': float(np.percentile(isi, 75)),
        'isi_min': float(np.min(isi)),
        'isi_max': float(np.max(isi)),
    }


def aktiflik_pencere(spike_zaman, sure_sn, pencere_sn=30,
                      min_spike=5):
    """
    Pencere bazlı aktiflik oranı.
    Pencere içinde min_spike kadar spike varsa "aktif" sayılır.

    Yorum yok — sadece "kaç pencere min eşiği geçti" der.

    Dönen değer:
        aktif_oran : 0-1 arası
        aktif_sayi : int
        toplam_pencere : int
    """
    if sure_sn < pencere_sn:
        return 0.0, 0, 0
    aktif = 0
    toplam = 0
    t = 0.0
    while t + pencere_sn <= sure_sn:
        n = int(np.sum((spike_zaman >= t) & (spike_zaman < t + pencere_sn)))
        if n >= min_spike:
            aktif += 1
        toplam += 1
        t += pencere_sn
    return aktif / toplam if toplam > 0 else 0.0, aktif, toplam


def cv_zaman_serisi(spike_zaman, sure_sn, pencere_sn=30):
    """
    Pencere bazlı CV zaman serisi.

    Dönen liste: her pencere için (zaman, cv, spike_sayisi)
    """
    seri = []
    t = 0.0
    while t + pencere_sn <= sure_sn:
        p = spike_zaman[(spike_zaman >= t) & (spike_zaman < t + pencere_sn)]
        if len(p) >= 2:
            isi = np.diff(p)
            if np.mean(isi) > 0:
                cv = float(np.std(isi) / np.mean(isi))
                seri.append((float(t), cv, int(len(p))))
        t += pencere_sn
    return seri


# ─────────────────────────────────────────────────────────────
# SPEKTRAL ANALİZ
# ─────────────────────────────────────────────────────────────

def fft_periyotlar(cv_serisi_degerleri, pencere_sn=30, top_n=3):
    """
    CV zaman serisi üzerinde FFT — dominant periyotları bulur.

    Dönen liste: [(periyot_sn, guc), ...] — top_n tanesi
    """
    if len(cv_serisi_degerleri) < 4:
        return []
    cv_arr = np.array(cv_serisi_degerleri)
    fft = np.abs(np.fft.rfft(cv_arr))
    freqler = np.fft.rfftfreq(len(cv_arr), d=pencere_sn)

    # En güçlü top_n+1 freq, ilki DC (0 Hz) hariç
    sirali = np.argsort(fft[1:])[-top_n:][::-1] + 1
    sonuc = []
    for i in sirali:
        if freqler[i] > 0:
            sonuc.append((float(1.0 / freqler[i]), float(fft[i])))
    return sonuc


def welch_periyotlar(cv_serisi_degerleri, pencere_sn=30, top_n=3):
    """
    Welch metodu — dominant periyotları bulur.

    Dönen liste: [(periyot_sn, psd), ...] — top_n tanesi
    """
    if len(cv_serisi_degerleri) < 4:
        return []
    cv_arr = np.array(cv_serisi_degerleri)
    nperseg = min(len(cv_arr), 8)
    noverlap = nperseg // 2
    try:
        freqler, psd = welch(cv_arr, fs=1.0/pencere_sn,
                              nperseg=nperseg, noverlap=noverlap)
    except Exception:
        return []
    sirali = np.argsort(psd[1:])[-top_n:][::-1] + 1
    sonuc = []
    for i in sirali:
        if i < len(freqler) and freqler[i] > 0:
            sonuc.append((float(1.0 / freqler[i]), float(psd[i])))
    return sonuc


def isi_shuffle_permutasyon(spike_zaman, sure_sn, pencere_sn=30,
                              n_permutasyon=1000, seed=42):
    """
    ISI shuffle permütasyon testi.

    Gerçek ISI'ları karıştırır (otokorelasyonu korur), her karıştırılmış
    sekansta dominant FFT gücünü hesaplar. Gerçek değer rastgele
    dağılımın neresinde duruyor — bunu sayısal olarak verir.

    Yorum yok — ham p değeri ve permütasyon sayısı döndürür.

    Dönen sözlük:
        p_deger        : float (raw, manipüle edilmemiş)
        permutasyon_n  : int
        gercek_guc     : float
        rastgele_guc_medyan : float
        rastgele_guc_q95    : float
        eslesen_permutasyon : int (kaç permütasyonda eşit veya büyük güç bulundu)
    """
    if len(spike_zaman) < 4:
        return None

    # Gerçek değer
    cv_seri = cv_zaman_serisi(spike_zaman, sure_sn, pencere_sn)
    if len(cv_seri) < 4:
        return None

    cv_arr = np.array([x[1] for x in cv_seri])
    fft = np.abs(np.fft.rfft(cv_arr))
    if len(fft) < 2:
        return None
    gercek_guc = float(np.max(fft[1:]))

    # Permütasyonlar
    isi = np.diff(np.sort(spike_zaman))
    if len(isi) < 2:
        return None

    rng = np.random.default_rng(seed)
    rastgele_gucler = []

    for _ in range(n_permutasyon):
        karisik = rng.permutation(isi)
        karisik_zaman = np.concatenate(
            ([spike_zaman[0]], spike_zaman[0] + np.cumsum(karisik)))
        karisik_zaman = karisik_zaman[karisik_zaman <= sure_sn]

        r_seri = cv_zaman_serisi(karisik_zaman, sure_sn, pencere_sn)
        if len(r_seri) < 4:
            continue
        r_cv = np.array([x[1] for x in r_seri])
        r_fft = np.abs(np.fft.rfft(r_cv))
        if len(r_fft) < 2:
            continue
        rastgele_gucler.append(float(np.max(r_fft[1:])))

    if len(rastgele_gucler) == 0:
        return None

    rastgele_arr = np.array(rastgele_gucler)
    eslesen = int(np.sum(rastgele_arr >= gercek_guc))
    p_deger = eslesen / len(rastgele_arr)

    return {
        'p_deger': p_deger,
        'permutasyon_n': len(rastgele_arr),
        'gercek_guc': gercek_guc,
        'rastgele_guc_medyan': float(np.median(rastgele_arr)),
        'rastgele_guc_q95': float(np.percentile(rastgele_arr, 95)),
        'eslesen_permutasyon': eslesen,
        'yontem': 'ISI shuffling',
        'seed': seed,
    }


# ─────────────────────────────────────────────────────────────
# İSTATİSTİKSEL TESTLER
# ─────────────────────────────────────────────────────────────

def shapiro_wilk(veri, max_n=5000):
    """
    Shapiro-Wilk normallik testi.

    Spike rate ve ISI verileri genelde log-normal dağılır.
    Bu test, t-test'in ön koşulunu kontrol eder.

    Dönen sözlük:
        w_stat, p_deger, n, normal_mi
        normal_mi : p > 0.05 ise True (yorum değil, eşik kontrolü)
                    Shapiro varsayılan eşiği bilim camiasında 0.05.

    Sınırlılık: n > 5000 için Shapiro yavaşlar, otomatik
    örnekleme yapılır (rastgele 5000 değer).
    """
    veri = np.asarray(veri, dtype=np.float64)
    veri = veri[np.isfinite(veri)]
    n = len(veri)
    if n < 3:
        return None
    if n > max_n:
        rng = np.random.default_rng(42)
        veri = rng.choice(veri, size=max_n, replace=False)
    try:
        w, p = stats.shapiro(veri)
    except Exception:
        return None
    return {
        'test': 'Shapiro-Wilk',
        'w_stat': float(w),
        'p_deger': float(p),
        'n': n,
        'normal_mi': bool(p > 0.05),
    }


def ttest_iki_grup(grup1, grup2, normalite_kontrol=True, welch=True):
    """
    Bağımsız iki grup t-testi.

    Varsayılan: Welch's t-test (equal_var=False) — eşit olmayan varyans
    varsayımı klasik t-test'ten daha sağlamdır. İki grubun varyansları
    farklı olabilir ve bu test bunu kabul eder.

    NOT: t-test verinin normal dağılım göstermesini varsayar.
    Spike rate / ISI gibi nörobilim verileri genelde log-normal
    dağılır. Bu fonksiyon Shapiro-Wilk sonucunu da döndürür —
    eğer normal_mi=False ise Mann-Whitney U tercih edilmeli.

    welch=False yaparsan klasik Student's t-test çalışır.

    Dönen sözlük:
        t_stat, p_deger, n1, n2, mean1, mean2, std1, std2,
        normalite_g1, normalite_g2, normalite_uyari, test_tipi
    """
    g1 = np.asarray(grup1, dtype=np.float64)
    g2 = np.asarray(grup2, dtype=np.float64)
    if len(g1) < 2 or len(g2) < 2:
        return None
    t_stat, p_val = stats.ttest_ind(g1, g2, equal_var=not welch)

    sonuc = {
        'test': "Welch's t-test" if welch else "Student's t-test",
        'test_tipi': 'welch' if welch else 'student',
        't_stat': float(t_stat),
        'p_deger': float(p_val),
        'n1': len(g1),
        'n2': len(g2),
        'mean1': float(np.mean(g1)),
        'mean2': float(np.mean(g2)),
        'std1': float(np.std(g1)),
        'std2': float(np.std(g2)),
    }

    if normalite_kontrol:
        sonuc['normalite_g1'] = shapiro_wilk(g1)
        sonuc['normalite_g2'] = shapiro_wilk(g2)
        n1_ok = sonuc['normalite_g1'] and sonuc['normalite_g1']['normal_mi']
        n2_ok = sonuc['normalite_g2'] and sonuc['normalite_g2']['normal_mi']
        sonuc['normalite_uyari'] = not (n1_ok and n2_ok)

    return sonuc


def mannwhitney_iki_grup(grup1, grup2):
    """
    Mann-Whitney U testi (non-parametrik).
    """
    g1 = np.asarray(grup1, dtype=np.float64)
    g2 = np.asarray(grup2, dtype=np.float64)
    if len(g1) < 2 or len(g2) < 2:
        return None
    u_stat, p_val = stats.mannwhitneyu(g1, g2, alternative='two-sided')
    return {
        'test': 'Mann-Whitney U',
        'u_stat': float(u_stat),
        'p_deger': float(p_val),
        'n1': len(g1),
        'n2': len(g2),
        'medyan1': float(np.median(g1)),
        'medyan2': float(np.median(g2)),
    }


def anova_n_grup(*gruplar):
    """
    Tek yönlü ANOVA. İki veya daha fazla grup.
    """
    g_listesi = [np.asarray(g, dtype=np.float64) for g in gruplar]
    if any(len(g) < 2 for g in g_listesi) or len(g_listesi) < 2:
        return None
    f_stat, p_val = stats.f_oneway(*g_listesi)
    return {
        'test': 'one-way ANOVA',
        'f_stat': float(f_stat),
        'p_deger': float(p_val),
        'n_gruplar': [len(g) for g in g_listesi],
        'mean_gruplar': [float(np.mean(g)) for g in g_listesi],
    }


def trial_stim_karsilastir(spike_zaman, trials):
    """
    Trial bazlı stim VAR vs YOK karşılaştırması.

    Bu örnekteki trial yapısı: trial['s'] iki elemanlı liste,
    [s1, s2] = [E1_aktif, E2_aktif]. (1,1) = stim var, (0,0) = stim yok.

    Dönen sözlük:
        ttest, mannwhitney, n_var, n_yok, mean_hz_var, mean_hz_yok
    """
    if not trials:
        return None

    hz_var = []
    hz_yok = []

    for t in trials:
        sure = t['stop_time'] - t['start_time']
        if sure <= 0:
            continue
        spk = spike_zaman[(spike_zaman >= t['start_time']) &
                          (spike_zaman < t['stop_time'])]
        hz = len(spk) / sure
        if 's' in t:
            s = t['s']
            if s[0] == 1 and len(s) > 1 and s[1] == 1:
                hz_var.append(hz)
            elif s[0] == 0 and len(s) > 1 and s[1] == 0:
                hz_yok.append(hz)

    if len(hz_var) < 2 or len(hz_yok) < 2:
        return None

    return {
        'ttest': ttest_iki_grup(hz_var, hz_yok),
        'mannwhitney': mannwhitney_iki_grup(hz_var, hz_yok),
        'n_var': len(hz_var),
        'n_yok': len(hz_yok),
        'mean_hz_var': float(np.mean(hz_var)),
        'mean_hz_yok': float(np.mean(hz_yok)),
        'std_hz_var': float(np.std(hz_var)),
        'std_hz_yok': float(np.std(hz_yok)),
    }


# ─────────────────────────────────────────────────────────────
# ÇOKLU KANAL SENKRON METRİKLERİ
# ─────────────────────────────────────────────────────────────

def sttc_iki_kanal(spike_t1, spike_t2, sure_sn, dt=0.01, sirali=True):
    """
    Spike-Time Tiling Coefficient (Cutts & Eglen 2014).

    İki spike treni arasındaki senkronizasyon — pencere bazlı
    cross-correlation'a göre çok daha sağlam (firing rate'e
    duyarsız, +1 ile -1 arası).

    Parametreler:
        spike_t1, spike_t2 : iki kanal/unit'in spike zamanları (sn)
        sure_sn            : kayıt süresi
        dt                 : tiling pencere yarı-genişliği (sn, varsayılan 10 ms)
        sirali             : True ise spike trenler zaten sıralı kabul edilir
                             (NWB ve sorting'den gelen veriler için)
                             False ise içeride sort yapılır

    Dönen değer: float (-1 ile 1 arası, None yetersiz veride)
    """
    if len(spike_t1) == 0 or len(spike_t2) == 0 or sure_sn <= 0:
        return None

    # Sort sadece gerekirse
    if not sirali:
        spike_t1 = np.sort(spike_t1)
        spike_t2 = np.sort(spike_t2)

    def TA(spike_train, sure, dt):
        """Spike train çevresindeki ±dt pencerelerin toplam oranı."""
        if len(spike_train) == 0:
            return 0.0
        # spike_train zaten sıralı (sirali=True veya yukarıda sıralandı)
        bas = np.maximum(spike_train - dt, 0)
        bit = np.minimum(spike_train + dt, sure)
        # Sıralı dizide bas zaten sıralı; örtüşmeleri lineer birleştir
        toplam = 0.0
        cur_bas = bas[0]
        cur_bit = bit[0]
        for i in range(1, len(bas)):
            if bas[i] <= cur_bit:
                cur_bit = max(cur_bit, bit[i])
            else:
                toplam += cur_bit - cur_bas
                cur_bas = bas[i]
                cur_bit = bit[i]
        toplam += cur_bit - cur_bas
        return toplam / sure

    def PA(spike_a, spike_b, dt):
        """A'nın spike'larından kaç tanesi B'nin ±dt penceresinde?"""
        if len(spike_a) == 0 or len(spike_b) == 0:
            return 0.0
        # spike_b zaten sıralı, searchsorted çalışır
        idx = np.searchsorted(spike_b, spike_a)
        sayac = 0
        for i, t in enumerate(spike_a):
            if idx[i] > 0 and abs(t - spike_b[idx[i]-1]) <= dt:
                sayac += 1
            elif idx[i] < len(spike_b) and abs(t - spike_b[idx[i]]) <= dt:
                sayac += 1
        return sayac / len(spike_a)

    TA1 = TA(spike_t1, sure_sn, dt)
    TA2 = TA(spike_t2, sure_sn, dt)
    PA12 = PA(spike_t1, spike_t2, dt)
    PA21 = PA(spike_t2, spike_t1, dt)

    if TA2 >= 1.0 or TA1 >= 1.0:
        return None

    pay1 = (PA12 - TA2) / (1 - PA12 * TA2) if (1 - PA12 * TA2) != 0 else 0
    pay2 = (PA21 - TA1) / (1 - PA21 * TA1) if (1 - PA21 * TA1) != 0 else 0

    return 0.5 * (pay1 + pay2)


def coklu_kanal_senkron(spike_treni_listesi, sure_sn, dt=0.01):
    """
    Tüm kanal/unit çiftleri için STTC matrisi.

    Parametreler:
        spike_treni_listesi : [np.array, np.array, ...] — her kanal/unit
        sure_sn             : kayıt süresi
        dt                  : STTC tiling penceresi

    Dönen sözlük:
        sttc_matrix : N x N simetrik matris (np.array)
        ortalama_sttc : tüm çiftlerin ortalama STTC'si
        max_sttc, min_sttc
        n_kanal
    """
    n = len(spike_treni_listesi)
    if n < 2:
        return None
    matris = np.full((n, n), np.nan)
    cift_sttc = []
    for i in range(n):
        matris[i, i] = 1.0
        for j in range(i+1, n):
            v = sttc_iki_kanal(spike_treni_listesi[i],
                                spike_treni_listesi[j],
                                sure_sn, dt)
            if v is not None:
                matris[i, j] = v
                matris[j, i] = v
                cift_sttc.append(v)
    if not cift_sttc:
        return None
    return {
        'sttc_matrix': matris,
        'ortalama_sttc': float(np.mean(cift_sttc)),
        'medyan_sttc': float(np.median(cift_sttc)),
        'max_sttc': float(np.max(cift_sttc)),
        'min_sttc': float(np.min(cift_sttc)),
        'n_kanal': n,
        'n_cift': len(cift_sttc),
        'dt_sn': dt,
    }


def stim_oncesi_sirasi_sonrasi(spike_zaman, sure_sn, stim_bas, stim_bit):
    """
    ES_STIM kanaalına göre stim öncesi/sırası/sonrası spike sayıları
    ve oranları.

    Yorum yok — etiket yok ("aktivasyon" demez), sadece sayı verir.

    Dönen sözlük:
        once_n, sirasinda_n, sonra_n,
        once_hz, sirasinda_hz, sonra_hz,
        sonra_once_orani (None yoksa)
    """
    if stim_bas is None or stim_bit is None:
        return None

    once = spike_zaman[spike_zaman < stim_bas]
    sira = spike_zaman[(spike_zaman >= stim_bas) & (spike_zaman < stim_bit)]
    sonra = spike_zaman[spike_zaman >= stim_bit]

    sure_once = stim_bas
    sure_sira = stim_bit - stim_bas
    sure_sonra = sure_sn - stim_bit

    once_hz = len(once) / sure_once if sure_once > 0 else None
    sira_hz = len(sira) / sure_sira if sure_sira > 0 else None
    sonra_hz = len(sonra) / sure_sonra if sure_sonra > 0 else None

    oran = None
    if once_hz is not None and once_hz > 0 and sonra_hz is not None:
        oran = sonra_hz / once_hz

    return {
        'stim_bas': stim_bas,
        'stim_bit': stim_bit,
        'once_n': len(once),
        'sirasinda_n': len(sira),
        'sonra_n': len(sonra),
        'once_hz': once_hz,
        'sirasinda_hz': sira_hz,
        'sonra_hz': sonra_hz,
        'sonra_once_orani': oran,
    }
