"""
organoid_signal.py — Sinyal İşleme
====================================

Görev: Ham sinyali filtrele, spike tespit et, SNR'ı hesapla.
Yorum yok — sadece sayısal işleme.

Spike tespiti zorunlu SNR kontrolünden geçer.
"""

import numpy as np
from scipy.signal import butter, filtfilt


def bandpass_filtre(veri, sr, low_hz=300, high_hz=3000, derece=4):
    """
    Butterworth band-pass filtre.

    Varsayılan 300-3000 Hz: organoid MEA kayıtları için standart aralık.
    Düşük: hareket artefaktı, baseline drift
    Yüksek: elektronik gürültü
    """
    nyq = sr / 2.0
    low_n = low_hz / nyq
    high_n = min(high_hz / nyq, 0.99)
    b, a = butter(derece, [low_n, high_n], btype='band')
    return filtfilt(b, a, veri).astype(np.float32)


def mad_std(veri):
    """
    Median Absolute Deviation tabanlı std tahmini.

    Quiroga 2004'ün yaklaşımı:
        sigma_hat = MAD / 0.6745

    Klasik np.std yüksek aktiviteli kayıtlarda spike'ları gürültü
    sayar — MAD spike'lara duyarsız, daha sağlam tahmin verir.
    """
    veri = np.asarray(veri)
    if len(veri) == 0:
        return 0.0
    med = np.median(veri)
    mad = np.median(np.abs(veri - med))
    return float(mad / 0.6745)


def gurultu_std_ornekle(io_modul, dosya_yolu, sr, low_hz, high_hz,
                         n_ornek=10, ornek_sn=30, kanal=0,
                         dc_cikar=False, mad_kullan=False):
    """
    Kayıt boyunca rastgele örneklerden gürültü standart sapması hesaplar.
    Tek bir chunk yerine medyan kullanılır — outlier'lara karşı robust.

    dc_cikar : True ise medyan çıkarılır (yavaş ama DC kayan kayıtlarda gerekli).
    mad_kullan : True ise MAD/0.6745 (Quiroga 2004), False ise klasik np.std.
                 Yüksek aktiviteli kayıtlarda True önerilir — spike'lar
                 gürültü tahminini şişirir, MAD bundan etkilenmez.
    """
    import pynwb
    rng = np.random.default_rng(42)
    std_listesi = []

    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        es = nwb.acquisition['ES']
        n_sample = es.data.shape[0]
        chunk_n = int(ornek_sn * sr)
        baslangic = int(5 * sr)

        if n_sample - baslangic < chunk_n:
            data = es.data[baslangic:, kanal].astype(np.float32)
            if dc_cikar:
                data = data - np.median(data)
            data_f = bandpass_filtre(data, sr, low_hz, high_hz)
            return mad_std(data_f) if mad_kullan else float(np.std(data_f))

        for _ in range(n_ornek):
            bas = int(rng.integers(baslangic, n_sample - chunk_n))
            bit = bas + chunk_n
            data = es.data[bas:bit, kanal].astype(np.float32)
            if dc_cikar:
                data = data - np.median(data)
            data_f = bandpass_filtre(data, sr, low_hz, high_hz)
            if mad_kullan:
                std_listesi.append(mad_std(data_f))
            else:
                std_listesi.append(float(np.std(data_f)))

    return float(np.median(std_listesi))


def spike_tespit_chunk(chunk_data, sr, low_hz=300, high_hz=3000,
                        sigma_esik=5, gurultu_std=None,
                        refractory_ms=1.0, dc_cikar=False,
                        sadece_negatif=False):
    """
    Tek bir chunk'tan spike tespit eder.

    Default: çift yönlü tespit (|signal| > sigma * std).
    Çoğu kayıtta ana spike negatif olsa da, elektrot pozisyonuna ve
    organoid yönelimine bağlı olarak pozitif spike'lar da görülebilir.

    sadece_negatif=True : sadece negatif eşik (signal < -sigma * std).
                          Yüksek SNR'lı sorting öncesi MUA için
                          tercih edilebilir.

    Parametreler:
        chunk_data : 1D array — ham sinyal parçası
        sr         : sampling rate
        sigma_esik : kaç sigma eşik (varsayılan 5)
        gurultu_std: önceden hesaplanan gürültü std
        refractory_ms: minimum spike-arası süre (varsayılan 1 ms)
        dc_cikar   : True ise medyan çıkarılır (DC drift için, yavaş)
        sadece_negatif: True ise sadece negatif eşik (opt-in)

    Dönen değer:
        spike_idx — chunk içindeki spike indeksleri
        chunk_filtreli — filtrelenmiş veri (SNR hesabı için)
    """
    if dc_cikar:
        chunk = chunk_data - np.median(chunk_data)
    else:
        chunk = chunk_data
    chunk_f = bandpass_filtre(chunk, sr, low_hz, high_hz)

    if gurultu_std is None:
        gurultu_std = float(np.std(chunk_f))

    if gurultu_std == 0:
        return np.array([], dtype=np.int64), chunk_f

    if sadece_negatif:
        spike_mask = chunk_f < -sigma_esik * gurultu_std
    else:
        spike_mask = np.abs(chunk_f) > sigma_esik * gurultu_std

    spike_idx = np.where(spike_mask)[0]

    # Dinamik refractory — SR'a bağlı
    min_spike_aralik = max(1, int((refractory_ms / 1000.0) * sr))

    if len(spike_idx) > 1:
        diff = np.diff(spike_idx)
        spike_idx = spike_idx[np.concatenate(
            ([True], diff > min_spike_aralik))]

    return spike_idx.astype(np.int64), chunk_f


def snr_hesapla(chunk_filtreli, spike_idx, sr, gurultu_std,
                pencere_ms=2, max_orneklem=200):
    """
    Spike SNR'ını hesaplar — tepe genlik / gürültü std.

    Yorum yok, sayı döndürür.

    Dönen değer:
        snr           : float
        ort_genlik_uv : float (uV cinsinden ortalama spike tepe genliği)
        n_spike       : int (kullanılan spike sayısı)
    """
    if len(spike_idx) == 0 or gurultu_std == 0:
        return 0.0, 0.0, 0

    pencere = int(sr * pencere_ms / 1000.0)
    n_kullan = min(max_orneklem, len(spike_idx))
    genlikler = []

    for idx in spike_idx[:n_kullan]:
        bas = max(0, idx - pencere)
        bit = min(len(chunk_filtreli), idx + pencere)
        if bit > bas:
            genlikler.append(float(np.max(np.abs(chunk_filtreli[bas:bit]))))

    if len(genlikler) == 0:
        return 0.0, 0.0, 0

    ort_genlik = float(np.mean(genlikler))
    snr = ort_genlik / gurultu_std
    return snr, ort_genlik * 1e6, len(genlikler)


def doygunluk_tespit(chunk_data, doygunluk_orani_esik=0.05):
    """
    Sinyal doygunluğu tespit eder — ADC'nin tepe değerine yakın
    örneklerin oranı.

    Dönen değer:
        doygun_oran : 0-1 arası, doygun örnek oranı
    """
    if len(chunk_data) == 0:
        return 0.0

    maks_deger = float(np.max(np.abs(chunk_data)))
    if maks_deger == 0:
        return 0.0

    # Tepe değerin %95'inden büyük örnekleri "doygun" say
    doygun_esik = maks_deger * 0.95
    doygun_n = int(np.sum(np.abs(chunk_data) >= doygun_esik))
    return doygun_n / len(chunk_data)
