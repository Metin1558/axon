"""
organoid_signal.py — Sinyal Operasyonları
==========================================

Görev: Ham sinyali filtrele, spike tespit et, SNR hesapla.
Yorum yok — sadece sayısal işlemler.

Spike detection, zorunlu SNR kontrolünden geçer.
"""

import numpy as np
from scipy.signal import butter, filtfilt


def bandpass_filter(seti, sr, low_hz=300, high_hz=3000, derece=4):
    """
    Butterworth band-pass filter.

    Default 300-3000 Hz: organoid MEA kayıtları için standart aralık.
    Low-cut kaldırır: hareket artefaktları, baseline drift
    High-cut kaldırır: elektronik gürültü
    """
    nyq = sr / 2.0
    low_n = low_hz / nyq
    high_n = min(high_hz / nyq, 0.99)
    b, a = butter(derece, [low_n, high_n], btype='band')
    return filtfilt(b, a, seti).astype(np.float32)


def mad_std(seti):
    """
    Median Absolute Deviation tabanlı std tahmini.

    Quiroga 2004 yaklaşımı:
        sigma_hat = MAD / 0.6745

    Klasik np.std, yüksek aktiviteli kayıtlarda spike'ları gürültü
    sayar — MAD spike'lara duyarsızdır, daha robust bir tahmin verir.
    """
    seti = np.asarray(seti)
    if len(seti) == 0:
        return 0.0
    med = np.median(seti)
    mad = np.median(np.abs(seti - med))
    return float(mad / 0.6745)


def gurultu_std_ornekle(io_modul, dosya_yolu, sr, low_hz, high_hz,
                         n_ornek=10, ornek_sn=30, channel=0,
                         dc_cikar=False, mad_kullan=False):
    """
    Kayıt boyunca rastgele chunk'lar örnekleyerek gürültü std'sini tahmin eder.
    Tek chunk yerine chunk'lar arası median kullanır — outlier'lara robust.

    dc_cikar   : True ise median çıkarılır (yavaş ama DC-drift'li kayıtlar
                 için gerekli).
    mad_kullan : True ise MAD/0.6745 (Quiroga 2004), False ise klasik np.std.
                 Yüksek aktiviteli kayıtlar için True önerilir — spike'lar
                 gürültü tahminini şişirir; MAD bundan etkilenmez.

    [DÜZELTME — bkz. README] Bu fonksiyon önceden acquisition anahtarını
    sabit 'ES' string'i ile arıyordu (nwb.acquisition['ES']); gerçek NWB
    dosyalarında ElectricalSeries objesi farklı bir isimle kayıtlı olabilir,
    bu durumda KeyError fırlatıyordu. Şimdi organoid_io.py'deki ile aynı
    generic (isimden bağımsız) ElectricalSeries arama mantığını kullanıyor.
    """
    import pynwb
    rng = np.random.default_rng(42)
    std_listesi = []

    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        es = next(
            (obj for obj in nwb.acquisition.values()
             if hasattr(obj, 'data') and hasattr(obj, 'rate')),
            None
        )
        if es is None:
            raise ValueError('NWB dosyası ElectricalSeries içermiyor')
        n_sample = es.data.shape[0]
        chunk_n = int(ornek_sn * sr)
        baslangic = int(5 * sr)

        if n_sample - baslangic < chunk_n:
            data = es.data[baslangic:, channel].astype(np.float32)
            if dc_cikar:
                data = data - np.median(data)
            data_f = bandpass_filter(data, sr, low_hz, high_hz)
            return mad_std(data_f) if mad_kullan else float(np.std(data_f))

        for _ in range(n_ornek):
            bas = int(rng.integers(baslangic, n_sample - chunk_n))
            bit = bas + chunk_n
            data = es.data[bas:bit, channel].astype(np.float32)
            if dc_cikar:
                data = data - np.median(data)
            data_f = bandpass_filter(data, sr, low_hz, high_hz)
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

    Default: bidirectional detection (|signal| > sigma * std).
    Birincil spike çoğu kayıtta negatif olsa da, elektrot pozisyonuna ve
    organoid oryantasyonuna bağlı olarak pozitif spike'lar da bulunabilir.

    sadece_negatif=True : sadece negatif threshold (signal < -sigma * std).
                          Yüksek SNR'lı sorting öncesi MUA için
                          tercih edilebilir.

    Parametreler:
        chunk_data : 1D array — ham sinyal parçası
        sr         : sampling rate
        sigma_esik : kaç sigma threshold (default 5)
        gurultu_std: önceden hesaplanan gürültü std
        refractory_ms: minimum spike-arası süre (default 1 ms)
        dc_cikar   : True ise median çıkarılır (DC drift için, yavaş)
        sadece_negatif: True ise sadece negatif threshold (opt-in)

    Returns:
        spike_idx — chunk içindeki spike indeksleri
        chunk_filterli — filtrelenmiş data (SNR hesaplaması için)
    """
    if dc_cikar:
        chunk = chunk_data - np.median(chunk_data)
    else:
        chunk = chunk_data
    chunk_f = bandpass_filter(chunk, sr, low_hz, high_hz)

    if gurultu_std is None:
        gurultu_std = float(np.std(chunk_f))

    if gurultu_std == 0:
        return np.array([], dtype=np.int64), chunk_f

    if sadece_negatif:
        spike_mask = chunk_f < -sigma_esik * gurultu_std
    else:
        spike_mask = np.abs(chunk_f) > sigma_esik * gurultu_std

    spike_idx = np.where(spike_mask)[0]

    # Dinamik refractory — sample rate'e bağlı
    min_spike_araligi = max(1, int((refractory_ms / 1000.0) * sr))

    if len(spike_idx) > 1:
        diff = np.diff(spike_idx)
        spike_idx = spike_idx[np.concatenate(
            ([True], diff > min_spike_araligi))]

    return spike_idx.astype(np.int64), chunk_f


def snr_compute(chunk_filterli, spike_idx, sr, gurultu_std,
                window_ms=2, max_ornek=200):
    """
    Spike SNR'ını hesaplar — tepe genlik / gürültü std.

    Yorum yok, sayı döner.

    Returns:
        snr           : float
        ort_genlik_uv : float (uV cinsinden ortalama spike tepe genliği)
        n_spike       : int (kullanılan spike sayısı)
    """
    if len(spike_idx) == 0 or gurultu_std == 0:
        return 0.0, 0.0, 0

    window = int(sr * window_ms / 1000.0)
    n_kullan = min(max_ornek, len(spike_idx))
    genlikler = []

    for idx in spike_idx[:n_kullan]:
        bas = max(0, idx - window)
        bit = min(len(chunk_filterli), idx + window)
        if bit > bas:
            genlikler.append(float(np.max(np.abs(chunk_filterli[bas:bit]))))

    if len(genlikler) == 0:
        return 0.0, 0.0, 0

    ort_genlik = float(np.mean(genlikler))
    snr = ort_genlik / gurultu_std
    return snr, ort_genlik * 1e6, len(genlikler)


def doygunluk_tespit(chunk_data, doygunluk_orani_esik=0.05):
    """
    Sinyal doygunluğu tespit eder — ADC'nin tepe değerine yakın
    örneklerin oranı.

    Returns:
        doygun_oran : 0-1 arasında, doygun örnek oranı
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
