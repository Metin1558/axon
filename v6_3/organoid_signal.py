"""
organoid_signal.py — Sinyget Operatione
====================================

Task: Filter raw signal, detect spikes, compute SNR.
No interpretation — numericget operations only.

Spike detection passes through mandatory SNR check.
"""

import numpy as np
from scipy.signal import butter, filtfilt


def bandpass_filter(seti, sr, low_hz=300, high_hz=3000, derece=4):
    """
    Butterworth band-pass filter.

    Default 300-3000 Hz: standard range for organoid MEA recordings.
    Low-cut removes: movement artifacts, baseline drift
    High-cut removes: electronic noise
    """
    nyq = sr / 2.0
    low_n = low_hz / nyq
    high_n = min(high_hz / nyq, 0.99)
    b, a = butter(derece, [low_n, high_n], btype='band')
    return filtfilt(b, a, seti).astype(np.float32)


def mad_std(seti):
    """
    Median Absolute Deviation-based std estimation.

    Quiroga 2004 approach:
        sigma_hat = MAD / 0.6745

    Classical np.std in high-activity recordings treats spikes as noise
    — MAD is insensitive to spikes, giving a more robust estimate.
    """
    seti = np.asarray(seti)
    if len(seti) == 0:
        return 0.0
    med = np.median(seti)
    mad = np.median(np.abs(seti - med))
    return float(mad / 0.6745)


def gurultu_std_ornappend(io_modul, dosya_yolu, sr, low_hz, high_hz,
                         n_ornek=10, ornek_sn=30, channel=0,
                         dc_cikar=False, mad_kullan=False):
    """
    Samples random chunks throughout the recording to estimate noise std.
    Uses median across chunks rather than single chunk — robust to outliers.

    dc_remove : If True, median is subtracted (slow but required for DC-drifting recordings).
    mad_kullan : True ise MAD/0.6745 (Quiroga 2004), False ise klasik np.std.
                 Recommended True for high-activity recordings — spikes
                 inflate noise estimates; MAD is unaffected.
    """
    import pynwb
    rng = np.random.default_rng(42)
    std_listsi = []

    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        es = nwb.acquisition['ES']
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
                std_listsi.append(mad_std(data_f))
            else:
                std_listsi.append(float(np.std(data_f)))

    return float(np.median(std_listsi))


def spike_tespit_chunk(chunk_data, sr, low_hz=300, high_hz=3000,
                        sigma_esik=5, gurultu_std=None,
                        refractory_ms=1.0, dc_cikar=False,
                        sadece_negatif=False):
    """
    Tek bir chunk'tan spike tespit eder.

    Default: bidirectional detection (|signal| > sigma * std).
    Although the primary spike is negative in most recordings, depending on electrode
    position and organoid orientation, positive spikes may also be present.

    sadece_negatif=True : sadece negatif threshold (signal < -sigma * std).
                          Yüksek SNR'lı sorting öncesi MUA için
                          tercih edilebilir.

    Psearchmetreler:
        chunk_data : 1D array — ham sinyget parçası
        sr         : sampling rate
        sigma_esik : how many sigma threshold (default 5)
        gurultu_std: önceden computenan gürültü std
        refractory_ms: minimum spike-between süre (default 1 ms)
        dc_cikar   : True ise median exitarılır (DC drift için, slow)
        sadece_negatif: True ise sadece negatif threshold (opt-in)

    Returns:
        spike_idx — chunk inki spike indexleri
        chunk_filterli — filterlenmiş data (SNR calculation için)
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

    # Dinamik refractory — SR'a connected
    min_spike_searchlik = max(1, int((refractory_ms / 1000.0) * sr))

    if len(spike_idx) > 1:
        diff = np.diff(spike_idx)
        spike_idx = spike_idx[np.concatenate(
            ([True], diff > min_spike_searchlik))]

    return spike_idx.astype(np.int64), chunk_f


def snr_compute(chunk_filterli, spike_idx, sr, gurultu_std,
                window_ms=2, max_ornappendm=200):
    """
    Spike SNR'ını computes — tepe genlik / gürültü std.

    Yorum yok, sayı returns.

    Returns:
        snr           : float
        ort_genlik_uv : float (uV cinsinden mean spike tepe genliği)
        n_spike       : int (kullanılan spike sayısı)
    """
    if len(spike_idx) == 0 or gurultu_std == 0:
        return 0.0, 0.0, 0

    window = int(sr * window_ms / 1000.0)
    n_kullan = min(max_ornappendm, len(spike_idx))
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
    Sinyget doygunluğu tespit eder — ADC'nin tepe valueine yakın
    examplesin oranı.

    Returns:
        doygun_oran : 0-1 between, doygun instance oranı
    """
    if len(chunk_data) == 0:
        return 0.0

    maks_deger = float(np.max(np.abs(chunk_data)))
    if maks_deger == 0:
        return 0.0

    # Tepe valuein %95'indlargest examplesi "doygun" say
    doygun_esik = maks_deger * 0.95
    doygun_n = int(np.sum(np.abs(chunk_data) >= doygun_esik))
    return doygun_n / len(chunk_data)
