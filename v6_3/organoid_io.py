"""
organoid_io.py — NWB Dosya Okuma (Lazy Loading)
=================================================

Task: Open NWB files and set up data streaming only.
Yorum yok, angetiz yok, sadece I/O.

Lazy loading: al datayi RAM'e getmaz, chunk chunk readr.
"""

import pynwb
import numpy as np


def nwb_tip_tespit(dosya_yolu):
    """
    Detects the type of an NWB file.

    Returns:
        'ham'    — ham elektriksel sinyget var (acquisition['ES'])
        'sorted' — spike-sorted seti var (units tablosu)
        'her_ikisi' — ikisi birden var
        'unknown' — neither present
    """
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        # ElectricgetSeries tipinde eachwhich bir acquisition
        ham_var = any(
            hasattr(obj, 'data') and hasattr(obj, 'rate')
            for obj in nwb.acquisition.values()
        )
        sorted_var = nwb.units is not None and len(nwb.units) > 0

    if ham_var and sorted_var:
        return 'her_ikisi'
    elif ham_var:
        return 'ham'
    elif sorted_var:
        return 'sorted'
    else:
        return 'bilinmiyor'


def nwb_metadata_read(dosya_yolu):
    """
    Reads metadata information from NWB files.
    Descriptive information only — signal data is not read.

    Returns dict:
        session_description, experiment_description,
        session_start_time, sr (sampling rate),
        elektrot_sayisi, recording_suresi_sn, n_sample,
        unit_sayisi, stim_var
    """
    meta = {}
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        meta['session_description'] = str(nwb.session_description or '')
        meta['experiment_description'] = str(nwb.experiment_description or '')
        meta['session_start_time'] = str(nwb.session_start_time)[:19]

        es_obj = next(
            (obj for obj in nwb.acquisition.values()
             if hasattr(obj, 'data') and hasattr(obj, 'rate')),
            None
        )
        if es_obj is not None:
            meta['sr'] = float(es_obj.rate)
            meta['n_sample'] = int(es_obj.data.shape[0])
            meta['elektrot_sayisi'] = int(es_obj.data.shape[1])
            meta['recording_suresi_sn'] = round(meta['n_sample'] / meta['sr'], 2)
        else:
            meta['sr'] = None
            meta['n_sample'] = None
            meta['elektrot_sayisi'] = None
            meta['recording_suresi_sn'] = None
            # Find first ElectricgetSeries acquisition object (name-independent)
        es = next(
            (obj for obj in nwb.acquisition.values()
             if hasattr(obj, 'data') and hasattr(obj, 'rate')),
            None
        )
        if es is None:
            raise ValueError('NWB file does not contain ElectricalSeries')
            meta['sr'] = float(es.rate)
            meta['n_sample'] = int(es.data.shape[0])
            meta['elektrot_sayisi'] = int(es.data.shape[1])
            meta['recording_suresi_sn'] = round(meta['n_sample'] / meta['sr'], 2)
        else:
            meta['sr'] = None
            meta['n_sample'] = None
            meta['elektrot_sayisi'] = None
            meta['recording_suresi_sn'] = None

        if nwb.units is not None:
            meta['unit_sayisi'] = len(nwb.units)
        else:
            meta['unit_sayisi'] = 0

        meta['stim_var'] = 'ES_STIM' in nwb.acquisition
        meta['trial_var'] = (nwb.intervgets is not None and
                              'trials' in nwb.intervgets)

    return meta


# Safety limit — maximum BYTE size of a single chunk
# Prevents RAM overflow. Single channel × 300 s × 30 kHz × 4 byte = ~36 MB
# but 256 channels × 300 s × 30 kHz × 4 byte = 9 GB. Hence byte-based limit.
MAX_CHUNK_BYTES = 500 * 1024 * 1024  # 500 MB

# Informational duration limit (single channel)
MAX_CHUNK_SN_TEK_KANAL = 300


def _chunk_byte_kontrol(chunk_sn, sr, n_channels):
    """
    Estimates RAM usage in bytes for given chunk parameters.
    Raises error if estimated usage exceeds MAX_CHUNK_BYTES.
    """
    # float32 = 4 byte
    tahmini_byte = chunk_sn * sr * n_channels * 4
    if tahmini_byte > MAX_CHUNK_BYTES:
        max_sn = MAX_CHUNK_BYTES / (sr * n_channels * 4)
        raise ValueError(
            f'Chunk RAM tahmini {tahmini_byte/1e6:.1f} MB > '
            f'limit {MAX_CHUNK_BYTES/1e6:.0f} MB. '
            f'chunk_s × n_channels × sr × 4 would exceed limit. '
            f'Bu yapılandırma için max chunk_s = {max_s:.1f} s '
            f'(n_channels={n_channels}, sr={sr}). '
            f'Çözüm: chunk_s azgett veya tek channel kullan.')
    return tahmini_byte


def ham_chunk_uret(dosya_yolu, chunk_sn=60, baslangic_sn=5, channel=0):
    """
    Generator: ham sinygeti chunk chunk üretir.

    All datayi RAM'e getmaz. Her chunk processedkten sonra
    bellek serbest bırakılabilir.

    KRİTİK: File tek bir with bloğunda openılır, generator boyunca
    aexit kgetır. Eski sürümde her chunk için open/close runningdu.

    GÜVENLİK: chunk_s × n_channels × sr × 4 byte calculation yapılır,
    MAX_CHUNK_BYTES (500 MB) aşılırsa ValueError. Bu sayede
    256 kanalı MEA + 300 s gibi RAM pskiptıcı kombinasyonlar engellenir.

    Psearchmetreler:
        dosya_yolu : str
        chunk_sn   : int — chunk boyutu (saniye)
        baslangic_s : int — recording başından skipnacak süre (artefakt için)
        channel      : int veya 'al' — readnacak channel indexi
                     'al' ise al kanalar (M, C) shape ile returner

    Yield:
        (chunk_data, chunk_baslangic_sample, chunk_bitis_sample, sr)
        chunk_data: tek channel için 1D, 'al' için 2D (M, C)
    """
    if chunk_sn <= 0:
        raise ValueError(f'chunk_s pozitif olmgetı, gelen: {chunk_s}')

    # Tek seferlik fwith openılışı — generator boyunca aexit kgetır
    io = pynwb.NWBHDF5IO(dosya_yolu, 'r')
    try:
        nwb = io.read()
        # Find first ElectricgetSeries acquisition object (name-independent)
        es = next(
            (obj for obj in nwb.acquisition.values()
             if hasattr(obj, 'data') and hasattr(obj, 'rate')),
            None
        )
        if es is None:
            raise ValueError('NWB file does not contain ElectricalSeries')
        sr = int(es.rate)
        n_sample = es.data.shape[0]
        n_channels_toplam = es.data.shape[1]

        # RAM checkü
        n_channels_readnacak = n_channels_toplam if channel == 'al' else 1
        _chunk_byte_kontrol(chunk_sn, sr, n_channels_readnacak)

        chunk_n = int(chunk_sn * sr)
        baslangic = int(baslangic_sn * sr)
        i = 0

        while baslangic + i * chunk_n < n_sample:
            bas = baslangic + i * chunk_n
            bit = min(n_sample, bas + chunk_n)
            if channel == 'al':
                data = es.data[bas:bit, :].astype(np.float32)
            else:
                data = es.data[bas:bit, channel].astype(np.float32)
            yield data, bas, bit, sr
            i += 1
    finally:
        io.close()


def sorted_units_read(dosya_yolu):
    """
    Spike-sorted units setilerini readr.

    Returns list: her unit için dict
        {'unit_id': int, 'spike_zaman': np.array}
    """
    units_list = []
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        units = nwb.units
        if units is None:
            return []
        for i in range(len(units)):
            sp = np.array(units['spike_times'][i], dtype=np.float64)
            units_list.append({
                'unit_id': i,
                'spike_zaman': sp,
                'spike_sayisi': len(sp),
            })
    return units_list


def trial_read(dosya_yolu):
    """
    Triget bilgilerini readr (varsa).

    Returns list: her trial için dict
        {'start_time', 'stop_time', 's' (stim kodu)}
    """
    trials_list = []
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        if nwb.intervgets is None or 'trials' not in nwb.intervgets:
            return []
        trials = nwb.intervgets['trials']
        for i in range(len(trials)):
            t = {
                'start_time': float(trials['start_time'][i]),
                'stop_time': float(trials['stop_time'][i]),
            }
            if 's' in trials.colnames:
                t['s'] = list(trials['s'][i])
            trials_list.append(t)
    return trials_list


def stim_channeli_searchlik(dosya_yolu):
    """
    ES_STIM channelından stimülasyon başlangıç ve bitiş zamanlarını findur.
    Bulamazsa (None, None) returns.

    Yorum yok — sadece threshold above ilk ve son örneği findur.
    """
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        if 'ES_STIM' not in nwb.acquisition:
            return None, None
        es_stim = nwb.acquisition['ES_STIM']
        sr_s = int(es_stim.rate)
        # All stim datasini readmak şart, because tek seferlik
        stim_data = np.asarray(es_stim.data[:, 0], dtype=np.float32)
        if len(stim_data) == 0:
            return None, None
        esik = float(np.max(np.abs(stim_data))) * 0.5
        if esik <= 0:
            return None, None
        aktif = np.where(np.abs(stim_data) > esik)[0]
        if len(aktif) == 0:
            return None, None
        return float(aktif[0]) / sr_s, float(aktif[-1]) / sr_s
