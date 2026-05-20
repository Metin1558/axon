"""
organoid_io.py — NWB Dosya Okuma (Lazy Loading)
=================================================

Görev: Sadece NWB dosyalarını açmak ve veri akışı sağlamak.
Yorum yok, analiz yok, sadece I/O.

Lazy loading: tüm veriyi RAM'e almaz, chunk chunk okur.
"""

import pynwb
import numpy as np


def nwb_tip_tespit(dosya_yolu):
    """
    NWB dosyasının tipini tespit eder.

    Dönen değerler:
        'ham'    — ham elektriksel sinyal var (acquisition['ES'])
        'sorted' — spike-sorted veri var (units tablosu)
        'her_ikisi' — ikisi birden var
        'bilinmiyor' — hiçbiri yok
    """
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        # ElectricalSeries tipinde herhangi bir acquisition
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


def nwb_metadata_oku(dosya_yolu):
    """
    NWB dosyasından metadata bilgilerini okur.
    Sadece tanımlayıcı bilgiler — sinyal verisi okunmaz.

    Dönen sözlük:
        session_description, experiment_description,
        session_start_time, sr (sampling rate),
        elektrot_sayisi, kayit_suresi_sn, n_sample,
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
            meta['kayit_suresi_sn'] = round(meta['n_sample'] / meta['sr'], 2)
        else:
            meta['sr'] = None
            meta['n_sample'] = None
            meta['elektrot_sayisi'] = None
            meta['kayit_suresi_sn'] = None
            # İlk ElectricalSeries acquisition nesnesini bul (isimden bağımsız)
        es = next(
            (obj for obj in nwb.acquisition.values()
             if hasattr(obj, 'data') and hasattr(obj, 'rate')),
            None
        )
        if es is None:
            raise ValueError('NWB dosyasında ElectricalSeries bulunamadı')
            meta['sr'] = float(es.rate)
            meta['n_sample'] = int(es.data.shape[0])
            meta['elektrot_sayisi'] = int(es.data.shape[1])
            meta['kayit_suresi_sn'] = round(meta['n_sample'] / meta['sr'], 2)
        else:
            meta['sr'] = None
            meta['n_sample'] = None
            meta['elektrot_sayisi'] = None
            meta['kayit_suresi_sn'] = None

        if nwb.units is not None:
            meta['unit_sayisi'] = len(nwb.units)
        else:
            meta['unit_sayisi'] = 0

        meta['stim_var'] = 'ES_STIM' in nwb.acquisition
        meta['trial_var'] = (nwb.intervals is not None and
                              'trials' in nwb.intervals)

    return meta


# Güvenlik sınırı — bir chunk'ın maksimum BYTE boyutu
# RAM patlamasını engelliyor. Tek kanal × 300 sn × 30 kHz × 4 byte = ~36 MB
# ama 256 kanal × 300 sn × 30 kHz × 4 byte = 9 GB. Bu yüzden byte tabanlı.
MAX_CHUNK_BYTES = 500 * 1024 * 1024  # 500 MB

# Bilgi amaçlı saniye sınırı (tek kanal için)
MAX_CHUNK_SN_TEK_KANAL = 300


def _chunk_byte_kontrol(chunk_sn, sr, n_kanal):
    """
    Verilen chunk parametrelerinin byte cinsinden RAM kullanımı.
    MAX_CHUNK_BYTES sınırını aşarsa hata fırlat.
    """
    # float32 = 4 byte
    tahmini_byte = chunk_sn * sr * n_kanal * 4
    if tahmini_byte > MAX_CHUNK_BYTES:
        max_sn = MAX_CHUNK_BYTES / (sr * n_kanal * 4)
        raise ValueError(
            f'Chunk RAM tahmini {tahmini_byte/1e6:.1f} MB > '
            f'limit {MAX_CHUNK_BYTES/1e6:.0f} MB. '
            f'chunk_sn × n_kanal × sr × 4 hesabı patlıyor. '
            f'Bu yapılandırma için max chunk_sn = {max_sn:.1f} sn '
            f'(n_kanal={n_kanal}, sr={sr}). '
            f'Çözüm: chunk_sn azalt veya tek kanal kullan.')
    return tahmini_byte


def ham_chunk_uret(dosya_yolu, chunk_sn=60, baslangic_sn=5, kanal=0):
    """
    Generator: ham sinyali chunk chunk üretir.

    Tüm veriyi RAM'e almaz. Her chunk işlendikten sonra
    bellek serbest bırakılabilir.

    KRİTİK: Dosya tek bir with bloğunda açılır, generator boyunca
    açık kalır. Eski sürümde her chunk için open/close yapılıyordu.

    GÜVENLİK: chunk_sn × n_kanal × sr × 4 byte hesabı yapılır,
    MAX_CHUNK_BYTES (500 MB) aşılırsa ValueError. Bu sayede
    256 kanallı MEA + 300 sn gibi RAM patlatıcı kombinasyonlar engellenir.

    Parametreler:
        dosya_yolu : str
        chunk_sn   : int — chunk boyutu (saniye)
        baslangic_sn : int — kayıt başından atlanacak süre (artefakt için)
        kanal      : int veya 'all' — okunacak kanal indeksi
                     'all' ise tüm kanallar (M, C) shape ile döner

    Yield:
        (chunk_data, chunk_baslangic_sample, chunk_bitis_sample, sr)
        chunk_data: tek kanal için 1D, 'all' için 2D (M, C)
    """
    if chunk_sn <= 0:
        raise ValueError(f'chunk_sn pozitif olmalı, gelen: {chunk_sn}')

    # Tek seferlik dosya açılışı — generator boyunca açık kalır
    io = pynwb.NWBHDF5IO(dosya_yolu, 'r')
    try:
        nwb = io.read()
        # İlk ElectricalSeries acquisition nesnesini bul (isimden bağımsız)
        es = next(
            (obj for obj in nwb.acquisition.values()
             if hasattr(obj, 'data') and hasattr(obj, 'rate')),
            None
        )
        if es is None:
            raise ValueError('NWB dosyasında ElectricalSeries bulunamadı')
        sr = int(es.rate)
        n_sample = es.data.shape[0]
        n_kanal_toplam = es.data.shape[1]

        # RAM kontrolü
        n_kanal_okunacak = n_kanal_toplam if kanal == 'all' else 1
        _chunk_byte_kontrol(chunk_sn, sr, n_kanal_okunacak)

        chunk_n = int(chunk_sn * sr)
        baslangic = int(baslangic_sn * sr)
        i = 0

        while baslangic + i * chunk_n < n_sample:
            bas = baslangic + i * chunk_n
            bit = min(n_sample, bas + chunk_n)
            if kanal == 'all':
                data = es.data[bas:bit, :].astype(np.float32)
            else:
                data = es.data[bas:bit, kanal].astype(np.float32)
            yield data, bas, bit, sr
            i += 1
    finally:
        io.close()


def sorted_units_oku(dosya_yolu):
    """
    Spike-sorted units verilerini okur.

    Dönen liste: her unit için sözlük
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


def trial_oku(dosya_yolu):
    """
    Trial bilgilerini okur (varsa).

    Dönen liste: her trial için sözlük
        {'start_time', 'stop_time', 's' (stim kodu)}
    """
    trials_list = []
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        if nwb.intervals is None or 'trials' not in nwb.intervals:
            return []
        trials = nwb.intervals['trials']
        for i in range(len(trials)):
            t = {
                'start_time': float(trials['start_time'][i]),
                'stop_time': float(trials['stop_time'][i]),
            }
            if 's' in trials.colnames:
                t['s'] = list(trials['s'][i])
            trials_list.append(t)
    return trials_list


def stim_kanali_aralik(dosya_yolu):
    """
    ES_STIM kanalından stimülasyon başlangıç ve bitiş zamanlarını bulur.
    Bulamazsa (None, None) döndürür.

    Yorum yok — sadece eşik üstü ilk ve son örneği bulur.
    """
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        if 'ES_STIM' not in nwb.acquisition:
            return None, None
        es_stim = nwb.acquisition['ES_STIM']
        sr_s = int(es_stim.rate)
        # Tüm stim verisini okumak şart, çünkü tek seferlik
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
