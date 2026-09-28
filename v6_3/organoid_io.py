"""
organoid_io.py — NWB Dosya Okuma (Lazy Loading)
=================================================

Görev: NWB dosyalarını aç ve sadece veri akışını kur.
Yorum yok, analiz yok, sadece I/O.

Lazy loading: tüm datayı RAM'e almaz, chunk chunk okur.
"""

import pynwb
import numpy as np


def nwb_tip_tespit(dosya_yolu):
    """
    Detects the type of an NWB file.

    Returns:
        'ham'    — ham elektriksel sinyal var (acquisition['ES'])
        'sorted' — spike-sorted seti var (units tablosu)
        'her_ikisi' — ikisi birden var
        'unknown' — neither present
    """
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        # ElectricalSeries tipinde herhangi bir acquisition var mı?
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
        unit_sayisi, stim_var, trial_var

    [DÜZELTME — bkz. README] Bu fonksiyon önceden, doğru hesaplanan
    sr/n_sample/elektrot_sayisi/recording_suresi_sn değerlerini hemen
    ardından tekrar hesaplayıp None'a ezen ölü/kopya bir kod bloğu
    içeriyordu; ayrıca ElectricalSeries içermeyen (yalnızca sorted
    units barındıran) NWB dosyalarında gereksiz yere ValueError
    fırlatıyordu. Her iki sorun da bu fonksiyonu, ham sinyal içeren
    dosyalarda metadata'yı sessizce None'a düşürüyor, sorted-units-only
    dosyalarda ise anında çökertiyordu. Şimdi tek, doğru bir hesaplama
    bloğu var: ElectricalSeries bulunamazsa ilgili alanlar None döner,
    hata fırlatılmaz — çünkü sorted-units-only dosyalar geçerli ve sık
    karşılaşılan bir tiptir (bkz. nwb_tip_tespit'in 'sorted' dönüşü).
    """
    meta = {}
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        meta['session_description'] = str(nwb.session_description or '')
        meta['experiment_description'] = str(nwb.experiment_description or '')
        meta['session_start_time'] = str(nwb.session_start_time)[:19]

        # Find first ElectricalSeries acquisition object (name-independent)
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

        if nwb.units is not None:
            meta['unit_sayisi'] = len(nwb.units)
        else:
            meta['unit_sayisi'] = 0

        meta['stim_var'] = 'ES_STIM' in nwb.acquisition
        # [DÜZELTME — bkz. README] nwb.intervgets -> nwb.intervals
        # (pynwb'nin gerçek attribute adı 'intervals'; eski hâli her NWB
        # dosyasında AttributeError fırlatıyordu).
        meta['trial_var'] = (nwb.intervals is not None and
                              'trials' in nwb.intervals)

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
        # [DÜZELTME — bkz. README] değişken adı max_s -> max_sn
        # (tanımsız isim yüzünden bu hata mesajının kendisi NameError
        # fırlatıp asıl RAM hatasını gizliyordu).
        max_sn = MAX_CHUNK_BYTES / (sr * n_channels * 4)
        raise ValueError(
            f'Chunk RAM tahmini {tahmini_byte/1e6:.1f} MB > '
            f'limit {MAX_CHUNK_BYTES/1e6:.0f} MB. '
            f'chunk_sn × n_channels × sr × 4 limiti aşıyor. '
            f'Bu yapılandırma için max chunk_sn = {max_sn:.1f} s '
            f'(n_channels={n_channels}, sr={sr}). '
            f'Çözüm: chunk_sn azalt veya tek channel kullan.')
    return tahmini_byte


def ham_chunk_uret(dosya_yolu, chunk_sn=60, baslangic_sn=5, channel=0,
                    max_sure_sn=None):
    """
    Generator: ham sinyali chunk chunk üretir.

    Tüm datayı RAM'e almaz. Her chunk işlendikten sonra
    bellek serbest bırakılabilir.

    KRİTİK: Dosya tek bir with bloğunda açılır, generator boyunca
    açık kalır. Eski sürümde her chunk için open/close çalışıyordu.

    GÜVENLİK: chunk_sn × n_channels × sr × 4 byte hesaplanır,
    MAX_CHUNK_BYTES (500 MB) aşılırsa ValueError. Bu sayede
    256 kanallı MEA + 300 sn gibi RAM patlatıcı kombinasyonlar engellenir.

    Parametreler:
        dosya_yolu   : str
        chunk_sn     : int — chunk boyutu (saniye)
        baslangic_sn : int — recording başından atlanacak süre (artefakt için)
        channel      : int veya 'al' — okunacak channel indeksi
                       'al' ise tüm kanallar (M, C) shape ile döner
        max_sure_sn  : float veya None — verilirse, recording başından
                       (baslangic_sn sonrası) bu kadar süre işlenince
                       üretim durur. [DÜZELTME — bkz. README] Bu parametre
                       önceden yoktu; organoid_units_analiz.py'nin ham
                       sinyal modu bunu var sayıp geçirdiği için
                       TypeError fırlatıyordu. Şimdi gerçekten destekleniyor.

    Yield:
        (chunk_data, chunk_baslangic_sample, chunk_bitis_sample, sr)
        chunk_data: tek channel için 1D, 'al' için 2D (M, C)
    """
    # [DÜZELTME — bkz. README] değişken adı chunk_s -> chunk_sn
    # (tanımsız isim yüzünden bu doğrulama NameError fırlatıp asıl
    # "chunk_sn pozitif olmalı" hatasını gizliyordu).
    if chunk_sn <= 0:
        raise ValueError(f'chunk_sn pozitif olmalı, gelen: {chunk_sn}')

    # Tek seferlik dosya açılışı — generator boyunca açık kalır
    io = pynwb.NWBHDF5IO(dosya_yolu, 'r')
    try:
        nwb = io.read()
        # Find first ElectricalSeries acquisition object (name-independent)
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

        # RAM kontrolü
        n_channels_okunacak = n_channels_toplam if channel == 'al' else 1
        _chunk_byte_kontrol(chunk_sn, sr, n_channels_okunacak)

        chunk_n = int(chunk_sn * sr)
        baslangic = int(baslangic_sn * sr)
        bitis_sinir = n_sample
        if max_sure_sn is not None:
            bitis_sinir = min(n_sample, baslangic + int(max_sure_sn * sr))
        i = 0

        while baslangic + i * chunk_n < bitis_sinir:
            bas = baslangic + i * chunk_n
            bit = min(bitis_sinir, bas + chunk_n)
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
    Spike-sorted units setlerini okur.

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
    Trial bilgilerini okur (varsa).

    Returns list: her trial için dict
        {'start_time', 'stop_time', 's' (stim kodu)}
    """
    trials_list = []
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        # [DÜZELTME — bkz. README] nwb.intervgets -> nwb.intervals
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


def stim_channeli_searchlik(dosya_yolu):
    """
    ES_STIM kanalından stimülasyon başlangıç ve bitiş zamanlarını bulur.
    Bulamazsa (None, None) döner.

    Yorum yok — sadece eşik üstü ilk ve son örneği bulur.
    """
    with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
        nwb = io.read()
        if 'ES_STIM' not in nwb.acquisition:
            return None, None
        es_stim = nwb.acquisition['ES_STIM']
        sr_s = int(es_stim.rate)
        # Tüm stim datasını okumak şart, çünkü tek seferlik
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
