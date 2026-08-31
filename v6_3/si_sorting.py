"""
si_sorting.py — Axon v6.3
==========================
SpikeInterface spike sorting modulu.
Kanal secimi + kalite metrikleri dahil.

Kullanim:
    import si_sorting
    result = si_sorting.interaktif_sort(dosya_yolu)
    # result: {
    #   'spike_dict': {unit_id: spike_times_s},
    #   'quality': {unit_id: {snr, isi_violation_ratio, firing_rate, ...}},
    #   'n_units': int,
    #   'n_kaliteli': int,
    # }
"""

import numpy as np
from pathlib import Path

DESTEKLENEN_ALGORITMALAR = [
    'mountainsort5',
    'spykingcircus2',
    'tridesclous2',
    'simple',
]

# Kalite eşikleri (ayarlanabilir)
KALITE_ESIK = {
    'snr_min': 3.0,               # minimum SNR
    'isi_violation_max': 0.10,    # max ISI violation ratio (%10)
    'firing_rate_min': 0.05,      # min ateşleme hızı (Hz)
}

# Kanal seçimi
KANAL_SECIM_N = 64        # kaç aktif kanal kullanılacak
KANAL_RMS_PENCERE = 10.0  # RMS hesabı için ilk N saniye


# ── Kanal seçimi ──────────────────────────────────────────

def aktif_kanallar_sec(recording, n_kanal=KANAL_SECIM_N,
                       pencere_sn=KANAL_RMS_PENCERE, verbose=True):
    """
    En yüksek RMS amplitüdlü N kanalı seçer.
    Gürültülü veya sessiz kanalları dışarıda bırakır.
    """
    sr = recording.get_sampling_frequency()
    n_toplam = recording.get_num_channels()

    if n_toplam <= n_kanal:
        if verbose:
            print(f"  Kanal seçimi: tüm {n_toplam} kanal kullanılıyor "
                  f"(<= {n_kanal})")
        return recording

    # İlk pencere_sn saniyelik veriyi oku
    n_ornek = min(int(pencere_sn * sr), recording.get_num_samples())
    try:
        traces = recording.get_traces(start_frame=0, end_frame=n_ornek)
    except Exception as e:
        if verbose:
            print(f"  Kanal seçimi atlandı (veri okunamadı): {e}")
        return recording

    # DC offset çıkar, sonra RMS hesapla
    traces_f = traces.astype(np.float32)
    dc_offset = np.mean(traces_f, axis=0)
    traces_centered = traces_f - dc_offset
    rms = np.sqrt(np.mean(traces_centered ** 2, axis=0))
    # Tüm nan/inf'leri sıfırla
    rms = np.nan_to_num(rms, nan=0.0, posinf=0.0, neginf=0.0)

    # En yüksek RMS'li N kanal
    top_idx = np.argsort(rms)[::-1][:n_kanal]
    top_idx_sorted = sorted(top_idx.tolist())

    if verbose:
        print(f"  Kanal seçimi: {n_toplam} → {n_kanal} kanal "
              f"(RMS bazlı, {pencere_sn:.0f}s pencere)")
        print(f"  RMS aralığı seçili: {rms[top_idx].min():.2f} – "
              f"{rms[top_idx].max():.2f} μV")

    # Kanalları seç
    channel_ids = recording.get_channel_ids()
    selected_ids = [channel_ids[i] for i in top_idx_sorted]
    return recording.select_channels(selected_ids)


# ── Kalite metrikleri ────────────────────────────────────

def kalite_metrikleri_hesapla(sorting, recording_f, sr, verbose=True):
    """
    Her unit için kalite metrikleri hesaplar.
    
    Döndürür:
        dict: {unit_id: {snr, isi_violation_ratio, firing_rate,
                         n_spikes, kaliteli}}
    """
    try:
        import spikeinterface.full as si
        from spikeinterface.postprocessing import compute_principal_components
        from spikeinterface.qualitymetrics import compute_quality_metrics
    except ImportError:
        if verbose:
            print("  Kalite metrikleri: spikeinterface eksik, atlanıyor.")
        return {}

    if verbose:
        print("  Kalite metrikleri hesaplanıyor...")

    try:
        # Waveform extractor
        we = si.create_sorting_analyzer(
            sorting, recording_f,
            format='memory',
            sparse=True,
        )
        we.compute('random_spikes')
        we.compute('waveforms')
        we.compute('templates')
        we.compute('noise_levels')

        # Kalite metrikleri
        qm = compute_quality_metrics(
            we,
            metric_names=['snr', 'isi_violation', 'firing_rate'],
            verbose=False
        )

        sonuc = {}
        for uid in sorting.get_unit_ids():
            spike_train = sorting.get_unit_spike_train(uid, segment_index=0)
            n_spikes = len(spike_train)
            firing_rate = n_spikes / (recording_f.get_num_samples() / sr)

            snr = float(qm.loc[uid, 'snr']) if 'snr' in qm.columns else None
            isi_viol = float(qm.loc[uid, 'isi_violations_ratio']) \
                if 'isi_violations_ratio' in qm.columns else None

            kaliteli = True
            if snr is not None and snr < KALITE_ESIK['snr_min']:
                kaliteli = False
            if isi_viol is not None and isi_viol > KALITE_ESIK['isi_violation_max']:
                kaliteli = False
            if firing_rate < KALITE_ESIK['firing_rate_min']:
                kaliteli = False

            sonuc[int(uid)] = {
                'snr': round(snr, 3) if snr is not None else None,
                'isi_violation_ratio': round(isi_viol, 4) if isi_viol is not None else None,
                'firing_rate': round(firing_rate, 4),
                'n_spikes': n_spikes,
                'kaliteli': kaliteli,
            }

        n_kaliteli = sum(1 for v in sonuc.values() if v['kaliteli'])
        if verbose:
            print(f"  Kalite: {n_kaliteli}/{len(sonuc)} unit kalite eşiğini geçti")
            print(f"  Eşikler: SNR>{KALITE_ESIK['snr_min']}, "
                  f"ISI<{KALITE_ESIK['isi_violation_max']:.0%}, "
                  f"Hz>{KALITE_ESIK['firing_rate_min']}")

        return sonuc

    except Exception as e:
        if verbose:
            print(f"  Kalite metrikleri hatası: {e}")
        # Fallback: sadece temel metrikler
        sonuc = {}
        for uid in sorting.get_unit_ids():
            spike_train = sorting.get_unit_spike_train(uid, segment_index=0)
            n_spikes = len(spike_train)
            firing_rate = n_spikes / (recording_f.get_num_samples() / sr)
            kaliteli = firing_rate >= KALITE_ESIK['firing_rate_min']
            sonuc[int(uid)] = {
                'snr': None,
                'isi_violation_ratio': None,
                'firing_rate': round(firing_rate, 4),
                'n_spikes': n_spikes,
                'kaliteli': kaliteli,
            }
        return sonuc


# ── Ana sorting fonksiyonu ────────────────────────────────

def sort_nwb(dosya_yolu, algorithm='mountainsort5',
             max_sure_sn=180.0, freq_min=300, freq_max=3000,
             n_kanal=KANAL_SECIM_N, verbose=True):
    """
    NWB dosyasındaki ham sinyali SpikeInterface ile spike sort eder.
    Kanal seçimi ve kalite metrikleri dahil.

    Returns:
        dict: {
            'spike_dict': {unit_id: spike_times_s},
            'quality': {unit_id: {snr, isi_violation_ratio, firing_rate, ...}},
            'n_units': int,
            'n_kaliteli': int,
            'sr': float,
        }
        Boş dict: başarısız
    """
    try:
        import spikeinterface.full as si
        from spikeinterface.preprocessing import (bandpass_filter,
                                                   unsigned_to_signed)
    except ImportError:
        print("  HATA: spikeinterface yüklü değil.")
        print("  Kur: pip install spikeinterface[full]")
        return {}

    dosya_yolu = str(dosya_yolu)

    if verbose:
        print(f"  SpikeInterface v{si.__version__}")
        print(f"  Algoritma: {algorithm}")
        print(f"  Dosya: {Path(dosya_yolu).name}")

    # 1. NWB okuma
    try:
        recording = si.read_nwb(dosya_yolu)
    except Exception:
        try:
            import pynwb
            with pynwb.NWBHDF5IO(dosya_yolu, 'r') as io:
                nwb = io.read()
                es_name = next(
                    (name for name, obj in nwb.acquisition.items()
                     if hasattr(obj, 'data') and hasattr(obj, 'rate')),
                    None
                )
            if es_name is None:
                print("  HATA: ElectricalSeries bulunamadı.")
                return {}
            recording = si.read_nwb(dosya_yolu,
                                     electrical_series_path=es_name)
        except Exception as e:
            print(f"  HATA: NWB okunamadı: {e}")
            return {}

    sr = recording.get_sampling_frequency()
    n_ch = recording.get_num_channels()
    n_s = recording.get_num_samples()
    sure = n_s / sr

    if verbose:
        print(f"  Sampling rate : {sr:.0f} Hz")
        print(f"  Kanal sayısı  : {n_ch}")
        print(f"  Kayıt süresi  : {sure:.1f} s")

    # 2. Süre sınırla
    n_max = int(max_sure_sn * sr)
    if n_s > n_max:
        recording = recording.frame_slice(0, n_max)
        if verbose:
            print(f"  Kayıt kısaltıldı: {max_sure_sn:.0f} s")

    # 3. Unsigned → signed
    try:
        recording = unsigned_to_signed(recording)
    except Exception:
        pass

    # 4. Kanal seçimi (YENİ)
    recording = aktif_kanallar_sec(recording, n_kanal=n_kanal, verbose=verbose)

    # 5. Bandpass
    if verbose:
        print(f"  Bandpass: {freq_min}–{freq_max} Hz...")
    recording_f = bandpass_filter(recording,
                                   freq_min=freq_min, freq_max=freq_max)

    # 6. Sorting
    if verbose:
        print(f"  Sorting başlıyor ({algorithm})...")

    import shutil, os
    output_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        f'{algorithm}_output'
    )
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)

    try:
        if algorithm == 'mountainsort5':
            sorting = si.run_sorter('mountainsort5', recording_f,
                                     verbose=False,
                                     scheme='1',
                                     detect_threshold=5.0)
        elif algorithm == 'spykingcircus2':
            sorting = si.run_sorter('spykingcircus2', recording_f,
                                     verbose=False)
        elif algorithm == 'tridesclous2':
            sorting = si.run_sorter('tridesclous2', recording_f,
                                     verbose=False)
        elif algorithm == 'simple':
            sorting = si.run_sorter('simple', recording_f,
                                     verbose=False)
        else:
            print(f"  HATA: Bilinmeyen algoritma: {algorithm}")
            return {}
    except Exception as e:
        print(f"  HATA: Sorting başarısız: {e}")
        return {}

    unit_ids = sorting.get_unit_ids()
    if verbose:
        print(f"  Bulunan unit sayısı: {len(unit_ids)}")

    # 7. Spike times
    spike_dict = {}
    for uid in unit_ids:
        frames = sorting.get_unit_spike_train(uid, segment_index=0)
        times_s = frames / sr
        spike_dict[int(uid)] = times_s

    # 8. Kalite metrikleri (YENİ)
    quality = kalite_metrikleri_hesapla(sorting, recording_f, sr, verbose)

    n_kaliteli = sum(1 for v in quality.values() if v.get('kaliteli', True))

    return {
        'spike_dict': spike_dict,
        'quality': quality,
        'n_units': len(spike_dict),
        'n_kaliteli': n_kaliteli,
        'sr': sr,
    }


def interaktif_sort(dosya_yolu, verbose=True):
    """İnteraktif algoritma seçimi."""
    print()
    print("  SpikeInterface Spike Sorting")
    print("  " + "-"*40)
    for i, alg in enumerate(DESTEKLENEN_ALGORITMALAR, 1):
        print(f"  {i} : {alg}")
    print("  0 : İptal")
    print()

    # Kanal sayısı seçimi
    print(f"  Kanal seçimi: varsayılan {KANAL_SECIM_N} "
          f"(RMS bazlı en aktif kanallar)")
    try:
        n_kanal_str = input(
            f"  Kanal sayısı [Enter={KANAL_SECIM_N}, 0=hepsi]: ").strip()
        if n_kanal_str == '':
            n_kanal = KANAL_SECIM_N
        elif n_kanal_str == '0':
            n_kanal = 999999
        else:
            n_kanal = int(n_kanal_str)
    except (ValueError, EOFError):
        n_kanal = KANAL_SECIM_N

    print()
    try:
        secim = input("  Algoritma seçin (0-4): ").strip()
    except (EOFError, KeyboardInterrupt):
        return {}

    if secim == '0' or not secim:
        return {}

    try:
        idx = int(secim) - 1
        if 0 <= idx < len(DESTEKLENEN_ALGORITMALAR):
            algorithm = DESTEKLENEN_ALGORITMALAR[idx]
        else:
            print("  Geçersiz seçim.")
            return {}
    except ValueError:
        print("  Geçersiz seçim.")
        return {}

    return sort_nwb(dosya_yolu, algorithm=algorithm,
                    n_kanal=n_kanal, verbose=verbose)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        result = sort_nwb(sys.argv[1])
        if result:
            print(f"\nSonuç: {result['n_units']} unit "
                  f"({result['n_kaliteli']} kaliteli)")
            for uid, times in list(result['spike_dict'].items())[:5]:
                q = result['quality'].get(uid, {})
                snr_s = f"SNR={q.get('snr','N/A')}"
                print(f"  Unit {uid}: {len(times)} spike, "
                      f"Hz={q.get('firing_rate','N/A')}, {snr_s}, "
                      f"kaliteli={q.get('kaliteli','?')}")
    else:
        print("Kullanım: python si_sorting.py <dosya.nwb>")
