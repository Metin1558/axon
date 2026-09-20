"""
organoid_sorting.py — Spike Sorting Sarmalı
=============================================

Görev: SpikeInterface üzerinden Tridesclous2 ile spike sorting.
Opsiyonel — sadece --sort flag ile çağırılır.

Bu modül sadece sarmalı yapar. Sonuçları organoid_metrics.py
analiz eder. Burada hiçbir analiz veya yorum yapılmaz.
"""

import os
import numpy as np


def spike_sort_calistir(nwb_dosya, cikis_klasoru, sigma_esik=5,
                         freq_min=300, freq_max=3000):
    """
    SpikeInterface + Tridesclous2 ile spike sorting çalıştırır.

    Parametreler:
        nwb_dosya       : NWB dosya yolu
        cikis_klasoru   : sorting çıktıları için klasör
        sigma_esik      : detection sigma eşiği
        freq_min/max    : band-pass filtre sınırları

    Dönen sözlük:
        units      : list of dicts — {unit_id, spike_zaman, n_spike}
        n_unit     : int
        sr         : sampling rate
        kayit_suresi : float

    Hata olursa None döndürür ve nedeni stdout'a yazar.
    """
    try:
        import spikeinterface.full as si
        from spikeinterface.sorters import run_sorter
    except ImportError:
        print('SpikeInterface kurulu degil. pip install spikeinterface')
        return None

    try:
        # NWB recording extractor
        recording = si.read_nwb_recording(nwb_dosya)
        sr = recording.sampling_frequency

        # Preprocess
        recording_f = si.bandpass_filter(recording,
                                          freq_min=freq_min,
                                          freq_max=freq_max)
        recording_cmr = si.common_reference(recording_f,
                                              reference='global',
                                              operator='median')

        # Sort
        os.makedirs(cikis_klasoru, exist_ok=True)
        sorting = run_sorter(
            sorter_name='tridesclous2',
            recording=recording_cmr,
            output_folder=os.path.join(cikis_klasoru, 'tdc2_output'),
            detect_threshold=sigma_esik,
            remove_existing_folder=True,
        )

        unit_ids = sorting.get_unit_ids()
        units_list = []
        for uid in unit_ids:
            spike_train = sorting.get_unit_spike_train(uid)
            spike_zaman = spike_train / sr
            units_list.append({
                'unit_id': int(uid),
                'spike_zaman': spike_zaman.astype(np.float64),
                'spike_sayisi': len(spike_train),
            })

        return {
            'units': units_list,
            'n_unit': len(units_list),
            'sr': float(sr),
            'kayit_suresi': float(recording.get_total_duration()),
            'sorter': 'tridesclous2',
            'parametreler': {
                'sigma_esik': sigma_esik,
                'freq_min': freq_min,
                'freq_max': freq_max,
            },
        }

    except Exception as e:
        print(f'Spike sorting hatasi: {e}')
        return None


def units_birlestir(units_list):
    """
    Tüm unit'lerin spike zamanlarını birleştirir, sıralar.

    Dönen değer: tek bir np.array — tüm spike'lar zamanı sıralı.
    """
    if not units_list:
        return np.array([], dtype=np.float64)
    tum_spikeler = []
    for u in units_list:
        tum_spikeler.extend(u['spike_zaman'].tolist())
    return np.sort(np.array(tum_spikeler, dtype=np.float64))
