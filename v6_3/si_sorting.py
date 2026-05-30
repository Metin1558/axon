"""
si_sorting.py — SpikeInterface spike sorting modulu
Axon v6.3 ile entegre calisir.

Kullanim:
    import si_sorting
    spike_dict = si_sorting.sort_nwb(dosya_yolu, algorithm='spykingcircus2')
    # spike_dict: {unit_id: np.array(spike_times_seconds)}
"""

import numpy as np
from pathlib import Path


DESTEKLENEN_ALGORITMALAR = [
    'mountainsort5',
    'spykingcircus2', 
    'tridesclous2',
    'simple',
]


def sort_nwb(dosya_yolu, algorithm='mountainsort5', 
             max_sure_sn=180.0, freq_min=300, freq_max=3000,
             verbose=True):
    """
    NWB dosyasindaki ham sinyali SpikeInterface ile spike sort eder.
    
    Parametreler:
        dosya_yolu  : str — NWB dosya yolu
        algorithm   : str — sorting algoritmasi
        max_sure_sn : float — max kayit suresi (RAM koruması)
        freq_min    : int — bandpass alt sinir (Hz)
        freq_max    : int — bandpass ust sinir (Hz)
        verbose     : bool — ilerleme mesajlari
    
    Returns:
        dict: {unit_id (int): spike_times (np.array, saniye)}
        Bos dict: sorting basarisiz olursa
    """
    try:
        import spikeinterface.full as si
        from spikeinterface.preprocessing import bandpass_filter, zscore
    except ImportError:
        print("  HATA: spikeinterface yuklu degil.")
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
    except Exception as e:
        # electrical_series_name bulamazsa dinamik bul
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
                print(f"  HATA: ElectricalSeries bulunamadi.")
                return {}
            recording = si.read_nwb(dosya_yolu, electrical_series_path=es_name)
        except Exception as e2:
            print(f"  HATA: NWB okunamadi: {e2}")
            return {}

    if verbose:
        sr = recording.get_sampling_frequency()
        n_ch = recording.get_num_channels()
        n_s = recording.get_num_samples()
        sure = n_s / sr
        print(f"  Sampling rate : {sr:.0f} Hz")
        print(f"  Kanal sayisi  : {n_ch}")
        print(f"  Kayit suresi  : {sure:.1f} s")

    # 2. Kayit suresini sinirla
    sr = recording.get_sampling_frequency()
    n_samples_max = int(max_sure_sn * sr)
    if recording.get_num_samples() > n_samples_max:
        recording = recording.frame_slice(0, n_samples_max)
        if verbose:
            print(f"  Kayit kisaltildi: {max_sure_sn:.0f} s")

    # 3. Preprocessing
    if verbose:
        print(f"  Bandpass: {freq_min}-{freq_max} Hz...")
    # Unsigned type kontrolu ve donusumu
    try:
        from spikeinterface.preprocessing import unsigned_to_signed
        recording = unsigned_to_signed(recording)
    except Exception:
        pass
    recording_f = bandpass_filter(recording, freq_min=freq_min, freq_max=freq_max)
    
    # 4. Sorting
    if verbose:
        print(f"  Sorting basliyor ({algorithm})...")
    
    try:
        if algorithm == 'mountainsort5':
            import shutil, os
            output_dir = os.path.join(os.path.dirname(__file__), 'mountainsort5_output')
            if os.path.exists(output_dir):
                shutil.rmtree(output_dir)
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
            print(f"  Desteklenenler: {DESTEKLENEN_ALGORITMALAR}")
            return {}
    except Exception as e:
        print(f"  HATA: Sorting basarisiz: {e}")
        return {}

    # 5. Spike times al
    unit_ids = sorting.get_unit_ids()
    if verbose:
        print(f"  Bulunan unit sayisi: {len(unit_ids)}")

    spike_dict = {}
    for uid in unit_ids:
        frames = sorting.get_unit_spike_train(uid, segment_index=0)
        times_s = frames / sr
        spike_dict[int(uid)] = times_s

    return spike_dict


def interaktif_sort(dosya_yolu, verbose=True):
    """
    Interaktif algoritma secimi ile sorting.
    organoid_units_analiz.py icin.
    """
    print()
    print("  SpikeInterface Spike Sorting")
    print("  " + "-"*40)
    for i, alg in enumerate(DESTEKLENEN_ALGORITMALAR, 1):
        print(f"  {i} : {alg}")
    print("  0 : Iptal")
    print()
    
    secim = input("  Algoritma secin (0-4): ").strip()
    
    if secim == '0' or not secim:
        return {}
    
    try:
        idx = int(secim) - 1
        if 0 <= idx < len(DESTEKLENEN_ALGORITMALAR):
            algorithm = DESTEKLENEN_ALGORITMALAR[idx]
        else:
            print("  Gecersiz secim.")
            return {}
    except ValueError:
        print("  Gecersiz secim.")
        return {}
    
    return sort_nwb(dosya_yolu, algorithm=algorithm, verbose=verbose)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        result = sort_nwb(sys.argv[1])
        print(f"\nSonuc: {len(result)} unit")
        for uid, times in list(result.items())[:5]:
            print(f"  Unit {uid}: {len(times)} spike, {times[0]:.3f}-{times[-1]:.3f} s")
    else:
        print("Kullanim: python si_sorting.py <dosya.nwb>")
