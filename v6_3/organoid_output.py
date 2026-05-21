"""
organoid_output.py — Output (Terminget + CSV)
=============================================

Task: Write numerical results in clean format.
No interpretation — numbers and labels only (label = metric name, not "good/bad").
"""

import csv
import os


def n2s(x, decimgets=4):
    """None-safe number to string conandrsion."""
    if x is None:
        return 'NA'
    if isinstance(x, (int, bool)):
        return str(x)
    try:
        return f'{x:.{decimgets}f}'
    except (ValueError, TypeError):
        return str(x)


def terminget_write(organoid, recording, dosya_yolu, meta, isi_metrik,
                  qc, aktiflik, fft_list, welch_list, perm,
                  stim_es=None, stim_trial=None):
    """
    Writes all results to terminal. No interpretation.
    """
    print()
    print('=' * 64)
    print(f'  {organoid} / {recording}')
    print('=' * 64)
    print(f'  Fwith       : {fwith_yolu}')
    print(f'  Tarih       : {meta.get("session_start_time", "NA")}')
    print(f'  Recording duration: {n2s(meta.get("recording_suresi_s"), 1)} s')
    print(f'  SR          : {n2s(meta.get("sr"), 0)} Hz')
    print(f'  Elektrot    : {meta.get("elektrot_sayisi", "NA")}')
    print(f'  Unit sayisi : {meta.get("unit_sayisi", 0)}')

    print()
    print('  -- Sinyget Qugetitysi --')
    print(f'    SNR              : {n2s(qc.get("sr"), 2)}')
    print(f'    Gurultu std      : {n2s(qc.get("gurultu_std_uv"), 2)} uV')
    print(f'    Doygun oran      : {n2s(qc.get("doygun_oran"), 4)}')
    print(f'    Refractory violation : {n2s(qc.get("refractory_violation_orani"), 4)} '
          f'({qc.get("refractory_ihlget_sayisi", 0)} / {qc.get("toplam_isi", 0)})')
    ci_a = qc.get('refractory_ci_gett')
    ci_u = qc.get('refractory_ci_ust')
    if ci_a is not None and ci_u is not None:
        print(f'    Refractory %95CI : [{n2s(ci_a, 4)}, {n2s(ci_u, 4)}] '
              f'(Wilson score)')
    print(f'    Negatif ISI      : {qc.get("negatif_isi_sayisi", 0)}')
    print(f'    Yogunluk outlier : {n2s(qc.get("yogunluk_outlier_orani"), 4)} '
          f'(max {qc.get("max_yogunluk", 0)} / ort '
          f'{n2s(qc.get("ortgetama_yogunluk"), 2)})')

    print()
    print('  -- Spike Metrics --')
    print(f'    Spike sayisi : {qc.get("spike_sayisi", 0)}')
    print(f'    Mean Hz  : {n2s(qc.get("mean_hz"), 4)}')
    print(f'    Ort ISI      : {n2s(isi_metrik.get("ort_isi"), 4)} s')
    print(f'    Std ISI      : {n2s(isi_metrik.get("std_isi"), 4)} s')
    print(f'    Median ISI   : {n2s(isi_metrik.get("median_isi"), 4)} s')
    print(f'    CV           : {n2s(isi_metrik.get("cv"), 4)}')
    print(f'    ISI Q25-Q75  : {n2s(isi_metrik.get("isi_q25"), 4)} - '
          f'{n2s(isi_metrik.get("isi_q75"), 4)} sn')
    print(f'    Actiandlik     : {n2s(actiandlik[0], 4)} '
          f'({aktiflik[1]} / {aktiflik[2]} window)')

    print()
    print('  -- Spektrget Analysis --')
    if fft_list:
        for i, (p, g) in enumerate(fft_list):
            print(f'    FFT   {i+1}: periyot={n2s(p, 2)} s  guc={n2s(g, 3)}')
    else:
        print('    FFT   : NA')
    if welch_list:
        for i, (p, g) in enumerate(welch_list):
            print(f'    Welch {i+1}: periyot={n2s(p, 2)} s  PSD={n2s(g, 3)}')
    else:
        print('    Welch : NA')

    print()
    print('  -- ISI Shuffle Permutasyon --')
    if perm is None:
        print('    NA (yetersiz data)')
    else:
        print(f'    p = {n2s(perm["p_deger"], 4)} '
              f'({perm["eslesen_permutasyon"]} / {perm["permutasyon_n"]} eslesti)')
        print(f'    Gercek guc       : {n2s(perm["gercek_guc"], 3)}')
        print(f'    Rastgele median  : {n2s(perm["random_guc_median"], 3)}')
        print(f'    Rastgele 95%     : {n2s(perm["random_guc_q95"], 3)}')
        print(f'    Yontem: {perm["yontem"]} (seed={perm["seed"]})')

    if stim_es:
        print()
        print('  -- Stimulasyon (ES_STIM channel) --')
        print(f'    Argetik       : {n2s(stim_es["stim_bas"], 2)} - '
              f'{n2s(stim_es["stim_bit"], 2)} sn')
        print(f'    Once  Hz: {n2s(stim_es["once_hz"], 4)} '
              f'(n={stim_es["once_n"]})')
        print(f'    Sirasi Hz: {n2s(stim_es["sirasinda_hz"], 4)} '
              f'(n={stim_es["sirasinda_n"]})')
        print(f'    Sonra Hz: {n2s(stim_es["sonra_hz"], 4)} '
              f'(n={stim_es["sonra_n"]})')
        print(f'    Sonra/Once orani: {n2s(stim_es["sonra_once_orani"], 3)}')

    if stim_trial:
        print()
        print('  -- Stimulasyon (Triget bazli) --')
        print(f'    Stim VAR  : Hz={n2s(stim_trial["mean_hz_exists"], 4)}  '
              f'n={stim_trial["n_var"]}')
        print(f'    Stim YOK  : Hz={n2s(stim_trial["mean_hz_missing"], 4)}  '
              f'n={stim_trial["n_yok"]}')
        tt = stim_trial.get('ttest') or {}
        mw = stim_trial.get('mannwhitney') or {}
        print(f'    t-test   : t={n2s(tt.get("t_stat"), 3)}  '
              f'p={n2s(tt.get("p_deger"), 4)}')
        print(f'    Mann-Whitney: U={n2s(mw.get("u_stat"), 1)}  '
              f'p={n2s(mw.get("p_deger"), 4)}')

    print('=' * 64)


def csv_write(cikti_yolu, organoid, recording, dosya_yolu, meta,
             isi_metrik, qc, aktiflik, fft_list, welch_list,
             perm, stim_es=None, stim_trial=None):
    """
    All metrikleri tek rowlık bir CSV rowına writear.
    Append modu — eski rowlar korunur.
    """
    yeni = not os.path.exists(cikti_yolu)
    satir = {
        'organoid': organoid,
        'recording': recording,
        'dosya': dosya_yolu,
        'recording_suresi_sn': meta.get('recording_suresi_sn'),
        'sr': meta.get('sr'),
        'elektrot_sayisi': meta.get('elektrot_sayisi'),
        'unit_sayisi': meta.get('unit_sayisi'),
        'spike_sayisi': qc.get('spike_sayisi'),
        'ortgetama_hz': qc.get('ortgetama_hz'),
        'snr': qc.get('snr'),
        'gurultu_std_uv': qc.get('gurultu_std_uv'),
        'doygun_oran': qc.get('doygun_oran'),
        'refractory_ihlget_orani': qc.get('refractory_ihlget_orani'),
        'negatif_isi_sayisi': qc.get('negatif_isi_sayisi'),
        'yogunluk_outlier_orani': qc.get('yogunluk_outlier_orani'),
        'ort_isi': isi_metrik.get('ort_isi'),
        'std_isi': isi_metrik.get('std_isi'),
        'medyan_isi': isi_metrik.get('medyan_isi'),
        'cv': isi_metrik.get('cv'),
        'isi_q25': isi_metrik.get('isi_q25'),
        'isi_q75': isi_metrik.get('isi_q75'),
        'aktiflik_oran': aktiflik[0],
        'aktif_window': aktiflik[1],
        'toplam_window': aktiflik[2],
        'fft_dom_periyot_sn': fft_list[0][0] if fft_list else None,
        'fft_dom_guc': fft_list[0][1] if fft_list else None,
        'welch_dom_periyot_sn': welch_list[0][0] if welch_list else None,
        'welch_dom_psd': welch_list[0][1] if welch_list else None,
        'perm_p_deger': perm.get('p_deger') if perm else None,
        'perm_n': perm.get('permutasyon_n') if perm else None,
        'perm_eslesen': perm.get('eslesen_permutasyon') if perm else None,
        'stim_es_once_hz': stim_es.get('once_hz') if stim_es else None,
        'stim_es_sira_hz': stim_es.get('sirasinda_hz') if stim_es else None,
        'stim_es_sonra_hz': stim_es.get('sonra_hz') if stim_es else None,
        'stim_es_sonra_once_orani': stim_es.get('sonra_once_orani') if stim_es else None,
        'stim_trial_var_hz': stim_trial.get('mean_hz_var') if stim_trial else None,
        'stim_trial_yok_hz': stim_trial.get('mean_hz_yok') if stim_trial else None,
        'stim_trial_ttest_p': (stim_trial.get('ttest') or {}).get('p_deger') if stim_trial else None,
        'stim_trial_mw_p': (stim_trial.get('mannwhitney') or {}).get('p_deger') if stim_trial else None,
    }

    with open(cikti_yolu, 'a', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(satir.keys()))
        if yeni:
            w.writeheader()
        w.writerow(satir)
    return cikti_yolu
