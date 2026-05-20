"""
organoid_output.py — Çıktı (Terminal + CSV)
=============================================

Görev: Sayısal sonuçları temiz formatta yazdır.
Yorum yok — sadece sayı ve etiket (etiket = metrik adı, "iyi/kötü" değil).
"""

import csv
import os


def n2s(x, decimals=4):
    """None-safe sayı stringe çevir."""
    if x is None:
        return 'NA'
    if isinstance(x, (int, bool)):
        return str(x)
    try:
        return f'{x:.{decimals}f}'
    except (ValueError, TypeError):
        return str(x)


def terminal_yaz(organoid, kayit, dosya_yolu, meta, isi_metrik,
                  qc, aktiflik, fft_list, welch_list, perm,
                  stim_es=None, stim_trial=None):
    """
    Tüm sonuçları terminale yazar. Yorum yok.
    """
    print()
    print('=' * 64)
    print(f'  {organoid} / {kayit}')
    print('=' * 64)
    print(f'  Dosya       : {dosya_yolu}')
    print(f'  Tarih       : {meta.get("session_start_time", "NA")}')
    print(f'  Kayit suresi: {n2s(meta.get("kayit_suresi_sn"), 1)} sn')
    print(f'  SR          : {n2s(meta.get("sr"), 0)} Hz')
    print(f'  Elektrot    : {meta.get("elektrot_sayisi", "NA")}')
    print(f'  Unit sayisi : {meta.get("unit_sayisi", 0)}')

    print()
    print('  -- Sinyal Kalitesi --')
    print(f'    SNR              : {n2s(qc.get("snr"), 2)}')
    print(f'    Gurultu std      : {n2s(qc.get("gurultu_std_uv"), 2)} uV')
    print(f'    Doygun oran      : {n2s(qc.get("doygun_oran"), 4)}')
    print(f'    Refractory ihlal : {n2s(qc.get("refractory_ihlal_orani"), 4)} '
          f'({qc.get("refractory_ihlal_sayisi", 0)} / {qc.get("toplam_isi", 0)})')
    ci_a = qc.get('refractory_ci_alt')
    ci_u = qc.get('refractory_ci_ust')
    if ci_a is not None and ci_u is not None:
        print(f'    Refractory %95CI : [{n2s(ci_a, 4)}, {n2s(ci_u, 4)}] '
              f'(Wilson score)')
    print(f'    Negatif ISI      : {qc.get("negatif_isi_sayisi", 0)}')
    print(f'    Yogunluk outlier : {n2s(qc.get("yogunluk_outlier_orani"), 4)} '
          f'(max {qc.get("max_yogunluk", 0)} / ort '
          f'{n2s(qc.get("ortalama_yogunluk"), 2)})')

    print()
    print('  -- Spike Metrikleri --')
    print(f'    Spike sayisi : {qc.get("spike_sayisi", 0)}')
    print(f'    Ortalama Hz  : {n2s(qc.get("ortalama_hz"), 4)}')
    print(f'    Ort ISI      : {n2s(isi_metrik.get("ort_isi"), 4)} sn')
    print(f'    Std ISI      : {n2s(isi_metrik.get("std_isi"), 4)} sn')
    print(f'    Medyan ISI   : {n2s(isi_metrik.get("medyan_isi"), 4)} sn')
    print(f'    CV           : {n2s(isi_metrik.get("cv"), 4)}')
    print(f'    ISI Q25-Q75  : {n2s(isi_metrik.get("isi_q25"), 4)} - '
          f'{n2s(isi_metrik.get("isi_q75"), 4)} sn')
    print(f'    Aktiflik     : {n2s(aktiflik[0], 4)} '
          f'({aktiflik[1]} / {aktiflik[2]} pencere)')

    print()
    print('  -- Spektral Analiz --')
    if fft_list:
        for i, (p, g) in enumerate(fft_list):
            print(f'    FFT   {i+1}: periyot={n2s(p, 2)} sn  guc={n2s(g, 3)}')
    else:
        print('    FFT   : NA')
    if welch_list:
        for i, (p, g) in enumerate(welch_list):
            print(f'    Welch {i+1}: periyot={n2s(p, 2)} sn  PSD={n2s(g, 3)}')
    else:
        print('    Welch : NA')

    print()
    print('  -- ISI Shuffle Permutasyon --')
    if perm is None:
        print('    NA (yetersiz veri)')
    else:
        print(f'    p = {n2s(perm["p_deger"], 4)} '
              f'({perm["eslesen_permutasyon"]} / {perm["permutasyon_n"]} eslesti)')
        print(f'    Gercek guc       : {n2s(perm["gercek_guc"], 3)}')
        print(f'    Rastgele medyan  : {n2s(perm["rastgele_guc_medyan"], 3)}')
        print(f'    Rastgele 95%     : {n2s(perm["rastgele_guc_q95"], 3)}')
        print(f'    Yontem: {perm["yontem"]} (seed={perm["seed"]})')

    if stim_es:
        print()
        print('  -- Stimulasyon (ES_STIM kanal) --')
        print(f'    Aralik       : {n2s(stim_es["stim_bas"], 2)} - '
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
        print('  -- Stimulasyon (Trial bazli) --')
        print(f'    Stim VAR  : Hz={n2s(stim_trial["mean_hz_var"], 4)}  '
              f'n={stim_trial["n_var"]}')
        print(f'    Stim YOK  : Hz={n2s(stim_trial["mean_hz_yok"], 4)}  '
              f'n={stim_trial["n_yok"]}')
        tt = stim_trial.get('ttest') or {}
        mw = stim_trial.get('mannwhitney') or {}
        print(f'    t-test   : t={n2s(tt.get("t_stat"), 3)}  '
              f'p={n2s(tt.get("p_deger"), 4)}')
        print(f'    Mann-Whitney: U={n2s(mw.get("u_stat"), 1)}  '
              f'p={n2s(mw.get("p_deger"), 4)}')

    print('=' * 64)


def csv_yaz(cikti_yolu, organoid, kayit, dosya_yolu, meta,
             isi_metrik, qc, aktiflik, fft_list, welch_list,
             perm, stim_es=None, stim_trial=None):
    """
    Tüm metrikleri tek satırlık bir CSV satırına yazar.
    Append modu — eski satırlar korunur.
    """
    yeni = not os.path.exists(cikti_yolu)
    satir = {
        'organoid': organoid,
        'kayit': kayit,
        'dosya': dosya_yolu,
        'kayit_suresi_sn': meta.get('kayit_suresi_sn'),
        'sr': meta.get('sr'),
        'elektrot_sayisi': meta.get('elektrot_sayisi'),
        'unit_sayisi': meta.get('unit_sayisi'),
        'spike_sayisi': qc.get('spike_sayisi'),
        'ortalama_hz': qc.get('ortalama_hz'),
        'snr': qc.get('snr'),
        'gurultu_std_uv': qc.get('gurultu_std_uv'),
        'doygun_oran': qc.get('doygun_oran'),
        'refractory_ihlal_orani': qc.get('refractory_ihlal_orani'),
        'negatif_isi_sayisi': qc.get('negatif_isi_sayisi'),
        'yogunluk_outlier_orani': qc.get('yogunluk_outlier_orani'),
        'ort_isi': isi_metrik.get('ort_isi'),
        'std_isi': isi_metrik.get('std_isi'),
        'medyan_isi': isi_metrik.get('medyan_isi'),
        'cv': isi_metrik.get('cv'),
        'isi_q25': isi_metrik.get('isi_q25'),
        'isi_q75': isi_metrik.get('isi_q75'),
        'aktiflik_oran': aktiflik[0],
        'aktif_pencere': aktiflik[1],
        'toplam_pencere': aktiflik[2],
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
