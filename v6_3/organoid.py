"""
organoid.py — Ana Komut (v6)
==============================

Kullanim:
    python organoid.py <dosya.nwb> <organoid> <kayit> [--secenekler]
    python organoid.py --ozet [organoid]
    python organoid.py --listele

Secenekler:
    --kanal N        : Kayittan okunacak kanal (varsayilan: 0)
    --snr SIGMA      : Spike detection sigma esigi (varsayilan: 5)
    --bandpass L H   : Band-pass aralik Hz (varsayilan: 300 3000)
    --chunk SN       : Chunk boyutu sn (varsayilan: 60)
    --sort           : Tridesclous2 spike sorting calistir
    --csv DOSYA      : CSV ciktisi yaz (varsayilan: organoid_v6_results.csv)
    --grafik DIR     : Quick-look grafiklerini bu klasore kaydet
                        (varsayilan: graphs/)
    --no-grafik      : Grafik olusturma
    --metadata J     : Metadata JSON dosyasi (organoid yasi vs)
    --refractory MS  : Refractory esigi ms (varsayilan: 1.0)
    --pencere SN     : CV/aktiflik pencere boyutu (varsayilan: 30)

Felsefe:
    - Modul yorum yapmaz. Sadece sayi verir.
    - Kontrol grubu hardcoded yok. Karsilastirma icin organoid_compare.py.
    - SQLite olcum verisi icin, JSON metadata icin.
    - Quick-look grafik zorunlu (kapatilabilir: --no-grafik).
    - Spike sorting opsiyonel (--sort).
"""

import sys
import os
import argparse
import time
import json
import numpy as np

import organoid_io as oio
import organoid_signal as osig
import organoid_qc as oqc
import organoid_metrics as omet
import organoid_db as odb
import organoid_plot as oplt
import organoid_output as oout


def analiz_calistir(args):
    """Tek bir NWB dosyasını analiz eder."""
    t0 = time.time()
    dosya = args.dosya
    organoid = args.organoid
    kayit = args.kayit

    if not os.path.exists(dosya):
        print(f'HATA: Dosya yok: {dosya}')
        return 1

    print(f'\n>> {organoid}/{kayit}')
    print(f'   Dosya: {dosya}')

    # 1. Tip ve metadata
    tip = oio.nwb_tip_tespit(dosya)
    meta = oio.nwb_metadata_oku(dosya)
    print(f'   Tip: {tip}, sure: {meta.get("kayit_suresi_sn")} sn, '
          f'elektrot: {meta.get("elektrot_sayisi")}, '
          f'unit: {meta.get("unit_sayisi")}')

    # DB bağlantısı
    conn = odb.db_ac()
    kayit_id = odb.kayit_ekle_veya_guncelle(
        conn, organoid, kayit, dosya, meta,
        kanal_indeks=args.kanal)

    # Metadata JSON
    if args.metadata and os.path.exists(args.metadata):
        with open(args.metadata, 'r', encoding='utf-8') as f:
            md = json.load(f)
        odb.metadata_yaz(organoid, kayit, md)
        print(f'   Metadata yazildi')

    # 2. Spike tespiti
    spike_zaman = None
    gurultu_std = None
    snr = None
    snr_genlik_uv = None
    doygun_oran = 0.0

    if args.sort and tip in ('ham', 'her_ikisi'):
        # Spike sorting yolu
        try:
            import organoid_sorting as osort
        except ImportError:
            if args.sort_or_mua:
                print('   UYARI: SpikeInterface yok, MUA fallback aktif.')
                print('   NOT: Sonuçlar "MUA-fallback" olarak işaretlenecek.')
                args.sort = False  # MUA yoluna düş
            else:
                print('   HATA: SpikeInterface kurulu değil.')
                print('   Çözüm: pip install spikeinterface')
                print('   Veya --sort-or-mua kullan (MUA fallback).')
                conn.close()
                return 2

        if args.sort:
            print('   Spike sorting calistiriliyor (Tridesclous2)...')
            sort_klasor = f'sorting_{organoid}_{kayit}'
            sort_sonuc = osort.spike_sort_calistir(
                dosya, sort_klasor,
                sigma_esik=args.snr,
                freq_min=args.bandpass[0],
                freq_max=args.bandpass[1])

            if sort_sonuc is None:
                if args.sort_or_mua:
                    print('   UYARI: Sorting başarısız, MUA fallback aktif.')
                    print('   NOT: Sonuçlar "MUA-fallback" olarak işaretlenecek.')
                    args.sort = False  # MUA yoluna düş
                else:
                    print('   HATA: Spike sorting başarısız.')
                    print('   Analiz durduruldu. Fallback yok.')
                    print('   --sort-or-mua kullan (etiketli MUA fallback).')
                    conn.close()
                    return 3
            else:
                print(f'   {sort_sonuc["n_unit"]} unit bulundu')
                spike_zaman = osort.units_birlestir(sort_sonuc['units'])
                if len(spike_zaman) == 0:
                    if args.sort_or_mua:
                        print('   UYARI: Sorting tamamlandı ama unit yok, MUA fallback.')
                        args.sort = False
                    else:
                        print('   HATA: Sorting tamamlandı ama hiç unit/spike yok.')
                        conn.close()
                        return 3
                else:
                    for u in sort_sonuc['units']:
                        odb.spike_zamanlari_yaz(
                            conn, kayit_id, u['spike_zaman'],
                            unit_id=u['unit_id'])

    # MUA yolu (sort yok veya fallback)
    if (not args.sort) and tip in ('ham', 'her_ikisi') and spike_zaman is None:
        # Hangi modda olduğumuzu belirt
        mod_etiket = 'MUA-fallback' if hasattr(args, '_fallback_oldu') else 'MUA'
        gurultu_yontem = 'MAD' if args.mad_esik else 'std'
        print(f'   Mod: {mod_etiket} | Gurultu: {gurultu_yontem}')
        print(f'   Gurultu std hesaplaniyor (n=10 ornek)...')
        gurultu_std = osig.gurultu_std_ornekle(
            oio, dosya, meta['sr'],
            args.bandpass[0], args.bandpass[1],
            n_ornek=10, ornek_sn=30, kanal=args.kanal,
            mad_kullan=args.mad_esik)
        gurultu_std_uv = gurultu_std * 1e6
        print(f'   Gurultu std: {gurultu_std_uv:.2f} uV ({gurultu_yontem})')

        print(f'   Spike detection (sigma={args.snr}, sadece negatif)...')
        spike_t_listesi = []
        ilk_chunk_f = None
        ilk_spike_idx = None
        toplam_doygun = 0.0
        chunk_count = 0

        for chunk_data, bas, bit, sr in oio.ham_chunk_uret(
                dosya, chunk_sn=args.chunk, kanal=args.kanal):
            doygun = osig.doygunluk_tespit(chunk_data)
            toplam_doygun += doygun
            chunk_count += 1

            spike_idx, chunk_f = osig.spike_tespit_chunk(
                chunk_data, sr,
                low_hz=args.bandpass[0],
                high_hz=args.bandpass[1],
                sigma_esik=args.snr,
                gurultu_std=gurultu_std,
                refractory_ms=args.refractory,
                sadece_negatif=args.sadece_negatif)

            spike_t = (spike_idx + (bas - int(5 * sr))) / sr
            spike_t_listesi.extend(spike_t.tolist())

            if ilk_chunk_f is None and len(spike_idx) > 0:
                ilk_chunk_f = chunk_f
                ilk_spike_idx = spike_idx

        spike_zaman = np.array(spike_t_listesi, dtype=np.float64)
        doygun_oran = toplam_doygun / chunk_count if chunk_count > 0 else 0.0

        # SNR hesapla
        if ilk_chunk_f is not None and ilk_spike_idx is not None:
            snr, snr_genlik_uv, _ = osig.snr_hesapla(
                ilk_chunk_f, ilk_spike_idx, meta['sr'], gurultu_std)

        odb.spike_zamanlari_yaz(conn, kayit_id, spike_zaman, unit_id=-1)
        snr_str = f'{snr:.2f}' if snr else 'NA'
        print(f'   {len(spike_zaman)} spike (SNR: {snr_str})')

    elif tip == 'sorted':
        # NWB içinde zaten sorted veri var
        units = oio.sorted_units_oku(dosya)
        if units:
            tum = [u['spike_zaman'] for u in units]
            spike_zaman = np.sort(np.concatenate(tum)) if tum else \
                          np.array([])
            for u in units:
                odb.spike_zamanlari_yaz(
                    conn, kayit_id, u['spike_zaman'],
                    unit_id=u['unit_id'])
            print(f'   {len(units)} unit, toplam {len(spike_zaman)} spike')

    if spike_zaman is None or len(spike_zaman) == 0:
        print('   UYARI: Spike yok, analiz atlanıyor')
        conn.close()
        return 0

    sure_sn = meta.get('kayit_suresi_sn') or float(spike_zaman[-1])

    # 3. ISI metrikleri
    isi_metrik = omet.isi_metrikleri(spike_zaman)

    # 4. Aktiflik
    aktiflik = omet.aktiflik_pencere(
        spike_zaman, sure_sn, pencere_sn=args.pencere)

    # 5. CV zaman serisi
    cv_seri = omet.cv_zaman_serisi(
        spike_zaman, sure_sn, pencere_sn=args.pencere)
    cv_degerleri = [x[1] for x in cv_seri]

    # 6. Spektral
    fft_list = omet.fft_periyotlar(cv_degerleri, pencere_sn=args.pencere)
    welch_list = omet.welch_periyotlar(cv_degerleri, pencere_sn=args.pencere)

    # 7. Permütasyon
    perm = omet.isi_shuffle_permutasyon(
        spike_zaman, sure_sn, pencere_sn=args.pencere,
        n_permutasyon=1000, seed=42)

    # 8. Stim analizi
    stim_es = None
    stim_trial = None
    if meta.get('stim_var'):
        stim_bas, stim_bit = oio.stim_kanali_aralik(dosya)
        if stim_bas is not None:
            stim_es = omet.stim_oncesi_sirasi_sonrasi(
                spike_zaman, sure_sn, stim_bas, stim_bit)

    if meta.get('trial_var'):
        trials = oio.trial_oku(dosya)
        stim_trial = omet.trial_stim_karsilastir(spike_zaman, trials)

    # 9. QC metrikleri
    qc = oqc.kalite_metrikleri_topla(
        spike_zaman, sure_sn,
        snr=snr,
        gurultu_std_uv=(gurultu_std * 1e6) if gurultu_std else None,
        doygun_oran=doygun_oran,
        refractory_ms=args.refractory)
    qc['gurultu_yontem'] = 'MAD' if args.mad_esik else 'std'

    # 10. DB'ye yaz
    odb.metrikler_yaz(conn, kayit_id, isi_metrik,
                       aktiflik[0], aktiflik[1], aktiflik[2],
                       qc.get('spike_sayisi'), qc.get('ortalama_hz'))
    odb.qc_yaz(conn, kayit_id, qc)
    odb.spektral_yaz(conn, kayit_id, fft_list, welch_list)
    odb.permutasyon_yaz(conn, kayit_id, perm)
    if stim_es:
        odb.stim_yaz(conn, kayit_id, 'es_stim', stim_es)
    if stim_trial:
        # trial sonucunu stim_es formatına uydur
        stim_trial_for_db = {
            'stim_bas': None, 'stim_bit': None,
            'once_n': stim_trial['n_yok'],
            'sirasinda_n': stim_trial['n_var'],
            'sonra_n': None,
            'once_hz': stim_trial['mean_hz_yok'],
            'sirasinda_hz': stim_trial['mean_hz_var'],
            'sonra_hz': None,
            'sonra_once_orani': None,
            'ttest': stim_trial.get('ttest'),
            'mannwhitney': stim_trial.get('mannwhitney'),
            'n_var': stim_trial['n_var'],
            'n_yok': stim_trial['n_yok'],
        }
        odb.stim_yaz(conn, kayit_id, 'trial', stim_trial_for_db)

    # 11. CSV
    csv_yolu = args.csv or 'organoid_v6_results.csv'
    oout.csv_yaz(csv_yolu, organoid, kayit, dosya, meta,
                  isi_metrik, qc, aktiflik, fft_list, welch_list,
                  perm, stim_es, stim_trial)
    print(f'   CSV: {csv_yolu}')

    # 12. Grafik
    if not args.no_grafik:
        grafik_dir = args.grafik or 'graphs'
        os.makedirs(grafik_dir, exist_ok=True)
        backend = 'plotly' if args.plotly else 'matplotlib'
        ext = '.html' if args.plotly else '.png'
        ql_yolu = os.path.join(grafik_dir,
                                f'{organoid}_{kayit}_quicklook{ext}')
        oplt.quick_look(spike_zaman, sure_sn, ql_yolu,
                         baslik_ek=f'{organoid} / {kayit}',
                         backend=backend)
        if not args.plotly:
            cv_yolu = os.path.join(grafik_dir, f'{organoid}_{kayit}_cv.png')
            oplt.cv_plot(cv_seri, cv_yolu,
                          baslik_ek=f'{organoid} / {kayit} - CV ({args.pencere}s)')
        print(f'   Grafik: {ql_yolu}')

    # 13. Terminal
    oout.terminal_yaz(organoid, kayit, dosya, meta,
                       isi_metrik, qc, aktiflik, fft_list, welch_list,
                       perm, stim_es, stim_trial)

    print(f'\n   Sure: {time.time()-t0:.1f} sn')
    conn.close()
    return 0


def ozet_goster(args):
    """Veritabanındaki kayıtların özetini gösterir."""
    conn = odb.db_ac()
    organoid = args.organoid if hasattr(args, 'organoid') and args.organoid else None
    rows = odb.kayit_ozet(conn, organoid)
    conn.close()
    if not rows:
        print('Kayit yok.')
        return 0

    print(f'\n{"Organoid":<12} {"Kayit":<14} {"Spike":>7} {"Hz":>9} '
          f'{"CV":>7} {"SNR":>6} {"RefIhl":>7}')
    print('-' * 68)
    for r in rows:
        print(f'{r["organoid"]:<12} {r["kayit"]:<14} '
              f'{r["spike_sayisi"] or "":>7} '
              f'{oout.n2s(r["ortalama_hz"], 3):>9} '
              f'{oout.n2s(r["cv"], 3):>7} '
              f'{oout.n2s(r["snr"], 2):>6} '
              f'{oout.n2s(r["refractory_ihlal_orani"], 3):>7}')
    return 0


def main():
    # Geriye dönük kullanım: ilk argüman .nwb dosyası ise 'analiz' olarak yorumla
    argv = sys.argv[1:]
    if len(argv) >= 1:
        if argv[0].endswith('.nwb') and os.path.exists(argv[0]):
            argv = ['analiz'] + argv
        elif argv[0] == '--ozet':
            argv = ['ozet'] + argv[1:]
        elif argv[0] == '--listele':
            argv = ['listele']

    parser = argparse.ArgumentParser(
        description='Organoid analiz - v6',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest='cmd')

    # Analiz
    p_an = sub.add_parser('analiz', help='NWB dosyasini analiz et')
    p_an.add_argument('dosya', help='NWB dosya yolu')
    p_an.add_argument('organoid', help='Organoid kodu (orn: SO1)')
    p_an.add_argument('kayit', help='Kayit kodu (orn: 13_Haz)')
    p_an.add_argument('--kanal', type=int, default=0)
    p_an.add_argument('--snr', type=float, default=5.0,
                       help='Spike detection sigma esigi')
    p_an.add_argument('--bandpass', type=float, nargs=2,
                       default=[300.0, 3000.0],
                       metavar=('LOW', 'HIGH'))
    p_an.add_argument('--chunk', type=int, default=60,
                       help='Chunk boyutu sn')
    p_an.add_argument('--sort', action='store_true',
                       help='Tridesclous2 ile spike sorting (başarısızsa exit 3)')
    p_an.add_argument('--sort-or-mua', action='store_true', dest='sort_or_mua',
                       help='Sort başarısızsa MUA fallback (etiketli)')
    p_an.add_argument('--mad-esik', action='store_true', dest='mad_esik',
                       help='Gurultu için MAD/0.6745 (Quiroga 2004), '
                            'yüksek aktiviteli kayıtlarda önerilir')
    p_an.add_argument('--sadece-negatif', action='store_true',
                       dest='sadece_negatif',
                       help='Sadece negatif eşikle spike tespit (default: çift yönlü)')
    p_an.add_argument('--csv', default=None)
    p_an.add_argument('--grafik', default=None,
                       help='Grafik klasoru (varsayilan: graphs/)')
    p_an.add_argument('--no-grafik', action='store_true',
                       dest='no_grafik')
    p_an.add_argument('--plotly', action='store_true',
                       help='Plotly ile interaktif HTML grafik (default: matplotlib PNG)')
    p_an.add_argument('--metadata', default=None,
                       help='Metadata JSON dosyasi')
    p_an.add_argument('--refractory', type=float, default=1.0,
                       help='Refractory esigi ms')
    p_an.add_argument('--pencere', type=int, default=30,
                       help='CV/aktiflik pencere sn')

    # Ozet
    p_oz = sub.add_parser('ozet', help='Veritabani ozeti')
    p_oz.add_argument('organoid', nargs='?', default=None)

    # Listele
    sub.add_parser('listele', help='Kayitlari listele')

    args = parser.parse_args(argv)

    if args.cmd is None:
        parser.print_help()
        return 0

    if args.cmd == 'analiz':
        return analiz_calistir(args)
    elif args.cmd in ('ozet', 'listele'):
        return ozet_goster(args)
    else:
        parser.print_help()
        return 0


if __name__ == '__main__':
    sys.exit(main())
