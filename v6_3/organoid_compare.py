"""
organoid_compare.py — Karşılaştırma Komutu
============================================

Görev: İki veya daha fazla kayıt arasında istatistiksel karşılaştırma.
       Tüm istatistik testler organoid_metrics.py'den çağırılır.
       Bu modül kendi başına hesap yapmaz — sadece organize eder.

Kullanim:
    python organoid_compare.py <organoid1> <kayit1> <organoid2> <kayit2> [...]

Ornek:
    python organoid_compare.py SO1 13_Haz SO1 14_Haz
    python organoid_compare.py SO1 13_Haz SO1 14_Haz BO2 11_Tem

Cikti: terminal + CSV (compare_results.csv)
       Yorum yok. Sadece test sonuclari ve sayisal degerler.
"""

import sys
import csv
import argparse
import numpy as np

import organoid_db as odb
import organoid_metrics as omet
import organoid_output as oout


def kayit_metriklerini_al(conn, organoid, kayit):
    """
    Bir kaydın spike zamanlarını ve metriklerini DB'den çeker.
    """
    row = conn.execute("""
        SELECT k.id, k.kayit_suresi_sn,
               m.spike_sayisi, m.ortalama_hz, m.cv, m.ort_isi, m.medyan_isi,
               q.snr, q.refractory_ihlal_orani
        FROM kayitlar k
        LEFT JOIN metrikler m ON m.kayit_id = k.id
        LEFT JOIN qc_metrikleri q ON q.kayit_id = k.id
        WHERE k.organoid=? AND k.kayit=? AND k.kanal_indeks=0
    """, (organoid, kayit)).fetchone()

    if row is None:
        return None

    spike_zaman = odb.spike_zaman_oku(conn, organoid, kayit)

    return {
        'organoid': organoid,
        'kayit': kayit,
        'kayit_id': row['id'],
        'spike_zaman': spike_zaman,
        'sure_sn': row['kayit_suresi_sn'],
        'spike_sayisi': row['spike_sayisi'],
        'ortalama_hz': row['ortalama_hz'],
        'cv': row['cv'],
        'ort_isi': row['ort_isi'],
        'medyan_isi': row['medyan_isi'],
        'snr': row['snr'],
        'refractory_ihlal_orani': row['refractory_ihlal_orani'],
    }


def hz_pencere_dagilimi(spike_zaman, sure_sn, pencere_sn=10):
    """
    Pencere bazlı Hz dağılımı.
    İstatistiksel testler için bağımsız örnek seti üretir.
    """
    if not sure_sn or sure_sn < pencere_sn:
        return np.array([])
    rates = []
    t = 0.0
    while t + pencere_sn <= sure_sn:
        n = int(np.sum((spike_zaman >= t) & (spike_zaman < t + pencere_sn)))
        rates.append(n / pencere_sn)
        t += pencere_sn
    return np.array(rates)


def cv_pencere_dagilimi(spike_zaman, sure_sn, pencere_sn=30):
    """
    Pencere bazlı CV dağılımı.
    """
    seri = omet.cv_zaman_serisi(spike_zaman, sure_sn, pencere_sn)
    return np.array([x[1] for x in seri])


def autocorrelation_lag1(veri):
    """
    Lag-1 autocorrelation katsayısı.

    Yüksek değer (>0.3) ardışık örneklerin bağımsız olmadığını gösterir.
    Bu durumda standart t-test/MW p değerleri şişmiş olabilir.
    """
    veri = np.asarray(veri, dtype=np.float64)
    if len(veri) < 3:
        return None
    v1 = veri[:-1]
    v2 = veri[1:]
    if np.std(v1) == 0 or np.std(v2) == 0:
        return None
    return float(np.corrcoef(v1, v2)[0, 1])


def circular_shift_permutation(veri1, veri2, n_perm=1000, seed=42):
    """
    Circular shift permutation testi — autocorrelated zaman serisi
    için bağımsız örneklem varsayımına gerek olmayan p değeri.

    Yöntem: birleşik serinin rastgele dairesel kaydırmasıyla
    iki grubun ortalama farkı dağılımı oluşturulur. Gerçek fark
    bu dağılımın neresinde durur — ham p değeri.

    Dönen sözlük (None yetersiz veride):
        p_deger : iki yönlü
        n_perm : kullanılan permutasyon sayısı
        gercek_fark : ortalama1 - ortalama2
        rastgele_fark_medyan, rastgele_fark_std
    """
    v1 = np.asarray(veri1, dtype=np.float64)
    v2 = np.asarray(veri2, dtype=np.float64)
    if len(v1) < 3 or len(v2) < 3:
        return None

    gercek_fark = float(np.mean(v1) - np.mean(v2))
    birlesik = np.concatenate([v1, v2])
    n1 = len(v1)
    n_total = len(birlesik)

    rng = np.random.default_rng(seed)
    rastgele_farklar = []
    for _ in range(n_perm):
        kaydir = int(rng.integers(1, n_total))
        kaydirilmis = np.roll(birlesik, kaydir)
        g1 = kaydirilmis[:n1]
        g2 = kaydirilmis[n1:]
        rastgele_farklar.append(np.mean(g1) - np.mean(g2))

    rastgele_arr = np.array(rastgele_farklar)
    # İki yönlü p
    p = float(np.mean(np.abs(rastgele_arr) >= abs(gercek_fark)))

    return {
        'test': 'Circular shift permutation',
        'p_deger': p,
        'n_perm': n_perm,
        'gercek_fark': gercek_fark,
        'rastgele_fark_medyan': float(np.median(rastgele_arr)),
        'rastgele_fark_std': float(np.std(rastgele_arr)),
        'seed': seed,
    }


def ikili_karsilastir(k1, k2, pencere_hz=10, pencere_cv=30,
                      circular_perm=True):
    """
    İki kayıt arası ikili karşılaştırma.
    Tüm testler organoid_metrics.py'den çağırılır.

    Pencere bazlı veriler ardışıktır — autocorrelation bulunabilir.
    Bu nedenle:
      - Lag-1 autocorrelation hesaplanır ve raporlanır.
      - circular_perm=True ise circular shift permutation eklenir.

    Dönen sözlük:
        hz_ttest, hz_mw, hz_perm, hz_acf
        cv_ttest, cv_mw, cv_perm, cv_acf
        isi_ttest, isi_mw
        ozet : { özet sayısal farklar }
    """
    sonuc = {
        'kayit1': f'{k1["organoid"]}/{k1["kayit"]}',
        'kayit2': f'{k2["organoid"]}/{k2["kayit"]}',
    }

    # Pencere bazlı Hz
    hz1 = hz_pencere_dagilimi(k1['spike_zaman'], k1['sure_sn'], pencere_hz)
    hz2 = hz_pencere_dagilimi(k2['spike_zaman'], k2['sure_sn'], pencere_hz)

    sonuc['n_pencere_1'] = len(hz1)
    sonuc['n_pencere_2'] = len(hz2)

    sonuc['hz_ttest'] = omet.ttest_iki_grup(hz1, hz2)
    sonuc['hz_mw'] = omet.mannwhitney_iki_grup(hz1, hz2)
    sonuc['hz_acf_1'] = autocorrelation_lag1(hz1)
    sonuc['hz_acf_2'] = autocorrelation_lag1(hz2)
    if circular_perm:
        sonuc['hz_perm'] = circular_shift_permutation(hz1, hz2)
    else:
        sonuc['hz_perm'] = None

    # Pencere bazlı CV
    cv1 = cv_pencere_dagilimi(k1['spike_zaman'], k1['sure_sn'], pencere_cv)
    cv2 = cv_pencere_dagilimi(k2['spike_zaman'], k2['sure_sn'], pencere_cv)

    sonuc['cv_ttest'] = omet.ttest_iki_grup(cv1, cv2)
    sonuc['cv_mw'] = omet.mannwhitney_iki_grup(cv1, cv2)
    sonuc['cv_acf_1'] = autocorrelation_lag1(cv1)
    sonuc['cv_acf_2'] = autocorrelation_lag1(cv2)
    if circular_perm and len(cv1) >= 3 and len(cv2) >= 3:
        sonuc['cv_perm'] = circular_shift_permutation(cv1, cv2)
    else:
        sonuc['cv_perm'] = None

    # ISI (zaten bağımsız sayılır, ama büyük dağılım — sadece ttest+MW)
    isi1 = np.diff(np.sort(k1['spike_zaman']))
    isi2 = np.diff(np.sort(k2['spike_zaman']))

    sonuc['isi_ttest'] = omet.ttest_iki_grup(isi1, isi2)
    sonuc['isi_mw'] = omet.mannwhitney_iki_grup(isi1, isi2)

    # Özet sayılar
    sonuc['ozet'] = {
        'k1_hz': k1['ortalama_hz'],
        'k2_hz': k2['ortalama_hz'],
        'hz_fark': (k2['ortalama_hz'] or 0) - (k1['ortalama_hz'] or 0),
        'k1_cv': k1['cv'],
        'k2_cv': k2['cv'],
        'cv_fark': (k2['cv'] or 0) - (k1['cv'] or 0),
        'k1_n': k1['spike_sayisi'],
        'k2_n': k2['spike_sayisi'],
    }

    return sonuc


def coklu_karsilastir(kayit_listesi):
    """
    İkiden fazla kayıt için ANOVA + ikili karşılaştırmalar.
    """
    sonuc = {'n_kayit': len(kayit_listesi),
             'kayit_etiketleri': [f'{k["organoid"]}/{k["kayit"]}'
                                   for k in kayit_listesi]}

    # ANOVA: pencere bazlı Hz dağılımları
    hz_dagilimlari = []
    for k in kayit_listesi:
        hz = hz_pencere_dagilimi(k['spike_zaman'], k['sure_sn'])
        hz_dagilimlari.append(hz)

    sonuc['hz_anova'] = omet.anova_n_grup(*hz_dagilimlari)

    # ANOVA: CV
    cv_dagilimlari = []
    for k in kayit_listesi:
        cv = cv_pencere_dagilimi(k['spike_zaman'], k['sure_sn'])
        cv_dagilimlari.append(cv)
    sonuc['cv_anova'] = omet.anova_n_grup(*cv_dagilimlari)

    # İkili karşılaştırmalar (post-hoc niyetinde değil, ham veri)
    ikiliyer = []
    for i in range(len(kayit_listesi)):
        for j in range(i+1, len(kayit_listesi)):
            ikiliyer.append(ikili_karsilastir(
                kayit_listesi[i], kayit_listesi[j]))
    sonuc['ikili'] = ikiliyer

    return sonuc


def yazdir_ikili(sonuc, acf_esik=0.3):
    """İkili karşılaştırma sonucunu terminale basar."""
    print()
    print('=' * 64)
    print(f'  {sonuc["kayit1"]}  vs  {sonuc["kayit2"]}')
    print('=' * 64)

    oz = sonuc['ozet']
    print(f'  k1: n={oz["k1_n"]}  Hz={oout.n2s(oz["k1_hz"], 4)}  '
          f'CV={oout.n2s(oz["k1_cv"], 4)}')
    print(f'  k2: n={oz["k2_n"]}  Hz={oout.n2s(oz["k2_hz"], 4)}  '
          f'CV={oout.n2s(oz["k2_cv"], 4)}')
    print(f'  Hz fark: {oout.n2s(oz["hz_fark"], 4)}')
    print(f'  CV fark: {oout.n2s(oz["cv_fark"], 4)}')

    print(f'  Pencere: n1={sonuc["n_pencere_1"]}, n2={sonuc["n_pencere_2"]}')

    def yaz(label, t):
        if t is None:
            print(f'  {label}: NA (yetersiz veri)')
            return
        if 't_stat' in t:
            uyari = ''
            if t.get('normalite_uyari'):
                uyari = '  [DIKKAT: veri normal dagilmiyor, MW tercih]'
            print(f'  {label}: t={oout.n2s(t["t_stat"], 3)}  '
                  f'p={oout.n2s(t["p_deger"], 4)}{uyari}')
            n1 = t.get('normalite_g1')
            n2 = t.get('normalite_g2')
            if n1:
                print(f'    Shapiro g1: W={oout.n2s(n1["w_stat"], 3)}  '
                      f'p={oout.n2s(n1["p_deger"], 4)}  '
                      f'normal={n1["normal_mi"]}')
            if n2:
                print(f'    Shapiro g2: W={oout.n2s(n2["w_stat"], 3)}  '
                      f'p={oout.n2s(n2["p_deger"], 4)}  '
                      f'normal={n2["normal_mi"]}')
        else:
            print(f'  {label}: U={oout.n2s(t["u_stat"], 1)}  '
                  f'p={oout.n2s(t["p_deger"], 4)}')

    print()
    print('  -- Pencere bazli Hz --')
    if sonuc.get('hz_acf_1') is not None and sonuc.get('hz_acf_2') is not None:
        a1, a2 = sonuc['hz_acf_1'], sonuc['hz_acf_2']
        uyari = ''
        if abs(a1) > acf_esik or abs(a2) > acf_esik:
            uyari = (f'  [DIKKAT: |ACF| > {acf_esik}, t-test/MW p şişmiş olabilir]')
        print(f'  Lag-1 ACF  : k1={oout.n2s(a1, 3)}, k2={oout.n2s(a2, 3)}{uyari}')
    yaz('t-test     ', sonuc['hz_ttest'])
    yaz('Mann-Whitney', sonuc['hz_mw'])
    if sonuc.get('hz_perm'):
        p = sonuc['hz_perm']
        print(f'  Circular shift perm: p={oout.n2s(p["p_deger"], 4)}  '
              f'n={p["n_perm"]}  fark={oout.n2s(p["gercek_fark"], 4)}')

    print()
    print('  -- Pencere bazli CV --')
    if sonuc.get('cv_acf_1') is not None and sonuc.get('cv_acf_2') is not None:
        a1, a2 = sonuc['cv_acf_1'], sonuc['cv_acf_2']
        uyari = ''
        if abs(a1) > acf_esik or abs(a2) > acf_esik:
            uyari = f'  [DIKKAT: |ACF| > {acf_esik}]'
        print(f'  Lag-1 ACF  : k1={oout.n2s(a1, 3)}, k2={oout.n2s(a2, 3)}{uyari}')
    yaz('t-test     ', sonuc['cv_ttest'])
    yaz('Mann-Whitney', sonuc['cv_mw'])
    if sonuc.get('cv_perm'):
        p = sonuc['cv_perm']
        print(f'  Circular shift perm: p={oout.n2s(p["p_deger"], 4)}  '
              f'n={p["n_perm"]}  fark={oout.n2s(p["gercek_fark"], 4)}')

    print()
    print('  -- ISI dagilimi --')
    yaz('t-test     ', sonuc['isi_ttest'])
    yaz('Mann-Whitney', sonuc['isi_mw'])

    print('=' * 64)


def yazdir_coklu(coklu):
    """Çoklu karşılaştırma sonucunu terminale basar."""
    print()
    print('=' * 64)
    print(f'  COKLU KARSILASTIRMA  ({coklu["n_kayit"]} kayit)')
    print('=' * 64)
    print(f'  Kayitlar: {", ".join(coklu["kayit_etiketleri"])}')

    print()
    print('  -- ANOVA: pencere bazli Hz --')
    a = coklu['hz_anova']
    if a:
        print(f'    F={oout.n2s(a["f_stat"], 3)}  '
              f'p={oout.n2s(a["p_deger"], 4)}')
        print(f'    n_gruplar={a["n_gruplar"]}')
        print(f'    mean_gruplar={[oout.n2s(x, 4) for x in a["mean_gruplar"]]}')

    print()
    print('  -- ANOVA: pencere bazli CV --')
    a = coklu['cv_anova']
    if a:
        print(f'    F={oout.n2s(a["f_stat"], 3)}  '
              f'p={oout.n2s(a["p_deger"], 4)}')

    print()
    print('  -- Ikili karsilastirmalar --')
    for s in coklu['ikili']:
        yazdir_ikili(s)


def csv_yaz_compare(yol, sonuclar):
    """Karşılaştırma sonuçlarını CSV'ye yazar."""
    import os
    yeni = not os.path.exists(yol)
    with open(yol, 'a', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        if yeni:
            w.writerow([
                'kayit1', 'kayit2', 'k1_n', 'k2_n',
                'k1_hz', 'k2_hz', 'k1_cv', 'k2_cv',
                'hz_ttest_t', 'hz_ttest_p', 'hz_mw_u', 'hz_mw_p',
                'cv_ttest_t', 'cv_ttest_p', 'cv_mw_u', 'cv_mw_p',
                'isi_ttest_t', 'isi_ttest_p', 'isi_mw_u', 'isi_mw_p',
            ])
        for s in sonuclar:
            oz = s['ozet']
            w.writerow([
                s['kayit1'], s['kayit2'], oz['k1_n'], oz['k2_n'],
                oz['k1_hz'], oz['k2_hz'], oz['k1_cv'], oz['k2_cv'],
                (s['hz_ttest'] or {}).get('t_stat'),
                (s['hz_ttest'] or {}).get('p_deger'),
                (s['hz_mw'] or {}).get('u_stat'),
                (s['hz_mw'] or {}).get('p_deger'),
                (s['cv_ttest'] or {}).get('t_stat'),
                (s['cv_ttest'] or {}).get('p_deger'),
                (s['cv_mw'] or {}).get('u_stat'),
                (s['cv_mw'] or {}).get('p_deger'),
                (s['isi_ttest'] or {}).get('t_stat'),
                (s['isi_ttest'] or {}).get('p_deger'),
                (s['isi_mw'] or {}).get('u_stat'),
                (s['isi_mw'] or {}).get('p_deger'),
            ])


def main():
    parser = argparse.ArgumentParser(description='Organoid karsilastirma')
    parser.add_argument('kayitlar', nargs='+',
                         help='Cift olarak: organoid1 kayit1 organoid2 kayit2 ...')
    parser.add_argument('--csv', default='compare_results.csv')
    parser.add_argument('--acf-esik', type=float, default=0.3,
                         dest='acf_esik',
                         help='Lag-1 ACF uyarı eşiği (default 0.3, '
                              'sıkı analiz için 0.1)')
    parser.add_argument('--no-perm', action='store_true',
                         help='Circular shift permutation atla')
    args = parser.parse_args()

    if len(args.kayitlar) % 2 != 0 or len(args.kayitlar) < 4:
        print('HATA: Cift sayida arguman ver. En az iki kayit (4 arguman).')
        print('Ornek: python organoid_compare.py SO1 13_Haz SO1 14_Haz')
        return 1

    # (organoid, kayit) çiftleri
    ciftler = []
    for i in range(0, len(args.kayitlar), 2):
        ciftler.append((args.kayitlar[i], args.kayitlar[i+1]))

    conn = odb.db_ac()

    # Tüm kayıtları yükle
    kayit_listesi = []
    for org, kay in ciftler:
        k = kayit_metriklerini_al(conn, org, kay)
        if k is None:
            print(f'HATA: {org}/{kay} veritabaninda yok.')
            conn.close()
            return 1
        if len(k['spike_zaman']) == 0:
            print(f'HATA: {org}/{kay} icin spike yok.')
            conn.close()
            return 1
        kayit_listesi.append(k)

    if len(kayit_listesi) == 2:
        sonuc = ikili_karsilastir(kayit_listesi[0], kayit_listesi[1],
                                    circular_perm=not args.no_perm)
        yazdir_ikili(sonuc, acf_esik=args.acf_esik)
        csv_yaz_compare(args.csv, [sonuc])
    else:
        coklu = coklu_karsilastir(kayit_listesi)
        yazdir_coklu(coklu)
        csv_yaz_compare(args.csv, coklu['ikili'])

    print(f'\n  CSV: {args.csv}')
    conn.close()
    return 0


if __name__ == '__main__':
    sys.exit(main())
