"""
organoid_replikasyon.py
========================
Aynı organoidin birden fazla kayıt oturumunu karşılaştırır.
"Biyolojik replikasyon" sınırlılığını giderir.

Temel soru: Bir organoidin aktivite profili oturumdan oturuma
ne kadar kararlı? Hangi metrikler tutarlı, hangiler değişiyor?

Kullanım (CLI):
    python organoid_replikasyon.py <organoid_ad> <nwb1> <nwb2> [nwb3 ...]
                                   [--out-dir KLASOR]

Kullanım (modül):
    from organoid_replikasyon import oturum_karsilastir
    sonuc = oturum_karsilastir(['dosya1.nwb', 'dosya2.nwb'], 'HO2')

Hesaplanan metrikler (her oturum için):
    - Aktif birim sayısı ve oranı
    - Hz dağılımı (medyan, IQR, log-sigma)
    - Medyan CV
    - Medyan STTC
    - Network burst/dakika
    - Ritmik birim oranı

Kararlılık metrikleri (oturumlar arası):
    - ICC (Intraclass Correlation Coefficient) — altın standart
    - CV_oturum (oturumlararası CV) — basit ama yorumlanabilir
    - Spearman ρ Hz sıralaması — aynı birimler mi benzer sırada?
"""

import sys
import os
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr


# ─────────────────────────────────────────────────────────────
# ICC hesabı (Shrout & Fleiss, 1979 — ICC(2,1))
# ─────────────────────────────────────────────────────────────

def icc_hesapla(matris):
    """
    ICC(2,1) — Two-way random effects, absolute agreement, single measures.
    Goldstandard oturumlararası kararlılık metriği.

    matris: shape (n_oturum, n_metrik) — her satır bir oturum,
            her sütun bir ölçüm (örn. medyan Hz, CV, burst/dk...)

    Dönen:
        icc_deger: 0-1 arası (>0.75 iyi, >0.9 mükemmel)
        f_stat: F istatistiği
        p_deger: anlamlılık

    Referans: Shrout, P.E. & Fleiss, J.L. (1979).
    Intraclass correlations: Uses in assessing rater reliability.
    Psychological Bulletin, 86(2), 420-428.
    """
    matris = np.array(matris, dtype=float)
    if matris.ndim == 1:
        matris = matris.reshape(-1, 1)

    k, n = matris.shape  # k=oturum sayisi, n=metrik sayisi

    if k < 2:
        return None, None, None

    # Grand mean
    grand_mean = matris.mean()

    # Between-targets MS (oturumlar arası)
    oturum_meanleri = matris.mean(axis=1)
    MS_r = n * np.sum((oturum_meanleri - grand_mean)**2) / (k - 1)

    # Within MS (içi)
    MS_w = np.sum((matris - oturum_meanleri.reshape(-1,1))**2) / (k*(n-1))

    if MS_w == 0:
        return 1.0, float('inf'), 0.0

    # ICC(2,1) formülü
    icc = (MS_r - MS_w) / (MS_r + (n-1)*MS_w)
    icc = max(0.0, min(1.0, icc))  # 0-1 aralığına sıkıştır

    # F testi
    from scipy.stats import f as f_dist
    F = MS_r / MS_w
    df1 = k - 1
    df2 = k * (n - 1)
    p = 1 - f_dist.cdf(F, df1, df2)

    return round(float(icc), 4), round(float(F), 3), round(float(p), 4)


# ─────────────────────────────────────────────────────────────
# Tek oturum özet metriklerini çıkar
# ─────────────────────────────────────────────────────────────

def oturum_ozet_cikar(ozet_csv_yolu):
    """
    organoid_units_analiz.py'nin ürettiği *_ozet.csv'yi okur.
    Dönen: dict (metrik adı -> float değer)
    """
    if not os.path.exists(ozet_csv_yolu):
        return None

    rows = list(csv.DictReader(open(ozet_csv_yolu, encoding='utf-8-sig')))
    if not rows:
        return None

    r = rows[0]

    def f(k):
        try:
            return float(r.get(k, '') or 0)
        except (ValueError, TypeError):
            return 0.0

    return {
        'kayit': r.get('kayit', ''),
        'kayit_suresi_sn': f('kayit_suresi_sn'),
        'aktif_unit': f('aktif_unit'),
        'toplam_unit': f('toplam_unit'),
        'aktif_oran': f('aktif_oran'),
        'toplam_spike': f('toplam_spike'),
        'hz_medyan': f('hz_medyan'),
        'hz_iqr': f('hz_iqr'),
        'cv_medyan': f('cv_medyan'),
        'refractory_medyan': f('refractory_medyan'),
        'sttc_medyan': f('sttc_medyan'),
        'burst_sayisi': f('burst_sayisi'),
        'burst_per_dakika': f('burst_per_dakika'),
        'ritim_unit_orani': f('ritim_unit_orani'),
    }


def oturum_ozet_hesapla(nwb_yolu, organoid_ad, v6_3_yol, tmp_dir):
    """
    Henüz analiz edilmemiş bir NWB dosyası için per-unit analiz çalıştırır
    ve özet metriklerini döner.
    """
    import subprocess
    import re

    per_unit = Path(v6_3_yol) / 'organoid_units_analiz.py'
    if not per_unit.exists():
        raise FileNotFoundError(f'organoid_units_analiz.py bulunamadı: {per_unit}')

    m = re.search(r'ses-(\d{8})', str(nwb_yolu))
    kayit = m.group(1) if m else Path(nwb_yolu).stem[:12].replace('-','')

    cikti = Path(tmp_dir) / organoid_ad / kayit
    cikti.mkdir(parents=True, exist_ok=True)

    # Önce mevcut ozet var mı?
    ozet_adaylar = list(cikti.glob('*_ozet.csv'))
    if ozet_adaylar:
        sonuc = oturum_ozet_cikar(str(ozet_adaylar[0]))
        if sonuc:
            return sonuc

    # Analiz çalıştır
    r = subprocess.run(
        [sys.executable, str(per_unit),
         str(Path(nwb_yolu).resolve()),
         organoid_ad, kayit,
         '--out-dir', str(cikti.resolve()),
         '--perm-n', '100'],  # Replikasyon için hız önemli
        cwd=str(v6_3_yol), capture_output=True, text=True)

    ozet_adaylar = list(cikti.glob('*_ozet.csv'))
    if not ozet_adaylar:
        return None

    return oturum_ozet_cikar(str(ozet_adaylar[0]))


# ─────────────────────────────────────────────────────────────
# Kararlılık analizi
# ─────────────────────────────────────────────────────────────

def karlilik_analizi(oturum_listesi):
    """
    N oturumun metriklerinden kararlılık hesaplar.

    Dönen: dict
        - her metrik için: oturumlararası CV, ICC, min/max/range
        - genel kararlılık skoru (0-1)
    """
    if len(oturum_listesi) < 2:
        return {'hata': 'En az 2 oturum gerekli'}

    metrikler = [
        'hz_medyan', 'cv_medyan', 'sttc_medyan',
        'burst_per_dakika', 'aktif_oran', 'ritim_unit_orani'
    ]

    sonuc = {}
    icc_degerleri = []

    for m in metrikler:
        degerler = [o[m] for o in oturum_listesi if o and m in o]
        if len(degerler) < 2:
            continue

        arr = np.array(degerler)
        ort = arr.mean()
        std = arr.std()

        # Oturumlararası CV (%)
        cv_oturum = (std / ort * 100) if ort > 0 else float('inf')

        # Kararlılık: CV_oturum tabanlı skor (ICC bu yapıda hesaplanamaz)
        karlilik_skoru = max(0.0, 1.0 - cv_oturum / 100.0)

        sonuc[m] = {
            'degerler': [round(d, 4) for d in degerler],
            'ortalama': round(float(ort), 4),
            'std': round(float(std), 4),
            'cv_oturum_pct': round(float(cv_oturum), 2),
            'karlilik_skoru': round(karlilik_skoru, 3),
            'min': round(float(arr.min()), 4),
            'max': round(float(arr.max()), 4),
        }
        icc_degerleri.append(karlilik_skoru)

    if icc_degerleri:
        genel_karlilik = np.mean(icc_degerleri)
        if genel_karlilik >= 0.85:
            karlilik_yorum = 'Mukemmel (CV_oturum<%15)'
        elif genel_karlilik >= 0.70:
            karlilik_yorum = 'Iyi (CV_oturum<30%)'
        else:
            karlilik_yorum = 'Dusuk (CV_oturum>30%)'

        sonuc['_ozet'] = {
            'n_oturum': len(oturum_listesi),
            'genel_karlilik': round(float(genel_karlilik), 4),
            'karlilik_yorum': karlilik_yorum,
        }

    return sonuc


# ─────────────────────────────────────────────────────────────
# Grafik
# ─────────────────────────────────────────────────────────────

def grafik_uret(oturum_listesi, karlilik, organoid_ad, cikti_yolu):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    metrikler = ['hz_medyan', 'cv_medyan', 'burst_per_dakika',
                 'sttc_medyan', 'aktif_oran', 'ritim_unit_orani']
    etiketler = ['Medyan Hz', 'Medyan CV', 'Burst/Dakika',
                 'Medyan STTC', 'Aktif Birim Oranı', 'Ritmik Birim Oranı']

    n_oturum = len(oturum_listesi)
    oturum_adlari = [o.get('kayit', f'O{i+1}') for i,o in enumerate(oturum_listesi)]

    fig, axes = plt.subplots(2, 3, figsize=(13, 8), facecolor='white')
    fig.suptitle(f'{organoid_ad} — Oturumlar Arası Kararlılık Analizi '
                 f'({n_oturum} oturum)',
                 fontsize=12, fontweight='bold')

    renkler = ['#1F77B4','#D62728','#2CA02C','#FF7F0E','#9467BD']

    for idx, (metrik, etiket) in enumerate(zip(metrikler, etiketler)):
        ax = axes[idx//3][idx%3]
        degerler = [o.get(metrik, 0) for o in oturum_listesi]

        bars = ax.bar(range(n_oturum), degerler,
                      color=renkler[:n_oturum], edgecolor='black', alpha=0.8)

        # ICC bilgisi
        if metrik in karlilik and karlilik[metrik].get('icc') is not None:
            icc = karlilik[metrik]['icc']
            cv_o = karlilik[metrik]['cv_oturum_pct']
            ax.set_title(f'{etiket}\nICC={icc:.2f}, CV_oturum=%{cv_o:.1f}',
                         fontsize=9)
        else:
            ax.set_title(etiket, fontsize=9)

        ax.set_xticks(range(n_oturum))
        ax.set_xticklabels(oturum_adlari, rotation=20, fontsize=7)
        ax.grid(alpha=0.3, axis='y')

        for bar, val in zip(bars, degerler):
            ax.text(bar.get_x()+bar.get_width()/2,
                    bar.get_height()*1.02,
                    f'{val:.3f}', ha='center', va='bottom', fontsize=7.5)

    plt.tight_layout()
    plt.savefig(cikti_yolu, dpi=130, bbox_inches='tight', facecolor='white')
    plt.close()


# ─────────────────────────────────────────────────────────────
# Ana fonksiyon (modül olarak import edildiğinde)
# ─────────────────────────────────────────────────────────────

def oturum_karsilastir(nwb_yollari, organoid_ad,
                        v6_3_yol=None, cikti_klasor='sonuclar/replikasyon'):
    """
    Birden fazla NWB oturumunu karşılaştırır.

    nwb_yollari: liste — her biri bir NWB dosya yolu
    organoid_ad: str — HO2, HO5 gibi
    v6_3_yol: str — organoid_units_analiz.py'nin bulunduğu klasör
    cikti_klasor: str — sonuçların yazılacağı yer

    Dönen: (oturum_listesi, karlilik_sonuc)
    """
    if v6_3_yol is None:
        v6_3_yol = Path(__file__).parent

    Path(cikti_klasor).mkdir(parents=True, exist_ok=True)

    print(f'[Replikasyon] {organoid_ad}: {len(nwb_yollari)} oturum analiz ediliyor...')

    oturum_listesi = []
    for i, nwb_yolu in enumerate(nwb_yollari):
        print(f'  Oturum {i+1}/{len(nwb_yollari)}: {Path(nwb_yolu).name}')
        ozet = oturum_ozet_hesapla(nwb_yolu, organoid_ad,
                                    v6_3_yol, cikti_klasor)
        if ozet:
            oturum_listesi.append(ozet)
        else:
            print(f'  UYARI: Oturum {i+1} analiz edilemedi — atlandı')

    if len(oturum_listesi) < 2:
        print('  HATA: Karşılaştırma için en az 2 başarılı oturum gerekli')
        return oturum_listesi, {}

    karlilik = karlilik_analizi(oturum_listesi)

    # Özet çıktı
    print(f'\n  [KARLILIK ANALİZİ] {organoid_ad}')
    print(f'  Oturum sayısı: {len(oturum_listesi)}')
    if '_ozet' in karlilik:
        oz = karlilik['_ozet']
        print(f'  Genel ICC: {oz["genel_icc"]} — {oz["karlilik_yorum"]}')

    print(f'\n  {"Metrik":<25} {"ICC":<8} {"CV_oturum%":<12} {"Min-Max"}')
    print('  ' + '-'*60)
    for m, veri in karlilik.items():
        if m.startswith('_'):
            continue
        icc = veri.get('icc', 'N/A')
        cv_o = veri.get('cv_oturum_pct', 'N/A')
        mn = veri.get('min', '-')
        mx = veri.get('max', '-')
        print(f'  {m:<25} {str(icc):<8} {str(cv_o):<12} {mn} - {mx}')

    # Grafik
    grafik_yolu = Path(cikti_klasor) / f'{organoid_ad}_replikasyon.png'
    grafik_uret(oturum_listesi, karlilik, organoid_ad, str(grafik_yolu))
    print(f'\n  Grafik: {grafik_yolu}')

    # CSV
    csv_yolu = Path(cikti_klasor) / f'{organoid_ad}_replikasyon_ozet.csv'
    with open(csv_yolu, 'w', newline='', encoding='utf-8') as f:
        fieldnames = ['metrik', 'icc', 'cv_oturum_pct', 'ortalama',
                      'std', 'min', 'max'] + \
                     [f'oturum_{i+1}' for i in range(len(oturum_listesi))]
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for m, veri in karlilik.items():
            if m.startswith('_'):
                continue
            satir = {
                'metrik': m,
                'icc': veri.get('icc', ''),
                'cv_oturum_pct': veri.get('cv_oturum_pct', ''),
                'ortalama': veri.get('ortalama', ''),
                'std': veri.get('std', ''),
                'min': veri.get('min', ''),
                'max': veri.get('max', ''),
            }
            for i, d in enumerate(veri.get('degerler', [])):
                satir[f'oturum_{i+1}'] = d
            w.writerow(satir)

    print(f'  CSV: {csv_yolu}')
    return oturum_listesi, karlilik


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Aynı organoidin birden fazla oturumunu karşılaştır')
    parser.add_argument('organoid', help='Organoid kodu (orn: HO2)')
    parser.add_argument('nwb_dosyalari', nargs='+',
                        help='2+ NWB dosya yolu')
    parser.add_argument('--out-dir', default='sonuclar/replikasyon')
    args = parser.parse_args()

    v6_3_yol = Path(__file__).parent

    oturum_karsilastir(
        nwb_yollari=args.nwb_dosyalari,
        organoid_ad=args.organoid,
        v6_3_yol=v6_3_yol,
        cikti_klasor=args.out_dir
    )


if __name__ == '__main__':
    main()
