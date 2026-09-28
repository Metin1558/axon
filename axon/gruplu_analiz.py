"""
gruplu_analiz.py — Axon v6.3
==============================
sonuclar/ klasöründeki tüm özet CSV'leri otomatik bulur,
kayıt süresine ve diğer kriterlere göre gruplar,
grup içi variabilite ve gruplar arası karşılaştırma yapar.

Kullanım:
    python gruplu_analiz.py                    # hepsini analiz et
    python gruplu_analiz.py --metric hz_medyan # belirli metrik
    python gruplu_analiz.py --sure-esik 300    # süre eşiği (sn)
    python gruplu_analiz.py --csv              # CSV kaydet

Otomatik gruplama:
    - Kayıt süresine göre (kısa: <300 sn, uzun: >=300 sn)
    - Türe göre (species alanı varsa)
    - İsme göre (HO*, MO*, M* prefix)
"""

import sys
import csv
import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict


# ── CSV okuma ─────────────────────────────────────────────

def ozet_csv_listesi_bul(base_dir='sonuclar'):
    """sonuclar/ altındaki tüm *_ozet.csv dosyalarını bulur."""
    base = Path(base_dir)
    if not base.exists():
        return []
    return sorted(base.rglob('*_ozet.csv'))


def csv_oku(dosya):
    """Özet CSV'yi dict olarak okur."""
    try:
        with open(dosya, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        if rows:
            r = rows[0]
            # sayısal dönüşüm
            for k, v in r.items():
                try:
                    r[k] = float(v)
                except (ValueError, TypeError):
                    pass
            # subject tespiti: klasör adı
            r['subject'] = dosya.parent.name
            return r
    except Exception as e:
        print(f"  UYARI: {dosya.name} okunamadı: {e}")
    return None


# ── Otomatik gruplama ─────────────────────────────────────

def gruplara_ayir(kayitlar, sure_esik=300):
    """
    Kayıtları otomatik gruplara ayırır.
    Öncelik sırası: prefix (HO/MO/M) > süre > tek grup
    """
    gruplar = defaultdict(list)

    for r in kayitlar:
        subj = r.get('subject', '')
        sure = float(r.get('kayit_suresi_sn', 0))

        # Prefix bazlı gruplama
        if subj.startswith('HO'):
            # İnsan organoidi — süreye göre alt grup
            alt_grup = 'HO_kisa' if sure < sure_esik else 'HO_uzun'
            gruplar[alt_grup].append(r)
        elif subj.startswith('MO'):
            gruplar['MO'].append(r)
        elif subj.upper().startswith('M') and not subj.startswith('MO'):
            gruplar['ExVivo'].append(r)
        else:
            # Bilinmeyen — süreye göre grupla
            alt_grup = 'Kisa' if sure < sure_esik else 'Uzun'
            gruplar[alt_grup].append(r)

    return dict(gruplar)


# ── İstatistik ────────────────────────────────────────────

def grup_istatistik(kayitlar, metrik):
    """Bir gruptaki kayıtlar için metrik istatistiği."""
    degerler = []
    for r in kayitlar:
        v = r.get(metrik)
        if v is not None:
            try:
                degerler.append(float(v))
            except (ValueError, TypeError):
                pass

    if not degerler:
        return None

    return {
        'n': len(degerler),
        'ortalama': np.mean(degerler),
        'medyan': np.median(degerler),
        'std': np.std(degerler),
        'min': np.min(degerler),
        'max': np.max(degerler),
        'cv_bireyler': np.std(degerler) / np.mean(degerler) if np.mean(degerler) != 0 else 0,
    }


def mann_whitney(a, b):
    """Basit Mann-Whitney U testi."""
    try:
        from scipy.stats import mannwhitneyu
        stat, p = mannwhitneyu(a, b, alternative='two-sided')
        return float(p)
    except Exception:
        return None


# ── Ana program ───────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Axon v6.3 — Gruplu Analiz')
    parser.add_argument('--sonuclar', type=str, default='sonuclar',
                        help='Sonuçlar klasörü (default: sonuclar)')
    parser.add_argument('--metric', type=str, default=None,
                        help='Analiz edilecek metrik (default: hepsi)')
    parser.add_argument('--sure-esik', type=float, default=300,
                        help='Kısa/uzun kayıt eşiği saniye (default: 300)')
    parser.add_argument('--csv', action='store_true',
                        help='Sonuçları gruplu_analiz.csv kaydet')
    args = parser.parse_args()

    print("=" * 65)
    print("  Axon v6.3 — Gruplu Analiz")
    print("=" * 65)

    csv_listesi = ozet_csv_listesi_bul(args.sonuclar)
    if not csv_listesi:
        print(f"\n  UYARI: {args.sonuclar}/ içinde özet CSV bulunamadı.")
        print("  Önce organoid_cli.py ile analiz yap.")
        return 1

    kayitlar = [csv_oku(f) for f in csv_listesi]
    kayitlar = [r for r in kayitlar if r is not None]

    print(f"\n  {len(kayitlar)} subject bulundu.")

    gruplar = gruplara_ayir(kayitlar, args.sure_esik)
    print(f"  {len(gruplar)} grup oluşturuldu: {', '.join(gruplar.keys())}\n")

    # Analiz edilecek metrikler
    if args.metric:
        metrikler = [args.metric]
    else:
        metrikler = ['hz_medyan', 'cv_medyan', 'sttc_medyan',
                     'burst_per_dakika', 'aktif_oran']

    # Her metrik için grup karşılaştırması
    cikti_satirlar = []

    for metrik in metrikler:
        print(f"  {'─'*55}")
        print(f"  Metrik: {metrik}")
        print(f"  {'─'*55}")
        print(f"  {'Grup':<15} {'n':>4} {'Medyan':>10} {'Ort±Std':>18} {'CV_bireyler':>12}")

        grup_degerler = {}
        for grup_adi, kayit_listesi in sorted(gruplar.items()):
            istat = grup_istatistik(kayit_listesi, metrik)
            if not istat:
                continue
            grup_degerler[grup_adi] = [
                float(r.get(metrik, 0)) for r in kayit_listesi
                if r.get(metrik) is not None
            ]
            print(f"  {grup_adi:<15} {istat['n']:>4} "
                  f"{istat['medyan']:>10.4f} "
                  f"{istat['ortalama']:>8.4f}±{istat['std']:<8.4f} "
                  f"{istat['cv_bireyler']:>12.3f}")

            for r in kayit_listesi:
                cikti_satirlar.append({
                    'metrik': metrik,
                    'grup': grup_adi,
                    'subject': r.get('subject', ''),
                    'deger': r.get(metrik, ''),
                })

        # Gruplar arası karşılaştırma (en fazla 2 grup)
        grup_adlari = list(grup_degerler.keys())
        if len(grup_adlari) == 2:
            a, b = grup_degerler[grup_adlari[0]], grup_degerler[grup_adlari[1]]
            p = mann_whitney(a, b)
            if p is not None:
                anlamli = "* p<0.05" if p < 0.05 else "ns"
                print(f"\n  Mann-Whitney [{grup_adlari[0]} vs {grup_adlari[1]}]: "
                      f"p={p:.4f} {anlamli}")
        print()

    if args.csv:
        cikti = Path('gruplu_analiz_sonuc.csv')
        with open(cikti, 'w', newline='', encoding='utf-8') as f:
            w = csv.DictWriter(f, fieldnames=['metrik', 'grup', 'subject', 'deger'])
            w.writeheader()
            w.writerows(cikti_satirlar)
        print(f"  CSV kaydedildi: {cikti}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
