"""
graph_analiz.py — Axon v6.3
============================
STTC matrislerinden graph theory metrikleri hesaplar.
sonuclar/ klasöründeki tüm analiz sonuçlarını otomatik bulur.

Kullanım:
    python graph_analiz.py                    # sonuclar/ içindeki hepsini analiz et
    python graph_analiz.py --subject HO2      # sadece tek subject
    python graph_analiz.py --dataset 001603   # sadece belirli dataset
    python graph_analiz.py --csv              # sonuçları CSV'e yaz

Gereksinimler:
    pip install networkx numpy scipy
"""

import sys
import os
import csv
import argparse
import subprocess
from pathlib import Path

# networkx kontrol
try:
    import networkx as nx
except ImportError:
    print("networkx kuruluyor...")
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--quiet', 'networkx'],
                   check=True)
    import networkx as nx

import numpy as np

sys.path.insert(0, str(Path(__file__).parent / 'v6_3'))
import organoid_io as oio
import organoid_metrics as omet


# ── Graph theory fonksiyonları ────────────────────────────

def sttc_matrisi_hesapla(spike_listesi, sure_sn, dt=0.01, max_unit=20):
    """En aktif N unitin STTC matrisini hesaplar."""
    spike_sayisi = [len(sp) for sp in spike_listesi]
    idx = sorted(range(len(spike_sayisi)),
                 key=lambda i: spike_sayisi[i], reverse=True)[:max_unit]
    secili = [spike_listesi[i] for i in idx]
    n = len(secili)
    matris = np.eye(n)
    for i in range(n):
        for j in range(i+1, n):
            v = omet.sttc_hesapla(secili[i], secili[j], sure_sn, dt=dt)
            matris[i, j] = v
            matris[j, i] = v
    return matris, idx


def graph_metrikleri_hesapla(matris, esik=0.1):
    """STTC matrisinden graph theory metrikleri."""
    n = len(matris)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i+1, n):
            if matris[i, j] > esik:
                G.add_edge(i, j, weight=matris[i, j])

    yoğunluk = nx.density(G)

    if G.number_of_edges() == 0:
        return {
            'yogunluk': 0.0,
            'ortalama_kumeleme': 0.0,
            'ortalama_yol_uzunlugu': None,
            'kucuk_dunya': None,
            'sigma': None,
            'baglantili': False,
            'kenar_sayisi': 0,
        }

    C = nx.average_clustering(G)
    baglantili = nx.is_connected(G)

    if baglantili and n > 2:
        L = nx.average_shortest_path_length(G)
        # Random graph referans
        if yoğunluk > 0 and yoğunluk < 1:
            C_rand = yoğunluk
            L_rand = (np.log(n) / np.log(max(2, n * yoğunluk)))
            sigma = (C / C_rand) / (L / L_rand) if L_rand > 0 else None
        else:
            sigma = None
        kucuk_dunya = sigma is not None and sigma > 1
    else:
        L = None
        sigma = None
        kucuk_dunya = None

    return {
        'yogunluk': round(yoğunluk, 4),
        'ortalama_kumeleme': round(C, 4),
        'ortalama_yol_uzunlugu': round(L, 4) if L else None,
        'kucuk_dunya': kucuk_dunya,
        'sigma': round(sigma, 4) if sigma else None,
        'baglantili': baglantili,
        'kenar_sayisi': G.number_of_edges(),
    }


# ── NWB dosyalarını analiz et ─────────────────────────────

def nwb_analiz(dosya_yolu, subject_id):
    """NWB dosyasından spike okur, graph metrikleri hesaplar."""
    try:
        units = oio.sorted_units_oku(dosya_yolu)
    except Exception:
        try:
            units = oio.sorted_units_read(dosya_yolu)
        except Exception as e:
            print(f"  HATA okuma: {e}")
            return None

    if not units:
        print(f"  Sorted units yok: {subject_id}")
        return None

    spike_listesi = [np.array(u['spike_zaman']) for u in units
                     if len(u['spike_zaman']) >= 10]

    if len(spike_listesi) < 2:
        print(f"  Yetersiz aktif unit: {subject_id}")
        return None

    # Kayıt süresi
    sure_sn = max(sp.max() for sp in spike_listesi if len(sp) > 0)

    matris, _ = sttc_matrisi_hesapla(spike_listesi, sure_sn)
    metriks = graph_metrikleri_hesapla(matris)
    metriks['subject'] = subject_id
    metriks['n_aktif_unit'] = len(spike_listesi)
    metriks['sure_sn'] = round(sure_sn, 1)
    return metriks


def ozet_csv_oku(csv_yolu):
    """Mevcut özet CSV'den spike listesi olmadan basit metrikleri okur."""
    try:
        with open(csv_yolu, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            satirlar = list(reader)
        if satirlar:
            return satirlar[0]
    except Exception:
        pass
    return None


# ── Otomatik subject tespiti ──────────────────────────────

def subject_listesi_bul(base_dir='sonuclar', subject_filtresi=None,
                        dataset_filtresi=None, indir_dir='dandi_indirilmis'):
    """sonuclar/ ve dandi_indirilmis/ klasörlerinden subject listesini toplar."""
    base = Path(base_dir)
    indir = Path(indir_dir)
    subjects = {}

    # sonuclar/ klasöründen ozet CSV'leri bul
    if base.exists():
        for ozet in base.rglob('*_ozet.csv'):
            subj = ozet.parent.name
            if subject_filtresi and subj != subject_filtresi:
                continue
            subjects[subj] = {'ozet_csv': ozet, 'nwb': None}

    # dandi_indirilmis/ klasöründen NWB dosyaları bul
    if indir.exists():
        for nwb in indir.glob('*.nwb'):
            # sub-HO2_ses-... -> HO2
            parts = nwb.stem.split('_')[0].replace('sub-', '')
            subj = parts
            if subject_filtresi and subj != subject_filtresi:
                continue
            if subj not in subjects:
                subjects[subj] = {'ozet_csv': None, 'nwb': nwb}
            else:
                subjects[subj]['nwb'] = nwb

    return subjects


# ── Ana program ───────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Axon v6.3 — Graph Theory Analizi')
    parser.add_argument('--subject', type=str, default=None,
                        help='Sadece bu subject (örn: HO2)')
    parser.add_argument('--dataset', type=str, default=None,
                        help='Dataset filtresi (şimdilik bilgi amaçlı)')
    parser.add_argument('--esik', type=float, default=0.1,
                        help='STTC eşik değeri (default: 0.1)')
    parser.add_argument('--csv', action='store_true',
                        help='Sonuçları graph_metrikleri.csv olarak kaydet')
    parser.add_argument('--sonuclar', type=str, default='sonuclar',
                        help='Sonuçlar klasörü (default: sonuclar)')
    parser.add_argument('--indir', type=str, default='dandi_indirilmis',
                        help='İndirilen NWB klasörü (default: dandi_indirilmis)')
    args = parser.parse_args()

    print("=" * 65)
    print("  Axon v6.3 — Graph Theory Analizi")
    print("=" * 65)

    subjects = subject_listesi_bul(
        args.sonuclar, args.subject, args.dataset, args.indir)

    if not subjects:
        print(f"\n  UYARI: {args.sonuclar}/ veya {args.indir}/ içinde")
        print("  analiz edilebilir subject bulunamadı.")
        print("\n  Önce organoid_cli.py ile analiz yap.")
        return 1

    print(f"\n  {len(subjects)} subject bulundu: {', '.join(sorted(subjects.keys()))}")
    print(f"  STTC eşiği: {args.esik}\n")

    tum_metriks = []

    for subj in sorted(subjects.keys()):
        info = subjects[subj]
        print(f"  [{subj}]")

        metriks = None

        # NWB varsa direkt analiz et
        if info['nwb'] and info['nwb'].exists():
            print(f"    NWB: {info['nwb'].name}")
            metriks = nwb_analiz(str(info['nwb']), subj)
        elif info['ozet_csv'] and info['ozet_csv'].exists():
            # Sadece özet CSV varsa — NWB olmadan kısmi bilgi
            print(f"    Ozet CSV: {info['ozet_csv'].name}")
            print(f"    (NWB yok — graph analizi için NWB gerekli)")
            metriks = None

        if metriks:
            sigma_str = f"{metriks['sigma']}" if metriks['sigma'] else "N/A"
            sw_str = "Evet" if metriks['kucuk_dunya'] else "Hayır"
            print(f"    Yoğunluk  : {metriks['yogunluk']:.3f}")
            print(f"    Kümeleme  : {metriks['ortalama_kumeleme']:.3f}")
            if metriks['ortalama_yol_uzunlugu']:
                print(f"    Yol uzunl : {metriks['ortalama_yol_uzunlugu']:.3f}")
            print(f"    σ (sigma) : {sigma_str}")
            print(f"    Küçük dünya: {sw_str}")
            print(f"    Kenar     : {metriks['kenar_sayisi']}")
            tum_metriks.append(metriks)
        else:
            print("    Graph analizi yapılamadı.")
        print()

    # Özet tablo
    if tum_metriks:
        print("=" * 65)
        print("  ÖZET")
        print("=" * 65)
        print(f"  {'Subject':<10} {'Yoğunluk':>10} {'Kümeleme':>10} "
              f"{'σ':>8} {'Küçük Dünya':>12}")
        print("  " + "-" * 53)
        for m in tum_metriks:
            sigma_s = f"{m['sigma']:.3f}" if m['sigma'] else "N/A"
            sw_s = "✓" if m['kucuk_dunya'] else "—"
            print(f"  {m['subject']:<10} {m['yogunluk']:>10.3f} "
                  f"{m['ortalama_kumeleme']:>10.3f} {sigma_s:>8} {sw_s:>12}")

        if args.csv:
            cikti = Path('graph_metrikleri.csv')
            fieldnames = ['subject', 'n_aktif_unit', 'sure_sn', 'yogunluk',
                          'ortalama_kumeleme', 'ortalama_yol_uzunlugu',
                          'sigma', 'kucuk_dunya', 'kenar_sayisi', 'baglantili']
            with open(cikti, 'w', newline='', encoding='utf-8') as f:
                w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
                w.writeheader()
                w.writerows(tum_metriks)
            print(f"\n  CSV kaydedildi: {cikti}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
