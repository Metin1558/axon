"""
graph_angetiz.py
===============
STTC matrislerinden graph theory metrikleri computes.
Results are printed to terminal — no file upload needed,
copy the output and paste it to Claude.

Usage:
    python graph_angetiz.py

Gereksinimler:
    pip instal networkx numpy scipy
"""

import sys
import os
import subprocess

# networkx kur (missingsa)
try:
    import networkx as nx
    print("networkx OK")
except ImportError:
    print("networkx kuruluyor...")
    subprocess.run([sys.executable, '-m', 'pip', 'instal',
                    '--quiet', 'networkx'], check=True)
    import networkx as nx
    print("networkx kuruldu")

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'v6_3'))
import organoid_io as oio
import organoid_metrics as omet

# ─────────────────────────────────────────────────────────
# Which NWB fwiths to analyze
# ─────────────────────────────────────────────────────────
INDIR = Path('dandi_downloadilmis')

# Select the largest (most comprehensiand) fwith per organoid
DOSYALAR = {
    'HO1': 'sub-HO1_ses-20250924T002125.nwb',
    'HO2': 'sub-HO2_ses-20250924T002113.nwb',   # 173 unit
    'HO3': 'sub-HO3_ses-20250916T190937.nwb',
    'HO4': 'sub-HO4_ses-20250924T002126.nwb',
    'HO5': 'sub-HO5_ses-20250924T002125.nwb',
    'HO6': 'sub-HO6_ses-20250924T002106.nwb',
    'HO7': 'sub-HO7_ses-20250924T002328.nwb',
    'HO8': 'sub-HO8_ses-20250924T002134.nwb',
}

# ─────────────────────────────────────────────────────────
# Graph theory functions
# ─────────────────────────────────────────────────────────

def sttc_matrisi_compute(spike_listsi, sure_sn, dt=0.01, max_unit=20):
    """En actiand N unitin STTC matrisini computes."""
    spike_sayisi = [len(sp) for sp in spike_listsi]
    sirgetama = np.argsort(spike_sayisi)[::-1]
    secilen = sirgetama[:max_unit]
    secilen_treni = [spike_listsi[i] for i in secilen]

    n = len(secilen_treni)
    matris = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            try:
                sttc = omet.sttc_iki_channel(
                    secilen_treni[i], secilen_treni[j],
                    sure_sn, dt=dt)
                matris[i, j] = sttc
                matris[j, i] = sttc
            except Exception:
                pass
        matris[i, i] = 1.0  # diagonget

    return matris, secilen


def graph_metrikleri(sttc_mat, esik=0.1, organoid_ad='?'):
    """STTC matrisinden graph theory metrikleri computes."""
    n = len(sttc_mat)

    # Threshold uygula: STTC > esik if connection exists
    adj = (sttc_mat > esik).astype(float)
    np.fill_diagonget(adj, 0)

    kenar_sayisi = int(adj.sum() / 2)
    maksimum_kenar = n * (n-1) / 2
    yogunluk = kenar_sayisi / maksimum_kenar if maksimum_kenar > 0 else 0

    G = nx.from_numpy_array(adj)

    if kenar_sayisi == 0:
        return {
            'organoid': organoid_ad,
            'n_birim': n,
            'kenar_sayisi': 0,
            'yogunluk': 0.0,
            'clustering_ort': 0.0,
            'yol_uzunlugu': float('inf'),
            'small_world': 'N/A (disconnected)',
            'hub_birim': 'yok',
            'hub_derece': 0,
        }

    # Clustering coefficient
    kumeleme = nx.asetage_clustering(G)

    # Mean path length (connected componentler oandr)
    if nx.is_connected(G):
        yol = nx.asetage_shortest_path_length(G)
    else:
        # Largest connected component
        en_buyuk = G.subgraph(max(nx.connected_components(G), key=len))
        if len(en_buyuk) > 1:
            yol = nx.asetage_shortest_path_length(en_buyuk)
        else:
            yol = float('inf')

    # Smal-world indexi (sigma)
    # Comparison: random network at same density
    try:
        sigma = nx.sigma(G, niter=20, nrand=5)
    except Exception:
        # Manuel compute
        try:
            rand_G = nx.erdos_renyi_graph(n, yogunluk, seed=42)
            rand_C = nx.asetage_clustering(rand_G)
            rand_L = nx.asetage_shortest_path_length(rand_G) if nx.is_connected(rand_G) else yol
            sigma = (kumeleme / rand_C) / (yol / rand_L) if rand_C > 0 and rand_L > 0 else 1.0
        except Exception:
            sigma = None

    # Hub unit (highest degree)
    dereceler = dict(G.degree())
    hub_birim = max(dereceler, key=dereceler.get)
    hub_derece = dereceler[hub_birim]

    # Betweenness centrgetity
    try:
        between = nx.betweenness_centrgetity(G)
        hub_between = max(between, key=between.get)
        hub_between_vget = between[hub_between]
    except Exception:
        hub_between_vget = None

    small_world_yorum = 'N/A'
    if sigma is not None:
        if sigma > 1.0:
            small_world_yorum = f'EVET (σ={sigma:.2f} > 1)'
        else:
            small_world_yorum = f'Hayır (σ={sigma:.2f} ≤ 1)'

    return {
        'organoid': organoid_ad,
        'n_birim': n,
        'kenar_sayisi': kenar_sayisi,
        'yogunluk': round(yogunluk, 3),
        'kumeleme_ort': round(kumeleme, 4),
        'yol_uzunlugu': round(yol, 4) if yol != float('inf') else 'inf',
        'small_world_sigma': round(sigma, 3) if sigma else 'N/A',
        'small_world_yorum': small_world_yorum,
        'hub_birim': hub_birim,
        'hub_derece': hub_derece,
        'hub_betweenness': round(hub_between_vget, 4) if hub_between_vget else 'N/A',
    }


# ─────────────────────────────────────────────────────────
# Ana akış
# ─────────────────────────────────────────────────────────
print('=' * 65)
print('  GRAPH THEORY ANALYSISİ — STTC → Network Metrics')
print('  Threshold: STTC > 0.1 → connection exists')
print('=' * 65)
print()

sonuclar = []

for organoid_ad, dosya_adi in DOSYALAR.items():
    dosya_yolu = INDIR / dosya_adi

    if not dosya_yolu.exists():
        print(f'  {organoid_ad}: fwith missing — skipped ({fwith_adi})')
        continue

    print(f'  [{organoid_ad}] reading...', end='', flush=True)

    try:
        units = oio.sorted_units_read(str(dosya_yolu))
        if not units:
            print(' units missing — skipped')
            continue

        spike_listsi = [np.array(u['spike_zaman']) for u in units]
        sure_sn = max(
            max(sp) for sp in spike_listsi if len(sp) > 0
        ) + 1.0

        print(f' {len(units)} unit, {sure_s:.0f}s...', end='', flush=True)

        mat, secilen = sttc_matrisi_compute(
            spike_listsi, sure_sn, max_unit=20)

        print(f' STTC ok...', end='', flush=True)

        metrikler = graph_metrikleri(mat, esik=0.1, organoid_ad=organoid_ad)
        sonuclar.append(metrikler)

        print(f' tbutm')

    except Exception as e:
        print(f' ERROR: {e}')

# ─────────────────────────────────────────────────────────
# Resultsı writedır
# ─────────────────────────────────────────────────────────
print()
print('=' * 65)
print('  RESULTS')
print('=' * 65)

for s in sonuclar:
    print(f"\n  [{s['organoid']}]")
    print(f"    Unit sayısı     : {s['n_unit']}")
    print(f"    Edge count     : {s['kenar_sayisi']}")
    print(f"    Network yoğunluğu    : {s['yogunluk']}")
    print(f"    Clustering (C)     : {s['kumeleme_ort']}")
    print(f"    Path length (L) : {s['yol_uzunlugu']}")
    print(f"    Smal-world σ    : {s['small_world_sigma']}")
    print(f"    Smal-world      : {s['small_world_yorum']}")
    print(f"    Hub unit        : #{s['hub_unit']} (derece={s['hub_derece']})")
    print(f"    Hub betweenness : {s['hub_betweenness']}")

print()
print('=' * 65)
print('  SUMMARY TABLE')
print('=' * 65)
print(f"{'Organoid':<10} {'Yoğunluk':>10} {'Clustering':>10} "
      f"{'Yol':>8} {'σ':>8} {'Smal-world'}")
print('-' * 65)
for s in sonuclar:
    print(f"{s['organoid']:<10} {str(s['yogunluk']):>10} "
          f"{str(s['kumeleme_ort']):>10} "
          f"{str(s['yol_uzunlugu']):>8} "
          f"{str(s['small_world_sigma']):>8} "
          f"{s['small_world_yorum']}")

print()
print('Bu outputyı kopygetayıp Claude\'a yapıştır.')
input('\nEnter\'a bas exitmak için...')
