"""
graph_analiz.py
===============
STTC matrislerinden graph theory metrikleri hesaplar.
Sonuçları terminale yazar — dosya atmana gerek yok,
çıktıyı kopyalayıp Claude'a yapıştır.

Kullanım:
    python graph_analiz.py

Gereksinimler:
    pip install networkx numpy scipy
"""

import sys
import os
import subprocess

# networkx kur (yoksa)
try:
    import networkx as nx
    print("networkx OK")
except ImportError:
    print("networkx kuruluyor...")
    subprocess.run([sys.executable, '-m', 'pip', 'install',
                    '--quiet', 'networkx'], check=True)
    import networkx as nx
    print("networkx kuruldu")

import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'v6_3'))
import organoid_io as oio
import organoid_metrics as omet

# ─────────────────────────────────────────────────────────
# Hangi NWB dosyaları analiz edilecek
# ─────────────────────────────────────────────────────────
INDIR = Path('dandi_indirilmis')

# Her organoidin en büyük (en kapsamlı) dosyasını seç
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
# Graph theory fonksiyonları
# ─────────────────────────────────────────────────────────

def sttc_matrisi_hesapla(spike_listesi, sure_sn, dt=0.01, max_unit=20):
    """En aktif N birimin STTC matrisini hesaplar."""
    spike_sayisi = [len(sp) for sp in spike_listesi]
    siralama = np.argsort(spike_sayisi)[::-1]
    secilen = siralama[:max_unit]
    secilen_treni = [spike_listesi[i] for i in secilen]

    n = len(secilen_treni)
    matris = np.zeros((n, n))

    for i in range(n):
        for j in range(i+1, n):
            try:
                sttc = omet.sttc_iki_kanal(
                    secilen_treni[i], secilen_treni[j],
                    sure_sn, dt=dt)
                matris[i, j] = sttc
                matris[j, i] = sttc
            except Exception:
                pass
        matris[i, i] = 1.0  # diagonal

    return matris, secilen


def graph_metrikleri(sttc_mat, esik=0.1, organoid_ad='?'):
    """STTC matrisinden graph theory metrikleri hesaplar."""
    n = len(sttc_mat)

    # Eşik uygula: STTC > esik ise bağlantı var
    adj = (sttc_mat > esik).astype(float)
    np.fill_diagonal(adj, 0)

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
            'kümeleme_ort': 0.0,
            'yol_uzunlugu': float('inf'),
            'small_world': 'N/A (bağlantısız)',
            'hub_birim': 'yok',
            'hub_derece': 0,
        }

    # Kümeleme katsayısı
    kumeleme = nx.average_clustering(G)

    # Ortalama yol uzunluğu (bağlı bileşenler üzerinden)
    if nx.is_connected(G):
        yol = nx.average_shortest_path_length(G)
    else:
        # En büyük bağlı bileşen
        en_buyuk = G.subgraph(max(nx.connected_components(G), key=len))
        if len(en_buyuk) > 1:
            yol = nx.average_shortest_path_length(en_buyuk)
        else:
            yol = float('inf')

    # Small-world indeksi (sigma)
    # Karşılaştırma: aynı yoğunlukta rastgele ağ
    try:
        sigma = nx.sigma(G, niter=20, nrand=5)
    except Exception:
        # Manuel hesapla
        try:
            rand_G = nx.erdos_renyi_graph(n, yogunluk, seed=42)
            rand_C = nx.average_clustering(rand_G)
            rand_L = nx.average_shortest_path_length(rand_G) if nx.is_connected(rand_G) else yol
            sigma = (kumeleme / rand_C) / (yol / rand_L) if rand_C > 0 and rand_L > 0 else 1.0
        except Exception:
            sigma = None

    # Hub birim (en yüksek derece)
    dereceler = dict(G.degree())
    hub_birim = max(dereceler, key=dereceler.get)
    hub_derece = dereceler[hub_birim]

    # Betweenness centrality
    try:
        between = nx.betweenness_centrality(G)
        hub_between = max(between, key=between.get)
        hub_between_val = between[hub_between]
    except Exception:
        hub_between_val = None

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
        'hub_betweenness': round(hub_between_val, 4) if hub_between_val else 'N/A',
    }


# ─────────────────────────────────────────────────────────
# Ana akış
# ─────────────────────────────────────────────────────────
print('=' * 65)
print('  GRAPH THEORY ANALİZİ — STTC → Ağ Metrikleri')
print('  Eşik: STTC > 0.1 → bağlantı var')
print('=' * 65)
print()

sonuclar = []

for organoid_ad, dosya_adi in DOSYALAR.items():
    dosya_yolu = INDIR / dosya_adi

    if not dosya_yolu.exists():
        print(f'  {organoid_ad}: dosya yok — atlandı ({dosya_adi})')
        continue

    print(f'  [{organoid_ad}] okunuyor...', end='', flush=True)

    try:
        units = oio.sorted_units_oku(str(dosya_yolu))
        if not units:
            print(' units yok — atlandı')
            continue

        spike_listesi = [np.array(u['spike_zaman']) for u in units]
        sure_sn = max(
            max(sp) for sp in spike_listesi if len(sp) > 0
        ) + 1.0

        print(f' {len(units)} unit, {sure_sn:.0f}sn...', end='', flush=True)

        mat, secilen = sttc_matrisi_hesapla(
            spike_listesi, sure_sn, max_unit=20)

        print(f' STTC ok...', end='', flush=True)

        metrikler = graph_metrikleri(mat, esik=0.1, organoid_ad=organoid_ad)
        sonuclar.append(metrikler)

        print(f' tamam')

    except Exception as e:
        print(f' HATA: {e}')

# ─────────────────────────────────────────────────────────
# Sonuçları yazdır
# ─────────────────────────────────────────────────────────
print()
print('=' * 65)
print('  SONUÇLAR')
print('=' * 65)

for s in sonuclar:
    print(f"\n  [{s['organoid']}]")
    print(f"    Birim sayısı     : {s['n_birim']}")
    print(f"    Kenar sayısı     : {s['kenar_sayisi']}")
    print(f"    Ağ yoğunluğu    : {s['yogunluk']}")
    print(f"    Kümeleme (C)     : {s['kumeleme_ort']}")
    print(f"    Yol uzunluğu (L) : {s['yol_uzunlugu']}")
    print(f"    Small-world σ    : {s['small_world_sigma']}")
    print(f"    Small-world      : {s['small_world_yorum']}")
    print(f"    Hub birim        : #{s['hub_birim']} (derece={s['hub_derece']})")
    print(f"    Hub betweenness  : {s['hub_betweenness']}")

print()
print('=' * 65)
print('  ÖZET TABLO')
print('=' * 65)
print(f"{'Organoid':<10} {'Yoğunluk':>10} {'Kümeleme':>10} "
      f"{'Yol':>8} {'σ':>8} {'Small-world'}")
print('-' * 65)
for s in sonuclar:
    print(f"{s['organoid']:<10} {str(s['yogunluk']):>10} "
          f"{str(s['kumeleme_ort']):>10} "
          f"{str(s['yol_uzunlugu']):>8} "
          f"{str(s['small_world_sigma']):>8} "
          f"{s['small_world_yorum']}")

print()
print('Bu çıktıyı kopyalayıp Claude\'a yapıştır.')
input('\nEnter\'a bas çıkmak için...')
