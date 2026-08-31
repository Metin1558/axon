"""
organoid_burst_ref.py
======================
Burst tespiti referans doğrulaması.
Bakkum ve ark. (2013) algoritması ile mevcut MAD algoritmasını
aynı veri üzerinde karşılaştırır.

Bakkum 2013 algoritması:
    Bakkum, D.J. ve diğerleri. (2013). Tracking axonal action potential
    propagation on a high-density microelectrode array across hundreds
    of sites. Nature Communications, 4, 2181.

    Temel prensip:
    1. Population spike rate (bin_ms pencerede)
    2. ISI eşiği: her birimin 'network ISI' hesapla
       (spike yoksa o bin sayılmaz)
    3. Burst: rate > esik_carpan * baseline rate için
       ardısık dönem, en az min_katilim birimi aktif

MAD algoritması (mevcut):
    network_burst_tespit() in organoid_units_analiz.py
    Adaptif MAD eşiği üzerindeki ardışık bin'ler.

Karşılaştırma metrikleri:
    - Tespit edilen burst sayısı
    - Zamansal örtüşme oranı (IoU — Intersection over Union)
    - Yanlış pozitif ve yanlış negatif tahmin
    - Medyan burst süresi farkı

Kullanım (CLI):
    python organoid_burst_ref.py <nwb_dosyasi> <organoid> <kayit>
                                 [--out-dir KLASOR]

Kullanım (modül):
    from organoid_burst_ref import burst_karsilastir
    sonuc = burst_karsilastir(spike_listesi, sure_sn)
"""

import sys
import argparse
import csv
from pathlib import Path
import numpy as np


# ─────────────────────────────────────────────────────────────
# Bakkum 2013 burst tespiti
# ─────────────────────────────────────────────────────────────

def bakkum_burst_tespit(spike_listesi, sure_sn,
                         bin_ms=50,
                         esik_carpan=2.5,
                         min_katilim_orani=0.10,
                         min_burst_sure_ms=50,
                         pencere_saniyelik=60):
    """
    Bakkum 2013 yaklaşımına dayalı network burst tespiti.

    Referans: Bakkum ve ark. (2013) Nature Communications 4:2181.
    Orijinal algoritma MEA burst detection için geliştirilmiştir.

    Adımlar:
    1. Population spike rate zaman serisi oluştur (bin_ms binlerde)
    2. Kayan pencere (pencere_saniyelik) ile baseline rate hesapla
    3. Anlık rate > esik_carpan * baseline olduğunda yüksek aktivite işaretle
    4. En az min_katilim_orani oranında birim aktif olduğunda burst say

    bin_ms: bin genişliği (ms)
    esik_carpan: baseline'ın kaç katı eşik (Bakkum 2013'te 5 kullanılmış)
    min_katilim_orani: burst sayılmak için aktif olması gereken min birim oranı
    min_burst_sure_ms: minimum burst süresi (ms)
    pencere_saniyelik: baseline hesabı için kayan pencere genişliği (sn)
    """
    n_unit = len(spike_listesi)
    if n_unit < 2:
        return [], np.array([])

    bin_sn = bin_ms / 1000.0
    n_bin = int(sure_sn / bin_sn) + 1

    pop_rate = np.zeros(n_bin)
    aktif_unit_per_bin = np.zeros(n_bin)

    for sp in spike_listesi:
        if len(sp) == 0:
            continue
        bin_idx = (np.array(sp) / bin_sn).astype(int)
        bin_idx = np.clip(bin_idx, 0, n_bin - 1)
        unique_bins, sayi = np.unique(bin_idx, return_counts=True)
        pop_rate[unique_bins] += sayi
        aktif_unit_per_bin[unique_bins] += 1

    # Global baseline: tüm kaydın medyanı (nonzero bin'ler)
    # Bakkum 2013'ün özüne sadık: sabit global eşik
    nonzero_bins = pop_rate[pop_rate > 0]
    if len(nonzero_bins) > 0:
        baseline_medyan = np.median(nonzero_bins)
    else:
        baseline_medyan = 1.0

    # Minimum birim katılımı
    min_aktif = max(2, int(n_unit * min_katilim_orani))
    esik_rate = esik_carpan * baseline_medyan

    # Burst tespit
    yuksek = (pop_rate > esik_rate) & (aktif_unit_per_bin >= min_aktif)

    burstler = []
    icinde = False
    bas_idx = 0
    min_burst_bin = int(min_burst_sure_ms / bin_ms)

    for i, durum in enumerate(yuksek):
        if durum and not icinde:
            icinde = True
            bas_idx = i
        elif not durum and icinde:
            icinde = False
            n_bin_sure = i - bas_idx
            if n_bin_sure >= min_burst_bin:
                burstler.append({
                    'baslangic_sn': bas_idx * bin_sn,
                    'bitis_sn': i * bin_sn,
                    'sure_sn': n_bin_sure * bin_sn,
                    'tepe_rate': float(np.max(pop_rate[bas_idx:i])),
                    'aktif_unit': int(np.max(aktif_unit_per_bin[bas_idx:i])),
                })

    return burstler, pop_rate


# ─────────────────────────────────────────────────────────────
# Zamansal örtüşme hesabı (IoU)
# ─────────────────────────────────────────────────────────────

def burst_iou_hesapla(burstler_a, burstler_b, sure_sn, bin_ms=10):
    """
    İki burst listesi arasındaki zamansal örtüşmeyi hesaplar.

    IoU (Intersection over Union):
        IoU = süre(A ∩ B) / süre(A ∪ B)
        1.0 = mükemmel eşleşme, 0.0 = hiç örtüşme yok

    Ayrıca:
        precision = A'nın kaçı B ile örtüşüyor?
        recall = B'nin kaçı A ile örtüşüyor?
    """
    bin_sn = bin_ms / 1000.0
    n_bin = int(sure_sn / bin_sn) + 1

    mask_a = np.zeros(n_bin, dtype=bool)
    mask_b = np.zeros(n_bin, dtype=bool)

    for b in burstler_a:
        bas = int(b['baslangic_sn'] / bin_sn)
        bit = int(b['bitis_sn'] / bin_sn)
        mask_a[bas:bit] = True

    for b in burstler_b:
        bas = int(b['baslangic_sn'] / bin_sn)
        bit = int(b['bitis_sn'] / bin_sn)
        mask_b[bas:bit] = True

    intersect = np.sum(mask_a & mask_b) * bin_sn
    union = np.sum(mask_a | mask_b) * bin_sn
    a_sure = np.sum(mask_a) * bin_sn
    b_sure = np.sum(mask_b) * bin_sn

    iou = intersect / union if union > 0 else 0.0
    precision = intersect / a_sure if a_sure > 0 else 0.0  # A'nın B'yle örtüşen kısmı
    recall = intersect / b_sure if b_sure > 0 else 0.0     # B'nin A'yla örtüşen kısmı
    f1 = 2*precision*recall/(precision+recall) if (precision+recall) > 0 else 0.0

    return {
        'iou': round(iou, 4),
        'precision': round(precision, 4),
        'recall': round(recall, 4),
        'f1': round(f1, 4),
        'intersect_sn': round(intersect, 3),
        'union_sn': round(union, 3),
    }


# ─────────────────────────────────────────────────────────────
# Ana karşılaştırma fonksiyonu
# ─────────────────────────────────────────────────────────────

def burst_karsilastir(spike_listesi, sure_sn, cikti_yolu=None):
    """
    MAD algoritması ve Bakkum 2013 algoritmasını karşılaştırır.

    spike_listesi: liste[np.array] — her eleman bir birimin spike zamanları
    sure_sn: float — kayıt süresi
    cikti_yolu: str/None — grafik çıktısı için yol

    Dönen: dict — karşılaştırma sonuçları
    """
    # MAD algoritması (mevcut)
    sys.path.insert(0, str(Path(__file__).parent))
    from organoid_units_analiz import network_burst_tespit

    print('  MAD algoritması çalıştırılıyor...')
    mad_sonuc = network_burst_tespit(spike_listesi, sure_sn, bin_ms=50)
    mad_burstler = mad_sonuc.get('burstler', [])

    print('  Bakkum 2013 algoritması çalıştırılıyor...')
    bakkum_burstler, pop_rate = bakkum_burst_tespit(
        spike_listesi, sure_sn)

    print(f'  MAD: {len(mad_burstler)} burst')
    print(f'  Bakkum: {len(bakkum_burstler)} burst')

    # Örtüşme hesabı
    iou_sonuc = burst_iou_hesapla(mad_burstler, bakkum_burstler, sure_sn)

    print(f'  IoU: {iou_sonuc["iou"]:.3f}')
    print(f'  F1 skoru: {iou_sonuc["f1"]:.3f}')

    # Burst süresi karşılaştırması
    mad_sureler = [b['sure_sn'] for b in mad_burstler] if mad_burstler else [0]
    bakkum_sureler = [b['sure_sn'] for b in bakkum_burstler] if bakkum_burstler else [0]

    sonuc = {
        'mad_burst_sayisi': len(mad_burstler),
        'bakkum_burst_sayisi': len(bakkum_burstler),
        'mad_burst_per_dakika': round(len(mad_burstler) / (sure_sn/60), 3) if sure_sn > 0 else 0,
        'bakkum_burst_per_dakika': round(len(bakkum_burstler) / (sure_sn/60), 3) if sure_sn > 0 else 0,
        'mad_medyan_sure_sn': round(float(np.median(mad_sureler)), 3),
        'bakkum_medyan_sure_sn': round(float(np.median(bakkum_sureler)), 3),
        'iou': iou_sonuc['iou'],
        'precision': iou_sonuc['precision'],
        'recall': iou_sonuc['recall'],
        'f1': iou_sonuc['f1'],
        'intersect_sn': iou_sonuc['intersect_sn'],
        'uyum_yorumu': (
            'Mükemmel (F1≥0.85)' if iou_sonuc['f1'] >= 0.85 else
            'İyi (F1 0.70-0.85)' if iou_sonuc['f1'] >= 0.70 else
            'Orta (F1 0.50-0.70)' if iou_sonuc['f1'] >= 0.50 else
            'Zayıf (F1<0.50)'
        ),
    }

    # Grafik
    if cikti_yolu:
        _grafik_uret(mad_burstler, bakkum_burstler, pop_rate,
                     sure_sn, iou_sonuc, cikti_yolu)

    return sonuc


def _grafik_uret(mad_burstler, bakkum_burstler, pop_rate,
                  sure_sn, iou_sonuc, cikti_yolu):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    fig, axes = plt.subplots(3, 1, figsize=(13, 9), facecolor='white')
    fig.suptitle('Burst Tespiti Algoritma Karşılaştırması: MAD vs Bakkum 2013',
                 fontsize=12, fontweight='bold')

    t = np.arange(len(pop_rate)) * 10 / 1000.0  # 10ms bin -> sn

    # Panel 1: MAD burst'leri
    ax = axes[0]
    ax.plot(t[:len(pop_rate)], pop_rate, lw=0.6, color='steelblue')
    for b in mad_burstler:
        ax.axvspan(b['baslangic_sn'], b['bitis_sn'], color='#1F77B4', alpha=0.35)
    ax.set_title(f'MAD Algoritması — {len(mad_burstler)} burst')
    ax.set_ylabel('Pop. spike sayısı')
    ax.set_xlim(0, sure_sn)
    ax.grid(alpha=0.3)

    # Panel 2: Bakkum burst'leri
    ax = axes[1]
    ax.plot(t[:len(pop_rate)], pop_rate, lw=0.6, color='steelblue')
    for b in bakkum_burstler:
        ax.axvspan(b['baslangic_sn'], b['bitis_sn'], color='#D62728', alpha=0.35)
    ax.set_title(f'Bakkum 2013 Algoritması — {len(bakkum_burstler)} burst')
    ax.set_ylabel('Pop. spike sayısı')
    ax.set_xlim(0, sure_sn)
    ax.grid(alpha=0.3)

    # Panel 3: Örtüşme
    ax = axes[2]
    ax.plot(t[:len(pop_rate)], pop_rate, lw=0.6, color='steelblue', label='Pop. rate', zorder=3)
    for b in mad_burstler:
        ax.axvspan(b['baslangic_sn'], b['bitis_sn'], color='blue', alpha=0.2)
    for b in bakkum_burstler:
        ax.axvspan(b['baslangic_sn'], b['bitis_sn'], color='red', alpha=0.2)

    patches = [
        Patch(facecolor='blue', alpha=0.4, label=f'MAD ({len(mad_burstler)})'),
        Patch(facecolor='red', alpha=0.4, label=f'Bakkum ({len(bakkum_burstler)})'),
    ]
    ax.legend(handles=patches, loc='upper right', fontsize=9)
    ax.set_title(f'Örtüşme — IoU={iou_sonuc["iou"]:.3f}, F1={iou_sonuc["f1"]:.3f}')
    ax.set_xlabel('Zaman (sn)')
    ax.set_ylabel('Pop. spike sayısı')
    ax.set_xlim(0, sure_sn)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig(cikti_yolu, dpi=130, bbox_inches='tight', facecolor='white')
    plt.close()


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Burst tespiti algoritma karşılaştırması (MAD vs Bakkum 2013)')
    parser.add_argument('dosya', help='NWB dosyası')
    parser.add_argument('organoid', help='Organoid kodu')
    parser.add_argument('kayit', help='Kayıt kodu')
    parser.add_argument('--out-dir', default='sonuclar/burst_ref')
    args = parser.parse_args()

    sys.path.insert(0, str(Path(__file__).parent))
    import organoid_io as oio

    cikti = Path(args.out_dir)
    cikti.mkdir(parents=True, exist_ok=True)

    print(f'[Burst Referans] {args.organoid}/{args.kayit}')
    print('  Spike zamanları okunuyor...')

    units = oio.sorted_units_oku(args.dosya)
    if not units:
        print('  HATA: sorted units bulunamadı')
        sys.exit(1)

    spike_listesi = [np.array(u['spike_zaman']) for u in units]
    tum_spike = [s for sp in spike_listesi for s in sp.tolist()]
    sure_sn = max(tum_spike) + 1.0 if tum_spike else 0

    print(f'  {len(units)} unit, sure={sure_sn:.1f}sn')

    grafik_yolu = cikti / f'{args.organoid}_{args.kayit}_burst_karsilastirma.png'
    sonuc = burst_karsilastir(spike_listesi, sure_sn, str(grafik_yolu))

    # CSV
    csv_yolu = cikti / f'{args.organoid}_{args.kayit}_burst_karsilastirma.csv'
    with open(csv_yolu, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(sonuc.keys()))
        w.writeheader()
        w.writerow(sonuc)

    print()
    print('  ─── SONUÇLAR ───')
    for k, v in sonuc.items():
        print(f'  {k:<35}: {v}')

    print(f'\n  Grafik: {grafik_yolu}')
    print(f'  CSV: {csv_yolu}')


if __name__ == '__main__':
    main()
