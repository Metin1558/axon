"""
organoid_lfp.py
================
LFP (Local Field Potential) analiz modülü.
Ham elektriksel sinyal içeren NWB dosyaları için.

NOT: DANDI:001603 spike-sorted dosyaları LFP içermez.
Bu modül ham sinyal (acquisition/ElectricalSeries) olan
NWB dosyaları için tasarlanmıştır. Örneğin DANDI:000034
(Gracias Lab ham kayıtlar) bu modülle kullanılabilir.

Hesaplanan metrikler:
    - Band gücü (delta, theta, alpha, beta, low-gamma)
    - Bant oranları (theta/delta, gamma/delta)
    - Welch PSD (Power Spectral Density)
    - Frekans tepe değeri
    - Dominant osilasyon periyodu
    - Band gücünün zamana göre değişimi (pencereli)

Frekans bantları (nörobiyoloji standardı):
    delta:      0.5 - 4   Hz  (derin uyku, kortikal yavaş dalgalar)
    theta:      4   - 8   Hz  (hafıza, navigasyon)
    alpha:      8   - 12  Hz  (relaks durum)
    beta:       12  - 30  Hz  (aktif düşünme)
    low-gamma:  30  - 80  Hz  (yüksek kognitif işlem)
    high-gamma: 80  - 150 Hz  (yerel ağ aktivitesi)

Kullanım (CLI):
    python organoid_lfp.py <nwb_dosyasi> <organoid> <kayit>
                           [--kanal N] [--out-dir KLASOR]

Kullanım (modül):
    from organoid_lfp import lfp_analiz
    sonuc = lfp_analiz('dosya.nwb', kanal=0)
"""

import sys
import argparse
import csv
from pathlib import Path
import numpy as np
from scipy import signal as scipy_signal


# ─────────────────────────────────────────────────────────────
# Bant tanımları
# ─────────────────────────────────────────────────────────────

BANTLAR = {
    'delta':      (0.5,  4.0),
    'theta':      (4.0,  8.0),
    'alpha':      (8.0, 12.0),
    'beta':       (12.0, 30.0),
    'low_gamma':  (30.0, 80.0),
    'high_gamma': (80.0, 150.0),
}

# ─────────────────────────────────────────────────────────────
# Filtreleme
# ─────────────────────────────────────────────────────────────

def bant_filtre(veri, sr, f_low, f_high, order=4):
    """
    Zero-phase Butterworth bandpass filtresi.
    Zero-phase: sosfiltfilt ile faz gecikmesi sıfır.
    """
    nyq = sr / 2.0
    low = max(f_low / nyq, 1e-6)
    high = min(f_high / nyq, 0.9999)

    if low >= high:
        return np.zeros_like(veri)

    try:
        sos = scipy_signal.butter(order, [low, high], btype='band', output='sos')
        return scipy_signal.sosfiltfilt(sos, veri)
    except Exception:
        return np.zeros_like(veri)


def lfp_filtre(ham_veri, sr, f_max=300.0):
    """
    Ham sinyali LFP aralığına (0.5-300 Hz) düşürür.
    Ham sinyaller genelde 30kHz örnekleme — bu aralık yeterli.
    """
    return bant_filtre(ham_veri, sr, 0.5, min(f_max, sr/2 - 1))


# ─────────────────────────────────────────────────────────────
# Band gücü hesabı
# ─────────────────────────────────────────────────────────────

def band_gucu_hesapla(veri, sr, bantlar=None):
    """
    Welch metodu ile her bant için güç hesaplar.

    Dönen: dict {bant_adı: güç_değeri}
    """
    if bantlar is None:
        bantlar = BANTLAR

    # Welch PSD
    nperseg = min(int(sr * 4), len(veri))  # 4 saniyelik segment
    freqs, psd = scipy_signal.welch(veri, fs=sr, nperseg=nperseg,
                                     scaling='density')

    band_gucler = {}
    for bant_adi, (f_low, f_high) in bantlar.items():
        mask = (freqs >= f_low) & (freqs <= f_high)
        if np.any(mask):
            # Trapez integral ile alan altı = band gücü
            guc = np.trapezoid(psd[mask], freqs[mask])
            band_gucler[bant_adi] = float(guc)
        else:
            band_gucler[bant_adi] = 0.0

    return band_gucler, freqs, psd


def band_oranlari_hesapla(band_gucler):
    """
    Diagnostik bant oranları.
    """
    d = band_gucler.get('delta', 1e-10)
    t = band_gucler.get('theta', 0)
    lg = band_gucler.get('low_gamma', 0)
    b = band_gucler.get('beta', 0)
    toplam = sum(band_gucler.values()) or 1e-10

    return {
        'theta_delta_orani': round(t / d if d > 0 else 0, 4),
        'gamma_delta_orani': round(lg / d if d > 0 else 0, 4),
        'beta_delta_orani': round(b / d if d > 0 else 0, 4),
        'delta_toplam_pct': round(d / toplam * 100, 2),
        'theta_toplam_pct': round(t / toplam * 100, 2),
        'beta_toplam_pct': round(b / toplam * 100, 2),
        'gamma_toplam_pct': round(lg / toplam * 100, 2),
    }


# ─────────────────────────────────────────────────────────────
# Dominant frekans
# ─────────────────────────────────────────────────────────────

def dominant_frekans(freqs, psd, f_min=0.5, f_max=100.0):
    """
    0.5-100 Hz aralığında en yüksek güçlü frekansı bulur.
    """
    mask = (freqs >= f_min) & (freqs <= f_max)
    if not np.any(mask):
        return 0.0, 0.0
    peak_idx = np.argmax(psd[mask])
    dom_freq = freqs[mask][peak_idx]
    dom_guc = psd[mask][peak_idx]
    return float(dom_freq), float(dom_guc)


# ─────────────────────────────────────────────────────────────
# Pencereli band güç zaman serisi
# ─────────────────────────────────────────────────────────────

def pencereli_band_gucu(veri, sr, bant_adi='delta', pencere_sn=10.0):
    """
    Zaman içinde band gücünün nasıl değiştiğini gösterir.
    Bir bandın gücü kayıt boyunca sabit mi yoksa değişiyor mu?

    Dönen: (zaman_dizisi, guc_dizisi)
    """
    f_low, f_high = BANTLAR[bant_adi]
    filtrelenmis = bant_filtre(veri, sr, f_low, f_high)

    pencere_ornek = int(pencere_sn * sr)
    adim = pencere_ornek // 2  # %50 örtüşme
    n_pencere = (len(veri) - pencere_ornek) // adim + 1

    zaman = []
    guc = []
    for i in range(n_pencere):
        bas = i * adim
        bit = bas + pencere_ornek
        if bit > len(veri):
            break
        pencere = filtrelenmis[bas:bit]
        guc.append(float(np.mean(pencere**2)))  # RMS güç
        zaman.append(float((bas + pencere_ornek/2) / sr))

    return np.array(zaman), np.array(guc)


# ─────────────────────────────────────────────────────────────
# Ana analiz fonksiyonu
# ─────────────────────────────────────────────────────────────

def lfp_analiz(nwb_yolu, kanal=0, max_sure_sn=300.0):
    """
    Ham NWB dosyasından LFP analizi yapar.

    nwb_yolu: str — NWB dosyası
    kanal: int — hangi kanal analiz edilecek
    max_sure_sn: float — maksimum analiz süresi (büyük dosyalar için)

    Dönen: dict veya None (LFP verisi yoksa)
    """
    sys.path.insert(0, str(Path(__file__).parent))
    import organoid_io as oio

    # Tip kontrolü
    tip = oio.nwb_tip_tespit(nwb_yolu)
    if tip not in ('ham', 'her_ikisi'):
        return {
            'durum': 'lfp_mevcut_degil',
            'sebep': f'NWB dosyası ham sinyal içermiyor (tip: {tip}). '
                     f'LFP analizi için ham elektriksel sinyal gereklidir.',
        }

    print(f'  Ham sinyal okunuyor (kanal {kanal}, max {max_sure_sn:.0f}sn)...')

    # Chunk okuma — sadece ilk max_sure_sn'yi al
    veri_parcalar = []
    sr = None
    toplam_ornek = 0
    max_ornek = None

    for chunk_data, bas, bit, chunk_sr in oio.ham_chunk_uret(
            nwb_yolu, chunk_sn=60, kanal=kanal):
        if sr is None:
            sr = chunk_sr
            max_ornek = int(max_sure_sn * sr)
        veri_parcalar.append(chunk_data)
        toplam_ornek += len(chunk_data)
        if toplam_ornek >= max_ornek:
            break

    if not veri_parcalar:
        return {'durum': 'hata', 'sebep': 'Veri okunamadı'}

    veri = np.concatenate(veri_parcalar)[:max_ornek]
    sure_sn = len(veri) / sr

    print(f'  {sure_sn:.1f}sn, SR={sr}Hz, {len(veri)} örnek')

    # uV dönüşümü (gerekirse)
    if np.std(veri) < 1e-3:
        veri = veri * 1e6

    # LFP filtresi (0.5-300 Hz)
    veri_lfp = lfp_filtre(veri, sr, f_max=min(300.0, sr/2 - 1))

    # Band güçleri
    print('  Band güçleri hesaplanıyor...')
    band_gucler, freqs, psd = band_gucu_hesapla(veri_lfp, sr)

    # Band oranları
    band_oranlari = band_oranlari_hesapla(band_gucler)

    # Dominant frekans
    dom_freq, dom_guc = dominant_frekans(freqs, psd)

    # En güçlü bant
    en_guclu_bant = max(band_gucler, key=lambda k: band_gucler[k])

    # Pencereli delta güç (zamana göre değişim)
    zaman_d, guc_d = pencereli_band_gucu(veri_lfp, sr, 'delta', pencere_sn=10.0)
    zaman_t, guc_t = pencereli_band_gucu(veri_lfp, sr, 'theta', pencere_sn=10.0)

    return {
        'durum': 'tamam',
        'kanal': kanal,
        'sr': sr,
        'sure_sn': round(sure_sn, 2),
        'dominant_frekans_hz': round(dom_freq, 3),
        'en_guclu_bant': en_guclu_bant,
        'band_gucler': {k: round(v, 6) for k, v in band_gucler.items()},
        'band_oranlari': band_oranlari,
        'freqs': freqs.tolist(),
        'psd': psd.tolist(),
        'zaman_delta': zaman_d.tolist(),
        'guc_delta': guc_d.tolist(),
        'zaman_theta': zaman_t.tolist(),
        'guc_theta': guc_t.tolist(),
    }


# ─────────────────────────────────────────────────────────────
# Grafik
# ─────────────────────────────────────────────────────────────

def lfp_grafik_uret(sonuc, organoid, kayit, cikti_yolu):
    """LFP analiz sonuçlarını görselleştirir."""
    if sonuc.get('durum') != 'tamam':
        return

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), facecolor='white')
    fig.suptitle(f'{organoid}/{kayit} — LFP Analizi (Kanal {sonuc["kanal"]})',
                 fontsize=12, fontweight='bold')

    band_renkler = {
        'delta': '#1F4E79',
        'theta': '#2E75B6',
        'alpha': '#5BA3D9',
        'beta': '#E67E22',
        'low_gamma': '#E74C3C',
        'high_gamma': '#8E44AD',
    }

    # Panel 1: Welch PSD
    ax = axes[0, 0]
    freqs = np.array(sonuc['freqs'])
    psd = np.array(sonuc['psd'])
    mask = (freqs >= 0.5) & (freqs <= 150)
    ax.semilogy(freqs[mask], psd[mask], color='black', lw=1.2)

    for bant_adi, (f_low, f_high) in BANTLAR.items():
        bant_mask = (freqs >= f_low) & (freqs <= f_high) & mask
        if np.any(bant_mask):
            ax.fill_between(freqs[bant_mask], psd[bant_mask],
                             alpha=0.4, color=band_renkler.get(bant_adi, 'gray'),
                             label=bant_adi)

    dom = sonuc['dominant_frekans_hz']
    ax.axvline(dom, color='red', ls='--', lw=1.2,
               label=f'Dominant: {dom:.2f} Hz')
    ax.set_xlabel('Frekans (Hz)')
    ax.set_ylabel('PSD (µV²/Hz)')
    ax.set_title('Welch Güç Spektral Yoğunluğu')
    ax.legend(fontsize=7, ncol=2)
    ax.grid(alpha=0.3)

    # Panel 2: Band güçleri çubuk
    ax = axes[0, 1]
    bantlar_sira = ['delta', 'theta', 'alpha', 'beta', 'low_gamma', 'high_gamma']
    gucler = [sonuc['band_gucler'].get(b, 0) for b in bantlar_sira]
    renkler = [band_renkler.get(b, 'gray') for b in bantlar_sira]
    bars = ax.bar(bantlar_sira, gucler, color=renkler, edgecolor='black', alpha=0.8)
    ax.set_yscale('log')
    ax.set_xlabel('Frekans bandı')
    ax.set_ylabel('Band gücü (µV²)')
    ax.set_title(f'Band Güçleri\nEn güçlü: {sonuc["en_guclu_bant"]}')
    ax.tick_params(axis='x', rotation=30)
    ax.grid(alpha=0.3, axis='y')

    # Panel 3: Delta güç zaman serisi
    ax = axes[1, 0]
    zaman_d = np.array(sonuc['zaman_delta'])
    guc_d = np.array(sonuc['guc_delta'])
    if len(zaman_d) > 0:
        ax.plot(zaman_d, guc_d, color=band_renkler['delta'], lw=1.2)
        ax.fill_between(zaman_d, 0, guc_d, alpha=0.3,
                         color=band_renkler['delta'])
    ax.set_xlabel('Zaman (sn)')
    ax.set_ylabel('Delta güç (µV² RMS)')
    ax.set_title('Delta Bandı (0.5-4 Hz) Zaman Serisi')
    ax.grid(alpha=0.3)

    # Panel 4: Band oranları
    ax = axes[1, 1]
    oranlar = sonuc['band_oranlari']
    oran_adlari = ['θ/δ', 'γ/δ', 'β/δ', 'δ%', 'θ%', 'β%', 'γ%']
    oran_degerler = [
        oranlar.get('theta_delta_orani', 0),
        oranlar.get('gamma_delta_orani', 0),
        oranlar.get('beta_delta_orani', 0),
        oranlar.get('delta_toplam_pct', 0),
        oranlar.get('theta_toplam_pct', 0),
        oranlar.get('beta_toplam_pct', 0),
        oranlar.get('gamma_toplam_pct', 0),
    ]
    bars = ax.barh(oran_adlari, oran_degerler,
                   color=['steelblue']*3 + ['darkred']*4,
                   edgecolor='black', alpha=0.8)
    ax.set_xlabel('Değer')
    ax.set_title('Diagnostik Band Oranları ve Yüzdeler')
    ax.grid(alpha=0.3, axis='x')
    for bar, val in zip(bars, oran_degerler):
        ax.text(bar.get_width() + 0.01 * max(oran_degerler),
                bar.get_y() + bar.get_height()/2,
                f'{val:.3f}', va='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(cikti_yolu, dpi=130, bbox_inches='tight', facecolor='white')
    plt.close()


# ─────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='LFP analizi — ham sinyal içeren NWB dosyaları için')
    parser.add_argument('dosya', help='NWB dosyası (ham sinyal içermeli)')
    parser.add_argument('organoid', help='Organoid kodu')
    parser.add_argument('kayit', help='Kayıt kodu')
    parser.add_argument('--kanal', type=int, default=0,
                        help='Kanal numarası (varsayılan: 0)')
    parser.add_argument('--max-sure', type=float, default=300.0,
                        help='Maksimum analiz süresi sn (varsayılan: 300)')
    parser.add_argument('--out-dir', default='sonuclar/lfp')
    args = parser.parse_args()

    cikti = Path(args.out_dir)
    cikti.mkdir(parents=True, exist_ok=True)

    print(f'[LFP Analizi] {args.organoid}/{args.kayit}')
    print(f'  Dosya: {args.dosya}')
    print(f'  Kanal: {args.kanal}')

    sonuc = lfp_analiz(args.dosya, kanal=args.kanal,
                        max_sure_sn=args.max_sure)

    if sonuc.get('durum') == 'lfp_mevcut_degil':
        print(f'\n  [BİLGİ] {sonuc["sebep"]}')
        print('  LFP analizi bu dosya için yapılamaz.')
        print('  Ham sinyal içeren dosyalar için kullanın.')
        sys.exit(0)

    if sonuc.get('durum') != 'tamam':
        print(f'  HATA: {sonuc.get("sebep", "Bilinmeyen hata")}')
        sys.exit(1)

    # Grafik
    grafik_yolu = cikti / f'{args.organoid}_{args.kayit}_lfp.png'
    lfp_grafik_uret(sonuc, args.organoid, args.kayit, str(grafik_yolu))

    # CSV (ağır verileri çıkar)
    csv_sonuc = {
        'organoid': args.organoid,
        'kayit': args.kayit,
        'kanal': sonuc['kanal'],
        'sr': sonuc['sr'],
        'sure_sn': sonuc['sure_sn'],
        'dominant_frekans_hz': sonuc['dominant_frekans_hz'],
        'en_guclu_bant': sonuc['en_guclu_bant'],
    }
    csv_sonuc.update({f'band_{k}': v
                       for k, v in sonuc['band_gucler'].items()})
    csv_sonuc.update(sonuc['band_oranlari'])

    csv_yolu = cikti / f'{args.organoid}_{args.kayit}_lfp.csv'
    with open(csv_yolu, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(csv_sonuc.keys()))
        w.writeheader()
        w.writerow(csv_sonuc)

    # Konsol özeti
    print()
    print('  ─── LFP SONUÇLARI ───')
    print(f'  Dominant frekans : {sonuc["dominant_frekans_hz"]:.3f} Hz')
    print(f'  En güçlü bant    : {sonuc["en_guclu_bant"]}')
    print()
    print('  Band güçleri:')
    for bant, guc in sonuc['band_gucler'].items():
        print(f'    {bant:<15}: {guc:.6f} µV²')
    print()
    print('  Band oranları:')
    for oran, val in sonuc['band_oranlari'].items():
        print(f'    {oran:<25}: {val}')
    print()
    print(f'  Grafik: {grafik_yolu}')
    print(f'  CSV: {csv_yolu}')


if __name__ == '__main__':
    main()
