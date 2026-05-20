"""
organoid_units_analiz.py
=========================
Sorted units (NWB units tablosu) icin per-unit analiz.

ESKI organoid.py problemi: Tum unit'lerin spike zamanlarini birlestirip
tek liste olarak analiz ediyordu. Sonuc: refractory ihlal=%52, CV=8.7
gibi ANLAMSIZ degerler. Cunku farkli noronlardan gelen spike'lar
"ayni noron" gibi muamele gordu.

YENI cozum: Her unit AYRI analiz edilir. Sonuclar:
- Her unit icin ayri satir CSV
- Ozet istatistikler (median Hz, median CV, kac aktif unit vs)
- Cross-unit STTC senkron matrisi
- Network-level burst tespiti
- Tum unit'lerin raster + populasyon Hz grafigi

Kullanim:
    python3 organoid_units_analiz.py <dosya.nwb> <organoid> <kayit>
                                     [--out-dir KLASOR]
                                     [--min-spike N]
                                     [--top-units N]
"""

import sys
import os
import argparse
import time
import csv
from pathlib import Path
from datetime import datetime

import numpy as np

# v6.3 modüllerini yükle
sys.path.insert(0, str(Path(__file__).parent))
import organoid_io as oio
import organoid_qc as oqc
import organoid_metrics as omet


# ============================================================
# UNIT BASINA ANALIZ
# ============================================================

def unit_analiz(unit_id, spike_zaman, sure_sn, refractory_ms=1.0,
                pencere_sn=30, perm_n=200, perm_seed=42):
    """
    Tek bir unit (noron) icin tum metrikleri hesaplar.
    
    Donen: dict (CSV satiri olur)
    """
    n_spike = len(spike_zaman)
    
    # Unit cok az spike'a sahipse atla
    if n_spike < 10:
        return {
            'unit_id': unit_id,
            'n_spike': n_spike,
            'durum': 'cok_az_spike',
            # Diger alanlar bos
        }
    
    # Bos hesap icin sirali ve unique olsun
    spike_zaman = np.sort(np.unique(spike_zaman))
    
    # ISI metrikleri
    isi = omet.isi_metrikleri(spike_zaman)
    
    # Refractory ihlal (Wilson CI)
    refr = oqc.refractory_ihlal_orani(spike_zaman, refractory_ms=refractory_ms)
    
    # Negatif ISI (sirali kontrol)
    neg_isi = oqc.negatif_isi_sayisi(spike_zaman)
    
    # Aktiflik (kac pencerede aktif) - tuple doner
    aktiflik_orani, aktif_pencere, toplam_pencere = omet.aktiflik_pencere(
        spike_zaman, sure_sn, pencere_sn=pencere_sn)
    
    # Ortalama Hz
    ort_hz = n_spike / sure_sn if sure_sn > 0 else 0
    
    # Permutasyon - cok pahali, sadece yeterince spike varsa
    perm_p = None
    if n_spike >= 50 and aktif_pencere >= 3:
        try:
            perm = omet.isi_shuffle_permutasyon(
                spike_zaman, sure_sn, pencere_sn=pencere_sn,
                n_permutasyon=perm_n, seed=perm_seed)
            if perm:
                perm_p = perm['p_deger']
        except Exception:
            pass
    
    # Sonuc dict
    return {
        'unit_id': unit_id,
        'n_spike': n_spike,
        'durum': 'ok',
        'ort_hz': round(ort_hz, 4),
        'ort_isi_sn': round(isi['ort_isi'], 6) if isi.get('ort_isi') is not None else None,
        'medyan_isi_sn': round(isi['medyan_isi'], 6) if isi.get('medyan_isi') is not None else None,
        'cv': round(isi['cv'], 4) if isi.get('cv') is not None else None,
        'isi_q25': round(isi['isi_q25'], 6) if isi.get('isi_q25') is not None else None,
        'isi_q75': round(isi['isi_q75'], 6) if isi.get('isi_q75') is not None else None,
        'refractory_ihlal_orani': round(refr['ihlal_orani'], 4),
        'refractory_ihlal_sayisi': refr['ihlal_sayisi'],
        'refractory_ci_alt': round(refr['ci_alt'], 4),
        'refractory_ci_ust': round(refr['ci_ust'], 4),
        'negatif_isi': neg_isi,
        'aktif_pencere': int(aktif_pencere),
        'toplam_pencere': int(toplam_pencere),
        'aktiflik_orani': round(aktiflik_orani, 4),
        'perm_p_deger': round(perm_p, 4) if perm_p is not None else None,
    }


def populasyon_metrikleri(unit_listesi, sure_sn, pencere_sn=1.0):
    """
    Tum unit'lerden cikarilacak ag-duzeyi metrikler.
    """
    aktif_unitler = [u for u in unit_listesi if u['durum'] == 'ok']
    n_aktif = len(aktif_unitler)
    n_toplam = len(unit_listesi)
    
    metrikler = {
        'toplam_unit': n_toplam,
        'aktif_unit': n_aktif,
        'aktif_oran': round(n_aktif / n_toplam, 4) if n_toplam > 0 else 0,
    }
    
    if n_aktif == 0:
        return metrikler
    
    # Hz dagilimi
    hz_listesi = [u['ort_hz'] for u in aktif_unitler]
    metrikler['hz_medyan'] = round(np.median(hz_listesi), 4)
    metrikler['hz_ort'] = round(np.mean(hz_listesi), 4)
    metrikler['hz_iqr'] = round(np.percentile(hz_listesi, 75) -
                                  np.percentile(hz_listesi, 25), 4)
    metrikler['hz_min'] = round(min(hz_listesi), 4)
    metrikler['hz_max'] = round(max(hz_listesi), 4)
    
    # CV dagilimi
    cv_listesi = [u['cv'] for u in aktif_unitler if u.get('cv') is not None]
    if cv_listesi:
        metrikler['cv_medyan'] = round(np.median(cv_listesi), 4)
        metrikler['cv_ort'] = round(np.mean(cv_listesi), 4)
    
    # Refractory genel saglik
    refr_listesi = [u['refractory_ihlal_orani'] for u in aktif_unitler]
    metrikler['refractory_medyan'] = round(np.median(refr_listesi), 4)
    metrikler['unit_temiz_orani'] = round(
        sum(1 for r in refr_listesi if r < 0.05) / n_aktif, 4)
    
    # Anlamli ritim (perm p<0.05) gosteren unit sayisi
    perm_listesi = [u['perm_p_deger'] for u in aktif_unitler
                     if u.get('perm_p_deger') is not None]
    if perm_listesi:
        anlamli = sum(1 for p in perm_listesi if p < 0.05)
        metrikler['ritim_unit_sayisi'] = anlamli
        metrikler['ritim_unit_orani'] = round(anlamli / len(perm_listesi), 4)
    
    return metrikler


# ============================================================
# CROSS-UNIT SENKRON
# ============================================================

def cross_unit_senkron(spike_listesi, sure_sn, dt=0.01, max_unit=20):
    """
    En aktif N unit icin STTC matrisi hesaplar.
    Cok unit varsa tumune yapmak yavas - en aktif olanlar alinir.
    """
    # En cok spike'i olan unit'leri sec
    n_unit = len(spike_listesi)
    if n_unit < 2:
        return None
    
    spike_sayisi = [len(sp) for sp in spike_listesi]
    siralama = np.argsort(spike_sayisi)[::-1]  # Cogundan aza
    secilen_idx = siralama[:max_unit]
    secilen_treni = [spike_listesi[i] for i in secilen_idx]
    
    print(f'  STTC matrisi: en aktif {len(secilen_treni)} unit (toplam {n_unit})...')
    
    sonuc = omet.coklu_kanal_senkron(secilen_treni, sure_sn, dt=dt)
    
    if sonuc is None:
        return None
    
    return {
        'matris': sonuc.get('sttc_matrix'),
        'ort_sttc': round(sonuc.get('ortalama_sttc', 0), 4),
        'medyan_sttc': round(sonuc.get('medyan_sttc', 0), 4),
        'max_sttc': round(sonuc.get('max_sttc', 0), 4),
        'min_sttc': round(sonuc.get('min_sttc', 0), 4),
        'kullanilan_unit': len(secilen_treni),
        'secilen_idx': secilen_idx.tolist(),
    }


# ============================================================
# NETWORK BURST TESPITI
# ============================================================

def network_burst_tespit(spike_listesi, sure_sn, bin_ms=50,
                          esik_carpan=3.0, min_aktif_unit=3,
                          min_burst_sure_ms=100):
    """
    Tum unit'lerin spike'lari toplu sayilir (population spike rate).
    Yuksek aktivite donemlerini bul (network bursts).
    
    Esik: medyan + esik_carpan * MAD (adaptif, robust)
    
    Donen: dict
    """
    n_unit = len(spike_listesi)
    if n_unit < min_aktif_unit:
        return {
            'burst_sayisi': 0,
            'burst_orani_per_dakika': 0,
            'durum': f'cok_az_unit ({n_unit})',
            'pop_rate': np.array([]),
            'bin_ms': bin_ms,
        }
    
    # Bin'leme
    bin_sn = bin_ms / 1000.0
    n_bin = int(sure_sn / bin_sn)
    
    pop_rate = np.zeros(n_bin)
    aktif_unit_per_bin = np.zeros(n_bin)
    
    for sp in spike_listesi:
        if len(sp) == 0:
            continue
        bin_idx = (np.array(sp) / bin_sn).astype(int)
        bin_idx = bin_idx[(bin_idx >= 0) & (bin_idx < n_bin)]
        # Population spike rate (her spike sayilir)
        unique_bins, sayi = np.unique(bin_idx, return_counts=True)
        pop_rate[unique_bins] += sayi
        # Aktif unit sayaci (bu unit'in DOKUNDUGU bin'lere +1)
        aktif_unit_per_bin[unique_bins] += 1
    
    # Adaptif esik: medyan + carpan * MAD
    # Sparse veride (cogu bin 0) medyan/MAD calismaz - sifir olmayanlari kullan
    nonzero = pop_rate[pop_rate > 0]
    if len(nonzero) > 10:
        # Yeteri kadar veri varsa MAD ile
        medyan = np.median(nonzero)
        mad = np.median(np.abs(nonzero - medyan))
        if mad < 1:
            mad = max(1, np.std(nonzero))
        esik = medyan + esik_carpan * mad
    else:
        # Cok sparse - basitce sabit esik
        esik = max(2, np.percentile(pop_rate, 95))
    
    # Esik en az min_aktif_unit (mantigen burst en az min_aktif_unit'i icermeli)
    esik = max(esik, min_aktif_unit)
    
    yuksek = pop_rate > esik
    min_burst_sure_sn = min_burst_sure_ms / 1000.0
    
    # Ardisik yuksek bin'leri burst say
    burstler = []
    icinde = False
    bas_idx = 0
    for i, durum in enumerate(yuksek):
        if durum and not icinde:
            icinde = True
            bas_idx = i
        elif not durum and icinde:
            icinde = False
            sure = (i - bas_idx) * bin_sn
            if sure >= min_burst_sure_sn:
                aktif_u = max(aktif_unit_per_bin[bas_idx:i])
                if aktif_u >= min_aktif_unit:
                    burstler.append({
                        'baslangic_sn': bas_idx * bin_sn,
                        'bitis_sn': i * bin_sn,
                        'sure_sn': sure,
                        'aktif_unit': int(aktif_u),
                        'tepe_rate': float(np.max(pop_rate[bas_idx:i])),
                    })
    
    return {
        'burst_sayisi': len(burstler),
        'burst_orani_per_dakika': round(len(burstler) / (sure_sn / 60), 3) if sure_sn > 0 else 0,
        'medyan_burst_sure_sn': round(np.median([b['sure_sn'] for b in burstler]), 3) if burstler else 0,
        'medyan_aktif_unit': int(np.median([b['aktif_unit'] for b in burstler])) if burstler else 0,
        'medyan_tepe_rate': round(np.median([b['tepe_rate'] for b in burstler]), 1) if burstler else 0,
        'esik_kullanilan': round(esik, 2),
        'pop_rate_medyan': round(float(np.median(pop_rate)), 2),
        'pop_rate_max': round(float(np.max(pop_rate)) if len(pop_rate) > 0 else 0, 2),
        'burstler': burstler[:20],
        'pop_rate': pop_rate,
        'bin_ms': bin_ms,
    }


# ============================================================
# GRAFIK URET
# ============================================================

def grafik_uret(spike_listesi, unit_sonuclari, pop_metrik, burst_sonuc,
                 sttc_sonuc, sure_sn, organoid, kayit, cikti_yolu):
    """
    4-panelli grafik:
    1. Raster (tum unit'ler)
    2. Population spike rate + burst tepecikleri
    3. Unit Hz dagilimi histogram
    4. STTC matrisi heatmap
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(14, 10), facecolor='white')
    fig.suptitle(f'{organoid} / {kayit}  —  '
                 f'{pop_metrik["aktif_unit"]}/{pop_metrik["toplam_unit"]} aktif unit',
                 fontsize=12, fontweight='bold')
    
    # Panel 1: Raster (en aktif 30 unit, yoksa hepsi)
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2)
    n_unit = len(spike_listesi)
    spike_sayisi = [len(sp) for sp in spike_listesi]
    siralama = np.argsort(spike_sayisi)[::-1]
    n_goster = min(30, n_unit)
    
    for plot_idx, unit_idx in enumerate(siralama[:n_goster]):
        sp = spike_listesi[unit_idx]
        if len(sp) == 0:
            continue
        ax1.scatter(sp, [plot_idx] * len(sp), s=0.5, c='black', marker='|',
                    alpha=0.7)
    
    ax1.set_xlim(0, sure_sn)
    ax1.set_ylim(-0.5, n_goster - 0.5)
    ax1.set_xlabel('Zaman (sn)')
    ax1.set_ylabel('Unit (en aktif)')
    ax1.set_title(f'Raster - en aktif {n_goster} unit')
    ax1.grid(alpha=0.3)
    
    # Panel 2: Population rate + burst
    ax2 = plt.subplot2grid((3, 2), (1, 0), colspan=2)
    if burst_sonuc and 'pop_rate' in burst_sonuc:
        pop_rate = burst_sonuc['pop_rate']
        bin_ms = burst_sonuc['bin_ms']
        t_bin = np.arange(len(pop_rate)) * bin_ms / 1000.0
        ax2.plot(t_bin, pop_rate, color='steelblue', lw=0.7)
        ax2.fill_between(t_bin, 0, pop_rate, alpha=0.3, color='steelblue')
        # Burst markerlari
        for b in burst_sonuc.get('burstler', []):
            ax2.axvspan(b['baslangic_sn'], b['bitis_sn'],
                        color='red', alpha=0.2)
        ax2.set_xlabel('Zaman (sn)')
        ax2.set_ylabel(f'Pop. spike sayisi ({bin_ms}ms bin)')
        ax2.set_title(f'Population rate + {burst_sonuc["burst_sayisi"]} network burst')
        ax2.grid(alpha=0.3)
    
    # Panel 3: Unit Hz dagilimi
    ax3 = plt.subplot2grid((3, 2), (2, 0))
    aktif_hz = [u['ort_hz'] for u in unit_sonuclari if u['durum'] == 'ok']
    if aktif_hz:
        ax3.hist(aktif_hz, bins=30, color='steelblue', edgecolor='black',
                 alpha=0.7)
        ax3.axvline(np.median(aktif_hz), color='red', ls='--', lw=1,
                    label=f'medyan={np.median(aktif_hz):.2f} Hz')
        ax3.set_xlabel('Unit ortalama Hz')
        ax3.set_ylabel('Unit sayisi')
        ax3.set_title(f'Hz dagilimi ({len(aktif_hz)} aktif unit)')
        ax3.set_yscale('log')
        ax3.legend()
        ax3.grid(alpha=0.3)
    
    # Panel 4: STTC matrisi
    ax4 = plt.subplot2grid((3, 2), (2, 1))
    if sttc_sonuc and sttc_sonuc.get('matris') is not None:
        matris = np.array(sttc_sonuc['matris'])
        im = ax4.imshow(matris, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        ax4.set_xlabel('Unit (sirali)')
        ax4.set_ylabel('Unit (sirali)')
        ax4.set_title(f'STTC matrisi ({sttc_sonuc["kullanilan_unit"]} unit)\n'
                      f'medyan={sttc_sonuc["medyan_sttc"]:.3f}')
        plt.colorbar(im, ax=ax4, fraction=0.046, pad=0.04)
    else:
        ax4.text(0.5, 0.5, 'STTC hesaplanamadi', ha='center', va='center',
                 transform=ax4.transAxes)
        ax4.axis('off')
    
    plt.tight_layout()
    plt.savefig(cikti_yolu, dpi=120, bbox_inches='tight', facecolor='white')
    plt.close()


# ============================================================
# ANA AKIS
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='NWB sorted units icin per-unit analiz')
    parser.add_argument('dosya', help='NWB dosyasi')
    parser.add_argument('organoid', help='Organoid kodu (orn: HO2)')
    parser.add_argument('kayit', help='Kayit kodu (orn: 20250916)')
    parser.add_argument('--out-dir', default='unit_sonuclari',
                        help='Cikti klasoru (varsayilan: unit_sonuclari)')
    parser.add_argument('--min-spike', type=int, default=10,
                        help='Bir unit en az bu kadar spike olmali (varsayilan: 10)')
    parser.add_argument('--top-units', type=int, default=20,
                        help='STTC icin kullanilan unit sayisi (varsayilan: 20)')
    parser.add_argument('--perm-n', type=int, default=200,
                        help='Permutasyon iterasyon (varsayilan: 200)')
    args = parser.parse_args()
    
    if not os.path.exists(args.dosya):
        print(f'HATA: Dosya yok: {args.dosya}')
        return 1
    
    cikti_klasor = Path(args.out_dir)
    cikti_klasor.mkdir(exist_ok=True)
    
    print('=' * 60)
    print(f'  PER-UNIT ANALIZ')
    print('=' * 60)
    print(f'  Dosya    : {args.dosya}')
    print(f'  Organoid : {args.organoid}')
    print(f'  Kayit    : {args.kayit}')
    print()
    
    t0 = time.time()
    
    # ========================================
    # 1. Veri oku
    # ========================================
    print('[1/5] NWB dosyasi okunuyor...')
    
    units = oio.sorted_units_oku(args.dosya)
    if not units:
        print('  UYARI: NWB dosyasinda sorted units tablosu bulunamadi.')
        print()
        print('  Bu dosya ham elektriksel sinyal iceriyor olabilir.')
        print('  Ham sinyal analizi farkli bir islem akisi gerektirir:')
        print('  bandpass filtre + spike detection + spike sorting.')
        print()
        print('  Ham veriyi islemek ister misiniz?')
        print('  1 : Evet — ham sinyal modunu baslat (deneysel)')
        print('  0 : Hayir — cik')
        print()
        
        try:
            secim = input('  Seciminiz (0/1): ').strip()
        except (EOFError, KeyboardInterrupt):
            secim = '0'
        
        if secim == '1':
            print()
            print('  Ham sinyal modu baslatiliyor...')
            print('  organoid_signal.py modulu kullanilacak.')
            print()
            
            # Ham sinyal analizi
            try:
                import organoid_signal as osig
                import organoid_io as oio2
                
                tip = oio2.nwb_tip_tespit(args.dosya)
                if tip not in ('ham', 'her_ikisi'):
                    print('  HATA: Bu dosyada ham sinyal de bulunamadi.')
                    print('  Dosya ne sorted units ne de ham sinyal iceriyor.')
                    return 1
                
                print('  Ham sinyal bulundu. Isleniyor...')
                print('  (Bu islem birkac dakika surebilir)')
                print()
                
                tum_spike_ham = []
                sure_sn = 0.0
                
                for chunk, bas, bit, sr in oio2.ham_chunk_uret(args.dosya, max_sure_sn=180.0):
                    sure_sn = bit
                    for kanal in range(min(4, chunk.shape[1])):
                        sinyal = chunk[:, kanal]
                        spike_idx = osig.spike_tespit_chunk(sinyal, sr)
                        spike_zaman = (spike_idx / sr) + bas
                        tum_spike_ham.extend(spike_zaman.tolist())
                
                if not tum_spike_ham:
                    print('  HATA: Ham sinyalden spike tespit edilemedi.')
                    return 1
                
                print(f'  Tespit edilen spike: {len(tum_spike_ham)}')
                print(f'  Kayit suresi: {sure_sn:.1f} sn')
                print()
                print('  NOT: Ham sinyal modu spike sorting YAPMAZ.')
                print('  Tum spike-lar tek populasyon olarak raporlanir.')
                print('  Tam analiz icin SpikeInterface ile sorting gerekir.')
                print()
                
                # Basit populasyon metrikleri
                tum_spike_ham = sorted(tum_spike_ham)
                spike_hz = len(tum_spike_ham) / sure_sn if sure_sn > 0 else 0
                
                print(f'  Populasyon atisle hizi: {spike_hz:.2f} Hz')
                print(f'  Toplam spike: {len(tum_spike_ham)}')
                print(f'  Kayit suresi: {sure_sn:.1f} sn ({sure_sn/60:.1f} dk)')
                print()
                print('  Ham sinyal analizi tamamlandi.')
                print('  Per-unit analiz icin spike sorting gereklidir.')
                return 0
                
            except Exception as e:
                print(f'  Ham sinyal analizi hatasi: {e}')
                print('  Sorted units olmadan tam analiz yapilamaz.')
                return 1
        else:
            print('  Cikiliyor.')
            return 1
    
    print(f'  {len(units)} unit bulundu')
    
    # Spike listelerini ve kayit suresini cikar
    spike_listesi = [np.array(u['spike_zaman']) for u in units]
    
    # Kayit suresini tahmin et (en buyuk spike zamani)
    tum_spike = []
    for sp in spike_listesi:
        if len(sp) > 0:
            tum_spike.extend(sp.tolist())
    
    if not tum_spike:
        print('  HATA: Hicbir unit\'te spike yok.')
        return 1

    # NaN ve gecersiz degerleri temizle
    import math
    tum_spike_temiz = [s for s in tum_spike
                       if not math.isnan(s) and 0 < s < 1e9]
    if not tum_spike_temiz:
        print('  HATA: Tum spike zamanlari NaN veya gecersiz.')
        return 1

    sure_sn = max(tum_spike_temiz) + 1.0  # Son spike + 1 sn buffer

    # spike_listesi icindeki NaN ve gecersiz degerleri temizle
    spike_listesi = [
        np.array([s for s in sp
                  if not math.isnan(float(s)) and 0 < float(s) < 1e9])
        for sp in spike_listesi
    ]
    print(f'  Toplam spike: {len(tum_spike)}')
    print()
    
    # ========================================
    # 2. Her unit icin ayri analiz
    # ========================================
    print('[2/5] Her unit icin metrik hesabi...')
    unit_sonuclari = []
    for i, sp in enumerate(spike_listesi):
        if (i + 1) % 20 == 0 or i + 1 == len(spike_listesi):
            print(f'  {i+1}/{len(spike_listesi)} unit islendi...', end='\r')
        
        sonuc = unit_analiz(
            unit_id=i,
            spike_zaman=sp,
            sure_sn=sure_sn,
            perm_n=args.perm_n,
        )
        unit_sonuclari.append(sonuc)
    
    print(f'  {len(unit_sonuclari)}/{len(spike_listesi)} unit islendi.' + ' '*20)
    
    aktif_unit_sayi = sum(1 for u in unit_sonuclari if u['durum'] == 'ok')
    print(f'  Aktif (>={args.min_spike} spike): {aktif_unit_sayi} unit')
    print()
    
    # ========================================
    # 3. Cross-unit senkron (STTC)
    # ========================================
    print('[3/5] Cross-unit senkron analizi (STTC)...')
    sttc_sonuc = cross_unit_senkron(spike_listesi, sure_sn,
                                       max_unit=args.top_units)
    if sttc_sonuc:
        print(f'  Medyan STTC: {sttc_sonuc["medyan_sttc"]:.3f}')
        print(f'  Max STTC: {sttc_sonuc["max_sttc"]:.3f}')
        print(f'  Min STTC: {sttc_sonuc["min_sttc"]:.3f}')
    print()
    
    # ========================================
    # 4. Network burst tespiti
    # ========================================
    print('[4/5] Network burst tespiti...')
    burst_sonuc = network_burst_tespit(spike_listesi, sure_sn)
    print(f'  Burst sayisi: {burst_sonuc.get("burst_sayisi", 0)}')
    print(f'  Burst/dakika: {burst_sonuc.get("burst_orani_per_dakika", 0)}')
    if burst_sonuc.get('burst_sayisi', 0) > 0:
        print(f'  Medyan burst suresi: '
              f'{burst_sonuc.get("medyan_burst_sure_sn", 0):.3f} sn')
    print()
    
    # ========================================
    # 5. Populasyon metrikleri + grafik + CSV
    # ========================================
    print('[5/5] Populasyon ozeti, CSV ve grafik...')
    pop_metrik = populasyon_metrikleri(unit_sonuclari, sure_sn)
    
    # CSV: per-unit
    csv_unit_yolu = cikti_klasor / f'{args.organoid}_{args.kayit}_per_unit.csv'
    if unit_sonuclari:
        keys = list(unit_sonuclari[0].keys())
        with open(csv_unit_yolu, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            for u in unit_sonuclari:
                # Eksik alanlari None yap
                satir = {k: u.get(k, '') for k in keys}
                writer.writerow(satir)
        print(f'  Per-unit CSV: {csv_unit_yolu}')
    
    # CSV: ozet
    csv_ozet_yolu = cikti_klasor / f'{args.organoid}_{args.kayit}_ozet.csv'
    ozet = {
        'organoid': args.organoid,
        'kayit': args.kayit,
        'kayit_suresi_sn': round(sure_sn, 2),
        'toplam_unit': len(units),
        'aktif_unit': aktif_unit_sayi,
        'toplam_spike': len(tum_spike),
        **pop_metrik,
        'sttc_medyan': sttc_sonuc.get('medyan_sttc') if sttc_sonuc else None,
        'sttc_ort': sttc_sonuc.get('ort_sttc') if sttc_sonuc else None,
        'sttc_max': sttc_sonuc.get('max_sttc') if sttc_sonuc else None,
        'burst_sayisi': burst_sonuc.get('burst_sayisi', 0),
        'burst_per_dakika': burst_sonuc.get('burst_orani_per_dakika', 0),
        'medyan_burst_sure_sn': burst_sonuc.get('medyan_burst_sure_sn', 0),
    }
    
    with open(csv_ozet_yolu, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=ozet.keys())
        writer.writeheader()
        writer.writerow(ozet)
    print(f'  Ozet CSV: {csv_ozet_yolu}')
    
    # Grafik
    grafik_yolu = cikti_klasor / f'{args.organoid}_{args.kayit}_quicklook.png'
    grafik_uret(spike_listesi, unit_sonuclari, pop_metrik, burst_sonuc,
                 sttc_sonuc, sure_sn, args.organoid, args.kayit, grafik_yolu)
    print(f'  Grafik: {grafik_yolu}')
    
    # Konsol ozeti
    print()
    print('=' * 60)
    print(f'  ORNEK ORGANOID: {args.organoid} / {args.kayit}')
    print('=' * 60)
    print(f'  Kayit suresi    : {sure_sn:.1f} sn ({sure_sn/60:.1f} dakika)')
    print(f'  Toplam unit     : {len(units)}')
    print(f'  Aktif unit      : {aktif_unit_sayi} '
          f'({100*aktif_unit_sayi/len(units):.0f}%)')
    print(f'  Toplam spike    : {len(tum_spike)}')
    print()
    print('  -- Hz Dagilimi (per-unit) --')
    print(f'    Medyan       : {pop_metrik.get("hz_medyan", "N/A")} Hz')
    print(f'    IQR          : {pop_metrik.get("hz_iqr", "N/A")} Hz')
    print(f'    Min - Max    : {pop_metrik.get("hz_min", "N/A")} - '
          f'{pop_metrik.get("hz_max", "N/A")} Hz')
    print()
    print('  -- Kalite (per-unit) --')
    print(f'    Medyan refractory ihlal : {pop_metrik.get("refractory_medyan", "N/A")}')
    print(f'    Temiz unit orani (<%5)   : '
          f'{pop_metrik.get("unit_temiz_orani", "N/A")}')
    print()
    print('  -- Salinim Tespiti --')
    print(f'    Medyan CV    : {pop_metrik.get("cv_medyan", "N/A")}')
    if 'ritim_unit_sayisi' in pop_metrik:
        print(f'    Ritim unit   : {pop_metrik["ritim_unit_sayisi"]} '
              f'({pop_metrik.get("ritim_unit_orani", 0)*100:.0f}%)')
    print()
    print('  -- Cross-Unit Senkron (STTC) --')
    if sttc_sonuc:
        print(f'    Medyan STTC  : {sttc_sonuc["medyan_sttc"]}')
        print(f'    Max STTC     : {sttc_sonuc["max_sttc"]}')
        print(f'    Min STTC     : {sttc_sonuc["min_sttc"]}')
    print()
    print('  -- Network Bursts --')
    print(f'    Burst sayisi      : {burst_sonuc.get("burst_sayisi", 0)}')
    print(f'    Burst/dakika      : {burst_sonuc.get("burst_orani_per_dakika", 0)}')
    if burst_sonuc.get('burst_sayisi', 0) > 0:
        print(f'    Medyan burst sure : '
              f'{burst_sonuc.get("medyan_burst_sure_sn", 0):.3f} sn')
        print(f'    Medyan aktif unit : '
              f'{burst_sonuc.get("medyan_aktif_unit", 0)}')
    print('=' * 60)
    
    sure = time.time() - t0
    print(f'\n  Toplam sure: {sure:.1f} sn')
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
