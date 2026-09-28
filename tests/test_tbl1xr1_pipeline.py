"""
test_tbl1xr1_pipeline.py — Uçtan Uca TBL1XR1 Pipeline Testi
=================================================================

Beş genotipin BEŞİ de test ediliyor, ama üçü BİLEREK "hipotez yok"
sonucu vermeli (G70N, L282P, Y446C) — bunlar PASS/FAIL değil,
"pipeline'ın kendini kandırmadığının" testi.
"""
import sys
from pathlib import Path

BASE = Path(__file__).parent.parent
sys.path.insert(0, str(BASE / 'tbl1xr1'))
sys.path.insert(0, str(BASE / 'v6_3'))  # GERÇEK Axon kodu (paket-içi)

import numpy as np
from tbl1xr1_types import Tbl1xr1Genotip
from tbl1xr1_metrics import (ortalama_unit_ateşleme_hizi, genotip_hipotezini_test_et,
                              sure_confound_kontrolu)
from tbl1xr1_synthetic import tbl1xr1_sentetik_organoid_uret

passed, failed = 0, 0
sonuclar = []


def check(isim, kosul, detay=""):
    global passed, failed
    if kosul:
        passed += 1
        sonuclar.append(f"PASS: {isim}")
    else:
        failed += 1
        sonuclar.append(f"FAIL: {isim}  {detay}")


N_CIFT, SURE_SN, N_UNIT = 10, 60.0, 19  # Mastrototaro'nun kendi n'i: 19 hücre

# ─────────────────────────────────────────────────────────────
# BÖLÜM 1: LOF — Mastrototaro'nun DOĞRUDAN sayılarıyla kalibre,
# hipotez testi UYUMLU çıkmalı (bu bir sanity check, gerçek keşif değil)
# ─────────────────────────────────────────────────────────────

kontrol_hz, mutant_hz = [], []
kontrol_sureler, mutant_sureler = [], []

for i in range(N_CIFT):
    k_spikes = tbl1xr1_sentetik_organoid_uret(N_UNIT, SURE_SN, genotip_artis_mi=False, seed=100+i)
    m_spikes = tbl1xr1_sentetik_organoid_uret(N_UNIT, SURE_SN, genotip_artis_mi=True, seed=200+i)
    kontrol_hz.append(ortalama_unit_ateşleme_hizi(k_spikes, SURE_SN)['ortalama_hz'])
    mutant_hz.append(ortalama_unit_ateşleme_hizi(m_spikes, SURE_SN)['ortalama_hz'])
    kontrol_sureler.append(SURE_SN)
    mutant_sureler.append(SURE_SN)

sonuc_lof = genotip_hipotezini_test_et(Tbl1xr1Genotip.LOF, mutant_hz, kontrol_hz)
print(f"LOF: kontrol={np.median(kontrol_hz):.2f} Hz, mutant={np.median(mutant_hz):.2f} Hz, "
      f"uyum={sonuc_lof['uyum']}")
check("LOF: Mastrototaro-kalibreli veri hipotezle UYUMLU çıkıyor (sanity check)",
      sonuc_lof['uyum'] == 'imza_ile_uyumlu', sonuc_lof)
check("LOF: Wilcoxon anlamlı", sonuc_lof['wilcoxon']['p_deger'] < 0.05, sonuc_lof['wilcoxon'])

confound_sonuc = sure_confound_kontrolu(mutant_sureler, kontrol_sureler)
check("LOF: eşit sürelerde confound YOK", confound_sonuc['confound_riski_var'] == False)

# ─────────────────────────────────────────────────────────────
# BÖLÜM 2: BELİRSİZ genotipler (G70N, L282P, Y446C) — sentetik veri
# ÜRETİLEMEMELİ (literatürde sayı yok), ve hipotez testi çağrılırsa
# "kesif" modunda kalmalı, sahte bir uyum/uyumsuzluk İDDİA ETMEMELİ.
# ─────────────────────────────────────────────────────────────

for genotip in (Tbl1xr1Genotip.G70N, Tbl1xr1Genotip.L282P, Tbl1xr1Genotip.Y446C):
    try:
        tbl1xr1_sentetik_organoid_uret(N_UNIT, SURE_SN, genotip_artis_mi=None, seed=1)
        check(f"{genotip.value}: sahte sentetik kalibrasyon REDDEDİLDİ", False,
              "beklenen NotImplementedError fırlamadı")
    except NotImplementedError:
        check(f"{genotip.value}: sahte sentetik kalibrasyon doğru şekilde REDDEDİLDİ", True)

    # Gerçek/varsayımsal ölçüm verisi geldiğini simüle edelim (rastgele,
    # kalibrasyonsuz) ve hipotez testinin "belirsiz" dediğini doğrulayalım
    rastgele_kontrol = np.random.default_rng(5).normal(4.0, 1.0, 8)
    rastgele_mutant = np.random.default_rng(6).normal(4.5, 1.0, 8)
    sonuc = genotip_hipotezini_test_et(genotip, rastgele_mutant, rastgele_kontrol)
    check(f"{genotip.value}: hipotez testi 'belirsiz' diyor (yanlış güven iddia ETMİYOR)",
          sonuc['uyum'] == 'imza_belirsiz_karsilastirilamiyor', sonuc)

# ─────────────────────────────────────────────────────────────
# BÖLÜM 3: BILINMEYEN_MISSENSE — tam keşif modu, hiçbir hipotez iddiası olmamalı
# ─────────────────────────────────────────────────────────────

sonuc_bilinmeyen = genotip_hipotezini_test_et(
    Tbl1xr1Genotip.BILINMEYEN_MISSENSE, mutant_hz, kontrol_hz)
check("BILINMEYEN_MISSENSE: keşif modunda, 'uyum' anahtarı YOK",
      sonuc_bilinmeyen['mod'] == 'kesif' and 'uyum' not in sonuc_bilinmeyen,
      sonuc_bilinmeyen)

# ─────────────────────────────────────────────────────────────
# BÖLÜM 4: Yanlış-yönlü veri testi — pipeline kendini KANDIRMAMALI
# ─────────────────────────────────────────────────────────────

ters_mutant_hz = list(np.random.default_rng(9).normal(2.0, 0.3, N_CIFT))  # LOF beklentisinin tersi
sonuc_ters = genotip_hipotezini_test_et(Tbl1xr1Genotip.LOF, ters_mutant_hz, kontrol_hz)
check("LOF: ters-yönlü veri dürüstçe UYUMSUZ raporlanıyor (zorla 'doğrulanmadı')",
      sonuc_ters['uyum'] == 'imza_ile_uyumsuz', sonuc_ters)

# ─────────────────────────────────────────────────────────────
print("\n" + "\n".join(sonuclar))
print(f"\n{'='*60}\nTOPLAM: {passed}/{passed+failed} test geçti\n{'='*60}")
sys.exit(0 if failed == 0 else 1)
