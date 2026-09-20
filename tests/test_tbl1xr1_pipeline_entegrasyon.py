"""
test_tbl1xr1_pipeline_entegrasyon.py — Gerçek Axon Modülleriyle Entegrasyon Testi
=====================================================================================
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / 'tbl1xr1'))

import numpy as np
from tbl1xr1_synthetic import tbl1xr1_sentetik_organoid_uret
from tbl1xr1_types import Tbl1xr1Genotip
from tbl1xr1_pipeline import (organoid_tam_analiz, kalite_kontrolu_gecti_mi,
                               uclu_metrik_ailesi_hesapla, boylamsal_karsilastir)

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


kontrol = tbl1xr1_sentetik_organoid_uret(19, 60.0, genotip_artis_mi=False, seed=1)
mutant = tbl1xr1_sentetik_organoid_uret(19, 60.0, genotip_artis_mi=True, seed=2)

# 1) QC gercek Axon kodundan
qc_k = kalite_kontrolu_gecti_mi(kontrol)
qc_m = kalite_kontrolu_gecti_mi(mutant)
check("QC (gerçek organoid_qc.refractory_ihlget_orani) ikisi de geçiyor",
      qc_k['gecti'] and qc_m['gecti'], (qc_k, qc_m))

# 2) Uc metrik ailesi de dogru hesaplaniyor ve etiketleniyor
m_k = uclu_metrik_ailesi_hesapla(kontrol, 60.0)
m_m = uclu_metrik_ailesi_hesapla(mutant, 60.0)
check("Birincil metrik (birim-ateşleme hızı) doğru yönde",
      m_m['birincil_ateşleme_hizi']['ortalama_hz'] > m_k['birincil_ateşleme_hizi']['ortalama_hz'])
check("İkincil/üçüncül metrikler 'hipotez testine dahil değil' notuyla dönüyor",
      'DAHİL DEĞİL' in m_k['ikincil_network_burst']['not']
      and 'DAHİL DEĞİL' in m_k['ucuncul_ag_topolojisi']['not'])

# 3) Tam zincir
sonuc_tam = organoid_tam_analiz(mutant, 60.0, Tbl1xr1Genotip.LOF)
check("organoid_tam_analiz uçtan uca çalışıyor", sonuc_tam['guvenilir_mi'] == True)

# 4) KRİTİK: boylamsal_karsilastir split-half NEGATİF kontrolü
#    (aynı kaydın kendi içi -> fark OLMAMALI)
SURE_TAM = 1200.0
tek_kayit = tbl1xr1_sentetik_organoid_uret(19, SURE_TAM, genotip_artis_mi=False, seed=42)
ilk_yari = [sp[sp < 600.0] for sp in tek_kayit]
ikinci_yari = [sp[sp >= 600.0] - 600.0 for sp in tek_kayit]
split_half_sonuc = boylamsal_karsilastir('org1', 'ilkyari', ilk_yari, 600.0,
                                          'ikinciyari', ikinci_yari, 600.0)
check("Split-half NEGATİF kontrol: aynı kaydın iki yarısı arasında fark YOK (p>0.05)",
      split_half_sonuc['hz_mw']['p_deger'] > 0.05, split_half_sonuc['hz_mw'])

# 5) KRİTİK: iki BAĞIMSIZ sentetik "oturum" (aynı dağılımdan farklı seed)
#    -> BİLEREK SAHTE-POZİTİF üretiyor, bunu GİZLEMİYORUZ, açıkça raporluyoruz
hafta10 = tbl1xr1_sentetik_organoid_uret(19, SURE_TAM, genotip_artis_mi=False, seed=10)
hafta14 = tbl1xr1_sentetik_organoid_uret(19, SURE_TAM, genotip_artis_mi=False, seed=14)
bagimsiz_sonuc = boylamsal_karsilastir('org1', 10, hafta10, SURE_TAM, 14, hafta14, SURE_TAM)
check("BİLİNEN SINIRLILIK doğrulandı: bağımsız oturumlar arası sahte-pozitif "
      "(p<<0.05, gerçek fark YOKKEN) — split-half kontrolsüz asla güvenilmemeli",
      bagimsiz_sonuc['hz_mw']['p_deger'] < 0.05, bagimsiz_sonuc['hz_mw'])

print("\n".join(sonuclar))
print(f"\n{'='*70}\nTOPLAM: {passed}/{passed+failed} test geçti "
      f"(5. test 'geçmesi', sınırlılığın gerçekten var olduğunu doğruluyor)\n{'='*70}")
sys.exit(0 if failed == 0 else 1)
