"""
tbl1xr1_pipeline.py — Tam Entegrasyon Katmanı
====================================================

ÖNCEKİ 4 dosya (types, metrics, synthetic, test) SADECE hipotez-testi
istatistiğini kapsıyordu — haklı eleştiri buydu. Bu dosya, gerçek
Axon'un (v6_3/ ve axon/) 24 dosyalık gerçek işlevselliğini KOPYALAMADAN,
doğrudan import ederek TBL1XR1 analizine bağlıyor.

Hangi gerçek Axon modülü, ne için, TBL1XR1 bağlamında:

    organoid_qc.refractory_ihlget_orani
        -> KALİTE KONTROL KAPISI. Genotip karşılaştırmasından ÖNCE
           çalışmalı — eğer mutant kanallarda refractory-ihlal oranı
           sistematik farklıysa, bu "genotip etkisi" değil "veri
           kalitesi farkı" anlamına gelebilir. axon-ndd'deki
           süre-confound kontrolüyle AYNI mantık, farklı bir gizli
           değişken için.

    graph_analiz.sttc_matrisi_hesapla + graph_metrikleri_hesapla
        -> ÜÇÜNCÜL metrik ailesi (birincil: birim-ateşleme hızı,
           ikincil: network burst). Axon'un kendi DANDI:001603
           validasyonunda YAŞA BAĞLI küçük-dünya topolojisi (σ)
           farkı bulunmuştu — bu, TBL1XR1'de de bağımsız bir sinyal
           taşıyabilir (nörogelişimsel hastalıklarda ağ topolojisi
           bozukluğu genel bir bulgu paterni).

    organoid_compare.ikili_karsilastir (permütasyon testi)
        -> BOYLAMSAL karşılaştırma — aynı organoidin farklı
           haftalardaki kayıtlarını karşılaştırmak için. Nagy 2024'ün
           bulduğu "nöbet başlangıcı yaşa göre değişiyor" ve
           MECP2'de gördüğümüz "erken/geç farklı yön" deseni göz
           önüne alınırsa, TEK bir zaman noktasına güvenmek riskli —
           bu fonksiyon "fenotip hangi haftada ortaya çıkıyor"
           sorusuna cevap arar.

    yas_metadata_cek.metadata_oku
        -> NWB dosyasından gerçek yaş bilgisini okur (organoid_types_ext
           gibi elle girilen bir alan değil, kayıttan doğrudan çekilen).

Bilerek BURAYA DAHIL ETMEDİKLERİM (nedeniyle birlikte):
    - organoid_signal.py / organoid_lfp.py: ham voltaj izine ihtiyaç
      duyuyor, elimizde sentetik ham voltaj yok (sadece spike zamanı
      sentezledik). Arayüzleri burada dokümante edildi ama gerçek ham
      veri gelmeden test edilemez — sahte voltaj izi uydurmak, daha
      önce sahte sayı uydurmamaya verdiğimiz önemle çelişirdi.
    - organoid_sorting.py / si_sorting.py: ACURARE'nin ön-sıralanmış
      spike verisi mi yoksa ham dalga formu mu göndereceği belli
      olmadan hangi moda ihtiyaç duyulacağı bilinmiyor.
    - organoid_db.py / organoid_output.py / batch_analiz.py /
      organoid_cli.py: bunlar TBL1XR1'e özel bilim değil, çok-organoidli
      gerçek bir çalışma başladığında OLDUĞU GİBİ (değiştirilmeden)
      kullanılacak altyapı — yeniden yazılmaya gerek yok.
"""
import sys
from pathlib import Path

# Paket-içi göreli yollar — bu dosya tbl1xr1/ altında, gerçek Axon
# kodu bir üst dizindeki v6_3/ ve axon/ klasörlerinde duruyor.
_KOK = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_KOK / 'v6_3'))
sys.path.insert(0, str(_KOK / 'axon'))
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
from organoid_qc import refractory_ihlget_orani          # GERÇEK Axon, dokunulmadı
from graph_analiz import sttc_matrisi_hesapla, graph_metrikleri_hesapla  # GERÇEK Axon
from organoid_units_analiz import network_burst_tespit    # GERÇEK Axon
from organoid_compare import ikili_karsilastir             # GERÇEK Axon
from tbl1xr1_metrics import ortalama_unit_ateşleme_hizi, genotip_hipotezini_test_et
from tbl1xr1_types import Tbl1xr1Genotip


def kalite_kontrolu_gecti_mi(spike_listesi, refractory_esik=0.02):
    """
    Genotip karşılaştırmasından ÖNCE her organoid için QC.

    refractory_esik: kabul edilebilir max refractory-ihlal oranı.
    Bu eşik keyfi — literatürden gelmiyor, ama organoid_qc'nin
    kendisi de "iyi/kötü" demiyor, sadece sayı veriyor; eşiği
    burada AÇIKÇA biz koyuyoruz, gizli değil.
    """
    ihlal_oranlari = []
    for sp in spike_listesi:
        r = refractory_ihlget_orani(sp)
        ihlal_oranlari.append(r['ihlal_orani'])
    ortalama_ihlal = float(np.mean(ihlal_oranlari)) if ihlal_oranlari else None
    return {
        'gecti': ortalama_ihlal is not None and ortalama_ihlal <= refractory_esik,
        'ortalama_ihlal_orani': ortalama_ihlal,
        'unit_bazli_ihlal_oranlari': ihlal_oranlari,
        'esik': refractory_esik,
    }


def uclu_metrik_ailesi_hesapla(spike_listesi, sure_sn):
    """
    Üç metrik ailesini BİRLİKTE, ama AÇIKÇA ETİKETLENMİŞ olarak
    hesaplar — hangisinin hipotez testine gireceği (birincil),
    hangisinin sadece bilgi amaçlı olduğu (ikincil, üçüncül)
    karıştırılmasın diye.
    """
    # Birincil (hipotez testine giren) — Mastrototaro sEPSC çevirisi
    birincil = ortalama_unit_ateşleme_hizi(spike_listesi, sure_sn)

    # İkincil (keşifsel) — axon-ndd'de kör nokta olduğu kanıtlanmış,
    # yine de bilgi amaçlı hesaplanıyor
    ikincil = network_burst_tespit(spike_listesi, sure_sn)

    # Üçüncül (keşifsel) — Axon'un kendi doğrulanmış yaş-bağlı bulgusu
    # (küçük-dünya σ) buraya da taşınabilir mi, bağımsız bir soru
    matris, idx = sttc_matrisi_hesapla(spike_listesi, sure_sn)
    uclu = graph_metrikleri_hesapla(matris)

    return {
        'birincil_ateşleme_hizi': birincil,
        'ikincil_network_burst': {
            'burst_orani_per_dakika': ikincil['burst_orani_per_dakika'],
            'not': 'Hipotez testine DAHİL DEĞİL — bkz. axon-ndd kör-nokta bulgusu.',
        },
        'ucuncul_ag_topolojisi': {
            'sigma': uclu.get('sigma'),
            'kucuk_dunya': uclu.get('kucuk_dunya'),
            'yogunluk': uclu.get('yogunluk'),
            'not': ('Hipotez testine DAHİL DEĞİL — TBL1XR1 için literatürde '
                    'network-topoloji beklentisi yok, sadece Axon\'un kendi '
                    'DANDI validasyonunda yaşa duyarlı olduğu bilinen bir '
                    'metrik, keşifsel amaçlı taşınıyor.'),
        },
    }


def organoid_tam_analiz(spike_listesi, sure_sn, genotip: Tbl1xr1Genotip,
                         refractory_esik=0.02):
    """
    ANA GİRİŞ NOKTASI — tek bir organoid için tam analiz zinciri:
    QC -> üç metrik ailesi -> (eğer QC geçildiyse) hipotez testine hazır çıktı.

    Hipotez testi burada YAPILMIYOR (o, izogenik ÇİFT gerektiriyor,
    tek organoid değil) — bu fonksiyon sadece tek-organoid seviyesinde
    QC + metrik hesaplama zincirini tamamlıyor.
    """
    qc = kalite_kontrolu_gecti_mi(spike_listesi, refractory_esik)
    metrikler = uclu_metrik_ailesi_hesapla(spike_listesi, sure_sn)
    return {
        'genotip': genotip.value,
        'qc': qc,
        'metrikler': metrikler,
        'guvenilir_mi': qc['gecti'],
    }


def _oturum_sozlugu_olustur(organoid_id, kayit_adi, spike_listesi, sure_sn):
    """
    organoid_compare.ikili_karsilastir tek bir POPÜLASYON spike
    treni bekliyor (çoklu-unit listesi değil) — çok-unit veriyi
    tek bir havuzlanmış (pooled) popülasyon treninde birleştiriyoruz.
    Bu bilinçli bir basitleştirme: birim-kimliği bilgisini bu
    karşılaştırma için feda ediyoruz, sadece "popülasyon aktivite
    örüntüsü zamanla nasıl değişiyor" sorusuna cevap arıyoruz.
    """
    havuz = np.sort(np.concatenate(spike_listesi)) if spike_listesi else np.array([])
    isi = np.diff(havuz) if len(havuz) > 1 else np.array([])
    ortalama_hz = len(havuz) / sure_sn if sure_sn > 0 else None
    cv = float(np.std(isi) / np.mean(isi)) if len(isi) > 0 and np.mean(isi) > 0 else None
    return {
        'organoid': organoid_id, 'kayit': kayit_adi, 'spike_zaman': havuz,
        'sure_sn': sure_sn, 'ortalama_hz': ortalama_hz, 'cv': cv,
        'spike_sayisi': len(havuz),
    }


def boylamsal_karsilastir(organoid_id, hafta1, spike_listesi1, sure_sn1,
                           hafta2, spike_listesi2, sure_sn2):
    """
    AYNI organoidin iki farklı haftadaki kaydını karşılaştırır —
    "fenotip hangi haftada ortaya çıkıyor" sorusu için (Nagy 2024'ün
    yaşa-bağlı nöbet başlangıcı bulgusu + MECP2'de gördüğümüz
    erken/geç yön-değişimi deseni bunu gerekli kılıyor).

    GERÇEK organoid_compare.ikili_karsilastir'ı çağırır, dokunulmadı.

    !!! GERÇEK BİR TEST SONUCU — KRİTİK SINIRLILIK !!!
    Bu fonksiyonu iki BAĞIMSIZ üretilmiş sentetik "oturum" ile test
    ettiğimde (aynı dağılımdan, farklı seed) p=4.8e-20 gibi SAHTE bir
    anlamlılık çıktı. Aynı kaydı gerçekten ikiye bölüp (aynı hücreler,
    aynı süreç) test ettiğimde p=0.70 — doğru şekilde anlamsız çıktı.

    Sonuç: `ikili_karsilastir`, TEK bir kayıt İÇİNDEKİ zamansal
    değişimi (drift/trend) tespit etmek için doğru çalışıyor, ama
    İKİ AYRI oturumu (gerçekte de elektrot kayması, farklı yakalanan
    birimler gibi oturumlar-arası varyasyon kaynakları olacağından)
    karşılaştırırken YANLIŞ-POZİTİF üretebilir — çünkü pencereler
    arası otokorelasyon, oturum-içi varyansı doğru modelliyor ama
    oturumlar-arası ek varyans kaynağını hesaba katmıyor.

    ZORUNLU KULLANIM KURALI: Gerçek hafta10-vs-hafta14 karşılaştırması
    yapılmadan ÖNCE, HER organoid için kendi split-half negatif
    kontrolü (bkz. tests/test_tbl1xr1_pipeline.py) çalıştırılmalı —
    bu gürültü tabanını göstermeden çıkan hiçbir "haftalar arası fark
    bulundu" iddiasına güvenilmemeli.
    """
    k1 = _oturum_sozlugu_olustur(organoid_id, f'hafta{hafta1}', spike_listesi1, sure_sn1)
    k2 = _oturum_sozlugu_olustur(organoid_id, f'hafta{hafta2}', spike_listesi2, sure_sn2)
    sonuc = ikili_karsilastir(k1, k2)
    sonuc['not'] = (
        'UYARI: Bu karşılaştırma, oturumlar-arası varyansı hesaba katmıyor — '
        'güvenmeden önce split-half negatif kontrolü çalıştırın (bkz. docstring). '
        'Popülasyon-havuzlanmış spike treni üzerinde çalışıyor, birim-kimliği '
        'feda edildi. Hipotez testi (genotip_hipotezini_test_et) ile '
        'KARIŞTIRILMAMALI, o izogenik ÇİFT arası, bu AYNI organoidin zaman-içi.'
    )
    return sonuc
