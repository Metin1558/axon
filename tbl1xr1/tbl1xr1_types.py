"""
tbl1xr1_types.py — TBL1XR1'e Özel Genotip Yapıları
========================================================

Bu, signatures/registry.py'nin (axon-ndd, genel çerçeve) TERSİ bir
tasarım kararı: orada "herhangi bir gen/hastalık" için genel bir
kayıt yapısı vardı. Burada YALNIZCA TBL1XR1'in bilinen 5 mutasyon
sınıfı, literatürden doğrudan aktarılan gerçek bulgularla, SABİT
olarak kodlanıyor. Yeni bir gen eklemek için bu dosyanın yeniden
yazılması gerekir — bu bilinçli bir tercih, genel bir dictionary'ye
gen ismi anahtarlamak değil.

Kaynak (hepsi bu proje kapsamında tam metni okunmuş 5 makaleden):
    - Heinen et al. 2016, J Med Genet (Y446C / Pierpont)
    - Mastrototaro et al. 2021, Front Cell Dev Biol (fare KO + complementation)
    - Nishi et al. 2017, Sci Rep (F10L / şizofreni)
    - Wei et al. 2025, BMC Med Genomics (genotip-fenotip heterojenite)
    - Nagy et al. 2024, Orphanet J Rare Dis (n=41 doğal seyir)
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional


class Tbl1xr1Genotip(Enum):
    """
    TBL1XR1'in literatürde FONKSİYONEL OLARAK test edilmiş 4 nokta
    mutasyonu + gerçek LOF (delesyon/frameshift) + bilinmeyen missense.
    Mastrototaro'nun NSC complementation deneyi bu 4 mutasyonu
    doğrudan test etti — bu enum keyfi değil, o deneyin kendi
    kategorileri.
    """
    LOF = 'lof'                      # gerçek delesyon/frameshift — tam fonksiyon kaybı
    F10L = 'f10l'                    # p.Phe10Leu — şizofreni (Nishi 2017)
    G70N = 'g70n'                    # p.Gly70Asp/Asn — West sendromu (Saitsu 2014)
    L282P = 'l282p'                  # p.Leu282Pro — ASD (O'Roak 2012)
    Y446C = 'y446c'                  # p.Tyr446Cys — Pierpont sendromu (Heinen 2016)
    BILINMEYEN_MISSENSE = 'bilinmeyen_missense'  # WD40 içinde ama yukarıdaki 4'ten biri değil
    KONTROL = 'kontrol'


class KanitSeviyesi(Enum):
    FARE_ELEKTROFIZYOLOJI = 'fare_elektrofizyoloji'      # Mastrototaro'nun patch-clamp verisi
    FARE_COMPLEMENTASYON = 'fare_complementasyon'         # NSC kurtarma deneyi (elektrofizyoloji değil)
    HUCRE_KULTURU_BIYOKIMYA = 'hucre_kulturu_biyokimya'   # Nishi'nin co-IP/TOPFlash verisi
    KLINIK_EEG = 'klinik_eeg'                             # Saitsu/Wei'nin hasta EEG bulguları
    YOK = 'insan_organoid_kaniti_yok'                     # dürüstçe: hiçbiri organoid değil


@dataclass
class MutasyonProfili:
    """Tek bir genotip sınıfı için literatürden çıkan tam profil."""
    genotip: Tbl1xr1Genotip
    klinik_karsilik: str
    npc_proliferasyon_farklilasma_mekanizmasi_var_mi: Optional[bool]
    # ^ Mastrototaro'nun complementation deneyinde KO fenotipini
    #   TEKRAR ÜRETIYOR mu (kurtarmıyor mu)? None = test edilmedi.
    beklenen_hucre_seviyesi_ateşleme_yonu: Optional[str]
    # ^ 'artis'/'azalis'/'belirsiz' — sEPSC frekansına (Mastrototaro)
    #   dayalı, SADECE npc_mekanizmasi_var ise doğrudan destekli
    kanit_seviyesi: KanitSeviyesi
    kaynak: str
    notlar: str


# ─────────────────────────────────────────────────────────────
# 5 PROFİL — hepsi bu projede tam metni okunan makalelerden
# ─────────────────────────────────────────────────────────────

PROFIL_LOF = MutasyonProfili(
    genotip=Tbl1xr1Genotip.LOF,
    klinik_karsilik='MRD41 (otizm, intellektüel yetersizlik, nöbetler) — otistik özellikler VAR',
    npc_proliferasyon_farklilasma_mekanizmasi_var_mi=True,
    beklenen_hucre_seviyesi_ateşleme_yonu='artis',
    kanit_seviyesi=KanitSeviyesi.FARE_ELEKTROFIZYOLOJI,
    kaynak='Mastrototaro et al. 2021 — Tbl1xr1 KO fare, P30 korteks, whole-cell '
           'voltage-clamp: sEPSC frekansı WT 4.0±0.3 Hz vs KO 6.3±0.4 Hz (p<0.001). '
           'Amplitüd/rise-time/decay-time değişmedi. sIPSC farksız.',
    notlar=(
        'BU DOĞRUDAN İNSAN VERİSİ DEĞİL — fare korteks dilimi, akut patch-clamp, '
        'tek hücre seviyesi. İnsan organoidinde MEA-seviyesi (popülasyon) karşılığı '
        'DOĞRULANMAMIŞ bir çeviri varsayımı. En yakın MEA karşılığı network burst '
        'SAYISI değil, birim-başı (per-unit) ortalama ateşleme hızıdır — bkz. '
        'tbl1xr1_metrics.ortalama_unit_ateşleme_hizi.'
    ),
)

PROFIL_F10L = MutasyonProfili(
    genotip=Tbl1xr1Genotip.F10L,
    klinik_karsilik='Şizofreni (sporadik, tek Japon hasta vakası)',
    npc_proliferasyon_farklilasma_mekanizmasi_var_mi=True,  # KO'yu KURTARAMADI — LOF gibi davranıyor
    beklenen_hucre_seviyesi_ateşleme_yonu='artis',  # LOF ile aynı mekanizma varsayımı
    kanit_seviyesi=KanitSeviyesi.FARE_COMPLEMENTASYON,
    kaynak='Mastrototaro 2021: F10L, Tbl1xr1 KO NSC proliferasyon/farklılaşma '
           'defektini KURTARAMADI (tek başarısız olan mutasyon) — LOF ile aynı '
           'kategoride. Nishi et al. 2017: F10L, N-CoR bağlanmasını AZALTIYOR, '
           'β-catenin bağlanmasını ARTIRIYOR (co-IP, 293FT + HT22 hücreleri), '
           'Wnt/β-catenin transkripsiyonel aktivitesi (TOPFlash) ARTIYOR.',
    notlar=(
        'Elektrofizyolojik yön DOĞRUDAN ölçülmedi — sadece "LOF fenotipini '
        'paylaşıyor" mekanizmasından ARTIŞ yönü çıkarsanıyor (complementation '
        'kurtaramama + Wnt aktivasyon artışı üzerinden dolaylı). Gerçek '
        'sEPSC/MEA ölçümü F10L için literatürde YOK.'
    ),
)

PROFIL_G70N = MutasyonProfili(
    genotip=Tbl1xr1Genotip.G70N,
    klinik_karsilik='West sendromu + otistik özellikler (Saitsu 2014, tek vaka), hipsaritmi',
    npc_proliferasyon_farklilasma_mekanizmasi_var_mi=False,  # KO'yu NORMAL KURTARDI
    beklenen_hucre_seviyesi_ateşleme_yonu='belirsiz',
    kanit_seviyesi=KanitSeviyesi.KLINIK_EEG,
    kaynak='Mastrototaro 2021: G70N, NSC proliferasyon/farklılaşma testinde KO '
           'fenotipini NORMAL şekilde kurtardı (yalnızca neurosphere çapında '
           'kısmi kurtarma). Saitsu et al. 2014: hastanın 7 aylıkken EEG '
           'hipsaritmi paterni gösterdiği, West sendromu tanısı aldığı bildirildi.',
    notlar=(
        'ÖNEMLİ: Mastrototaro"nun testi G70N"nin LOF ile AYNI mekanizmayı '
        'PAYLAŞMADIĞINI gösteriyor — yani buradaki elektrofizyolojik beklentiyi '
        'LOF"tan miras almak YANLIŞ olur. Hipsaritmi klinik bir EEG paterni '
        '(kaotik, düzensiz, yüksek genlikli) — insan organoidinde/farede hiç '
        'test edilmemiş, doğrudan bir hücresel mekanizma/yön iddia edilemez. '
        'Bilinçli olarak BELİRSİZ bırakıldı.'
    ),
)

PROFIL_L282P = MutasyonProfili(
    genotip=Tbl1xr1Genotip.L282P,
    klinik_karsilik='Otizm spektrum bozukluğu (O\'Roak 2012, sporadik ASD kohortu)',
    npc_proliferasyon_farklilasma_mekanizmasi_var_mi=False,  # KO'yu NORMAL KURTARDI
    beklenen_hucre_seviyesi_ateşleme_yonu='belirsiz',
    kanit_seviyesi=KanitSeviyesi.KLINIK_EEG,
    kaynak='Mastrototaro 2021: L282P, NSC proliferasyon/farklılaşma testinde KO '
           'fenotipini normal kurtardı. O\'Roak et al. 2012: ASD kohortunda de '
           'novo bulundu, spesifik elektrofizyolojik karakterizasyon yapılmadı.',
    notlar=(
        'G70N ile aynı durum: LOF mekanizmasını paylaşmıyor, elektrofizyolojik '
        'yön için doğrudan hiçbir kanıt (klinik ya da deneysel) yok. BELİRSİZ.'
    ),
)

PROFIL_Y446C = MutasyonProfili(
    genotip=Tbl1xr1Genotip.Y446C,
    klinik_karsilik='Pierpont sendromu — otizm YOK, farklı dismorfik özellikler',
    npc_proliferasyon_farklilasma_mekanizmasi_var_mi=False,  # KO'yu NORMAL KURTARDI
    beklenen_hucre_seviyesi_ateşleme_yonu='belirsiz',
    kanit_seviyesi=KanitSeviyesi.HUCRE_KULTURU_BIYOKIMYA,
    kaynak='Mastrototaro 2021: Y446C, NSC proliferasyon/farklılaşma testinde KO '
           'fenotipini normal kurtardı. Heinen et al. 2016: mutant protein '
           'NCoR/SMRT/HDAC3 kompleksine DOĞRU şekilde katılıyor (yapısal bozulma '
           'yok) — patoloji, henüz tanımlanmamış bir protein-protein '
           'etkileşiminin bozulmasından kaynaklanıyor, dominant-negatif mekanizma.',
    notlar=(
        'Pierpont hastalarında otizm/nöbet profili LOF hastalarından FARKLI — bu '
        'yüzden LOF\'un elektrofizyolojik yönünü buraya taşımak metodolojik hata '
        'olur. BELİRSİZ bırakıldı, veri gelince ilk test edilecek/düzeltilecek.'
    ),
)

PROFIL_KONTROL = MutasyonProfili(
    genotip=Tbl1xr1Genotip.KONTROL,
    klinik_karsilik='İzogenik kontrol / sağlıklı donör',
    npc_proliferasyon_farklilasma_mekanizmasi_var_mi=None,
    beklenen_hucre_seviyesi_ateşleme_yonu=None,
    kanit_seviyesi=KanitSeviyesi.YOK,
    kaynak='—',
    notlar='Referans grup, kendi başına bir beklenti taşımaz.',
)

TUM_PROFILLER = {
    Tbl1xr1Genotip.LOF: PROFIL_LOF,
    Tbl1xr1Genotip.F10L: PROFIL_F10L,
    Tbl1xr1Genotip.G70N: PROFIL_G70N,
    Tbl1xr1Genotip.L282P: PROFIL_L282P,
    Tbl1xr1Genotip.Y446C: PROFIL_Y446C,
    Tbl1xr1Genotip.KONTROL: PROFIL_KONTROL,
}


def profil_getir(genotip: Tbl1xr1Genotip) -> MutasyonProfili:
    """
    BILINMEYEN_MISSENSE için özel davranış: TUM_PROFILLER'da yok,
    çünkü literatür açıkça gösteriyor ki (Wei 2025, Nagy 2024) WD40
    içindeki bilinmeyen bir missense'in hangi profile benzeyeceği
    TAHMİN EDİLEMEZ (aynı domain, hatta aynı rezidü farklı fenotip
    üretebiliyor — S447N vs S447R). Bu yüzden BILINMEYEN_MISSENSE
    için sahte bir profil UYDURMAK yerine hata fırlatılıyor —
    çağıran, bunu açıkça "keşif modu" olarak ele almalı.
    """
    if genotip == Tbl1xr1Genotip.BILINMEYEN_MISSENSE:
        raise ValueError(
            'BILINMEYEN_MISSENSE için literatürden türetilmiş bir profil YOK '
            '(bkz. Wei 2025: aynı WD40 rezidüsü bile zıt fenotip üretebiliyor). '
            'Bu genotip için önceden bir yön varsaymayın — keşif modunda, '
            'hipotez-testi yapmadan, sadece tanımlayıcı istatistik raporlayın.'
        )
    return TUM_PROFILLER[genotip]


@dataclass
class Tbl1xr1Kayit:
    """Tek bir organoid/kayıt için TBL1XR1-özel metadata."""
    organoid_id: str
    genotip: Tbl1xr1Genotip
    izogenik_cift_id: Optional[str] = None
    yas_hafta: Optional[float] = None
    kayit_suresi_sn: Optional[float] = None
    elektrot_sayisi: Optional[int] = None

    def __post_init__(self):
        if self.yas_hafta is not None and self.yas_hafta < 0:
            raise ValueError(f'yas_hafta negatif olamaz: {self.organoid_id}')
        if self.kayit_suresi_sn is not None and self.kayit_suresi_sn <= 0:
            raise ValueError(f'kayit_suresi_sn pozitif olmalı: {self.organoid_id}')
