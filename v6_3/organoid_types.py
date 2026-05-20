"""
organoid_types.py — Veri Taşıyıcı Dataclass'lar
=================================================

Görev: Modüller arası veri aktarımı için tip-güvenli taşıyıcılar.

Önemli: Bu dataclass'lar SADECE veri taşır, davranış içermez.
Modüller hâlâ bağımsız çalışır — sınıf metodu yok, business logic yok.
Sadece "5 parametre yerine 1 dataclass" sadeleştirmesi.

Kullanım örneği:
    kayit = KayitBilgisi(organoid='SO1', kayit='13_Haz', sr=30000.0)
    metric = ttest_iki_grup(g1, g2)  # eskiden 5 arg, dataclass değil
    qc = SinyalKalitesi(snr=6.92, gurultu_std_uv=20.5, ...)

Modüller bu dataclass'ları opsiyonel olarak alabilir veya hâlâ
ayrı parametrelerle çalışabilir. Geriye uyumluluk korunur.
"""

from dataclasses import dataclass, field
from typing import Optional, List
import numpy as np


@dataclass
class KayitBilgisi:
    """
    Bir kaydın temel meta bilgileri.
    NWB metadata + kullanıcı tanımlamaları.
    """
    organoid: str
    kayit: str
    dosya_yolu: Optional[str] = None
    sr: Optional[float] = None
    kayit_suresi_sn: Optional[float] = None
    elektrot_sayisi: Optional[int] = None
    unit_sayisi: Optional[int] = 0
    kanal_indeks: int = 0
    session_start_time: Optional[str] = None
    stim_var: bool = False
    trial_var: bool = False

    @property
    def etiket(self):
        """Standart etiket: 'organoid/kayit'."""
        return f'{self.organoid}/{self.kayit}'


@dataclass
class SpikeVerisi:
    """
    Spike zamanları + nasıl üretildiği.

    Kaynak (mod):
        'mua'           — threshold-based MUA
        'sort'          — Tridesclous2 sorting
        'sort_or_mua'   — sort denenip MUA fallback (etiketli)
        'nwb_units'     — NWB içindeki sorted units
    """
    spike_zaman: np.ndarray
    sure_sn: float
    mod: str = 'mua'  # 'mua', 'sort', 'sort_or_mua', 'nwb_units'
    sigma_esik: Optional[float] = None
    refractory_ms: Optional[float] = None
    sadece_negatif: bool = False
    mad_esik: bool = False
    n_unit: int = 1

    @property
    def n_spike(self):
        return int(len(self.spike_zaman))

    @property
    def ortalama_hz(self):
        if self.sure_sn <= 0:
            return 0.0
        return self.n_spike / self.sure_sn


@dataclass
class SinyalKalitesi:
    """
    QC metrikleri — sayısal, etiketsiz.
    """
    snr: Optional[float] = None
    gurultu_std_uv: Optional[float] = None
    gurultu_yontem: str = 'std'  # 'std' veya 'mad'
    doygun_oran: float = 0.0
    refractory_ihlal_orani: float = 0.0
    refractory_ihlal_sayisi: int = 0
    toplam_isi: int = 0
    refractory_ci_alt: Optional[float] = None
    refractory_ci_ust: Optional[float] = None
    negatif_isi_sayisi: int = 0
    yogunluk_outlier_orani: float = 0.0
    yogunluk_outlier_sayisi: int = 0


@dataclass
class AnalizParametre:
    """
    Bir analiz çalışmasının ayar bilgileri.
    Reproducibility için saklamakta fayda var.
    """
    snr_sigma: float = 5.0
    bandpass_low: float = 300.0
    bandpass_high: float = 3000.0
    chunk_sn: int = 60
    refractory_ms: float = 1.0
    pencere_sn: int = 30
    sadece_negatif: bool = False
    mad_esik: bool = False
    sort: bool = False
    sort_or_mua: bool = False
    versiyon: str = 'organoid v6.3'
