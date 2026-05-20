"""
organoid_db.py — Veri Depolama
================================

İki ayrı amaç, iki ayrı format:
  - SQLite: ölçülen sayısal veriler (metrikler, test sonuçları)
            Spike zamanları artık BLOB olarak (numpy serialize)
  - JSON  : metadata (organoid yaşı, besiyeri, elektrot tipi, deney koşulları)

v6.3 değişikliği: spike_zamanlari tablosu yerine spike_arrays tablosu.
Her unit/kayıt tek bir satır, spike zamanları BLOB olarak.
np.frombuffer ile direkt numpy array'e açılır — milyonlarca spike için
hem disk hem hız avantajı.
"""

import io as ioo
import os
import json
import sqlite3
from datetime import datetime
import numpy as np


# ─────────────────────────────────────────────────────────────
# SQLITE — ÖLÇÜM VERİSİ
# ─────────────────────────────────────────────────────────────

def db_ac(db_yolu='organoid_v6.db'):
    """
    SQLite bağlantısı aç, tabloları oluştur.

    PRAGMA ayarları:
        journal_mode = WAL  : Yoğun yazma için
        synchronous = NORMAL: Güvenli + hızlı
        temp_store = MEMORY : Geçici tablolar RAM'de
    """
    conn = sqlite3.connect(db_yolu)
    conn.row_factory = sqlite3.Row

    conn.execute('PRAGMA journal_mode = WAL')
    conn.execute('PRAGMA synchronous = NORMAL')
    conn.execute('PRAGMA temp_store = MEMORY')
    conn.execute('PRAGMA cache_size = -32000')

    conn.executescript("""
        CREATE TABLE IF NOT EXISTS kayitlar (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            organoid TEXT NOT NULL,
            kayit TEXT NOT NULL,
            dosya_yolu TEXT,
            tarih_analiz TEXT,
            versiyon TEXT,
            kayit_suresi_sn REAL,
            sr REAL,
            elektrot_sayisi INTEGER,
            unit_sayisi INTEGER,
            kanal_indeks INTEGER,
            UNIQUE(organoid, kayit, kanal_indeks)
        );

        CREATE TABLE IF NOT EXISTS spike_arrays (
            kayit_id INTEGER REFERENCES kayitlar(id),
            unit_id INTEGER,
            n_spike INTEGER,
            dtype TEXT,
            spike_blob BLOB,
            PRIMARY KEY (kayit_id, unit_id)
        );

        CREATE TABLE IF NOT EXISTS metrikler (
            kayit_id INTEGER PRIMARY KEY REFERENCES kayitlar(id),
            spike_sayisi INTEGER,
            ortalama_hz REAL,
            ort_isi REAL,
            std_isi REAL,
            medyan_isi REAL,
            cv REAL,
            isi_q25 REAL,
            isi_q75 REAL,
            isi_min REAL,
            isi_max REAL,
            aktiflik_oran REAL,
            aktif_pencere INTEGER,
            toplam_pencere INTEGER
        );

        CREATE TABLE IF NOT EXISTS qc_metrikleri (
            kayit_id INTEGER PRIMARY KEY REFERENCES kayitlar(id),
            snr REAL,
            gurultu_std_uv REAL,
            gurultu_yontem TEXT,
            doygun_oran REAL,
            refractory_ihlal_orani REAL,
            refractory_ihlal_sayisi INTEGER,
            refractory_ci_alt REAL,
            refractory_ci_ust REAL,
            negatif_isi_sayisi INTEGER,
            yogunluk_outlier_orani REAL,
            yogunluk_outlier_sayisi INTEGER,
            ortalama_yogunluk REAL,
            max_yogunluk INTEGER
        );

        CREATE TABLE IF NOT EXISTS spektral_sonuc (
            kayit_id INTEGER REFERENCES kayitlar(id),
            yontem TEXT,
            siralama INTEGER,
            periyot_sn REAL,
            guc REAL
        );

        CREATE TABLE IF NOT EXISTS permutasyon_sonuc (
            kayit_id INTEGER PRIMARY KEY REFERENCES kayitlar(id),
            yontem TEXT,
            permutasyon_n INTEGER,
            seed INTEGER,
            p_deger REAL,
            gercek_guc REAL,
            rastgele_guc_medyan REAL,
            rastgele_guc_q95 REAL,
            eslesen_permutasyon INTEGER
        );

        CREATE TABLE IF NOT EXISTS stim_sonuc (
            kayit_id INTEGER PRIMARY KEY REFERENCES kayitlar(id),
            stim_tipi TEXT,
            stim_bas REAL,
            stim_bit REAL,
            once_n INTEGER,
            sirasinda_n INTEGER,
            sonra_n INTEGER,
            once_hz REAL,
            sirasinda_hz REAL,
            sonra_hz REAL,
            sonra_once_orani REAL,
            ttest_p REAL,
            ttest_t REAL,
            mw_p REAL,
            mw_u REAL,
            n_var INTEGER,
            n_yok INTEGER
        );

        CREATE INDEX IF NOT EXISTS idx_organoid ON kayitlar(organoid);
        CREATE INDEX IF NOT EXISTS idx_spike_kayit ON spike_arrays(kayit_id);
    """)
    conn.commit()
    return conn


# ─────────────────────────────────────────────────────────────
# SPIKE BLOB DEPOLAMA
# ─────────────────────────────────────────────────────────────

def spike_zamanlari_yaz(conn, kayit_id, spike_zaman, unit_id=-1):
    """
    Spike zamanlarını BLOB olarak kaydeder (np.tobytes).

    v6.3 değişikliği: her spike için ayrı satır yerine tek satır BLOB.
    Milyon spike için satır sayısı 1.000.000'dan 1'e düşer.

    Okuma için spike_zamanlari_oku() kullan, np.frombuffer ile açılır.
    """
    spike_zaman = np.asarray(spike_zaman, dtype=np.float64)
    blob = spike_zaman.tobytes()
    n = len(spike_zaman)

    with conn:
        conn.execute("""
            INSERT INTO spike_arrays (kayit_id, unit_id, n_spike, dtype, spike_blob)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT(kayit_id, unit_id) DO UPDATE SET
                n_spike=excluded.n_spike,
                dtype=excluded.dtype,
                spike_blob=excluded.spike_blob
        """, (kayit_id, unit_id, n, 'float64', blob))


def spike_zamanlari_oku_blob(conn, kayit_id, unit_id=-1):
    """
    BLOB'dan spike zamanlarını numpy array olarak okur.
    np.frombuffer — kopya yapmadan açar.
    """
    row = conn.execute(
        "SELECT spike_blob, dtype, n_spike FROM spike_arrays "
        "WHERE kayit_id=? AND unit_id=?",
        (kayit_id, unit_id)).fetchone()
    if row is None:
        return np.array([], dtype=np.float64)
    arr = np.frombuffer(row['spike_blob'], dtype=row['dtype'] or 'float64')
    return arr.copy()  # kopya — yazılabilir olsun


def vacuum(db_yolu='organoid_v6.db'):
    """
    Veritabanını VACUUM ile sıkıştırır.
    BLOB depolama ile zaten az gerekli, ama şişme ihtimali için.
    """
    conn = sqlite3.connect(db_yolu)
    conn.execute('VACUUM')
    conn.close()


def kayit_ekle_veya_guncelle(conn, organoid, kayit, dosya_yolu,
                                meta, kanal_indeks=0,
                                versiyon='organoid v6'):
    """Ana kayıt satırını ekler/günceller, kayit_id döndürür."""
    tarih = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    conn.execute("""
        INSERT INTO kayitlar
            (organoid, kayit, dosya_yolu, tarih_analiz, versiyon,
             kayit_suresi_sn, sr, elektrot_sayisi, unit_sayisi, kanal_indeks)
        VALUES (?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(organoid, kayit, kanal_indeks) DO UPDATE SET
            dosya_yolu=excluded.dosya_yolu,
            tarih_analiz=excluded.tarih_analiz,
            versiyon=excluded.versiyon,
            kayit_suresi_sn=excluded.kayit_suresi_sn,
            sr=excluded.sr,
            elektrot_sayisi=excluded.elektrot_sayisi,
            unit_sayisi=excluded.unit_sayisi
    """, (organoid, kayit, dosya_yolu, tarih, versiyon,
          meta.get('kayit_suresi_sn'),
          meta.get('sr'),
          meta.get('elektrot_sayisi'),
          meta.get('unit_sayisi'),
          kanal_indeks))
    conn.commit()
    row = conn.execute(
        "SELECT id FROM kayitlar WHERE organoid=? AND kayit=? AND kanal_indeks=?",
        (organoid, kayit, kanal_indeks)).fetchone()
    return row['id']


def metrikler_yaz(conn, kayit_id, isi_metrik, aktiflik_oran,
                   aktif_pencere, toplam_pencere, spike_sayisi,
                   ortalama_hz):
    """ISI metriklerini ve aktiflik oranını yazar."""
    conn.execute("""
        INSERT INTO metrikler
            (kayit_id, spike_sayisi, ortalama_hz, ort_isi, std_isi, medyan_isi,
             cv, isi_q25, isi_q75, isi_min, isi_max,
             aktiflik_oran, aktif_pencere, toplam_pencere)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(kayit_id) DO UPDATE SET
            spike_sayisi=excluded.spike_sayisi,
            ortalama_hz=excluded.ortalama_hz,
            ort_isi=excluded.ort_isi,
            std_isi=excluded.std_isi,
            medyan_isi=excluded.medyan_isi,
            cv=excluded.cv,
            isi_q25=excluded.isi_q25,
            isi_q75=excluded.isi_q75,
            isi_min=excluded.isi_min,
            isi_max=excluded.isi_max,
            aktiflik_oran=excluded.aktiflik_oran,
            aktif_pencere=excluded.aktif_pencere,
            toplam_pencere=excluded.toplam_pencere
    """, (kayit_id, spike_sayisi, ortalama_hz,
          isi_metrik.get('ort_isi'), isi_metrik.get('std_isi'),
          isi_metrik.get('medyan_isi'), isi_metrik.get('cv'),
          isi_metrik.get('isi_q25'), isi_metrik.get('isi_q75'),
          isi_metrik.get('isi_min'), isi_metrik.get('isi_max'),
          aktiflik_oran, aktif_pencere, toplam_pencere))
    conn.commit()


def qc_yaz(conn, kayit_id, qc):
    """QC metriklerini yazar (CI alanları + gurultu_yontem dahil)."""
    conn.execute("""
        INSERT INTO qc_metrikleri
            (kayit_id, snr, gurultu_std_uv, gurultu_yontem, doygun_oran,
             refractory_ihlal_orani, refractory_ihlal_sayisi,
             refractory_ci_alt, refractory_ci_ust,
             negatif_isi_sayisi,
             yogunluk_outlier_orani, yogunluk_outlier_sayisi,
             ortalama_yogunluk, max_yogunluk)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(kayit_id) DO UPDATE SET
            snr=excluded.snr,
            gurultu_std_uv=excluded.gurultu_std_uv,
            gurultu_yontem=excluded.gurultu_yontem,
            doygun_oran=excluded.doygun_oran,
            refractory_ihlal_orani=excluded.refractory_ihlal_orani,
            refractory_ihlal_sayisi=excluded.refractory_ihlal_sayisi,
            refractory_ci_alt=excluded.refractory_ci_alt,
            refractory_ci_ust=excluded.refractory_ci_ust,
            negatif_isi_sayisi=excluded.negatif_isi_sayisi,
            yogunluk_outlier_orani=excluded.yogunluk_outlier_orani,
            yogunluk_outlier_sayisi=excluded.yogunluk_outlier_sayisi,
            ortalama_yogunluk=excluded.ortalama_yogunluk,
            max_yogunluk=excluded.max_yogunluk
    """, (kayit_id, qc.get('snr'), qc.get('gurultu_std_uv'),
          qc.get('gurultu_yontem', 'std'),
          qc.get('doygun_oran'),
          qc.get('refractory_ihlal_orani'),
          qc.get('refractory_ihlal_sayisi'),
          qc.get('refractory_ci_alt'),
          qc.get('refractory_ci_ust'),
          qc.get('negatif_isi_sayisi'),
          qc.get('yogunluk_outlier_orani'),
          qc.get('yogunluk_outlier_sayisi'),
          qc.get('ortalama_yogunluk'),
          qc.get('max_yogunluk')))
    conn.commit()


def spektral_yaz(conn, kayit_id, fft_list, welch_list):
    """Spektral analiz sonuçlarını yazar."""
    conn.execute("DELETE FROM spektral_sonuc WHERE kayit_id=?", (kayit_id,))
    rows = []
    for i, (p, g) in enumerate(fft_list):
        rows.append((kayit_id, 'fft', i, p, g))
    for i, (p, g) in enumerate(welch_list):
        rows.append((kayit_id, 'welch', i, p, g))
    if rows:
        conn.executemany(
            "INSERT INTO spektral_sonuc "
            "(kayit_id, yontem, siralama, periyot_sn, guc) VALUES (?,?,?,?,?)",
            rows)
    conn.commit()


def permutasyon_yaz(conn, kayit_id, perm):
    """ISI shuffle permütasyon sonucunu yazar."""
    if perm is None:
        return
    conn.execute("""
        INSERT INTO permutasyon_sonuc
            (kayit_id, yontem, permutasyon_n, seed, p_deger,
             gercek_guc, rastgele_guc_medyan, rastgele_guc_q95,
             eslesen_permutasyon)
        VALUES (?,?,?,?,?,?,?,?,?)
        ON CONFLICT(kayit_id) DO UPDATE SET
            yontem=excluded.yontem,
            permutasyon_n=excluded.permutasyon_n,
            seed=excluded.seed,
            p_deger=excluded.p_deger,
            gercek_guc=excluded.gercek_guc,
            rastgele_guc_medyan=excluded.rastgele_guc_medyan,
            rastgele_guc_q95=excluded.rastgele_guc_q95,
            eslesen_permutasyon=excluded.eslesen_permutasyon
    """, (kayit_id, perm.get('yontem'), perm.get('permutasyon_n'),
          perm.get('seed'), perm.get('p_deger'),
          perm.get('gercek_guc'), perm.get('rastgele_guc_medyan'),
          perm.get('rastgele_guc_q95'),
          perm.get('eslesen_permutasyon')))
    conn.commit()


def stim_yaz(conn, kayit_id, stim_tipi, stim):
    """Stimülasyon analiz sonucunu yazar."""
    if stim is None:
        return
    ttest = stim.get('ttest') or {}
    mw = stim.get('mannwhitney') or {}
    conn.execute("""
        INSERT INTO stim_sonuc
            (kayit_id, stim_tipi, stim_bas, stim_bit,
             once_n, sirasinda_n, sonra_n,
             once_hz, sirasinda_hz, sonra_hz, sonra_once_orani,
             ttest_p, ttest_t, mw_p, mw_u, n_var, n_yok)
        VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
        ON CONFLICT(kayit_id) DO UPDATE SET
            stim_tipi=excluded.stim_tipi,
            stim_bas=excluded.stim_bas,
            stim_bit=excluded.stim_bit,
            once_n=excluded.once_n,
            sirasinda_n=excluded.sirasinda_n,
            sonra_n=excluded.sonra_n,
            once_hz=excluded.once_hz,
            sirasinda_hz=excluded.sirasinda_hz,
            sonra_hz=excluded.sonra_hz,
            sonra_once_orani=excluded.sonra_once_orani,
            ttest_p=excluded.ttest_p,
            ttest_t=excluded.ttest_t,
            mw_p=excluded.mw_p,
            mw_u=excluded.mw_u,
            n_var=excluded.n_var,
            n_yok=excluded.n_yok
    """, (kayit_id, stim_tipi, stim.get('stim_bas'), stim.get('stim_bit'),
          stim.get('once_n'), stim.get('sirasinda_n'), stim.get('sonra_n'),
          stim.get('once_hz'), stim.get('sirasinda_hz'),
          stim.get('sonra_hz'), stim.get('sonra_once_orani'),
          ttest.get('p_deger'), ttest.get('t_stat'),
          mw.get('p_deger'), mw.get('u_stat'),
          stim.get('n_var'), stim.get('n_yok')))
    conn.commit()


# ─────────────────────────────────────────────────────────────
# JSON / YAML — METADATA
# ─────────────────────────────────────────────────────────────

def metadata_yaz(organoid, kayit, metadata, klasor='metadata'):
    """
    Hiyerarşik / esnek metadata için JSON dosyası.
    Organoid yaşı, besiyeri, elektrot tipi, deney koşulları vb.

    Yapı: metadata/<organoid>/<kayit>.json
    """
    os.makedirs(os.path.join(klasor, organoid), exist_ok=True)
    yol = os.path.join(klasor, organoid, f'{kayit}.json')
    with open(yol, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2, ensure_ascii=False)
    return yol


def metadata_oku(organoid, kayit, klasor='metadata'):
    """Metadata JSON oku."""
    yol = os.path.join(klasor, organoid, f'{kayit}.json')
    if not os.path.exists(yol):
        return {}
    with open(yol, 'r', encoding='utf-8') as f:
        return json.load(f)


# ─────────────────────────────────────────────────────────────
# SORGULAR
# ─────────────────────────────────────────────────────────────

def kayit_ozet(conn, organoid=None):
    """Veritabanındaki kayıtların özet listesini döndürür."""
    if organoid:
        rows = conn.execute("""
            SELECT k.organoid, k.kayit, k.tarih_analiz,
                   k.kayit_suresi_sn, m.spike_sayisi, m.ortalama_hz, m.cv,
                   q.snr, q.refractory_ihlal_orani
            FROM kayitlar k
            LEFT JOIN metrikler m ON m.kayit_id = k.id
            LEFT JOIN qc_metrikleri q ON q.kayit_id = k.id
            WHERE k.organoid=?
            ORDER BY k.kayit
        """, (organoid,)).fetchall()
    else:
        rows = conn.execute("""
            SELECT k.organoid, k.kayit, k.tarih_analiz,
                   k.kayit_suresi_sn, m.spike_sayisi, m.ortalama_hz, m.cv,
                   q.snr, q.refractory_ihlal_orani
            FROM kayitlar k
            LEFT JOIN metrikler m ON m.kayit_id = k.id
            LEFT JOIN qc_metrikleri q ON q.kayit_id = k.id
            ORDER BY k.organoid, k.kayit
        """).fetchall()
    return [dict(r) for r in rows]


def spike_zaman_oku(conn, organoid, kayit, kanal_indeks=0, unit_id=-1):
    """
    Bir kaydın spike zamanlarını DB'den okur.
    BLOB'dan np.frombuffer ile direkt numpy array'e açılır.
    """
    row = conn.execute(
        "SELECT id FROM kayitlar WHERE organoid=? AND kayit=? AND kanal_indeks=?",
        (organoid, kayit, kanal_indeks)).fetchone()
    if row is None:
        return np.array([], dtype=np.float64)
    return spike_zamanlari_oku_blob(conn, row['id'], unit_id)
