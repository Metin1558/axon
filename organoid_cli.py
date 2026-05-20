"""
organoid_cli.py — Organoid Analiz Komut Satırı Aracı
======================================================
Kullanım:
    python organoid_cli.py <DANDISET_ID> <SUBJECT_PREFIX> [seçenekler]

Örnekler:
    python organoid_cli.py 001603 sub-HO2
    python organoid_cli.py 001603 sub-HO5 --max-mb 200
    python organoid_cli.py 001603 sub-MO14 --max-mb 50
    python organoid_cli.py 001603 sub-HO9 --kayit-sec 2
    python organoid_cli.py 001603 sub-HO2 --listele

Ne yapar:
    1. DANDI'de subject'i arar
    2. En küçük (veya seçilen) NWB dosyasını bulur
    3. Zaten indirildiyse atlar, yoksa indirir
    4. organoid_units_analiz.py ile per-unit analiz yapar
    5. sonuclar/<SUBJECT>/ altına yazar

Seçenekler:
    --max-mb N       Max dosya boyutu MB (varsayılan: 300)
    --listele        Sadece uygun dosyaları listele, indirme/analiz yapma
    --kayit-sec N    N. dosyayı seç (varsayılan: 1 = en küçük)
    --perm-n N       Permutasyon iterasyon sayısı (varsayılan: 200)
    --yeniden-indir  Zaten varsa bile yeniden indir
    --sadece-analiz  İndirme atla, sadece mevcut dosyayı analiz et
"""

import sys
import os
import argparse
import subprocess
import time
import re
import urllib.request
from pathlib import Path


# ─────────────────────────────────────────────────────────────
# Paket kontrolü
# ─────────────────────────────────────────────────────────────

def paket_kontrol():
    eksik = []
    for p in ['dandi', 'pynwb', 'h5py', 'numpy', 'scipy', 'matplotlib']:
        try:
            __import__(p)
        except ImportError:
            eksik.append(p)
    if eksik:
        print(f'Eksik paketler: {", ".join(eksik)}')
        print(f'Kurulum: pip install {" ".join(eksik)}')
        sys.exit(1)


# ─────────────────────────────────────────────────────────────
# DANDI arama
# ─────────────────────────────────────────────────────────────

def dandi_ara(dandiset_id, subject_prefix, max_mb=300, sadece_listele=False):
    """
    DANDI'de subject_prefix ile eşleşen NWB dosyalarını arar.
    Boyuta göre sıralar, listeler.
    """
    from dandi.dandiapi import DandiAPIClient

    print(f'  DANDI:{dandiset_id} taraniyor — "{subject_prefix}"...')

    adaylar = []
    with DandiAPIClient() as client:
        dandiset = client.get_dandiset(dandiset_id, 'draft')
        sayac = 0
        for asset in dandiset.get_assets():
            sayac += 1
            if sayac % 30 == 0:
                print(f'    {sayac} dosya tarandı...', end='\r')

            path = asset.path
            size_mb = asset.size / 1e6

            if not path.endswith('.nwb'):
                continue
            if subject_prefix not in path:
                continue
            if size_mb > max_mb:
                continue

            adaylar.append({
                'path': path,
                'size_mb': size_mb,
                'asset': asset,
            })

    print(f'  Tarama tamamlandı ({sayac} dosya).     ')

    if not adaylar:
        print(f'  SONUÇ YOK: "{subject_prefix}" için {max_mb} MB altında NWB bulunamadı.')
        print(f'  İpucu: --max-mb değerini artır veya subject adını kontrol et.')
        return []

    # Boyuta göre sırala
    adaylar.sort(key=lambda x: x['size_mb'])

    print(f'\n  {len(adaylar)} uygun dosya bulundu:')
    for i, a in enumerate(adaylar):
        print(f'    {i+1}. ({a["size_mb"]:.1f} MB) {a["path"]}')

    return adaylar


# ─────────────────────────────────────────────────────────────
# İndirme
# ─────────────────────────────────────────────────────────────

def dosya_indir(dandiset_id, asset_path, hedef_klasor, yeniden=False):
    """
    DANDI'den dosyayı indirir. Zaten varsa atlar.
    Dönen: Path (yerel dosya yolu) veya None
    """
    from dandi.dandiapi import DandiAPIClient

    hedef_klasor = Path(hedef_klasor)
    hedef_klasor.mkdir(exist_ok=True)

    dosya_adi = Path(asset_path).name
    yerel_yol = hedef_klasor / dosya_adi

    if yerel_yol.exists() and not yeniden:
        boyut_mb = yerel_yol.stat().st_size / 1e6
        print(f'  Zaten mevcut: {dosya_adi} ({boyut_mb:.1f} MB) — atlandı.')
        print(f'  (Yeniden indirmek için --yeniden-indir kullan)')
        return yerel_yol

    print(f'  İndiriliyor: {dosya_adi}')

    with DandiAPIClient() as client:
        asset = client.get_dandiset(dandiset_id, 'draft').get_asset_by_path(asset_path)
        s3_url = asset.get_content_url(follow_redirects=1, strip_query=True)

    print(f'  URL: {s3_url[:65]}...')

    def progress(blok, blok_boyut, toplam):
        if toplam > 0:
            yuzde = min(100, blok * blok_boyut * 100 / toplam)
            mb = blok * blok_boyut / 1e6
            print(f'    {yuzde:.0f}% — {mb:.1f}/{toplam/1e6:.1f} MB', end='\r')

    try:
        urllib.request.urlretrieve(s3_url, str(yerel_yol), reporthook=progress)
        boyut = yerel_yol.stat().st_size / 1e6
        print(f'  İndirme tamam: {boyut:.1f} MB        ')
        return yerel_yol
    except Exception as e:
        print(f'\n  HATA: {e}')
        if yerel_yol.exists():
            yerel_yol.unlink()  # Yarım dosyayı sil
        return None


# ─────────────────────────────────────────────────────────────
# Analiz
# ─────────────────────────────────────────────────────────────

def analiz_et(dosya_yolu, organoid_ad, v6_3_yol, sonuc_klasor, perm_n=200):
    """
    organoid_units_analiz.py çalıştırır.
    Dönen: (basarili: bool, sure_sn: float)
    """
    dosya_yolu = Path(dosya_yolu)
    v6_3_yol = Path(v6_3_yol)
    per_unit_script = v6_3_yol / 'organoid_units_analiz.py'

    if not per_unit_script.exists():
        print(f'  HATA: organoid_units_analiz.py bulunamadı — {per_unit_script}')
        return False, 0

    m = re.search(r'ses-(\d{8})', dosya_yolu.name)
    kayit = m.group(1) if m else 'kayit'

    cikti = Path(sonuc_klasor) / organoid_ad
    cikti.mkdir(parents=True, exist_ok=True)

    komut = [
        sys.executable, str(per_unit_script),
        str(dosya_yolu.resolve()),
        organoid_ad, kayit,
        '--out-dir', str(cikti.resolve()),
        '--perm-n', str(perm_n),
    ]

    print(f'\n  Analiz başlıyor...')
    print(f'  Komut: {" ".join(str(x) for x in komut[2:])}')
    print()

    t0 = time.time()
    r = subprocess.run(komut, cwd=str(v6_3_yol), text=True)
    sure = time.time() - t0

    basarili = r.returncode == 0
    if not basarili:
        print(f'  UYARI: exit code {r.returncode}')

    return basarili, sure


# ─────────────────────────────────────────────────────────────
# Ana akış
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='DANDI organoid arama, indirme ve analiz CLI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('dandiset_id',
                        help='DANDI dataset ID (orn: 001603)')
    parser.add_argument('subject',
                        help='Subject prefix (orn: sub-HO2, sub-MO14)')
    parser.add_argument('--max-mb', type=float, default=300,
                        help='Max dosya boyutu MB (varsayılan: 300)')
    parser.add_argument('--listele', action='store_true',
                        help='Sadece dosyaları listele')
    parser.add_argument('--kayit-sec', type=int, default=1,
                        help='Hangi dosya seçilsin? 1=en küçük (varsayılan)')
    parser.add_argument('--perm-n', type=int, default=200,
                        help='Permutasyon sayısı (varsayılan: 200)')
    parser.add_argument('--yeniden-indir', action='store_true',
                        help='Zaten varsa bile yeniden indir')
    parser.add_argument('--sadece-analiz', action='store_true',
                        help='İndirme atla, mevcut dosyayı analiz et')
    parser.add_argument('--indir-klasor', default='dandi_indirilmis',
                        help='İndirme klasörü (varsayılan: dandi_indirilmis)')
    parser.add_argument('--sonuc-klasor', default='sonuclar',
                        help='Sonuç klasörü (varsayılan: sonuclar)')

    args = parser.parse_args()

    # v6.3 yolunu belirle
    v6_3_yol = Path(__file__).parent / 'v6_3'
    if not v6_3_yol.exists():
        print(f'HATA: v6_3 klasörü bulunamadı: {v6_3_yol}')
        sys.exit(1)

    print('=' * 60)
    print(f'  ORGANOİD CLI')
    print(f'  Dataset: DANDI:{args.dandiset_id}')
    print(f'  Subject: {args.subject}')
    print(f'  Max boyut: {args.max_mb} MB')
    print('=' * 60)
    print()

    # Paket kontrolü
    paket_kontrol()

    # ── 1. DANDI'de ara ──────────────────────────────────────
    print('[1] Dosya arama...')

    if args.sadece_analiz:
        # Yerel klasörde ara
        yerel_klasor = Path(args.indir_klasor)
        prefix = args.subject.replace('sub-', '')
        mevcut = sorted(yerel_klasor.glob(f'*{prefix}*.nwb'))
        if not mevcut:
            print(f'  HATA: {yerel_klasor}/ içinde {prefix}*.nwb bulunamadı.')
            sys.exit(1)
        seçilen_yerel = mevcut[0]
        print(f'  Mevcut dosya: {seçilen_yerel.name}')
        organoid_ad = args.subject.replace('sub-', '')
        basarili, sure = analiz_et(
            seçilen_yerel, organoid_ad, v6_3_yol,
            args.sonuc_klasor, args.perm_n)
        sys.exit(0 if basarili else 1)

    adaylar = dandi_ara(
        args.dandiset_id, args.subject,
        max_mb=args.max_mb, sadece_listele=args.listele)

    if args.listele or not adaylar:
        sys.exit(0)

    # Dosya seç
    secim_idx = min(args.kayit_sec - 1, len(adaylar) - 1)
    secilen = adaylar[secim_idx]
    print(f'\n  Seçilen (#{secim_idx+1}): {secilen["path"]} '
          f'({secilen["size_mb"]:.1f} MB)')

    # ── 2. İndir ─────────────────────────────────────────────
    print(f'\n[2] İndirme...')
    yerel_dosya = dosya_indir(
        args.dandiset_id,
        secilen['path'],
        args.indir_klasor,
        yeniden=args.yeniden_indir)

    if yerel_dosya is None:
        print('  İndirme başarısız — çıkılıyor.')
        sys.exit(1)

    # ── 3. Analiz ────────────────────────────────────────────
    print(f'\n[3] Per-unit analiz...')
    organoid_ad = args.subject.replace('sub-', '')
    basarili, sure = analiz_et(
        yerel_dosya, organoid_ad, v6_3_yol,
        args.sonuc_klasor, args.perm_n)

    # ── 4. Özet ──────────────────────────────────────────────
    print()
    print('=' * 60)
    durum = 'TAMAMLANDI' if basarili else 'UYARILARLA BİTTİ'
    print(f'  {durum} ({sure:.1f} sn)')
    print(f'  Sonuçlar: {args.sonuc_klasor}/{organoid_ad}/')
    print('=' * 60)

    sys.exit(0 if basarili else 1)


if __name__ == '__main__':
    main()
