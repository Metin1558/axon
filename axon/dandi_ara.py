"""
dandi_ara.py
============
DANDI Archive'deki datasetleri ve dosyaları listeler.

Kullanım:
    python dandi_ara.py                          # Tüm datasetleri listele
    python dandi_ara.py organoid                 # "organoid" içeren datasetler
    python dandi_ara.py organoid --max-mb 100    # 100MB altı dosyaları göster
    python dandi_ara.py MEA --dataset 001603     # Belirli datasette ara
    python dandi_ara.py brain --max-ds 20        # En fazla 20 dataset göster

Seçenekler:
    kelime          Anahtar kelime (dataset adı/açıklamada arar)
    --dataset ID    Belirli bir dataset içinde dosya listeler
    --max-mb N      Sadece N MB altı dosyaları göster
    --max-ds N      En fazla N dataset göster (varsayılan: 50)
    --tur ISIM      Tür filtresi (örn: "Homo sapiens", "Mus musculus")
    --sorted        Sadece sorted units içerenleri göster (NWB kontrolü)
"""

import sys
import argparse
from dandi.dandiapi import DandiAPIClient

def boyut_str(byte):
    mb = byte / 1024 / 1024
    if mb >= 1000:
        return f"{mb/1024:.1f} GB"
    return f"{mb:.1f} MB"

def dataset_listele(kelime=None, max_ds=50):
    """DANDI'deki datasetleri listeler, opsiyonel anahtar kelime filtresi."""
    client = DandiAPIClient()
    
    print(f"\n{'='*65}")
    print(f"  DANDI DATASET ARAMA")
    if kelime:
        print(f"  Anahtar kelime: '{kelime}'")
    print(f"{'='*65}\n")
    
    bulunan = 0
    gosterilen = 0
    
    try:
        dandisets = client.get_dandisets()
        
        for ds in dandisets:
            bulunan += 1
            
            try:
                meta = ds.get_raw_metadata()
                isim = meta.get('name', 'İsimsiz')
                aciklama = meta.get('description', '')[:100]
                identifier = meta.get('identifier', ds.identifier)
                
                # Anahtar kelime filtresi
                if kelime:
                    arama_metni = (isim + ' ' + aciklama).lower()
                    if kelime.lower() not in arama_metni:
                        continue
                
                gosterilen += 1
                print(f"  [{identifier}] {isim[:60]}")
                if aciklama:
                    print(f"           {aciklama[:70]}...")
                
                # Katkıda bulunanlar
                kontrib = meta.get('contributor', [])
                if kontrib:
                    isimler = []
                    for k in kontrib[:3]:
                        if isinstance(k, dict):
                            n = k.get('name', '')
                            if n:
                                isimler.append(n)
                    if isimler:
                        print(f"           Yazarlar: {', '.join(isimler)}")
                
                print()
                
                if gosterilen >= max_ds:
                    print(f"  ... {max_ds} dataset gösterildi. --max-ds ile artır.")
                    break
                    
            except Exception:
                continue
        
        print(f"  Toplam taranan: {bulunan} | Gösterilen: {gosterilen}")
        
    except Exception as e:
        print(f"  HATA: {e}")

def dosya_listele(dataset_id, kelime=None, max_mb=None, tur=None):
    """Belirli bir datasetteki dosyaları listeler."""
    client = DandiAPIClient()
    
    print(f"\n{'='*65}")
    print(f"  DANDI:{dataset_id} — Dosya Listesi")
    if kelime:
        print(f"  Filtre: '{kelime}'")
    if max_mb:
        print(f"  Max boyut: {max_mb} MB")
    print(f"{'='*65}\n")
    
    try:
        ds = client.get_dandiset(dataset_id)
        assets = list(ds.get_assets())
        print(f"  Toplam dosya: {len(assets)}")
        
        # Filtrele ve sırala
        uygun = []
        for a in assets:
            boyut_mb = a.size / 1024 / 1024
            
            # Boyut filtresi
            if max_mb and boyut_mb > max_mb:
                continue
            
            # Kelime filtresi
            if kelime and kelime.lower() not in a.path.lower():
                continue
                
            uygun.append((boyut_mb, a.path))
        
        uygun.sort()  # Boyuta göre sırala
        
        print(f"  Uygun dosya: {len(uygun)}\n")
        
        # Subject bazında grupla
        from collections import defaultdict
        subjectler = defaultdict(list)
        for mb, yol in uygun:
            parts = yol.split('/')
            sub = parts[0] if parts else yol
            subjectler[sub].append((mb, yol))
        
        print(f"  {'Subject':<15} {'Dosya':>6} {'Min MB':>8} {'Max MB':>8} {'Toplam MB':>10}")
        print(f"  {'-'*55}")
        
        for sub in sorted(subjectler.keys()):
            dosyalar = subjectler[sub]
            min_mb = min(d[0] for d in dosyalar)
            max_mb_val = max(d[0] for d in dosyalar)
            toplam = sum(d[0] for d in dosyalar)
            print(f"  {sub:<15} {len(dosyalar):>6} {min_mb:>8.1f} {max_mb_val:>8.1f} {toplam:>10.1f}")
        
        print(f"\n  Toplam: {sum(d[0] for d in uygun):.1f} MB")
        
    except Exception as e:
        print(f"  HATA: {e}")

def main():
    parser = argparse.ArgumentParser(
        description='DANDI Archive dataset ve dosya arama aracı',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('kelime', nargs='?', default=None,
                        help='Arama anahtar kelimesi')
    parser.add_argument('--dataset', '-d', default=None,
                        help='Belirli dataset ID (örn: 001603)')
    parser.add_argument('--max-mb', '-m', type=float, default=None,
                        help='Maksimum dosya boyutu (MB)')
    parser.add_argument('--max-ds', type=int, default=50,
                        help='Gösterilecek max dataset sayısı (varsayılan: 50)')
    parser.add_argument('--tur', '-t', default=None,
                        help='Tür filtresi (örn: "Homo sapiens")')
    
    args = parser.parse_args()
    
    if args.dataset:
        # Belirli datasette dosya listele
        dosya_listele(
            dataset_id=args.dataset,
            kelime=args.kelime,
            max_mb=args.max_mb,
            tur=args.tur
        )
    else:
        # Dataset listele
        dataset_listele(
            kelime=args.kelime,
            max_ds=args.max_ds
        )
    
    print()

if __name__ == '__main__':
    main()
