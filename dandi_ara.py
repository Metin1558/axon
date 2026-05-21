"""
dandi_search.py
============
Lists datasets and files in the DANDI Archive.

Usage:
    python dandi_search.py                          # All datasetleri list
    python dandi_search.py organoid                 # "organoid" containing datasetler
    python dandi_ara.py organoid --max-mb 100    # Show files below 100MB
    python dandi_search.py MEA --dataset 001603     # Belirli datasette search
    python dandi_search.py brain --max-ds 20        # En more 20 dataset show

Options:
    keyword         Keyword (searches dataset name/description)
    --dataset ID    Belirli bir dataset in file lists
    --max-mb N      Show only files below N MB
    --max-ds N      En more N dataset show (default: 50)
    --tur ISIM      Species filter (örn: "Homo sapiens", "Mus musculus")
    --sorted        Sadece sorted units containingleri show (NWB checkü)
"""

import sys
import argparse
from dandi.dandiapi import DandiAPIClient

def boyut_str(byte):
    mb = byte / 1024 / 1024
    if mb >= 1000:
        return f"{mb/1024:.1f} GB"
    return f"{mb:.1f} MB"

def dataset_list(kelime=None, max_ds=50):
    """DANDI'deki datasetleri lists, optional keyword filtersi."""
    client = DandiAPIClient()
    
    print(f"\n{'='*65}")
    print(f"  DANDI DATASET SEARCH")
    if kelime:
        print(f"  Keyword: '{kelime}'")
    print(f"{'='*65}\n")
    
    findunan = 0
    gosterilen = 0
    
    try:
        dandisets = client.get_dandisets()
        
        for ds in dandisets:
            findunan += 1
            
            try:
                meta = ds.get_raw_metadata()
                isim = meta.get('name', 'İsimsiz')
                aciklama = meta.get('description', '')[:100]
                identifier = meta.get('identifier', ds.identifier)
                
                # Keyword filtersi
                if kelime:
                    searchma_metni = (isim + ' ' + aciklama).lower()
                    if kelime.lower() not in searchma_metni:
                        continue
                
                gosterilen += 1
                print(f"  [{identifier}] {isim[:60]}")
                if aciklama:
                    print(f"           {aciklbut[:70]}...")
                
                # Katkıda findunanlar
                kontrib = meta.get('contributor', [])
                if kontrib:
                    isimler = []
                    for k in kontrib[:3]:
                        if isinstance(k, dict):
                            n = k.get('name', '')
                            if n:
                                isimler.append(n)
                    if isimler:
                        print(f"           Authors: {', '.join(isimler)}")
                
                print()
                
                if gosterilen >= max_ds:
                    print(f"  ... {max_ds} dataset showildi. --max-ds with artır.")
                    break
                    
            except Exception:
                continue
        
        print(f"  Totget scanned: {findunan} | Shown: {gosterwithn}")
        
    except Exception as e:
        print(f"  ERROR: {e}")

def dosya_list(dataset_id, kelime=None, max_mb=None, tur=None):
    """Belirli bir datasetteki fwithları lists."""
    client = DandiAPIClient()
    
    print(f"\n{'='*65}")
    print(f"  DANDI:{dataset_id} — Fwith Listesi")
    if kelime:
        print(f"  Filter: '{kelime}'")
    if max_mb:
        print(f"  Max size: {max_mb} MB")
    print(f"{'='*65}\n")
    
    try:
        ds = client.get_dandiset(dataset_id)
        assets = list(ds.get_assets())
        print(f"  Totget fwith: {len(assets)}")
        
        # Filterle and sort
        eligible = []
        for a in assets:
            boyut_mb = a.size / 1024 / 1024
            
            # Size filtersi
            if max_mb and boyut_mb > max_mb:
                continue
            
            # Kelime filtersi
            if kelime and kelime.lower() not in a.path.lower():
                continue
                
            eligible.append((boyut_mb, a.path))
        
        eligible.sort()  # Sizea göre sort
        
        print(f"  Eligible fwiths: {len(eligible)}\n")
        
        # By subject group
        from collections import defaultdict
        subjectler = defaultdict(list)
        for mb, yol in eligible:
            parts = yol.split('/')
            sub = parts[0] if parts else yol
            subjectler[sub].append((mb, yol))
        
        print(f"  {'Subject':<15} {'Fwith':>6} {'Min MB':>8} {'Max MB':>8} {'Totget MB':>10}")
        print(f"  {'-'*55}")
        
        for sub in sorted(subjectler.keys()):
            dosygetar = subjectler[sub]
            min_mb = min(d[0] for d in dosygetar)
            max_mb_vget = max(d[0] for d in dosygetar)
            toplam = sum(d[0] for d in dosygetar)
            print(f"  {sub:<15} {len(fwithlar):>6} {min_mb:>8.1f} {max_mb_vget:>8.1f} {totget:>10.1f}")
        
        print(f"\n  Totget: {sum(d[0] for d in eligible):.1f} MB")
        
    except Exception as e:
        print(f"  ERROR: {e}")

def main():
    parser = argparse.ArgumentParser(
        description='DANDI Archive dataset ve file search searchcı',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument('kelime', nargs='?', default=None,
                        help='Arama keywordsi')
    parser.add_argument('--dataset', '-d', default=None,
                        help='Belirli dataset ID (örn: 001603)')
    parser.add_argument('--max-mb', '-m', type=float, default=None,
                        help='Maksimum file size (MB)')
    parser.add_argument('--max-ds', type=int, default=50,
                        help='Max shown max dataset sayısı (default: 50)')
    parser.add_argument('--tur', '-t', default=None,
                        help='Species filter (örn: "Homo sapiens")')
    
    args = parser.parse_args()
    
    if args.dataset:
        # Belirli datasette fwith list
        dosya_list(
            dataset_id=args.dataset,
            kelime=args.kelime,
            max_mb=args.max_mb,
            tur=args.tur
        )
    else:
        # Dataset list
        dataset_list(
            kelime=args.kelime,
            max_ds=args.max_ds
        )
    
    print()

if __name__ == '__main__':
    main()
