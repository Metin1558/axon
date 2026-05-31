"""
batch_analiz.py — Axon v6.3
============================
Bir dataset'ten birden fazla subject'i sırayla indirir ve analiz eder.
Sorted units varsa direkt analiz, yoksa SpikeInterface ile sorting.

Kullanım:
    # 5 subject otomatik seç ve analiz et
    python batch_analiz.py --dataset 001611 --n 5 --max-mb 30

    # Belirli subject'leri analiz et
    python batch_analiz.py --dataset 001132 --subjects sub-3C-1 sub-5C-1 sub-7D-7

    # Ham sinyal için sorting modu (SpikeInterface)
    python batch_analiz.py --dataset 001132 --n 5 --max-mb 2000 --sorting 1

Parametreler:
    --dataset   : DANDI dataset ID (örn: 001611)
    --n         : Kaç subject analiz edilecek (default: 5)
    --max-mb    : Maksimum dosya boyutu MB (default: 100)
    --subjects  : Manuel subject listesi (opsiyonel)
    --sorting   : Ham sinyal modu (0=hayır, 1=threshold, 2=SpikeInterface)
    --n-kanal   : SpikeInterface kanal sayısı (default: 64)
    --algoritma : Sorting algoritması (default: mountainsort5)
    --perm-n    : STTC permütasyon sayısı (default: 200)
"""

import sys
import os
import subprocess
import argparse
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / 'v6_3'))

try:
    from dandi.dandiapi import DandiAPIClient
except ImportError:
    print("HATA: dandi kurulu değil. pip install dandi")
    sys.exit(1)


def dataset_subject_listesi(dataset_id, max_mb=100, n=5, secili_subjects=None):
    """Dataset'ten uygun subject listesini çeker."""
    print(f"  DANDI:{dataset_id} taranıyor...")
    
    with DandiAPIClient() as client:
        dandiset = client.get_dandiset(dataset_id)
        assets = list(dandiset.get_assets())
    
    # NWB filtrele
    nwb_assets = [a for a in assets if a.path.endswith('.nwb')]
    
    # Subject -> dosya mapping
    subject_dosyalar = {}
    for asset in nwb_assets:
        parts = asset.path.split('/')
        subj = parts[0] if len(parts) > 1 else 'unknown'
        size_mb = asset.size / (1024 * 1024)
        
        if size_mb > max_mb:
            continue
        
        if subj not in subject_dosyalar:
            subject_dosyalar[subj] = []
        subject_dosyalar[subj].append((asset.path, size_mb))
    
    if secili_subjects:
        # Manuel seçim
        secili = []
        for subj in secili_subjects:
            if subj in subject_dosyalar:
                dosyalar = sorted(subject_dosyalar[subj], key=lambda x: x[1])
                secili.append((subj, dosyalar[0][0], dosyalar[0][1]))
            else:
                print(f"  UYARI: {subj} bulunamadı.")
        return secili
    else:
        # Otomatik seç — her subject'ten en küçük dosya
        secili = []
        for subj, dosyalar in sorted(subject_dosyalar.items()):
            dosyalar_sirali = sorted(dosyalar, key=lambda x: x[1])
            secili.append((subj, dosyalar_sirali[0][0], dosyalar_sirali[0][1]))
        
        return secili[:n]


def analiz_et(dataset_id, subject, dosya_yolu, max_mb, sorting_mode,
              n_kanal, algoritma, perm_n):
    """Tek subject'i indirir ve analiz eder."""
    print(f"\n{'='*60}")
    print(f"  Subject: {subject}")
    print(f"  Dosya  : {dosya_yolu}")
    print(f"  Boyut  : {max_mb:.1f} MB")
    print(f"{'='*60}")
    
    # organoid_cli.py çağır
    cmd = [
        sys.executable, 'organoid_cli.py',
        dataset_id, subject,
        '--max-mb', str(int(max_mb) + 50),
    ]
    
    if sorting_mode == 0:
        # Sorted units — interaktif input yok
        env = os.environ.copy()
        proc = subprocess.run(cmd, capture_output=False)
        return proc.returncode == 0
    else:
        # Ham sinyal — interaktif input simüle et
        # sorting_mode: 1=threshold, 2=SpikeInterface
        if sorting_mode == 2:
            # SpikeInterface: "2\n{n_kanal}\n1\n" (2=SI, kanal, 1=mountainsort5)
            alg_idx = {
                'mountainsort5': '1',
                'spykingcircus2': '2',
                'tridesclous2': '3',
                'simple': '4',
            }.get(algoritma, '1')
            input_str = f"2\n{n_kanal}\n{alg_idx}\n"
        else:
            # Threshold: "1\n"
            input_str = "1\n"
        
        proc = subprocess.run(
            cmd,
            input=input_str,
            text=True,
            capture_output=False,
        )
        return proc.returncode == 0


def main():
    parser = argparse.ArgumentParser(
        description='Axon v6.3 — Batch Analiz',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    
    parser.add_argument('--dataset', required=True, help='DANDI dataset ID')
    parser.add_argument('--n', type=int, default=5, help='Subject sayısı')
    parser.add_argument('--max-mb', type=float, default=100, dest='max_mb')
    parser.add_argument('--subjects', nargs='+', default=None,
                        help='Manuel subject listesi')
    parser.add_argument('--sorting', type=int, default=0,
                        choices=[0, 1, 2],
                        help='0=sorted units, 1=threshold, 2=SpikeInterface')
    parser.add_argument('--n-kanal', type=int, default=64, dest='n_kanal')
    parser.add_argument('--algoritma', type=str, default='mountainsort5')
    parser.add_argument('--perm-n', type=int, default=200, dest='perm_n')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("  Axon v6.3 — Batch Analiz")
    print("=" * 60)
    print(f"  Dataset  : DANDI:{args.dataset}")
    print(f"  Max boyut: {args.max_mb} MB")
    print(f"  Sorting  : {['sorted units', 'threshold', 'SpikeInterface'][args.sorting]}")
    if args.sorting == 2:
        print(f"  Kanal    : {args.n_kanal}")
        print(f"  Algoritma: {args.algoritma}")
    print()
    
    # Subject listesi
    try:
        subjects = dataset_subject_listesi(
            args.dataset, args.max_mb, args.n, args.subjects)
    except Exception as e:
        print(f"HATA: Dataset okunamadı: {e}")
        return 1
    
    if not subjects:
        print("UYARI: Uygun subject bulunamadı.")
        return 1
    
    print(f"  {len(subjects)} subject analiz edilecek:")
    for subj, dosya, boyut in subjects:
        print(f"    {subj:<20} {boyut:.1f} MB")
    print()
    
    # Analiz
    basarili = 0
    basarisiz = []
    toplam_sure = 0
    
    for i, (subj, dosya, boyut) in enumerate(subjects, 1):
        print(f"\n[{i}/{len(subjects)}] {subj}")
        t0 = time.time()
        
        ok = analiz_et(
            args.dataset, subj, dosya, boyut,
            args.sorting, args.n_kanal, args.algoritma, args.perm_n
        )
        
        sure = time.time() - t0
        toplam_sure += sure
        
        if ok:
            basarili += 1
            print(f"  ✓ {subj} tamamlandı ({sure:.0f}s)")
        else:
            basarisiz.append(subj)
            print(f"  ✗ {subj} başarısız")
    
    # Özet
    print()
    print("=" * 60)
    print("  BATCH ANALİZ TAMAMLANDI")
    print("=" * 60)
    print(f"  Başarılı : {basarili}/{len(subjects)}")
    print(f"  Toplam süre: {toplam_sure:.0f}s ({toplam_sure/60:.1f} dakika)")
    if basarisiz:
        print(f"  Başarısız: {', '.join(basarisiz)}")
    print(f"  Sonuçlar : sonuclar/ klasöründe")
    
    return 0 if not basarisiz else 1


if __name__ == '__main__':
    sys.exit(main())
