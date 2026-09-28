"""
yas_metadata_cek.py — Axon v6.3
=================================
dandi_indirilmis/ klasöründeki tüm NWB dosyalarının
subject metadata'sını okur. Yaş, tür, cinsiyet, subject ID
gibi bilgileri çeker ve CSV olarak kaydeder.

Kullanım:
    python yas_metadata_cek.py                         # dandi_indirilmis/
    python yas_metadata_cek.py --indir baska_klasor    # farklı klasör
    python yas_metadata_cek.py --cikti meta.csv        # farklı çıktı
"""

import os
import csv
import sys
import argparse
import pynwb
from pathlib import Path


def metadata_oku(dosya_yolu):
    """Tek NWB dosyasından metadata çeker."""
    kayit = {
        'dosya': Path(dosya_yolu).name,
        'subject_id': '',
        'species': '',
        'age': '',
        'age_reference': '',
        'sex': '',
        'description': '',
        'session_description': '',
        'institution': '',
        'lab': '',
        'recording_date': '',
        'hata': '',
    }

    try:
        with pynwb.NWBHDF5IO(str(dosya_yolu), 'r') as io:
            nwb = io.read()

            # Session bilgileri
            kayit['session_description'] = str(
                getattr(nwb, 'session_description', ''))[:200]
            kayit['institution'] = str(getattr(nwb, 'institution', ''))
            kayit['lab'] = str(getattr(nwb, 'lab', ''))

            # Tarih
            ts = getattr(nwb, 'session_start_time', None)
            if ts:
                kayit['recording_date'] = str(ts.date())

            # Subject
            subj = getattr(nwb, 'subject', None)
            if subj:
                kayit['subject_id'] = str(getattr(subj, 'subject_id', ''))
                kayit['species'] = str(getattr(subj, 'species', ''))
                kayit['age'] = str(getattr(subj, 'age', ''))
                kayit['age_reference'] = str(getattr(subj, 'age__reference', ''))
                kayit['sex'] = str(getattr(subj, 'sex', ''))
                kayit['description'] = str(
                    getattr(subj, 'description', ''))[:200]

            # Acquisition bilgisi
            acq_bilgi = []
            for name, obj in nwb.acquisition.items():
                tip = type(obj).__name__
                if hasattr(obj, 'data'):
                    try:
                        shape = obj.data.shape
                        sr = getattr(obj, 'rate', 'N/A')
                        acq_bilgi.append(f"{name}:{tip}:{shape}:{sr}Hz")
                    except Exception:
                        acq_bilgi.append(f"{name}:{tip}")
                else:
                    acq_bilgi.append(f"{name}:{tip}")
            kayit['acquisition'] = ' | '.join(acq_bilgi)

            # Units
            if nwb.units is not None:
                kayit['n_units'] = len(nwb.units)
                kayit['sorted_units'] = 'evet'
            else:
                kayit['n_units'] = 0
                kayit['sorted_units'] = 'hayır'

    except Exception as e:
        kayit['hata'] = str(e)[:200]

    return kayit


def main():
    parser = argparse.ArgumentParser(
        description='Axon v6.3 — NWB Metadata Okuyucu')
    parser.add_argument('--indir', type=str, default='dandi_indirilmis',
                        help='NWB dosyalarının bulunduğu klasör')
    parser.add_argument('--cikti', type=str,
                        default='sonuclar/metadata_ozet.csv',
                        help='Çıktı CSV dosyası')
    parser.add_argument('--verbose', action='store_true',
                        help='Her dosya için detaylı çıktı')
    args = parser.parse_args()

    indir = Path(args.indir)
    cikti = Path(args.cikti)

    print('=' * 60)
    print('  Axon v6.3 — NWB Metadata Okuyucu')
    print('=' * 60)
    print()

    if not indir.exists():
        print(f'HATA: {indir} klasörü bulunamadı.')
        print('dandi_indirilmis/ içinde NWB dosyaları olmalı.')
        return 1

    nwb_dosyalari = sorted(indir.glob('*.nwb'))
    if not nwb_dosyalari:
        print(f'  {indir} içinde NWB dosyası yok.')
        return 1

    print(f'  {len(nwb_dosyalari)} NWB dosyası bulundu.\n')

    sonuclar = []
    for dosya in nwb_dosyalari:
        print(f'  Okunuyor: {dosya.name}')
        kayit = metadata_oku(dosya)
        sonuclar.append(kayit)

        if args.verbose:
            print(f"    Subject    : {kayit['subject_id']}")
            print(f"    Species    : {kayit['species']}")
            print(f"    Age        : {kayit['age']} ({kayit['age_reference']})")
            print(f"    Sex        : {kayit['sex']}")
            print(f"    Units      : {kayit['n_units']} ({kayit['sorted_units']})")
            print(f"    Tarih      : {kayit['recording_date']}")
            if kayit['hata']:
                print(f"    HATA       : {kayit['hata']}")
            print()

    # CSV kaydet
    cikti.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ['dosya', 'subject_id', 'species', 'age', 'age_reference',
                  'sex', 'n_units', 'sorted_units', 'recording_date',
                  'institution', 'lab', 'session_description',
                  'acquisition', 'description', 'hata']

    with open(cikti, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sonuclar)

    print(f'\n  {len(sonuclar)} dosya işlendi.')
    print(f'  CSV kaydedildi: {cikti}')

    # Özet
    sorted_count = sum(1 for r in sonuclar if r['sorted_units'] == 'evet')
    print(f'\n  Özet:')
    print(f'    Sorted units olan: {sorted_count}/{len(sonuclar)}')

    species_list = [r['species'] for r in sonuclar if r['species']]
    if species_list:
        from collections import Counter
        for sp, n in Counter(species_list).most_common():
            print(f'    {sp}: {n} dosya')

    return 0


if __name__ == '__main__':
    sys.exit(main())
