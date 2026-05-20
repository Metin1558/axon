"""
yas_metadata_cek.py
====================
dandi_indirilmis/ klasöründeki tüm NWB dosyalarının
subject metadata'sını okur. Yaş, tür, iPSC hattı gibi
bilgileri çeker ve CSV olarak kaydeder.

Kullanım:
    python yas_metadata_cek.py
"""

import os
import csv
import pynwb
from pathlib import Path

INDIR_KLASOR = Path('dandi_indirilmis')
CIKTI = Path('sonuclar/metadata_ozet.csv')

print('=' * 55)
print('  NWB Metadata Okuyucu')
print('=' * 55)
print()

if not INDIR_KLASOR.exists():
    print(f'HATA: {INDIR_KLASOR} klasörü bulunamadı.')
    print('dandi_indirilmis/ içinde NWB dosyaları olmalı.')
    input('Enter\'a bas çıkmak için...')
    exit(1)

nwb_dosyalari = sorted(INDIR_KLASOR.glob('*.nwb'))
print(f'{len(nwb_dosyalari)} NWB dosyası bulundu.\n')

sonuclar = []

for dosya in nwb_dosyalari:
    print(f'  Okunuyor: {dosya.name}')

    kayit = {
        'dosya': dosya.name,
        'organoid': '',
        'identifier': '',
        'session_description': '',
        'subject_id': '',
        'subject_age': '',
        'subject_species': '',
        'subject_description': '',
        'subject_sex': '',
        'kayit_suresi_sn': '',
        'unit_sayisi': '',
    }

    try:
        with pynwb.NWBHDF5IO(str(dosya), 'r') as io:
            nwb = io.read()

            # Genel
            kayit['identifier'] = str(nwb.identifier or '')
            kayit['session_description'] = str(
                nwb.session_description or '')[:80]

            # Subject metadata
            if nwb.subject is not None:
                s = nwb.subject
                kayit['subject_id'] = str(getattr(s, 'subject_id', '') or '')
                kayit['subject_age'] = str(getattr(s, 'age', '') or '')
                kayit['subject_species'] = str(
                    getattr(s, 'species', '') or '')
                kayit['subject_description'] = str(
                    getattr(s, 'description', '') or '')[:100]
                kayit['subject_sex'] = str(getattr(s, 'sex', '') or '')

            # Organoid adını subject_id veya identifier'dan çıkar
            # "HO2_7M_sorted" -> HO2, 7M
            import re
            for pattern in [r'(HO\d+)', r'(MO\d+)']:
                m = re.search(pattern, kayit['identifier'] + ' ' +
                              kayit['subject_id'])
                if m:
                    kayit['organoid'] = m.group(1)
                    break

            # Kayıt süresi
            if nwb.units is not None and len(nwb.units) > 0:
                kayit['unit_sayisi'] = str(len(nwb.units))
                # Tüm spike zamanlarından max değer = kayıt süresi tahmini
                import numpy as np
                maks = 0
                for i in range(min(10, len(nwb.units))):
                    sp = nwb.units['spike_times'][i]
                    if len(sp) > 0:
                        maks = max(maks, float(np.max(sp)))
                kayit['kayit_suresi_sn'] = f'{maks:.1f}'

        print(f'    Organoid: {kayit["organoid"]}, '
              f'Yaş: {kayit["subject_age"]}, '
              f'Tür: {kayit["subject_species"]}, '
              f'Unit: {kayit["unit_sayisi"]}')

    except Exception as e:
        print(f'    HATA: {e}')

    sonuclar.append(kayit)

# CSV kaydet
CIKTI.parent.mkdir(exist_ok=True)
if sonuclar:
    with open(CIKTI, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=sonuclar[0].keys())
        writer.writeheader()
        writer.writerows(sonuclar)
    print(f'\n  CSV kaydedildi: {CIKTI}')

# Ekrana da yaz
print()
print('=' * 55)
print('  ÖZET')
print('=' * 55)
print(f'{"Organoid":<10} {"Yaş":<12} {"Tür":<20} {"Unit":<8} {"Süre(sn)"}')
print('-' * 65)
for r in sonuclar:
    print(f'{r["organoid"]:<10} {r["subject_age"]:<12} '
          f'{r["subject_species"][:18]:<20} {r["unit_sayisi"]:<8} '
          f'{r["kayit_suresi_sn"]}')

print()
print('Bu CSV\'yi Claude\'a yapıştır veya yükle.')
input('\nEnter\'a bas çıkmak için...')
