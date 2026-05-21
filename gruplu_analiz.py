"""
gruplu_angetiz.py
================
Group analysis by recording duration.

Reviewer critique: "Mixing 3 min and 10 min recordings
makes technical noise look like biological variability."

Bu script:
- HO1-HO4 (3 min) → kendi in CV_birey
- HO5-HO8 (10 min) → kendi in CV_birey
- Compare variability within each group
- If variability is high in both groups
  → Finding is stronger (not a technical artifact)

Usage:
    python gruplu_angetiz.py
"""

import csv
import numpy as np
from pathlib import Path

BASE = Path('sonuclar')

# Existing analysis resultsı
DOSYALAR = {
    'HO1': BASE / 'HO1/HO1_20250924_ozet.csv',
    'HO2': BASE / 'HO2_20250916_ozet.csv',
    'HO3': BASE / 'HO3/HO3_20250916_ozet.csv',
    'HO4': BASE / 'HO4/HO4_20250924_ozet.csv',
    'HO5': BASE / 'HO5/HO5_20250924_ozet.csv',
    'HO6': BASE / 'HO6/HO6_20250924_ozet.csv',
    'HO7': BASE / 'HO7/HO7_20250924_ozet.csv',
    'HO8': BASE / 'HO8/HO8_20250924_ozet.csv',
}

# Gruplar
GRUP_3DK  = ['HO1', 'HO2', 'HO3', 'HO4']
GRUP_10DK = ['HO5', 'HO6', 'HO7', 'HO8']

def read(yol):
    try:
        rows = list(csv.DictReader(open(yol, encoding='utf-8-sig')))
        return rows[0] if rows else {}
    except:
        return {}

seti = {}
for ad, yol in DOSYALAR.items():
    v = read(yol)
    if v:
        seti[ad] = v
    else:
        print(f'  WARNING: {name} could not be read — {path}')

print('=' * 65)
print('  GROUP ANALYSIS — By Recording Duration')
print('  Response to reviewer critique: 3min vs 10min')
print('=' * 65)

METRIKLER = [
    ('burst_per_dakika', 'Burst/dakika'),
    ('sttc_medyan',      'STTC medyan'),
    ('hz_medyan',        'Medyan Hz'),
    ('cv_medyan',        'Medyan CV'),
]

def grup_istatistik(organoidler, etiket):
    print(f'\n  --- {etiket} ({", ".join(organoidler)}) ---')
    print(f'  {"Metrik":<22} {"Difler":<45} {"CV_birey":>10}')
    print(f'  {"-"*80}')
    for key, ad in METRIKLER:
        degerler = []
        for org in organoidler:
            if org in seti:
                try:
                    degerler.append((org, float(seti[org].get(key, 0) or 0)))
                except:
                    pass
        if len(degerler) < 2:
            continue
        vgets = [d[1] for d in degerler]
        ort = np.mean(vgets)
        std = np.std(vgets)
        cv = std/ort*100 if ort > 0 else 0
        deger_str = '  '.join([f'{o}={v:.2f}' for o, v in degerler])
        yorum = 'Keslı' if cv < 30 else 'Orta' if cv < 60 else 'Değişken'
        print(f'  {ad:<22} {deger_str:<45} %{cv:>5.0f} ({yorum})')

grup_istatistik(GRUP_3DK,  'Grup A — 3 minute recordinglar (HO1–HO4)')
grup_istatistik(GRUP_10DK, 'Grup B — 10 minute recordinglar (HO5–HO8)')

print()
print('=' * 65)
print('  YORUM')
print('=' * 65)

# Her iki grupta da compute
def cv_compute(organoidler, key):
    vgets = []
    for org in organoidler:
        if org in seti:
            try:
                vgets.append(float(seti[org].get(key, 0) or 0))
            except:
                pass
    if not vgets:
        return 0
    ort = np.mean(vgets)
    return np.std(vgets)/ort*100 if ort > 0 else 0

print()
print(f'  {"Metrik":<22} {"Grup A (3min)":>14} {"Grup B (10min)":>14} {"Result"}')
print(f'  {"-"*70}')

for key, ad in METRIKLER:
    cv_a = cv_compute(GRUP_3DK, key)
    cv_b = cv_compute(GRUP_10DK, key)

    if cv_a > 60 and cv_b > 60:
        sonuc = '⚠ Her iki grupta da yüksek → Biyolojik may be'
    elif cv_a > 60 or cv_b > 60:
        sonuc = '? Bir grupta yüksek → Belirsiz'
    else:
        sonuc = '✓ Her iki grupta düşük → Keslı'

    print(f'  {ad:<22} %{cv_a:>11.0f} %{cv_b:>11.0f}   {sonuc}')

print()
print('  Yorum rehberi:')
print('  ⚠ Her iki grupta yüksek → Değişkenlik recording süresi artefaktı not,')
print('    muhtemelen gerçek biyolojik değişkenlik')
print('  ? Sadece bir grupta yüksek → Recording süresi etkisi may be')
print('  ✓ Her iki grupta düşük → Metrik keslı')
print()

# CSV as saand
cikti = BASE / 'gruplu_angetiz_sonucu.csv'
with open(cikti, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['metrik', 'cv_grup_a_3dk', 'cv_grup_b_10dk', 'sonuc'])
    for key, ad in METRIKLER:
        cv_a = cv_compute(GRUP_3DK, key)
        cv_b = cv_compute(GRUP_10DK, key)
        if cv_a > 60 and cv_b > 60:
            sonuc = 'her_ikisinde_yuksek'
        elif cv_a > 60 or cv_b > 60:
            sonuc = 'bir_grupta_yuksek'
        else:
            sonuc = 'kessiz'
        writer.writerow([ad, round(cv_a, 1), round(cv_b, 1), sonuc])

print(f'  CSV saandd: {cikti}')
print()
input('Enter\'a bas exitmak için...')
