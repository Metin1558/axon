"""
gruplu_analiz.py
================
Kayıt sürelerine göre gruplu analiz.

Hakem eleştirisi: "3 dk ve 10 dk kayıtları karıştırmak
teknik gürültüyü biyolojik değişkenlik gibi gösterir."

Bu script:
- HO1-HO4 (3 dk) → kendi içinde CV_birey
- HO5-HO8 (10 dk) → kendi içinde CV_birey
- İki gruptaki değişkenliği karşılaştır
- Eğer her iki grupta da değişkenlik yüksekse
  → Bulgu daha güçlü (teknik artefakt değil)

Kullanım:
    python gruplu_analiz.py
"""

import csv
import numpy as np
from pathlib import Path

BASE = Path('sonuclar')

# Mevcut analiz sonuçları
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

def oku(yol):
    try:
        rows = list(csv.DictReader(open(yol, encoding='utf-8-sig')))
        return rows[0] if rows else {}
    except:
        return {}

veri = {}
for ad, yol in DOSYALAR.items():
    v = oku(yol)
    if v:
        veri[ad] = v
    else:
        print(f'  UYARI: {ad} okunamadı — {yol}')

print('=' * 65)
print('  GRUPLU ANALİZ — Kayıt Süresine Göre')
print('  Hakem eleştirisine yanıt: 3dk vs 10dk')
print('=' * 65)

METRIKLER = [
    ('burst_per_dakika', 'Burst/dakika'),
    ('sttc_medyan',      'STTC medyan'),
    ('hz_medyan',        'Medyan Hz'),
    ('cv_medyan',        'Medyan CV'),
]

def grup_istatistik(organoidler, etiket):
    print(f'\n  --- {etiket} ({", ".join(organoidler)}) ---')
    print(f'  {"Metrik":<22} {"Değerler":<45} {"CV_birey":>10}')
    print(f'  {"-"*80}')
    for key, ad in METRIKLER:
        degerler = []
        for org in organoidler:
            if org in veri:
                try:
                    degerler.append((org, float(veri[org].get(key, 0) or 0)))
                except:
                    pass
        if len(degerler) < 2:
            continue
        vals = [d[1] for d in degerler]
        ort = np.mean(vals)
        std = np.std(vals)
        cv = std/ort*100 if ort > 0 else 0
        deger_str = '  '.join([f'{o}={v:.2f}' for o, v in degerler])
        yorum = 'Kararlı' if cv < 30 else 'Orta' if cv < 60 else 'Değişken'
        print(f'  {ad:<22} {deger_str:<45} %{cv:>5.0f} ({yorum})')

grup_istatistik(GRUP_3DK,  'Grup A — 3 dakika kayıtlar (HO1–HO4)')
grup_istatistik(GRUP_10DK, 'Grup B — 10 dakika kayıtlar (HO5–HO8)')

print()
print('=' * 65)
print('  YORUM')
print('=' * 65)

# Her iki grupta da hesapla
def cv_hesapla(organoidler, key):
    vals = []
    for org in organoidler:
        if org in veri:
            try:
                vals.append(float(veri[org].get(key, 0) or 0))
            except:
                pass
    if not vals:
        return 0
    ort = np.mean(vals)
    return np.std(vals)/ort*100 if ort > 0 else 0

print()
print(f'  {"Metrik":<22} {"Grup A (3dk)":>14} {"Grup B (10dk)":>14} {"Sonuç"}')
print(f'  {"-"*70}')

for key, ad in METRIKLER:
    cv_a = cv_hesapla(GRUP_3DK, key)
    cv_b = cv_hesapla(GRUP_10DK, key)

    if cv_a > 60 and cv_b > 60:
        sonuc = '⚠ Her iki grupta da yüksek → Biyolojik olabilir'
    elif cv_a > 60 or cv_b > 60:
        sonuc = '? Bir grupta yüksek → Belirsiz'
    else:
        sonuc = '✓ Her iki grupta düşük → Kararlı'

    print(f'  {ad:<22} %{cv_a:>11.0f} %{cv_b:>11.0f}   {sonuc}')

print()
print('  Yorum rehberi:')
print('  ⚠ Her iki grupta yüksek → Değişkenlik kayıt süresi artefaktı değil,')
print('    muhtemelen gerçek biyolojik değişkenlik')
print('  ? Sadece bir grupta yüksek → Kayıt süresi etkisi olabilir')
print('  ✓ Her iki grupta düşük → Metrik kararlı')
print()

# CSV olarak kaydet
cikti = BASE / 'gruplu_analiz_sonucu.csv'
with open(cikti, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['metrik', 'cv_grup_a_3dk', 'cv_grup_b_10dk', 'sonuc'])
    for key, ad in METRIKLER:
        cv_a = cv_hesapla(GRUP_3DK, key)
        cv_b = cv_hesapla(GRUP_10DK, key)
        if cv_a > 60 and cv_b > 60:
            sonuc = 'her_ikisinde_yuksek'
        elif cv_a > 60 or cv_b > 60:
            sonuc = 'bir_grupta_yuksek'
        else:
            sonuc = 'kararsiz'
        writer.writerow([ad, round(cv_a, 1), round(cv_b, 1), sonuc])

print(f'  CSV kaydedildi: {cikti}')
print()
input('Enter\'a bas çıkmak için...')
