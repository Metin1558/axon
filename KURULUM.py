"""
KURULUM.py
==========
Yeni bilgisayarda bir kez çalıştır, hepsi hazır olur.
Sonrasında sadece organoid_cli.py kullanırsın.

Çalıştırma:
    python KURULUM.py
"""

import sys
import os
import subprocess
import platform

print("=" * 55)
print("  organoid v6.3 — Yeni Bilgisayar Kurulumu")
print("=" * 55)
print()

# ── 1. Python versiyonu kontrol ──────────────────────────
print("[1/4] Python versiyonu kontrol ediliyor...")
major = sys.version_info.major
minor = sys.version_info.minor
print(f"  Python {major}.{minor} bulundu.")

if major < 3 or (major == 3 and minor < 8):
    print()
    print("  HATA: Python 3.8 veya üzeri gerekli.")
    print("  Şu an yüklü versiyon çok eski.")
    print()
    print("  Çözüm:")
    print("  1. https://www.python.org/downloads/ adresine git")
    print("  2. En güncel Python 3.x sürümünü indir ve kur")
    print("  3. Kurulumda 'Add Python to PATH' kutusunu işaretle")
    print("  4. Bu dosyayı tekrar çalıştır")
    input("\nEnter'a bas çıkmak için...")
    sys.exit(1)

print("  OK — versiyon uyumlu.\n")

# ── 2. pip kontrol ───────────────────────────────────────
print("[2/4] pip kontrol ediliyor...")
try:
    subprocess.run([sys.executable, "-m", "pip", "--version"],
                   check=True, capture_output=True)
    print("  OK — pip mevcut.\n")
except subprocess.CalledProcessError:
    print("  HATA: pip bulunamadı.")
    print("  Çözüm: Python'u yeniden kur, 'Add to PATH' seç.")
    input("\nEnter'a bas çıkmak için...")
    sys.exit(1)

# ── 3. Gerekli kütüphaneleri kur ─────────────────────────
print("[3/4] Gerekli kütüphaneler kuruluyor...")
print("  (İnternet bağlantısı gerekli, 2-5 dakika sürebilir)\n")

kutuphaneler = [
    # paket adı      # import adı    # açıklama
    ("dandi",         "dandi",        "DANDI Archive bağlantısı"),
    ("pynwb",         "pynwb",        "NWB dosya okuma"),
    ("h5py",          "h5py",         "HDF5 dosya desteği"),
    ("numpy",         "numpy",        "Sayısal hesaplama"),
    ("scipy",         "scipy",        "İstatistik fonksiyonları"),
    ("matplotlib",    "matplotlib",   "Grafik üretimi"),
    ("remfile",       "remfile",      "DANDI streaming"),
    ("pypdf",         "pypdf",        "PDF okuma (opsiyonel)"),
]

basarisiz = []

for pip_adi, import_adi, aciklama in kutuphaneler:
    # Önce zaten kurulu mu kontrol et
    try:
        __import__(import_adi)
        print(f"  ✓ {pip_adi:<15} zaten kurulu  ({aciklama})")
        continue
    except ImportError:
        pass

    # Kurulu değil, kur
    print(f"  → {pip_adi:<15} kuruluyor...  ({aciklama})", end="", flush=True)
    sonuc = subprocess.run(
        [sys.executable, "-m", "pip", "install", "--quiet", pip_adi],
        capture_output=True, text=True
    )

    if sonuc.returncode == 0:
        print(" TAMAM")
    else:
        print(" HATA")
        basarisiz.append(pip_adi)

print()

if basarisiz:
    print(f"  UYARI: Şu paketler kurulamadı: {', '.join(basarisiz)}")
    print("  İnternet bağlantını kontrol et ve tekrar dene.")
    print()
else:
    print("  Tüm kütüphaneler hazır.\n")

# ── 4. v6_3 klasörü kontrol ──────────────────────────────
print("[4/4] Klasör yapısı kontrol ediliyor...")

klasor_durumu = {
    "v6_3/":                       os.path.isdir("v6_3"),
    "v6_3/organoid.py":            os.path.isfile("v6_3/organoid.py"),
    "v6_3/organoid_units_analiz.py": os.path.isfile("v6_3/organoid_units_analiz.py"),
    "v6_3/organoid_io.py":         os.path.isfile("v6_3/organoid_io.py"),
    "v6_3/organoid_metrics.py":    os.path.isfile("v6_3/organoid_metrics.py"),
    "organoid_cli.py":             os.path.isfile("organoid_cli.py"),
}

eksik = []
for yol, var_mi in klasor_durumu.items():
    durum = "✓" if var_mi else "✗ EKSİK"
    print(f"  {durum}  {yol}")
    if not var_mi:
        eksik.append(yol)

print()

if eksik:
    print("  UYARI: Bazı dosyalar eksik:")
    for e in eksik:
        print(f"    - {e}")
    print()
    print("  Bu dosyalar calistirma_klasoru içinde olmalı.")
    print("  Flash'tan veya kaynaktan eksik dosyaları kopyala.")
else:
    print("  Tüm dosyalar yerli yerinde.\n")

# ── Sonuç ─────────────────────────────────────────────────
print("=" * 55)

if not basarisiz and not eksik:
    print()
    print("  KURULUM TAMAMLANDI!")
    print()
    print("  Artık organoid_cli.py kullanabilirsin.")
    print()
    print("  Örnek ilk komut:")
    print("  python organoid_cli.py 001603 sub-HO2 --listele")
    print()
else:
    print()
    print("  Kurulum eksik veya hatalı bitti.")
    print("  Yukarıdaki UYARI mesajlarını çöz, tekrar çalıştır.")
    print()

print("=" * 55)
input("\nEnter'a bas çıkmak için...")
