"""
KURULUM.py — Axon v6.3
=======================
Yeni bilgisayarda bir kez çalıştır, hepsi hazır olur.
Sonrasında sadece organoid_cli.py kullanırsın.

Çalıştırma:
    python KURULUM.py
"""

import sys
import os
import subprocess
import platform

print("=" * 60)
print("  Axon v6.3 — Kurulum")
print("=" * 60)
print()

# ── 1. Python versiyonu ───────────────────────────────────
print("[1/5] Python versiyonu kontrol ediliyor...")
major = sys.version_info.major
minor = sys.version_info.minor
print(f"  Python {major}.{minor} bulundu.")

if major < 3 or (major == 3 and minor < 9):
    print()
    print("  HATA: Python 3.9 veya üzeri gerekli.")
    print("  https://www.python.org/downloads/ adresinden güncelle.")
    input("\nEnter'a bas çıkmak için...")
    sys.exit(1)

print("  OK\n")

# ── 2. pip ───────────────────────────────────────────────
print("[2/5] pip kontrol ediliyor...")
try:
    subprocess.run([sys.executable, "-m", "pip", "--version"],
                   check=True, capture_output=True)
    print("  OK\n")
except subprocess.CalledProcessError:
    print("  HATA: pip bulunamadı.")
    input("\nEnter'a bas çıkmak için...")
    sys.exit(1)

# ── 3. Temel kütüphaneler ────────────────────────────────
print("[3/5] Temel kütüphaneler kuruluyor...")
print("  (2-3 dakika sürebilir)\n")

temel = [
    ("dandi",        "dandi",        "DANDI Archive bağlantısı"),
    ("pynwb",        "pynwb",        "NWB dosya okuma"),
    ("h5py",         "h5py",         "HDF5 desteği"),
    ("numpy",        "numpy",        "Sayısal hesaplama"),
    ("scipy",        "scipy",        "İstatistik"),
    ("matplotlib",   "matplotlib",   "Grafik"),
    ("networkx",     "networkx",     "Graph theory"),
    ("remfile",      "remfile",      "DANDI streaming"),
    ("pandas",       "pandas",       "CSV işleme"),
]

basarisiz = []
for pip_adi, import_adi, aciklama in temel:
    try:
        __import__(import_adi)
        print(f"  ✓ {pip_adi:<15} zaten kurulu  ({aciklama})")
        continue
    except ImportError:
        pass
    print(f"  → {pip_adi:<15} kuruluyor...", end="", flush=True)
    r = subprocess.run(
        [sys.executable, "-m", "pip", "install", "--quiet", pip_adi],
        capture_output=True, text=True
    )
    if r.returncode == 0:
        print(" OK")
    else:
        print(" HATA")
        basarisiz.append(pip_adi)

print()

# ── 4. SpikeInterface ────────────────────────────────────
print("[4/5] SpikeInterface kuruluyor (spike sorting)...")
print("  (5-10 dakika sürebilir, büyük paket)\n")

try:
    import spikeinterface
    print(f"  ✓ spikeinterface {spikeinterface.__version__} zaten kurulu")
except ImportError:
    print("  → spikeinterface[full] kuruluyor...", end="", flush=True)
    r = subprocess.run(
        [sys.executable, "-m", "pip", "install", "--quiet", "spikeinterface[full]"],
        capture_output=True, text=True
    )
    if r.returncode == 0:
        print(" OK")
    else:
        print(" HATA")
        basarisiz.append("spikeinterface")

# Spike sorting algoritmaları
print("\n  Spike sorting algoritmaları kontrol ediliyor...")
algoritmalar = [
    ("mountainsort5", "mountainsort5", "MountainSort5 (önerilen)"),
    ("tridesclous2",  "tridesclous",   "Tridesclous2"),
    ("spykingcircus2","spyking_circus", "SpykingCircus2"),
]

for pip_adi, import_adi, aciklama in algoritmalar:
    try:
        __import__(import_adi)
        print(f"  ✓ {pip_adi:<20} kurulu  ({aciklama})")
    except ImportError:
        print(f"  → {pip_adi:<20} kuruluyor...", end="", flush=True)
        r = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--quiet", pip_adi],
            capture_output=True, text=True
        )
        print(" OK" if r.returncode == 0 else " HATA (opsiyonel)")

print()

# ── 5. Dosya yapısı ──────────────────────────────────────
print("[5/5] Dosya yapısı kontrol ediliyor...")

dosyalar = {
    "v6_3/":                          os.path.isdir("v6_3"),
    "v6_3/organoid_units_analiz.py":  os.path.isfile("v6_3/organoid_units_analiz.py"),
    "v6_3/organoid_io.py":            os.path.isfile("v6_3/organoid_io.py"),
    "v6_3/organoid_metrics.py":       os.path.isfile("v6_3/organoid_metrics.py"),
    "v6_3/organoid_signal.py":        os.path.isfile("v6_3/organoid_signal.py"),
    "organoid_cli.py":                os.path.isfile("organoid_cli.py"),
    "dandi_ara.py":                   os.path.isfile("dandi_ara.py"),
    "graph_analiz.py":                os.path.isfile("graph_analiz.py"),
    "gruplu_analiz.py":               os.path.isfile("gruplu_analiz.py"),
    "si_sorting.py":                  os.path.isfile("si_sorting.py"),
    "yas_metadata_cek.py":            os.path.isfile("yas_metadata_cek.py"),
}

eksik = []
for yol, var_mi in dosyalar.items():
    durum = "✓" if var_mi else "✗ EKSİK"
    print(f"  {durum}  {yol}")
    if not var_mi:
        eksik.append(yol)

print()

# ── Sonuç ─────────────────────────────────────────────────
print("=" * 60)

if not basarisiz and not eksik:
    print()
    print("  KURULUM TAMAMLANDI!")
    print()
    print("  Kullanım:")
    print("  python organoid_cli.py 001603 sub-HO2")
    print("  python dandi_ara.py organoid")
    print("  python graph_analiz.py")
    print()
    print("  Spike sorting (ham sinyal için):")
    print("  python organoid_cli.py 001132 sub-5C-1 --max-mb 750")
    print("  → Seçenek 2 (SpikeInterface)")
    print()
else:
    if basarisiz:
        print(f"\n  UYARI: Kurulamayan paketler: {', '.join(basarisiz)}")
    if eksik:
        print(f"\n  UYARI: Eksik dosyalar:")
        for e in eksik:
            print(f"    - {e}")
    print("\n  Yukarıdaki sorunları çöz, tekrar çalıştır.")

print("=" * 60)
input("\nEnter'a bas çıkmak için...")
