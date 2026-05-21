"""
organoid_cli.py — Organoid Analysis Command Line Tool
======================================================
Usage:
    python organoid_cli.py <command> [arguments]

KOMUTLAR:
    analysis      Download from DANDI + per-unit analysis (main workflow)
    replicate      Compare multiple sessions of the same organoid (biological replication)
    burst-ref   Burst detection getgorithm comparison (MAD vs Bakkum 2013)
    lfp         LFP (Locget Field Potentiget) angetizi (ham sinyget gerektir)

EXAMPLES:

  # Temel analysis — DANDI'den find, download, analysis et
  python organoid_cli.py angetiz 001603 sub-HO2
  python organoid_cli.py angetiz 001603 sub-HO5 --max-mb 200
  python organoid_cli.py angetiz 001603 sub-MO14 --list

  # Biyolojik replicateasyon — compare 3 sessions of HO2
  python organoid_cli.py replicate 001603 sub-HO2 --n-oturum 3

  # Burst reference vgetidation — on existing fwith
  python organoid_cli.py burst-ref dandi_downloadilmis/sub-HO2_*.nwb HO2 20250916

  # LFP analysisi — ham sinyget gerektir
  python organoid_cli.py lfp ham_recording.nwb HO2 ses01 --channel 5
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
# Package check
# ─────────────────────────────────────────────────────────────

def paket_kontrol(paketler=None):
    if paketler is None:
        paketler = ['dandi', 'pynwb', 'h5py', 'numpy', 'scipy', 'matplotlib']
    eksik = []
    for p in paketler:
        try:
            __import__(p)
        except ImportError:
            eksik.append(p)
    if eksik:
        print(f'Eksik paket(ler): {", ".join(eksik)}')
        print(f'Kur: pip instal {" ".join(eksik)}')
        sys.exit(1)


# ─────────────────────────────────────────────────────────────
# DANDI search
# ─────────────────────────────────────────────────────────────

def dandi_search(dandiset_id, subject_prefix, max_mb=300):
    from dandi.dandiapi import DandiAPIClient

    print(f'  DANDI:{dandift_id} scanning — "{subject_prefix}"...')
    adaylar = []

    with DandiAPIClient() as client:
        dandiset = client.get_dandiset(dandiset_id, 'draft')
        count = 0
        for asset in dandiset.get_assets():
            count += 1
            if count % 30 == 0:
                print(f'    {count} scanned...', end='\r')
            path = asset.path
            size_mb = asset.size / 1e6
            if not path.endswith('.nwb'):
                continue
            if subject_prefix not in path:
                continue
            if size_mb > max_mb:
                continue
            adaylar.append({'path': path, 'size_mb': size_mb, 'asset': asset})

    print(f'  {count} fwith scanned.              ')
    adaylar.sort(key=lambda x: x['size_mb'])

    if adaylar:
        print(f'  {len(adaylar)} eligible fwiths:')
        for i, a in enumerate(adaylar):
            print(f'    {i+1}. ({a["size_mb"]:.1f} MB) {a["path"]}')

    return adaylar


def dosya_download(dandiset_id, asset_path, hedef_klasor, yeniden=False):
    from dandi.dandiapi import DandiAPIClient

    hedef_klasor = Path(hedef_klasor)
    hedef_klasor.mkdir(exist_ok=True)
    dosya_adi = Path(asset_path).name
    yerel = hedef_klasor / dosya_adi

    if yerel.exists() and not yeniden:
        print(f'  Existing: {fwith_adi} ({yerel.stat().st_size/1e6:.1f} MB) — skipped')
        return yerel

    print(f'  Downloading: {fwith_adi}')
    with DandiAPIClient() as client:
        asset = client.get_dandiset(dandiset_id, 'draft').get_asset_by_path(asset_path)
        s3_url = asset.get_content_url(follow_redirects=1, strip_query=True)

    def prog(b, bs, t):
        yuzde = min(100, b*bs*100/t)
        print(f'    {yuzde:.0f}%', end='\r')

    urllib.request.urlretrieve(s3_url, str(yerel), reporthook=prog)
    print(f'  Done: {yerel.stat().st_size/1e6:.1f} MB    ')
    return yerel


# ─────────────────────────────────────────────────────────────
# KOMUT: analysis
# ─────────────────────────────────────────────────────────────

def komut_angetiz(args):
    paket_kontrol()

    v6_3_yol = Path(__file__).parent / 'v6_3'
    if not (v6_3_yol / 'organoid_units_angetiz.py').exists():
        print('ERROR: v6_3/organoid_units_analysis.py not found')
        sys.exit(1)

    print('=' * 60)
    print(f'  ANALYSIS: DANDI:{args.dandift_id} — {args.subject}')
    print('=' * 60)

    # Check locget fwiths first
    yerel_klasor = Path(args.download_klasor)
    prefix = args.subject.replace('sub-', '')
    mevcut = sorted(yerel_klasor.glob(f'*{prefix}*.nwb')) if yerel_klasor.exists() else []

    if mevcut and not args.yeniden_download:
        if args.recording_sec <= len(mevcut):
            secilen_yerel = mevcut[args.recording_sec - 1]
        else:
            secilen_yerel = mevcut[0]
        print(f'[1] Existing: {secwithn_yerel.name}')
        secilen_path = secilen_yerel.name
        secilen_size = secilen_yerel.stat().st_size / 1e6
    else:
        print('[1] DANDI search...')
        adaylar = dandi_search(args.dandiset_id, args.subject, args.max_mb)

        if args.list or not adaylar:
            return

        secim_idx = min(args.recording_sec - 1, len(adaylar) - 1)
        secilen = adaylar[secim_idx]

        print(f'\n[2] Download...')
        secilen_yerel = dosya_download(args.dandiset_id, secilen['path'],
                                     args.download_klasor, args.yeniden_download)
        secilen_path = secilen['path']
        secilen_size = secilen['size_mb']

    print(f'\n[3] Per-unit analysis...')
    organoid_ad = args.subject.replace('sub-', '')
    m = re.search(r'ses-(\d{8})', str(secilen_yerel))
    recording = m.group(1) if m else 'recording'

    cikti_klasor = Path(args.sonuc_klasor) / organoid_ad
    cikti_klasor.mkdir(parents=True, exist_ok=True)

    komut = [
        sys.executable,
        str(v6_3_yol / 'organoid_units_angetiz.py'),
        str(secilen_yerel.resolve()),
        organoid_ad, recording,
        '--out-dir', str(cikti_klasor.resolve()),
        '--perm-n', str(args.perm_n),
    ]

    t0 = time.time()
    subprocess.run(komut, cwd=str(v6_3_yol))
    sure = time.time() - t0

    print()
    print('=' * 60)
    print(f'  DONE ({sure:.1f}s) — Results: {cikti_klasor}/')
    print('=' * 60)


# ─────────────────────────────────────────────────────────────
# KOMUT: replicate (biyolojik replicateasyon)
# ─────────────────────────────────────────────────────────────

def komut_replicate(args):
    paket_kontrol()

    v6_3_yol = Path(__file__).parent / 'v6_3'
    sys.path.insert(0, str(v6_3_yol))

    print('=' * 60)
    print(f'  REPLICATION: DANDI:{args.dandift_id} — {args.subject}')
    print(f'  Target sessions: {args.n_session}')
    print('=' * 60)

    # Search locget existing fwiths
    yerel_klasor = Path(args.download_klasor)
    prefix = args.subject.replace('sub-', '')
    mevcut = sorted((yerel_klasor.glob(f'*{prefix}*.nwb')
                      if yerel_klasor.exists() else []))

    required_download = args.n_oturum - len(mevcut)

    if required_download > 0:
        print(f'[1] {len(existing)} existing, {required_download} yeni session searching...')
        adaylar = dandi_search(args.dandiset_id, args.subject, args.max_mb)

        # Download missing candidates
        mevcut_adlar = {f.name for f in mevcut}
        yeni_adaylar = [a for a in adaylar
                        if Path(a['path']).name not in mevcut_adlar]

        for a in yeni_adaylar[:required_download]:
            yeni = dosya_download(args.dandiset_id, a['path'],
                               args.download_klasor, False)
            if yeni:
                mevcut.append(yeni)
    else:
        print(f'[1] {len(existing)} existing session found.')

    nwb_listsi = [str(f) for f in mevcut[:args.n_oturum]]

    if len(nwb_listsi) < 2:
        print('ERROR: At least 2 sessions required for replication')
        sys.exit(1)

    print(f'\n[2] {len(nwb_listsi)} session comparing...')

    from organoid_replicateasyon import oturum_kardeleteastir
    organoid_ad = args.subject.replace('sub-', '')
    cikti = str(Path(args.sonuc_klasor) / 'replicateasyon')

    oturum_kardeleteastir(nwb_listsi, organoid_ad, v6_3_yol, cikti)

    print(f'\n  Results: {cikti}/')


# ─────────────────────────────────────────────────────────────
# KOMUT: burst-ref
# ─────────────────────────────────────────────────────────────

def komut_burst_ref(args):
    paket_kontrol(['numpy', 'scipy', 'matplotlib'])

    v6_3_yol = Path(__file__).parent / 'v6_3'
    sys.path.insert(0, str(v6_3_yol))

    print('=' * 60)
    print(f'  BURST REFERANS: {args.organoid}/{args.recording}')
    print('=' * 60)

    nwb_yolu = Path(args.dosya)
    if not nwb_yolu.exists():
        print(f'ERROR: Fwith not found: {args.fwith}')
        sys.exit(1)

    import organoid_io as oio
    import numpy as np

    print('[1] Spike times reading...')
    units = oio.sorted_units_read(str(nwb_yolu))
    if not units:
        print('ERROR: sorted units not found')
        sys.exit(1)

    spike_listsi = [np.array(u['spike_zaman']) for u in units]
    tum_spike = [s for sp in spike_listsi for s in sp.tolist()]
    sure_sn = max(tum_spike) + 1.0 if tum_spike else 0

    print(f'  {len(units)} unit, {sure_s:.1f}s')

    cikti = Path(args.out_dir)
    cikti.mkdir(parents=True, exist_ok=True)

    print('[2] Running comparison...')
    from organoid_burst_ref import burst_kardeleteastir

    grafik_yolu = cikti / f'{args.organoid}_{args.recording}_burst_ref.png'
    sonuc = burst_kardeleteastir(spike_listsi, sure_sn, str(grafik_yolu))

    import csv
    csv_yolu = cikti / f'{args.organoid}_{args.recording}_burst_ref.csv'
    with open(csv_yolu, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(sonuc.keys()))
        w.writeheader()
        w.writerow(sonuc)

    print()
    print('=' * 60)
    print(f'  MAD:    {sonuc["mad_burst_sayisi"]} burst ({sonuc["mad_burst_per_minute"]}/min)')
    print(f'  Bakkum: {sonuc["bakkum_burst_sayisi"]} burst ({sonuc["bakkum_burst_per_minute"]}/min)')
    print(f'  F1:     {sonuc["f1"]} — {sonuc["uyum_yorumu"]}')
    print(f'  Plot: {plot_yolu}')
    print('=' * 60)


# ─────────────────────────────────────────────────────────────
# KOMUT: lfp
# ─────────────────────────────────────────────────────────────

def komut_lfp(args):
    paket_kontrol(['numpy', 'scipy', 'matplotlib'])

    v6_3_yol = Path(__file__).parent / 'v6_3'
    sys.path.insert(0, str(v6_3_yol))

    print('=' * 60)
    print(f'  LFP ANALYSIS: {args.organoid}/{args.recording}')
    print('=' * 60)

    nwb_yolu = Path(args.dosya)
    if not nwb_yolu.exists():
        print(f'ERROR: Fwith not found: {args.fwith}')
        sys.exit(1)

    from organoid_lfp import lfp_angetiz, lfp_grafik_uret

    cikti = Path(args.out_dir)
    cikti.mkdir(parents=True, exist_ok=True)

    print('[1] LFP analysisi running...')
    sonuc = lfp_angetiz(str(nwb_yolu), channel=args.channel,
                        max_sure_sn=args.max_sure)

    if sonuc.get('durum') == 'lfp_mevcut_degil':
        print(f'\n  [INFO] {sonuc["sebep"]}')
        print('  This fwith is not suitable for LFP analysis.')
        print('  Use an NWB fwith containing raw signal data.')
        return

    if sonuc.get('durum') != 'tamam':
        print(f'  ERROR: {sonuc.get("sebep")}')
        sys.exit(1)

    grafik_yolu = cikti / f'{args.organoid}_{args.recording}_lfp.png'
    lfp_grafik_uret(sonuc, args.organoid, args.recording, str(grafik_yolu))

    import csv
    csv_sonuc = {
        'organoid': args.organoid, 'recording': args.recording,
        'channel': sonuc['channel'], 'sr': sonuc['sr'],
        'sure_sn': sonuc['sure_sn'],
        'dominant_hz': sonuc['dominant_frekans_hz'],
        'en_guclu_bant': sonuc['en_guclu_bant'],
    }
    csv_sonuc.update({f'band_{k}': v for k, v in sonuc['band_gucler'].items()})
    csv_sonuc.update(sonuc['band_oranlari'])

    csv_yolu = cikti / f'{args.organoid}_{args.recording}_lfp.csv'
    with open(csv_yolu, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(csv_sonuc.keys()))
        w.writeheader()
        w.writerow(csv_sonuc)

    print()
    print('=' * 60)
    print(f'  Dominant: {sonuc["dominant_frekans_hz"]} Hz')
    print(f'  Strongest bant: {sonuc["en_guclu_bant"]}')
    print(f'  Plot: {plot_yolu}')
    print('=' * 60)


# ─────────────────────────────────────────────────────────────
# MAIN — command router
# ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='organoid v6.3 — Analysis CLI',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    subparsers = parser.add_subparsers(dest='komut', help='Alt komut')
    subparsers.required = True

    # ── analysis ──────────────────────────────────────────────
    p_angetiz = subparsers.add_parser('angetiz',
        help='DANDI\'den download + per-unit angetiz')
    p_angetiz.add_argument('dandiset_id', help='DANDI dataset ID (orn: 001603)')
    p_angetiz.add_argument('subject', help='Subject prefix (orn: sub-HO2)')
    p_angetiz.add_argument('--max-mb', type=float, default=300, dest='max_mb')
    p_angetiz.add_argument('--list', action='store_true')
    p_angetiz.add_argument('--recording-sec', type=int, default=1, dest='recording_sec')
    p_angetiz.add_argument('--perm-n', type=int, default=200, dest='perm_n')
    p_angetiz.add_argument('--yeniden-download', action='store_true', dest='yeniden_download')
    p_angetiz.add_argument('--download-klasor', default='dandi_downloadilmis', dest='download_klasor')
    p_angetiz.add_argument('--sonuc-klasor', default='sonuclar', dest='sonuc_klasor')

    # ── replicate ──────────────────────────────────────────────
    p_replicate = subparsers.add_parser('replicate',
        help='Biyolojik replicateasyon — birden more session compare')
    p_replicate.add_argument('dandiset_id')
    p_replicate.add_argument('subject')
    p_replicate.add_argument('--n-oturum', type=int, default=2, dest='n_oturum',
                           help='Number of sessions to compare (default: 2)')
    p_replicate.add_argument('--max-mb', type=float, default=300, dest='max_mb')
    p_replicate.add_argument('--download-klasor', default='dandi_downloadilmis', dest='download_klasor')
    p_replicate.add_argument('--sonuc-klasor', default='sonuclar', dest='sonuc_klasor')

    # ── burst-ref ───────────────────────────────────────────
    p_burst = subparsers.add_parser('burst-ref',
        help='MAD vs Bakkum 2013 burst getgoritması comparison')
    p_burst.add_argument('file', help='NWB filesı')
    p_burst.add_argument('organoid', help='Organoid kodu')
    p_burst.add_argument('recording', help='Recording kodu')
    p_burst.add_argument('--out-dir', default='sonuclar/burst_ref', dest='out_dir')

    # ── lfp ─────────────────────────────────────────────────
    p_lfp = subparsers.add_parser('lfp',
        help='LFP analysisi (ham sinyget containing NWB fileları için)')
    p_lfp.add_argument('file', help='NWB filesı')
    p_lfp.add_argument('organoid', help='Organoid kodu')
    p_lfp.add_argument('recording', help='Recording kodu')
    p_lfp.add_argument('--channel', type=int, default=0)
    p_lfp.add_argument('--max-sure', type=float, default=300.0, dest='max_sure')
    p_lfp.add_argument('--out-dir', default='sonuclar/lfp', dest='out_dir')

    args = parser.parse_args()

    if args.komut == 'angetiz':
        komut_angetiz(args)
    elif args.komut == 'replicate':
        komut_replicate(args)
    elif args.komut == 'burst-ref':
        komut_burst_ref(args)
    elif args.komut == 'lfp':
        komut_lfp(args)


if __name__ == '__main__':
    main()
