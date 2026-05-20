# organoid v6.3

**An open-source electrophysiology analysis pipeline for NWB-formatted brain organoid recordings.**

Developed independently. Validated on [DANDI:001603](https://dandiarchive.org/dandiset/001603) (van der Molen et al., *Nature Neuroscience* 2026) — 18 subjects, ~1,610 individual neurons across human organoids, mouse organoids, and ex vivo mouse cortical slices.

---

## What it does

- **Per-unit analysis:** ISI statistics, CV, firing rate, refractory violation, Wilson score CI, ISI shuffle permutation (rhythmicity test)
- **Network synchrony:** STTC cross-unit synchrony matrices (Cutts & Eglen 2014)
- **Burst detection:** Adaptive MAD-threshold network burst detection
- **Graph theory:** Clustering coefficient, mean path length, small-world index (NetworkX)
- **Multi-session stability:** CV_session across recording sessions
- **DANDI integration:** Search, download, and analyze any NWB subject in one command
- **Raw signal fallback:** Interactive option for subjects without sorted units

---

## Installation

```bash
pip install pynwb dandi numpy scipy matplotlib pandas networkx
```

Python 3.10+ required.

---

## Quick Start

```bash
# Search DANDI for organoid datasets
python dandi_ara.py organoid

# List files in a specific dataset
python dandi_ara.py --dataset 001603

# Download and analyze a subject (sorted units mode)
python organoid_cli.py 001603 sub-HO2

# List available recordings for a subject
python organoid_cli.py 001603 sub-HO2 --listele

# Select a specific recording
python organoid_cli.py 001603 sub-HO2 --kayit-sec 3

# Run graph theory analysis on downloaded subjects
python graph_analiz.py
```

Output files are saved to `sonuclar/<subject>/`:
- `*_per_unit.csv` — per-neuron metrics
- `*_ozet.csv` — organoid-level summary
- `*_quicklook.png` — 4-panel visualization (raster, population rate, Hz histogram, STTC matrix)

---

## Pipeline Architecture

```
organoid_cli.py          ← CLI entry point: DANDI search + download + analysis
dandi_ara.py             ← DANDI-wide dataset and file search
graph_analiz.py          ← Graph theory from STTC matrices
gruplu_analiz.py         ← Multi-organoid group comparison
yas_metadata_cek.py      ← Extract age metadata from NWB files

v6_3/
├── organoid_io.py       ← NWB reading, chunk management, RAM safety
├── organoid_signal.py   ← Bandpass filter, MAD-threshold spike detection
├── organoid_qc.py       ← SNR, refractory violation, Wilson CI
├── organoid_metrics.py  ← ISI, CV, STTC, permutation test
├── organoid_units_analiz.py  ← Per-unit analysis, STTC matrix, burst detection
├── organoid_burst_ref.py     ← Bakkum 2013 reference burst comparison
├── organoid_lfp.py      ← LFP band power analysis
├── organoid_replikasyon.py   ← Multi-session stability
├── organoid_compare.py  ← Two-recording statistical comparison
├── organoid_db.py       ← SQLite + BLOB storage
├── organoid_output.py   ← Console output, CSV writing
├── organoid_plot.py     ← Quicklook PNG generation
├── organoid_sorting.py  ← SpikeInterface integration (optional)
└── organoid_types.py    ← Shared type definitions
```

---

## Validation

15 synthetic tests — 100% pass rate. Key tests:

| Test | Expected | Result |
|------|----------|--------|
| Spike detection | 300 spikes ±30% | 250 (σ=5 threshold) |
| ISI CV (Poisson) | 0.85–1.15 | 0.944 |
| STTC (identical trains) | ≈1.0 | 1.000 |
| STTC (independent) | ≈0.0 | −0.015 |
| MAD robustness | std +50%, MAD <10% | std +403%, MAD +0.8% |
| RAM safety | ValueError on overflow | 18 GB > 500 MB → caught |

---

## Applied to DANDI:001603

All 18 subjects with sorted units analyzed:

| Category | n | Median Hz | Median STTC | Burst/min |
|----------|---|-----------|-------------|-----------|
| Human organoid (HO1–HO8) | 8 | 0.42 | 0.159 | 9.84 |
| Mouse organoid (MO10–MO14) | 5 | 4.42 | 0.014 | 5.40 |
| Ex vivo mouse cortex (M1S1–M3S2) | 5 | 0.39 | 0.144 | 10.10 |

Age group comparison (human organoids):
- G1 (6–7 months): STTC = 0.295 ± 0.096
- G2 (~3.3 months): STTC = 0.023 ± 0.008
- Mann-Whitney p = 0.029

---

## Known Limitations

- Sorted units mode only (NWB `units` table required). Raw signal mode is experimental.
- Recording duration and organoid age are confounded in DANDI:001603 (G1: 3 min, G2: 10 min).
- Graph theory sigma values unreliable for small networks (n=20); σ > 2 results more trustworthy.
- n = 4–8 per group — indicative, not conclusive.

---

## References

- van der Molen et al. (2026). Preconfigured neuronal firing sequences in human brain organoids. *Nature Neuroscience*, 29(1), 123–135. DANDI:001603.
- Cutts & Eglen (2014). Detecting pairwise correlations in spike trains. *Journal of Neuroscience*, 34(43), 14288–14303.
- Buzsáki & Mizuseki (2014). The log-dynamic brain. *Nature Reviews Neuroscience*, 15(4), 264–278.
- Rübel et al. (2022). The Neurodata Without Borders ecosystem. *eLife*, 11, e78362.

---

*Data: DANDI:001603, CC-BY-4.0. Pipeline: independent development. Not peer reviewed.*
