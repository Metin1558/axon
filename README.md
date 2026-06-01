# Axon

**Open-source pipeline for brain organoid MEA electrophysiology analysis.**

Axon processes NWB-formatted brain organoid MEA recordings from the DANDI Archive end-to-end. Per-unit firing rate, ISI statistics, coefficient of variation, STTC cross-unit synchrony matrices, adaptive network burst detection, graph-theoretic network topology, and SpikeInterface spike sorting — via a single CLI command.

Validated on 29 subjects from four public DANDI datasets: DANDI:001603 (van der Molen et al., *Nature Neuroscience* 2026), DANDI:001132 (Andrews et al.), DANDI:000774, and DANDI:001611.

**Preprint:** [bioRxiv — link TBA]

---

## Repository Structure

```
axon/
├── KURULUM.py                   # Setup: installs all dependencies including SpikeInterface
├── organoid_cli.py              # Main CLI: search, download, and analyze
├── dandi_ara.py                 # DANDI Archive dataset search and file listing
├── graph_analiz.py              # Graph theory analysis (auto-detects subjects)
├── gruplu_analiz.py             # Group comparison with Mann-Whitney U
├── yas_metadata_cek.py          # NWB metadata extraction (age, species, acquisition)
├── si_sorting.py                # SpikeInterface spike sorting module
├── batch_analiz.py              # Batch analysis across multiple subjects
│
├── v6_3/                        # Core analysis modules
│   ├── organoid_units_analiz.py # Per-unit orchestrator: STTC, burst, plot
│   ├── organoid_io.py           # NWB reading, file type detection, RAM safety
│   ├── organoid_signal.py       # Bandpass filter, MAD-threshold spike detection
│   ├── organoid_metrics.py      # ISI, CV, STTC (Cutts & Eglen 2014)
│   ├── organoid_qc.py           # Refractory violation, Wilson CI, SNR
│   ├── organoid_burst_ref.py    # Bakkum 2013 reference burst comparison
│   ├── organoid_replikasyon.py  # Multi-session stability analysis
│   ├── organoid_compare.py      # Two-recording statistical comparison
│   ├── organoid_plot.py         # 4-panel quicklook visualization
│   └── organoid_types.py        # Data types
│
└── paper/
    ├── axon_preprint.tex        # LaTeX source
    ├── references.bib           # BibTeX references
    └── Axon_Preprint_LaTeX.pdf  # Compiled PDF
```

---

## Installation

```bash
git clone https://github.com/metin/axon
cd axon
python KURULUM.py
```

**Requirements:** Python 3.9+, internet connection

**Automatically installed:** `dandi`, `pynwb`, `h5py`, `numpy`, `scipy`, `matplotlib`, `networkx`, `remfile`, `pandas`, `spikeinterface[full]`, `mountainsort5`, `tridesclous2`, `spykingcircus2`

---

## Quick Start

```bash
# Search DANDI datasets
python dandi_ara.py organoid

# List files in a dataset
python dandi_ara.py --dataset 001603

# Download and analyze (sorted units)
python organoid_cli.py 001603 sub-HO2

# Download and analyze (raw signal + SpikeInterface sorting)
python organoid_cli.py 001132 sub-5C-1 --max-mb 750
# → Select option 2 → Select algorithm (1: mountainsort5 recommended)

# Batch analysis across multiple subjects
python batch_analiz.py --dataset 001611 --n 5 --max-mb 30

# Graph theory analysis (auto-reads results/ directory)
python graph_analiz.py

# Group comparison
python gruplu_analiz.py

# NWB metadata summary
python yas_metadata_cek.py
```

---

## Outputs

Each analysis writes to `results/<subject>/`:

| File | Content |
|------|---------|
| `*_per_unit.csv` | Per-neuron metrics: Hz, CV, STTC, ISI, refractory violation, SNR, quality flag |
| `*_summary.csv` | Population summary: n_units, burst rate, median STTC |
| `*_quicklook.png` | 4-panel figure: raster, population rate, Hz distribution, STTC matrix |

---

## SpikeInterface Integration

For raw electrophysiology data (no pre-sorted units), Axon integrates SpikeInterface with:

- **DC-corrected RMS channel selection** — removes uint16→int16 DC offset before RMS computation; selects top N active channels (default N=64)
- **Four sorting algorithms** — mountainsort5 (recommended), spykingcircus2, tridesclous2, simple
- **Per-unit quality metrics** — SNR (threshold >3.0), ISI violation ratio (<10%), firing rate (>0.05 Hz)

| Algorithm | Notes |
|-----------|-------|
| mountainsort5 | Recommended — fast, stable, scheme1 |
| spykingcircus2 | High channel-count HD-MEA |
| tridesclous2 | Alternative template matching |
| simple | Lowest RAM usage |

---

## Validation

**Synthetic tests:** 15/15 pass

**DANDI:001603** (18 subjects, ~1,610 neurons, sorted units):

| Category | n | Median Hz | Median CV | Median STTC | Burst/min |
|----------|---|-----------|-----------|-------------|-----------|
| Human organoid (HO1–HO8) | 8 | 0.42 | 1.92 | 0.159 | 9.84 |
| Mouse organoid (MO10–MO14) | 5 | 4.42 | 1.15 | 0.014 | 5.40 |
| Ex vivo mouse cortex (M1S1–M3S2) | 5 | 0.39 | 2.83 | 0.144 | 10.10 |

Age-group difference (human organoids): Mann-Whitney p=0.029 for STTC and CV (n=4 per group).

Small-world topology: HO6 and HO7 (P100D) σ ≈ 4.3; older organoids σ ≈ 1.0.

**DANDI:001132** (5 subjects, HD-MEA raw signal, SpikeInterface):
- 745–972 channels, DC-corrected RMS → 64 channels selected
- mountainsort5: 2–15 units sorted per subject, 67–87% qualified

**DANDI:001611** (1 subject, scalability test):
- 950 units, 974,562 spikes, 261,569 s (~3 days) recording
- Full pipeline completed without modification

---

## Citation

If you use Axon in your work, please cite:

```
Metin (2026). Axon: A Standardized Open-Source Pipeline...
Zenodo: https://doi.org/10.5281/zenodo.20492811
```

---

## Data Availability

All datasets are publicly available on the DANDI Archive (CC-BY-4.0):
- [DANDI:001603](https://dandiarchive.org/dandiset/001603)
- [DANDI:001132](https://dandiarchive.org/dandiset/001132)
- [DANDI:000774](https://dandiarchive.org/dandiset/000774)
- [DANDI:001611](https://dandiarchive.org/dandiset/001611)

---

*Axon v6.3 · Independent Research · Not peer reviewed · May 2026*
