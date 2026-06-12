---
title: 'Axon: A Standardized Open-Source Pipeline for Brain Organoid MEA Electrophysiology'
tags:
  - Python
  - neuroscience
  - brain organoid
  - MEA
  - electrophysiology
  - NWB
  - DANDI
  - spike sorting
  - STTC
  - organoid intelligence
authors:
  - name: Metin
    affiliation: 1
affiliations:
  - name: Independent Researcher
    index: 1
date: 12 June 2026
bibliography: paper.bib
---

# Summary

Brain organoid microelectrode array (MEA) electrophysiology is a rapidly growing field with no established analysis standards. Different groups use different spike detection thresholds, burst detection algorithms, and synchrony metrics, making cross-study comparison difficult. `Axon` (v6.3) is an open-source Python pipeline for end-to-end analysis of Neurodata Without Borders (NWB)-formatted brain organoid MEA recordings with direct DANDI Archive integration. The pipeline computes per-unit firing rate, inter-spike interval (ISI) statistics, coefficient of variation (CV), Spike Time Tiling Coefficient (STTC) cross-unit synchrony matrices, adaptive network burst detection, and graph-theoretic network topology via a single command-line interface. Version 6.3 adds full SpikeInterface integration for spike sorting of raw high-density MEA data, supporting four algorithms with DC-corrected RMS-based channel selection and per-unit quality metrics.

# Statement of Need

The Neurodata Without Borders (NWB) format [@ruebel2022nwb] and the DANDI Archive [@dandi] provide a common data standard for neurophysiology, but the analytical layer for brain organoid MEA data remains fragmented. Existing general-purpose spike analysis tools (e.g., SpikeInterface [@buccino2020spikeinterface], MountainSort [@chung2017mountainsort]) do not provide organoid-specific metrics such as STTC-based synchrony matrices, adaptive burst detection calibrated for organoid firing statistics, or graph-theoretic topology analysis. `Axon` fills this gap by providing a reproducible, dataset-agnostic pipeline specifically designed for NWB-formatted organoid MEA data with DANDI integration.

A critical technical contribution is the DC-offset correction step for HD-MEA raw signal data. Without mean subtraction before RMS computation, uint16-to-int16 converted recordings exhibit near-identical RMS values (~32,000 ADU) across all channels, rendering channel selection meaningless. After DC correction, RMS reflects genuine signal amplitude variation (13–68 μV across channels in DANDI:001132), enabling meaningful selection of the most active electrode subset.

# Functionality

`Axon` automatically detects NWB file type: if a units table is present, per-unit analysis proceeds directly; if only raw ElectricalSeries data is present, the user is offered threshold detection or SpikeInterface spike sorting. All outputs are written to a subject-specific results directory as CSV files and a 4-panel quicklook figure.

**Key algorithms:**

- **Spike Time Tiling Coefficient (STTC)** [@cutts2014sttc]: Cross-unit pairwise synchrony, independent of firing rate, with coincidence window dt = 10 ms.
- **Adaptive burst detection**: Population firing rate binned at 50 ms; burst threshold set using Median Absolute Deviation (MAD), requiring ≥3 co-firing units.
- **Graph-theoretic topology**: Clustering coefficient, mean path length, and small-world index σ computed via NetworkX [@hagberg2008networkx].
- **SpikeInterface integration**: Four sorting algorithms (mountainsort5, spykingcircus2, tridesclous2, simple) with DC-corrected RMS channel selection (default: top 64 of up to 942 channels).

# Validation

`Axon` was validated on 29 subjects from four public DANDI datasets spanning human and mouse brain organoids, ex vivo cortical tissue, and chronic dissociated cultures (Table 1). Synthetic unit tests (15/15 pass) cover spike detection accuracy, STTC correctness on known inputs, MAD robustness, RAM overflow protection, and Wilson CI correctness.

| Dataset | Description | Subjects | Data type |
|---------|-------------|----------|-----------|
| DANDI:001603 | Human and mouse brain organoids, ex vivo cortex | 18 | Sorted units |
| DANDI:001132 | Human hippocampal organoids, HD-MEA | 5 | Raw signal |
| DANDI:000774 | In vitro MEA | 5 | Sorted units |
| DANDI:001611 | Chronic dissociated cortical culture | 1 | Sorted units |

: Datasets used for validation. \label{tab:datasets}

A key finding from DANDI:001603 was a statistically significant age-dependent difference in synchrony between early- and late-stage human organoids (Mann-Whitney U, p = 0.029, n = 4 per group), demonstrating the pipeline's sensitivity to biologically relevant differences.

# Acknowledgements

The author thanks the DANDI Archive team and the NWB community for open data infrastructure, and the SpikeInterface developers for the spike sorting framework.

# References
