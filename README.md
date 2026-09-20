# Axon v6.3 + TBL1XR1 Extension

This repository bundles two things together, unmodified where possible and
clearly separated where new:

1. **Axon** — Metin's MEA (microelectrode array) electrophysiology analysis
   toolkit, validated across 4 DANDI Archive datasets and 29 subjects
   (human and mouse brain organoids). This is the same codebase that lives
   on GitHub today, plus three bug fixes described in detail below.
2. **The TBL1XR1 extension** — a new, gene-specific hypothesis-testing layer
   built on top of Axon, designed as a gene-specific hypothesis-testing
   layer for TBL1XR1 brain organoid electrophysiology.
   It does not duplicate any of Axon's logic; it imports and orchestrates it.

Everything in this document reflects what was actually built, actually
tested, and actually found while building it — including a real bug found
in this same update cycle, and a real methodological limitation discovered
while testing the new longitudinal-comparison wrapper. Nothing here is
aspirational; every number quoted below came from an actual test run of
the code in this archive.

---

# PART 1 — Axon (standard usage, unchanged from the 29-subject analysis)

## 1.1 What Axon does

Axon is a command-line tool and Python library that extracts standard
electrophysiology metrics from MEA recordings stored in NWB (Neurodata
Without Borders) format — whether sourced from the DANDI Archive or any
other NWB-producing pipeline. It was built to be disease-agnostic: none of
its core metrics assume anything about which organism, genotype, or
condition produced the recording.

Its capabilities, grouped by what question they answer:

**"How active is this tissue?"**
- Inter-spike interval (ISI) statistics: mean, median ISI
- Coefficient of variation (CV) of the ISI distribution — a measure of
  firing regularity (CV=0 is perfectly regular, CV→1 approaches Poisson)
- Per-unit and population-level firing rate (Hz)

**"How synchronized is this tissue?"**
- STTC (Spike Time Tiling Coefficient), following Cutts & Eglen (2014) —
  a rate-independent measure of pairwise spike-train correlation that
  avoids the firing-rate confound that plagues cross-correlation-based
  synchrony measures
- Full pairwise STTC matrices across a population of units

**"Is this tissue producing coordinated network events?"**
- Adaptive, MAD-based (median + k × median absolute deviation) network
  burst detection on the population spike-rate time series — this
  threshold adapts to each recording's own baseline rather than using a
  fixed cutoff, which matters because baseline activity varies enormously
  across organoids, ages, and recording conditions

**"Is the network topologically organized?"**
- Graph-theoretic metrics computed from the STTC adjacency matrix:
  density, average clustering coefficient, average shortest path length,
  and small-world sigma (σ) — the ratio of clustering-to-path-length
  relative to a random graph of the same density, a standard measure of
  whether a network has "small-world" organization

**"Can I trust this data before I compare groups?"**
- Quality-control metrics that report numbers, not verdicts: refractory-
  period violation rate (with a Wilson score 95% confidence interval),
  negative-ISI count, firing-rate outlier density within sliding windows,
  and channel silence detection. The module's own internal philosophy,
  stated in its docstring, is explicit: it does not label anything
  "good" or "bad" — it reports "X.Y% refractory violation, 95% CI [a, b]"
  and leaves interpretation to the analyst.

**"Are these two groups/recordings different?"**
- Mann-Whitney U, unpaired t-test, one-way ANOVA for n-group comparisons
- Circular-shift permutation testing for autocorrelated, windowed data
  (i.e., comparing two recordings by their within-recording windowed
  firing-rate/CV distributions rather than treating each window as an
  independent sample)
- Lag-1 autocorrelation reporting alongside every windowed comparison, so
  the analyst can see whether the permutation-based correction was
  actually necessary for that particular dataset

**"Is my measurement itself reliable, session to session?"**
- Intraclass correlation coefficient (ICC) computation across repeated
  sessions of the same preparation

**"How do I go from raw voltage to spike times?"**
- Bandpass filtering (default 300–3000 Hz, order-4 Butterworth),
  MAD-based noise-floor estimation, chunked spike detection with
  signal-to-noise ratio computation per detected event, and saturation
  (clipping) detection on the raw trace
- An entry point into SpikeInterface for full spike sorting when the
  incoming data is raw waveforms rather than pre-sorted units

**"What about local field potentials?"**
- Band-power decomposition into delta/theta/alpha/beta/gamma bands,
  band-power ratios, dominant-frequency extraction via Welch PSD, and
  windowed band-power time courses

**"How do I compare a reference burst-detection method against the one Axon uses?"**
- An implementation of the Bakkum et al. (2013, *Nature Communications*)
  burst-detection algorithm, run side-by-side against Axon's own MAD
  algorithm on the same data, with IoU (Intersection-over-Union),
  precision, recall, and F1 agreement scoring, plus a diagnostic plot

**"How do I run this at scale, and keep results somewhere?"**
- SQLite-backed result storage (`organoid_db.py`), CSV/terminal report
  generation (`organoid_output.py`), batch analysis across many subjects
  (`batch_analiz.py`), a DANDI Archive search-and-download front end
  (`dandi_ara.py`), and a unified CLI (`organoid_cli.py`)

This toolkit was run and validated across DANDI:001603 and three other
DANDI datasets, covering 29 subjects total, spanning both human and mouse
organoid preparations at ages ranging from P100D (~14 weeks) to 6–7
months. The full validation write-up, including every statistical test
and its result, lives in `paper/axon_preprint.tex` (with the compiled PDF
alongside it) — that document is unchanged by this update.

## 1.2 Directory structure, file by file

```
v6_3/                        Core analysis library (17 files)
├── organoid_io.py             Reads NWB files; extracts sorted-unit spike times
├── organoid_signal.py         Raw voltage → bandpass filter → spike detection
│                               (bandpass_filter, mad_std, spike_tespit_chunk,
│                               snr_compute, doygunluk_tespit/saturation check)
├── organoid_sorting.py        Orchestrates the spike-sorting workflow
├── si_sorting.py               SpikeInterface integration for full sorting
├── organoid_qc.py             Quality metrics: refractory_ihlget_orani (refractory
│                               violation rate + Wilson CI), negatif_isi_sayisi,
│                               yogunluk_outlier_orani, channel_sessizligi
├── organoid_metrics.py        Core metrics: isi_metrikleri, sttc_iki_channel,
│                               coklu_channel_senkron, mannwhitney_iki_grup,
│                               ttest_iki_grup, anova_n_grup, aktiflik_window
├── organoid_units_analiz.py   network_burst_tespit — the adaptive MAD burst detector
├── organoid_burst_ref.py      Bakkum 2013 reference comparison  [FIXED — see 1.3]
├── organoid_lfp.py            LFP band-power analysis (band_gucu_hesapla,
│                               dominant_frekans, pencereli_band_gucu)
├── organoid_compare.py        Pairwise recording comparison (ikili_karsilastir,
│                               circular_shift_permutation, autocorrelation_lag1)
├── organoid_replikasyon.py    Session-to-session reliability (icc_hesapla)
├── organoid_db.py             SQLite storage for spike times, metrics, QC results
├── organoid_output.py         Terminal and CSV report writers
├── organoid_plot.py           Visualization helpers
├── organoid_types.py          KayitBilgisi dataclass (recording metadata)
├── organoid.py                 Single-NWB-file analysis entry point
└── test_organoid_metrics.py    15-assertion synthetic-data validation suite  [FIXED — see 1.3]

axon/                         CLI and batch tooling (7 files)
├── organoid_cli.py             Main CLI: search, download, and dispatch analysis
├── dandi_ara.py                 DANDI Archive search
├── batch_analiz.py              Multi-subject batch analysis
├── gruplu_analiz.py             Automatic grouping + group comparison
│                                (gruplara_ayir — see the documented limitation
│                                in 1.3)
├── graph_analiz.py              Network topology: sttc_matrisi_hesapla,
│                                graph_metrikleri_hesapla (density, clustering,
│                                small-world sigma)
├── yas_metadata_cek.py          Extracts age metadata from NWB files
└── KURULUM.py                   Installation helper script

paper/                        Axon's own preprint (PDF, LaTeX source, bibliography)
├── Axon_Preprint_LaTeX.pdf
├── axon_preprint.tex
└── references.bib

LICENSE, requirements.txt     Unchanged from the original repository
```

## 1.3 Bug fixes applied in this update

Three concrete, verifiable defects were found and fixed. After every fix,
the original 15-assertion test suite was re-run and confirmed to still
pass in full — **no regression was introduced**, and the existing
29-subject analysis workflow is unaffected by any of these changes.

### Fix 1 — `organoid_burst_ref.py`: threshold contradicted its own citation

**Before:** the function `bakkum_burst_tespit()` carried a docstring
explicitly stating *"Bakkum 2013'te 5 kullanılmış"* (Bakkum 2013 used a
value of 5 for the threshold multiplier), yet the function's own default
parameter was `esik_carpan=2.5` — exactly half the cited reference value,
with no explanation anywhere in the code for the deviation.

**Why this mattered:** the entire purpose of this module is to compare
Axon's own MAD-based burst detector against Bakkum et al.'s published
method as an external reference check. If the "Bakkum" implementation
silently uses a threshold that is not Bakkum's own threshold, the
resulting IoU/F1 agreement scores are not actually measuring agreement
with the cited method — they are measuring agreement with an
undocumented, unrelated variant of it. That undermines the validity of
using this module as a reference-method sanity check at all.

**Fix applied:** the default was changed to `esik_carpan=5`, matching the
cited value exactly. The docstring was updated to note the correction and
explain the reasoning, so a future reader understands why the value is
what it is rather than re-discovering the discrepancy independently.

### Fix 2 — `organoid_burst_ref.py`: diagnostic plot's time axis was wrong by a factor of 5

**Before:** the internal plotting function `_grafik_uret()` computed its
time axis as:
```python
t = np.arange(len(pop_rate)) * 10 / 1000.0  # hardcoded: assumes 10 ms bins
```
regardless of what bin width was actually used to compute `pop_rate`. In
practice, `pop_rate` is produced by `bakkum_burst_tespit()`, which is
called elsewhere in the same file with `bin_ms=50` (the function's own
default). So the plotted population-rate curve was being drawn on a time
axis compressed by a factor of 5 relative to the recording's real
duration — while the burst-interval overlay rectangles were drawn using
the correct, real-second timestamps stored in each burst dictionary.

**Effect:** the diagnostic PNG this function produces would show the
population-rate trace squeezed into roughly one-fifth of the plot's true
time span, with the burst-highlight rectangles (correctly scaled)
appearing badly misaligned against it — a plot that visually looks broken
even though the underlying detection was fine.

**What was *not* affected:** every numeric output of this module — burst
counts, IoU, precision, recall, F1, median burst duration — is computed
directly from the burst dictionaries and the raw `pop_rate` array, never
from the plot's x-axis. So this was a purely visual bug; no downstream
statistic was ever wrong because of it. It's still worth fixing, because
a visual sanity-check tool that itself needs a sanity check defeats its
own purpose.

**Fix applied:** `bin_ms` is now threaded explicitly from `burst_karsilastir()`
into `_grafik_uret()` as a parameter and used to build the time axis
(`t = np.arange(len(pop_rate)) * bin_ms / 1000.0`), so the axis always
reflects whatever bin width was actually used, whatever that value is
changed to in the future.

### Fix 3 — `test_organoid_metrics.py`: a hardcoded path to the author's own machine

**Before:**
```python
sys.path.insert(0, '/tmp/axon_final/v6_3')
```
This is an absolute path specific to whichever machine originally
developed this test file. On any other machine — including, most likely,
the machine of anyone cloning this repository from GitHub — this line
either silently inserts a non-existent directory into `sys.path` (harmless
but pointless) or, worse, could silently pick up a stale, unrelated
`organoid_metrics.py` if such a path happened to exist and contain one.
The tests happened to still pass on most machines only because Python
already searches a script's own directory by default — this line was
never actually necessary for the tests to work, it was simply a leftover
that could not do anything correct on anyone else's machine.

**Fix applied:**
```python
sys.path.insert(0, str(Path(__file__).resolve().parent))
```
The path is now derived from the test file's own location, so it is
correct on every machine, including whichever one this is eventually run on.

### A documented (not fixed) limitation: `gruplu_analiz.py`'s automatic grouping

This is **not** counted among the three fixes above, because it is not a
bug in the sense of "produces a wrong answer where a right answer was
expected" — it is a design choice with a real, demonstrated failure mode,
and fixing it properly would require adding genuine genotype/age metadata
fields to the pipeline (a larger design change than a bug fix, and one
that risks changing behavior for the existing 29-subject workflow).

`gruplara_ayir()`, the automatic grouping function in `gruplu_analiz.py`,
partitions recordings into groups using only two signals: the subject's
name prefix (e.g., "HO" for human organoid) and a recording-duration
threshold (default 300 seconds — recordings shorter than this go into one
bucket, longer ones into another). This was empirically tested during
this project: synthetic recordings were constructed with **no age
information whatsoever**, only a name prefix and a duration matching the
real durations used in Axon's own published age-group comparison (3-minute
vs. 10-minute recordings, corresponding to the paper's "G1" vs. "G2" age
groups). Running `gruplara_ayir()` on this age-blind synthetic input
reproduced the exact same subject partition — {HO1, HO2, HO3, HO4} versus
{HO5, HO6, HO7, HO8} — that the published age-group comparison reports.

This does not necessarily mean the *published* age comparison in the
paper was computed via this specific heuristic rather than genuine age
metadata (that could not be fully confirmed from the code alone) — but it
does mean that **recording duration and subject age are perfectly
collinear in at least this one real dataset**, to the point where an
age-blind heuristic based purely on duration can fully recover the "age"
grouping. Anyone using `gruplara_ayir()` going forward, on new data where
duration and true group identity are *not* collinear by coincidence, risks
silently grouping by an accidental proxy variable instead of the variable
they actually intend to compare.

**No code was changed for this.** The function still behaves exactly as
it did in the current GitHub version, so the 29-subject workflow is
unaffected. The recommendation, stated here rather than enforced in code,
is: for any new grouped comparison, prefer explicit age/genotype metadata
(as extracted by `yas_metadata_cek.py`, or supplied directly) over this
automatic heuristic, and if the heuristic must be used, verify that
recording duration is not itself the thing driving the apparent group
difference — precisely the kind of duration/genotype confound this whole
extension was built to guard against (see Part 2, §2.4).

## 1.4 Installation and usage (unchanged)

```bash
pip install -r requirements.txt
python3 axon/organoid_cli.py --help
python3 v6_3/test_organoid_metrics.py   # should print 15/15 passed
```

Nothing about how you invoke Axon has changed. Every command that worked
before this update works identically now, with the three corrected
behaviors described above.

---

# PART 2 — The TBL1XR1 Extension

## 2.1 Why this exists

TBL1XR1 is a monogenic neurodevelopmental disease gene for which brain
organoid modeling is an active area of research interest — organoid
electrophysiology is a natural readout for a gene whose loss-of-function
and missense variants produce autism, intellectual disability, seizures,
and other neurodevelopmental phenotypes in patients (see the literature
review in §2.2 for the full clinical picture).

Axon's measurement core — firing rate, ISI/CV, STTC, network bursts,
network topology — is already disease-agnostic; none of it assumes
anything about TBL1XR1 specifically. What was genuinely missing was the
gene-specific knowledge layer: for a given TBL1XR1 mutation class, which
metric should move, in which direction, and with how much confidence —
and, just as importantly, for which mutation classes does no such
prediction currently exist in the literature at all.

That knowledge layer could not be guessed or assumed. It required reading
the primary literature. Five papers were read in full (not abstract-only)
for this extension, and their findings are what the code below encodes.

## 2.2 The five papers, and exactly what each one contributed

**Heinen et al. 2016**, *Journal of Medical Genetics* — *"A specific
mutation in TBL1XR1 causes Pierpont syndrome."* Whole-exome sequencing of
six Pierpont syndrome patients identified a single recurrent de novo
missense variant, c.1337A>C (p.Tyr446Cys), located in the WD40 domain.
Critically, the mutant protein was purified and shown to assemble
*correctly* into the NCoR/SMRT/HDAC3 corepressor complex — the pathology
is not a folding or complex-assembly defect, but an as-yet-unidentified
disrupted protein–protein interaction, consistent with a dominant-negative
mechanism distinct from simple loss-of-function. The paper is explicit
that Pierpont patients with this specific mutation do **not** present with
autism, in contrast to patients carrying whole-gene deletions or other
loss-of-function variants, who frequently do.

**Mastrototaro et al. 2021**, *Frontiers in Cell and Developmental
Biology* — *"TBL1XR1 Ensures Balanced Neural Development Through NCOR
Complex-Mediated Regulation of the MAPK Pathway."* This is the single most
important paper for this extension, for two reasons.

First, it is the **only source of direct electrophysiological data**
identified anywhere in this literature review. In whole-cell voltage-clamp
recordings from P30 mouse cortical slices, *Tbl1xr1* knockout neurons
showed a significantly elevated spontaneous excitatory postsynaptic
current (sEPSC) frequency compared to wild-type littermates: 6.3 ± 0.4 Hz
versus 4.0 ± 0.3 Hz (n=19 cells per group, p<0.001), with amplitude,
10–90% rise time, and decay time constant all unchanged, and no difference
in inhibitory (sIPSC) parameters. This is a hyperexcitability phenotype at
the single-cell, synaptic-current level — not a network-level burst
phenotype, and not measured in organoid tissue at all, but it is the most
direct, quantitative electrophysiological number this literature offers,
and it is what calibrates this extension's synthetic data (§2.5) and its
choice of primary metric (§2.4).

Second, the paper performed a **complementation experiment** that turned
out to be decisive for this extension's design: *Tbl1xr1* knockout neural
stem cells were transduced with lentiviral constructs re-expressing either
wild-type TBL1XR1 or one of four disease-associated point mutations —
F10L (schizophrenia), G70N (West syndrome), L282P (autism spectrum
disorder), and Y446C (Pierpont syndrome) — and assayed for rescue of the
knockout's proliferation/differentiation defect. Wild-type rescued fully.
F10L was the **only** mutation that failed to rescue — it behaved
identically to the null allele. G70N, L282P, and Y446C all rescued the
defect essentially normally. In other words: whatever mechanism produces
the West syndrome, ASD, and Pierpont phenotypes in humans, it is *not* the
same NCoR-complex/NPC-proliferation mechanism disrupted by true
loss-of-function or by F10L. This single result is why this extension
treats those three mutation classes as having **no predicted direction**
rather than inheriting the loss-of-function prediction, and it is
discussed further in §2.4.

The paper additionally found that TBL1XR1 loss disrupts MAPK/ERK signaling
(via loss of NCOR1/HDAC3 protein levels, not just mislocalization),
produces a reduced-proliferation/premature-differentiation phenotype in
neural progenitors (matching what earlier conversation summaries in this
project already knew about), and produces behavioral deficits in adult
knockout mice (motor coordination, spatial working memory, sociability)
that recapitulate aspects of the human phenotype.

**Nishi et al. 2017**, *Scientific Reports* — *"De novo non-synonymous
TBL1XR1 mutation alters Wnt signaling activity."* Whole-exome sequencing
of 18 Japanese schizophrenia trios identified a de novo TBL1XR1 mutation,
c.30C>G (p.Phe10Leu), absent from 1,191 independent schizophrenia patients
and 1,986 controls (i.e., a private, ultra-rare variant in this one
patient). Structural modeling predicted decreased stability of the
N-terminal domain. Co-immunoprecipitation in 293FT and HT22 (mouse
hippocampal) cells confirmed the predicted mechanism directly: F10L
**decreases** TBL1XR1's interaction with N-CoR and **increases** its
interaction with β-catenin, and a TOPFlash Wnt-reporter assay confirmed
that F10L correspondingly **increases** Wnt/β-catenin transcriptional
activity relative to wild-type TBL1XR1 overexpression. This is the
mechanistic explanation for why F10L behaves like a loss-of-function
allele in Mastrototaro's NCoR-dependent assay despite being a missense
variant rather than a null.

**Wei et al. 2025**, *BMC Medical Genomics* — *"Different mutations in
TBL1XR1 lead to diverse phenotypes of neurodevelopmental disorder: two
case reports."* Reports two novel missense variants (p.Val314Phe and
p.Asp463Tyr, both in the WD40 domain) with clearly distinct clinical
presentations (one with global developmental delay, intellectual
disability, and seizures consistent with a developmental epileptic
encephalopathy pattern on EEG; the other with a phenotype more consistent
with Pierpont syndrome). The paper's broader literature synthesis is the
source of a critical fact for this extension's design: of 48 reported
missense variants, 43 (89.5%) fall within the WD40 domain — and *even a
single residue* can produce opposite phenotypic classifications depending
on the exact substitution (the paper cites S447N producing a partial
Pierpont phenotype while S447R at the same position does not). This is
direct evidence against ever assuming a default phenotype/direction for
an uncharacterized WD40 missense variant purely by its location in the
protein.

**Nagy et al. 2024**, *Orphanet Journal of Rare Diseases* — *"The spectrum
of neurological presentation in individuals affected by TBL1XR1 gene
defects."* A cross-sectional caregiver survey (n=41 patients, recruited
through two disease-specific Facebook support groups, REDCap-based data
collection, IRB-approved through Massachusetts General Hospital) providing
the field's first systematic natural-history data. Seizures occurred in
34.1% of patients (median onset age 2.5 years, range 0.4–23.2 years,
types including absence, tonic-clonic, focal, tonic, and epileptic
spasms), with the large majority occurring before age 10. Language was the
most severely affected developmental domain; a minority of patients
experienced developmental regression, most commonly in language, extending
into the second decade of life. Genetically, every reported missense
variant and in-frame deletion fell within the WD40 repeat region, and the
paper explicitly concludes that genotype does not reliably predict
phenotype — patients diagnosed with Pierpont syndrome carried a
substantially wider range of variants than the single canonical Y446C
mutation historically associated with it. The paper also reports that
TBL1XR1 is highly intolerant of true loss-of-function variants in the
general population (LOEUF score 0.11, zero LOF alleles observed in
gnomAD), consistent with most clinically observed patients carrying
missense or partial-function alleles rather than true null alleles.

## 2.3 Directory structure

```
tbl1xr1/
├── tbl1xr1_types.py       Five hardcoded genotype profiles — LOF, F10L
│                          (schizophrenia), G70N (West syndrome), L282P
│                          (ASD), Y446C (Pierpont) — each carrying its
│                          literature source, evidence tier, and predicted
│                          direction (or explicit "unknown")
├── tbl1xr1_metrics.py     Primary metric (mean per-unit firing rate — see
│                          §2.4 for why this and not network burst count)
│                          and isogenic-pair hypothesis comparison
│                          (Wilcoxon signed-rank + duration-confound check)
├── tbl1xr1_pipeline.py    Orchestrator: wires the above into REAL,
│                          unmodified Axon functions (QC, network topology,
│                          network bursts, longitudinal comparison) rather
│                          than reimplementing any of them
└── tbl1xr1_synthetic.py   Synthetic spike-train generator calibrated to
                           Mastrototaro 2021's actual published numbers
                           (WT 4.0 ± 0.3 Hz, KO 6.3 ± 0.4 Hz, n=19 cells)

tests/
├── test_tbl1xr1_pipeline.py               11 assertions — hypothesis-test logic
└── test_tbl1xr1_pipeline_entegrasyon.py    6 assertions — real-Axon integration,
                                            including the discovered limitation
                                            described in §2.6
```

## 2.4 Design decision: why the primary metric is per-unit firing rate, not network burst count

Mastrototaro's only direct electrophysiological finding is a **cell-level,
synaptic-current-frequency** measurement (sEPSC frequency via patch-clamp).
The nearest conceptual equivalent obtainable from MEA data is **mean
per-unit spike rate**, not population-level network burst count — those
are different quantities measuring different things (a burst is a
population-synchrony event; sEPSC frequency is a per-cell excitability
measure).

This was not left as a theoretical argument — it was tested directly.
Synthetic data calibrated to Mastrototaro's real numbers was passed through
**both** metrics: through the new `ortalama_unit_ateşleme_hizi()` function,
and through Axon's real, unmodified `organoid_units_analiz.network_burst_tespit()`.
The per-unit rate metric correctly recovered the calibrated signal (≈4.0 Hz
control vs. ≈6.3 Hz mutant). The network burst count returned **0.0 bursts
per minute for both groups**, because Axon's burst detector uses a
threshold that is adaptive *relative to each recording's own median
activity level* — it is designed to detect discrete, punctuated bursts
against a quiet baseline, and it is structurally unable to detect a
uniformly, continuously elevated firing rate that never dips to a quiet
baseline within the same recording. Had network burst count been chosen as
the primary metric without this check, the entire calibrated signal would
have been silently invisible to the pipeline.

Network burst count is still computed by the pipeline (as a labeled,
secondary/exploratory metric, via the real Axon function, untouched), and
so is a network-topology metric (STTC-matrix-derived small-world sigma,
via the real `graph_analiz.py`, also untouched, as a labeled
tertiary/exploratory metric) — but neither enters the hypothesis test.

## 2.5 Design decision: five genotypes, three of them deliberately "unknown"

`tbl1xr1_types.py` hardcodes five mutation classes rather than using a
generic, open-ended registry, because — as of this literature review —
we actually know the relevant facts for exactly five classes, and
guessing beyond them would misrepresent the state of the evidence:

| Genotype | Clinical correlate | Shares LOF's NPC mechanism? (Mastrototaro) | Predicted firing-rate direction | Evidence tier |
|---|---|---|---|---|
| LOF (true deletion/frameshift) | MRD41 — autism, ID, seizures present | Yes (defining case) | **Increase** | Direct mouse patch-clamp electrophysiology |
| F10L | Sporadic schizophrenia | Yes (only rescue failure) | **Increase** (inferred via shared mechanism + Wnt activation) | Complementation rescue-failure + biochemistry |
| G70N | West syndrome, autistic features | **No** — rescued normally | Unknown | Clinical EEG only, no functional data |
| L282P | Autism spectrum disorder | **No** — rescued normally | Unknown | Clinical genetics only, no functional data |
| Y446C | Pierpont syndrome, autism absent | **No** — rescued normally | Unknown | Cell-culture biochemistry (complex assembly), no electrophysiology |

For any as-yet-uncharacterized WD40 missense variant not matching one of
the four specific mutations above, `profil_getir()` deliberately **raises
an exception rather than returning a guessed profile** — directly because
Wei et al. 2025 documented that even a single residue's identity can flip
the phenotype (S447N vs. S447R). Any caller reaching this genotype is
required to treat it as an open discovery question, not a
confirmation/rejection test.

## 2.6 A real limitation, found by testing, not assumed

`organoid_compare.ikili_karsilastir()` — Axon's real, unmodified pairwise
comparison function — was wired into a new longitudinal-comparison
wrapper (`boylamsal_karsilastir()`), intended to answer "at what
developmental week does a phenotype emerge in the same organoid," a
question directly motivated by Nagy et al. 2024's finding that seizure
onset age varies substantially across patients.

Testing this wrapper surfaced something worth documenting rather than
hiding. Two independently generated synthetic "sessions," both drawn from
the *identical* wild-type distribution (same mean, same generative
process, only a different random seed), were compared — and returned a
spurious p-value of 4.8 × 10⁻²⁰, a false positive of an extreme order of
magnitude for what should have been a null result. Splitting a single
continuous recording in half and comparing the two halves against each
other (same underlying cells, same continuous process, genuinely the
correct kind of "no difference" negative control) instead returned
p=0.70 — correctly non-significant.

The interpretation: `ikili_karsilastir()`'s windowed, autocorrelation-aware
permutation test is well-suited to detecting drift or trend *within* one
continuous recording, but it does not model the *additional* variance that
exists between two genuinely separate recording sessions of nominally "the
same" preparation (different electrode contact, different captured units,
etc.). Applying it naively across sessions can manufacture significance
from nothing.

`boylamsal_karsilastir()` now carries this finding directly in its
docstring and return value, and mandates that any real use of it first be
validated against a split-half negative control on each session
individually — a requirement enforced by convention (documented) rather
than by code (since Axon's own comparison function was left untouched),
and demonstrated in `tests/test_tbl1xr1_pipeline_entegrasyon.py` with both
the correct negative control (p=0.70, passing) and the known false-positive
case (p<<0.05, also passing, as an explicit regression guard against ever
silently losing track of this limitation).

## 2.7 What this extension is actually useful for

1. **It will be the first dataset to test whether Mastrototaro's mouse
   finding translates to human tissue at all.** No prior study has
   measured TBL1XR1-related electrophysiology in a human organoid.
2. **For a specific patient's variant, this may produce the first
   functional evidence in existence for that variant**, if — as is
   statistically likely given that most reported TBL1XR1 variants are
   private, single-patient findings (Nagy et al. found 26 of 32 genotyped
   patients carried previously unreported variants) — that mutation does
   not match one of the four literature-characterized classes. Functional
   evidence of this kind is the category of data used under ACMG
   variant-classification criteria (specifically PS3) to strengthen a
   variant's pathogenicity classification — this is not only
   scientifically informative but potentially clinically relevant to the
   family involved.
3. **It produces a candidate outcome measure for future treatment
   testing**, useful for any research program investigating possible
   treatment approaches for TBL1XR1-related disease.
4. **It directly answers an open question the natural-history literature
   itself poses.** Nagy et al. 2024 explicitly conclude that "it is not
   known whether functional differences caused by different variants in
   TBL1XR1 explain the phenotypic differences" and call for exactly the
   kind of systematic functional characterization ("edgotyping") this
   pipeline is built to support, however small a first contribution.

It is equally important to be explicit about what it is **not**: with a
likely sample size of one patient (absent an isogenic control line, this
project cannot even claim n=1 versus a matched control), this is a data
point, not a statistically powered study; a measured electrophysiological
difference would not by itself explain *why* it occurs (that requires the
kind of molecular assays — western blot, RNA-seq — that Mastrototaro
performed and this pipeline does not); and this is one step in a treatment
development process, not a treatment in itself.

## 2.8 Test results (all obtained by actually running this exact code)

| Suite | Result |
|---|---|
| `v6_3/test_organoid_metrics.py` (original Axon suite, post bug-fix) | 15/15 |
| `tests/test_tbl1xr1_pipeline.py` | 11/11 |
| `tests/test_tbl1xr1_pipeline_entegrasyon.py` | 6/6 |

## 2.9 Usage

```python
import sys
sys.path.insert(0, 'tbl1xr1')
from tbl1xr1_types import Tbl1xr1Genotip
from tbl1xr1_pipeline import organoid_tam_analiz, boylamsal_karsilastir
from tbl1xr1_metrics import genotip_hipotezini_test_et

# Per-organoid: QC gate + primary/secondary/tertiary metric families
sonuc = organoid_tam_analiz(spike_listesi, sure_sn, Tbl1xr1Genotip.LOF)

# Once an isogenic mutant/control pair is available:
karsilastirma = genotip_hipotezini_test_et(
    Tbl1xr1Genotip.LOF, mutant_hizlari, kontrol_hizlari)

# Longitudinal, same-organoid comparison across weeks — ONLY after running
# a split-half negative control on each session (see §2.6):
boylamsal = boylamsal_karsilastir('org1', 10, spikes_week10, sure10,
                                   14, spikes_week14, sure14)
```

## 2.10 Deliberately left undone (and why)

- **`organoid_signal.py` / `organoid_lfp.py` are not wired in.** Both
  require raw voltage traces. No raw voltage data exists yet in this
  project (only synthetic spike-time data was available to test against),
  and fabricating synthetic raw voltage traces to force a test would have
  contradicted the same principle applied everywhere else in this
  extension: never invent a number where real data or real literature
  should decide it.
- **No electrophysiological hypothesis exists for G70N, L282P, or Y446C**,
  because none exists in the literature (§2.5). This was not an oversight;
  it is the correct state of current knowledge, and the code enforces it
  by refusing to guess.
- **`organoid_db.py`, `organoid_output.py`, `batch_analiz.py`, and
  `organoid_cli.py` were left untouched and unwrapped.** These are
  general-purpose infrastructure for running Axon at scale across many
  subjects — not TBL1XR1-specific science — and will be used as-is once a
  real multi-organoid dataset exists. Rewriting them here would have
  added risk without adding scientific value.
