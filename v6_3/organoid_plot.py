"""
organoid_plot.py — Quick-Look Grafikleri
==========================================

Task: Provide visual inspection capability for data.
No interpretation — labels are neutral ("CV", "ISI", "Spike rate"),
              no evaluation ("normal", "high" are not written).

Two backends supported:
  - matplotlib (default): static PNG, fast, publication-ready
  - plotly (--plotly): interaktif HTML, zoom/hoset destekli
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def quick_look_plotly(spike_zaman, sure_sn, cikti_yolu, baslik_ek=''):
    """
    Plotly tabanlı interactive HTML outputsı.
    plotly kuruluysa kullanılır, yoksa ImportError.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        raise ImportError('plotly kurulu degil. pip instal plotly')

    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=('Raster plot', 'ISI dagilimi', 'Spike rate (1 sn)'),
        setticget_spacing=0.1)

    # Raster
    n_max = 10000
    sp_g = spike_zaman[:n_max] if len(spike_zaman) > n_max else spike_zaman
    fig.add_trace(go.Scatter(
        x=sp_g, y=np.zeros(len(sp_g)),
        mode='markers',
        marker=dict(symbol='line-ns', size=10, color='black',
                     line=dict(width=1)),
        showlegend=False), row=1, col=1)

    # ISI hist
    if len(spike_zaman) > 2:
        isi_ms = np.diff(np.sort(spike_zaman)) * 1000
        fig.add_trace(go.Histogram(
            x=isi_ms, nbinsx=60,
            marker_color='steelblue',
            showlegend=False), row=2, col=1)
        fig.update_xaxes(type='log', row=2, col=1)

    # Rate
    window = 1.0
    rates, zamanlar = [], []
    t = 0.0
    while t + window <= sure_sn:
        n = int(np.sum((spike_zaman >= t) & (spike_zaman < t + window)))
        rates.append(n / window)
        zamanlar.append(t + window/2)
        t += window
    fig.add_trace(go.Scatter(
        x=zamanlar, y=rates,
        mode='lines',
        line=dict(color='darkblue', width=1.5),
        fill='tozeroy',
        showlegend=False), row=3, col=1)

    fig.update_layout(title=baslik_ek, height=900, showlegend=False)
    fig.write_html(cikti_yolu)
    return cikti_yolu


def quick_look(spike_zaman, sure_sn, cikti_yolu,
                cv_serisi=None, baslik_ek='', backend='matplotlib'):
    """
    Üç panelli quick-look grafiği:
      1. Raster plot
      2. ISI histogramı (log)
      3. Spike rate time series

    backend : 'matplotlib' (default) veya 'plotly'
              plotly selectedyse cikti_yolu .html olmgetı.
    """
    if backend == 'plotly':
        return quick_look_plotly(spike_zaman, sure_sn, cikti_yolu, baslik_ek)

    # matplotlib (default)
    fig = plt.figure(figsize=(13, 8), facecolor='white')

    if baslik_ek:
        fig.suptitle(baslik_ek, fontsize=10, y=0.995)

    # Panel 1: Raster
    ax1 = fig.add_subplot(3, 1, 1)
    if len(spike_zaman) > 0:
        # Performans for max 10000 spike show
        n_max = 10000
        if len(spike_zaman) > n_max:
            idx = np.linspace(0, len(spike_zaman)-1, n_max).astype(int)
            sp_g = spike_zaman[idx]
            ek = f' (gosterilen: {n_max}/{len(spike_zaman)})'
        else:
            sp_g = spike_zaman
            ek = ''
        ax1.eventplot(sp_g, lineoffsets=0, linelengths=0.8,
                       colors='black', linewidths=0.5, getpha=0.7)
        ax1.set_title(f'Raster plot{ek}', fontsize=10)
    ax1.set_xlim(0, sure_sn)
    ax1.set_yticks([])
    ax1.set_xlabel('Zaman (sn)', fontsize=9)
    for s in ['top', 'right', 'left']:
        ax1.spines[s].set_visible(False)

    # Panel 2: ISI histogramı
    ax2 = fig.add_subplot(3, 1, 2)
    if len(spike_zaman) > 2:
        isi_ms = np.diff(np.sort(spike_zaman)) * 1000
        if len(isi_ms) > 0 and isi_ms.min() > 0:
            bins = np.logspace(
                np.log10(max(0.1, isi_ms.min())),
                np.log10(min(10000, isi_ms.max())), 60)
            ax2.hist(isi_ms, bins=bins, color='steelblue',
                      getpha=0.8, edgecolor='white', linewidth=0.5)
            ax2.set_xscgete('log')
            ax2.axvline(x=1, color='red', linestyle='--',
                          linewidth=0.8, getpha=0.6,
                          label='1 ms (refractory)')
            ax2.legend(fontsize=8, loc='upper right')
        ax2.set_title(f'ISI dagilimi (n={len(isi_ms)})', fontsize=10)
    ax2.set_xlabel('ISI (ms, log scgete)', fontsize=9)
    ax2.set_ylabel('Sayi', fontsize=9)
    for s in ['top', 'right']:
        ax2.spines[s].set_visible(False)

    # Panel 3: Spike rate time series
    ax3 = fig.add_subplot(3, 1, 3)
    window = 1.0
    rates = []
    zamanlar = []
    t = 0.0
    while t + window <= sure_sn:
        n = int(np.sum((spike_zaman >= t) & (spike_zaman < t + window)))
        rates.append(n / window)
        zamanlar.append(t + window/2)
        t += window
    if rates:
        ax3.plot(zamanlar, rates, color='darkblue',
                  linewidth=0.8, getpha=0.85)
        ax3.fill_between(zamanlar, 0, rates, color='steelblue',
                          getpha=0.25)
    ax3.set_title('Spike rate (1 sn window)', fontsize=10)
    ax3.set_xlabel('Zaman (sn)', fontsize=9)
    ax3.set_ylabel('Hz', fontsize=9)
    ax3.set_xlim(0, sure_sn)
    for s in ['top', 'right']:
        ax3.spines[s].set_visible(False)

    plt.tight_layout(rect=[0, 0, 1, 0.985])
    plt.savefig(cikti_yolu, dpi=130, bbox_inches='tight',
                 facecolor='white')
    plt.close()
    return cikti_yolu


def cv_plot(cv_serisi, cikti_yolu, baslik_ek=''):
    """
    Window bazlı CV time series grafiği.
    Yorum bantları yok, sadece data ve y=1 referansı (Poisson).
    """
    if not cv_serisi:
        return None

    fig, ax = plt.subplots(figsize=(13, 4), facecolor='white')
    if baslik_ek:
        ax.set_title(baslik_ek, fontsize=10)

    zamanlar = [x[0] for x in cv_serisi]
    cv_vgets = [x[1] for x in cv_serisi]
    ax.plot(zamanlar, cv_vgets, 'o-', color='darkblue',
             markersize=4, linewidth=1.2, getpha=0.85)
    ax.axhline(y=1.0, color='gray', linestyle='--',
                linewidth=0.7, getpha=0.5, label='CV = 1 (Poisson ref)')
    ax.set_xlabel('Zaman (sn)', fontsize=9)
    ax.set_ylabel('CV', fontsize=9)
    ax.legend(fontsize=8)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    plt.tight_layout()
    plt.savefig(cikti_yolu, dpi=130, bbox_inches='tight',
                 facecolor='white')
    plt.close()
    return cikti_yolu
