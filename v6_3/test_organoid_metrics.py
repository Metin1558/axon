"""
test_organoid_metrics.py
=========================
Axon v6.3 için gerçek, çalıştırılabilir doğrulama testleri.
Bilinen sentetik veri + bilinen beklenen sonuç mantığıyla yazıldı.
"""
import sys
sys.path.insert(0, '/tmp/axon_final/v6_3')
import numpy as np
from organoid_metrics import (
    isi_metrikleri, aktiflik_window, sttc_iki_channel,
    coklu_channel_senkron, mannwhitney_iki_grup, anova_n_grup,
    stim_oncesi_sirasi_sonrasi,
)

passed = 0
failed = 0
results = []

def check(name, condition, detail=""):
    global passed, failed
    if condition:
        passed += 1
        results.append(f"PASS: {name}")
    else:
        failed += 1
        results.append(f"FAIL: {name}  {detail}")

# ── T1: isi_metrikleri — perfectly regular spike train ──
# spikes at 0,1,2,3,4 -> ISIs all = 1.0 -> mean=1, std=0, cv=0
sp = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
r = isi_metrikleri(sp)
check("T1a: mean ISI == 1.0", abs(r['ort_isi'] - 1.0) < 1e-9, f"got {r['ort_isi']}")
check("T1b: std ISI == 0.0 (perfectly regular)", abs(r['std_isi']) < 1e-9, f"got {r['std_isi']}")
check("T1c: CV == 0.0 (regular train)", abs(r['cv']) < 1e-9, f"got {r['cv']}")

# ── T2: isi_metrikleri — insufficient data (n<2) returns None dict ──
r2 = isi_metrikleri(np.array([1.0]))
check("T2: insufficient data returns None fields", r2['ort_isi'] is None)

# ── T3: aktiflik_window — spikes concentrated in first window only ──
# 30s windows, min_spike=5. Put 10 spikes in [0,30), 0 spikes in [30,60)
sp3 = np.concatenate([np.linspace(0, 29, 10), np.array([])])
ratio, active, total = aktiflik_window(sp3, sure_sn=60, window_sn=30, min_spike=5)
check("T3: exactly 1 of 2 windows active", active == 1 and total == 2, f"got active={active} total={total}")

# ── T4: sttc_iki_channel — IDENTICAL spike trains -> STTC should be high (~1) ──
sp_a = np.sort(np.array([1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]))
sttc_identical = sttc_iki_channel(sp_a, sp_a.copy(), sure_sn=35, dt=0.01)
check("T4: identical spike trains -> STTC > 0.9", sttc_identical > 0.9, f"got {sttc_identical}")

# ── T5: sttc_iki_channel — completely disjoint, far-apart trains -> STTC near 0 or negative ──
sp_b1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
sp_b2 = np.array([100.0, 101.0, 102.0, 103.0, 104.0])
sttc_disjoint = sttc_iki_channel(sp_b1, sp_b2, sure_sn=200, dt=0.01)
check("T5: disjoint far-apart trains -> STTC near 0", abs(sttc_disjoint) < 0.15, f"got {sttc_disjoint}")

# ── T6: coklu_channel_senkron — 3 channels, matrix diagonal should be 1.0 ──
chans = [sp_a, sp_a.copy(), sp_b1]
multi = coklu_channel_senkron(chans, sure_sn=35, dt=0.01)
diag_ok = all(abs(multi['sttc_matrix'][i, i] - 1.0) < 1e-9 for i in range(3))
check("T6: STTC matrix diagonal == 1.0", diag_ok)
check("T6b: n_channels correct", multi['n_channels'] == 3)

# ── T7: mannwhitney_iki_grup — clearly different groups -> significant p-value ──
g1 = [1, 2, 1, 2, 1, 2, 1, 2]
g2 = [100, 101, 100, 101, 100, 101, 100, 101]
mw = mannwhitney_iki_grup(g1, g2)
check("T7: clearly different groups -> p < 0.05", mw['p_deger'] < 0.05, f"got p={mw['p_deger']}")

# ── T8: mannwhitney_iki_grup — insufficient data -> None ──
mw2 = mannwhitney_iki_grup([1], [2])
check("T8: insufficient data (n=1) returns None", mw2 is None)

# ── T9: anova_n_grup — 3 clearly different groups -> significant p-value ──
a1 = [1, 2, 1, 2]
a2 = [50, 51, 50, 51]
a3 = [200, 201, 200, 201]
av = anova_n_grup(a1, a2, a3)
check("T9: 3 clearly different groups -> p < 0.05", av['p_deger'] < 0.05, f"got p={av['p_deger']}")

# ── T10: stim_oncesi_sirasi_sonrasi — known before/during/after counts ──
sp10 = np.array([1.0, 2.0, 5.0, 5.5, 6.0, 12.0, 13.0])  # before=2, during(5-7)=3, after=2
r10 = stim_oncesi_sirasi_sonrasi(sp10, sure_sn=15, stim_bas=5.0, stim_bit=7.0)
check("T10a: before-count correct", r10['once_n'] == 2, f"got {r10['once_n']}")
check("T10b: during-count correct", r10['sirasinda_n'] == 3, f"got {r10['sirasinda_n']}")
check("T10c: after-count correct", r10['sonra_n'] == 2, f"got {r10['sonra_n']}")

# ── Summary ──
print("\n".join(results))
print(f"\n{'='*50}")
print(f"TOPLAM: {passed}/{passed+failed} test geçti")
print(f"{'='*50}")
sys.exit(0 if failed == 0 else 1)
