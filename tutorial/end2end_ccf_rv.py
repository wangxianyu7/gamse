#!/usr/bin/env python3
"""
End-to-end CCF radial velocities — worked example
=================================================

Extracted spectrum  ->  CCF against a line mask  ->  barycentric RV.

Data: the RV standard star **HD 9407** (G6.5 V) observed on three nights with
the Xinglong 2.16 m HRS (2026-10-18, 10-26, 10-27).  Because a standard star has
a constant true RV, the scatter of the measured RVs *is* the instrument's RV
stability.  This reproduces the three-timescale result:

    within a night        : a few - tens of m/s
    adjacent nights (1 d)  : ~ 2 m/s   (10-26 vs 10-27)
    across ~8 days         : ~ 300-500 m/s common zero-point drift (10-18)

Run:   python end2end_ccf_rv.py
"""
import glob
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

# gamse.rv — the module this tutorial demonstrates
from gamse.rv import ccf_rv_batch, read_mask

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

# HD 9407 coordinates (used for the barycentric correction)
coord = SkyCoord('01:34:33.0 +68:56:53', unit=(u.hourangle, u.deg))

files = sorted(glob.glob(os.path.join(DATA, 'hd9407_*.fits')))
mask = read_mask(os.path.join(DATA, 'mask_GK.txt'))   # 4300-6476 A, O2 excluded

# One call: extracted spectra -> barycentric RVs.
# Xinglong DATE-OBS is Beijing local time, so utc_offset_hours=8.
#
#   method='gauss' : Gaussian centroid of the mask CCF, absolute RV + BERV.
#                    Best for sharp-lined, decent-S/N stars (like this G dwarf);
#                    reproduces the LSD-Gaussian benchmark to a few m/s.
#   method='chunk' : robust chunk cross-correlation with a data-driven error;
#                    preferred for broad-lined / fast rotators / low S/N, where a
#                    single Gaussian on the full CCF is biased.
tab = ccf_rv_batch(files, mask, coord=coord, method='gauss', utc_offset_hours=8.0)

# Report RV relative to the overall median (standard star -> anomaly ~ 0)
tab['rv_ms'] -= np.median(tab['rv_ms'])

print('\n  file                          BJD_TDB        RV[m/s]   err[m/s]  BERV[km/s]')
night, nightkey = [], []
for r in tab:
    nk = r['file'].split('_')[1]          # night tag yyyymmdd
    nightkey.append(nk)
    print('  %-28s %.5f  %+8.1f  %7.1f  %+8.3f'
          % (r['file'], r['bjd'], r['rv_ms'], r['err_ms'], r['berv_kms']))
nightkey = np.array(nightkey)

print('\n  --- per-night mean RV (instrument zero-point) ---')
means = {}
for nk in ['20261018', '20261026', '20261027']:
    m = nightkey == nk
    if m.any():
        means[nk] = np.mean(tab['rv_ms'][m])
        print('  %s : %+8.1f m/s  (scatter %.0f m/s, N=%d)'
              % (nk, means[nk], np.std(tab['rv_ms'][m]), m.sum()))

if {'20261026', '20261027'} <= means.keys():
    print('\n  adjacent nights 10-26 vs 10-27 : %+.0f m/s  (rock-stable)'
          % (means['20261027'] - means['20261026']))
if {'20261018', '20261027'} <= means.keys():
    print('  ~8-day drift   10-18 vs 10-27 : %+.0f m/s  (needs standard-star tie)'
          % (means['20261027'] - means['20261018']))

# ---- plot ----
colors = {'20261018': 'tab:red', '20261026': 'tab:green', '20261027': 'tab:blue'}
labels = {'20261018': '2026-10-18', '20261026': '2026-10-26', '20261027': '2026-10-27'}
fig, ax = plt.subplots(figsize=(8, 5))
for nk in colors:
    m = nightkey == nk
    if not m.any():
        continue
    x = (tab['bjd'][m] - tab['bjd'].min())
    ax.errorbar(x, tab['rv_ms'][m], tab['err_ms'][m], fmt='o', ms=8, capsize=3,
                color=colors[nk], label=labels[nk])
ax.axhline(0, color='gray', lw=0.5)
ax.set_xlabel('BJD - %.0f (days)' % tab['bjd'].min())
ax.set_ylabel('CCF RV - median (m/s)')
ax.set_title('HD 9407 (RV standard) — gamse.rv CCF velocities over 3 nights\n'
             'adjacent nights stable to ~m/s; ~8-day zero-point drift ~hundreds m/s')
ax.legend()
out = os.path.join(HERE, 'hd9407_3night_rv.png')
fig.tight_layout()
fig.savefig(out, dpi=140)
print('\n  saved figure ->', out)
