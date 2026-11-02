#!/usr/bin/env python3
"""
Build a Teff-matched G-star line mask from a public solar atlas.
================================================================

HD 9407 is a G6.5 V dwarf, so the **Sun** (G2 V) is an excellent, freely
available template. This script downloads the Kurucz (2005) solar irradiance
atlas, converts it from vacuum to air wavelengths, and detects absorption lines
to produce `data/mask_solar.txt` — a matched mask that performs at least as well
as a generic hotter VALD mask (and better on some nights).

Key gotcha demonstrated here: external atlases are usually in **vacuum**, while
ThAr wavelength solutions and VALD "extract stellar" masks are in **air**.
Mixing them shifts the CCF by ~80 km/s near 5000 A. We convert with
`gamse.rv.vac2air` before building the mask.

Run:  python make_solar_mask.py        (needs network; writes data/mask_solar.txt)
"""
import os
import urllib.request
import numpy as np

from gamse.rv import build_mask, vac2air

HERE = os.path.dirname(os.path.abspath(__file__))
ATLAS_URL = 'http://kurucz.harvard.edu/sun/irradiance2005/irradthu.dat'
CACHE = os.path.join(HERE, 'data', '_kurucz_solar_irradiance.dat')
OUT = os.path.join(HERE, 'data', 'mask_solar.txt')

# 1. download the solar atlas (cached)
if not os.path.exists(CACHE):
    print('downloading Kurucz 2005 solar irradiance atlas ...')
    urllib.request.urlretrieve(ATLAS_URL, CACHE)

# 2. read (wavelength[nm, vacuum], irradiance); skip the text header robustly
wl, fl = [], []
for line in open(CACHE):
    p = line.split()
    if len(p) == 2:
        try:
            wl.append(float(p[0])); fl.append(float(p[1]))
        except ValueError:
            pass
wave_vac = np.array(wl) * 10.0        # nm -> Angstrom (vacuum)
flux = np.array(fl)

# 3. vacuum -> air (so it matches the ThAr/VALD air frame)
wave_air = vac2air(wave_vac)

# 4. detect lines -> mask (continuum via rolling 95th percentile inside build_mask)
m = (wave_air > 4290) & (wave_air < 6490)
mwl, mdepth = build_mask(wave_air[m], flux[m], depth_min=0.10)

hdr = ('Solar (G2V) line mask from Kurucz 2005 solar irradiance atlas\n'
       '%s ; Thuillier+2004 calibration\n'
       'vacuum->air converted (gamse.rv.vac2air); 4300-6476 A, O2 6270-6359 removed, depth>0.10\n'
       'wave_Angstrom_AIR  line_depth' % ATLAS_URL)
np.savetxt(OUT, np.column_stack([mwl, mdepth]), fmt='%.4f %.4f', header=hdr)
print('wrote %s  (%d lines)' % (OUT, len(mwl)))
