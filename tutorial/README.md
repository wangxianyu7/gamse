# Tutorial — end-to-end CCF radial velocities

This tutorial shows the new `gamse.rv` layer: **an extracted spectrum in, a
barycentric radial velocity out**, via cross-correlation against a line mask.
It sits on top of GAMSE's reduction (overscan → bias → trace → flat →
extraction); feed it the 1-D spectra GAMSE produces.

## The example: an RV standard star over three nights

`data/` holds the RV standard **HD 9407** (G6.5 V) extracted from Xinglong
2.16 m HRS frames on three nights (2026-10-18, 10-26, 10-27), plus a line mask
(`mask_GK.txt`, VALD, 4300–6476 Å with the 6270–6359 Å O₂ band removed).

Each spectrum is a small FITS file with a `WAVE / FLUX / ERR` table and the
header keywords the RV step needs (`DATE-OBS`, `EXPTIME`, `RA`, `DEC`). Note the
Xinglong `DATE-OBS` is **Beijing local time**, so the example passes
`utc_offset_hours=8`.

Because a standard star has a constant true RV, the scatter of the measured RVs
*is* the instrument's velocity stability.

## Run it

```bash
cd tutorial
python end2end_ccf_rv.py
```

Output (`hd9407_3night_rv.png` + a table):

```
per-night mean RV (instrument zero-point)
  2026-10-18 : +319 m/s   (scatter ~36 m/s)
  2026-10-26 :  -15 m/s   (scatter ~26 m/s)
  2026-10-27 :  -23 m/s   (scatter ~10 m/s)

  adjacent nights 10-26 vs 10-27 :  ~ few m/s   (rock-stable)
  ~8-day drift   10-18 vs 10-27  :  ~ -340 m/s  (common zero-point drift)
```

i.e. the three-timescale picture: within a night a few–tens of m/s; **adjacent
nights stable to ~m/s**; a **few-hundred-m/s common zero-point drift over ~8
days** that must be removed with standard stars before combining multi-night RVs.

## What the code does (see `gamse/rv/ccf_rv.py`)

1. `read_spectrum` — load wave/flux/err (+ header) from the FITS (or text).
2. `read_mask` — load the line mask inside the clean window (telluric band out).
3. `compute_ccf` — depth-weighted cross-correlation vs the mask.
4. measure the RV, two estimators:
   * **`gauss`** — Gaussian centroid of the mask CCF → absolute RV, then add
     BERV. Best for sharp-lined, decent-S/N stars (used here; matches the
     LSD-Gaussian benchmark to a few m/s).
   * **`chunk`** — split into wavelength chunks, cross-correlate each against a
     common (barycentric-shifted) master, σ-clip outlier chunks, robust-mean →
     RV with a *data-driven* error. Preferred for broad-lined / fast-rotating /
     low-S/N targets, where a single Gaussian on the full CCF is biased.
5. `berv_bjd` — barycentric correction and BJD_TDB (astropy), always applied.

One-call batch helper:

```python
from gamse.rv import ccf_rv_batch, read_mask
from astropy.coordinates import SkyCoord; import astropy.units as u

tab = ccf_rv_batch(sorted(glob.glob('data/hd9407_*.fits')),
                   read_mask('data/mask_GK.txt'),
                   coord=SkyCoord('01:34:33.0 +68:56:53', unit=(u.hourangle, u.deg)),
                   method='gauss', utc_offset_hours=8.0)
# fields: file, bjd, rv_ms, err_ms, berv_kms
```

The RVs are relative to the common template — subtract the median for a per-epoch
anomaly, or tie nights with a standard star (as above). This module stops at RV;
orbit / Rossiter-McLaughlin modelling is left to downstream tools.
