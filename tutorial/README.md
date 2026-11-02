# Tutorial — end-to-end CCF radial velocities

This tutorial shows the new `gamse.rv` layer: **an extracted spectrum in, a
barycentric radial velocity out**, via cross-correlation against a line mask.
It sits on top of GAMSE's reduction (overscan → bias → trace → flat →
extraction); feed it the 1-D spectra GAMSE produces.

## The example: an RV standard star over three nights

`data/` holds the RV standard **HD 9407** (G6.5 V) extracted from Xinglong
2.16 m HRS frames on three nights (2026-10-18, 10-26, 10-27), plus the official
**ESPRESSO G8 CCF mask** (`ESPRESSO_G8.fits`, from
[`j-faria/iCCF`](https://github.com/j-faria/iCCF); G2/G9 also included). The
ESPRESSO/HARPS masks are the community standard (also the basis of NEID-type CCF
pipelines) and are in **air** wavelengths; `read_mask` reads the FITS table
directly and keeps 4300–6476 Å, dropping the 6270–6359 Å O₂ band. G8 matches our
G6.5 V target.

Each spectrum is a small FITS file with a `WAVE / FLUX / ERR` table and the
header keywords the RV step needs (`DATE-OBS`, `EXPTIME`, `RA`, `DEC`). Note the
Xinglong `DATE-OBS` is **Beijing local time**, so the example passes
`utc_offset_hours=8`.

> **How these spectra were extracted, and important caveats** (double-fiber flat
> not applied → continuum via median filter; hard-wired aperture/pixel offsets;
> fiber-B at trace−14 pix) are documented in **[`REDUCTION_NOTES.md`](REDUCTION_NOTES.md)**.
> Short version: fine for *differential CCF RV* (demonstrated here), not for
> absolute flux / line depths until the flat is fixed.

Because a standard star has a constant true RV, the scatter of the measured RVs
*is* the instrument's velocity stability.

## Run it

Open the notebook [`end2end_ccf_rv.ipynb`](end2end_ccf_rv.ipynb) (it has the
executed outputs and the figure inline):

```bash
cd tutorial
jupyter lab end2end_ccf_rv.ipynb      # or: jupyter notebook
```

Output (a table + an inline figure, also saved as `hd9407_3night_rv.png`):

```
per-night mean RV (instrument zero-point, ESPRESSO-G8 mask)
  2026-10-18 : +344 m/s   (scatter ~39 m/s)
  2026-10-26 :   -3 m/s   (scatter ~7 m/s)
  2026-10-27 :   -6 m/s   (scatter ~6 m/s)

  adjacent nights 10-26 vs 10-27 :  ~ few m/s   (rock-stable)
  ~8-day drift   10-18 vs 10-27  :  ~ -349 m/s  (common zero-point drift)
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
                   read_mask('data/ESPRESSO_G8.fits'),
                   coord=SkyCoord('01:34:33.0 +68:56:53', unit=(u.hourangle, u.deg)),
                   method='gauss', utc_offset_hours=8.0)
# fields: file, bjd, rv_ms, err_ms, berv_kms
```

The RVs are relative to the common template — subtract the median for a per-epoch
anomaly, or tie nights with a standard star (as above). This module stops at RV;
orbit / Rossiter-McLaughlin modelling is left to downstream tools.

## Masks

The tutorial uses the ESPRESSO masks (`data/ESPRESSO_G2/G8/G9.fits`) — pick the
spectral type closest to your target. They are **air** wavelengths, so no
conversion is needed. `gamse.rv` also provides helpers to roll your own matched
mask: `build_mask(wave, flux, ...)` detects lines from any reference spectrum (a
public atlas, or a co-add of your own frames), and `vac2air`/`air2vac` align the
wavelength frame (an atlas in vacuum vs an air-based ThAr solution differ by
~80 km/s in the CCF).
