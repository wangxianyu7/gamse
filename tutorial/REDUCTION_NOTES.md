# How the tutorial spectra were produced — and important caveats

The `gamse.rv` layer takes **already-extracted** 1-D spectra and returns CCF
radial velocities. The HD 9407 FITS in `data/` were extracted from raw Xinglong
2.16 m HRS frames with the steps and **work-arounds** below. Read this before
trusting absolute fluxes or applying the recipe to new data.

## Extraction recipe used for `data/*.fits`

1. **Bias**: median of the night's `*001`–`*010` overscan-corrected frames.
2. **Trace**: reused a Sep-2026 trace (`trace_A.trc`); no per-night re-trace.
3. **Aperture**: `extract_aperset(..., lower_limit=4, upper_limit=4)` — a **±4 pix**
   window (narrower than the ±6 default; see "pixel offsets" below).
4. **Wavelength**: a **cached** `wlcalib_202412210039.fits`, mapped to the trace
   apertures with an offset of **`APOFF = +6`** (trace aperture *a* → wavecal
   aperture *a*+6). No per-night ThAr re-solve.
5. **Continuum**: per order, divide by a `median_filter(flux, 151)` response.
   This is the normalisation; see the flat caveat next.

## ⚠ Caveat 1 — the double-fiber flat is NOT applied

The current Xinglong double-fiber flat step (`get_doublefiber_flat`) does **not**
produce usable flats in this pipeline (it returns nothing end-to-end), so these
spectra were extracted **without flat-fielding / blaze correction**. The
continuum was instead set per order by a median filter.

What this means:
- **For relative/ CCF RV: little impact.** The per-order median-filter
  normalisation removes the blaze shape, and CCF RV is differential, so the
  missing flat largely cancels. The tutorial demonstrates this empirically —
  HD 9407 repeats to ~m/s between adjacent nights.
- **For absolute flux, line depths, S/N weighting, or chromatic indices: not
  reliable.** These need a real flat. Treat the extracted flux as
  continuum-normalised only.
- A proper fix (make the double-fiber flat work, divide by it before continuum
  normalisation) is left as future work and would mainly improve S/N weighting
  and chromatic diagnostics, not the first-order RV.

## ⚠ Caveat 2 — hard-wired pixel/aperture offsets

These are specific to this instrument + the cached calibrations and are **not**
auto-derived:
- `APOFF = 6` — trace-to-wavecal aperture-number offset. If you re-trace or use a
  different wavecal, re-check this (a wrong offset mis-labels every order's
  wavelength).
- **fiber-B is at `trace_A − 14 pix`** — the simultaneous-ThAr fiber. The ±4
  science aperture is deliberately narrow so it does not pick up fiber-B light; a
  wider aperture (±6/±8) dilutes the stellar CCF and adds ThAr contamination.
- The cached wavecal is from a **different date** than the science frames; it is
  adequate for RV at the ~tens-of-m/s level demonstrated here but a same-night
  ThAr solution would be better for absolute work.

## ⚠ Caveat 3 — `DATE-OBS` is Beijing local time

Xinglong headers store `DATE-OBS` in UTC+8. The tutorial passes
`utc_offset_hours=8` so the barycentric correction and BJD use the correct UTC.

## Bottom line

For the intended use — **differential CCF radial velocities** — these spectra
are fine, and the tutorial proves it against an RV standard. For anything that
depends on a true flat (absolute flux, precise line depths, per-order S/N
weighting, chromatic RV), fix the double-fiber flat first.
