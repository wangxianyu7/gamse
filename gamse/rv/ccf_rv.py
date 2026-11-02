"""
CCF-based radial velocities from extracted échelle spectra.
==========================================================

End-to-end: an **extracted spectrum** (wavelength + flux, optionally per-order)
is cross-correlated against a stellar line mask, and a radial velocity is
measured from the CCF and corrected to the barycentric frame.

This is the layer that sits on top of GAMSE's reduction (overscan -> bias ->
trace -> flat -> extraction). Feed it the 1-D extracted spectra GAMSE produces.

Design choices baked in (hard-won on Xinglong 2.16 m HRS; see tutorial/):
  * default line window 4300-6476 A, O2 telluric band 6270-6359 A excluded;
  * default RV estimator = 'chunk': split the spectrum into wavelength chunks,
    cross-correlate each against a master template, sigma-clip outlier chunks
    (tellurics / cosmics / blends), robust-mean -> RV, with a *data-driven*
    error from the chunk scatter.  This is stable on broad-lined / rotating
    stars where a single Gaussian centroid on the full CCF is biased.
  * 'gauss' (Gaussian centroid of the mask CCF) is offered for high-S/N,
    sharp-lined stars.
  * barycentric correction is always applied (astropy); report BJD_TDB.

No orbit / Rossiter-McLaughlin modelling here — this module stops at RV.
"""

import os
import gzip
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u

C_KMS = 299792.458  # speed of light, km/s

# Xinglong 2.16 m default site (override via obsloc=)
XINGLONG = dict(lat=40.3937, lon=117.5772, height=960.0)


# ---------------------------------------------------------------------------
# Line mask
# ---------------------------------------------------------------------------
def read_mask(path, wmin=4300.0, wmax=6476.0,
              exclude=((6270.0, 6359.0),), depth_min=0.1):
    """Read a stellar line mask and return (wave_A, depth) within the window.

    Supports two formats, auto-detected:
      * VALD "extract stellar" long-format list (optionally .gz): the central
        wavelength (nm) is the first number of each line entry and the depth is
        taken from the line-depth column.
      * plain two-column text:  wavelength[Angstrom]  depth

    Parameters
    ----------
    wmin, wmax : float
        Keep only lines inside [wmin, wmax] Angstrom.
    exclude : iterable of (lo, hi)
        Wavelength bands (Angstrom) to drop, e.g. the 6270-6359 O2 band.
    depth_min : float
        Minimum line depth to keep.
    """
    wl, dp = _read_mask_raw(path)
    wl = np.asarray(wl, float)
    dp = np.asarray(dp, float)
    keep = np.isfinite(wl) & np.isfinite(dp) & (wl >= wmin) & (wl <= wmax) & (dp >= depth_min)
    for lo, hi in exclude:
        keep &= ~((wl > lo) & (wl < hi))
    return wl[keep], dp[keep]


def vac2air(wave_vac):
    """Vacuum -> air wavelengths (Angstrom), Ciddor/Edlen (good in the optical).

    External templates/atlases are often in vacuum while ThAr-based wavelength
    solutions (and VALD "extract stellar" masks) are in air; mixing them shifts
    the CCF by ~80 km/s near 5000 A. Convert before building/using a mask.
    """
    wave_vac = np.asarray(wave_vac, float)
    s2 = (1e4 / wave_vac) ** 2
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    return wave_vac / n


def air2vac(wave_air):
    """Air -> vacuum wavelengths (Angstrom)."""
    wave_air = np.asarray(wave_air, float)
    s2 = (1e4 / wave_air) ** 2
    n = 1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) \
        + 0.0001599740894897 / (38.92568793293 - s2)
    return wave_air * n


def build_mask(wave, flux, depth_min=0.10, wmin=4300.0, wmax=6476.0,
               exclude=((6270.0, 6359.0),), order=5, cont_window=4000,
               normalized=False):
    """Build a line mask (wave_A, depth) by detecting absorption lines.

    Feed a high-S/N reference spectrum (a solar/stellar atlas, or a co-added
    template of your own frames). Lines are local flux minima deeper than
    ``depth_min`` after continuum normalisation.

    Parameters
    ----------
    normalized : bool
        If False (default) the continuum is estimated as a rolling 95th
        percentile (window ``cont_window`` points) and divided out first.
    order : int
        Neighbourhood (pixels) for the local-minimum test.
    Wavelengths must already be in the SAME frame (air/vacuum) as the spectra
    you will cross-correlate — see :func:`vac2air` / :func:`air2vac`.
    """
    from scipy.signal import argrelmin
    from scipy.ndimage import percentile_filter, uniform_filter1d
    wave = np.asarray(wave, float)
    flux = np.asarray(flux, float)
    if normalized:
        norm = flux
    else:
        cont = percentile_filter(flux, 95, size=cont_window)
        cont = uniform_filter1d(cont, max(1, cont_window // 2))
        norm = flux / cont
    loc = argrelmin(norm, order=order)[0]
    w0, dep = wave[loc], 1.0 - norm[loc]
    keep = (dep > depth_min) & (w0 >= wmin) & (w0 <= wmax)
    for lo, hi in exclude:
        keep &= ~((w0 > lo) & (w0 < hi))
    return w0[keep], np.clip(dep[keep], 0.0, 1.0)


def _read_mask_raw(path):
    """Return raw (wave_Angstrom, depth) arrays from a mask file.

    Supports a FITS binary table (ESPRESSO / HARPS / CORALIE / NIRPS style, with
    ``lambda``+``contrast`` columns — these are in AIR wavelengths, use as-is),
    a VALD "extract stellar" list (via specpolFlow, nm), or plain two-column text.
    """
    sp = str(path).lower()
    if sp.endswith(('.fits', '.fit')):
        with fits.open(path) as hdul:
            for h in hdul:
                if getattr(h, 'columns', None) is None:
                    continue
                names = {c.lower(): c for c in h.columns.names}
                wcol = next((names[k] for k in ('lambda', 'wave', 'wavelength', 'll')
                             if k in names), None)
                dcol = next((names[k] for k in ('contrast', 'weight', 'depth', 'w')
                             if k in names), None)
                if wcol and dcol:
                    return (np.asarray(h.data[wcol], float),
                            np.asarray(h.data[dcol], float))
        raise ValueError('no lambda/contrast table found in %s' % path)
    opener = gzip.open if sp.endswith('.gz') else open
    # Try specpolFlow first (handles VALD masks cleanly) if available.
    try:
        import specpolFlow as spf
        m = spf.make_mask(path, depthCutoff=0.0)
        return np.array(m.wl) * 10.0, np.array(m.depth)   # nm -> Angstrom
    except Exception:
        pass
    # Fallback: parse two columns of numbers (wavelength Angstrom, depth).
    wl, dp = [], []
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(('#', '!', '*')):
                continue
            parts = line.replace(',', ' ').split()
            try:
                wl.append(float(parts[0]))
                dp.append(float(parts[1]) if len(parts) > 1 else 1.0)
            except (ValueError, IndexError):
                continue
    return np.array(wl), np.array(dp)


# ---------------------------------------------------------------------------
# Spectrum I/O
# ---------------------------------------------------------------------------
def read_spectrum(path):
    """Read an extracted spectrum into a dict.

    Accepts:
      * FITS with a binary table HDU containing columns WAVE, FLUX [, ERR],
        and header keywords DATE-OBS, EXPTIME [, RA, DEC, OBJECT];
      * plain text with columns wave flux [err] (wavelength in Angstrom; if the
        median wavelength is < 1000 it is assumed to be nm and multiplied by 10).

    Returns dict with keys: wave (Angstrom), flux, err, header (dict-like).
    """
    if str(path).lower().endswith(('.fits', '.fit', '.fits.gz')):
        with fits.open(path) as hdul:
            hdr = hdul[0].header
            tab = None
            for h in hdul:
                if getattr(h, 'columns', None) is not None:
                    tab = h
                    break
            if tab is None:
                raise ValueError('no table HDU with WAVE/FLUX columns in %s' % path)
            cols = {c.upper(): c for c in tab.columns.names}
            wave = np.asarray(tab.data[cols['WAVE']], float)
            flux = np.asarray(tab.data[cols['FLUX']], float)
            err = np.asarray(tab.data[cols['ERR']], float) if 'ERR' in cols else np.full_like(flux, np.nan)
            header = dict(hdr)
    else:
        d = np.loadtxt(path)
        wave, flux = d[:, 0].astype(float), d[:, 1].astype(float)
        err = d[:, 2].astype(float) if d.shape[1] > 2 else np.full_like(flux, np.nan)
        if np.nanmedian(wave) < 1000:   # nm -> Angstrom
            wave = wave * 10.0
        header = {}
    order = np.argsort(wave)
    return dict(wave=wave[order], flux=flux[order], err=err[order], header=header)


# ---------------------------------------------------------------------------
# CCF
# ---------------------------------------------------------------------------
def compute_ccf(wave, flux, mask_wl, mask_wt, vgrid):
    """Mask-weighted cross-correlation function.

    For each velocity v in vgrid, the CCF is the depth-weighted mean of
    (1 - flux) interpolated at each mask line shifted to v.  The continuum of
    the returned CCF is subtracted using its edges.
    """
    wave = np.asarray(wave, float)
    flux = np.asarray(flux, float)
    num = np.zeros_like(vgrid, float)
    den = 0.0
    wlo, whi = wave.min() + 2, wave.max() - 2
    for w0, wt in zip(mask_wl, mask_wt):
        if w0 < wlo or w0 > whi:
            continue
        den += wt
        shifted = w0 * (1.0 + vgrid / C_KMS)
        fi = np.interp(shifted, wave, flux, left=np.nan, right=np.nan)
        num += wt * np.nan_to_num(1.0 - fi)
    if den == 0:
        raise ValueError('no mask lines fall inside the spectrum range')
    ccf = num / den
    edge = max(5, len(ccf) // 20)
    ccf = ccf - np.median(np.r_[ccf[:edge], ccf[-edge:]])
    return ccf


# ---------------------------------------------------------------------------
# RV estimators
# ---------------------------------------------------------------------------
def rv_gauss(vgrid, ccf, half_window=40.0):
    """RV = inverted-Gaussian centroid of the CCF peak (restricted window)."""
    from scipy.optimize import curve_fit
    g = lambda v, c0, amp, v0, sig: c0 + amp * np.exp(-0.5 * ((v - v0) / sig) ** 2)
    v0 = vgrid[np.argmax(ccf)]
    sel = np.abs(vgrid - v0) < half_window
    try:
        p, pc = curve_fit(g, vgrid[sel], ccf[sel],
                          p0=[0, ccf.max(), v0, 8.0], maxfev=10000)
        return p[2], float(np.sqrt(abs(pc[2, 2])))
    except Exception:
        return v0, np.nan


def _xcorr_shift(vgrid, ccf, template, max_lag=15.0, step=0.1):
    """Sub-pixel velocity lag of ccf relative to template (parabolic refine)."""
    lags = np.arange(-max_lag, max_lag + step, step)
    chi = np.array([np.sum((np.interp(vgrid, vgrid - L, ccf, 0, 0) - template) ** 2)
                    for L in lags])
    i = int(np.argmin(chi))
    if 0 < i < len(lags) - 1:
        dd = chi[i - 1] - 2 * chi[i] + chi[i + 1]
        return lags[i] + (0.5 * (chi[i - 1] - chi[i + 1]) / dd * step if dd else 0.0)
    return lags[i]


def rv_chunk(wave, flux, master_wave, master_flux, chunk_width=90.0,
             wmin=4300.0, wmax=6476.0, exclude=((6270.0, 6359.0),),
             max_lag=12.0, sigma_clip=2.5):
    """Robust chunk cross-correlation RV with a data-driven error.

    The spectrum is split into ``chunk_width`` (Angstrom) windows; each chunk is
    cross-correlated against the master template; outlier chunks are
    sigma-clipped (tellurics / cosmics / blends) and the surviving chunk RVs are
    averaged.  Returns (rv_kms, err_kms) with err = std(chunks)/sqrt(N).

    ``master_wave``/``master_flux`` should be a high-S/N template on a common
    grid (e.g. the median of all epochs; see :func:`build_master`).
    """
    edges = np.arange(wmin, wmax + chunk_width, chunk_width)
    chunks = [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]
    chunks = [c for c in chunks
              if not any(c[0] >= lo and c[1] <= hi for lo, hi in exclude)]
    step = 0.1
    lags = np.arange(-max_lag, max_lag + step, step)
    rvs = []
    for lo, hi in chunks:
        m = (master_wave >= lo) & (master_wave < hi)
        if m.sum() < 200:
            continue
        x = master_wave[m]
        b = master_flux[m] - np.nanmean(master_flux[m])
        a = np.interp(x, wave, flux, np.nan, np.nan)
        if not np.isfinite(a).all():
            continue
        a = a - np.nanmean(a)
        chi = np.array([np.sum((np.interp(x, x * (1 + L / C_KMS), a, 0, 0) - b) ** 2)
                        for L in lags])
        i = int(np.argmin(chi))
        if 0 < i < len(lags) - 1:
            dd = chi[i - 1] - 2 * chi[i] + chi[i + 1]
            rvs.append(lags[i] + (0.5 * (chi[i - 1] - chi[i + 1]) / dd * step if dd else 0.0))
        else:
            rvs.append(lags[i])
    rvs = np.array(rvs, float)
    rvs = rvs[np.isfinite(rvs)]
    if rvs.size == 0:
        return np.nan, np.nan
    for _ in range(3):
        mu, sd = np.median(rvs), np.std(rvs)
        keep = np.abs(rvs - mu) < sigma_clip * sd
        if keep.all():
            break
        rvs = rvs[keep]
    return float(np.mean(rvs)), float(np.std(rvs) / np.sqrt(len(rvs)))


def build_master(spectra, wmin=4300.0, wmax=6476.0, dwave=0.02):
    """Median template on a common Angstrom grid from a list of spectrum dicts."""
    wg = np.arange(wmin, wmax, dwave)
    stack = np.array([np.interp(wg, s['wave'], s['flux'], np.nan, np.nan)
                      for s in spectra])
    return wg, np.nanmedian(stack, axis=0)


# ---------------------------------------------------------------------------
# Barycentric correction
# ---------------------------------------------------------------------------
def berv_bjd(dateobs, coord, obsloc=XINGLONG, exptime=0.0, utc_offset_hours=0.0):
    """Barycentric RV correction (km/s) and BJD_TDB at mid-exposure.

    ``utc_offset_hours`` subtracts a local-time offset from DATE-OBS before use
    (Xinglong headers are Beijing local time -> pass 8).
    """
    loc = EarthLocation(lat=obsloc['lat'] * u.deg, lon=obsloc['lon'] * u.deg,
                        height=obsloc['height'] * u.m)
    t = Time(dateobs, scale='utc') - utc_offset_hours * u.hour + (exptime / 2.0) * u.s
    t = Time(t.iso, scale='utc', location=loc)
    berv = coord.radial_velocity_correction(obstime=t).to(u.km / u.s).value
    bjd = (t.tdb + t.light_travel_time(coord, 'barycentric')).jd
    return float(berv), float(bjd)


# ---------------------------------------------------------------------------
# High level
# ---------------------------------------------------------------------------
def measure_rv(spectrum, mask, coord=None, method='chunk', master=None,
               vspan=80.0, dv=0.25, obsloc=XINGLONG, utc_offset_hours=0.0,
               apply_berv=True):
    """Measure one barycentric RV from a single extracted spectrum.

    Parameters
    ----------
    spectrum : dict or path
        Output of :func:`read_spectrum`, or a path to read.
    mask : (wave_A, depth) tuple or path
        Line mask; a path is passed through :func:`read_mask`.
    coord : astropy SkyCoord, optional
        Target coordinates for BERV.  If None, taken from the spectrum header
        RA/DEC keywords.
    method : {'chunk', 'gauss'}
        RV estimator (see module docstring).  'chunk' needs ``master``.
    master : (wave_A, flux) tuple, optional
        Template for method='chunk'.
    """
    if not isinstance(spectrum, dict):
        spectrum = read_spectrum(spectrum)
    if not isinstance(mask, tuple):
        mask = read_mask(mask)
    mask_wl, mask_wt = mask
    hdr = spectrum['header']

    if coord is None:
        coord = _coord_from_header(hdr)

    if method == 'chunk':
        if master is None:
            raise ValueError("method='chunk' requires master=(wave, flux)")
        rv_rest, err = rv_chunk(spectrum['wave'], spectrum['flux'],
                                master[0], master[1])
    elif method == 'gauss':
        vgrid = np.arange(-vspan, vspan + dv, dv)
        ccf = compute_ccf(spectrum['wave'], spectrum['flux'], mask_wl, mask_wt, vgrid)
        rv_rest, err = rv_gauss(vgrid, ccf)
    else:
        raise ValueError('unknown method %r' % method)

    out = dict(rv_rest=rv_rest, err=err, berv=np.nan, bjd=np.nan)
    if apply_berv and coord is not None and hdr.get('DATE-OBS'):
        berv, bjd = berv_bjd(hdr['DATE-OBS'], coord, obsloc=obsloc,
                             exptime=float(hdr.get('EXPTIME', hdr.get('EXPOSURE', 0.0))),
                             utc_offset_hours=utc_offset_hours)
        out['berv'], out['bjd'] = berv, bjd
        out['rv'] = rv_rest + berv
    else:
        out['rv'] = rv_rest
    return out


def ccf_rv_batch(spec_files, mask, coord=None, method='chunk',
                 obsloc=XINGLONG, utc_offset_hours=0.0, **kw):
    """Measure barycentric RVs for a list of extracted spectra.

    Returns a structured array with fields (file, bjd, rv_ms, err_ms, berv_kms).

    For method='chunk' the spectra are first shifted to the **barycentric rest
    frame** (wavelength * (1 + BERV/c)); a common master template is built from
    these shifted spectra, and each epoch is cross-correlated against it.  The
    resulting RV is therefore already barycentric (no double counting), which is
    essential across nights where BERV changes by km/s.  RVs are relative to the
    template — subtract the median for a per-epoch anomaly, or tie nights with a
    standard star.

    For method='gauss' each epoch's absolute mask-rest RV is measured and BERV
    is added (absolute barycentric RV, directly comparable across nights).
    """
    if not isinstance(mask, tuple):
        mask = read_mask(mask)
    specs = [read_spectrum(f) for f in spec_files]
    rows = []

    if method == 'chunk':
        bary = []
        for s in specs:
            c = coord if coord is not None else _coord_from_header(s['header'])
            berv, bjd = berv_bjd(s['header']['DATE-OBS'], c, obsloc=obsloc,
                                 exptime=float(s['header'].get('EXPTIME',
                                               s['header'].get('EXPOSURE', 0.0))),
                                 utc_offset_hours=utc_offset_hours)
            bary.append(dict(wave=s['wave'] * (1.0 + berv / C_KMS),
                             flux=s['flux'], berv=berv, bjd=bjd))
        master = build_master(bary)
        for f, b in zip(spec_files, bary):
            rv, err = rv_chunk(b['wave'], b['flux'], master[0], master[1])
            rows.append((os.path.basename(f), b['bjd'], rv * 1e3,
                         err * 1e3 if np.isfinite(err) else np.nan, b['berv']))
    else:
        for f, s in zip(spec_files, specs):
            r = measure_rv(s, mask, coord=coord, method=method,
                           obsloc=obsloc, utc_offset_hours=utc_offset_hours, **kw)
            rows.append((os.path.basename(f), r['bjd'], r['rv'] * 1e3,
                         r['err'] * 1e3 if np.isfinite(r['err']) else np.nan,
                         r['berv']))
    return np.array(rows, dtype=[('file', 'U64'), ('bjd', float),
                                 ('rv_ms', float), ('err_ms', float),
                                 ('berv_kms', float)])


def _coord_from_header(hdr):
    ra, dec = hdr.get('RA'), hdr.get('DEC')
    if ra is None or dec is None:
        return None
    try:
        return SkyCoord(float(ra), float(dec), unit=u.deg)
    except (TypeError, ValueError):
        return SkyCoord(str(ra), str(dec), unit=(u.hourangle, u.deg))
