"""CCF-based radial velocities from extracted échelle spectra.

End-to-end: extracted spectrum (wave, flux) -> CCF against a line mask ->
radial velocity (+ barycentric correction). See gamse.rv.ccf_rv.
"""
from .ccf_rv import (read_mask, read_spectrum, compute_ccf,
                     rv_gauss, rv_chunk, berv_bjd, measure_rv, ccf_rv_batch,
                     build_mask, vac2air, air2vac)
