"""Cross-section computations for the JAX backend.

Pure Python/NumPy port of cpp/src/flux.cpp. These operate on small O(rmax)
arrays so plain loops are fine — no JAX JIT needed.
"""

import numpy as np


def _rmax_to_lmax(rmax):
    return int(-1 + (1 + rmax)**0.5)


def particle_cross_sections_jax(p_scat, p_inc, p_src, k):
    """Compute scattering, absorption, extinction cross-sections.

    Arguments:
        p_scat[2,rmax]  particle scattering coefficients
        p_inc[2,rmax]   particle incident coefficients
        p_src[2,rmax]   source coefficients at particle
        k               wavenumber

    Returns: (Cscat[2,lmax], Cabs[2,lmax], Cext[2,lmax])
    """
    rmax = p_scat.shape[1]
    lmax = _rmax_to_lmax(rmax)

    Cext = np.zeros((2, lmax), dtype=float)
    Cabs = np.zeros((2, lmax), dtype=float)

    factor = 4 * np.pi / k**2

    for n in range(1, lmax + 1):
        for m in range(-n, n + 1):
            r = n**2 + n + m - 1
            for a in range(2):
                idx = n - 1
                Cabs[a, idx] += factor * (np.real(np.conj(p_inc[a, r]) * p_scat[a, r])
                                          - np.abs(p_scat[a, r])**2)
                Cext[a, idx] += factor * np.real(np.conj(p_src[a, r]) * p_scat[a, r])

    Cscat = Cext - Cabs

    return Cscat, Cabs, Cext
