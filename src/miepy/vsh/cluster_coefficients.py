"""Defines function to calculate scattering coefficients of a cluster."""

import numpy as np

import miepy.coordinates as coordinates
from miepy.cpp.vsh_functions import vsh_mode

from .mode_indices import lmax_to_rmax, mode_indices, rmax_to_lmax
from .vsh_translation import vsh_translation


# TODO: equations for rmax, r, lmax (here and elsewhere) should be a function call
# TODO: iteration over (n,m,r) could be simplified through a generator call (see all interactions)
def cluster_coefficients(positions, p_scat, k, origin, lmax=None):
    """Solve for the cluster scattering coefficients of N particles around an origin.

    Arguments:
        positions[N,3]   particle positions
        p_scat[N,2,rmax]   scattering coefficients
        k                medium wavenumber
        origin           position around which to calculate the cluster coefficients
        lmax             (optional) compute scattering for up to lmax terms (default: lmax of input p/q)
    """
    Nparticles = positions.shape[0]
    rmax_in = p_scat.shape[-1]
    lmax_in = rmax_to_lmax(rmax_in)

    if lmax is None:
        lmax = lmax_in

    rmax = lmax_to_rmax(lmax)
    p_cluster = np.zeros([2, rmax], dtype=complex)

    for i in range(Nparticles):
        if np.all(positions[i] == origin):
            p_cluster[0, :rmax_in] += p_scat[i, 0, :rmax_in]
            p_cluster[1, :rmax_in] += p_scat[i, 1, :rmax_in]
            continue

        rij = origin - positions[i]
        rad, theta, phi = coordinates.cart_to_sph(*rij)

        for r, n, m in mode_indices(lmax):
            for rp, v, u in mode_indices(lmax_in):
                a = p_scat[i, 0, rp]
                b = p_scat[i, 1, rp]

                A, B = vsh_translation(m, n, u, v, rad, theta, phi, k, vsh_mode.incident)

                p_cluster[0, r] += a * A + b * B
                p_cluster[1, r] += a * B + b * A

    return p_cluster
