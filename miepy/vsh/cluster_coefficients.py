"""
Defines function to calculate scattering coefficients of a cluster
"""

import numpy as np
import miepy.coordinates as coordinates
from miepy import vsh

#TODO: equations for rmax, r, Lmax (here and elsewhere) should be a function call
#TODO: iteration over (n,m,r) could be simplified through a generator call (see all interactions)
def cluster_coefficients(positions, p_scat, k, origin, Lmax=None):
    """Solve for the cluster scattering coefficients of N particles around an origin

    Arguments:
        positions[N,3]   particle positions
        p_scat[N,2,rmax]   scattering coefficients 
        k                medium wavenumber
        origin           position around which to calculate the cluster coefficients
        Lmax             (optional) compute scattering for up to Lmax terms (default: Lmax of input p/q)
    """

    Nparticles = positions.shape[0]
    Lmax_in = vsh.rmax_to_Lmax(p_scat.shape[-1])

    if Lmax is None:
        Lmax = Lmax_in

    rmax = vsh.Lmax_to_rmax(Lmax)
    p_cluster = np.zeros([2,rmax], dtype=complex)

    for i in range(Nparticles):
        if np.all(positions[i] == origin):
            p_cluster[...] = p_scat[i]
            continue

        rij = origin - positions[i]
        rad, theta, phi = coordinates.cart_to_sph(*rij)
        
        for r,(n,m) in enumerate(vsh.mode_indices(Lmax)):
            for rp,(v,u) in enumerate(vsh.mode_indices(Lmax_in)):
                a = p_scat[i,0,rp]
                b = p_scat[i,1,rp]

                A, B = vsh.vsh_translation(m, n, u, v, rad, theta, phi, k, vsh.VSH_mode.incident)

                p_cluster[0,r] += a*A + b*B
                p_cluster[1,r] += a*B + b*A

    return p_cluster
