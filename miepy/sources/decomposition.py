"""
methods for perfomring decomposition of source into vsh expansion coefficients
"""

import numpy as np
import miepy

def sample_sphere_point_matching(position, radius, sampling):
    """Sample points on the surface of the sphere for the point matching method
       Returns points[3,N]

       Arguments:
           position[3]   position of sphere
           radius        radius of sphere
           sampling      angular points sampled per pi radians
    """

    theta = np.linspace(0, np.pi, sampling)
    phi = np.linspace(0, 2*np.pi, 2*sampling)[:-1]
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    X,Y,Z = miepy.coordinates.sph_to_cart(radius, THETA, PHI, origin=position)

    Nphi = phi.shape[0]
    return np.reshape(np.array([X,Y,Z]), [3, -1])[:,Nphi-1:-Nphi+1]

#TODO: create a root function: instead of source... pass in raw fields for more generality, (no position needed)
def point_matching(source, position, radius, k, Lmax, sampling=6):
    """Decompose a source into VSHs using the point matching method
       Returns p_src[2,rmax]
       
       Arguments:
           source      source object
           position    position around which to decompose
           radius      radius of sphere for sampling the source field
           k           medium wavenumber
           Lmax        maximum number of multipoles
           sampling    angular points sampled per pi radians (default: 5)
    """
    points = sample_sphere_point_matching(position, radius, sampling)
    Npoints = points.shape[1]
    X = points[0]
    Y = points[1]
    Z = points[2]
    _, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)

    rmax = miepy.vsh.Lmax_to_rmax(Lmax)
    E_vsh = np.zeros([3*Npoints, 2*rmax], dtype=complex)
    E_src = np.zeros([3*Npoints], dtype=complex)

    for N in range(Npoints):
        E = source.E_field(X[N], Y[N], Z[N], k)
        E = miepy.coordinates.vec_cart_to_sph(E, THETA[N], PHI[N])
        E_src[3*N:3*N+3] = E

    for i,(n,m) in enumerate(miepy.vsh.mode_indices(Lmax)):
        Nfunc, Mfunc = miepy.vsh.VSH(n, m, mode=miepy.vsh.VSH_mode.incident)
        Emn_val = miepy.vsh.Emn(m, n)
        for N in range(Npoints):
            E_vsh[3*N:3*N+3,i] = -1j*Emn_val*Nfunc(radius, THETA[N], PHI[N], k)
            E_vsh[3*N:3*N+3,rmax+i] = -1j*Emn_val*Mfunc(radius, THETA[N], PHI[N], k)

    sol = np.linalg.lstsq(E_vsh, E_src, rcond=None)
    p_src = sol[0]

    return np.reshape(p_src, [2,rmax])
