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
#TODO: sampling default based on value of Lmax
#TODO: radius default based on value of k
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
    _, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z, origin=position)

    rmax = miepy.vsh.Lmax_to_rmax(Lmax)

    E_src = source.E_field(X, Y, Z, k)
    H_src = source.H_field(X, Y, Z, k)
    E_src = miepy.coordinates.vec_cart_to_sph(E_src, THETA, PHI)
    H_src = miepy.coordinates.vec_cart_to_sph(H_src, THETA, PHI)

    E_vsh = np.zeros([3, Npoints, 2, rmax], dtype=complex)
    for i,(n,m) in enumerate(miepy.vsh.mode_indices(Lmax)):
        Nfunc, Mfunc = miepy.vsh.VSH(n, m, mode=miepy.vsh.VSH_mode.incident)
        Emn_val = miepy.vsh.Emn(m, n)
        E_vsh[...,0,i] = -1j*Emn_val*Nfunc(radius, THETA, PHI, k)
        E_vsh[...,1,i] = -1j*Emn_val*Mfunc(radius, THETA, PHI, k)

    H_vsh = -1j*E_vsh[...,::-1,:]

    column = np.concatenate([E_src.reshape([-1]), H_src.reshape([-1])])
    matrix = np.concatenate([E_vsh.reshape([-1, 2*rmax]), H_vsh.reshape([-1, 2*rmax])])
    sol = np.linalg.lstsq(matrix, column, rcond=None)
    p_src = sol[0]

    return np.reshape(p_src, [2,rmax])
