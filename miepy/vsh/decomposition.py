"""
Decomposition of electric fields and sources into VSH coefficients using:
    (1) A point matching method
    (2) An integral projection method
"""

import numpy as np
from my_pytools.my_numpy.integrate import simps_2d
import miepy.coordinates as coordinates
from miepy import vsh

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
    X,Y,Z = coordinates.sph_to_cart(radius, THETA, PHI, origin=position)

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
    _, THETA, PHI = coordinates.cart_to_sph(X, Y, Z, origin=position)

    rmax = vsh.Lmax_to_rmax(Lmax)

    E_src = source.E_field(X, Y, Z, k)
    H_src = source.H_field(X, Y, Z, k)
    E_src = coordinates.vec_cart_to_sph(E_src, THETA, PHI)
    H_src = coordinates.vec_cart_to_sph(H_src, THETA, PHI)

    E_vsh = np.zeros([3, Npoints, 2, rmax], dtype=complex)
    for i,(n,m) in enumerate(vsh.mode_indices(Lmax)):
        Nfunc, Mfunc = vsh.VSH(n, m, mode=vsh.VSH_mode.incident)
        Emn_val = vsh.Emn(m, n)
        E_vsh[...,0,i] = -1j*Emn_val*Nfunc(radius, THETA, PHI, k)
        E_vsh[...,1,i] = -1j*Emn_val*Mfunc(radius, THETA, PHI, k)

    H_vsh = -1j*E_vsh[...,::-1,:]

    column = np.concatenate([E_src.reshape([-1]), H_src.reshape([-1])])
    matrix = np.concatenate([E_vsh.reshape([-1, 2*rmax]), H_vsh.reshape([-1, 2*rmax])])
    sol = np.linalg.lstsq(matrix, column, rcond=None)
    p_src = sol[0]

    return np.reshape(p_src, [2,rmax])

def project_fields_onto(E, r, k, ftype, n, m, mode=vsh.VSH_mode.outgoing, spherical=False):
    """Project fields onto a given mode

    Arguments:
        E[3,Ntheta,Nphi]     electric field values on the surface of a sphere
        r                    radius
        k                    wavenumber
        ftype                'electric' or 'magnetic'
        n                    vsh order (1, 2, ...)
        m                    vsh orientation (-n, -n+1, ..., n)
        mode: VSH_mode       type of VSH (outgoing, incident) (default: outgoing)
        spherical            If true, E should be in spherical components (default: False (cartesian))
    """
    Ntheta, Nphi = E.shape[1:]
    sampling = Ntheta
    THETA, PHI = coordinates.sphere_mesh(sampling)

    tau = np.linspace(-1, 1, sampling)
    phi = np.linspace(0, 2*np.pi, 2*sampling)

    N,M = vsh.VSH(n, m, mode)
    if ftype == 'electric':
        base_function = N
    elif ftype == 'magnetic':
        base_function = M

    vsh_data = base_function(r,THETA,PHI,k).squeeze()

    if not spherical:
        E = coordinates.vec_cart_to_sph(E, THETA, PHI)

    Emn_val = vsh.Emn(m, n)

    if mode == vsh.VSH_mode.outgoing:
        factor = 1/(1j*Emn_val)
    elif mode in (vsh.VSH_mode.incident, vsh.VSH_mode.ingoing):
        factor = -1/(1j*Emn_val)
    else:
        raise ValueError(f'{mode} is not a valid type of mode')

    norm = vsh.vsh_normalization_values(mode, ftype, n, m, r, k)

    proj_data  = np.sum(E*np.conj(vsh_data), axis=0)
    integrated = simps_2d(tau, phi, proj_data)

    return factor*integrated/norm

def project_source_onto(src, k, ftype, n, m, origin=[0,0,0], sampling=30, mode=vsh.VSH_mode.incident):
    """Project source object onto a given mode

    Arguments:
        src        source object
        k          wavenumber
        ftype      'electric' or 'magnetic'
        n          vsh order (1, 2, ...)
        m          vsh orientation (-n, -n+1, ..., n)
        origin     origin around which to perform the expansion (default: [0,0,0])
        sampling   number of points to sample between 0 and pi (default: 30)
        mode: VSH_mode       type of VSH (outgoing, incident) (default: incident)
    """

    r = 2*np.pi/k   # choose radius to be a wavelength of the light

    THETA, PHI = coordinates.sphere_mesh(sampling)
    X,Y,Z = coordinates.sph_to_cart(r, THETA, PHI, origin=origin)
    E = src.E_field(X, Y, Z, k)

    return project_fields_onto(E, r, k, ftype, n, m, mode, spherical=False)

def decompose_fields(E, r, k, Lmax, mode=vsh.VSH_mode.outgoing, spherical=False):
    """Decompose fields into the VSHs
    Returns p[2,rmax]

    Arguments:
        E[3,Ntheta,Nphi]   electric field values on the surface of a sphere
        r                  radius
        k                  wavenumber
        Lmax               maximum number of multipoles
        mode: VSH_mode     type of VSH (outgoing, incident) (default: outgoing)
        spherical          If true, E should be in spherical components (default: False (cartesian))
    """

    rmax = Lmax_to_rmax(Lmax)
    p = np.zeros([2,rmax], dtype=complex)

    for i,(n,m) in enumerate(mode_indices(Lmax)):
        p[0,i] = project_fields_onto(E, r, k, 'electric', n, m, mode, spherical)
        p[1,i] = project_fields_onto(E, r, k, 'magnetic', n, m, mode, spherical)

    return p

def decompose_source(src, k, Lmax, origin=[0,0,0], sampling=30, mode=vsh.VSH_mode.incident):
    """Decompose a source object into VSHs
    Returns p[2,rmax]

    Arguments:
        src        source object
        k          wavenumber
        Lmax       maximum number of multipoles
        origin     origin around which to perform the expansion (default: [0,0,0])
        sampling   number of points to sample between 0 and pi (default: 30)
        mode: VSH_mode       type of VSH (outgoing, incident) (default: incident)
    """

    rmax = Lmax_to_rmax(Lmax)
    p = np.zeros([2,rmax], dtype=complex)

    for i,(n,m) in enumerate(mode_indices(Lmax)):
        p[0,i] = project_source_onto(src, k, 'electric', n, m, origin, sampling, mode)
        p[1,i] = project_source_onto(src, k, 'magnetic', n, m, origin, sampling, mode)

    return p

