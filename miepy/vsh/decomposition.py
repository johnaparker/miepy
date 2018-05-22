"""
Decomposition of electric fields and sources into VSH coefficients using:
    (1) A point matching method
    (2) An integral projection method
"""

import numpy as np
from miepy import vsh, coordinates

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

def sample_plane_point_matching(position, size, sampling):
    """Sample points on a planar surface (z-oriented) for the point matching method
       Returns points[3,N]

       Arguments:
           position[3]   center of plane
           size          length of the plane (square shape)
           sampling      number of sampling points along a dimension
    """
    x = position[0] + np.linspace(-size/2, size/2, sampling)
    y = position[1] + np.linspace(-size/2, size/2, sampling)
    X, Y = np.meshgrid(x, y)
    Z = positio[2]*np.ones_like(X)

    points = np.array([X, Y, Z])
    return np.reshape(points, [3,-1])

#TODO: create a root function: instead of source... pass in raw fields for more generality, (no position needed)
#TODO: sampling default based on value of Lmax
def far_field_point_matching(source, position, radius, k, Lmax, sampling=6):
    """Decompose a source into VSHs using the point matching method in the far field
       Returns p_src[2,rmax]
       
       Arguments:
           source      source object
           position    position around which to decompose
           radius      radius of sphere for sampling the source field
           k           medium wavenumber
           Lmax        maximum number of multipoles
           sampling    angular points sampled per pi radians (default: 5)
    """
    position = np.asarray(position)

    points = sample_sphere_point_matching(position, radius, sampling)
    Npoints = points.shape[1]
    X = points[0]
    Y = points[1]
    Z = points[2]
    _, THETA, PHI = coordinates.cart_to_sph(X, Y, Z, origin=position)

    rmax = vsh.Lmax_to_rmax(Lmax)

    E_src = np.zeros([2, Npoints], dtype=complex)
    E_vsh = np.zeros([2, Npoints, 2, rmax], dtype=complex)

    # fields in the upper-hemisphere are zero
    idx = THETA < np.pi/2
    E_src[:,idx] = source.spherical_ingoing(THETA[idx], PHI[idx], k)

    # phase correction for moving away from the center of the source
    rhat, *_ = miepy.coordinates.sph_basis_vectors(THETA, PHI)
    delta = source.center - position
    phase = k*np.einsum('ij,i', rhat, delta)
    E_src *= np.exp(1j*phase)

    for i,n,m in vsh.mode_indices(Lmax):
        Nfunc, Mfunc = miepy.VSH(n, m, mode=miepy.VSH_mode.ingoing)
        Emn_val = miepy.vsh.Emn(m, n)
        E_vsh[...,0,i] = -1j*Emn_val*Nfunc(radius, THETA, PHI, k)[1:]
        E_vsh[...,1,i] = -1j*Emn_val*Mfunc(radius, THETA, PHI, k)[1:]

    column = E_src.reshape([2*Npoints])
    matrix = E_vsh.reshape([2*Npoints, 2*rmax])
    sol = np.linalg.lstsq(matrix, column, rcond=None)
    p_src = sol[0]

    return np.reshape(p_src, [2,rmax])

def near_field_point_matching(source, position, size, k, Lmax, sampling):
    """Decompose a source into VSHs using the point matching method in the near field
       Returns p_src[2,rmax]
       
       Arguments:
           source      source object
           position    position around which to decompose
           size        size of xy planar region to perform point matching over
           k           medium wavenumber
           Lmax        maximum number of multipoles
           sampling    number of sampling points along a dimension
    """
    points = sample_plane_point_matching(position, size, sampling)
    Npoints = points.shape[1]
    X = points[0]
    Y = points[1]
    Z = points[2]
    RAD, THETA, PHI = coordinates.cart_to_sph(X, Y, Z, origin=position)

    rmax = vsh.Lmax_to_rmax(Lmax)

    E_src = source.E_field(X, Y, Z, k)[:2]
    H_src = source.H_field(X, Y, Z, k)[:2]
    # TODO: is this true?
    # H_src = E_src[::-1]
    E_vsh = np.zeros([2, Npoints, 2, rmax], dtype=complex)

    for i,n,m in vsh.mode_indices(Lmax):
        Nfunc, Mfunc = vsh.VSH(n, m, mode=vsh.VSH_mode.incident)
        Emn_val = vsh.Emn(m, n)
        E_vsh[...,0,i] = -1j*Emn_val*coordinates.vec_sph_to_cart(Nfunc(RAD, THETA, PHI, k), THETA, PHI)[:2]
        E_vsh[...,1,i] = -1j*Emn_val*coordinates.vec_sph_to_cart(Mfunc(RAD, THETA, PHI, k), THETA, PHI)[:2]

    H_vsh = -1j*E_vsh[...,::-1,:]

    column = np.concatenate([E_src.reshape([-1]), H_src.reshape([-1])])
    matrix = np.concatenate([E_vsh.reshape([-1, 2*rmax]), H_vsh.reshape([-1, 2*rmax])])
    sol = np.linalg.lstsq(matrix, column, rcond=None)
    p_src = sol[0]

    return np.reshape(p_src, [2,rmax])

def integral_project_fields_onto(E, r, k, ftype, n, m, mode=vsh.VSH_mode.outgoing, spherical=False):
    """Project fields onto a given mode using integral method

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
    integrated = vsh.misc.simps_2d(tau, phi, proj_data)

    return factor*integrated/norm

def integral_project_source_onto(src, k, ftype, n, m, origin=[0,0,0], sampling=30, mode=vsh.VSH_mode.incident):
    """Project source object onto a given mode using integral method

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
    r = 50e-9

    THETA, PHI = coordinates.sphere_mesh(sampling)
    X,Y,Z = coordinates.sph_to_cart(r, THETA, PHI, origin=origin)
    E = src.E_field(X, Y, Z, k)

    return project_fields_onto(E, r, k, ftype, n, m, mode, spherical=False)

def integral_project_fields(E, r, k, Lmax, mode=vsh.VSH_mode.outgoing, spherical=False):
    """Decompose fields into the VSHs using integral method
    Returns p[2,rmax]

    Arguments:
        E[3,Ntheta,Nphi]   electric field values on the surface of a sphere
        r                  radius
        k                  wavenumber
        Lmax               maximum number of multipoles
        mode: VSH_mode     type of VSH (outgoing, incident) (default: outgoing)
        spherical          If true, E should be in spherical components (default: False (cartesian))
    """

    rmax = vsh.Lmax_to_rmax(Lmax)
    p = np.zeros([2,rmax], dtype=complex)

    for i,n,m in vsh.mode_indices(Lmax):
        p[0,i] = project_fields_onto(E, r, k, 'electric', n, m, mode, spherical)
        p[1,i] = project_fields_onto(E, r, k, 'magnetic', n, m, mode, spherical)

    return p

def integral_project_source(src, k, Lmax, origin=[0,0,0], sampling=30, mode=vsh.VSH_mode.incident):
    """Decompose a source object into VSHs using integral method
    Returns p[2,rmax]

    Arguments:
        src        source object
        k          wavenumber
        Lmax       maximum number of multipoles
        origin     origin around which to perform the expansion (default: [0,0,0])
        sampling   number of points to sample between 0 and pi (default: 30)
        mode: VSH_mode       type of VSH (outgoing, incident) (default: incident)
    """

    rmax = vsh.Lmax_to_rmax(Lmax)
    p = np.zeros([2,rmax], dtype=complex)

    for i,n,m in vsh.mode_indices(Lmax):
        p[0,i] = project_source_onto(src, k, 'electric', n, m, origin, sampling, mode)
        p[1,i] = project_source_onto(src, k, 'magnetic', n, m, origin, sampling, mode)

    return p

