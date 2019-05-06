"""
Decomposition of electric fields and sources into VSH coefficients using:
    (1) A point matching method
    (2) An integral projection method
"""

import numpy as np
from miepy import vsh, coordinates

#TODO: this should be called by the point_matching methods below directly
def sampling_from_lmax(lmax, method):
    """Determine the required sampling from lmax for point matching"""
    rmax = vsh.lmax_to_rmax(lmax)
    N = max(3, int(np.ceil(rmax**0.5)))

    if method == 'far':
        return 3*N
    elif method == 'near':
        return 2*N
    else:
        raise ValueError("'{method}' is not a valid method. Use 'far' or 'near'".format(method=method))

def sample_sphere_point_matching(position, radius, sampling):
    """Sample points on the surface of the sphere for the point matching method
       Returns points[3,N]

       Arguments:
           position[3]   position of sphere
           radius        radius of sphere
           sampling      angular points sampled per pi radians
    """
    theta = np.linspace(0, np.pi, sampling+2)[1:-1]
    phi = np.linspace(0, 2*np.pi, 2*sampling + 1)[:-1]
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    X,Y,Z = coordinates.sph_to_cart(radius, THETA, PHI, origin=position)

    Nphi = phi.shape[0]
    return np.reshape(np.array([X,Y,Z]), [3, -1])

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
    Z = position[2]*np.ones_like(X)

    points = np.array([X, Y, Z])
    return np.reshape(points, [3,-1])

#TODO: create a root function: instead of source... pass in raw fields for more generality, (no position needed)
#TODO: sampling default based on value of lmax
#TODO: use far-field expressions for the vsh functions
def far_field_point_matching(source, position, radius, k, lmax, sampling=6):
    """Decompose a source into VSHs using the point matching method in the far field
       Returns p_src[2,rmax]
       
       Arguments:
           source      source object
           position    position around which to decompose
           radius      radius of sphere for sampling the source field
           k           medium wavenumber
           lmax        maximum number of multipoles
           sampling    angular points sampled per pi radians (default: 5)
    """
    position = np.asarray(position)

    points = sample_sphere_point_matching(position, radius, sampling)
    Npoints = points.shape[1]
    X = points[0]
    Y = points[1]
    Z = points[2]
    _, THETA, PHI = coordinates.cart_to_sph(X, Y, Z, origin=position)

    rmax = vsh.lmax_to_rmax(lmax)

    E_src = np.zeros([2, Npoints], dtype=complex)
    E_vsh = np.zeros([2, Npoints, 2, rmax], dtype=complex)

    # fields in the upper-hemisphere are zero
    idx = THETA >= np.pi/2
    E_src[:,idx] = source.spherical_ingoing(THETA[idx], PHI[idx], k)

    # phase correction for moving away from the center of the source
    rhat, *_ = coordinates.sph_basis_vectors(THETA, PHI)
    delta = source.center - position
    phase = k*np.einsum('ij,i', rhat, delta)
    E_src *= np.exp(1j*phase)

    for i,n,m in vsh.mode_indices(lmax):
        Nfunc, Mfunc = vsh.VSH(n, m, mode=vsh.vsh_mode.ingoing)
        Emn_val = vsh.Emn(m, n)
        E_vsh[...,0,i] = -1j*Emn_val*Nfunc(radius, THETA, PHI, k)[1:]
        E_vsh[...,1,i] = -1j*Emn_val*Mfunc(radius, THETA, PHI, k)[1:]

    column = E_src.reshape([2*Npoints])
    matrix = E_vsh.reshape([2*Npoints, 2*rmax])
    sol = np.linalg.lstsq(matrix, column, rcond=None)
    p_src = sol[0]

    return np.reshape(p_src, [2,rmax])

def near_field_point_matching(source, position, size, k, lmax, sampling):
    """Decompose a source into VSHs using the point matching method in the near field
       Returns p_src[2,rmax]
       
       Arguments:
           source      source object
           position    position around which to decompose
           size        size of xy planar region to perform point matching over
           k           medium wavenumber
           lmax        maximum number of multipoles
           sampling    number of sampling points along a dimension
    """
    points = sample_plane_point_matching(position, size, sampling)
    Npoints = points.shape[1]
    X = points[0]
    Y = points[1]
    Z = points[2]
    RAD, THETA, PHI = coordinates.cart_to_sph(X, Y, Z + 1e-9, origin=position)

    rmax = vsh.lmax_to_rmax(lmax)

    E_src = source.E_field(X, Y, Z, k)[:2]
    H_src = source.H_field(X, Y, Z, k)[:2]
    # TODO: is this true?
    # H_src = E_src[::-1]
    E_vsh = np.zeros([2, Npoints, 2, rmax], dtype=complex)

    for i,n,m in vsh.mode_indices(lmax):
        Nfunc, Mfunc = vsh.VSH(n, m, mode=vsh.vsh_mode.incident)
        Emn_val = vsh.Emn(m, n)
        E_vsh[...,0,i] = -1j*Emn_val*coordinates.vec_sph_to_cart(Nfunc(RAD, THETA, PHI, k), THETA, PHI)[:2]
        E_vsh[...,1,i] = -1j*Emn_val*coordinates.vec_sph_to_cart(Mfunc(RAD, THETA, PHI, k), THETA, PHI)[:2]

    H_vsh = -1j*E_vsh[...,::-1,:]

    column = np.concatenate([E_src.reshape([-1]), H_src.reshape([-1])])
    matrix = np.concatenate([E_vsh.reshape([-1, 2*rmax]), H_vsh.reshape([-1, 2*rmax])])
    sol = np.linalg.lstsq(matrix, column, rcond=None)
    p_src = sol[0]

    return np.reshape(p_src, [2,rmax])

def integral_project_fields_onto(E, r, k, ftype, n, m, mode=vsh.vsh_mode.outgoing, spherical=False):
    """Project fields onto a given mode using integral method

    Arguments:
        E[3,Ntheta,Nphi]     electric field values on the surface of a sphere
        r                    radius
        k                    wavenumber
        ftype                'electric' or 'magnetic'
        n                    vsh order (1, 2, ...)
        m                    vsh orientation (-n, -n+1, ..., n)
        mode: vsh_mode       type of VSH (outgoing, incident) (default: outgoing)
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

    if mode == vsh.vsh_mode.outgoing:
        factor = 1/(1j*Emn_val)
    elif mode in (vsh.vsh_mode.incident, vsh.vsh_mode.ingoing):
        factor = -1/(1j*Emn_val)
    else:
        raise ValueError('{mode} is not a valid type of mode'.format(mode=mode))

    norm = vsh.vsh_normalization_values(mode, ftype, n, m, r, k)

    proj_data  = np.sum(E*np.conj(vsh_data), axis=0)
    integrated = vsh.misc.trapz_2d(tau, phi, proj_data)

    return factor*integrated/norm

def integral_project_source_onto(src, k, ftype, n, m, origin=[0,0,0], sampling=30, mode=vsh.vsh_mode.incident):
    """Project source object onto a given mode using integral method

    Arguments:
        src        source object
        k          wavenumber
        ftype      'electric' or 'magnetic'
        n          vsh order (1, 2, ...)
        m          vsh orientation (-n, -n+1, ..., n)
        origin     origin around which to perform the expansion (default: [0,0,0])
        sampling   number of points to sample between 0 and pi (default: 30)
        mode: vsh_mode       type of VSH (outgoing, incident) (default: incident)
    """

    r = 2*np.pi/k   # choose radius to be a wavelength of the light
    r = 50e-9

    THETA, PHI = coordinates.sphere_mesh(sampling)
    X,Y,Z = coordinates.sph_to_cart(r, THETA, PHI, origin=origin)
    E = src.E_field(X, Y, Z, k)

    return project_fields_onto(E, r, k, ftype, n, m, mode, spherical=False)

def integral_project_fields(E, r, k, lmax, mode=vsh.vsh_mode.outgoing, spherical=False):
    """Decompose fields into the VSHs using integral method
    Returns p[2,rmax]

    Arguments:
        E[3,Ntheta,Nphi]   electric field values on the surface of a sphere
        r                  radius
        k                  wavenumber
        lmax               maximum number of multipoles
        mode: vsh_mode     type of VSH (outgoing, incident) (default: outgoing)
        spherical          If true, E should be in spherical components (default: False (cartesian))
    """

    rmax = vsh.lmax_to_rmax(lmax)
    p = np.zeros([2,rmax], dtype=complex)

    for i,n,m in vsh.mode_indices(lmax):
        p[0,i] = integral_project_fields_onto(E, r, k, 'electric', n, m, mode, spherical)
        p[1,i] = integral_project_fields_onto(E, r, k, 'magnetic', n, m, mode, spherical)

    return p

def integral_project_source(src, k, lmax, origin=[0,0,0], sampling=30, mode=vsh.vsh_mode.incident):
    """Decompose a source object into VSHs using integral method
    Returns p[2,rmax]

    Arguments:
        src        source object
        k          wavenumber
        lmax       maximum number of multipoles
        origin     origin around which to perform the expansion (default: [0,0,0])
        sampling   number of points to sample between 0 and pi (default: 30)
        mode: vsh_mode       type of VSH (outgoing, incident) (default: incident)
    """

    rmax = vsh.lmax_to_rmax(lmax)
    p = np.zeros([2,rmax], dtype=complex)

    for i,n,m in vsh.mode_indices(lmax):
        p[0,i] = project_source_onto(src, k, 'electric', n, m, origin, sampling, mode)
        p[1,i] = project_source_onto(src, k, 'magnetic', n, m, origin, sampling, mode)

    return p

#TODO: rad=1 .....
#TODO: (theta_min, theta_max) for integral bounds
def integral_project_source_far(src, k, lmax, sampling=20, theta_0=np.pi/2):
    """Decompose a source object into VSHs using integral method in the far-field
    Returns p[2,rmax]

    Arguments:
        src        source object
        k          wavenumber
        lmax       maximum number of multipoles
        sampling   number of points to sample between 0 and pi (default: 20)
        theta_0    integral performed from theta_0 to pi (default: pi/2)
    """
    rmax = vsh.lmax_to_rmax(lmax)

    theta = np.linspace(theta_0, np.pi, sampling)
    phi = np.linspace(0, 2*np.pi, 2*sampling)
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    rad = 1
    rhat, *_ = coordinates.sph_basis_vectors(THETA, PHI)
    Esrc = src.angular_spectrum(THETA, PHI, k)*np.exp(-1j*k)/rad

    p0 = np.zeros((2,rmax) + THETA.shape, dtype=complex)
    E0 = src.E0(k)*np.exp(1j*src.phase)

    for i,n,m in vsh.mode_indices(lmax):
        Emn_val = vsh.Emn(m, n)
        factor = E0*k**2*1j**(2-n)*np.abs(Emn_val)/(4*np.pi)
        N, M = vsh.VSH_far(n, m, vsh.vsh_mode.ingoing)

        E = N(rad, THETA, PHI, k)
        U = np.sum(Esrc*np.conjugate(E), axis=0)*np.sin(THETA)
        integrand = U*rad**2
        p0[0,i] = 2*factor*integrand

        E = M(rad, THETA, PHI, k)
        U = np.sum(Esrc*np.conjugate(E), axis=0)*np.sin(THETA)
        integrand = U*rad**2
        p0[1,i] = 2*factor*integrand

    def f(origin):
        p = np.zeros([2,rmax], dtype=complex)

        phase = k*np.einsum('i...,i', rhat, -origin)
        exp_phase = np.exp(1j*phase)

        for i in range(rmax):
            p[0,i] = vsh.misc.trapz_2d(theta, phi, p0[0,i]*exp_phase)
            p[1,i] = vsh.misc.trapz_2d(theta, phi, p0[1,i]*exp_phase)

        return p

    return f
