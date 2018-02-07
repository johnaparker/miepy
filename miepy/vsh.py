"""
Vector spherical harmonic related functions
"""

import numpy as np
import sympy
import scipy
from scipy import special, integrate, misc
from my_pytools.my_numpy.integrate import simps_2d
import enum
from functools import lru_cache
from math import factorial

def sph_to_cart(r, theta, phi, origin=[0,0,0]):
    """convert spherical coordinates (r, theta, phi) centered at origin to cartesian coordinates (x, y, z)"""
    x = origin[0] + r*np.sin(theta)*np.cos(phi)
    y = origin[1] + r*np.sin(theta)*np.sin(phi)
    z = origin[2] + r*np.cos(theta)

    return x,y,z

def cart_to_sph(x, y, z, origin=[0,0,0]):
    """convert cartesian coordinates (x, y, z) to spherical coordinates (r, theta, phi) centered at origin"""
    x0,y0,z0 = origin
    r = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2)**0.5
    theta = np.arccos((z - z0)/r)
    phi = np.arctan2(y - y0, x - x0)

    return r, theta, phi

def sph_basis_vectors(theta, phi):
    """obtain the spherical basis vectors (r_hat, theta_hat, phi_hat) for given theta, phi"""
    r_hat = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    theta_hat = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -1*np.sin(theta)])
    phi_hat = np.array([-1*np.sin(phi), np.cos(phi), np.zeros_like(phi)])

    return r_hat, theta_hat, phi_hat

def vec_cart_to_sph(F, theta, phi):
    """convert a vector field F from cartesian to spherical coordinates

    Arguments:
        F[3,...]     vector field values
        theta        theta coordinates
        phi          phi coordinates
    """
    Fsph = np.zeros_like(F)
    r_hat, theta_hat, phi_hat = sph_basis_vectors(theta, phi)
    Fsph[0] = np.sum(F*r_hat, axis=0)
    Fsph[1] = np.sum(F*theta_hat, axis=0)
    Fsph[2] = np.sum(F*phi_hat, axis=0)

    return Fsph

def vec_sph_to_cart(F, theta, phi):
    """convert a vector field F from spherical to cartesian coordinates

    Arguments:
        F[3,...]     vector field values
        theta        theta coordinates
        phi          phi coordinates
    """
    Fcart = np.zeros_like(F)
    r_hat, theta_hat, phi_hat = sph_basis_vectors(theta, phi)
    for i in range(3):
        Fcart[i] = F[0]*r_hat[i] + F[1]*theta_hat[i] + F[2]*phi_hat[i]

    return Fcart

def spherical_hn(n, z, derivative=False):
    """spherical hankel function of the first kind or its derivative

            n: int,array-like        order of the bessel function
            z: array[complex/float]  argument
            derivative: bool         If True, compute the derivative instead    """

    return special.spherical_jn(n,z,derivative) + 1j*special.spherical_yn(n,z,derivative)

@lru_cache(maxsize=None)
def associated_legendre(n,m, deriv=0):
    """associated legendre function of integer order and degree

            n: int         order
            m: int         degree
            deriv: int     derivative to take

        returns lpmv(x) function     """

    x = sympy.symbols('x')
    legfun_sym = (-1)**abs(m)*sympy.functions.special.polynomials.assoc_legendre(n,m,x)
    legfunc_sym_deriv = sympy.diff(legfun_sym, x, deriv)

    legfun_num_deriv = sympy.lambdify(x,legfunc_sym_deriv, modules='numpy')
    return legfun_num_deriv

@lru_cache(maxsize=None)
def pi_func(n, m):
    """pi special function that appears in the vector spherical harmonics (VSH)

            n: int         order
            m: int         degree

       returns pi(theta)                                                           """

    x = sympy.symbols('x')
    legfunc_sym = m*(-1)**abs(m)*sympy.functions.special.polynomials.assoc_legendre(n,m,x)
    legfunc_sym /= sympy.sqrt(1 - x**2)
    legfunc_sym = sympy.simplify(legfunc_sym)
    legfunc_sym = legfunc_sym.subs(x, sympy.cos(x))

    legfun_num = sympy.lambdify(x,legfunc_sym, modules='numpy')
    return legfun_num

@lru_cache(maxsize=None)
def tau_func(n, m):
    """pi special function that appears in the vector spherical harmonics (VSH)

            n: int         order
            m: int         degree

       returns tau(theta)                                                           """

    x = sympy.symbols('x')
    legfun_sym = (-1)**abs(m)*sympy.functions.special.polynomials.assoc_legendre(n,m,x)
    legfunc_sym_deriv = sympy.diff(legfun_sym, x, 1)
    legfunc_sym_deriv *= -1*sympy.sqrt(1 - x**2)
    legfunc_sym_deriv = sympy.expand(legfunc_sym_deriv)
    legfunc_sym_deriv = legfunc_sym_deriv.subs(x, sympy.cos(x))

    legfun_num_deriv = sympy.lambdify(x,legfunc_sym_deriv, modules='numpy')
    return legfun_num_deriv



@lru_cache(maxsize=None)
def a_func(m, n, u, v, p):
    """a function that appears in the VSH translation coefficients"""

    if np.abs(m) > n or np.abs(u) > v or np.abs(m+u) > p:
        return 0

    factor = (2*p+1)/2*factorial(p-m-u)/factorial(p+m+u)

    Pnm = associated_legendre(n,m)
    Pvu = associated_legendre(v,u)
    Ppmu = associated_legendre(p,m+u)
    integrand = lambda x: Pnm(x)*Pvu(x)*Ppmu(x)
    integral = integrate.quad(integrand, -1, 1)[0]

    return factor*integral

@lru_cache(maxsize=None)
def b_func(m, n, u, v, p):
    """b function that appears in the VSH translation coefficients"""

    b1 = (2*p+1)/(2*p-1)
    b2 = (v-u)*(v+u+1)*a_func(m,n,-u-1,v,p-1)
    b3 = -(p-m+u)*(p-m+u-1)*a_func(m,n,-u+1,v,p-1)
    b4 = 2*u*(p-m+u)*a_func(m,n,-u,v,p-1)

    return b1*(b2+b3+b4)

def Emn(m, n, E0):
    return E0*1j**n*(2*n+1)*factorial(n-m)/factorial(n+m)

def A_translation(m, n, u, v, r, theta, phi, k):
    factor = (-1.)**u *1.j**(v-n)*(2*v+1)/(2*v*(v+1))
    normalization = Emn(u,v,1)/Emn(m,n,1)

    sum_term = 0
    for p in range(np.abs(n-v), n+v+1):
        if np.abs(m-u) > p:
            continue
        else:
            Pnm = associated_legendre(p,m-u)
            sum_term += (-1j)**p *(n*(n+1) + v*(v+1) - p*(p+1))*a_func(m,n,-u,v,p)*spherical_hn(p, k*r)*Pnm(np.cos(theta))*np.exp(1j*(m-u)*phi)
    
    return normalization*factor*sum_term

def B_translation(m, n, u, v, r, theta, phi, k):
    factor = (-1.)**u *1.j**(v-n)*(2*v+1)/(2*v*(v+1))
    normalization = Emn(u,v,1)/Emn(m,n,1)

    sum_term = 0
    for p in range(np.abs(n-v), n+v+1):
        if np.abs(m-u) > p:
            continue
        else:
            Pnm = associated_legendre(p,m-u)
            sum_term += (-1j)**p *b_func(m,n,u,v,p)*spherical_hn(p, k*r)*Pnm(np.cos(theta))*np.exp(1j*(m-u)*phi)
    
    return normalization*factor*sum_term

class VSH_mode(enum.Enum):
    outgoing = enum.auto()
    ingoing  = enum.auto()
    incident = enum.auto()

def VSH(n, m, mode=VSH_mode.outgoing):
    """electric and magnetic vector spherical harmonic function

            n: int           order
            m: int           degree
            mode: VSH_mode   type of VSH (outgoing, incident)


       returns (N(r,θ,ϕ,k) -> [3,...], M(r,θ,ϕ,k) -> [3,...]), the 3 x,y,z components"""

    pi_f = pi_func(n,m)
    tau_f = tau_func(n,m)
    Pnm = associated_legendre(n,m)

    if mode is VSH_mode.outgoing:
        zn = spherical_hn
    elif mode is VSH_mode.incident:
        zn = special.spherical_jn
    else:
        raise TypeError('mode must be of enum type VSH_mode')
        
    def N(r, theta, phi, k):
        H = zn(n, k*r)
        Hp = zn(n, k*r, derivative=True)
        Pnm_val = Pnm(np.cos(theta))

        factor = (H + r*k*Hp)*np.exp(1j*m*phi)/(k*r)

        r_comp = n*(n+1)*Pnm_val*H/(k*r)*np.exp(1j*m*phi)
        theta_comp = tau_f(theta)*factor
        phi_comp = 1j*pi_f(theta)*factor

        return np.array([r_comp, theta_comp, phi_comp])

    def M(r, theta, phi, k):
        H = zn(n, k*r)
        factor = H*np.exp(1j*m*phi)

        theta_comp = 1j*pi_f(theta)*factor
        phi_comp = -1*tau_f(theta)*factor
        r_comp = np.zeros_like(theta_comp)

        return np.array([r_comp, theta_comp, phi_comp])

    return N,M

def vsh_normalization_values(mode, ftype, n, m, r, k):
    """Determine the norm of a given vsh mode
    
    Arguments:
        mode: VSH_mode    type of VSH (outgoing, incident)
        ftype             'electric' or 'magnetic'
        n                 vsh order (1, 2, ...)
        m                 vsh orientation (-n, -n+1, ..., n)
        r                 radius
        k                 wavenumber
    """
    if mode == VSH_mode.outgoing:
        zn = spherical_hn
    elif mode == VSH_mode.incident:
        zn = special.spherical_jn

    Emn_val = Emn(m, n, E0=1)
    zn_val = zn(n, k*r)
    angular_term = 4*np.pi*n*(n+1)/np.abs(Emn_val)

    if ftype == 'magnetic':
        radial_term = np.abs(zn_val)**2
        return angular_term*radial_term

    elif ftype == 'electric':
        znp_val = zn(n, k*r, derivative=True)
        radial_term = (np.abs(zn_val + k*r*znp_val)**2 + n*(n+1)*np.abs(zn_val)**2)/(k*r)**2
        return angular_term*radial_term

##### riccati needs to be verified with miepy #####

def riccati_1(n,z, derivative = False):
    jn = special.spherical_jn(n, z)

    if derivative:
        jn_p = special.spherical_jn(n, z, derivative=True)
        return z*jn_p + jn
    return z*jn

def riccati_2(n,z, derivative = False):
    yn = special.spherical_yn(n, z, derivative = derivative) 

    if derivative:
        yn_p = special.spherical_yn(n, z, derivative=True)
        return -z*yn_p - yn
    return -z*yn

def riccati_3(n,z, derivative = False):
    return riccati_2(n,z, derivative) - riccati_1(n,z, derivative)

###### below are pi,tau,VSH used in Mie theory, which may differ from those defined above ######

def pi_tau_func(n):
    # if np.sin(theta) == 0: return 0
    lpn = special.legendre(n)
    lpn_p = lpn.deriv()
    lpn_p2 = lpn_p.deriv()

    def pi_func(theta):
        return -1*lpn_p(np.cos(theta))
        # with np.errstate(divide='ignore', invalid='ignore'):
            # val = lpn(np.cos(theta))/np.sin(theta)
            # val[val == np.inf] = 0
            # val = np.nan_to_num(val)
            # return val

    def tau_func(theta):
        # val = -1*np.sin(theta)*lpn_p(np.cos(theta))
        val = -1*np.cos(theta)*lpn_p(np.cos(theta)) + np.sin(theta)**2*lpn_p2(np.cos(theta))
        return val

    return pi_func, tau_func 

class vector_spherical_harmonics:
    def __init__(self, n, superscript=3):
        self.pi_func, self.tau_func = pi_tau_func(n)
        self.n = n

        if superscript == 1:
            self.z_func = lambda x: special.spherical_jn(n,x)
            self.zp_func = lambda x: special.spherical_jn(n,x, derivative=True)
        elif superscript == 3:
            self.z_func = lambda x: spherical_hn(n,x)
            self.zp_func = lambda x: spherical_hn(n,x, derivative=True)

    def M_o1n(self, k):
        def f(r, theta, phi):
            theta_comp = np.cos(phi)*self.pi_func(theta)*self.z_func(k*r)
            phi_comp = -1*np.sin(phi)*self.tau_func(theta)*self.z_func(k*r)
            r_comp = np.zeros(shape = theta.shape, dtype=np.complex)
            return np.array([r_comp, theta_comp, phi_comp])
        return f

    def M_e1n(self, k):
        def f(r, theta, phi):
            theta_comp = -1*np.sin(phi)*self.pi_func(theta)*self.z_func(k*r)
            phi_comp = -1*np.cos(phi)*self.tau_func(theta)*self.z_func(k*r)
            r_comp = np.zeros(shape = theta.shape, dtype=np.complex)
            return np.array([r_comp, theta_comp, phi_comp])
        return f

    def N_o1n(self, k):
        def f(r, theta, phi):
            p = k*r
            theta_comp = np.sin(phi)*self.tau_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            phi_comp = np.cos(phi)*self.pi_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            r_comp = np.sin(phi)*self.n*(self.n+1)*np.sin(theta)*self.pi_func(theta)*self.z_func(p)/p
            return np.array([r_comp, theta_comp, phi_comp])
        return f

    def N_e1n(self, k):
        def f(r, theta, phi):
            p = k*r
            theta_comp = np.cos(phi)*self.tau_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            phi_comp = -1*np.sin(phi)*self.pi_func(theta)*(self.z_func(p) + p*self.zp_func(p))/p
            r_comp = np.cos(phi)*self.n*(self.n+1)*np.sin(theta)*self.pi_func(theta)*self.z_func(p)/p
            return np.array([r_comp, theta_comp, phi_comp])
        return f

def sphere_mesh(sampling):
    """
    Obtain a THETA,PHI mesh for discretizing the surface of the sphere, consistent
    with the format required by the project and decompose functions
    Returns (THETA,PHI) meshgrids

    Arguments:
        sampling   number of points to sample between 0 and pi
    """

    phi = np.linspace(0, 2*np.pi, 2*sampling)
    tau = np.linspace(-1, 1, sampling)
    theta = np.arccos(tau)

    THETA,PHI = np.meshgrid(theta, phi, indexing='ij')
    return THETA, PHI


def project_fields_onto(E, r, k, ftype, n, m, mode=VSH_mode.outgoing, spherical=False):
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
    THETA, PHI = sphere_mesh(sampling)

    tau = np.linspace(-1, 1, sampling)
    phi = np.linspace(0, 2*np.pi, 2*sampling)

    N,M = VSH(n, m, mode)
    if ftype == 'electric':
        base_function = N
    elif ftype == 'magnetic':
        base_function = M

    vsh_data = base_function(r,THETA,PHI,k).squeeze()

    if not spherical:
        E = vec_cart_to_sph(E, THETA, PHI)

    Emn_val = Emn(m, n, E0=1)

    if mode == VSH_mode.outgoing:
        factor = 1/(1j*Emn_val)
    elif mode == VSH_mode.incident:
        factor = -1/(1j*Emn_val)
    else:
        raise ValueError(f'{mode} is not a valid type of mode')

    norm = vsh_normalization_values(mode, ftype, n, m, r, k)

    proj_data  = np.sum(E*np.conj(vsh_data), axis=0)
    integrated = simps_2d(tau, phi, proj_data)

    return factor*integrated/norm

def project_source_onto(src, k, ftype, n, m, origin=[0,0,0], sampling=30, mode=VSH_mode.incident):
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

    THETA, PHI = sphere_mesh(sampling)
    X,Y,Z = sph_to_cart(r, THETA, PHI, origin=origin)
    E = src.E([X,Y,Z], k)

    return project_fields_onto(E, r, k, ftype, n, m, mode, spherical=False)

def decompose_fields(E, r, k, Nmax, mode=VSH_mode.outgoing, spherical=False):
    """Decompose fields into the VSHs
    Returns p[Nmax,2*Nmax+1], q[Nmax,2*Nmax+1]

    Arguments:
        E[3,Ntheta,Nphi]   electric field values on the surface of a sphere
        r                  radius
        k                  wavenumber
        Nmax               maximum number of multipoles
        mode: VSH_mode     type of VSH (outgoing, incident) (default: outgoing)
        spherical          If true, E should be in spherical components (default: False (cartesian))
    """

    p = np.zeros((Nmax,2*Nmax+1), dtype=np.complex)
    q = np.zeros((Nmax,2*Nmax+1), dtype=np.complex)
    for n in range(1, Nmax+1):
        for m in range(-n, n+1):
            p[n-1,m+n] = project_fields_onto(E, r, k, 'electric', n, m, mode, spherical)
            q[n-1,m+n] = project_fields_onto(E, r, k, 'magnetic', n, m, mode, spherical)
    return p,q

def decompose_source(src, k, Nmax, origin=[0,0,0], sampling=30, mode=VSH_mode.incident):
    """Decompose a source object into VSHs
    Returns p[Nmax,2*Nmax+1], q[Nmax,2*Nmax+1]

    Arguments:
        src        source object
        k          wavenumber
        Nmax       maximum number of multipoles
        origin     origin around which to perform the expansion (default: [0,0,0])
        sampling   number of points to sample between 0 and pi (default: 30)
        mode: VSH_mode       type of VSH (outgoing, incident) (default: incident)
    """

    p = np.zeros((Nmax,2*Nmax+1), dtype=np.complex)
    q = np.zeros((Nmax,2*Nmax+1), dtype=np.complex)
    for n in range(1, Nmax+1):
        for m in range(-n, n+1):
            p[n-1,m+n] = project_source_onto(src, k, 'electric', n, m, origin, sampling, mode)
            q[n-1,m+n] = project_source_onto(src, k, 'magnetic', n, m, origin, sampling, mode)
    return p,q

# TODO: implement origin
# TODO: E0 = 1 issue... absorb into p,q or pass it in
# TODO: implement correct factor for each mode
def expand(p, q, k, mode, origin=[0,0,0]):
    """Expand VSH coefficients to obtain an electric field function
    Returns E(x,y,z) function
    
    Arguments:
        p[Nmax,2*Nmax+1]   p coefficients 
        q[Nmax,2*Nmax+1]   q coefficients 
        k                  wavenumber
        mode: VSH_mode     type of VSH (outgoing, incident)
        origin             origin around which to perform the expansion (default: [0,0,0])
    """
    Nmax = p.shape[0]

    def f(x, y, z):
        r, theta, phi = cart_to_sph(x,y,z)
        rhat,that,phat = sph_basis_vectors(theta,phi)

        expanded_E = np.zeros(shape=(3,)+x.shape, dtype=complex)
        for n in range(1,Nmax+1):
            for m in range(-n, n+1):
                Nfunc,Mfunc = VSH(n, m, mode=mode)

                Emn_val = Emn(m, n, E0=1)

                N = Nfunc(r, theta, phi, k)
                M = Mfunc(r, theta, phi, k)

                N = vec_sph_to_cart(N, theta, phi)
                M = vec_sph_to_cart(M, theta, phi)

                expanded_E += -1j*Emn_val*(p[n-1,n+m]*N + q[n-1,n+m]*M)

        return expanded_E
    
    return f