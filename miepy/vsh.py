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
import miepy.coordinates as coordinates

def rmax_to_Lmax(rmax):
    """obtain Lmax from rmax"""
    Lmax = int(-1 + (1+rmax)**0.5)
    return Lmax

def Lmax_to_rmax(Lmax):
    """obtain rmax from Lmax"""
    rmax = Lmax*(Lmax + 2)
    return rmax

def mode_indices(Lmax):
    """generate (n,m) index pairs up to n=Lmax"""
    n = 1
    m = -1

    for n in range(1,Lmax+1):
        for m in range(-n,n+1):
            yield n, m

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

def wigner_3j(j1, j2, j3, m1, m2, m3):
    """wigner3j coefficients"""
    kmin = max(0, j2-j3-m1, j1-j3+m2)
    kmax = max(j1+j2-j3, j1-m1, j2+m2)

    if m1 + m2 + m3 != 0:
        return 0

    j3min = max(abs(j1-j2), abs(m1+m2))
    if not j3min <= j3 <= j1 + j2:
        return 0

    f = factorial
    binomial = special.binom
    numerator = f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2)*f(j3-m3)*f(j3+m3)
    denominator = f(j1+j2-j3)*f(j1-j2+j3)*f(-j1+j2+j3)*f(j1+j2+j3+1)
    factor = (-1)**(j1+j2+m3)*(numerator/denominator)**0.5

    sum_term = 0
    for k in range(kmin, kmax+1):
        b1 = binomial(j1+j2-j3, k)
        b2 = binomial(j1-j2+j3, j1-m1-k)
        b3 = binomial(-j1+j2+j3, j2+m2-k)
        sum_term += (-1)**k * b1 * b2 * b3
    
    return factor*sum_term

@lru_cache(None)
def a_func(m,n,u,v,p):
    """gaunt coefficient"""

    f = factorial
    numerator = f(n+m)*f(v+u)*f(p-m-u)
    denominator = f(n-m)*f(v-u)*f(p+m+u)
    factor = (-1.)**(m+u)*(2*p+1)*(numerator/denominator)**0.5

    w1 = wigner_3j(n,v,p,0,0,0)
    w2 = wigner_3j(n,v,p,m,u,-m-u)

    return factor*w1*w2

@lru_cache(None)
def b_func(m,n,u,v,p):
    """b function"""

    f = factorial
    numerator = f(n+m)*f(v+u)*f(p-m-u+1)
    denominator = f(n-m)*f(v-u)*f(p+m+u+1)
    factor = (-1.)**(m+u)*(2*p+3)*(numerator/denominator)**0.5

    w1 = wigner_3j(n,v,p,0,0,0)
    w2 = wigner_3j(n,v,p+1,m,u,-m-u)

    return factor*w1*w2

def Emn(m, n):
    return 1j**n*np.sqrt((2*n+1)*factorial(n-m)/(n*(n+1)*factorial(n+m)))

def vsh_translation(m, n, u, v, r, theta, phi, k, mode):
    """VSH translation coefficients"""
    m *= -1
    f = factorial
    zn = get_zn(mode)

    factor = 0.5 * (-1.)**m * np.sqrt((2*v+1)*(2*n+1)*f(v-u)*f(n-m)
            /(v*(v+1)*n*(n+1)*f(v+u)*f(n+m))) * np.exp(1j*(u+m)*phi)

    qmax = min(n, v, (n+v - abs(m+u))//2)
    sum_term = 0
    for q in range(0, qmax+1):
        p = n + v - 2*q
        aq = a_func(m,n,u,v,p)
        A = 1j**p*(n*(n+1) + v*(v+1) - p*(p+1))*aq

        Pnm = associated_legendre(p,u+m)
        sum_term += A*zn(p, k*r)*Pnm(np.cos(theta))

    A_translation = factor*sum_term

    qmax = min(n, v, (n+v+1 - abs(m+u))//2)
    sum_term = 0
    for q in range(1, qmax+1):
        p = n + v - 2*q
        bq = b_func(m,n,u,v,p)
        A = 1j**(p+1)*(((p+1)**2 - (n-v)**2)*((n+v+1)**2 - (p+1)**2))**0.5*bq

        Pnm = associated_legendre(p+1,u+m)
        sum_term += A*zn(p+1, k*r)*Pnm(np.cos(theta))

    B_translation = -factor*sum_term

    return A_translation, B_translation

class VSH_mode(enum.Enum):
    outgoing = enum.auto()
    ingoing  = enum.auto()
    incident = enum.auto()
    interior = enum.auto()

def get_zn(mode):
    """determine the zn function for a given mode"""
    if mode is VSH_mode.outgoing:
        return spherical_hn
    elif mode in (VSH_mode.incident, VSH_mode.ingoing, VSH_mode.interior):
        return special.spherical_jn
    else:
        raise TypeError(f'{mode} is not a valid type of mode')

#TODO: this whole interface could probably be nicer...
#TODO: specify spherical flag (either in VSH or the N/M functions themselves)
def VSH(n, m, mode=VSH_mode.outgoing):
    """electric and magnetic vector spherical harmonic function

            n: int           order
            m: int           degree
            mode: VSH_mode   type of VSH (outgoing, incident)


       returns (N(r,θ,ϕ,k) -> [3,...], M(r,θ,ϕ,k) -> [3,...]), the 3 x,y,z components"""

    pi_f = pi_func(n,m)
    tau_f = tau_func(n,m)
    Pnm = associated_legendre(n,m)

    zn = get_zn(mode)
        
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
    zn = get_zn(mode)

    norm = 1j**n*(2*n+1)*factorial(n-m)/factorial(n+m)
    zn_val = zn(n, k*r)
    angular_term = 4*np.pi*n*(n+1)/np.abs(norm)

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
    THETA, PHI = coordinates.sphere_mesh(sampling)

    tau = np.linspace(-1, 1, sampling)
    phi = np.linspace(0, 2*np.pi, 2*sampling)

    N,M = VSH(n, m, mode)
    if ftype == 'electric':
        base_function = N
    elif ftype == 'magnetic':
        base_function = M

    vsh_data = base_function(r,THETA,PHI,k).squeeze()

    if not spherical:
        E = coordinates.vec_cart_to_sph(E, THETA, PHI)

    Emn_val = Emn(m, n)

    if mode == VSH_mode.outgoing:
        factor = 1/(1j*Emn_val)
    elif mode == VSH_mode.incident or mode == VSH_mode.ingoing:
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

    THETA, PHI = coordinates.sphere_mesh(sampling)
    X,Y,Z = coordinates.sph_to_cart(r, THETA, PHI, origin=origin)
    E = src.E_field(X, Y, Z, k)

    return project_fields_onto(E, r, k, ftype, n, m, mode, spherical=False)

def decompose_fields(E, r, k, Lmax, mode=VSH_mode.outgoing, spherical=False):
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

def decompose_source(src, k, Lmax, origin=[0,0,0], sampling=30, mode=VSH_mode.incident):
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

def expand_E(p, k, mode):
    """Expand VSH coefficients to obtain an electric field function
    Returns E(r,θ,φ) function
    
    Arguments:
        p[2,rmax]          expansion coefficients 
        k                  wavenumber
        mode: VSH_mode     type of VSH (outgoing, incident, interior, ingoing)
    """
    Lmax = rmax_to_Lmax(p.shape[1])
    factor = 1j if mode == VSH_mode.outgoing else -1j

    #TODO: depends on theta.shape
    def f(rad, theta, phi):
        (rad, theta, phi) = map(lambda A: np.asarray(A, dtype=float), (rad, theta, phi))
        E_sph = np.zeros(shape=(3,) + theta.shape, dtype=complex)

        for i,(n,m) in enumerate(mode_indices(Lmax)):
            Nfunc,Mfunc = VSH(n, m, mode=mode)

            Emn_val = Emn(m, n)

            N = Nfunc(rad, theta, phi, k)
            M = Mfunc(rad, theta, phi, k)

            E_sph += factor*Emn_val*(p[0,i]*N + p[1,i]*M)

        return E_sph
    
    return f

def expand_E_far(p_scat, k):
    """Expand VSH scattering coefficients to obtain an electric field function for the far-field
    Returns E(r,θ,φ) function
    
    Arguments:
        p_scat[2,rmax]    scattering coefficients 
        k                 wavenumber
    """
    Lmax = rmax_to_Lmax(p_scat.shape[1])

    #TODO: depends on theta.shape
    def f(rad, theta, phi):
        (rad, theta, phi) = map(lambda A: np.asarray(A, dtype=float), (rad, theta, phi))

        E_sph = np.zeros(shape=(3,) + theta.shape, dtype=complex)
        factor = np.exp(1j*k*rad)/(k*rad)

        for i,(n,m) in enumerate(mode_indices(Lmax)):
            Emn_val = Emn(m, n)

            tau = tau_func(n,m)(theta)
            pi = pi_func(n,m)(theta)

            E_sph[1] += 1j*factor*Emn_val*(-1j)**(n)*(p_scat[0,i]*tau + p_scat[1,i]*pi)*np.exp(1j*m*phi)
            E_sph[2] += -factor*Emn_val*(-1j)**(n)*(p_scat[0,i]*pi + p_scat[1,i]*tau)*np.exp(1j*m*phi)

        return E_sph

    return f

def expand_H(p, k, mode, eps_b, mu_b):
    """Expand VSH coefficients to obtain a magnetic field function
    Returns H(r,θ,φ) function
    
    Arguments:
        p[2,rmax]       expansion coefficients 
        k               wavenumber
        mode: VSH_mode     type of VSH (outgoing, incident, interior, ingoing)
        eps_b     background permitiviity
        mu_b      background permeability
    """

    factor = -1j*np.sqrt(eps_b/mu_b)
    E_func = expand_E(p[::-1], k, mode)
    return lambda *args: factor*E_func(*args)

def expand_H_far(p_scat, k, eps_b, mu_b):
    """Expand VSH scattering coefficients to obtain a magnetic field function for the far-field
    Returns H(r,θ,φ) function
    
    Arguments:
        p_scat[2,rmax]    scattering coefficients 
        k                 wavenumber
        eps_b     background permitiviity
        mu_b      background permeability
    """
    factor = -1j*np.sqrt(eps_b/mu_b)
    E_func = expand_E_far(p[::-1], k)
    return lambda *args: factor*E_func(*args)

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
    Lmax_in = rmax_to_Lmax(p_scat.shape[-1])

    if Lmax is None:
        Lmax = Lmax_in

    rmax = Lmax_to_rmax(Lmax)
    p_cluster = np.zeros([2,rmax], dtype=complex)

    for i in range(Nparticles):
        if np.all(positions[i] == origin):
            p_cluster[...] = p_scat[i]
            continue

        rij = origin - positions[i]
        rad, theta, phi = coordinates.cart_to_sph(*rij)
        
        for r,(n,m) in enumerate(mode_indices(Lmax)):
            for rp,(v,u) in enumerate(mode_indices(Lmax_in)):
                a = p_scat[i,0,rp]
                b = p_scat[i,1,rp]

                A, B = vsh_translation(m, n, u, v, rad, theta, phi, k, VSH_mode.incident)

                p_cluster[0,r] += a*A + b*B
                p_cluster[1,r] += a*B + b*A

    return p_cluster
