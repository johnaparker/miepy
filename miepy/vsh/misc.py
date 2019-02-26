"""
Miscellaneous functions related to vsh
"""

import numpy as np
from scipy import special
from scipy.integrate import simps
from miepy.cpp.decomposition import trapz, trapz_2d
import miepy
from miepy.constants import Z0

def simps_2d(xd,yd,fd):
    """1d simpsons rule extended to 2d"""
    if np.iscomplexobj(fd):
        return simps_2d(xd, yd, fd.real) + 1j*simps_2d(xd, yd, fd.imag)

    xData = np.zeros(len(xd))
    for i,x in enumerate(xd):
        xData[i] = simps(fd[i,:], yd)

    return simps(xData, xd)

def power_through_aperature(source, center, radius, k, sampling=150):
    """Calculate the power through an aperature for a source

       Arguments:
           source      source object
           center      center position of the aperature
           radius      radius of the aperature
           k           wavenumber
           sampling    numer of samples used along the diamter of the aperature (default: 150)
    """
    x = center[0] + np.linspace(-radius, radius, sampling)
    y = center[1] + np.linspace(-radius, radius, sampling)
    X, Y = np.meshgrid(x, y)
    Z = center[2] + np.zeros_like(X)

    lmax = int(np.ceil(k*radius))
    lmax = 8
    p_src = source.structure(center, k, lmax)
    mask = (X**2 + Y**2 < radius **2)

    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)
    Efunc = miepy.expand_E(p_src, k, miepy.vsh_mode.incident)
    Hfunc = miepy.expand_H(p_src, k, miepy.vsh_mode.incident, 1, 1)

    E = np.zeros((3,) + X.shape, dtype=complex)
    H = np.zeros((3,) + X.shape, dtype=complex)
    E[:,mask] = Efunc(R[mask], THETA[mask], PHI[mask])
    H[:,mask] = Hfunc(R[mask], THETA[mask], PHI[mask])
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
    H = miepy.coordinates.vec_sph_to_cart(H, THETA, PHI)/Z0

    S = 0.5*np.cross(E, np.conjugate(H), axis=0)[2]

    P = trapz_2d(x, y, S).real
    return P

###### below are pi,tau,VSH used in Mie theory, which may differ from those defined in GMT ######

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
