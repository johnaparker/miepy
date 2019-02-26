"""
Defines the vsh wave functions and related functions
"""

import numpy as np
from scipy import special
import enum
from math import factorial
from miepy import vsh
from miepy.cpp.vsh_functions import vsh_mode, Emn

def get_zn(mode):
    """determine the zn function for a given mode"""
    if mode is vsh_mode.outgoing:
        return vsh.special.spherical_hn
    elif mode is vsh_mode.ingoing:
        return vsh.special.spherical_hn_2
    elif mode in (vsh_mode.incident, vsh_mode.interior):
        return vsh.special.spherical_jn
    else:
        raise TypeError('{mode} is not a valid type of mode'.format(mode=mode))

def get_zn_far(mode):
    """determine the zn function for a given mode, in the far-field limit"""
    if mode is vsh_mode.outgoing:
        return lambda n, z: np.exp(1j*(z - (n+1)*np.pi/2))/z
    elif mode is vsh_mode.ingoing:
        return lambda n, z: np.exp(-1j*(z - (n+1)*np.pi/2))/z
    elif mode in (vsh_mode.incident, vsh_mode.interior):
        return lambda n, z: np.cos(z - (n+1)*np.pi/2)/z
    else:
        raise TypeError('{mode} is not a valid type of mode'.format(mode=mode))

#TODO: this whole interface could probably be nicer...
#TODO: specify spherical flag (either in VSH or the N/M functions themselves)
#TODO: expansion issues at origin (r=0) for incident modes
def VSH(n, m, mode=vsh_mode.outgoing):
    """electric and magnetic vector spherical harmonic function

            n: int           order
            m: int           degree
            mode: vsh_mode   type of VSH (outgoing, incident)

       returns (N(r,θ,ϕ,k) -> [3,...], M(r,θ,ϕ,k) -> [3,...]), the 3 spherical components"""

    pi_f = vsh.special.pi_func
    tau_f = vsh.special.tau_func
    Pnm = vsh.special.associated_legendre

    zn = get_zn(mode)
        
    def N(r, theta, phi, k):
        H = zn(n, k*r)
        Hp = zn(n, k*r, derivative=True)
        Pnm_val = Pnm(n, m, np.cos(theta))

        factor = (H + r*k*Hp)*np.exp(1j*m*phi)/(k*r)

        r_comp = n*(n+1)*Pnm_val*H/(k*r)*np.exp(1j*m*phi)
        theta_comp = tau_f(n, m, theta)*factor
        phi_comp = 1j*pi_f(n, m, theta)*factor

        return np.array([r_comp, theta_comp, phi_comp])

    def M(r, theta, phi, k):
        H = zn(n, k*r)
        factor = H*np.exp(1j*m*phi)

        theta_comp = 1j*pi_f(n, m, theta)*factor
        phi_comp = -1*tau_f(n, m, theta)*factor
        r_comp = np.zeros_like(theta_comp)

        return np.array([r_comp, theta_comp, phi_comp])

    return N,M

def VSH_far(n, m, mode=vsh_mode.outgoing):
    """electric and magnetic vector spherical harmonic function in the far field

            n: int           order
            m: int           degree
            mode: vsh_mode   type of VSH (outgoing, incident)

       returns (N(r,θ,ϕ,k) -> [2,...], M(r,θ,ϕ,k) -> [2,...]), the 2 theta/phi components"""

    pi_f = vsh.special.pi_func
    tau_f = vsh.special.tau_func
    zn = get_zn_far(mode)
    sign = -1 if mode is vsh.vsh_mode.ingoing else 1
        
    def N(r, theta, phi, k):
        factor = sign*zn(n, k*r)*np.exp(1j*m*phi)
        theta_comp = 1j*tau_f(n, m, theta)*factor
        phi_comp = -pi_f(n, m, theta)*factor

        return np.array([theta_comp, phi_comp])

    def M(r, theta, phi, k):
        factor = zn(n, k*r)*np.exp(1j*m*phi)
        theta_comp = 1j*pi_f(n, m, theta)*factor
        phi_comp = -tau_f(n, m, theta)*factor

        return np.array([theta_comp, phi_comp])

    return N,M

def vsh_normalization_values(mode, ftype, n, m, r, k):
    """Determine the norm of a given vsh mode
    
    Arguments:
        mode: vsh_mode    type of VSH (outgoing, incident)
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

def vsh_normalization_values_far(n, m):
    """Determine the norm of a given vsh mode in the far-field
    
    Arguments:
        n                 vsh order (1, 2, ...)
        m                 vsh orientation (-n, -n+1, ..., n)
    """

    return 4*np.pi*n*(n+1)*factorial(n+m)/((2*n+1)*factorial(n-m))
