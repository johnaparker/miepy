"""
Expansion of electric and magnetic fields given expansion coefficients
"""

import numpy as np
from miepy import vsh

#TODO: move k argument to field function for consistency
def expand_E(p, k, mode):
    """Expand VSH coefficients to obtain an electric field function
    Returns E(r,θ,φ) function
    
    Arguments:
        p[2,rmax]          expansion coefficients 
        k                  wavenumber
        mode: vsh_mode     type of VSH (outgoing, incident, interior, ingoing)
    """
    lmax = vsh.rmax_to_lmax(p.shape[1])
    factor = 1j if mode == vsh.vsh_mode.outgoing else -1j

    #TODO: depends on theta.shape
    def f(rad, theta, phi):
        (rad, theta, phi) = map(lambda A: np.asarray(A, dtype=float), (rad, theta, phi))
        E_sph = np.zeros(shape=(3,) + theta.shape, dtype=complex)

        for i,n,m in vsh.mode_indices(lmax):
            Nfunc,Mfunc = vsh.VSH(n, m, mode=mode)

            Emn_val = vsh.Emn(m, n)

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
    lmax = vsh.rmax_to_lmax(p_scat.shape[1])

    #TODO: depends on theta.shape
    def f(rad, theta, phi):
        (rad, theta, phi) = map(lambda A: np.asarray(A, dtype=float), (rad, theta, phi))

        E_sph = np.zeros(shape=(3,) + theta.shape, dtype=complex)
        factor = np.exp(1j*k*rad)/(k*rad)

        for i,n,m in vsh.mode_indices(lmax):
            Emn_val = vsh.Emn(m, n)

            pi = vsh.special.pi_func(n, m, theta)
            tau = vsh.special.tau_func(n, m, theta)

            E_sph[1] += 1j*factor*Emn_val*(-1j)**(n)*(p_scat[0,i]*tau + p_scat[1,i]*pi)*np.exp(1j*m*phi)
            E_sph[2] += -factor*Emn_val*(-1j)**(n)*(p_scat[0,i]*pi + p_scat[1,i]*tau)*np.exp(1j*m*phi)

        return E_sph

    return f

def expand_H(p, k, mode, eps, mu):
    """Expand VSH coefficients to obtain a magnetic field function
    Returns H(r,θ,φ) function
    
    Arguments:
        p[2,rmax]       expansion coefficients 
        k               wavenumber
        mode: vsh_mode     type of VSH (outgoing, incident, interior, ingoing)
        eps     medium permitiviity
        mb      medium permeability
    """

    factor = -1j*np.sqrt(eps/mu)
    E_func = expand_E(p[::-1], k, mode)
    return lambda *args: factor*E_func(*args)

def expand_H_far(p_scat, k, eps, mu):
    """Expand VSH scattering coefficients to obtain a magnetic field function for the far-field
    Returns H(r,θ,φ) function
    
    Arguments:
        p_scat[2,rmax]    scattering coefficients 
        k                 wavenumber
        eps     medium permitiviity
        mu      medium permeability
    """
    factor = -1j*np.sqrt(eps/mu)
    E_func = expand_E_far(p_scat[::-1], k)
    return lambda *args: factor*E_func(*args)

