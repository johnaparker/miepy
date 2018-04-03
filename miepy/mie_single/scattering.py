"""
Scattering defines all functions that make use of the scattering coefficients an, bn
Calculations include scattering, absorbption, and electric and magnetic field computations
Mie sphere and Mie core shell both contain an, bn as part of their solution
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
import scipy.constants as constants
from miepy.special_functions import riccati_1,riccati_2,vector_spherical_harmonics

def scattering_per_multipole(an, bn, k):
    """Scattering cross-section per multipole. Returns scat[Nfreq,2,Lmax].
            an[N]    an scattering coefficients
            bn[N]    bn scattering coefficients
            k[N]     wavenumbers
    """
    Nfreq, Lmax = an.shape
    flux = np.zeros([Nfreq,2,Lmax])
    nvals = np.arange(1, Lmax+1)
    flux[:,0,:] = 2*np.pi*(2*nvals+1)*np.abs(an)**2/k[:,np.newaxis]**2
    flux[:,1,:] = 2*np.pi*(2*nvals+1)*np.abs(bn)**2/k[:,np.newaxis]**2

    return flux

def extinction_per_multipole(an, bn, k):
    """Extinction cross-section per multipole. Returns extinct[Nfreq,2,Lmax].
            an[N]    an scattering coefficients
            bn[N]    bn scattering coefficients
            k[N]     wavenumbers
    """
    Nfreq, Lmax = an.shape
    flux = np.zeros([Nfreq,2,Lmax])
    nvals = np.arange(1, Lmax+1)
    flux[:,0,:] = 2*np.pi*(2*nvals+1)*np.real(an)/k[:,np.newaxis]**2
    flux[:,1,:] = 2*np.pi*(2*nvals+1)*np.real(bn)/k[:,np.newaxis]**2

    return flux

def absorbption_per_multipole(an, bn, k):
    """Absorbption cross-section per multipole. Returns absorb[Nfreq,2,Lmax].
            an[N]    an scattering coefficients
            bn[N]    bn scattering coefficients
            k[N]     wavenumbers
    """
    return extinction_per_multipole(an, bn, k) - scattering_per_multipole(an, bn, k)

def cross_sections(an, bn, k):
    """Return the 3 cross-sections, (Scattering, Absorbption, Extinction)
            an[N]    an scattering coefficients
            bn[N]    bn scattering coefficients
            k[N]     wavenumbers
    """
    scat_flux = scattering_per_multipole(an, bn, k)
    extinct_flux = extinction_per_multipole(an, bn, k)
    abs_flux = extinct_flux - scat_flux

    return map(lambda arr: np.sum(arr, axis=(1,2)), [scat_flux, abs_flux, extinct_flux])

def multipole_label(T,L):
    """Get multipole label.
            T = 0 (electric), 1(magnetic)
            L = 0,1,2... (order)
    """
    first = ['e', 'm'][T]
    if L <= 3:
        last = ['D', 'Q', 'O', 'H'][L]
    else:
        last = f" (L = {L})" 
    return first + last

def scattered_E(an, bn, k):
    """For a given an, bn, k, return the scattered electric field function E(r,theta,phi)
                an[L]       an coefficients
                an[L]       bn coefficients
                k           wavenumber in the medium
    """
    Lmax = an.shape[0]
    def E_func(r, theta, phi):
        E = np.zeros(shape = [3] + list(r.shape), dtype=np.complex)
        for L in range(1,Lmax+1):
            En = 1j**L*(2*L+1)/(L*(L+1))

            VSH = vector_spherical_harmonics(L,3)
            E += En*(1j*an[L-1]*VSH.N_e1n(k)(r,theta,phi)  \
                        - bn[L-1]*VSH.M_o1n(k)(r,theta,phi))
        return -E
    return E_func

def interior_E(cn, dn, k):
    """For a given cn, dn, k, return the interior electric field function E(r,theta,phi) for a sphere
                cn[L]       cn coefficients
                dn[L]       dn coefficients
                k           wavenumber inside the sphere
    """
    Lmax = cn.shape[0]
    def E_func(r, theta, phi):
        E = np.zeros(shape = [3] + list(r.shape), dtype=np.complex)
        for L in range(1,Lmax+1):
            En = 1j**L*(2*L+1)/(L*(L+1))

            VSH = vector_spherical_harmonics(L,1)
            E += En*(cn[L-1]*VSH.M_o1n(k)(r,theta,phi)  \
                     - 1j*dn[L-1]*VSH.N_e1n(k)(r,theta,phi))
        return -E
    return E_func

def scattered_H(an, bn, k, n_b, mu_b):
    """For a given an, bn, k, return the scattered electric field function H(r,theta,phi)
                an[L]       an coefficients
                an[L]       bn coefficients
                k           wavenumber in the medium
                n_b         index of refraction of the medium
                mu_b        permeability of the medium
    """
    Lmax = an.shape[0]
    def H_func(r, theta, phi):
        H = np.zeros(shape = [3] + list(r.shape), dtype=np.complex)
        for L in range(1,Lmax+1):
            En = 1j**L*(2*L+1)/(L*(L+1))

            VSH = vector_spherical_harmonics(L,3)
            H += n_b*En/mu_b*(1j*bn[L-1]*VSH.N_o1n(k)(r,theta,phi)  \
                            + an[L-1]*VSH.M_e1n(k)(r,theta,phi))
        return -H
    return H_func

def interior_H(cn, dn, k, n, mu):
    """For a given cn, dn, k, return the interior electric field function H(r,theta,phi) for a sphere
                cn[L]       cn coefficients
                dn[L]       dn coefficients
                k           wavenumber inside the sphere
                n           index of refraction of the sphere
                mu          permeability of the sphere
    """
    Lmax = cn.shape[0]
    def H_func(r, theta, phi):
        H = np.zeros(shape = [3] + list(r.shape), dtype=np.complex)
        for L in range(1,Lmax+1):
            En = 1j**L*(2*L+1)/(L*(L+1))

            VSH = vector_spherical_harmonics(L,1)
            H += -n*En/mu*(dn[L-1]*VSH.M_e1n(k)(r,theta,phi)  \
                            + 1j*cn[L-1]*VSH.N_o1n(k)(r,theta,phi))
        return -H
    return H_func
