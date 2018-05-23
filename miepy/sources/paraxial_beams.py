"""
beam sources
"""

import numpy as np
import miepy
from miepy.sources.source_base import source, combined_source
from math import factorial
from scipy.special import eval_genlaguerre, eval_hermite

def fields_from_potential(U, polarization):
    """Determine cartesian fields from the scalar potential.
    Returns E[2,...], the x & y components
    
    Arguments:
        U[...]            values of the scalar potential
        polarization[2]   polarization direction
    """
    Ex = U*polarization[0]
    Ey = U*polarization[1]

    return np.array([Ex, Ey])

def far_fields_from_potential(U, polarization, phi):
    """Determine far-field spherical fields from the scalar potential.
    Returns E[2,...], the theta & phi components
    
    Arguments:
        U[...]            values of the scalar potential
        polarization[2]   polarization direction
        phi[...]          phi angles
    """
    Ex = U*polarization[0]
    Ey = U*polarization[1]
    Etheta = -Ex*np.cos(phi) - Ey*np.sin(phi)
    Ephi   = -Ex*np.sin(phi) + Ey*np.cos(phi)

    return np.array([Etheta, Ephi])

def sampling_from_Lmax(Lmax):
    """Determine the required sampling from Lmax for point matching"""
    rmax = miepy.vsh.Lmax_to_rmax(Lmax)
    np.ceil(rmax**0.5)
    return max(3, rmax)

def zr(w0, wav):
    return np.pi*w0**2/wav

def w(z, w0, wav):
    return w0*np.sqrt(1 + (z/zr(w0,wav))**2)
    
def Rinv(z, w0, wav):
    return z/(z**2 + zr(w0,wav)**2)

def gouy(z, w0, wav):
    return np.arctan2(z, zr(w0,wav))

class gaussian_beam(source):
    def __init__(self, width, polarization, amplitude=1, center=np.zeros(3)):
        super().__init__(amplitude)
        self.width = width
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
        self.center = np.asarray(center)
    
    def E_field(self, x, y, z, k):
        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k
        amp = self.amplitude*self.width/w(rp[2], self.width, wav) * np.exp(-rho_sq/w(rp[2],self.width,wav)**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - gouy(rp[2],self.width,wav)
        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))


    def H_field(self, x, y, z, k):
        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k
        amp = self.amplitude*self.width/w(rp[2], self.width, wav) * np.exp(-rho_sq/w(rp[2],self.width,wav)**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - gouy(rp[2],self.width,wav)
        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])

        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))

    def is_paraxial(self, k):
        return 2*np.pi/k > self.width

    def structure(self, position, k, Lmax, radius):
        sampling = sampling_from_Lmax(Lmax)

        if self.is_paraxial(k):
            return miepy.vsh.decomposition.near_field_point_matching(self, 
                              position, 2*radius, k, Lmax, sampling)
        else:
            r = 1e6*(2*np.pi/k)
            return miepy.vsh.decomposition.far_field_point_matching(self, 
                              position, r, k, Lmax, sampling)

    def spherical_ingoing(self, theta, phi, k):
        U = np.exp(-(k*self.width*np.tan(theta)/2)**2)
        return far_fields_from_potential(U, self.polarization, phi)

class hermite_gaussian_beam(source):
    def __init__(self, l, m, width, polarization, amplitude=1, center=np.zeros(3)):
        super().__init__(amplitude)
        self.l = l
        self.m = m
        self.width = width
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
        self.center = np.asarray(center)
    
    def E_field(self, x, y, z, k):
        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k

        wz = w(rp[2], self.width, wav)
        HG_l = eval_hermite(self.l, np.sqrt(2)*rp[0]/wz)
        HG_m = eval_hermite(self.m, np.sqrt(2)*rp[1]/wz)
        N = self.l + self.m

        amp = self.amplitude*self.width/wz * HG_l * HG_m * np.exp(-rho_sq/wz**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - (N+1)*gouy(rp[2],self.width,wav)

        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))


    def H_field(self, x, y, z, k):
        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k

        wz = w(rp[2], self.width, wav)
        HG_l = eval_hermite(self.l, np.sqrt(2)*rp[0]/wz)
        HG_m = eval_hermite(self.m, np.sqrt(2)*rp[1]/wz)
        N = self.l + self.m

        amp = self.amplitude*self.width/wz * HG_l * HG_m * np.exp(-rho_sq/wz**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - (N+1)*gouy(rp[2],self.width,wav)

        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])

        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))

class laguerre_gaussian_beam(source):
    def __init__(self, p, l, width, polarization, amplitude=1, center=np.zeros(3)):
        super().__init__(amplitude)
        self.p = p
        self.l = l
        self.width = width
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
        self.center = np.asarray(center)
    
    def E_field(self, x, y, z, k):
        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        phi = np.arctan2(rp[1], rp[0])
        wav = 2*np.pi/k

        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        wz = w(rp[2], self.width, wav)

        Lpl = eval_genlaguerre(self.p, abs(self.l), 2*rho_sq/wz**2)
        N = abs(self.l) + 2*self.p

        amp = self.amplitude*C/wz * np.exp(-rho_sq/wz**2) * ((2*rho_sq)**0.5/wz)**abs(self.l) * Lpl
        phase = self.l*phi + k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - (N+1)*gouy(rp[2],self.width,wav)

        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))


    def H_field(self, x, y, z, k):
        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        phi = np.arctan2(rp[1], rp[0])
        wav = 2*np.pi/k

        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        wz = w(rp[2], self.width, wav)

        Lpl = eval_genlaguerre(self.p, abs(self.l), 2*rho_sq/wz**2)
        N = abs(self.l) + 2*self.p

        amp = self.amplitude*C/wz * np.exp(-rho_sq/wz**2) * ((2*rho_sq)**0.5/wz)**abs(self.l) * Lpl
        phase = self.l*phi + k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - (N+1)*gouy(rp[2],self.width,wav)

        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])

        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))

def azimuthal_beam(width, amplitude=1, center=np.zeros(3)):
    """azimuthally polarized beam"""
    HG_1 = hermite_gaussian_beam(1, 0, width, [0,1], amplitude, center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [-1,0], amplitude, center)
    return HG_1 + HG_2

def radial_beam(width, amplitude=1, center=np.zeros(3)):
    """radially polarized beam"""
    HG_1 = hermite_gaussian_beam(1, 0, width, [1,0], amplitude, center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [0,1], amplitude, center)
    return HG_1 + HG_2

def shear_beam(width, amplitude=1, center=np.zeros(3)):
    """shear polarized beam"""
    HG_1 = hermite_gaussian_beam(1, 0, width, [0,1], amplitude, center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [1,0], amplitude, center)
    return HG_1 + HG_2
