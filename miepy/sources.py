"""
Pre-defined sources that can be used with particle_system
"""

import numpy as np

from abc import ABCMeta, abstractmethod
from miepy.vsh import pi_func, tau_func, project_source_onto
from math import factorial

class source:
    """source interface base class"""
    __metaclass__ = ABCMeta

    def __init__(self, amplitude):
        self.amplitude = amplitude

    @abstractmethod
    def E(self, r, k): pass

    @abstractmethod
    def H(self, r, k): pass

    @abstractmethod
    def structure_of_mode(self, n, m, r, k): pass

    def structure(self, position, k, Nmax):
        p = np.zeros([Nmax, 2*Nmax+1], dtype=complex)
        q = np.zeros([Nmax, 2*Nmax+1], dtype=complex)
        for n in range(1, Nmax+1):
            for m in range(-n,n+1):
                p[n-1,m+n], q[n-1,m+n] = self.structure_of_mode(n, m, position, k)

        return p,q

class plane_wave(source):
    def __init__(self, polarization, amplitude=1):
        super().__init__(amplitude)
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
    
    def E(self, r, k):
        amp = self.amplitude*np.exp(1j*k*r[2])
        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp)

    def H(self, r, k):
        amp = self.amplitude*np.exp(1j*k*r[2])
        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])
        return np.einsum('i...,...->i...', pol, amp)

    def structure_of_mode(self, n, m, r, k):
        phase = np.exp(1j*k*r[2])

        alpha = 0
        pi_value = pi_func(n,m)(alpha)
        tau_value = tau_func(n,m)(alpha)

        p = phase/(n*(n+1))*(tau_value*self.polarization[0] - 1j*pi_value*self.polarization[1])
        q = m*phase/(n*(n+1))*(tau_value*self.polarization[0] - 1j*pi_value*self.polarization[1])
        return (p,q)

        # if m == 1:
        #     return (phase*1/2, phase*1/2)
        # elif m == -1:
        #     return (-phase*1/(2*n*(n+1)), phase*1/(2*n*(n+1)))
        # else:
        #     return (np.zeros_like(r[0]), np.zeros_like(r[0]))

def x_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[1,0], amplitude=amplitude)

def y_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[0,1], amplitude=amplitude)

def rhc_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[1,1j], amplitude=amplitude)

def lhc_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[1,-1j], amplitude=amplitude)



def zr(w0, wav):
    return np.pi*w0**2/wav

def w(z, w0, wav):
    return w0*np.sqrt(1 + (z/zr(w0,wav))**2)
    
def Rinv(z, w0, wav):
    return z/(z**2 + zr(w0,wav)**2)

def gouy(z, w0, wav):
    return np.arctan2(z, zr(w0,wav))

class gaussian_beam(source):
    def __init__(self, width, polarization, amplitude=1):
        super().__init__(amplitude)
        self.width = width
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
    
    def E(self, r, k):
        rho_sq = r[0]**2 + r[1]**2
        wav = 2*np.pi/k
        amp = self.amplitude*self.width/w(r[2], self.width, wav) * np.exp(-rho_sq/w(r[2],self.width,wav)**2)
        phase = k*r[2] + k*rho_sq*Rinv(r[2],self.width,wav)/2 - gouy(r[2],self.width,wav)
        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))


    def H(self, r, k):
        rho_sq = r[0]**2 + r[1]**2
        wav = 2*np.pi/k
        amp = self.amplitude*self.width/w(r[2], self.width, wav) * np.exp(-rho_sq/w(r[2],self.width,wav)**2)
        phase = k*r[2] + k*rho_sq*Rinv(r[2],self.width,wav)/2 - gouy(r[2],self.width,wav)
        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])

        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))

    def structure_of_mode(self, n, m, r, k):
        p = project_source_onto(self, k, 'electric', n, m, r)
        q = project_source_onto(self, k, 'magnetic', n, m, r)

        return (p,q)



class laguerre_gaussian_beam(source):
    def __init__(self, l, width, polarization, amplitude=1):
        super().__init__(amplitude)
        self.l = l
        self.p = 1
        self.width = width
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
    
    def E(self, r, k):
        rho_sq = r[0]**2 + r[1]**2
        phi = np.arctan2(r[1], r[0])
        wav = 2*np.pi/k
        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        N = abs(self.l) + 2*self.p

        amp = self.amplitude*C/w(r[2], self.width, wav) * np.exp(-rho_sq/w(r[2],self.width,wav)**2) * ((2*rho_sq)**0.5/(w(r[2], self.width, wav)))**abs(self.l) * (1 + self.l - 2*rho_sq/w(r[2],self.width,wav)**2)
        phase = self.l*phi + k*r[2] + k*rho_sq*Rinv(r[2],self.width,wav)/2 - (N+1)*gouy(r[2],self.width,wav)

        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))


    def H(self, r, k):
        rho_sq = r[0]**2 + r[1]**2
        phi = np.arctan2(r[1], r[0])
        wav = 2*np.pi/k
        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        N = abs(self.l) + 2*self.p

        amp = self.amplitude*C/w(r[2], self.width, wav) * np.exp(-rho_sq/w(r[2],self.width,wav)**2) * ((2*rho_sq)**0.5/(w(r[2], self.width, wav)))**abs(self.l) * (1 + self.l - 2*rho_sq/w(r[2],self.width,wav)**2)
        phase = self.l*phi + k*r[2] + k*rho_sq*Rinv(r[2],self.width,wav)/2 - (N+1)*gouy(r[2],self.width,wav)

        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])

        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))

    def structure_of_mode(self, n, m, r, k):
        p = project_source_onto(self, k, 'electric', n, m, r)
        q = project_source_onto(self, k, 'magnetic', n, m, r)

        return (p,q)
