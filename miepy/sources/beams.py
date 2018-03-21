"""
beam sources
"""

import numpy as np
from miepy.vsh import project_source_onto
from miepy.sources.source_base import source, combined_source
from math import factorial
from scipy.special import eval_genlaguerre, eval_hermite

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
    
    def E(self, r, k):
        rp = np.array([r[0] - self.center[0], r[1] - self.center[1], r[2] - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k
        amp = self.amplitude*self.width/w(rp[2], self.width, wav) * np.exp(-rho_sq/w(rp[2],self.width,wav)**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - gouy(rp[2],self.width,wav)
        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))


    def H(self, r, k):
        rp = np.array([r[0] - self.center[0], r[1] - self.center[1], r[2] - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k
        amp = self.amplitude*self.width/w(rp[2], self.width, wav) * np.exp(-rho_sq/w(rp[2],self.width,wav)**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - gouy(rp[2],self.width,wav)
        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])

        return np.einsum('i...,...->i...', pol, amp*np.exp(1j*phase))

    def structure_of_mode(self, n, m, r, k):
        p = project_source_onto(self, k, 'electric', n, m, r)
        q = project_source_onto(self, k, 'magnetic', n, m, r)

        return (p,q)

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
    
    def E(self, r, k):
        rp = np.array([r[0] - self.center[0], r[1] - self.center[1], r[2] - self.center[2]])
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


    def H(self, r, k):
        rp = np.array([r[0] - self.center[0], r[1] - self.center[1], r[2] - self.center[2]])
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

    def structure_of_mode(self, n, m, r, k):
        p = project_source_onto(self, k, 'electric', n, m, r)
        q = project_source_onto(self, k, 'magnetic', n, m, r)

        return (p,q)


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
    
    def E(self, r, k):
        rp = np.array([r[0] - self.center[0], r[1] - self.center[1], r[2] - self.center[2]])
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


    def H(self, r, k):
        rp = np.array([r[0] - self.center[0], r[1] - self.center[1], r[2] - self.center[2]])
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

    def structure_of_mode(self, n, m, r, k):
        p = project_source_onto(self, k, 'electric', n, m, r)
        q = project_source_onto(self, k, 'magnetic', n, m, r)

        return (p,q)

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
