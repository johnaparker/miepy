"""
plane wave sources
"""

import numpy as np
from miepy.vsh import pi_func, tau_func
from miepy.sources.source_base import source

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
