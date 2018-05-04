"""
plane wave sources
"""

import numpy as np
import miepy
from miepy.vsh import pi_func, tau_func
from miepy.sources.source_base import source

class plane_wave(source):
    def __init__(self, polarization, amplitude=1):
        super().__init__(amplitude)
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
    
    def E_field(self, x, y, z, k):
        amp = self.amplitude*np.exp(1j*k*z)
        pol = np.array([*self.polarization, 0])
        return np.einsum('i...,...->i...', pol, amp)

    def H_field(self, x, y, z, k):
        amp = self.amplitude*np.exp(1j*k*z)
        H0_x, H0_y = -self.polarization[1], self.polarization[0]
        pol = np.array([H0_x, H0_y, 0])
        return np.einsum('i...,...->i...', pol, amp)

    def structure(self, position, k, Lmax):
        rmax = miepy.vsh.Lmax_to_rmax(Lmax)
        p_src = np.zeros([2, rmax], dtype=complex)

        phase = 1j*k*position[2]
        alpha = 0

        for i, (n,m) in enumerate(miepy.vsh.mode_indices(Lmax)):
            pi_value = pi_func(n, m)(alpha)
            tau_value = tau_func(n, m)(alpha)

            p_src[0,i] = self.amplitude*np.exp(phase)*np.sqrt(2*n+1)*tau_value*(self.polarization[0] - 1j*m*self.polarization[1])
            p_src[1,i] = self.amplitude*np.exp(phase)*np.sqrt(2*n+1)*pi_value*(self.polarization[0] - 1j*m*self.polarization[1])

            if m == 1:
                p_src[:,i] /= (n*(n+1))

        return p_src

def x_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[1,0], amplitude=amplitude)

def y_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[0,1], amplitude=amplitude)

def rhc_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[1,1j], amplitude=amplitude)

def lhc_polarized_plane_wave(amplitude=1):
    return plane_wave(polarization=[1,-1j], amplitude=amplitude)
