"""
Paraxial reprentation of beams
"""

import numpy as np
import miepy
from miepy.sources.source_base import source, combined_source
from math import factorial
from scipy.special import eval_genlaguerre, eval_hermite, erfc
from scipy.constants import physical_constants
from copy import deepcopy

Z0 = physical_constants['characteristic impedance of vacuum'][0]

def zr(w0, wav):
    return np.pi*w0**2/wav

def w(z, w0, wav):
    return w0*np.sqrt(1 + (z/zr(w0,wav))**2)
    
def Rinv(z, w0, wav):
    return z/(z**2 + zr(w0,wav)**2)

def gouy(z, w0, wav):
    return np.arctan2(z, zr(w0,wav))

def gaussian_paraxial(self, x, y, z, k):
    if self.amplitude is None:
        E0 = 2/self.width*np.sqrt(Z0*self.power/np.pi)
    else:
        E0 = self.amplitude

    rho_sq = x**2 + y**2
    wav = 2*np.pi/k
    amp = E0*self.width/w(z, self.width, wav) * np.exp(-rho_sq/w(z,self.width,wav)**2)
    phase = k*z[2] + k*rho_sq*Rinv(z[2],self.width,wav)/2 - gouy(z[2],self.width,wav)

    return amp*np.exp(1j*phase)
