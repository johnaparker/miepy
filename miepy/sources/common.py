"""
commonly defined sources
"""

import numpy as np
import miepy
from miepy.sources import polarized_beam
from math import factorial
from scipy.special import eval_genlaguerre, eval_hermite, erfc
from copy import deepcopy
from miepy.constants import Z0

class gaussian_beam(polarized_beam):
    def __init__(self, width, polarization, power=1, theta_max=np.pi/2, phase=0, center=None,
                theta=0, phi=0, standing=False):
        super().__init__(polarization=polarization, power=power, theta_max=theta_max,
                phase=phase, center=center, theta=theta, phi=phi, standing=standing)
        self.width = width

    def __repr__(self):
        return f'gaussian_beam(width={self.width}, polarization={self.polarization}, power={self.power}, ' \
               f'center={self.center}, theta={self.theta}, phi={self.phi})'

    def E0(self, k):
        c = 0.5*(k*self.width)**2
        if c < 700:
            U0 = np.sqrt(Z0*self.power/(np.pi*(1 - np.sqrt(np.pi*c)*np.exp(c)*erfc(np.sqrt(c)))))
        else:
            U0 = np.sqrt(Z0*self.power/(np.pi*(1/(2*c) - 3/(4*c**2) + 15/(8*c**3))))

        return U0

    def scalar_angular_spectrum(self, theta, phi, k):
        return np.exp(-(k*self.width*np.tan(theta)/2)**2)

    def theta_cutoff(self, k, cutoff=1e-9, tol=None):
        return np.arctan(np.sqrt(-2*np.log(cutoff))/(k*self.width))

class bigaussian_beam(polarized_beam):
    def __init__(self, width_x, width_y, polarization, power=1, theta_max=np.pi/2, phase=0, center=None,
                theta=0, phi=0, standing=False):
        super().__init__(polarization=polarization, power=power, theta_max=theta_max,
                phase=phase, center=center, theta=theta, phi=phi, standing=standing)
        self.width_x = width_x
        self.width_y = width_y

    def __repr__(self):
        return f'bigaussian_beam(width_x={self.width_x}, width_y={self.width_y}, polarization={self.polarization}, ' \
               f'power={self.power}, center={self.center}, theta={self.theta}, phi={self.phi})'

    def scalar_angular_spectrum(self, theta, phi, k):
        return np.exp(-(k*np.tan(theta)/2)**2*((self.width_x*np.cos(phi))**2 + (self.width_y*np.sin(phi))**2))

class hermite_gaussian_beam(polarized_beam):
    def __init__(self, l, m, width, polarization, power=1, theta_max=np.pi/2, phase=0, center=None,
                theta=0, phi=0, standing=False):
        super().__init__(polarization=polarization, power=power, theta_max=theta_max,
                phase=phase, center=center, theta=theta, phi=phi, standing=standing)
        self.width = width
        self.l = l
        self.m = m

    def __repr__(self):
        return f'HG_beam(width={self.width}, l={self.l}, m={self.m}, polarization={self.polarization}, ' \
               f'power={self.power}, center={self.center}, theta={self.theta}, phi={self.phi})'

    def scalar_angular_spectrum(self, theta, phi, k):
        HG_l = eval_hermite(self.l, k*self.width/np.sqrt(2)*np.tan(theta)*np.cos(phi))
        HG_m = eval_hermite(self.m, k*self.width/np.sqrt(2)*np.tan(theta)*np.sin(phi))
        exp = np.exp(-(k*self.width*np.tan(theta)/2)**2)
        factor = (-1j)**(self.l + self.m)

        return factor * HG_l * HG_m * exp

class laguerre_gaussian_beam(polarized_beam):
    def __init__(self, p, l, width, polarization, power=1, theta_max=np.pi/2, phase=0, center=None,
                theta=0, phi=0, standing=False):
        super().__init__(polarization=polarization, power=power, theta_max=theta_max,
                phase=phase, center=center, theta=theta, phi=phi, standing=standing)
        self.width = width
        self.p = p
        self.l = l

    def __repr__(self):
        return f'LG_beam(width={self.width}, p={self.p}, l={self.l}, polarization={self.polarization}, ' \
               f'power={self.power}, center={self.center}, theta={self.theta}, phi={self.phi})'

    def scalar_angular_spectrum(self, theta, phi, k):
        Lpl = eval_genlaguerre(self.p, abs(self.l), 0.5*(k*self.width*np.tan(theta))**2)
        exp = np.exp(-(k*self.width*np.tan(theta)/2)**2)
        phase = np.exp(1j*self.l*phi)
        amp = 1j*(k*self.width*np.tan(theta)/np.sqrt(2))**abs(self.l)

        return amp*Lpl*exp*phase

def azimuthal_beam(width, theta=0, phi=0, power=None, phase=0, center=None, theta_max=np.pi/2):
    """azimuthally polarized beam"""
    if power is None:
        power = 1.0

    HG_1 = hermite_gaussian_beam(1, 0, width, [0,-1],  theta=theta, phi=phi, 
                  power=power/2, phase=phase, center=center, theta_max=theta_max)
    HG_2 = hermite_gaussian_beam(0, 1, width, [1,0], theta=theta, phi=phi,
                  power=power/2, phase=phase, center=center, theta_max=theta_max)
    return HG_1 + HG_2

def radial_beam(width, theta=0, phi=0, power=None, phase=0, center=None, theta_max=np.pi/2):
    """radially polarized beam"""
    if power is None:
        power = 1.0

    HG_1 = hermite_gaussian_beam(1, 0, width, [-1,0],  theta=theta, phi=phi, 
                  power=power/2, phase=phase, center=center, theta_max=theta_max)
    HG_2 = hermite_gaussian_beam(0, 1, width, [0,-1], theta=theta, phi=phi,
                  power=power/2, phase=phase, center=center, theta_max=theta_max)
    return HG_1 + HG_2

def shear_beam(width, theta=0, phi=0, power=None, phase=0, center=None, theta_max=np.pi/2):
    """shear polarized beam"""
    if power is None:
        power = 1.0

    HG_1 = hermite_gaussian_beam(1, 0, width, [0,-1],  theta=theta, phi=phi, 
                  power=power/2, phase=phase, center=center, theta_max=theta_max)
    HG_2 = hermite_gaussian_beam(0, 1, width, [-1,0], theta=theta, phi=phi,
                  power=power/2, phase=phase, center=center, theta_max=theta_max)
    return HG_1 + HG_2
