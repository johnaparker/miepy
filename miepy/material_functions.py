"""
materials defines the material class, containing eps and mu functions of wavelength
"""

import matplotlib.pyplot as plt
import numpy as np
import miepy
from scipy import constants
from scipy.interpolate import interp1d
import pandas as pd
from abc import ABCMeta, abstractmethod

class material:
    """material interface base class"""
    __metaclass__ = ABCMeta

    def __init__(self):
        pass

    @abstractmethod
    def eps(self, wavelength): pass

    @abstractmethod
    def mu(self, wavelength): pass

class constant_material(material):
    def __init__(self, eps, mu=1.0):
        """create a material with a constant eps and mu"""
        self.eps_value = eps
        self.mu_value = mu
    
    def eps(self, wavelength):
        f = np.vectorize(lambda wav: self.eps_value)
        return f(wavelength)

    def mu(self, wavelength):
        f = np.vectorize(lambda wav: self.mu_value)
        return f(wavelength)

class function_material(material):
    def __init__(self, eps_function, mu_function=None):
        """create a material with an eps and mu function"""
        self.eps_function = eps_function

        if mu_function is None:
            self.mu_function = lambda wavelength: 1.0
        else:
            self.mu_function = mu_function
    
    def eps(self, wavelength):
        f = np.vectorize(lambda wav: self.eps_function(wav))
        return f(wavelength)

    def mu(self, wavelength):
        f = np.vectorize(lambda wav: self.mu_function(wav))
        return f(wavelength)

class data_material(material):
    def __init__(self, wavelength, eps, mu=None):
        """create a material with raw wavelength, eps, and mu data"""
        if np.isscalar(wavelength):
            raise ValueError("Material has only 1 data point. Use constant_material instead")

        self.data = pd.DataFrame()
        
        self.data['wavelength'] = np.asarray(wavelength, dtype=float)
        self.data['eps'] = np.asarray(eps, dtype=complex)
        if mu is None:
            self.data['mu'] = np.ones_like(self.data['eps'])
        else:
            self.data['mu'] = np.asarray(mu, dtype=complex)

    def eps(self, wavelength):
        f = interp1d(self.data['wavelength'], self.data['eps'], kind='cubic')
        return f(wavelength)

    def mu(self, wavelength):
        f = interp1d(self.data['wavelength'], self.data['mu'], kind='cubic')
        return f(wavelength)


def wavelength_to_energy(wavelength):
    """return the wavelength in energy units (eV)"""
    return constants.h*constants.c/wavelength/constants.e

def wavelength_to_wavenumber(wavelength):
    """return the wavelength as a wavenumber"""
    return 2*np.pi/wavelength

#TODO implement
def plot_material(mat):
    """Plot material properties of 'mat' """
    pass
    
#TODO fix and test
def drude_lorentz(wp, sig, f, gam, magnetic_only=False, eps_inf=1):
    """Create a function_material using a Drude-Lorentz function.
       All arguments must be in eV units

       Arguments can be either 1D or 2D arrays: if 2D, eps/mu
       are both specified. If 1D, eps or mu is specified,
       depending on the variable magnetic_only.

            wp    =  plasma frequency      (eV), [2] array
            sig   =  strength factors      (eV), [2, #poles] array
            f     =  resonant frequencies  (eV), [2, #poles] array
            gam   =  damping factors       (eV), [2, #poles] array
            magnetic_only   =  True if 1D array is to specify mu
    """

    #convert iterable input to numpy arrays if necessary
    sig = np.asarray(sig)
    f = np.asarray(f)
    gam = np.asarray(gam)
    wav = np.asarray(wav)

    if len(sig.shape) == 1:
        size = len(sig)
        c = np.array([0,1]) if magnetic_only else np.array([1,0])
        wp = np.array([wp,wp])*c
        sig = np.array([sig,sig])
        f = np.array([f,f])
        gam = np.array([gam,gam])

    Nfreq = len(wav)
    size = len(sig[0])
    omega = 2*np.pi*constants.c*constants.hbar/(constants.e*wav*1e-9)

    eps = eps_inf*np.ones(Nfreq, dtype=np.complex)
    mu = np.ones(Nfreq, dtype=np.complex)
    for i in range(size):
        eps_add = sig[0,i]*wp[0]**2/(f[0,i]**2 - omega**2 -1j*omega*gam[0,i])
        mu_add = sig[1,i]*wp[1]**2/(f[1,i]**2 - omega**2 -1j*omega*gam[1,i])
        eps += eps_add
        mu += mu_add

    return data_material(wav,eps,mu)