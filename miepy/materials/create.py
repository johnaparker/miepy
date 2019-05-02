"""
creation new materials
"""

import numpy as np
import miepy
from scipy import constants
from scipy.interpolate import interp1d
from abc import ABCMeta, abstractmethod

def wavelength_to_energy(wavelength):
    """return the wavelength in energy units (eV)"""
    return constants.h*constants.c/wavelength/constants.e

def wavelength_to_wavenumber(wavelength):
    """return the wavelength as a wavenumber"""
    return 2*np.pi/wavelength

class material:
    """material interface base class"""
    __metaclass__ = ABCMeta

    def __init__(self, name=None):
        self.name = name

    @abstractmethod
    def eps(self, wavelength): pass

    @abstractmethod
    def mu(self, wavelength): pass

    def index(self, wavelength):
        return np.sqrt(self.eps(wavelength)*self.mu(wavelength))

    def __repr__(self):
        return self.name

class constant_material(material):
    def __init__(self, eps=1.0, mu=1.0, index=None, name=None):
        """create a material with a constant eps and mu"""
        super().__init__(name)

        self.eps_value = eps
        self.mu_value = mu

        self.f_eps = np.vectorize(lambda wav: self.eps_value)
        self.f_mu  = np.vectorize(lambda wav: self.mu_value)

        if index is not None:
            self.eps_value = index**2
    
    def eps(self, wavelength):
        return self.f_eps(wavelength)

    def mu(self, wavelength):
        return self.f_mu(wavelength)

    def __repr__(self):
        if self.name is None:
            return f'constant_material(eps={self.eps_value}, mu={self.mu_value})'
        else:
            return super().__repr__()

class function_material(material):
    def __init__(self, eps_function, mu_function=None, name=None):
        """create a material with an eps and mu function"""
        super().__init__(name)

        self.eps_function = eps_function

        if mu_function is None:
            self.mu_function = lambda wavelength: 1.0
        else:
            self.mu_function = mu_function

        self.f_eps = np.vectorize(lambda wav: self.eps_function(wav))
        self.f_mu = np.vectorize(lambda wav: self.mu_function(wav))
    
    def eps(self, wavelength):
        return self.f_eps(wavelength)

    def mu(self, wavelength):
        return self.f_mu(wavelength)

    def __repr__(self):
        if self.name is None:
            return f'function_material(eps={self.eps_function.__name__}, mu={self.mu_function.__name__})'
        else:
            return super().__repr__()

class data_material(material):
    def __init__(self, wavelength, eps, mu=None, name=None):
        """create a material with raw wavelength, eps, and mu data"""
        super().__init__(name)

        if np.isscalar(wavelength):
            raise ValueError("Material has only 1 data point. Use constant_material instead")

        self.wavelength = np.asarray(wavelength, dtype=float)
        self.eps_data = np.asarray(eps, dtype=complex)
        if mu is None:
            self.mu_data = np.ones_like(self.eps_data)
        else:
            self.mu_data = np.asarray(mu, dtype=complex)

        self.f_eps = interp1d(self.wavelength, self.eps_data, kind='linear')
        self.f_mu  = interp1d(self.wavelength, self.mu_data,  kind='linear')

    def eps(self, wavelength):
        return self.f_eps(wavelength)

    def mu(self, wavelength):
        return self.f_mu(wavelength)

    def __repr__(self):
        if self.name is None:
            return f'date_material(wavelengh_range=[{self.wavelength[0]:.2e}, {self.wavelength[-1]:.2e}])'
        else:
            return super().__repr__()

#TODO fix and test
def drude_lorentz(wp, sig, f, gam, magnetic_only=False, eps_inf=1, name=None):
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

    return data_material(wav,eps,mu, name=name)
