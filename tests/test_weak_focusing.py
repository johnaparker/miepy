"""
Test tighly focused beams by comparing the E field from two methods:
    (i) Directly integrating the far-field angular spectrum to obtain focal fields
    (ii) Using the expansion coefficients around the center of the beam
"""

import numpy as np
import miepy
from miepy.constants import Z0

nm = 1e-9
wav = 600*nm
k = 2*np.pi/wav

width = 10000*nm
polarization = [1,0]
power = 1

### grid points to evaluate field over
x = np.linspace(-.5*width, .5*width, 5)
y = np.linspace(-.5*width, .5*width, 5)
z = np.linspace(-.5*width, .5*width, 3)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')


def zr(w0, wav):
    return np.pi*w0**2/wav

def w(z, w0, wav):
    return w0*np.sqrt(1 + (z/zr(w0,wav))**2)
    
def Rinv(z, w0, wav):
    return z/(z**2 + zr(w0,wav)**2)

def gouy(z, w0, wav):
    return np.arctan2(z, zr(w0,wav))

def test_gaussian_beam_weak_focusing():
    """Weak focus test for Gaussian beam"""
    def gaussian_paraxial(x, y, z, k):
        E0 = 2/width*np.sqrt(Z0*power/np.pi)

        rho_sq = x**2 + y**2
        wav = 2*np.pi/k
        amp = E0*width/w(z, width, wav) * np.exp(-rho_sq/w(z, width, wav)**2)
        phase = k*z + k*rho_sq*Rinv(z, width, wav)/2 - gouy(z, width, wav)

        return amp*np.exp(1j*phase)

    source = miepy.sources.gaussian_beam(width=width, polarization=polarization, power=power)
    E1 = source.E_field(X, Y, Z, k, sampling=100)[0]
    E2 = gaussian_paraxial(X, Y, Z, k)

    assert np.allclose(E1, E2, rtol=8e-4, atol=0)

def test_gouy_phase():
    """Weak focus test for Gaussian beam fails without the Guoy phase term"""
    def gaussian_paraxial(x, y, z, k):
        E0 = 2/width*np.sqrt(Z0*power/np.pi)

        rho_sq = x**2 + y**2
        wav = 2*np.pi/k
        amp = E0*width/w(z, width, wav) * np.exp(-rho_sq/w(z, width, wav)**2)
        phase = k*z + k*rho_sq*Rinv(z, width, wav)/2

        return amp*np.exp(1j*phase)

    source = miepy.sources.gaussian_beam(width=width, polarization=polarization, power=power)
    E1 = source.E_field(X, Y, Z, k, sampling=100)[0]
    E2 = gaussian_paraxial(X, Y, Z, k)

    assert not np.allclose(E1, E2, rtol=8e-3, atol=0)

# def test_E_field_weak_focusing():
    # def hermite_gaussian_paraxial(x, y, z, k):
        # E0 = 2/width*np.sqrt(Z0*power/np.pi)

        # rho_sq = x**2 + y**2
        # wav = 2*np.pi/k
        # amp = E0*width/w(z, width, wav) * np.exp(-rho_sq/w(z, width, wav)**2)
        # phase = k*z + k*rho_sq*Rinv(z, width, wav)/2 - gouy(z, width, wav)

        # return amp*np.exp(1j*phase)

    # source = miepy.sources.gaussian_beam(width=width, polarization=polarization, power=power)
    # E1 = source.E_field(X, Y, Z, k, sampling=100)[0]
    # E2 = gaussian_paraxial(X, Y, Z, k)

    # assert np.allclose(E1, E2, rtol=8e-4, atol=0)
