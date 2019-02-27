"""
Test tighly focused beams by comparing the E field from two methods:
    (i) Directly integrating the far-field angular spectrum to obtain focal fields
    (ii) Using the expansion coefficients around the center of the beam
"""

import numpy as np
import miepy
from miepy.constants import Z0
from math import factorial
from scipy.special import eval_genlaguerre, eval_hermite

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

def test_hermite_gaussian_beam_weak_focusing():
    """Weak focus test for Hermite-Gaussian beam"""
    def hermite_gaussian_paraxial(src, x, y, z, k):
        factor = 1/src.width*np.sqrt(2/(np.pi*2**src.l*2**src.m*factorial(src.l)*factorial(src.m)))
        E0 = factor*np.sqrt((2*Z0*src.power))

        rho_sq = x**2 + y**2
        wav = 2*np.pi/k

        wz = w(z, src.width, wav)
        HG_l = eval_hermite(src.l, np.sqrt(2)*x/wz)
        HG_m = eval_hermite(src.m, np.sqrt(2)*y/wz)
        N = src.l + src.m

        amp = E0*src.width/wz * HG_l * HG_m * np.exp(-rho_sq/wz**2)
        phase = k*z + k*rho_sq*Rinv(z,src.width,wav)/2 - (N+1)*gouy(z,src.width,wav)

        return amp*np.exp(1j*phase)

    source = miepy.sources.hermite_gaussian_beam(1, 0, width=width, polarization=polarization, power=power,
            theta_max=.15)
    E1 = source.E_field(X, Y, Z, k, sampling=100)[0]
    E2 = hermite_gaussian_paraxial(source, X, Y, Z, k)

    assert np.allclose(E1, E2, rtol=2e-3, atol=1e-9)

def test_laguerre_gaussian_beam_weak_focusing():
    """Weak focus test for Laguerre-Gaussian beam"""
    def laguerre_gaussian_paraxial(src, x, y, z, k):
        E0 = np.sqrt((2*Z0*src.power))

        rho_sq = x**2 + y**2
        phi = np.arctan2(y, x)
        wav = 2*np.pi/k

        C = np.sqrt(2*factorial(src.p)/(np.pi*factorial(src.p + abs(src.l))))
        wz = w(z, src.width, wav)

        Lpl = eval_genlaguerre(src.p, abs(src.l), 2*rho_sq/wz**2)
        N = abs(src.l) + 2*src.p

        amp = E0*C/wz * np.exp(-rho_sq/wz**2) * ((2*rho_sq)**0.5/wz)**abs(src.l) * Lpl
        phase = src.l*phi + k*z + k*rho_sq*Rinv(z,src.width,wav)/2 - (N+1)*gouy(z,src.width,wav)

        return amp*np.exp(1j*phase)

    source = miepy.sources.laguerre_gaussian_beam(1, 1, width=width, polarization=polarization, power=power, theta_max=.1)
    E1 = source.E_field(X, Y, Z, k, sampling=100)[0]
    E2 = laguerre_gaussian_paraxial(source, X, Y, Z, k)

    assert np.allclose(E1, E2, rtol=3e-3, atol=1e-9)
