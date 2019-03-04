"""
Test weakly focused beams by comparing the E/H field from two methods:
    (i) Directly integrating the far-field angular spectrum to obtain focal fields
    (ii) Using the analytic paraxial expression for the beam
"""

import numpy as np
import miepy
from miepy.constants import Z0
from math import factorial
from scipy.special import eval_genlaguerre, eval_hermite
import pytest

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

def gaussian_paraxial(x, y, z, k):
    E0 = 2/width*np.sqrt(Z0*power/np.pi)

    rho_sq = x**2 + y**2
    wav = 2*np.pi/k
    amp = E0*width/w(z, width, wav) * np.exp(-rho_sq/w(z, width, wav)**2)
    phase = k*z + k*rho_sq*Rinv(z, width, wav)/2 - gouy(z, width, wav)

    return amp*np.exp(1j*phase)

def test_gaussian_beam_weak_focusing():
    """Weak focus test for Gaussian beam"""
    source = miepy.sources.gaussian_beam(width=width, polarization=polarization, power=power)
    E1 = source.E_field(X, Y, Z, k, sampling=100)[0]
    E2 = gaussian_paraxial(X, Y, Z, k)

    H1 = source.H_field(X, Y, Z, k, sampling=100)[1]
    H2 = E2

    assert np.allclose(E1, E2, rtol=8e-4, atol=0), 'electric field'
    assert np.allclose(H1, H2, rtol=8e-4, atol=0), 'magnetic field'

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

def test_gaussian_beam_z_component():
    """Check z-component of weakly focused Gaussian beam in xy-plane"""
    source = miepy.sources.gaussian_beam(width=width, polarization=polarization, power=power)

    x = np.linspace(-1.5*width, 1.5*width, 5)
    y = np.linspace(-1.5*width, 1.5*width, 5)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    
    Ez1 = source.E_field(X, Y, Z, k)[2]

    Ex = gaussian_paraxial(X, Y, Z, k)
    Ez2 = -1j*2*X/(k*width**2)*Ex

    Hz1 = source.H_field(X, Y, Z, k)[2]

    assert np.allclose(Ez1, Ez2, rtol=2e-3, atol=1e-9), 'Ez component matches'
    assert np.allclose(Hz1, Ez2.T, rtol=2e-3, atol=1e-9), 'Hz component matches (Hz(x,y) = Ez(y,x))'

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

    source = miepy.sources.hermite_gaussian_beam(1, 0, width=width, polarization=polarization, power=power)
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

    source = miepy.sources.laguerre_gaussian_beam(1, 1, width=width, polarization=polarization, power=power)
    E1 = source.E_field(X, Y, Z, k, sampling=100)[0]
    E2 = laguerre_gaussian_paraxial(source, X, Y, Z, k)

    assert np.allclose(E1, E2, rtol=4e-3, atol=1e-9)

def test_stiching_by_expansion_gaussian_beam():
    """reconstruct the paraxial E-field by expanding over p_src in a 'stitched' together fashion"""

    source = miepy.sources.gaussian_beam(width=width, polarization=polarization, power=power)
    E1 = gaussian_paraxial(X, Y, Z, k)

    E2 = np.zeros_like(E1)
    for xi,xval in enumerate(x):
        for yi,yval in enumerate(y):
            for zi,zval in enumerate(z):
                z0 = 50*nm
                p_src = source.structure([xval,yval,zval - z0], k, lmax=3)
                Efunc = miepy.expand_E(p_src, k, miepy.vsh_mode.incident)

                R, THETA, PHI =miepy.coordinates.cart_to_sph(0, 0, z0)
                E = Efunc(R, THETA, PHI)
                E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
                E2[xi,yi,zi] = E[0]

    assert np.allclose(E1, E2, rtol=8e-3, atol=0)
