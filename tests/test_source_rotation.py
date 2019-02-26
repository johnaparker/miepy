"""
Tests for rotation of incident source expansion coefficients
"""

import numpy as np
import miepy

nm = 1e-9
um = 1e-6

def test_plane_wave_rotation():
    """
    A rotated plane-wave: analytic compared to rotated expansion coefficients
    """
    Nx = 6
    Ny = 6
    Nz = 6
    x = np.linspace(-100*nm, 100*nm, Nx)
    y = np.linspace(-100*nm, 100*nm, Nx)
    z = np.linspace(-100*nm, 100*nm, Ny)

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    Y = np.zeros_like(X)

    wavelength = 600*nm
    k = 2*np.pi/wavelength

    source = miepy.sources.plane_wave(polarization=[1,0], theta=np.pi/3, phi=np.pi/2.7)
    E1 = source.E_field(X, Y, Z, k)

    lmax = 6
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)
    p_src = source.structure([0,0,0], k, lmax)
    Efunc = miepy.expand_E(p_src, k, mode=miepy.vsh_mode.incident)

    E2 = Efunc(R, THETA, PHI)
    E2 = miepy.coordinates.vec_sph_to_cart(E2, THETA, PHI)

    L2 = np.sqrt(np.sum(np.abs(E1 - E2)**2))/np.product(X.shape)
    avg = np.average(np.abs(E1) + np.abs(E2))/2
    assert L2 < 7e-6*avg
