"""
Tests for the power quantity of different incident beams
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
from scipy.integrate import nquad, dblquad
from scipy.special import erfc
from miepy.vsh.misc import trapz_2d
Z0 = 376.7

nm = 1e-9

#TODO implement tests
def test_numeric_power():
    """compute power of a Gaussian beam using various methods"""

    ### analytic expression
    k = 2*np.pi/(600*nm)
    width = 200*nm

    radius = 1e6*(2*np.pi/k)
    c = 0.5*(k*width)**2
    P = radius**2*np.pi*(1 - np.sqrt(np.pi*c)*np.exp(c)*erfc(np.sqrt(c)))
    P *= 6

    print(P)
    assert(abs(P-1) < 0.002)

    ### power through aperature
    lmax = 6
    source = miepy.sources.gaussian_beam(width, [1,0], power=1)
    p_src = source.structure([0,0,0], k, lmax, width)

    P = miepy.vsh.misc.power_through_aperature(source, [0,0,0], 3*width, k, sampling=20)

    print(P)
    assert(abs(P-1) < 0.03)

    ### integration over spherically ingoing waves
    theta = np.linspace(np.pi/2, np.pi, 20)
    phi = np.linspace(0, 2*np.pi, 20)
    THETA, PHI = np.meshgrid(theta, phi)
    R = radius*np.ones_like(THETA)

    Efunc = miepy.expand_E(p_src/2, k, miepy.vsh_mode.ingoing)
    E = Efunc(R, THETA, PHI)
    S = 0.5/Z0*np.sum(np.abs(E)**2, axis=0)*np.sin(THETA)
    P = radius**2*trapz_2d(theta, phi, S.T).real

    print(P)
    assert(abs(P-1) < 0.02)

    ### sum over expansion coefficients from analytic integration
    factor = 0.5/Z0*np.pi/k**2
    P = factor*np.sum(np.abs(p_src)**2)

    print(P)
    assert(abs(P-1) < 0.008)

    ### 2D integration of Poynting vector
    a = 3000*nm
    x = np.linspace(-a, a, 20)
    y = np.linspace(-a, a, 20)
    X, Y = np.meshgrid(x, y)
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, 0)
    Efunc = miepy.expand_E(p_src, k, miepy.vsh_mode.incident)
    Hfunc = miepy.expand_H(p_src, k, miepy.vsh_mode.incident, 1, 1)
    E = Efunc(R, THETA, PHI)
    H = Hfunc(R, THETA, PHI)
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
    H = miepy.coordinates.vec_sph_to_cart(H, THETA, PHI)/Z0

    S = 0.5*np.linalg.norm(np.cross(E, np.conjugate(H), axis=0), axis=0)
    S = 0.5*np.cross(E, np.conjugate(H), axis=0)[2]
    # S = 0.5*np.sum(np.abs(E)**2, axis=0)

    P = trapz_2d(x, y, S).real

    print(P)
    assert(abs(P-1) < 0.01)
