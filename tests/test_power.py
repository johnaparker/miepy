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
def todo():
    k = 2*np.pi/(600*nm)
    width = 200*nm

    def A():
        source = miepy.sources.gaussian_beam(width, [1,0], power=1)

        def f(x, y):
            return 1e9*0.5*np.sum(np.abs(source.E_field(x, y, 0, k))**2)

        a = 20000*nm
        P = dblquad(f, -a, a, lambda x: -a, lambda x: a)
        print(P)

        def g(theta, phi):
            r = 1
            x, y, z = miepy.coordinates.sph_to_cart(r, theta, phi)
            return 1e9*0.5*np.sum(np.abs(source.E_field(x, y, z, k))**2)*r**2*np.sin(theta)

        def h(theta, phi):
            return np.exp(-2*(k*width*np.tan(theta)/2)**2)

        P = dblquad(h, np.pi/2, np.pi, lambda x: 0, lambda x: 2*np.pi)
        print(P[0]*5.5e-4)

    N, M = miepy.VSH(1, 0, mode=miepy.VSH_mode.incident)

    def f(x, y):
        r, theta, phi = miepy.coordinates.cart_to_sph(x, y, 1*nm*np.ones_like(x))
        E = N(r, theta, phi, k)
        H = M(r, theta, phi, k)
        S = 0.5*np.linalg.norm(np.cross(E, np.conjugate(H), axis=0), axis=0)
        return S

    # a = 68000*nm
    # x = np.linspace(-a, a, 2000)
    # y = np.linspace(-a, a, 2000)
    # X, Y = np.meshgrid(x, y)

    # S = f(X, Y)
    # P = trapz_2d(x, y, S)

    # print(P)
    # print(2*np.pi/k**2)


    # def f(x, y):
        # return 0.5*np.sum(np.abs(source.E_field(x, y, 0, k))**2)

    # a = 2000*nm
    # P = dblquad(f, -a, a, lambda x: -a, lambda x: a)
    # print(P)

    radius = 1e6*(2*np.pi/k)
    c = 0.5*(k*width)**2
    P = radius**2*np.pi*(1 - np.sqrt(np.pi*c)*np.exp(c)*erfc(np.sqrt(c)))
    print(P)


    lmax = 8
    source = miepy.sources.gaussian_beam(width, [1,0], power=1)
    p_src = source.structure([0,0,0], k, lmax, width)

    P = miepy.vsh.misc.power_through_aperature(source, [0,0,0], 3*width, k)
    print(P)

    theta = np.linspace(np.pi/2, np.pi, 30)
    phi = np.linspace(0, 2*np.pi, 60)
    THETA, PHI = np.meshgrid(theta, phi)
    R = radius*np.ones_like(THETA)

    Efunc = miepy.expand_E(p_src/2, k, miepy.VSH_mode.ingoing)
    E = Efunc(R, THETA, PHI)
    S = 0.5/Z0*np.sum(np.abs(E)**2, axis=0)*np.sin(THETA)
    P = radius**2*trapz_2d(theta, phi, S.real.T)
    print(P)

    E = Efunc(radius*np.ones_like(theta), theta, np.zeros_like(theta))
    plt.plot(theta, np.sum(np.abs(E)**2, axis=0))
    y = np.exp(-k**2*width**2*np.tan(theta)**2/2)
    y[theta<np.pi/2] = 0
    plt.plot(theta, y, 'o')

    factor = 0.5/Z0*np.pi/k**2
    P = factor*np.sum(np.abs(p_src)**2)
    print(P)

    a = 6000*nm
    x = np.linspace(-a, a, 150)
    y = np.linspace(-a, a, 150)
    X, Y = np.meshgrid(x, y)
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, 0)
    Efunc = miepy.expand_E(p_src, k, miepy.VSH_mode.incident)
    Hfunc = miepy.expand_H(p_src, k, miepy.VSH_mode.incident, 1, 1)
    E = Efunc(R, THETA, PHI)
    H = Hfunc(R, THETA, PHI)
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
    H = miepy.coordinates.vec_sph_to_cart(H, THETA, PHI)/Z0

    S = 0.5*np.linalg.norm(np.cross(E, np.conjugate(H), axis=0), axis=0)
    S = 0.5*np.cross(E, np.conjugate(H), axis=0)[2]
    # S = 0.5*np.sum(np.abs(E)**2, axis=0)

    P = trapz_2d(x, y, S.real)
    print(P)

    plt.show()
