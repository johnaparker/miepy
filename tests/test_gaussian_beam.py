import numpy as np
import miepy
import matplotlib.pyplot as plt

nm = 1e-9

wav = 600*nm
k = 2*np.pi/wav


def sim():
    s = miepy.sources.gaussian_beam(width=1900*nm, polarization=[1,1j], center=[100*nm, 0, 200*nm], theta=np.pi/4)

    x = np.linspace(-2900*nm, 2900*nm, 100)
    y = np.linspace(-2900*nm, 2900*nm, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    ### E by angular spectrum integration
    E = s.E_field(X, Z, Y, k)
    I = np.sum(np.abs(E)**2, axis=0)
    print(E[:,50,50])
    a1 = np.max(I)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(X/nm, Y/nm, I, shading='gouraud')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax)

    ### E by field expansion
    lmax = 6
    p_src = s.structure([0, 0, 0], k, lmax)
    Efunc = miepy.vsh.expand_E(p_src, k, miepy.vsh_mode.incident)
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Z, Y)
    E = Efunc(R, THETA, PHI)
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
    I = np.sum(np.abs(E)**2, axis=0)
    a2 = np.max(I)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(X/nm, Y/nm, I, shading='gouraud')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax)

    print(E[:,50,50])
    print(a2/a1)

    plt.show()
