import numpy as np
import miepy
import matplotlib.pyplot as plt

nm = 1e-9

wav = 600*nm
k = 2*np.pi/wav

### E by angular spectrum integration

def sim():
    s = miepy.sources.gaussian_beam(width=30*nm, polarization=[1,1])

    x = np.linspace(-900*nm, 900*nm, 100)
    y = np.linspace(-900*nm, 900*nm, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)

    E = s.E_field(X, Y, Z, k)
    I = np.sum(np.abs(E)**2, axis=0)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(X/nm, Y/nm, I, shading='gouraud')
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax)

    ### E by field expansion

    lmax = 6
    p_src = s.structure([0, 0, 0], k, lmax)
    Efunc = miepy.vsh.expand_E(p_src, k, miepy.vsh_mode.incident)
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z + 10*nm)
    E = Efunc(R, THETA, PHI)
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
    I = np.sum(np.abs(E)**2, axis=0)

    fig, ax = plt.subplots()
    im = ax.pcolormesh(X/nm, Y/nm, I)
    ax.set_aspect('equal')
    plt.colorbar(im, ax=ax)

    plt.show()
