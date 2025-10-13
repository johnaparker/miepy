"""Scattering, absroption, and extinction cross-sections of an Au dimer."""

import matplotlib.pyplot as plt
import numpy as np

import miepy

nm = 1e-9
um = 1e-6

Nx = 150
Ny = 150
x = np.linspace(-1500 * nm, 1500 * nm, Nx)
y = np.linspace(-1500 * nm, 1500 * nm, Ny)
X, Y = np.meshgrid(x, y, indexing="ij")

R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, 0)

source = miepy.sources.azimuthal_beam(600 * nm)
lmax = 7

wavelength = 900 * nm
k = 2 * np.pi / wavelength

# idx_center = (X**2 + Y**2 < (20*nm)**2)


def plot(ax, E):
    I = np.sum(np.abs(E) ** 2, axis=0)
    # I[idx_center] = 0
    ax.pcolormesh(X / nm, Y / nm, I, rasterized=True)

    skip = 8
    np.s_[::skip, ::skip]
    # ax.quiver(X[idx]/nm, Y[idx]/nm, E[0][idx].real, E[1][idx].real, color='black', pivot='mid')


fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(14, 8))
E = source.E([X, Y, 0], k)
plot(axes[0, 0], E)

for i in range(1, lmax + 1):
    p, q = miepy.vsh.decompose_source(source, k, i)
    Efunc = miepy.vsh.expand(p, q, k, miepy.vsh.vsh_mode.incident)
    E = Efunc(X, Y, 0)
    print(E.shape)
    plot(axes.flatten()[i], E)

for ax in axes.flatten():
    ax.set_aspect("equal")

for ax in axes[1]:
    ax.set(xlabel="x (nm)")

for ax in axes[:, 0]:
    ax.set(ylabel="y (nm)")

axes[0, 0].set_title("azimuthal beam")
for i, ax in enumerate(axes.flatten()[1:]):
    ax.set_title(f"lmax = {i + 1}")

# plt.savefig('beam_expansion.pdf')
plt.show()
