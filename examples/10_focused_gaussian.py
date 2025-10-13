import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

import miepy

nm = 1e-9

width = 10 * nm
size = 600 * nm
wavelength = 600 * nm
k = 2 * np.pi / wavelength
lmax = 6

source = miepy.sources.gaussian_beam(width, [1, 1j])

Nx = 50
x = np.linspace(-size / 2, size / 2, Nx)
y = np.linspace(-size / 2, size / 2, Nx)
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

E = source.E_field(X, Y, Z, k)
I = np.sum(E.real**2, axis=0)
vmax = np.max(np.sum(np.abs(E) ** 2, axis=0)) / 2

fig, axes = plt.subplots(ncols=3, figsize=(12, 4.4), sharey=True)

ax = axes[0]
ax.pcolormesh(X / nm, Y / nm, np.sum(np.abs(E) ** 2, axis=0), shading="gouraud", rasterized=True, vmin=0)
ax.set_aspect("equal")
ax.set_title("time averaged energy", weight="bold")

ax = axes[1]
im = ax.pcolormesh(X / nm, Y / nm, I.T, shading="gouraud", rasterized=True, vmax=vmax, vmin=0)
ax.set_aspect("equal")
ax.set_title("instantaneous energy", weight="bold")


def update(phase):
    Enow = (E * (np.exp(1j * phase))).real
    I = np.sum(Enow**2, axis=0)
    im.set_array(np.ravel(I.T))
    return [im]


ani = animation.FuncAnimation(fig, update, np.linspace(0, 2 * np.pi, 120), interval=15, blit=True)

ax = axes[2]
H = source.H_field(X, Y, Z, k)

S = np.real(np.cross(E, np.conj(H), axis=0))
skip = 4
idx = np.s_[::skip, ::skip]
arrows = ax.quiver(X[idx] / nm, Y[idx] / nm, S[0][idx], S[1][idx], pivot="mid")
ax.set_aspect("equal")
ax.set_title("time averaged Poynting vector", weight="bold")

for ax in axes:
    ax.set(xlim=[-size / 2 / nm, size / 2 / nm], ylim=[-size / 2 / nm, size / 2 / nm], xlabel="x (nm)")
axes[0].set(ylabel="y (nm)")
fig.suptitle("Orbital angular momentum in a tightly focused RHC polarized Gaussian beam", fontsize=16)

plt.show()
