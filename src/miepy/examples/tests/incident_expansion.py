"""Decompose an x-polarized plane wave into the VSHs."""

import matplotlib.pyplot as plt
import numpy as np

import miepy

### source definition
source = miepy.sources.plane_wave.from_string(polarization="x")
k = 2 * np.pi / 1
Nmax = 5

### grid plane
x = np.linspace(-0.3, 0.3, 30)
y = np.linspace(-0.3, 0.3, 30)
z = 0
X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
coords = np.array([X, Y, Z])


### exact fields
def exact():
    E = source.E(coords, k).squeeze()

    fig, ax = plt.subplots()
    I = np.sum(np.abs(E) ** 2, axis=0)
    im = ax.pcolormesh(X.squeeze(), Y.squeeze(), I, vmin=0.9, vmax=1.1)
    plt.colorbar(im)
    ax.quiver(X.squeeze(), Y.squeeze(), E[0].real, E[1].real)


### analytic expansion
def analytic():
    p, q = source.structure([0, 0, 0], k, Nmax)
    f = miepy.vsh.expand(p, q, k, miepy.vsh.vsh_mode.incident)
    expanded_E = f(X, Y, Z).squeeze()

    fig, ax = plt.subplots()
    I = np.sum(np.abs(expanded_E) ** 2, axis=0)
    im = ax.pcolormesh(X.squeeze(), Y.squeeze(), I)
    plt.colorbar(im)
    ax.quiver(X.squeeze(), Y.squeeze(), expanded_E[0].real, expanded_E[1].real)

    return p, q


### numerical decomposition and expansion
def numeric():
    p, q = miepy.vsh.decompose_source(source, k, Nmax, sampling=31)
    f = miepy.vsh.expand(p, q, k, miepy.vsh.vsh_mode.incident)
    expanded_E = f(X, Y, Z).squeeze()

    fig, ax = plt.subplots()
    I = np.sum(np.abs(expanded_E) ** 2, axis=0)
    im = ax.pcolormesh(X.squeeze(), Y.squeeze(), I)
    plt.colorbar(im)
    ax.quiver(X.squeeze(), Y.squeeze(), expanded_E[0].real, expanded_E[1].real)

    return p, q


exact()
p1, q1 = analytic()
p2, q2 = numeric()

err = np.sum(np.abs(p2 - p1) ** 2 + np.abs(q2 - q1) ** 2)
print(err)

plt.show()
