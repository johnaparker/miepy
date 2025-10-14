import matplotlib.pyplot as plt
import numpy as np

import miepy


def hexagonal_lattice_layers(L):
    """Return a hexagonal lattice with unit spacing and L layers."""
    k1 = np.array([1, 0, 0], dtype=float)
    k2 = np.array([np.cos(np.pi / 3), np.sin(np.pi / 3), 0], dtype=float)

    lattice = [np.zeros(3)]

    for step in range(1, L + 1):
        lattice.append(step * k1)

        for direc in [k2 - k1, -k1, -k2, -k2 + k1, k1, k2]:
            for _ in range(step):
                lattice.append(lattice[-1] + direc)
        lattice.pop()

    return np.asarray(lattice, dtype=float)


nm = 1e-9
um = 1e-6

theta_obj = 90 / 180 * np.pi

Ag = miepy.materials.Ag()
radius = 75 * nm
source = miepy.sources.plane_wave.from_string(polarization="rhc", direction="-z")
lmax = 2
water = miepy.materials.water()

wavelength = 800 * nm
separation = 600 * nm
L = 1
lattice = hexagonal_lattice_layers(L)

cluster = miepy.sphere_cluster(
    position=separation * lattice,
    radius=radius,
    material=Ag,
    source=source,
    wavelength=wavelength,
    lmax=lmax,
    medium=water,
)
xmax = L * separation + 8 * radius

x = np.linspace(-xmax, xmax, 100)
y = np.linspace(-xmax, xmax, 100)
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

E = cluster.E_field(X, Y, Z)
I = np.sum(np.abs(E) ** 2, axis=0)

fig, axes = plt.subplots(ncols=3, figsize=(15, 6))
ax = axes[0]
ax.pcolormesh(X / nm, Y / nm, I**0.5, shading="gouraud", rasterized=True)
ax.set_aspect("equal")
ax.set_title("near fields", weight="bold")

theta = np.linspace(0, theta_obj, 50)
phi = np.linspace(0, 2 * np.pi, 100)
THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
E_far = cluster.E_angular(THETA, PHI)

I_far = np.sum(np.abs(E_far) ** 2, axis=0)
ax = axes[1]

X = np.sin(THETA) * np.cos(PHI)
Y = np.sin(THETA) * np.sin(PHI)
ax.pcolormesh(X, Y, I_far, shading="gouraud", rasterized=True)
ax.set_aspect("equal")
ax.set_title("angular far fields", weight="bold")

ax = axes[2]
x = np.linspace(-xmax, xmax, 50)
y = np.linspace(-xmax, xmax, 50)
X, Y = np.meshgrid(x, y)

im = cluster.microscope(X, Y, theta_obj=theta_obj)
I = np.sum(np.abs(im) ** 2, axis=0)

ax.pcolormesh(X / nm, Y / nm, I, shading="gouraud", cmap="gray", rasterized=True)
ax.set_aspect("equal")
ax.set_title("image", weight="bold")

plt.show()
