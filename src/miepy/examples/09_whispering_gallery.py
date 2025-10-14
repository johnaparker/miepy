"""demonstration of a whispering gallery mode in a single dielectric sphere excited by a point-dipole near-field."""

import matplotlib.pyplot as plt
import numpy as np

import miepy

nm = 1e-9

f = miepy.sources.point_dipole([700 * nm, 1100 * nm, 0], direction=[0, 0, 1])
x = np.linspace(-500 * nm, 500 * nm, 50)
y = np.linspace(200 * nm, 1500 * nm, 50)
X, Y = np.meshgrid(x, y, indexing="ij")
Z = np.zeros_like(X)
k = 2 * np.pi / (400 * nm)

E = f.E_field(X, Y, Z, k).squeeze()

I = np.sum(np.abs(E) ** 2, axis=0)

fig, ax = plt.subplots()
vmax = np.max(E[2].real)
ax.pcolormesh(X / nm, Y / nm, E[2].real, cmap="bwr", shading="gouraud", vmax=vmax, vmin=-vmax)
ax.set_aspect("equal")

radius = 250 * nm
source = miepy.sources.point_dipole([-radius - 10 * nm, 0, 0], direction=[0, 0, 1])
# source = miepy.sources.plane_wave.from_string(polarization='x')
sphere = miepy.sphere_cluster(
    position=[0, 0, 0],
    material=miepy.constant_material(4**2),
    # material=miepy.materials. Au(),
    radius=radius,
    source=source,
    lmax=4,
    wavelength=400 * nm,
)

x = np.linspace(-radius - 50 * nm, radius + 50 * nm, 150)
z = np.linspace(-radius - 50 * nm, radius + 50 * nm, 150)
X, Z = np.meshgrid(x, z, indexing="ij")
Y = np.zeros_like(X)

E = sphere.E_field(X, Y, Z)
# E = sphere.source.E([X,Y,Z], sphere.material_data.k_b).squeeze()
I = np.sum(np.abs(E) ** 2, axis=0)

fig, ax = plt.subplots()
# im = ax.pcolormesh(X/nm, Z/nm, I, shading='gouraud')
# plt.colorbar(im)
vmax = np.max(E[0].real) / 30000
ax.pcolormesh(X / nm, Z / nm, E[0].real, cmap="bwr", shading="gouraud", vmax=vmax, vmin=-vmax)
ax.set_aspect("equal")

plt.show()
