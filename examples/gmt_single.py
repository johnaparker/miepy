import matplotlib.pyplot as plt
import numpy as np

import miepy

nm = 1e-9
r = 20 * nm

fig, axes = plt.subplots(ncols=2, figsize=plt.figaspect(1 / 2))

Nx = 50
Ny = 50
x = np.linspace(-5 * r, 5 * r, Nx)
y = np.linspace(-5 * r, 5 * r, Ny)
z = np.array([0])
X, Y, Z = np.meshgrid(x, y, z, indexing="xy")
R = (X**2 + Y**2 + Z**2) ** 0.5
THETA = np.arccos(Z / R)
PHI = np.arctan2(Y, X)

system = miepy.sphere_cluster(
    position=[0, 0, 0],
    radius=r,
    material=miepy.constant_material(1.3),
    source=miepy.sources.plane_wave.from_string(polarization="x"),
    wavelength=600 * nm,
    lmax=2,
    interactions=False,
)

E = np.squeeze(system.E_field(X, Y, Z, source=False))
I = np.sum(np.abs(E) ** 2, axis=0)
print(E[:, 0, 0])

mask = np.zeros((Nx, Ny), dtype=bool)
mask[(np.squeeze(Y)) ** 2 + np.squeeze(X) ** 2 < 3.5 * r**2] = True
I[mask] = 0

im = axes[0].pcolormesh(np.squeeze(X) / nm, np.squeeze(Y) / nm, I, shading="gouraud", rasterized=True)
arrow = E[...]
arrow[0][mask] = 0
arrow[1][mask] = 0
skip = 2
axes[0].quiver(
    np.squeeze(X)[::skip, ::skip] / nm,
    np.squeeze(Y)[::skip, ::skip] / nm,
    np.real(arrow[0])[::skip, ::skip],
    np.real(arrow[1])[::skip, ::skip],
    pivot="mid",
)
plt.colorbar(im, ax=axes[0], label="Intensity")
axes[0].set(aspect="equal", xlabel="x (nm)", ylabel="y (nm)")

system = miepy.single_mie_sphere(r, miepy.constant_material(1.3), 600 * nm, 2)

E = np.squeeze(system.E_field(index=0)(R, THETA, PHI))

E_cart = np.array(
    [
        E[0] * np.sin(THETA.squeeze()) * np.cos(PHI.squeeze())
        + E[1] * np.cos(THETA.squeeze()) * np.cos(PHI.squeeze())
        - E[2] * np.sin(PHI.squeeze()),
        E[0] * np.sin(THETA.squeeze()) * np.sin(PHI.squeeze())
        + E[1] * np.cos(THETA.squeeze()) * np.sin(PHI.squeeze())
        + E[2] * np.cos(PHI.squeeze()),
        E[0] * np.cos(THETA.squeeze()) - E[1] * np.sin(THETA.squeeze()),
    ]
)
E = E_cart
# E[0] += 1

I = np.sum(np.abs(E) ** 2, axis=0)
print(E[:, 0, 0])
I[mask] = 0

im2 = axes[1].pcolormesh(np.squeeze(X) / nm, np.squeeze(Y) / nm, I, shading="gouraud", rasterized=True)
arrow = E[...]
arrow[0][mask] = 0
arrow[1][mask] = 0
axes[1].quiver(
    np.squeeze(X)[::skip, ::skip] / nm,
    np.squeeze(Y)[::skip, ::skip] / nm,
    np.real(arrow[0])[::skip, ::skip],
    np.real(arrow[1])[::skip, ::skip],
    pivot="mid",
)
plt.colorbar(im2, ax=axes[1], label="Intensity")
axes[1].set(aspect="equal", xlabel="x (nm)", ylabel="y (nm)")

plt.show()
