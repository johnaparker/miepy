import matplotlib.pyplot as plt
import numpy as np

import miepy

nm = 1e-9
sep = 50 * nm
r = 20 * nm

fig, axes = plt.subplots(nrows=2, figsize=plt.figaspect(2))

Nx = 180
Ny = 181
x = np.linspace(-sep / 2 - 2 * r, sep / 2 + 2 * r, Nx)
y = np.linspace(-2 * r, 2 * r, Ny)
z = 0
X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

for i, sim in enumerate([True, False]):
    system = miepy.sphere_cluster(
        position=[[-sep / 2, 0, 0], [sep / 2, 0, 0]],
        radius=r,
        material=miepy.materials.Ag(),
        source=miepy.sources.plane_wave.from_string(polarization="x"),
        wavelength=600 * nm,
        lmax=2,
        interactions=sim,
    )
    E = np.squeeze(system.E_field(X, Y, Z, source=True))
    I = np.sum(np.abs(E) ** 2, axis=0)
    print(I.shape)

    mask = np.zeros((Nx, Ny), dtype=bool)
    mask[(np.squeeze(X) - sep / 2) ** 2 + np.squeeze(Y) ** 2 < r**2] = True
    mask[(np.squeeze(X) + sep / 2) ** 2 + np.squeeze(Y) ** 2 < r**2] = True
    I[mask] = 0

    im = axes[i].pcolormesh(np.squeeze(X) / nm, np.squeeze(Y) / nm, I, shading="gouraud", rasterized=True)
    plt.colorbar(im, ax=axes[i], label="Intensity")
    axes[i].set(aspect="equal", xlabel="x (nm)", ylabel="y (nm)")


plt.show()
