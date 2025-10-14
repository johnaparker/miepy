import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import miepy


def sph_to_cart(r, theta, phi, origin=None):
    """Convert spherical coordinates (r, theta, phi) centered at origin to cartesian coordinates (x, y, z)."""
    if origin is None:
        origin = [0, 0, 0]
    x = origin[0] + r * np.sin(theta) * np.cos(phi)
    y = origin[1] + r * np.sin(theta) * np.sin(phi)
    z = origin[2] + r * np.cos(theta)

    return x, y, z


def sphere_mesh(sampling):
    """Obtain a THETA,PHI mesh for discretizing the surface of the sphere, consistent
    with the format required by the project and decompose functions
    Returns (THETA,PHI) meshgrids.

    Arguments:
        sampling   number of points to sample between 0 and pi
    """
    phi = np.linspace(0, 2 * np.pi, 2 * sampling)
    tau = np.linspace(-1, 1, sampling)
    theta = np.arccos(tau)

    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    return THETA, PHI


nm = 1e-9
sampling = 40
Nwav = 200

Ag = miepy.materials.Ag()
radius = 75 * nm
source = miepy.sources.plane_wave.from_string(polarization="x", amplitude=1)
# separations = np.linspace(2*radius+10e-9,2*radius+700e-9, 50)
separation = 153 * nm
wavelengths = np.linspace(330 * nm, 1000 * nm, Nwav)

separations = np.linspace(153 * nm, 300 * nm, 10)

plt.figure(figsize=(5, 10))

r = 10000 * nm
THETA, PHI = sphere_mesh(sampling)
X, Y, Z = sph_to_cart(r, THETA, PHI)
tau = np.linspace(-1, 1, sampling)
phi = np.linspace(0, 2 * np.pi, 2 * sampling)

for sep_idx, separation in enumerate(tqdm(separations)):
    flux = np.zeros_like(wavelengths)

    for i, wavelength in enumerate(wavelengths):
        sol = miepy.sphere_cluster(
            position=[[separation / 2, 0, 0], [-separation / 2, 0, 0]],
            radius=radius,
            material=Ag,
            source=source,
            wavelength=wavelength,
            lmax=3,
            interactions=True,
        )

        E = sol.E_field(X, Y, Z, source=False)
        I = np.sum(np.abs(E) ** 2, axis=0)
        flux[i] = miepy.vsh.misc.simps_2d(tau, phi, I)

    plt.plot(wavelengths / nm, flux - 0.0004 * sep_idx, color="C0")

plt.xlabel("wavelength (nm)")
plt.ylabel("scattering")
plt.savefig("out.pdf", bbox_inches="tight")
plt.show()
