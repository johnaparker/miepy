"""Test to ensure that scattering cross-section computed via 2 different methods yields identical results:

1. Numerical integration of Poynting vector over sphere
2. Sum over a,b coefficients of the cluster
"""

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
Nwav = 60

Au = miepy.materials.Au()
radius = 40 * nm
source = miepy.sources.x_polarized_plane_wave(amplitude=1)
# separations = np.linspace(2*radius+10e-9,2*radius+700e-9, 50)
wavelengths = np.linspace(300 * nm, 800 * nm, Nwav)
separation = 83 * nm
flux = np.zeros_like(wavelengths)
flux2 = np.zeros_like(wavelengths)

plt.figure()
for wav_idx, wav in enumerate(tqdm(wavelengths)):
    spheres = miepy.spheres([[separation / 2, 0, 0], [-separation / 2, 0, 0]], radius, Au)
    sol = miepy.gmt(spheres, source, wav, 2, interactions=True)

    scat = sol.scattering_cross_section()
    flux2[wav_idx] = np.sum(scat)

    # sol.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))

    r = 10000 * nm
    THETA, PHI = sphere_mesh(sampling)
    X, Y, Z = sph_to_cart(r, THETA, PHI)

    E = sol.E_field(X, Y, Z, inc=False)
    # H = sol.H_field(X,Y,Z, inc=False)
    I = np.sum(np.abs(E) ** 2, axis=0).squeeze()

    dA = r**2
    tau = np.linspace(-1, 1, sampling)
    phi = np.linspace(0, 2 * np.pi, 2 * sampling)

    flux[wav_idx] = miepy.vsh.misc.simps_2d(tau, phi, I) * dA

plt.plot(wavelengths / nm, flux, color="C0")
plt.plot(wavelengths / nm, flux2, color="C1")

plt.xlabel("wavelength (nm)")
plt.ylabel("scattering")

plt.show()
