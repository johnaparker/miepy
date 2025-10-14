"""Scattering, absroption, and extinction cross-sections of an Au dimer."""

import matplotlib.pyplot as plt
import numpy as np

import miepy

nm = 1e-9
um = 1e-6
Nwav = 60

Au = miepy.materials.Au()
radius = 40 * nm
source = miepy.sources.plane_wave.from_string(polarization="x")
lmax = 3

wavelengths = np.linspace(300 * nm, 800 * nm, Nwav)
separation = 83 * nm

scat = np.zeros_like(wavelengths)
absorb = np.zeros_like(wavelengths)
extinct = np.zeros_like(wavelengths)

for i, wavelength in enumerate(wavelengths):
    sol = miepy.sphere_cluster(
        position=[[separation / 2, 0, 0], [-separation / 2, 0, 0]],
        radius=radius,
        material=Au,
        source=source,
        wavelength=wavelength,
        lmax=lmax,
    )
    scat[i], absorb[i], extinct[i] = sol.cross_sections()

plt.figure(figsize=(8, 6))
plt.plot(wavelengths / nm, scat / um**2, label="scattering (dimer)", color="C0")
plt.plot(wavelengths / nm, absorb / um**2, label="absorption (dimer)", color="C1")
plt.plot(wavelengths / nm, extinct / um**2, label="extinction (dimer)", color="C2")

sphere = miepy.single_mie_sphere(radius, Au, wavelengths, lmax)
scat, absorb, extinct = sphere.cross_sections()

plt.plot(wavelengths / nm, 2 * scat / um**2, label="scattering (single x 2)", color="C0", linestyle="--")
plt.plot(wavelengths / nm, 2 * absorb / um**2, label="absorption (single x 2)", color="C1", linestyle="--")
plt.plot(wavelengths / nm, 2 * extinct / um**2, label="extinction (single x 2)", color="C2", linestyle="--")

plt.xlabel("wavelength (nm)")
plt.ylabel(r"cross-section ($\mu$m$^2$)")
plt.legend()

plt.show()
