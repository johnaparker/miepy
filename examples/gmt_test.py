"""
Comparison of single Mie theory and GMT
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
from tqdm import tqdm

# wavelength from 400nm to 1000nm
wavelength = np.linspace(400e-9,1000e-9,1000)

# create a material with n = 3.7 (eps = n^2) at all wavelengths
dielectric = miepy.material_functions.constant_material(3.7**2)

# calculate scattering coefficients
radius = 100e-9       # 100 nm radius

# Single Mie Theory
Lmax = 10       # Use up to 10 multipoles
sphere = miepy.single_mie_sphere(radius, dielectric, wavelength, Lmax)
S,*_ = sphere.cross_sections()
plt.plot(wavelength*1e9,S,label="Single Mie theory", linewidth=2)


# Generalized Mie Theory (GMT)
particles = miepy.spheres(position=[[0,0,0]], radius=radius, material=dielectric)
# alternative call
# particles = miepy.spheres(positions=[[0,0,0]], radius=radius, material=dielectric)
wavelength = np.linspace(400e-9,1000e-9,100)
source = miepy.sources.x_polarized_plane_wave()
medium = None

system = miepy.gmt(particles, source, wavelength, 2, medium, interactions=True)
flux = system.flux_from_particle(0).squeeze()
plt.plot(wavelength*1e9, flux, label='GMT')

# Plot labels
plt.xlabel("wavelength (nm)")
plt.ylabel("Scattering cross-section")
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.legend()
plt.show()
