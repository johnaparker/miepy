"""
Example of how to make a dielectric material and plot scattering,
absorption, and scattering per multipole
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy

nm = 1e-9
# wavelength from 400nm to 1000nm
wavelengths = np.array([500*nm])

# create a material with n = 3.7 (eps = n^2) at all wavelengths
dielectric = miepy.constant_material(3.7**2)

# Calculate scattering coefficients
radius = 100e-9    # 100 nm radius
Lmax = 2          # Use up to 10 multipoles
sphere = miepy.single_mie_sphere(radius, dielectric, wavelengths, Lmax)

print(sphere.solve_exterior())
a,b = miepy.mie_single.mie_sphere_scattering_coefficients(radius, 2, 3.7**2, 1, 1, 1, 2*np.pi/(500*nm))
print(a,b)
