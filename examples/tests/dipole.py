"""Displaying the fields in an xy cross section of the sphere (x polarized light, z-propagating)."""

import numpy as np

from miepy import sphere
from miepy.materials import material

# wavelength from 400nm to 1000nm
wav = np.linspace(300, 1100, 1000)

# create a material with n = 3.7 (eps = n^2) at all wavelengths
eps = 1.7**2 * np.ones(1000)
mu = 1 * np.ones(1000)
dielectric = material(wav, eps, mu)  # material object

# calculate scattering coefficients
rad = 20  # 100 nm radius
Nmax = 1  # Use up to 10 multipoles
m = sphere(Nmax, dielectric, rad)

E_func = m.E_field(999)
H_func = m.H_field(999)

E = -1 * E_func(R, THETA, PHI)
H = -H_func(R, THETA, PHI)

E = np.squeeze(E)
Ex = E[0] * np.sin(THETA) * np.cos(PHI) + E[1] * np.cos(THETA) * np.cos(PHI) - E[2] * np.sin(PHI)
Ey = E[0] * np.sin(THETA) * np.sin(PHI) + E[1] * np.cos(THETA) * np.sin(PHI) + E[2] * np.cos(PHI)
Ez = E[0] * np.cos(THETA) - E[1] * np.sin(THETA)

H = np.squeeze(H)
Hx = H[0] * np.sin(THETA) * np.cos(PHI) + H[1] * np.cos(THETA) * np.cos(PHI) - H[2] * np.sin(PHI)
Hy = H[0] * np.sin(THETA) * np.sin(PHI) + H[1] * np.cos(THETA) * np.sin(PHI) + H[2] * np.cos(PHI)
Hz = H[0] * np.cos(THETA) - H[1] * np.sin(THETA)
