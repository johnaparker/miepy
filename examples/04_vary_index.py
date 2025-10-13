"""Scattering intensity of a dielectric sphere for variable wavelength and dielectric constant."""

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import miepy

# Variable wavelengths and core index of refraction
N_index = 250
N_wav = 250

data = np.zeros(shape=(N_index, N_wav))
indices = np.linspace(1, 4, N_index)
wavelengths = np.linspace(400e-9, 1000e-9, N_wav)

# Calculate scattering coefficients
radius = 165e-9  # 165 nm radius
lmax = 10  # Use up to 10 multipoles

for i, index in enumerate(tqdm(indices)):
    dielectric = miepy.constant_material(index**2)
    sphere = miepy.single_mie_sphere(radius, dielectric, wavelengths, lmax)

    sphere.solve_exterior()
    S = sphere.cross_sections().scattering

    data[i] = S

# Plot results
plt.figure(1)
data /= np.max(data)

X, Y = np.meshgrid(indices, wavelengths * 1e9, indexing="ij")

plt.pcolormesh(X, Y, data, shading="gouraud")
plt.xlabel("index of refraction")
plt.ylabel("wavelength (nm)")
plt.colorbar(label="scattering cross-section")

plt.show()
