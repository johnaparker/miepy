"""Example of how to make a Ag sphere and plot scattering,
absorption, and scattering per multipole.
"""

import matplotlib.pyplot as plt
import numpy as np

import miepy

# wavelength from 400nm to 1000nm
wavelengths = np.linspace(400e-9, 800e-9, 1000)

# Ag material
Ag = miepy.materials.Ag()
water = miepy.materials.water()

# Calculate scattering coefficients
radius = 100e-9  # 100 nm radius
lmax = 10  # Use up to 10 multipoles
sphere = miepy.single_mie_sphere(radius, Ag, wavelengths, lmax, medium=water)

# Figure 1: Scattering and Absorption
fig, ax1 = plt.subplots()
C = sphere.cross_sections()
plt.plot(wavelengths * 1e9, C.scattering / C.scattering[-1], label="Scattering", linewidth=2)
plt.plot(wavelengths * 1e9, C.absorption / C.absorption[-1], label="Absorption", linewidth=2)

# Figure 2: Scattering per multipole
fig, ax2 = plt.subplots()
S = sphere.cross_sections_per_multipole().scattering
plt.plot(wavelengths * 1e9, C.scattering, label="Total", linewidth=2)

for i in range(2):
    for j in range(2):
        plt.plot(wavelengths * 1e9, S[:, i, j], label=miepy.multipole_label(i, j))

# Set figure properties
for ax in (ax1, ax2):
    ax.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    ax.legend()
    ax.set(xlabel="wavelength (nm)", ylabel=r"cross section (m$^2$)")
ax1.set_title("Total scattering")
ax2.set_title("Scattering per multipole")

plt.show()
