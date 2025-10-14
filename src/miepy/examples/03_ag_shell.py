"""Example of how to make a core-shell and plot scattering,
absorption, and scattering per multipole.
"""

import matplotlib.pyplot as plt
import numpy as np

import miepy

# wavelength from 400nm to 1000nm
wavelengths = np.linspace(400e-9, 1000e-9, 1000)

# Ag shell and dielectric core
Ag = miepy.materials.Ag()
dielectric = miepy.constant_material(1.46**2)

# Calculate scattering coefficients
radius_in = 135e-9
radius_out = 145e-9
lmax = 10  # Use up to 10 multipoles
core_shell = miepy.single_mie_core_shell(radius_in, radius_out, dielectric, Ag, wavelengths, lmax)

# Figure 1: Scattering and Absorption
fig, ax1 = plt.subplots()
C = core_shell.cross_sections()
plt.plot(wavelengths * 1e9, C.scattering, label="Scattering", linewidth=2)
plt.plot(wavelengths * 1e9, C.absorption, label="Absorption", linewidth=2)

# Figure 2: Scattering per multipole
fig, ax2 = plt.subplots()
plt.plot(wavelengths * 1e9, C.scattering, label="Total", linewidth=2)
S = core_shell.cross_sections_per_multipole().scattering

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
