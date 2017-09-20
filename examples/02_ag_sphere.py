"""
Example of how to make a Ag sphere and plot scattering,
absorption, and scattering per multipole
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy

# wavelength from 400nm to 1000nm
wavelengths = np.linspace(400e-9,1000e-9,1000)

# Ag material
Ag = miepy.materials.predefined.Ag()

# Calculate scattering coefficients
radius = 100e-9    # 100 nm radius
Lmax = 10          # Use up to 10 multipoles
sphere = miepy.single_mie_sphere(radius, Ag, wavelengths, Lmax)

# Figure 1: Scattering and Absorption
fig, ax1 = plt.subplots()
S,A,_ = sphere.cross_sections()
plt.plot(wavelengths*1e9, S, label="Scattering", linewidth=2)
plt.plot(wavelengths*1e9, A, label="Absorption", linewidth=2)

# Figure 2: Scattering per multipole
fig, ax2 = plt.subplots()
plt.plot(wavelengths*1e9, S, label="Total", linewidth=2)
S,*_ = sphere.cross_sections_per_multipole()

for i in range(2):
    for j in range(2):
        plt.plot(wavelengths*1e9, S[:,i,j], label=miepy.scattering.multipole_label(i,j))

# Set figure properties
for ax in (ax1,ax2):
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.legend()
    ax.set(xlabel="wavelength (nm)", ylabel=r"cross section (m$^2$)")
ax1.set_title("Total scattering")
ax2.set_title("Scattering per multipole")

plt.show()
