"""
Scattering, absroption, and extinction cross-sections of an Au dimer
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
import matplotlib.patheffects as path_effects


nm = 1e-9
um = 1e-6
Nwav = 30

Au = miepy.materials.predefined.Au()
radius = 50*nm

Lmax = 2
nb = 1.0
medium = miepy.constant_material(nb**2)

wavelengths = np.linspace(500*nm, 1000*nm, Nwav)
separation = 2*radius + 40*nm
source = miepy.sources.x_polarized_plane_wave()


spheres = miepy.spheres([[separation/2,0,0], [-separation/2,0,0]], radius, Au)
sol = miepy.gmt(spheres, source, wavelengths, Lmax, medium=medium)
scat, absorb, extinct = sol.cross_sections()

THETA, PHI = miepy.coordinates.sphere_mesh(30)
R = (separation/2 + radius + 20*nm)*np.ones_like(THETA)
X, Y, Z = miepy.coordinates.sph_to_cart(R, THETA, PHI)

E = sol.E_field(X,Y,Z, source=False)
H = sol.H_field(X,Y,Z, source=False)
C = np.array([miepy.flux.flux_from_poynting_sphere(E[:,i], H[:,i], R, eps=nb**2) for i in range(Nwav)])

E = sol.E_field(X,Y,Z, source=True)
H = sol.H_field(X,Y,Z, source=True)
A = np.array([-miepy.flux.flux_from_poynting_sphere(E[:,i], H[:,i], R, eps=nb**2) for i in range(Nwav)])

line, = plt.plot(wavelengths/nm, C, 'o',   color='C0', label='scattering')
line, = plt.plot(wavelengths/nm, A, 'o',   color='C1', label='absorbption')
line, = plt.plot(wavelengths/nm, C + A, 'o',   color='C2', label='extinction')

plt.axhline(color='black', linestyle='--')
line, = plt.plot(wavelengths/nm, scat,    color='C0', label='scattering')
line, = plt.plot(wavelengths/nm, absorb,  color='C1', label='absorbption')
line, = plt.plot(wavelengths/nm, extinct, color='C2', label='extinction')

plt.xlabel('wavelength (nm)')
plt.ylabel(r'cross-section ($\mu$m$^2$)')
plt.legend()

plt.show()
