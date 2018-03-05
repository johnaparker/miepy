"""
Comparison of single Mie theory and GMT
"""

import numpy as np
import miepy
from tqdm import tqdm

nm = 1e-9

# wavelength from 400nm to 1000nm
wavelength = np.linspace(400*nm,1000*nm,10)

# create a material with n = 3.7 (eps = n^2) at all wavelengths
dielectric = miepy.constant_material(3.7**2)

# calculate scattering coefficients
radius = 100*nm       # 100 nm radius

# Single Mie Theory
Lmax = 10       # Use up to 10 multipoles
sphere = miepy.single_mie_sphere(radius, dielectric, wavelength, Lmax)
S,*_ = sphere.cross_sections()

# Generalized Mie Theory (GMT)
particles = miepy.spheres(position=[[0,0,0]], radius=radius, material=dielectric)
source = miepy.sources.x_polarized_plane_wave()

system = miepy.gmt(particles, source, wavelength, 2)
scat,*_ = system.cross_sections()

def test_scattering():
    """compare scattering cross-section of GMT and single Mie theory"""
    L2 = np.linalg.norm(scat - S)/scat.shape[0]
    avg = np.average(np.abs(S) + np.abs(scat))/2

    assert np.all(L2 < 2e-3*avg)

test_scattering()
if __name__ == '__main__':
    import matplotlib.pyplot as plt

    plt.plot(wavelength/nm, S, 'o', label="Single Mie theory", linewidth=2)
    plt.plot(wavelength/nm, scat, label='GMT')

    plt.xlabel("wavelength (nm)")
    plt.ylabel("Scattering cross-section")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.show()
