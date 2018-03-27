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
dielectric = miepy.constant_material(3.7**2 + .1j)

# calculate scattering coefficients
radius = 100*nm       # 100 nm radius

# water medium
medium = miepy.materials.predefined.water()

# Single Mie Theory
Lmax = 5       # Use up to 5 multipoles
sphere = miepy.single_mie_sphere(radius, dielectric, wavelength, Lmax, medium=medium)
S,A,E = sphere.cross_sections()
Fz = sphere.radiation_force()

# Generalized Mie Theory (GMT)
particles = miepy.spheres(position=[[0,0,0]], radius=radius, material=dielectric)
source = miepy.sources.x_polarized_plane_wave()

system = miepy.gmt(particles, source, wavelength, Lmax, medium=medium)
scat,absorb,extinct = system.cross_sections()
force = system.force().squeeze()

def test_scattering():
    """compare scattering cross-section of GMT and single Mie theory"""
    L2 = np.linalg.norm(scat - S)/scat.shape[0]
    avg = np.average(np.abs(S) + np.abs(scat))/2

    assert np.all(L2 < 1e-15*avg)

def test_absoprtion():
    """compare absoprtion cross-section of GMT and single Mie theory"""
    L2 = np.linalg.norm(absorb - A)/scat.shape[0]
    avg = np.average(np.abs(A) + np.abs(absorb))/2

    assert np.all(L2 < 1e-15*avg)

def test_force():
    """compare radiation force of GMT and single Mie theory"""
    L2 = np.linalg.norm(force[2] - Fz)/Fz.shape[0]
    avg = np.average(np.abs(force[2]) + np.abs(Fz))/2

    assert np.all(L2 < 2e-8*avg)

if __name__ == '__main__':
    import matplotlib.pyplot as plt

    plt.figure()
    plt.plot(wavelength/nm, scat, color='C0', label='GMT scattering')
    plt.plot(wavelength/nm, S, 'o', color='C0', label="Single Mie theory", linewidth=2)

    plt.plot(wavelength/nm, absorb, color='C1', label='GMT absorption')
    plt.plot(wavelength/nm, A, 'o', color='C1', label="Single Mie theory", linewidth=2)

    plt.plot(wavelength/nm, extinct, color='C2', label='GMT extinction')
    plt.plot(wavelength/nm, E, 'o', color='C2', label="Single Mie theory", linewidth=2)

    plt.xlabel("wavelength (nm)")
    plt.ylabel("Scattering cross-section")
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()

    plt.figure()
    plt.plot(wavelength/nm, force[2], color='C1', label='GMT force')
    plt.plot(wavelength/nm, Fz, 'o', color='C1', label="Single Mie theory", linewidth=2)

    plt.show()
