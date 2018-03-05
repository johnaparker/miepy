"""
Test properties of a cluster of particles
"""

import numpy as np
import h5py
import miepy
from tqdm import tqdm

nm = 1e-9

Ag = miepy.materials.predefined.Ag()
radius = 75*nm
source = miepy.sources.y_polarized_plane_wave()
Lmax = 8

sep = 600*nm
wavelength = np.linspace(400*nm, 1000*nm, 10)

spheres = miepy.spheres([[sep/2,0,0], [-sep/2,0,0]], radius, Ag)
mie = miepy.gmt(spheres, source, wavelength, Lmax, interactions=False)

C1,A1,E1 = mie.cross_sections()

sphere = miepy.single_mie_sphere(radius, Ag, wavelength, Lmax)
Cs,As,Es, *_ = sphere.cross_sections()

#TODO impolement
def test_particle_scattering():
    """compare scattering of cluster to single Mie scattering with interactions disabled"""
    pass

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    ax.plot(wavelength/nm, C1, color='C0', label='scattering')
    ax.plot(wavelength/nm, 2*Cs, 'o', color='C0')

    ax.plot(wavelength/nm, A1, color='C1', label='absorption')
    ax.plot(wavelength/nm, 2*As, 'o', color='C1')

    ax.plot(wavelength/nm, E1, color='C2', label='extinction')
    ax.plot(wavelength/nm, 2*Es, 'o', color='C2')

    ax.legend()
    ax.set(xlabel='wavelength (nm)')
    plt.show()
