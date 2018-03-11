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
Lmax = 1

wavelength = np.linspace(600*nm, 1000*nm, 5)


def test_off_center_particle(plot=False):
    """make sure scattering by a sphere away from the origin is equal to scattering of a particle at the origin"""

    # at the origin
    sphere = miepy.spheres([0,0,0], radius, Ag)
    mie = miepy.gmt(sphere, source, wavelength, Lmax)
    C1,A1,E1 = mie.cross_sections()

    # displaced sphere
    mie.update_position([40*nm,50*nm,60*nm])
    C2,A2,E2 = mie.cross_sections(4)

    if not plot:
        for a,b,tol in [(C1,C2,1e-5), (A1,A2,1e-1), (E1,E2,1e-3)]:
            L2 = np.linalg.norm(a - b)/a.shape[0]
            avg = np.average(a + b)/2
            assert L2 < tol*avg

    else:
        fig, ax = plt.subplots()
        ax.plot(wavelength/nm, C1, color='C0', label='scattering (origin)')
        ax.plot(wavelength/nm, A1, color='C1', label='absorption (origin)')
        ax.plot(wavelength/nm, E1, color='C2', label='extinction (origin)')

        ax.plot(wavelength/nm, C2, 'o', color='C0', label='scattering (displaced)')
        ax.plot(wavelength/nm, A2, 'o', color='C1', label='absorption (displaced)')
        ax.plot(wavelength/nm, E2, 'o', color='C2', label='extinction (displaced)')

        ax.set(xlabel='wavelength (nm)', ylabel='cross-section', title='test_off_center_particle')
        ax.legend()

def test_interactions_off(plot=False):
    """compare scattering of cluster to single Mie scattering with interactions disabled
            expect: C(cluster) = N*C(single) 
    """
    sep = 200*nm

    # cluster
    spheres = miepy.spheres([[sep/2,0,0], [-sep/2,0,0]], radius, Ag)
    cluster = miepy.gmt(spheres, source, wavelength, Lmax, interactions=True)
    C1,A1,E1 = cluster.cross_sections(4)

    # single
    sphere = miepy.single_mie_sphere(radius, Ag, wavelength, Lmax)
    Cs,As,Es = sphere.cross_sections()

    if not plot:
        for a,b,tol in [(C1,Cs,1e-1), (A1,As,1e-1), (E1,Es,1e-1)]:
            L2 = np.linalg.norm(a - b)/a.shape[0]
            avg = np.average(a + b)/2
            assert L2 < tol*avg

    else:
        fig, ax = plt.subplots()
        ax.plot(wavelength/nm, C1, color='C0', label='scattering (cluster)')
        ax.plot(wavelength/nm, A1, color='C1', label='absorption (cluster)')
        ax.plot(wavelength/nm, E1, color='C2', label='extinction (cluster)')

        ax.plot(wavelength/nm, 2*Cs, 'o', color='C0', label='scattering (N x single)')
        ax.plot(wavelength/nm, 2*As, 'o', color='C1', label='absorption (N x single)')
        ax.plot(wavelength/nm, 2*Es, 'o', color='C2', label='extinction (N x single)')

        ax.set(xlabel='wavelength (nm)', ylabel='cross-section', title='test_interactions_off')
        ax.legend()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_off_center_particle(plot=True)
    test_interactions_off(plot=True)
    plt.show()
