"""
Test properties of a cluster of particles
"""

import numpy as np
import miepy
from tqdm import tqdm

nm = 1e-9

Ag = miepy.materials. Ag()
radius = 75*nm
source = miepy.sources.plane_wave.from_string(polarization='y')
lmax = 1

wavelengths = np.linspace(600*nm, 1000*nm, 5)


def test_off_center_particle(plot=False):
    """make sure scattering by a sphere away from the origin is equal to scattering of a particle at the origin"""

    # at the origin
    C1,A1,E1,C2,A2,E2 = [np.zeros_like(wavelengths, dtype=float) for i in range(6)]
    for i,wavelength in enumerate(wavelengths):
        sol = miepy.sphere_cluster(position=[0,0,0],
                                   radius=radius,
                                   material=Ag,
                                   source=source,
                                   wavelength=wavelength,
                                   lmax=lmax)

        C1[i], A1[i], E1[i] = sol.cross_sections()

        # displace sphere
        sol.update_position([40*nm,50*nm,60*nm])
        C2[i], A2[i], E2[i] = sol.cross_sections()

    if not plot:
        for a,b,tol in [(C1,C2,1e-15), (A1,A2,1e-14), (E1,E2,1e-15)]:
            L2 = np.linalg.norm(a - b)/a.shape[0]
            avg = np.average(a + b)/2
            print(L2, avg)
            assert L2 < tol*avg

    else:
        fig, ax = plt.subplots()
        ax.plot(wavelengths/nm, C1, color='C0', label='scattering (origin)')
        ax.plot(wavelengths/nm, A1, color='C1', label='absorption (origin)')
        ax.plot(wavelengths/nm, E1, color='C2', label='extinction (origin)')

        ax.plot(wavelengths/nm, C2, 'o', color='C0', label='scattering (displaced)')
        ax.plot(wavelengths/nm, A2, 'o', color='C1', label='absorption (displaced)')
        ax.plot(wavelengths/nm, E2, 'o', color='C2', label='extinction (displaced)')

        ax.set(xlabel='wavelength (nm)', ylabel='cross-section', title='test_off_center_particle')
        ax.legend()

def test_interactions_off(plot=False):
    """compare extinction of cluster to single Mie extinction with interactions disabled
            expect: E(cluster) = N*E(single) 
    """
    sep = 200*nm

    C1,A1,E1 = [np.zeros_like(wavelengths, dtype=float) for i in range(3)]
    for i,wavelength in enumerate(wavelengths):
        sol = miepy.sphere_cluster(position=[[sep/2,0,0], [-sep/2,0,0]],
                                   radius=radius,
                                   material=Ag,
                                   source=source,
                                   wavelength=wavelength,
                                   lmax=lmax,
                                   interactions=False)

        C1[i], A1[i], E1[i] = sol.cross_sections()

    # single
    sphere = miepy.single_mie_sphere(radius, Ag, wavelengths, lmax)
    Cs,As,Es = sphere.cross_sections()

    if not plot:
        for a,b,tol in [(E1,2*Es,1e-15)]:
            L2 = np.linalg.norm(a - b)/a.shape[0]
            avg = np.average(a + b)/2
            print(L2, avg)
            assert L2 < tol*avg

    else:
        fig, ax = plt.subplots()
        ax.plot(wavelengths/nm, C1, color='C0', label='scattering (cluster)')
        ax.plot(wavelengths/nm, A1, color='C1', label='absorption (cluster)')
        ax.plot(wavelengths/nm, E1, color='C2', label='extinction (cluster)')

        ax.plot(wavelengths/nm, 2*Cs, 'o', color='C0', label='scattering (N x single)')
        ax.plot(wavelengths/nm, 2*As, 'o', color='C1', label='absorption (N x single)')
        ax.plot(wavelengths/nm, 2*Es, 'o', color='C2', label='extinction (N x single)')

        ax.set(xlabel='wavelength (nm)', ylabel='cross-section', title='test_interactions_off')
        ax.legend()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_off_center_particle(plot=True)
    test_interactions_off(plot=True)
    plt.show()
