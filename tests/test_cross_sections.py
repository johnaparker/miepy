"""
Compare analytic cross-sections to integrated Poynting vector for an Ag dimer in x-polarized light
"""

import numpy as np
import miepy
from tqdm import tqdm
from scipy import constants

nm = 1e-9

Ag = miepy.materials. Ag()
radius = 75*nm
source = miepy.sources.plane_wave.from_string(polarization='x')
separations = np.linspace(2*radius + 10*nm, 2*radius + 200*nm, 5)

analytic_scattering = np.zeros(separations.shape)
analytic_absorption = np.zeros(separations.shape)

poynting_scattering = np.zeros(separations.shape)
poynting_absorption = np.zeros(separations.shape)

for i, separation in enumerate(tqdm(separations)):
    mie = miepy.sphere_cluster(position=[[separation/2,0,0], [-separation/2,0,0]],
                               radius=radius,
                               material=Ag,
                               source=source,
                               wavelength=800*nm,
                               lmax=1)

    analytic_scattering[i], analytic_absorption[i], _ = mie.cross_sections()

    poynting_scattering[i], poynting_absorption[i], _ = miepy.flux._gmt_cross_sections_from_poynting(mie, radius=separation+radius, sampling=60)

def test_scattering():
    """comapre analytic scattering to numerical Poynting vector approach"""
    L2 = np.linalg.norm(analytic_scattering - poynting_scattering)/poynting_scattering.shape[0]
    avg = np.average(np.abs(analytic_scattering) + np.abs(poynting_scattering))/2

    assert L2 < 1e-5*avg

def test_absorption():
    """comapre analytic absorption to numerical Poynting vector approach"""
    L2 = np.linalg.norm(analytic_absorption - poynting_absorption)/poynting_absorption.shape[0]
    avg = np.average(np.abs(analytic_absorption) + np.abs(poynting_absorption))/2

    assert L2 < 2e-3*avg

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(ncols=2, figsize=plt.figaspect(1/2))

    axes[0].plot(separations/nm, analytic_scattering, color='C0', label='analytic')
    axes[1].plot(separations/nm, analytic_absorption, color='C0', label='analytic')

    axes[0].plot(separations/nm, poynting_scattering, 'o', color='C1', label='numerical poynting vector')
    axes[1].plot(separations/nm, poynting_absorption, 'o', color='C1', label='numerical poynting vector')

    axes[0].set_title('scattering', weight='bold')
    axes[1].set_title('absorption', weight='bold')

    for ax in axes:
        ax.legend()
        ax.set_xlabel('separation (nm)')

    axes[0].set_ylabel('cross-section')

    plt.tight_layout()
    plt.show()
