"""
Comapre analytic force and torque expressions to integrated Maxwell stress tensor
"""

import numpy as np
import miepy
from math import factorial
from tqdm import tqdm
from scipy import constants

nm = 1e-9

Ag = miepy.materials. Ag()
radius = 75*nm
source = miepy.sources.plane_wave.from_string(polarization='rhc')
separations = np.linspace(2*radius + 10*nm, 2*radius + 700*nm, 5)

analytic_force = np.zeros((3,) + separations.shape)
analytic_torque = np.zeros_like(analytic_force)

mst_force = np.zeros_like(analytic_force)
mst_torque = np.zeros_like(analytic_force)

for i, separation in enumerate(tqdm(separations)):
    mie = miepy.sphere_cluster(position=[[separation/2,0,0], [-separation/2,0,0]],
                               radius=radius,
                               material=Ag,
                               source=source,
                               wavelength=600*nm,
                               lmax=2)

    analytic_force[:,i] = mie.force_on_particle(0).squeeze()
    analytic_torque[:,i] = mie.torque_on_particle(0).squeeze()

    mst_force[:,i], mst_torque[:,i] = map(np.squeeze, 
            miepy.forces._gmt_force_and_torque_from_mst(mie, 0, sampling=30))

def test_force():
    """comapre MST force to analytic force"""
    L2 = np.linalg.norm(analytic_force - mst_force, axis=1)/analytic_force.shape[1]
    avg = np.average(np.abs(analytic_force) + np.abs(mst_force), axis=1)/2

    assert np.all(L2 < 1e-3*avg)

def test_torque():
    """comapre MST torque to analytic torque"""
    L2 = np.linalg.norm(analytic_torque - mst_torque, axis=1)/analytic_torque.shape[1]
    avg = np.average(np.abs(analytic_torque) + np.abs(mst_torque), axis=1)/2

    assert np.all(L2 < 1e-1*avg)



if __name__ == '__main__':
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=plt.figaspect(2/3)*2)

    for i in range(3):
        comp = ['x', 'y', 'z'][i]

        axes[0,i].plot(separations/nm, mst_force[i], 'o', color='C{}'.format(i), label='Numerical Stress Tensor')
        axes[1,i].plot(separations/nm, mst_torque[i], 'o', color='C{}'.format(i), label='Numerical Stress Tensor')

        axes[0,i].plot(separations/nm, analytic_force[i], color='C{}'.format(i), label='Analytic Equation')
        axes[1,i].plot(separations/nm, analytic_torque[i], color='C{}'.format(i), label='Analytic Equation')

        axes[0,i].set_title(label='F{comp}'.format(comp=comp), weight='bold')
        axes[1,i].set_title(label='T{comp}'.format(comp=comp), weight='bold')

    for ax in axes.flatten():
        ax.axhline(y=0, color='black', linestyle='--')
        ax.legend()

    for ax in axes[1]:
        ax.set_xlabel('separation (nm)')

    axes[0,0].set_ylabel('force (N)')
    axes[1,0].set_ylabel('torque (mN)')

    plt.tight_layout()
    plt.show()
