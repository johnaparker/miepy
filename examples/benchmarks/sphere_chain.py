"""
performance test
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
from tqdm import tqdm
from functools import partial
from scipy.sparse.linalg import bicg, bicgstab
from miepy.interactions import solve_linear_system
from topics.photonic_clusters.create_lattice import hexagonal_lattice_particles
from timer import time_function

nm = 1e-9

Ag = miepy.materials.Ag()
radius = 75*nm
source = miepy.sources.plane_wave.from_string(polarization='rhc')
source = miepy.sources.gaussian_beam(2500*nm, [1, 1j], power=.002)
separation = 600*nm

def tests(Nmax, step=1):
    Nparticles = np.arange(1, Nmax+1, step)
    t_force, t_flux, t_build, t_solve, t_source = [np.zeros_like(Nparticles, dtype=float) for i in range(5)]
    for i,N in enumerate(Nparticles):
        print(N, Nmax)
        # positions = [[n*separation, 0, 0] for n in range(N)]
        positions = hexagonal_lattice_particles(N)*separation
        mie = miepy.sphere_cluster(position=positions,
                                   radius=radius,
                                   material=Ag,
                                   source=source,
                                   wavelength=600*nm,
                                   lmax=2)
        
        t_force[i] = time_function(mie.force)
        t_flux[i] = time_function(mie.cross_sections)
        t_build[i] = time_function(partial(miepy.interactions.sphere_aggregate_tmatrix, 
                              mie.position, mie.mie_scat, mie.material_data.k_b))

        A = miepy.interactions.sphere_aggregate_tmatrix(mie.position, mie.mie_scat, k=mie.material_data.k_b)
        t_solve[i] = time_function(partial(solve_linear_system, A, mie.p_src, method=miepy.solver.bicgstab))
        
        x = np.linspace(0, N*separation, 1)
        y = 2*radius*np.ones_like(x)
        z = np.zeros_like(x)

        t_source[i] = time_function(mie._solve_source_decomposition)

    fig, ax = plt.subplots()

    ax.plot(Nparticles, t_force, '-o', label='force')
    ax.plot(Nparticles, t_flux,'-o',  label='flux')
    ax.plot(Nparticles, t_build, '-o', label='build')
    ax.plot(Nparticles, t_solve, '-o', label='solve')
    ax.plot(Nparticles, t_source, '-o', label='source')

    ax.legend()
    ax.set(xlabel='number of particles', ylabel='runtime (s)')

    plt.show()

tests(10, step=1)
