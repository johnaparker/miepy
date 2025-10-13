"""performance test."""

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from timer import time_function
from topics.photonic_clusters.create_lattice import hexagonal_lattice_particles

import miepy
from miepy.interactions import solve_linear_system

nm = 1e-9

Ag = miepy.materials.Ag()
radius = 75 * nm
source = miepy.sources.plane_wave.from_string(polarization="rhc")
source = miepy.sources.gaussian_beam(2500 * nm, [1, 1j], power=0.002)

x = np.linspace(-4000 * nm, 4000 * nm, 100)
y = np.linspace(-4000 * nm, 4000 * nm, 100)
source = miepy.sources.grid_interpolate_source(source, [x, y, 0])
separation = 600 * nm

x = np.arange(-1500 * nm, 1500 * nm, 50 * nm)
y = np.arange(-1500 * nm, 1500 * nm, 50 * nm)
# source = miepy.sources.grid_interp_source(source, grid=(x,y,0))


def tests(Nmax, step=1):
    Nparticles = np.arange(1, Nmax + 1, step)
    t_force, t_flux, t_build, t_solve, t_source = [np.zeros_like(Nparticles, dtype=float) for i in range(5)]
    for i, N in enumerate(Nparticles):
        print(N, Nmax)
        # positions = [[n*separation, 0, 0] for n in range(N)]
        positions = hexagonal_lattice_particles(N) * separation
        mie = miepy.sphere_cluster(
            position=positions, radius=radius, material=Ag, source=source, wavelength=600 * nm, lmax=2
        )

        t_force[i] = time_function(mie.force)
        t_flux[i] = time_function(mie.cross_sections)
        t_build[i] = time_function(
            partial(miepy.interactions.sphere_aggregate_tmatrix, mie.position, mie.mie_scat, mie.material_data.k_b)
        )

        A = miepy.interactions.sphere_aggregate_tmatrix(mie.position, mie.mie_scat, k=mie.material_data.k_b)
        t_solve[i] = time_function(partial(solve_linear_system, A, mie.p_src, method=miepy.solver.bicgstab))

        x = np.linspace(0, N * separation, 1)
        2 * radius * np.ones_like(x)
        np.zeros_like(x)

        t_source[i] = time_function(mie._solve_source_decomposition)

    fig, ax = plt.subplots()

    ax.plot(Nparticles, t_force * 1e3, "-o", label="force")
    ax.plot(Nparticles, t_flux * 1e3, "-o", label="flux")
    ax.plot(Nparticles, t_build * 1e3, "-o", label="build")
    ax.plot(Nparticles, t_solve * 1e3, "-o", label="solve")
    ax.plot(Nparticles, t_source * 1e3, "-o", label="source")

    ax.legend()
    ax.set(xlabel="number of particles", ylabel="runtime (ms)")

    plt.show()


tests(10, step=1)
