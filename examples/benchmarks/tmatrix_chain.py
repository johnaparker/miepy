"""performance test."""

from functools import partial

import matplotlib.pyplot as plt
import numpy as np
from timer import time_function

import miepy
from miepy.interactions import solve_linear_system

nm = 1e-9

Ag = miepy.materials.Ag()
radius = 75 * nm
source = miepy.sources.plane_wave.from_string(polarization="rhc")
separation = 165 * nm


def tmatrix_time(gmt):
    tmatrices = {}
    for i in range(gmt.Nparticles):
        gmt.position[i] = gmt.particles[i].position
        gmt.material[i] = gmt.particles[i].material

        key = gmt.particles[i]._dict_key(gmt.wavelength)
        if key in tmatrices:
            gmt.particles[i].tmatrix_fixed = tmatrices[key]
            gmt.particles[i]._rotate_fixed_tmatrix()
            gmt.tmatrix[i] = gmt.particles[i].tmatrix
        else:
            gmt.tmatrix[i] = gmt.particles[i].compute_tmatrix(gmt.lmax, gmt.wavelength, gmt.medium.eps(gmt.wavelength))
            tmatrices[key] = gmt.particles[i].tmatrix_fixed


def tests(Nmax, step=1):
    Nparticles = np.arange(1, Nmax + 1, step)
    t_force, t_flux, t_build, t_solve, t_expand, t_tmatrix = [np.zeros_like(Nparticles, dtype=float) for i in range(6)]
    for i, N in enumerate(Nparticles):
        print(N, Nmax)
        positions = [[n * separation, 0, 0] for n in range(N)]
        particles = [miepy.spheroid(pos, radius, 0.5 * radius, Ag) for pos in positions]
        mie = miepy.cluster(particles=particles, source=source, wavelength=600 * nm, lmax=2)

        t_force[i] = time_function(mie.force)
        t_flux[i] = time_function(mie.cross_sections)
        t_build[i] = time_function(
            partial(miepy.interactions.particle_aggregate_tmatrix, mie.position, mie.tmatrix, mie.material_data.k_b)
        )

        A = miepy.interactions.particle_aggregate_tmatrix(mie.position, mie.tmatrix, k=mie.material_data.k_b)
        t_solve[i] = time_function(partial(solve_linear_system, A, mie.p_src, method=miepy.solver.bicgstab))
        t_tmatrix[i] = time_function(partial(tmatrix_time, mie))

        x = np.linspace(0, N * separation, 1)
        y = 2 * radius * np.ones_like(x)
        z = np.zeros_like(x)

        t_expand[i] = time_function(partial(mie.E_field, x, y, z))

    fig, ax = plt.subplots()

    ax.plot(Nparticles, t_force, label="force")
    ax.plot(Nparticles, t_flux, label="flux")
    ax.plot(Nparticles, t_build, label="build")
    ax.plot(Nparticles, t_solve, label="solve")
    ax.plot(Nparticles, t_expand, label="expand")
    ax.plot(Nparticles, t_tmatrix, label="tmatrix")

    ax.legend()
    ax.set(xlabel="number of particles", ylabel="runtime (s)")

    plt.show()


tests(50, step=5)
