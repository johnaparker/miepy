"""performance test."""

from functools import partial

import numpy as np
from timer import time_function

import miepy
from miepy.interactions import solve_linear_system

nm = 1e-9

Ag = miepy.materials.Ag()
radius = 75 * nm
source = miepy.sources.plane_wave.from_string(polarization="rhc")
separation = 250 * nm


def tests(Nmax, step=1):
    Nparticles = np.arange(1, Nmax + 1, step)
    names = ["build", "solve", "flux", "force"]
    ftimer = {name: np.zeros_like(Nparticles, dtype=float) for name in names}
    for i, N in enumerate(Nparticles):
        print(N, Nmax)
        positions = [[n * separation, 0, 0] for n in range(N)]
        mie = miepy.sphere_cluster(
            position=positions, radius=radius, material=Ag, source=source, wavelength=600 * nm, lmax=2
        )

        ftimer["force"][i] = time_function(mie.force)
        ftimer["flux"][i] = time_function(mie.cross_sections)
        ftimer["build"][i] = time_function(
            partial(miepy.interactions.sphere_aggregate_tmatrix, mie.position, mie.mie_scat, mie.material_data.k_b)
        )

        A = miepy.interactions.sphere_aggregate_tmatrix(mie.position, mie.mie_scat, k=mie.material_data.k_b)
        ftimer["solve"][i] = time_function(partial(solve_linear_system, A, mie.p_src, method=miepy.solver.bicgstab))

    return ftimer, Nparticles


import os

import h5py

OMP = os.environ["OMP_NUM_THREADS"]
ftimer, Nparticles = tests(300, step=10)

with h5py.File(f"out_{OMP}.h5", "w") as f:
    for name, time in ftimer.items():
        f[name] = time
    f.attrs["Nparticles"] = Nparticles
