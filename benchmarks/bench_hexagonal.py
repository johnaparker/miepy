"""Benchmark tightly-coupled hexagonal cluster solving."""

import numpy as np

import miepy

from .bench_utils import time_function_safe

nm = 1e-9


def hexagonal_lattice_layers(L):
    """Return a hexagonal lattice with unit spacing and L layers."""
    k1 = np.array([1, 0, 0], dtype=float)
    k2 = np.array([np.cos(np.pi / 3), np.sin(np.pi / 3), 0], dtype=float)

    lattice = [np.zeros(3)]

    for step in range(1, L + 1):
        lattice.append(step * k1)

        for direc in [k2 - k1, -k1, -k2, -k2 + k1, k1, k2]:
            for _ in range(step):
                lattice.append(lattice[-1] + direc)
        lattice.pop()

    return np.asarray(lattice, dtype=float)


# Layer counts that give particle counts up to ~100
# L=0: 1, L=1: 7, L=2: 19, L=3: 37, L=4: 61, L=5: 91
LAYER_COUNTS = [1, 2, 3, 4, 5]


def _build_hex_cluster(L, lmax, wavelength=600 * nm):
    """Build a tightly-coupled hexagonal cluster.

    Interparticle spacing = 1 wavelength (center-to-center).
    """
    Ag = miepy.materials.Ag()
    source = miepy.sources.plane_wave.from_string(polarization="rhc")
    positions = hexagonal_lattice_layers(L) * wavelength

    cluster = miepy.sphere_cluster(
        position=positions,
        radius=75 * nm,
        material=Ag,
        source=source,
        wavelength=wavelength,
        lmax=lmax,
    )
    return cluster


def run(preset="standard"):
    """Run hexagonal cluster benchmarks."""
    lmax_values = [2, 3]

    results = []
    for lmax in lmax_values:
        for L in LAYER_COUNTS:
            positions = hexagonal_lattice_layers(L)
            N = len(positions)
            print(f"  hex L={L} (N={N}), lmax={lmax}...", end=" ", flush=True)

            try:
                cluster = _build_hex_cluster(L, lmax)

                # Time aggregate T-matrix
                t_agg, err_agg = time_function_safe(
                    lambda: miepy.interactions.sphere_aggregate_tmatrix(
                        cluster.position, cluster.mie_scat, cluster.material_data.k_b
                    )
                )

                # Time solve (call C++ directly to work across code versions)
                A = miepy.interactions.sphere_aggregate_tmatrix(
                    cluster.position, cluster.mie_scat, k=cluster.material_data.k_b
                )
                rmax = cluster.p_src.shape[-1]
                N = len(cluster.position)
                size = N * 2 * rmax
                A_flat = A.reshape(size, size)
                p_src_flat = cluster.p_src.reshape(-1)
                t_solve, err_solve = time_function_safe(
                    lambda: miepy.cpp.interactions.solve_linear_system(
                        A_flat, p_src_flat, miepy.solver.bicgstab
                    )
                )

                entry = {
                    "L": L,
                    "N": N,
                    "lmax": lmax,
                    "aggregate_tmatrix": {"time": t_agg, "error": err_agg},
                    "solve_linear_system": {"time": t_solve, "error": err_solve},
                }
                results.append(entry)
                print(f"agg={t_agg*1e3:.1f}ms, solve={t_solve*1e3:.1f}ms")
            except Exception as e:
                print(f"FAILED: {e}")
                results.append({"L": L, "N": N, "lmax": lmax, "error": str(e)})

    return {"hexagonal": results}
