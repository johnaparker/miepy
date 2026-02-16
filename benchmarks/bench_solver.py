"""Benchmark BiCGSTAB solver convergence and performance."""

from functools import partial
from timeit import default_timer as timer

import numpy as np

import miepy

from .bench_utils import linear_chain_positions, time_function_safe

nm = 1e-9

PRESETS = {
    "quick": {
        "configs": [
            {"N": 10, "lmax": 2, "sep_factors": [3, 10]},
            {"N": 20, "lmax": 2, "sep_factors": [3, 10]},
        ]
    },
    "standard": {
        "configs": [
            {"N": 10, "lmax": 2, "sep_factors": [2, 3, 5, 10]},
            {"N": 20, "lmax": 2, "sep_factors": [2, 3, 5, 10]},
            {"N": 20, "lmax": 4, "sep_factors": [2, 3, 5, 10]},
            {"N": 50, "lmax": 2, "sep_factors": [3, 5, 10]},
        ]
    },
    "thorough": {
        "configs": [
            {"N": 10, "lmax": 2, "sep_factors": [2, 3, 5, 10]},
            {"N": 20, "lmax": 2, "sep_factors": [2, 3, 5, 10]},
            {"N": 20, "lmax": 4, "sep_factors": [2, 3, 5, 10]},
            {"N": 50, "lmax": 2, "sep_factors": [2, 3, 5, 10]},
            {"N": 50, "lmax": 4, "sep_factors": [3, 5, 10]},
            {"N": 100, "lmax": 2, "sep_factors": [3, 5, 10]},
        ]
    },
}


def _profile_solver(N, lmax, sep_factor, radius=75 * nm):
    """Profile the BiCGSTAB solver for a given configuration.

    sep_factor: separation as multiple of diameter (2*radius).
    """
    separation = sep_factor * 2 * radius
    Ag = miepy.materials.Ag()
    source = miepy.sources.plane_wave.from_string(polarization="rhc")
    positions = linear_chain_positions(N, separation)

    cluster = miepy.sphere_cluster(
        position=positions,
        radius=radius,
        material=Ag,
        source=source,
        wavelength=600 * nm,
        lmax=lmax,
    )

    # Build aggregate T-matrix at C++ level (flat 2D matrix)
    Nparticles = len(cluster.position)
    mie_flat = cluster.mie_scat.reshape([Nparticles, -1])
    agg_tmatrix_flat = miepy.cpp.interactions.sphere_aggregate_tmatrix(
        cluster.position, mie_flat, cluster.material_data.k_b
    )

    # Build interaction matrix (add identity, same as C++ solve_linear_system)
    interaction_matrix = agg_tmatrix_flat.copy()
    for i in range(interaction_matrix.shape[0]):
        interaction_matrix[i, i] += 1

    # Flatten p_src for the solver
    p_src_flat = cluster.p_src.flatten()

    # Profile with bicgstab_profiled
    t0 = timer()
    solution, iterations, residual = miepy.cpp.interactions.bicgstab_profiled(
        interaction_matrix, p_src_flat
    )
    solve_time = timer() - t0

    # Matrix properties
    matrix_size = interaction_matrix.shape[0]
    frobenius_norm = float(np.linalg.norm(interaction_matrix, "fro"))

    return {
        "N": N,
        "lmax": lmax,
        "sep_factor": sep_factor,
        "separation_nm": separation / nm,
        "matrix_size": matrix_size,
        "iterations": iterations,
        "residual": residual,
        "solve_time": solve_time,
        "time_per_iteration": solve_time / max(iterations, 1),
        "frobenius_norm": frobenius_norm,
    }


def run(preset="standard"):
    """Run solver convergence benchmarks."""
    configs = PRESETS[preset]["configs"]

    print(f"Solver benchmark ({preset} preset)")
    results = []

    for cfg in configs:
        N = cfg["N"]
        lmax = cfg["lmax"]
        for sep_factor in cfg["sep_factors"]:
            print(f"  N={N}, lmax={lmax}, sep={sep_factor}x...", end=" ", flush=True)
            try:
                entry = _profile_solver(N, lmax, sep_factor)
                results.append(entry)
                print(f"iters={entry['iterations']}, "
                      f"residual={entry['residual']:.2e}, "
                      f"time={entry['solve_time']*1e3:.1f}ms")
            except Exception as e:
                print(f"FAILED: {e}")
                results.append({
                    "N": N, "lmax": lmax, "sep_factor": sep_factor,
                    "error": str(e),
                })

    return {"benchmark": "solver", "data": results}


if __name__ == "__main__":
    import json
    from .bench_utils import NumpyEncoder

    results = run("quick")
    print(json.dumps(results, indent=2, cls=NumpyEncoder))
