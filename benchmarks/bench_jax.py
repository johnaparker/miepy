"""Benchmark JAX backend vs C++ backend for sphere_cluster operations.

Usage:
    python -m benchmarks.bench_jax [--preset quick|standard|thorough]
"""

from functools import partial

import numpy as np

import miepy
from miepy.interactions import solve_linear_system, sphere_aggregate_tmatrix

from .bench_utils import linear_chain_positions, time_function, time_function_safe

nm = 1e-9

PRESETS = {
    "quick": {
        "N_values": [2, 5, 20],
        "lmax_values": [2, 4],
    },
    "standard": {
        "N_values": [2, 5, 10, 20, 50],
        "lmax_values": [2, 3, 4, 6],
    },
    "thorough": {
        "N_values": [2, 5, 10, 20, 50, 100],
        "lmax_values": [2, 3, 4, 5, 6, 8],
    },
}


def _build_cluster(N, lmax, separation=250 * nm, backend='cpp'):
    """Build a sphere_cluster with the specified backend."""
    Ag = miepy.materials.Ag()
    source = miepy.sources.plane_wave.from_string(polarization="rhc")
    positions = linear_chain_positions(N, separation)

    with miepy.backends.backend(backend):
        cluster = miepy.sphere_cluster(
            position=positions,
            radius=75 * nm,
            material=Ag,
            source=source,
            wavelength=600 * nm,
            lmax=lmax,
        )
    return cluster


def _time_backend_operations(cluster, backend_name):
    """Time aggregate_tmatrix and solve_linear_system with a given backend."""
    results = {}
    N = len(cluster.position)

    with miepy.backends.backend(backend_name):
        # T-matrix build
        t, err = time_function_safe(
            partial(
                sphere_aggregate_tmatrix,
                cluster.position,
                cluster.mie_scat,
                cluster.material_data.k_b,
            )
        )
        results["aggregate_tmatrix"] = {"time": t, "error": err}

        # Solve linear system (exact)
        if N > 1:
            A = sphere_aggregate_tmatrix(
                cluster.position, cluster.mie_scat, k=cluster.material_data.k_b
            )
            t, err = time_function_safe(
                partial(solve_linear_system, A, cluster.p_src, method=miepy.solver.exact)
            )
            results["solve_exact"] = {"time": t, "error": err}

            # Solve linear system (BiCGSTAB)
            t, err = time_function_safe(
                partial(solve_linear_system, A, cluster.p_src, method=miepy.solver.bicgstab)
            )
            results["solve_bicgstab"] = {"time": t, "error": err}
        else:
            results["solve_exact"] = {"time": 0.0, "error": None}
            results["solve_bicgstab"] = {"time": 0.0, "error": None}

        # Cross-sections
        t, err = time_function_safe(cluster.cross_sections)
        results["cross_sections"] = {"time": t, "error": err}

    return results


def _fmt_time(t):
    """Format time in ms, or 'FAIL' if None."""
    if t is None:
        return "FAIL"
    return f"{t*1e3:.2f}"


def bench_comparison(N_values=None, lmax_values=None, preset="standard"):
    """Compare JAX vs C++ for each (N, lmax) combination."""
    if N_values is None:
        N_values = PRESETS[preset]["N_values"]
    if lmax_values is None:
        lmax_values = PRESETS[preset]["lmax_values"]

    print(f"JAX vs C++ comparison: N={N_values}, lmax={lmax_values}")
    results = []

    for lmax in lmax_values:
        for N in N_values:
            print(f"  N={N}, lmax={lmax}...", end=" ", flush=True)
            try:
                # Build cluster once (uses C++ to ensure correctness)
                cluster = _build_cluster(N, lmax, backend='cpp')

                # Time C++ backend
                cpp_timings = _time_backend_operations(cluster, 'cpp')

                # Time JAX backend
                jax_timings = _time_backend_operations(cluster, 'jax')

                entry = {
                    "N": N,
                    "lmax": lmax,
                    "cpp": cpp_timings,
                    "jax": jax_timings,
                }

                # Print summary
                ops = ["aggregate_tmatrix", "solve_exact", "solve_bicgstab"]
                parts = []
                for op in ops:
                    ct = cpp_timings.get(op, {}).get("time")
                    jt = jax_timings.get(op, {}).get("time")
                    if ct and jt and ct > 0:
                        ratio = jt / ct
                        parts.append(f"{op}={ratio:.2f}x")
                print(", ".join(parts) if parts else "done")

                results.append(entry)
            except Exception as e:
                print(f"FAILED: {e}")
                results.append({"N": N, "lmax": lmax, "error": str(e)})

    return {"comparison": results}


def run(preset="standard"):
    """Run the JAX benchmark suite."""
    return bench_comparison(preset=preset)


if __name__ == "__main__":
    import argparse
    import json
    from .bench_utils import NumpyEncoder, system_metadata

    parser = argparse.ArgumentParser(description="JAX vs C++ benchmarks")
    parser.add_argument("--preset", choices=["quick", "standard", "thorough"],
                        default="quick")
    args = parser.parse_args()

    results = {
        "metadata": system_metadata(),
        "preset": args.preset,
        "jax_comparison": run(args.preset),
    }
    print(json.dumps(results, indent=2, cls=NumpyEncoder))
