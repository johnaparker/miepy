"""Benchmark N-particle and lmax scaling for sphere_cluster."""

from functools import partial

import numpy as np

import miepy
from miepy.interactions import solve_linear_system

from .bench_utils import linear_chain_positions, time_function, time_function_safe

nm = 1e-9
BYTES_PER_COMPLEX128 = 16  # numpy complex128

PRESETS = {
    "quick": {
        "N_values": [1, 5, 20, 50],
        "lmax_values": [2, 4],
    },
    "standard": {
        "N_values": [1, 2, 5, 10, 20, 50, 100, 200],
        "lmax_values": [2, 3, 4, 5, 6, 8],
    },
    "thorough": {
        "N_values": [1, 2, 5, 10, 20, 50, 100, 200, 300, 500],
        "lmax_values": [2, 3, 4, 5, 6, 8, 10],
    },
}


def _build_cluster(N, lmax, separation=250 * nm):
    """Build a sphere_cluster for benchmarking."""
    Ag = miepy.materials.Ag()
    source = miepy.sources.plane_wave.from_string(polarization="rhc")
    positions = linear_chain_positions(N, separation)

    cluster = miepy.sphere_cluster(
        position=positions,
        radius=75 * nm,
        material=Ag,
        source=source,
        wavelength=600 * nm,
        lmax=lmax,
    )
    return cluster


def _memory_info(N, lmax):
    """Compute theoretical memory usage for the aggregate T-matrix.

    The matrix is (2 * rmax * N)^2 complex128 entries.
    """
    rmax = lmax * (lmax + 2)
    dim = 2 * rmax * N
    matrix_bytes = dim * dim * BYTES_PER_COMPLEX128
    return {
        "matrix_dim": dim,
        "matrix_bytes": matrix_bytes,
        "matrix_mb": matrix_bytes / (1024 * 1024),
    }


def _measure_rss_delta(func):
    """Measure RSS increase during func() using resource.getrusage.

    Returns RSS delta in bytes, or None if resource module is unavailable.
    Note: this captures all process memory (C++ and Python), but is noisy
    due to fragmentation and other concurrent allocations.
    """
    try:
        import resource
        import platform as _platform
    except ImportError:
        return None

    # ru_maxrss is in bytes on macOS, kilobytes on Linux
    scale = 1 if _platform.system() == "Darwin" else 1024

    before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale
    func()
    after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * scale

    delta = after - before
    return delta if delta > 0 else None


def _time_operations(cluster, label=""):
    """Time the main operations on a cluster. Returns dict of timings."""
    results = {}
    N = len(cluster.position)

    # T-matrix build
    t, err = time_function_safe(
        partial(
            miepy.interactions.sphere_aggregate_tmatrix,
            cluster.position,
            cluster.mie_scat,
            cluster.material_data.k_b,
        )
    )
    results["aggregate_tmatrix"] = {"time": t, "error": err}

    # Solve linear system
    if N > 1:
        A = miepy.interactions.sphere_aggregate_tmatrix(
            cluster.position, cluster.mie_scat, k=cluster.material_data.k_b
        )
        t, err = time_function_safe(
            partial(solve_linear_system, A, cluster.p_src, method=miepy.solver.bicgstab)
        )
        results["solve_linear_system"] = {"time": t, "error": err}
    else:
        results["solve_linear_system"] = {"time": 0.0, "error": None}

    # Cross-sections
    t, err = time_function_safe(cluster.cross_sections)
    results["cross_sections"] = {"time": t, "error": err}

    # E-field at a single point
    def eval_field():
        cluster.E_field(0, 0, 5000 * nm)

    t, err = time_function_safe(eval_field)
    results["E_field"] = {"time": t, "error": err}

    return results


def bench_n_sweep(N_values=None, lmax=2, preset="standard"):
    """Benchmark scaling with number of particles at fixed lmax."""
    if N_values is None:
        N_values = PRESETS[preset]["N_values"]

    print(f"N-sweep benchmark: N={N_values}, lmax={lmax}")
    results = []

    for N in N_values:
        print(f"  N={N}...", end=" ", flush=True)
        try:
            cluster = _build_cluster(N, lmax)
            timings = _time_operations(cluster)
            mem = _memory_info(N, lmax)
            entry = {"N": N, "lmax": lmax, "timings": timings, "memory": mem}

            # Measure RSS delta during T-matrix build (captures C++ allocations)
            rss_delta = _measure_rss_delta(
                partial(
                    miepy.interactions.sphere_aggregate_tmatrix,
                    cluster.position,
                    cluster.mie_scat,
                    cluster.material_data.k_b,
                )
            )
            if rss_delta is not None:
                entry["memory"]["rss_delta_bytes"] = rss_delta

            results.append(entry)
            total = sum(
                v["time"] for v in timings.values() if v["time"] is not None
            )
            print(f"total={total*1e3:.1f}ms, matrix={mem['matrix_mb']:.1f}MB")
        except Exception as e:
            print(f"FAILED: {e}")
            results.append({"N": N, "lmax": lmax, "error": str(e)})

    return {"sweep": "N", "fixed_lmax": lmax, "data": results}


def bench_lmax_sweep(lmax_values=None, N=20, preset="standard"):
    """Benchmark scaling with lmax at fixed N."""
    if lmax_values is None:
        lmax_values = PRESETS[preset]["lmax_values"]

    print(f"lmax-sweep benchmark: lmax={lmax_values}, N={N}")
    results = []

    for lmax in lmax_values:
        print(f"  lmax={lmax}...", end=" ", flush=True)
        try:
            cluster = _build_cluster(N, lmax)
            timings = _time_operations(cluster)
            mem = _memory_info(N, lmax)
            entry = {"N": N, "lmax": lmax, "timings": timings, "memory": mem}

            rss_delta = _measure_rss_delta(
                partial(
                    miepy.interactions.sphere_aggregate_tmatrix,
                    cluster.position,
                    cluster.mie_scat,
                    cluster.material_data.k_b,
                )
            )
            if rss_delta is not None:
                entry["memory"]["rss_delta_bytes"] = rss_delta

            results.append(entry)
            total = sum(
                v["time"] for v in timings.values() if v["time"] is not None
            )
            print(f"total={total*1e3:.1f}ms, matrix={mem['matrix_mb']:.1f}MB")
        except Exception as e:
            print(f"FAILED: {e}")
            results.append({"N": N, "lmax": lmax, "error": str(e)})

    return {"sweep": "lmax", "fixed_N": N, "data": results}


def run(preset="standard"):
    """Run all scaling benchmarks."""
    cfg = PRESETS[preset]
    return {
        "n_sweep": bench_n_sweep(cfg["N_values"], lmax=2, preset=preset),
        "lmax_sweep": bench_lmax_sweep(cfg["lmax_values"], N=20, preset=preset),
    }


if __name__ == "__main__":
    import json
    from .bench_utils import NumpyEncoder

    results = run("quick")
    print(json.dumps(results, indent=2, cls=NumpyEncoder))
