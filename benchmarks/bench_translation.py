"""Benchmark VSH translation coefficient computation."""

from timeit import default_timer as timer

import numpy as np

import miepy
from miepy.cpp.vsh_functions import vsh_mode

from .bench_utils import time_function, time_function_safe

PRESETS = {
    "quick": {"lmax_values": [2, 4, 6]},
    "standard": {"lmax_values": [2, 3, 4, 5, 6, 8]},
    "thorough": {"lmax_values": [2, 3, 4, 5, 6, 8, 10, 12]},
}


def bench_single_translation(lmax_values=None, preset="standard"):
    """Benchmark individual VSH translation coefficient calls.

    For each lmax, iterates all (n,m,v,u) mode pairs and times the total.
    """
    if lmax_values is None:
        lmax_values = PRESETS[preset]["lmax_values"]

    print(f"Single translation benchmark: lmax={lmax_values}")
    results = []

    # Fixed geometry: two particles separated by 250nm along x
    rad = 250e-9
    theta = np.pi / 2
    phi = 0.0
    k = 2 * np.pi / 600e-9

    for lmax in lmax_values:
        rmax = lmax * (lmax + 2)
        n_calls = 0

        # Count total calls
        for n in range(1, lmax + 1):
            for m in range(-n, n + 1):
                for v in range(1, lmax + 1):
                    for u in range(-v, v + 1):
                        n_calls += 1

        def run_all_translations():
            for n in range(1, lmax + 1):
                for m in range(-n, n + 1):
                    for v in range(1, lmax + 1):
                        for u in range(-v, v + 1):
                            miepy.cpp.vsh_translation.vsh_translation(
                                m, n, u, v, rad, theta, phi, k,
                                vsh_mode.outgoing,
                            )

        t, err = time_function_safe(run_all_translations, min_runtime=0.05)

        entry = {
            "lmax": lmax,
            "rmax": rmax,
            "n_calls": n_calls,
            "total_time": t,
            "per_call_time": t / n_calls if t is not None else None,
            "error": err,
        }
        results.append(entry)
        if t is not None:
            print(f"  lmax={lmax}: {t*1e3:.2f}ms total, "
                  f"{t/n_calls*1e6:.2f}us/call ({n_calls} calls)")
        else:
            print(f"  lmax={lmax}: FAILED ({err})")

    return {"benchmark": "single_translation", "data": results}


def bench_cache_creation(lmax_values=None, preset="standard"):
    """Benchmark VSH cache map creation (used in aggregate T-matrix build)."""
    if lmax_values is None:
        lmax_values = PRESETS[preset]["lmax_values"]

    print(f"Cache creation benchmark: lmax={lmax_values}")
    results = []

    for lmax in lmax_values:
        def create_cache(lmax=lmax):
            miepy.cpp.vsh_translation.create_vsh_cache_map(lmax)

        t, err = time_function_safe(create_cache, min_runtime=0.05)

        entry = {"lmax": lmax, "time": t, "error": err}
        results.append(entry)
        if t is not None:
            print(f"  lmax={lmax}: {t*1e3:.2f}ms")
        else:
            print(f"  lmax={lmax}: FAILED ({err})")

    return {"benchmark": "cache_creation", "data": results}


def compute_scaling_exponents(translation_results):
    """Fit power-law scaling exponents from translation benchmark data.

    Returns exponent from log-log linear fit. Theoretical is O(lmax^4) for
    total translation time (rmax^2 calls, each O(lmax)).
    """
    data = translation_results["data"]
    valid = [(d["lmax"], d["total_time"]) for d in data
             if d["total_time"] is not None and d["lmax"] >= 2]

    if len(valid) < 2:
        return {"exponent": None, "error": "insufficient data"}

    lmax_arr = np.array([v[0] for v in valid], dtype=float)
    time_arr = np.array([v[1] for v in valid])

    coeffs = np.polyfit(np.log(lmax_arr), np.log(time_arr), 1)
    exponent = coeffs[0]

    print(f"  Scaling exponent: {exponent:.2f} (theoretical ~4-5)")
    return {"exponent": float(exponent), "theoretical": "O(lmax^4 to lmax^5)"}


def run(preset="standard"):
    """Run all translation benchmarks."""
    cfg = PRESETS[preset]
    lmax_values = cfg["lmax_values"]

    single = bench_single_translation(lmax_values, preset=preset)
    cache = bench_cache_creation(lmax_values, preset=preset)
    scaling = compute_scaling_exponents(single)

    return {
        "single_translation": single,
        "cache_creation": cache,
        "scaling_analysis": scaling,
    }


if __name__ == "__main__":
    import json
    from .bench_utils import NumpyEncoder

    results = run("quick")
    print(json.dumps(results, indent=2, cls=NumpyEncoder))
