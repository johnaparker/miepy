"""Benchmark parallel scaling (strong scaling with OMP_NUM_THREADS).

Uses subprocess-based design because OMP_NUM_THREADS must be set before
the C++ library is loaded.

Orchestrator mode: spawns child processes with different thread counts.
Worker mode (--worker): runs benchmarks at current thread count.
"""

import json
import os
import subprocess
import sys
from functools import partial
from timeit import default_timer as timer

import numpy as np

PRESETS = {
    "quick": {
        "thread_counts": [1, 2],
        "configs": [{"N": 20, "lmax": 2}],
    },
    "standard": {
        "thread_counts": [1, 2, 4],
        "configs": [
            {"N": 50, "lmax": 2},
            {"N": 20, "lmax": 4},
        ],
    },
    "thorough": {
        "thread_counts": [1, 2, 4, 8],
        "configs": [
            {"N": 50, "lmax": 2},
            {"N": 20, "lmax": 4},
            {"N": 100, "lmax": 2},
            {"N": 50, "lmax": 4},
        ],
    },
}


def _worker_benchmark(configs):
    """Run benchmarks in worker mode (called in subprocess)."""
    # Import miepy here so OMP_NUM_THREADS is already set
    import miepy
    from miepy.interactions import solve_linear_system

    nm = 1e-9
    results = []

    for cfg in configs:
        N = cfg["N"]
        lmax = cfg["lmax"]
        separation = 250 * nm

        positions = np.zeros((N, 3))
        positions[:, 0] = np.arange(N) * separation

        Ag = miepy.materials.Ag()
        source = miepy.sources.plane_wave.from_string(polarization="rhc")

        cluster = miepy.sphere_cluster(
            position=positions,
            radius=75 * nm,
            material=Ag,
            source=source,
            wavelength=600 * nm,
            lmax=lmax,
        )

        timings = {}

        # aggregate T-matrix build
        def build():
            miepy.interactions.sphere_aggregate_tmatrix(
                cluster.position, cluster.mie_scat, cluster.material_data.k_b
            )

        t0 = timer()
        build()  # warmup
        t_warmup = timer() - t0

        # Time with repeated runs
        runs = max(1, int(0.2 / max(t_warmup, 1e-6)))
        t0 = timer()
        for _ in range(runs):
            build()
        timings["aggregate_tmatrix"] = (timer() - t0) / runs

        # Solve
        if N > 1:
            A = miepy.interactions.sphere_aggregate_tmatrix(
                cluster.position, cluster.mie_scat, k=cluster.material_data.k_b
            )

            def solve():
                solve_linear_system(A, cluster.p_src, method=miepy.solver.bicgstab)

            t0 = timer()
            solve()
            t_warmup = timer() - t0
            runs = max(1, int(0.2 / max(t_warmup, 1e-6)))
            t0 = timer()
            for _ in range(runs):
                solve()
            timings["solve_linear_system"] = (timer() - t0) / runs
        else:
            timings["solve_linear_system"] = 0.0

        # Cross-sections
        t0 = timer()
        cluster.cross_sections()
        t_warmup = timer() - t0
        runs = max(1, int(0.2 / max(t_warmup, 1e-6)))
        t0 = timer()
        for _ in range(runs):
            cluster.cross_sections()
        timings["cross_sections"] = (timer() - t0) / runs

        results.append({"N": N, "lmax": lmax, "timings": timings})

    return results


def run_worker():
    """Entry point for worker subprocess. Reads config from stdin, writes JSON to stdout."""
    configs = json.loads(sys.stdin.read())
    results = _worker_benchmark(configs)
    print(json.dumps(results))


def run(preset="standard"):
    """Orchestrate parallel scaling benchmarks by spawning subprocesses."""
    cfg = PRESETS[preset]
    thread_counts = cfg["thread_counts"]
    configs = cfg["configs"]

    # Cap thread counts at available CPUs
    n_cpus = os.cpu_count() or 1
    thread_counts = [t for t in thread_counts if t <= n_cpus]
    if not thread_counts:
        thread_counts = [1]

    print(f"Parallel scaling benchmark ({preset} preset)")
    print(f"  Thread counts: {thread_counts}")
    print(f"  Configs: {configs}")

    all_results = {}

    for n_threads in thread_counts:
        print(f"\n  OMP_NUM_THREADS={n_threads}...", flush=True)

        env = os.environ.copy()
        env["OMP_NUM_THREADS"] = str(n_threads)

        try:
            proc = subprocess.run(
                [sys.executable, "-m", "benchmarks.bench_parallel", "--worker"],
                input=json.dumps(configs),
                capture_output=True,
                text=True,
                env=env,
                timeout=300,
            )

            if proc.returncode != 0:
                print(f"    Worker failed: {proc.stderr[:200]}")
                all_results[str(n_threads)] = {"error": proc.stderr[:500]}
                continue

            worker_results = json.loads(proc.stdout)
            all_results[str(n_threads)] = worker_results

            for r in worker_results:
                timings_str = ", ".join(
                    f"{k}={v*1e3:.1f}ms" for k, v in r["timings"].items()
                )
                print(f"    N={r['N']}, lmax={r['lmax']}: {timings_str}")

        except subprocess.TimeoutExpired:
            print(f"    Worker timed out")
            all_results[str(n_threads)] = {"error": "timeout"}
        except Exception as e:
            print(f"    Error: {e}")
            all_results[str(n_threads)] = {"error": str(e)}

    # Compute speedups relative to single-threaded
    speedups = _compute_speedups(all_results)

    return {
        "benchmark": "parallel",
        "thread_counts": thread_counts,
        "configs": configs,
        "raw_results": all_results,
        "speedups": speedups,
    }


def _compute_speedups(all_results):
    """Compute speedup factors relative to single-threaded results."""
    if "1" not in all_results or isinstance(all_results["1"], dict) and "error" in all_results["1"]:
        return {}

    baseline = all_results["1"]
    if not isinstance(baseline, list):
        return {}

    speedups = {}
    for n_threads_str, results in all_results.items():
        if not isinstance(results, list):
            continue
        n_threads = int(n_threads_str)
        thread_speedups = []
        for i, r in enumerate(results):
            if i >= len(baseline):
                break
            base = baseline[i]
            entry = {"N": r["N"], "lmax": r["lmax"]}
            for op in r.get("timings", {}):
                base_time = base.get("timings", {}).get(op, 0)
                cur_time = r["timings"][op]
                if base_time > 0 and cur_time > 0:
                    entry[op] = base_time / cur_time
            thread_speedups.append(entry)
        speedups[n_threads_str] = thread_speedups

    return speedups


if __name__ == "__main__":
    if "--worker" in sys.argv:
        run_worker()
    else:
        from .bench_utils import NumpyEncoder

        results = run("quick")
        print(json.dumps(results, indent=2, cls=NumpyEncoder))
