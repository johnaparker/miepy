"""Benchmark JAX GPU vs C++ CPU on hexagonal clusters.

Usage:
    python benchmarks/bench_jax_hex.py [--lmax 2] [--layers 0 1 2 3 4 5]
"""

import argparse
import numpy as np
from timeit import default_timer as timer

import miepy

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


def time_fn(func, warmup=2, repeats=5):
    for _ in range(warmup):
        func()
    times = []
    for _ in range(repeats):
        t0 = timer()
        func()
        times.append(timer() - t0)
    return np.mean(times)


def build_cluster(L, lmax, backend='cpp', wavelength=600 * nm):
    positions = hexagonal_lattice_layers(L) * wavelength
    Ag = miepy.materials.Ag()
    source = miepy.sources.plane_wave.from_string(polarization="rhc")

    with miepy.backends.backend(backend):
        cluster = miepy.sphere_cluster(
            position=positions,
            radius=75 * nm,
            material=Ag,
            source=source,
            wavelength=wavelength,
            lmax=lmax,
        )
    return cluster


def bench(layers, lmax):
    import jax
    print(f"JAX backend: {jax.default_backend()}, device: {jax.devices()[0]}")
    print(f"lmax={lmax}, layers={layers}")
    print()

    header = (
        f"{'L':>3} {'N':>4} {'dim':>6} | "
        f"{'C++ solve':>10} {'C++ force':>10} {'C++ total':>10} | "
        f"{'JAX solve':>10} {'JAX force':>10} {'JAX total':>10} | "
        f"{'solve ratio':>12} {'total ratio':>12}"
    )
    print(header)
    print("-" * len(header))

    for L in layers:
        positions = hexagonal_lattice_layers(L)
        N = len(positions)
        rmax = lmax * (lmax + 2)
        dim = 2 * rmax * N

        # --- C++ backend ---
        with miepy.backends.backend('cpp'):
            cluster_cpp = build_cluster(L, lmax, backend='cpp')

            # T-matrix + solve separately
            A_cpp = miepy.interactions.sphere_aggregate_tmatrix(
                cluster_cpp.position, cluster_cpp.mie_scat, cluster_cpp.material_data.k_b
            )
            if N > 1:
                cpp_solve = time_fn(
                    lambda: miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.exact
                    )
                )
            else:
                cpp_solve = 0.0

            # Force only (after solve)
            cluster_cpp.solve()
            cpp_force = time_fn(cluster_cpp.force)

            # Full solve + force
            def cpp_solve_and_force():
                cluster_cpp.solve()
                cluster_cpp.force()
            cpp_total = time_fn(cpp_solve_and_force)

        # --- JAX backend ---
        with miepy.backends.backend('jax'):
            # Use C++ T-matrix, only JAX solver
            if N > 1:
                # Warmup JIT
                for _ in range(3):
                    miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.exact
                    )
                jax_solve = time_fn(
                    lambda: miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.exact
                    )
                )
            else:
                jax_solve = 0.0

            # Full solve + force (includes JAX T-matrix assembly)
            cluster_jax = build_cluster(L, lmax, backend='jax')
            # Warmup
            cluster_jax.solve()
            cluster_jax.force()

            # Force only (after solve)
            cluster_jax.solve()
            jax_force = time_fn(cluster_jax.force)

            # Full solve + force
            def jax_solve_and_force():
                cluster_jax.solve()
                cluster_jax.force()
            jax_total = time_fn(jax_solve_and_force)

        solve_ratio = jax_solve / cpp_solve if cpp_solve > 0 else float('nan')
        total_ratio = jax_total / cpp_total if cpp_total > 0 else float('nan')

        print(
            f"{L:3d} {N:4d} {dim:6d} | "
            f"{cpp_solve*1e3:9.2f}ms {cpp_force*1e3:9.2f}ms {cpp_total*1e3:9.2f}ms | "
            f"{jax_solve*1e3:9.2f}ms {jax_force*1e3:9.2f}ms {jax_total*1e3:9.2f}ms | "
            f"{solve_ratio:11.2f}x {total_ratio:11.2f}x"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--lmax", type=int, default=2)
    parser.add_argument("--layers", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5])
    args = parser.parse_args()

    bench(args.layers, args.lmax)
