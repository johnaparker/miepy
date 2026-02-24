"""Benchmark JAX GPU vs C++ CPU on hexagonal clusters.

Times each pipeline stage separately:
  tmat   — aggregate T-matrix assembly
  src    — source decomposition
  solve  — linear system solve
  force  — force computation
  total  — full solve() + force()

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


def build_cluster(L, lmax, backend='cpu', wavelength=600 * nm):
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
        f"{'C++ tmat':>9} {'C++ src':>8} {'C++ solve':>10} {'C++ force':>10} {'C++ total':>10} | "
        f"{'JAX tmat':>9} {'JAX LU':>9} {'JAX iter':>9} {'JAX force':>10} {'JAX total':>10} | "
        f"{'LU/C++':>7} {'iter/C++':>8}"
    )
    print(header)
    print("-" * len(header))

    for L in layers:
        positions = hexagonal_lattice_layers(L)
        N = len(positions)
        rmax = lmax * (lmax + 2)
        dim = 2 * rmax * N

        # --- C++ backend ---
        with miepy.backends.backend('cpu'):
            cluster_cpp = build_cluster(L, lmax, backend='cpu')

            # T-matrix assembly
            cpp_tmat = time_fn(
                lambda: miepy.interactions.sphere_aggregate_tmatrix(
                    cluster_cpp.position, cluster_cpp.mie_scat, cluster_cpp.material_data.k_b
                )
            )

            # Source decomposition
            source = cluster_cpp.source
            cpp_src = time_fn(
                lambda: source.structure(cluster_cpp.position, cluster_cpp.material_data.k_b, cluster_cpp.lmax)
            )

            # Solve only (excludes T-matrix assembly)
            A_cpp = miepy.interactions.sphere_aggregate_tmatrix(
                cluster_cpp.position, cluster_cpp.mie_scat, cluster_cpp.material_data.k_b
            )
            if N > 1:
                cpp_solve = time_fn(
                    lambda: miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.bicgstab
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
        with miepy.backends.backend('gpu'):
            # JAX T-matrix assembly
            if N > 1:
                # Warmup
                for _ in range(3):
                    miepy.interactions.sphere_aggregate_tmatrix(
                        cluster_cpp.position, cluster_cpp.mie_scat, cluster_cpp.material_data.k_b
                    )
                jax_tmat = time_fn(
                    lambda: miepy.interactions.sphere_aggregate_tmatrix(
                        cluster_cpp.position, cluster_cpp.mie_scat, cluster_cpp.material_data.k_b
                    )
                )
            else:
                jax_tmat = 0.0

            # JAX solve — LU (using C++ T-matrix for fair comparison)
            if N > 1:
                for _ in range(3):
                    miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.exact
                    )
                jax_solve_lu = time_fn(
                    lambda: miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.exact
                    )
                )
            else:
                jax_solve_lu = 0.0

            # JAX solve — BiCGSTAB
            if N > 1:
                for _ in range(3):
                    miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.bicgstab
                    )
                jax_solve_iter = time_fn(
                    lambda: miepy.interactions.solve_linear_system(
                        A_cpp, cluster_cpp.p_src, method=miepy.solver.bicgstab
                    )
                )
            else:
                jax_solve_iter = 0.0

            # Full solve + force (includes JAX T-matrix assembly)
            cluster_jax = build_cluster(L, lmax, backend='gpu')
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

        lu_ratio = jax_solve_lu / cpp_solve if cpp_solve > 0 else float('nan')
        iter_ratio = jax_solve_iter / cpp_solve if cpp_solve > 0 else float('nan')

        print(
            f"{L:3d} {N:4d} {dim:6d} | "
            f"{cpp_tmat*1e3:8.2f}ms {cpp_src*1e3:7.2f}ms {cpp_solve*1e3:9.2f}ms {cpp_force*1e3:9.2f}ms {cpp_total*1e3:9.2f}ms | "
            f"{jax_tmat*1e3:8.2f}ms {jax_solve_lu*1e3:8.2f}ms {jax_solve_iter*1e3:8.2f}ms {jax_force*1e3:9.2f}ms {jax_total*1e3:9.2f}ms | "
            f"{lu_ratio:6.2f}x {iter_ratio:7.2f}x"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--lmax", type=int, default=2)
    parser.add_argument("--layers", type=int, nargs="+", default=[0, 1, 2, 3, 4, 5])
    args = parser.parse_args()

    bench(args.layers, args.lmax)
