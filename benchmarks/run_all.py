"""Orchestrate all benchmarks and provide comparison tools.

Usage:
    python -m benchmarks.run_all run [--preset quick|standard|thorough] [--skip-parallel] [--output FILE]
    python -m benchmarks.run_all compare FILE_A FILE_B
"""

import argparse
import datetime
import os
import sys

from .bench_utils import NumpyEncoder, load_results, save_results, system_metadata


def run_benchmarks(preset="standard", skip_parallel=False, output=None):
    """Run all benchmark suites and save consolidated results."""
    from . import bench_scaling, bench_solver, bench_translation

    print(f"=== MiePy Benchmark Suite ({preset} preset) ===\n")

    results = {
        "metadata": system_metadata(),
        "preset": preset,
    }

    print("--- Scaling Benchmarks ---")
    results["scaling"] = bench_scaling.run(preset)
    print()

    print("--- Translation Benchmarks ---")
    results["translation"] = bench_translation.run(preset)
    print()

    print("--- Solver Benchmarks ---")
    results["solver"] = bench_solver.run(preset)
    print()

    if not skip_parallel:
        print("--- Parallel Benchmarks ---")
        from . import bench_parallel
        results["parallel"] = bench_parallel.run(preset)
        print()
    else:
        print("--- Parallel Benchmarks (skipped) ---\n")

    # Save results
    if output is None:
        date_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        output = os.path.join(
            os.path.dirname(__file__), "results", f"baseline_{date_str}.json"
        )

    save_results(output, results)
    return output


def compare_results(file_a, file_b):
    """Compare two benchmark result files and print percentage changes."""
    a = load_results(file_a)
    b = load_results(file_b)

    print(f"Comparing:")
    print(f"  A (baseline): {file_a}")
    print(f"    {a.get('metadata', {}).get('timestamp', 'unknown')}")
    print(f"  B (current):  {file_b}")
    print(f"    {b.get('metadata', {}).get('timestamp', 'unknown')}")
    print()

    # Compare scaling benchmarks
    _compare_scaling(a, b, "n_sweep", "N")
    _compare_scaling(a, b, "lmax_sweep", "lmax")

    # Compare solver benchmarks
    _compare_solver(a, b)

    # Compare translation benchmarks
    _compare_translation(a, b)


def _compare_scaling(a, b, sweep_key, param_key):
    """Compare scaling benchmark results."""
    a_data = a.get("scaling", {}).get(sweep_key, {}).get("data", [])
    b_data = b.get("scaling", {}).get(sweep_key, {}).get("data", [])

    if not a_data or not b_data:
        return

    print(f"=== Scaling: {sweep_key} ===")
    print(f"{'Param':>8} {'Operation':>20} {'A (ms)':>10} {'B (ms)':>10} {'Change':>10}")
    print("-" * 62)

    a_by_param = {d[param_key]: d for d in a_data if "timings" in d}
    b_by_param = {d[param_key]: d for d in b_data if "timings" in d}

    for param in sorted(set(a_by_param) & set(b_by_param)):
        a_timings = a_by_param[param]["timings"]
        b_timings = b_by_param[param]["timings"]

        for op in ["aggregate_tmatrix", "solve_linear_system", "cross_sections", "E_field"]:
            a_t = a_timings.get(op, {}).get("time")
            b_t = b_timings.get(op, {}).get("time")

            if a_t is not None and b_t is not None and a_t > 0:
                change = (b_t - a_t) / a_t * 100
                sign = "+" if change > 0 else ""
                print(f"{param:>8} {op:>20} {a_t*1e3:>10.2f} {b_t*1e3:>10.2f} {sign}{change:>8.1f}%")

    print()


def _compare_solver(a, b):
    """Compare solver benchmark results."""
    a_data = a.get("solver", {}).get("data", [])
    b_data = b.get("solver", {}).get("data", [])

    if not a_data or not b_data:
        return

    print("=== Solver ===")
    print(f"{'N':>5} {'lmax':>5} {'sep':>5} "
          f"{'A iters':>8} {'B iters':>8} "
          f"{'A (ms)':>10} {'B (ms)':>10} {'Change':>10}")
    print("-" * 72)

    def key(d):
        return (d.get("N", 0), d.get("lmax", 0), d.get("sep_factor", 0))

    a_by_key = {key(d): d for d in a_data if "error" not in d}
    b_by_key = {key(d): d for d in b_data if "error" not in d}

    for k in sorted(set(a_by_key) & set(b_by_key)):
        ad = a_by_key[k]
        bd = b_by_key[k]
        a_t = ad.get("solve_time", 0)
        b_t = bd.get("solve_time", 0)

        if a_t > 0:
            change = (b_t - a_t) / a_t * 100
            sign = "+" if change > 0 else ""
            print(f"{k[0]:>5} {k[1]:>5} {k[2]:>5}x "
                  f"{ad.get('iterations', '?'):>8} {bd.get('iterations', '?'):>8} "
                  f"{a_t*1e3:>10.2f} {b_t*1e3:>10.2f} {sign}{change:>8.1f}%")

    print()


def _compare_translation(a, b):
    """Compare translation benchmark results."""
    a_data = a.get("translation", {}).get("single_translation", {}).get("data", [])
    b_data = b.get("translation", {}).get("single_translation", {}).get("data", [])

    if not a_data or not b_data:
        return

    print("=== Translation ===")
    print(f"{'lmax':>6} {'A (ms)':>10} {'B (ms)':>10} {'Change':>10}")
    print("-" * 40)

    a_by_lmax = {d["lmax"]: d for d in a_data if d.get("total_time") is not None}
    b_by_lmax = {d["lmax"]: d for d in b_data if d.get("total_time") is not None}

    for lmax in sorted(set(a_by_lmax) & set(b_by_lmax)):
        a_t = a_by_lmax[lmax]["total_time"]
        b_t = b_by_lmax[lmax]["total_time"]

        if a_t > 0:
            change = (b_t - a_t) / a_t * 100
            sign = "+" if change > 0 else ""
            print(f"{lmax:>6} {a_t*1e3:>10.2f} {b_t*1e3:>10.2f} {sign}{change:>8.1f}%")

    print()


def main():
    parser = argparse.ArgumentParser(description="MiePy Benchmark Suite")
    subparsers = parser.add_subparsers(dest="command")

    # run subcommand
    run_parser = subparsers.add_parser("run", help="Run benchmarks")
    run_parser.add_argument(
        "--preset", choices=["quick", "standard", "thorough"], default="standard"
    )
    run_parser.add_argument("--skip-parallel", action="store_true")
    run_parser.add_argument("--output", "-o", help="Output JSON file path")

    # compare subcommand
    cmp_parser = subparsers.add_parser("compare", help="Compare two result files")
    cmp_parser.add_argument("file_a", help="Baseline results JSON")
    cmp_parser.add_argument("file_b", help="Current results JSON")

    args = parser.parse_args()

    if args.command == "run":
        output = run_benchmarks(
            preset=args.preset,
            skip_parallel=args.skip_parallel,
            output=args.output,
        )
        print(f"\nDone! Results at: {output}")
    elif args.command == "compare":
        compare_results(args.file_a, args.file_b)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
