"""Shared utilities for benchmarking."""

import datetime
import json
import os
import platform
from timeit import default_timer as timer

import numpy as np


def time_function(func, min_runtime=0.1):
    """Time a function by running it repeatedly for at least min_runtime seconds.

    Adapted from src/miepy/examples/benchmarks/timer.py.
    Returns average time per call in seconds.
    """
    t = 0
    count = 0

    while t < min_runtime:
        t0 = timer()
        func()
        tf = timer()
        t += tf - t0
        count += 1

    return t / count


def time_function_safe(func, min_runtime=0.1):
    """Like time_function but catches MemoryError and other exceptions.

    Returns (time_seconds, error_string_or_None).
    """
    try:
        result = time_function(func, min_runtime)
        return result, None
    except MemoryError:
        return None, "MemoryError"
    except Exception as e:
        return None, str(e)


def linear_chain_positions(N, separation):
    """Generate positions for N particles in a linear chain along x-axis."""
    positions = np.zeros((N, 3))
    positions[:, 0] = np.arange(N) * separation
    return positions


def system_metadata():
    """Capture system metadata for reproducibility."""
    meta = {
        "timestamp": datetime.datetime.now().isoformat(),
        "cpu": platform.processor() or platform.machine(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "not set"),
    }

    try:
        import miepy
        meta["miepy_version"] = getattr(miepy, "__version__", "unknown")
    except ImportError:
        meta["miepy_version"] = "not installed"

    return meta


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""

    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_results(filepath, results):
    """Save benchmark results as JSON with numpy type handling."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as f:
        json.dump(results, f, indent=2, cls=NumpyEncoder)
    print(f"Results saved to {filepath}")


def load_results(filepath):
    """Load benchmark results from JSON."""
    with open(filepath) as f:
        return json.load(f)
