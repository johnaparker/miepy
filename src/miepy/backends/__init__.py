"""Backend registry for MiePy compute backends.

Supports 'cpp' (default, C++ via pybind11) and 'jax' (JAX/XLA).
"""

import contextlib

_backend = 'cpp'


def get_backend():
    """Return the current active backend name ('cpp' or 'jax')."""
    return _backend


def set_backend(name):
    """Set the global compute backend.

    Arguments:
        name: 'cpp' or 'jax'
    """
    global _backend
    if name not in ('cpp', 'jax'):
        raise ValueError(f"Unknown backend '{name}'. Choose 'cpp' or 'jax'.")
    if name == 'jax':
        # Trigger the lazy import guard / x64 init
        from miepy.backends import jax as _jax_backend  # noqa: F401
    _backend = name


@contextlib.contextmanager
def backend(name):
    """Context manager to temporarily switch the compute backend.

    Usage:
        with miepy.backends.backend('jax'):
            cluster.solve()
    """
    global _backend
    old = _backend
    set_backend(name)
    try:
        yield
    finally:
        _backend = old
