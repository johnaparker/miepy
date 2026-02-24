"""Backend registry for MiePy compute backends.

Supports 'cpu' (default, C++ via pybind11) and 'gpu' (JAX/XLA).
"""

import contextlib

_backend = 'cpu'


def get_backend():
    """Return the current active backend name ('cpu' or 'gpu')."""
    return _backend


def set_backend(name):
    """Set the global compute backend.

    Arguments:
        name: 'cpu' or 'gpu'
    """
    global _backend
    if name not in ('cpu', 'gpu'):
        raise ValueError(f"Unknown backend '{name}'. Choose 'cpu' or 'gpu'.")
    if name == 'gpu':
        # Trigger the lazy import guard / x64 init
        from miepy.backends import jax as _jax_backend  # noqa: F401
    _backend = name


@contextlib.contextmanager
def backend(name):
    """Context manager to temporarily switch the compute backend.

    Usage:
        with miepy.backends.backend('gpu'):
            cluster.solve()
    """
    global _backend
    old = _backend
    set_backend(name)
    try:
        yield
    finally:
        _backend = old
