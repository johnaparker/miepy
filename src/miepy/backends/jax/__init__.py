"""JAX backend for MiePy.

Enables float64 precision on first import. All MiePy computations require
double precision complex128 arithmetic.
"""

try:
    import jax
    jax.config.update('jax_enable_x64', True)
    import jax.numpy as jnp  # noqa: F401
    HAS_JAX = True
except ImportError:
    HAS_JAX = False

if not HAS_JAX:
    raise ImportError(
        "JAX is required for the 'jax' backend. "
        "Install it with: pip install jax jaxlib"
    )
