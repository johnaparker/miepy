"""JAX implementations of Mie scattering coefficients for spheres.

Ports the functions from miepy.mie_single.mie_sphere to use JAX Riccati-Bessel
functions from miepy.backends.jax.special.
"""

import jax.numpy as jnp
from miepy.backends.jax.special import riccati_1_single, riccati_3_single


def mie_sphere_scattering_coefficients(radius, n, eps, mu, eps_b, mu_b, k, conducting=False):
    """Mie scattering coefficients (a_n, b_n) for a sphere.

    Arguments:
        radius      sphere radius
        n           coefficient order (n=1,2,...)
        eps         sphere permittivity
        mu          sphere permeability
        eps_b       medium permittivity
        mu_b        medium permeability
        k           medium wavenumber
        conducting  if True, calculate for conducting sphere (default: False)
    """
    xvals = k * radius
    m_ratio = (eps / eps_b) ** 0.5
    jn = riccati_1_single(n, xvals)
    yn = riccati_3_single(n, xvals)

    if conducting:
        a = jn[1] / yn[1]
        b = jn[0] / yn[0]
    else:
        jnm = riccati_1_single(n, m_ratio * xvals)
        mt = m_ratio * mu_b / mu

        a = (mt * jnm[0] * jn[1] - jn[0] * jnm[1]) / (mt * jnm[0] * yn[1] - yn[0] * jnm[1])
        b = (jnm[0] * jn[1] - mt * jn[0] * jnm[1]) / (jnm[0] * yn[1] - mt * yn[0] * jnm[1])

    return a, b


def mie_sphere_interior_coefficients(radius, n, eps, mu, eps_b, mu_b, k, conducting=False):
    """Mie interior coefficients (c_n, d_n) for a sphere.

    Arguments:
        radius      sphere radius
        n           coefficient order (n=1,2,...)
        eps         sphere permittivity
        mu          sphere permeability
        eps_b       medium permittivity
        mu_b        medium permeability
        k           medium wavenumber
        conducting  if True, calculate for conducting sphere (default: False)
    """
    xvals = k * radius
    m_ratio = (eps / eps_b) ** 0.5
    mt = m_ratio * mu_b / mu

    jn = riccati_1_single(n, xvals)
    yn = riccati_3_single(n, xvals)
    jnm = riccati_1_single(n, m_ratio * xvals)

    if conducting:
        c = jnp.zeros_like(xvals, dtype=jnp.complex128)
        d = jnp.zeros_like(xvals, dtype=jnp.complex128)
    else:
        c = (m_ratio * jn[0] * yn[1] - m_ratio * yn[0] * jn[1]) / (jnm[0] * yn[1] - mt * yn[0] * jnm[1])
        d = (m_ratio * jn[0] * yn[1] - m_ratio * yn[0] * jn[1]) / (mt * jnm[0] * yn[1] - yn[0] * jnm[1])

    return c, d
