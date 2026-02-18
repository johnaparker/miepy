"""JAX implementation of VSH translation coefficients."""

import numpy as np


def vsh_translation_jax(m, n, u, v, rad, theta, phi, k, mode):
    """Compute VSH translation coefficients A and B.

    This mirrors the C++ vsh_translation() function from vsh_translation.cpp.

    Arguments:
        m, n    source mode indices
        u, v    target mode indices
        rad     radial distance
        theta   polar angle
        phi     azimuthal angle
        k       wavenumber
        mode    vsh_mode (outgoing, incident, etc.)

    Returns:
        (A_translation, B_translation) complex scalars
    """
    from math import factorial
    from miepy.backends.jax.special import (
        associated_legendre,
        gaunt_batch,
        spherical_hn,
        spherical_jn,
    )
    from miepy.cpp.vsh_functions import vsh_mode

    # The C++ code negates m at the start
    m = -m

    # Select radial function based on mode
    if mode is vsh_mode.outgoing:
        zn = spherical_hn
    elif mode is vsh_mode.ingoing:
        from miepy.backends.jax.special import spherical_hn_2
        zn = spherical_hn_2
    elif mode in (vsh_mode.incident, vsh_mode.interior):
        zn = spherical_jn
    else:
        raise ValueError(f"Unknown vsh_mode: {mode}")

    factor = 0.5 * (-1)**m * np.sqrt(
        (2*v + 1) * (2*n + 1) * factorial(v - u) * factorial(n - m)
        / (v * (v + 1) * n * (n + 1) * factorial(v + u) * factorial(n + m))
    )

    cos_theta = np.cos(theta)
    exp_phi = np.exp(1j * (u + m) * phi)

    a_vals, b_vals = gaunt_batch(m, n, u, v)

    # A translation: sum over q
    qmax_A = min(n, v, (n + v - abs(m + u)) // 2)
    sum_A = 0.0 + 0.0j
    for q in range(qmax_A + 1):
        p = n + v - 2 * q
        aq = a_vals[q]
        A_coeff = aq * (1j)**p * (n*(n+1) + v*(v+1) - p*(p+1))
        Pnm_val = associated_legendre(p, u + m, cos_theta)
        zn_val = zn(p, k * rad)
        sum_A += A_coeff * Pnm_val * zn_val

    A_translation = factor * exp_phi * sum_A

    # B translation: sum over q (starts at 1)
    qmax_B = min(n, v, (n + v + 1 - abs(m + u)) // 2)
    sum_B = 0.0 + 0.0j
    for q in range(1, qmax_B + 1):
        p = n + v - 2 * q
        bq = b_vals[q]
        B_coeff = bq * (1j)**(p+1) * np.sqrt(
            ((p+1)**2 - (n-v)**2) * ((n+v+1)**2 - (p+1)**2)
        )
        Pnm_val = associated_legendre(p + 1, u + m, cos_theta)
        zn_val = zn(p + 1, k * rad)
        sum_B += B_coeff * Pnm_val * zn_val

    B_translation = -factor * exp_phi * sum_B

    return A_translation, B_translation
