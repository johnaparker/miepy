"""JAX implementations of special functions for MiePy.

Spherical Bessel/Hankel, Riccati-Bessel, associated Legendre,
angular pi/tau functions, Emn normalization, and Wigner 3j / Gaunt coefficients.

All array-valued functions use JAX arrays (complex128 where needed).
Wigner 3j and Gaunt coefficients use plain NumPy (integer-only precomputation).
"""

import numpy as np
import jax.numpy as jnp
from math import factorial as _factorial


# ---------------------------------------------------------------------------
# Spherical Bessel functions
# ---------------------------------------------------------------------------

def spherical_jn(n, z, derivative=False):
    """Spherical Bessel function of the first kind j_n(z).

    Forward recursion for complex z. n is a non-negative integer scalar;
    z may be an array.
    """
    z = jnp.asarray(z, dtype=jnp.complex128)
    sin_z = jnp.sin(z)
    cos_z = jnp.cos(z)

    jn_prev2 = sin_z / z          # j_0
    if n == 0 and not derivative:
        return jn_prev2

    jn_prev1 = jn_prev2 / z - cos_z / z  # j_1
    if n == 1 and not derivative:
        return jn_prev1

    if n >= 2:
        jn_cur = jn_prev1
        for i in range(2, n + 1):
            jn_cur = (2 * i - 1) / z * jn_prev1 - jn_prev2
            jn_prev2 = jn_prev1
            jn_prev1 = jn_cur

    if not derivative:
        return jn_prev1

    # Derivative: j_n'(z) = j_{n-1}(z) - (n+1)/z * j_n(z)
    jn_val = jn_prev1
    jn_minus1 = spherical_jn(n - 1, z)
    return jn_minus1 - (n + 1) / z * jn_val


def spherical_yn(n, z, derivative=False):
    """Spherical Bessel function of the second kind y_n(z) (Neumann).

    Forward recursion for real z. n is a non-negative integer scalar;
    z may be an array.
    """
    z = jnp.asarray(z, dtype=jnp.float64)
    sin_z = jnp.sin(z)
    cos_z = jnp.cos(z)

    yn_prev2 = -cos_z / z          # y_0
    if n == 0 and not derivative:
        return yn_prev2

    yn_prev1 = yn_prev2 / z - sin_z / z  # y_1  (= y_0/z - j_0)
    if n == 1 and not derivative:
        return yn_prev1

    if n >= 2:
        for i in range(2, n + 1):
            yn_cur = (2 * i - 1) / z * yn_prev1 - yn_prev2
            yn_prev2 = yn_prev1
            yn_prev1 = yn_cur

    if not derivative:
        return yn_prev1

    # Derivative: y_n'(z) = y_{n-1}(z) - (n+1)/z * y_n(z)
    yn_val = yn_prev1
    yn_minus1 = spherical_yn(n - 1, z)
    return yn_minus1 - (n + 1) / z * yn_val


def spherical_hn(n, z, derivative=False):
    """Spherical Hankel function of the first kind h_n^(1)(z) = j_n + i*y_n.

    For real z, uses simultaneous forward recursion of j and y.
    """
    z = jnp.asarray(z, dtype=jnp.float64)

    if not derivative:
        sin_z = jnp.sin(z)
        cos_z = jnp.cos(z)

        jn_prev2 = sin_z / z
        yn_prev2 = -cos_z / z

        if n == 0:
            return jn_prev2 + 1j * yn_prev2

        jn_prev1 = jn_prev2 / z - cos_z / z
        yn_prev1 = yn_prev2 / z - sin_z / z

        if n == 1:
            return jn_prev1 + 1j * yn_prev1

        for i in range(2, n + 1):
            jn_cur = (2 * i - 1) / z * jn_prev1 - jn_prev2
            jn_prev2 = jn_prev1
            jn_prev1 = jn_cur

            yn_cur = (2 * i - 1) / z * yn_prev1 - yn_prev2
            yn_prev2 = yn_prev1
            yn_prev1 = yn_cur

        return jn_prev1 + 1j * yn_prev1
    else:
        return spherical_hn(n - 1, z) - (n + 1) / z * spherical_hn(n, z)


def spherical_hn_2(n, z, derivative=False):
    """Spherical Hankel function of the second kind h_n^(2) = conj(h_n^(1))."""
    return jnp.conj(spherical_hn(n, z, derivative))


def spherical_hn_recursion(nmax, z):
    """All spherical Hankel h_n^(1)(z) for n = 0..nmax. Real z only.

    Returns complex array of shape [nmax+1, ...].
    """
    z = jnp.asarray(z, dtype=jnp.float64)
    sin_z = jnp.sin(z)
    cos_z = jnp.cos(z)

    jn = [None] * (nmax + 1)
    yn = [None] * (nmax + 1)

    jn[0] = sin_z / z
    yn[0] = -cos_z / z

    if nmax >= 1:
        jn[1] = jn[0] / z + yn[0]   # j_0/z - cos/z  but yn[0] = -cos/z so j_0/z + yn[0]
        yn[1] = yn[0] / z - jn[0]

    for i in range(2, nmax + 1):
        jn[i] = (2 * i - 1) / z * jn[i - 1] - jn[i - 2]
        yn[i] = (2 * i - 1) / z * yn[i - 1] - yn[i - 2]

    jn_arr = jnp.stack(jn, axis=0)
    yn_arr = jnp.stack(yn, axis=0)
    return jn_arr + 1j * yn_arr


# ---------------------------------------------------------------------------
# Riccati-Bessel functions
# ---------------------------------------------------------------------------

def riccati_1(nmax, x):
    r"""Riccati-Bessel function of the 1st kind.

    Returns shape [2, nmax, ...] where:
      result[0, n] = x * j_{n+1}(x)    (psi)
      result[1, n] = j_{n+1}(x) + x * j_{n+1}'(x)   (psi')
    for n = 0, 1, ..., nmax-1 representing orders 1..nmax.
    """
    x = jnp.asarray(x, dtype=jnp.complex128)
    vals = []
    derivs = []

    for n in range(nmax):
        jn_val = spherical_jn(n + 1, x)
        jn_der = spherical_jn(n + 1, x, derivative=True)
        vals.append(x * jn_val)
        derivs.append(jn_val + x * jn_der)

    vals = jnp.stack(vals, axis=0)
    derivs = jnp.stack(derivs, axis=0)
    return jnp.stack([vals, derivs], axis=0)


def riccati_3(nmax, x):
    r"""Riccati-Bessel function of the 3rd kind (using y_n).

    Returns shape [2, nmax, ...] where:
      result[0, n] = x * y_{n+1}(x)
      result[1, n] = y_{n+1}(x) + x * y_{n+1}'(x)
    """
    x = jnp.asarray(x, dtype=jnp.float64)
    vals = []
    derivs = []

    for n in range(nmax):
        yn_val = spherical_yn(n + 1, x)
        yn_der = spherical_yn(n + 1, x, derivative=True)
        vals.append(x * yn_val)
        derivs.append(yn_val + x * yn_der)

    vals = jnp.stack(vals, axis=0)
    derivs = jnp.stack(derivs, axis=0)
    return jnp.stack([vals, derivs], axis=0)


def riccati_2(nmax, x):
    r"""Riccati-Bessel function of the 2nd kind = riccati_1 + 1j * riccati_3."""
    return riccati_1(nmax, x) + 1j * riccati_3(nmax, x)


def riccati_1_single(n, x):
    r"""Riccati-Bessel function of the 1st kind for a single order n.

    Returns shape [2, ...]:
      result[0] = x * j_n(x)
      result[1] = j_n(x) + x * j_n'(x)
    """
    x = jnp.asarray(x, dtype=jnp.complex128)
    jn_val = spherical_jn(n, x)
    jn_der = spherical_jn(n, x, derivative=True)
    return jnp.stack([x * jn_val, jn_val + x * jn_der], axis=0)


def riccati_2_single(n, x):
    r"""Riccati-Bessel using y_n for a single order n.

    Returns shape [2, ...]:
      result[0] = x * y_n(x)
      result[1] = y_n(x) + x * y_n'(x)
    """
    x = jnp.asarray(x, dtype=jnp.float64)
    yn_val = spherical_yn(n, x)
    yn_der = spherical_yn(n, x, derivative=True)
    return jnp.stack([x * yn_val, yn_val + x * yn_der], axis=0)


def riccati_3_single(n, x):
    r"""Riccati-Bessel function of the 3rd kind for a single order.

    riccati_3 = riccati_1 + 1j * riccati_2_single  (matching the NumPy convention).
    """
    return riccati_1_single(n, x) + 1j * riccati_2_single(n, x)


# ---------------------------------------------------------------------------
# Associated Legendre functions
# ---------------------------------------------------------------------------

def associated_legendre(n, m, z):
    r"""Unnormalized associated Legendre polynomial P_n^m(z).

    Uses three-term recurrence. z = cos(theta), real.
    Matches the GSL_SF_LEGENDRE_NONE convention used in the C++ backend.
    """
    z = jnp.asarray(z, dtype=jnp.float64)
    abs_m = abs(m)

    if abs_m > n:
        return jnp.zeros_like(z)

    # Compute P_n^{|m|}(z) via recursion
    # GSL_SF_LEGENDRE_NONE convention: no Condon-Shortley phase
    # P_{|m|}^{|m|} = (2|m|-1)!! * (1-z^2)^{|m|/2}
    sin_theta_sq = 1.0 - z * z
    sin_theta_abs = jnp.sqrt(jnp.clip(sin_theta_sq, 0.0))

    # P_m^m (without Condon-Shortley phase, matching GSL)
    pmm = 1.0
    for i in range(1, abs_m + 1):
        pmm *= (2 * i - 1)
    pmm = pmm * sin_theta_abs ** abs_m

    if n == abs_m:
        leg = pmm
    else:
        # P_{m+1}^m = z * (2m+1) * P_m^m
        pmm1 = z * (2 * abs_m + 1) * pmm
        if n == abs_m + 1:
            leg = pmm1
        else:
            # (n-m)*P_n^m = (2n-1)*z*P_{n-1}^m - (n+m-1)*P_{n-2}^m
            p_prev2 = pmm
            p_prev1 = pmm1
            for nn in range(abs_m + 2, n + 1):
                p_cur = ((2 * nn - 1) * z * p_prev1 - (nn + abs_m - 1) * p_prev2) / (nn - abs_m)
                p_prev2 = p_prev1
                p_prev1 = p_cur
            leg = p_prev1

    if m < 0:
        # P_n^{-|m|} = (-1)^m * (n+m)!/(n-m)! * P_n^{|m|}
        # where m is negative, so (-1)^m = (-1)^{-|m|}
        factor = (-1) ** m * _factorial(n + m) / _factorial(n - m)
        return factor * leg
    else:
        return leg


def associated_legendre_recursion(nmax, z):
    r"""All associated Legendre polynomials P_n^m(z) for 0 <= n <= nmax.

    Returns array indexed by n*(n+2) - n + m for all (n, m) with -n <= m <= n.
    Total size: nmax*(nmax+2) + 1.
    """
    z = jnp.asarray(z, dtype=jnp.float64)
    total_size = nmax * (nmax + 2) + 1

    # Build all P_n^m(z) for m >= 0 via GSL-style recursion
    # First compute P_n^m for m >= 0 using the standard recursion
    sin_theta_sq = 1.0 - z * z
    sin_theta_abs = jnp.sqrt(jnp.clip(sin_theta_sq, 0.0))

    # Store P_n^m for all (n, m >= 0)
    # Use dict keyed by (n, m) for clarity; pack into output array at end
    pnm = {}
    pnm[(0, 0)] = jnp.ones_like(z)

    # Sectoral: P_m^m = (2m-1) * sqrt(1-z^2) * P_{m-1}^{m-1}  (no Condon-Shortley phase)
    for m in range(1, nmax + 1):
        pnm[(m, m)] = (2 * m - 1) * sin_theta_abs * pnm[(m - 1, m - 1)]

    # Semi-sectoral: P_{m+1}^m = z*(2m+1)*P_m^m
    for m in range(0, nmax):
        pnm[(m + 1, m)] = z * (2 * m + 1) * pnm[(m, m)]

    # General recursion: (n-m)*P_n^m = (2n-1)*z*P_{n-1}^m - (n+m-1)*P_{n-2}^m
    for m in range(0, nmax + 1):
        for nn in range(m + 2, nmax + 1):
            pnm[(nn, m)] = ((2 * nn - 1) * z * pnm[(nn - 1, m)] -
                            (nn + m - 1) * pnm[(nn - 2, m)]) / (nn - m)

    # Pack into output array
    parts = []
    indices = []
    for nn in range(0, nmax + 1):
        for m in range(-nn, nn + 1):
            idx = nn * (nn + 2) - nn + m
            abs_m = abs(m)
            leg = pnm[(nn, abs_m)]
            if m < 0:
                factor = (-1) ** m * _factorial(nn + m) / _factorial(nn - m)
                leg = factor * leg
            parts.append(leg)
            indices.append(idx)

    # Build result array
    # For scalar z, result is 1D of length total_size
    # For array z, result shape is [total_size, ...z.shape]
    result_list = [None] * total_size
    for i, idx in enumerate(indices):
        result_list[idx] = parts[i]

    return jnp.stack(result_list, axis=0)


# ---------------------------------------------------------------------------
# Angular pi and tau functions
# ---------------------------------------------------------------------------

def pi_func(n, m, theta):
    r"""Angular function pi_{n,m}(theta) = m/sin(theta) * P_n^m(cos(theta)).

    Special cases at theta=0 and theta=pi.
    """
    theta = jnp.asarray(theta, dtype=jnp.float64)

    # Handle scalar theta for edge cases
    at_zero = jnp.abs(theta) < 1e-30
    at_pi = jnp.abs(theta - jnp.pi) < 1e-30
    regular = ~at_zero & ~at_pi

    result = jnp.zeros_like(theta)

    # theta = 0 cases
    if m == 1:
        result = jnp.where(at_zero, n * (n + 1) / 2.0, result)
    elif m == -1:
        result = jnp.where(at_zero, 0.5, result)
    # else: already 0

    # theta = pi cases
    if m == 1:
        result = jnp.where(at_pi, (-1) ** (n + 1) * n * (n + 1) / 2.0, result)
    elif m == -1:
        result = jnp.where(at_pi, (-1) ** (n + 1) / 2.0, result)
    # else: already 0

    # General case
    z = jnp.cos(theta)
    sin_theta = jnp.sin(theta)
    # Avoid division by zero; use safe_sin that won't produce NaN
    safe_sin = jnp.where(regular, sin_theta, 1.0)
    leg = associated_legendre(n, m, z)
    general = m / safe_sin * leg
    result = jnp.where(regular, general, result)

    return result


def tau_func(n, m, theta):
    r"""Angular function tau_{n,m}(theta) = d/dtheta P_n^m(cos(theta)).

    Computed as -sin(theta) * dP_n^m/d(cos theta).
    Special cases at theta=0 and theta=pi.
    """
    theta = jnp.asarray(theta, dtype=jnp.float64)

    at_zero = jnp.abs(theta) < 1e-30
    at_pi = jnp.abs(theta - jnp.pi) < 1e-30
    regular = ~at_zero & ~at_pi

    result = jnp.zeros_like(theta)

    # theta = 0 cases
    if m == 1:
        result = jnp.where(at_zero, n * (n + 1) / 2.0, result)
    elif m == -1:
        result = jnp.where(at_zero, -0.5, result)

    # theta = pi cases
    if m == 1:
        result = jnp.where(at_pi, (-1) ** n * n * (n + 1) / 2.0, result)
    elif m == -1:
        result = jnp.where(at_pi, -(-1) ** n / 2.0, result)

    # General case: use finite difference on associated_legendre
    # tau = d/dtheta P_n^m(cos(theta)) = -sin(theta) * dP_n^m/dz
    # Compute dP_n^m/dz via the identity:
    # (1 - z^2) dP_n^m/dz = (n+m)(n-m+1) P_n^{m-1} - m*z*P_n^m   (for |m| <= n)
    # But more robustly: -sin(theta) * dP/dz = [n*z*P_n^m - (n+m)*P_{n-1}^m] / (z^2-1) * -(1-z^2)
    # Actually the standard relation:
    # -sin(theta) dP_n^m/dz = [ n*cos(theta)*P_n^m(cos theta) - (n+m)*P_{n-1}^m(cos theta) ] / sin(theta)
    #
    # Use the recurrence:  (careful with 1-z^2 = sin^2 theta)
    # -sin(theta) * dP_n^m/dz = [ n*z*P_n^m - (n+m)*P_{n-1}^m ] / sin(theta)
    # when sin(theta) != 0
    z = jnp.cos(theta)
    sin_theta = jnp.sin(theta)
    safe_sin = jnp.where(regular, sin_theta, 1.0)

    pnm = associated_legendre(n, m, z)
    if n >= 1:
        pn1m = associated_legendre(n - 1, m, z)
    else:
        pn1m = jnp.zeros_like(z)

    general = (n * z * pnm - (n + m) * pn1m) / safe_sin
    result = jnp.where(regular, general, result)

    return result


# ---------------------------------------------------------------------------
# Emn normalization factor
# ---------------------------------------------------------------------------

def Emn(m, n):
    r"""VSH normalization factor E_{m,n}.

    Emn(m, n) = i^n * sqrt((2n+1) * (n-m)! / (n*(n+1)*(n+m)!))
    """
    val = 1j ** n * np.sqrt(
        (2 * n + 1) * _factorial(n - m)
        / (n * (n + 1) * _factorial(n + m))
    )
    return complex(val)


# ---------------------------------------------------------------------------
# Wigner 3j symbols (pure NumPy -- integer arguments, precomputation only)
# ---------------------------------------------------------------------------

def wigner_3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol using the Racah formula.

    All arguments are integers. Returns a float.
    """
    # Selection rules
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    # Triangle condition
    if j1 + j2 + j3 < 0:
        return 0.0

    # Racah formula
    # w3j = (-1)^(j1-j2-m3) * sqrt(delta) * sum_t ...
    # delta = (j1+j2-j3)!(j1-j2+j3)!(-j1+j2+j3)! / (j1+j2+j3+1)!
    t1 = j1 + j2 - j3
    t2 = j1 - j2 + j3
    t3 = -j1 + j2 + j3
    t4 = j1 + j2 + j3 + 1

    if t1 < 0 or t2 < 0 or t3 < 0:
        return 0.0

    delta = _factorial(t1) * _factorial(t2) * _factorial(t3) / _factorial(t4)

    # Prefactor from m-dependent factorials
    pre = (_factorial(j1 + m1) * _factorial(j1 - m1) *
           _factorial(j2 + m2) * _factorial(j2 - m2) *
           _factorial(j3 + m3) * _factorial(j3 - m3))

    # Sum over t
    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

    s = 0.0
    for t in range(t_min, t_max + 1):
        num = (-1) ** t
        den = (_factorial(t) *
               _factorial(j1 + j2 - j3 - t) *
               _factorial(j1 - m1 - t) *
               _factorial(j2 + m2 - t) *
               _factorial(j3 - j2 + m1 + t) *
               _factorial(j3 - j1 - m2 + t))
        s += num / den

    sign = (-1) ** (j1 - j2 - m3)
    return float(sign * np.sqrt(delta * pre) * s)


def wigner_3j_batch(j2, j3, m1, m2, m3):
    """Batch Wigner 3j symbols via Schulten-Gordon three-term recurrence.

    Returns (jmin, jmax, values) where values[i] = wigner_3j(jmin+i, j2, j3, m1, m2, m3).
    """
    if m1 + m2 + m3 != 0:
        return (0, -1, np.array([]))

    jmin = max(abs(j2 - j3), abs(m1))
    jmax = j2 + j3
    size = jmax - jmin + 1

    if size <= 0:
        return (jmin, jmax, np.array([]))

    values = np.zeros(size)

    if size <= 2:
        for i in range(size):
            values[i] = wigner_3j(jmin + i, j2, j3, m1, m2, m3)
        return (jmin, jmax, values)

    # Seed with 2 direct calls
    values[0] = wigner_3j(jmin, j2, j3, m1, m2, m3)
    values[1] = wigner_3j(jmin + 1, j2, j3, m1, m2, m3)

    # Forward recursion (Schulten-Gordon)
    d23 = j2 - j3
    s23p1 = j2 + j3 + 1
    j2j3_term = float(j2) * (j2 + 1) - float(j3) * (j3 + 1)
    m_diff = m3 - m2

    for i in range(1, size - 1):
        j = jmin + i
        j_d = float(j)
        jp1_d = j + 1.0

        sj_sq = (j_d**2 - float(d23)**2) * (float(s23p1)**2 - j_d**2) * (j_d**2 - float(m1)**2)
        sj = np.sqrt(max(0.0, sj_sq))

        sjp1_sq = (jp1_d**2 - float(d23)**2) * (float(s23p1)**2 - jp1_d**2) * (jp1_d**2 - float(m1)**2)
        sjp1 = np.sqrt(max(0.0, sjp1_sq))

        A_coeff = jp1_d * sj
        B_coeff = (2.0 * j + 1) * (-float(m1) * j2j3_term + m_diff * j_d * jp1_d)
        C_coeff = j_d * sjp1

        if abs(C_coeff) < 1e-30:
            values[i + 1] = wigner_3j(j + 1, j2, j3, m1, m2, m3)
        else:
            values[i + 1] = (-A_coeff * values[i - 1] - B_coeff * values[i]) / C_coeff

    return (jmin, jmax, values)


def gaunt_batch(m, n, u, v):
    """Batch Gaunt coefficients a(q) and b(q).

    Returns (a_vals, b_vals) NumPy arrays.
    """
    qmax_A = min(n, v, (n + v - abs(m + u)) // 2)
    qmax_B = min(n, v, (n + v + 1 - abs(m + u)) // 2)

    a_vals = np.zeros(qmax_A + 1)
    b_vals = np.zeros(qmax_B + 1)
    # b_vals[0] is always 0

    # Batch Wigner 3j sequences
    # W3j(n,v,p, 0,0,0) = W3j(p,n,v, 0,0,0) by column permutation symmetry
    w3j_zero = wigner_3j_batch(n, v, 0, 0, 0)
    # W3j(n,v,p, m,u,-(m+u)) = W3j(p,n,v, -(m+u),m,u)
    w3j_mu = wigner_3j_batch(n, v, -(m + u), m, u)

    fac_nm = _factorial(n + m)
    fac_n_m = _factorial(n - m)
    fac_vu = _factorial(v + u)
    fac_v_u = _factorial(v - u)
    sign = (-1) ** (m + u)

    # a_func values
    for q in range(qmax_A + 1):
        p = n + v - 2 * q
        fac_ratio = fac_nm * fac_vu * _factorial(p - m - u) / (fac_n_m * fac_v_u * _factorial(p + m + u))
        factor = sign * (2 * p + 1) * np.sqrt(fac_ratio)

        w1 = w3j_zero[2][p - w3j_zero[0]]
        w2 = w3j_mu[2][p - w3j_mu[0]]
        a_vals[q] = factor * w1 * w2

    # b_func values
    for q in range(1, qmax_B + 1):
        p = n + v - 2 * q
        fac_ratio = (fac_nm * fac_vu * _factorial(p - m - u + 1) /
                     (fac_n_m * fac_v_u * _factorial(p + m + u + 1)))
        factor = sign * (2 * p + 3) * np.sqrt(fac_ratio)

        w1 = w3j_zero[2][p - w3j_zero[0]]
        w2 = w3j_mu[2][(p + 1) - w3j_mu[0]]
        b_vals[q] = factor * w1 * w2

    return (a_vals, b_vals)
