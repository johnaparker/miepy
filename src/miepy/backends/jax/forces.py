"""Force and torque computations for the JAX backend.

Pure Python/NumPy port of cpp/src/forces.cpp. These operate on small O(rmax)
arrays so plain loops are fine — no JAX JIT needed.
"""

import numpy as np
from scipy import constants


def _rmax_to_lmax(rmax):
    return int(-1 + (1 + rmax)**0.5)


def force_jax(p_scat, p_inc, k, eps_b, mu_b):
    """Compute the optical force on a particle from expansion coefficients.

    Arguments:
        p_scat[2*rmax]  scattered field coefficients (flattened)
        p_inc[2*rmax]   incident field coefficients (flattened)
        k               wavenumber
        eps_b           background permittivity
        mu_b            background permeability

    Returns: F[3]
    """
    eps_0 = constants.epsilon_0
    rmax = len(p_scat) // 2
    lmax = _rmax_to_lmax(rmax)

    p = p_scat[:rmax]
    q = p_scat[rmax:]
    pi = p_inc[:rmax]
    qi = p_inc[rmax:]

    Axy = np.pi / k**2 * eps_0 * eps_b
    Az = -2 * Axy

    Fxy = 0 + 0j
    Fz = 0 + 0j

    for n in range(1, lmax + 1):
        for m in range(-n, n + 1):
            r = n**2 + n + m - 1

            # Fxy, term 1/3
            if m != n:
                factor = Axy * np.sqrt((n + m + 1) * (n - m)) / (n * (n + 1))
                r1 = r + 1
                Fxy += factor * (2.0 * p[r] * np.conj(q[r1])
                                 - p[r] * np.conj(qi[r1])
                                 - pi[r] * np.conj(q[r1])
                                 + 2.0 * q[r] * np.conj(p[r1])
                                 - q[r] * np.conj(pi[r1])
                                 - qi[r] * np.conj(p[r1]))

            # Fz, term 1/2
            factor = Az * m / (n * (n + 1))
            Fz += factor * (2.0 * p[r] * np.conj(q[r])
                            - p[r] * np.conj(qi[r])
                            - pi[r] * np.conj(q[r]))

            if n < lmax:
                # Fxy, term 2/3
                factor = (-Axy * np.sqrt((n + m + 2) * (n + m + 1) * n * (n + 2)
                          / ((2 * n + 3) * (2 * n + 1))) / (n + 1))
                r1 = (n + 1)**2 + (n + 1) + m
                Fxy += factor * (2.0 * p[r] * np.conj(p[r1])
                                 - p[r] * np.conj(pi[r1])
                                 - pi[r] * np.conj(p[r1])
                                 + 2.0 * q[r] * np.conj(q[r1])
                                 - q[r] * np.conj(qi[r1])
                                 - qi[r] * np.conj(q[r1]))

                # Fxy, term 3/3
                factor = (Axy * np.sqrt((n - m + 1) * (n - m + 2) * n * (n + 2)
                          / ((2 * n + 3) * (2 * n + 1))) / (n + 1))
                r1 = (n + 1)**2 + (n + 1) + m - 2
                Fxy += factor * (2.0 * p[r1] * np.conj(p[r])
                                 - p[r1] * np.conj(pi[r])
                                 - pi[r1] * np.conj(p[r])
                                 + 2.0 * q[r1] * np.conj(q[r])
                                 - q[r1] * np.conj(qi[r])
                                 - qi[r1] * np.conj(q[r]))

                # Fz, term 2/2
                factor = (Az * np.sqrt((n - m + 1) * (n + m + 1) * n * (n + 2)
                          / ((2 * n + 3) * (2 * n + 1))) / (n + 1))
                r1 = (n + 1)**2 + (n + 1) + m - 1
                Fz += factor * (2.0 * p[r1] * np.conj(p[r])
                                - p[r1] * np.conj(pi[r])
                                - pi[r1] * np.conj(p[r])
                                + 2.0 * q[r1] * np.conj(q[r])
                                - q[r1] * np.conj(qi[r])
                                - qi[r1] * np.conj(q[r]))

    return np.array([np.real(Fxy), np.imag(Fxy), np.real(Fz)])


def torque_jax(p_scat, p_inc, k, eps_b, mu_b):
    """Compute the optical torque on a particle from expansion coefficients.

    Arguments:
        p_scat[2*rmax]  scattered field coefficients (flattened)
        p_inc[2*rmax]   incident field coefficients (flattened)
        k               wavenumber
        eps_b           background permittivity
        mu_b            background permeability

    Returns: T[3]
    """
    eps_0 = constants.epsilon_0
    rmax = len(p_scat) // 2
    lmax = _rmax_to_lmax(rmax)
    A = -2 * np.pi / k**3 * eps_0 * eps_b

    p = p_scat[:rmax]
    q = p_scat[rmax:]
    pi = p_inc[:rmax]
    qi = p_inc[rmax:]

    T = np.zeros(3, dtype=float)

    for n in range(1, lmax + 1):
        for m in range(-n, n + 1):
            r = n**2 + n + m - 1

            if m != n:
                factor = -A * np.sqrt((n - m) * (n + m + 1))
                r1 = r + 1
                T[0] += factor * np.real(
                    p[r] * np.conj(p[r1]) + q[r] * np.conj(q[r1])
                    - 0.5 * (p[r1] * np.conj(pi[r]) + p[r] * np.conj(pi[r1])
                             + q[r1] * np.conj(qi[r]) + q[r] * np.conj(qi[r1])))
                T[1] += factor * np.imag(
                    p[r] * np.conj(p[r1]) + q[r] * np.conj(q[r1])
                    + 0.5 * (p[r1] * np.conj(pi[r]) - p[r] * np.conj(pi[r1])
                             + q[r1] * np.conj(qi[r]) - q[r] * np.conj(qi[r1])))

            factor = A * m
            T[2] += factor * (np.abs(p[r])**2 + np.abs(q[r])**2
                              - np.real(p[r] * np.conj(pi[r]) + q[r] * np.conj(qi[r])))

    return T
