from functools import lru_cache

import numpy as np
from scipy import special
from sympy.physics.wigner import gaunt

import miepy

# w = miepy.vsh.wigner_3j(1,1,1, 1,-1,0)
# print(w)
# w = wigner_3j(1,1,1,1,-1,0)
# print(w)
# from IPython import embed; embed()


@lru_cache(None)
def gaunt(m, n, u, v, p):
    """Gaunt coefficient."""
    def f(n):
        return special.gamma(n + 1)
    numerator = f(n + m) * f(v + u) * f(p - m - u)
    denominator = f(n - m) * f(v - u) * f(p + m + u)
    factor = (-1) ** (m + u) * (2 * p + 1) * (numerator / denominator) ** 0.5

    w1 = miepy.vsh.wigner_3j(n, v, p, 0, 0, 0)
    w2 = miepy.vsh.wigner_3j(n, v, p, m, u, -m - u)

    return factor * w1 * w2


@lru_cache(None)
def b_func(m, n, u, v, p):
    """B function."""
    def f(n):
        return special.gamma(n + 1)
    numerator = f(n + m) * f(v + u) * f(p - m - u + 1)
    denominator = f(n - m) * f(v - u) * f(p + m + u + 1)
    factor = (-1) ** (m + u) * (2 * p + 3) * (numerator / denominator) ** 0.5

    w1 = miepy.vsh.wigner_3j(n, v, p, 0, 0, 0)
    w2 = miepy.vsh.wigner_3j(n, v, p + 1, m, u, -m - u)

    return factor * w1 * w2


# g = gaunt(1,2,0,1,1)
# print(g)
# g = miepy.vsh.a_func(1,2,0,1,1)
# print(g)
# from IPython import embed; embed()


def A_translation(m, n, u, v, r, theta, phi, k):
    m *= -1
    def f(n):
        return special.gamma(n + 1)
    numerator = (2 * v + 1) * f(n - m) * f(v - u)
    denominator = 2 * n * (n + 1) * f(n + m) * f(v + u)

    factor = (-1) ** m * numerator / denominator * np.exp(1j * (u + m) * phi)

    qmax = min(n, v, (n + v - abs(m + u)) // 2)
    sum_term = 0
    for q in range(0, qmax + 1):
        p = n + v - 2 * q
        aq = gaunt(m, n, u, v, p)
        A = 1j**p * (n * (n + 1) + v * (v + 1) - p * (p + 1)) * aq

        Pnm = miepy.vsh.associated_legendre(p, u + m)
        sum_term += A * miepy.vsh.spherical_hn(p, k * r) * Pnm(np.cos(theta))

    return factor * sum_term


def B_translation(m, n, u, v, r, theta, phi, k):
    m *= -1
    def f(n):
        return special.gamma(n + 1)
    numerator = (2 * v + 1) * f(n - m) * f(v - u)
    denominator = 2 * n * (n + 1) * f(n + m) * f(v + u)

    factor = (-1) ** (m + 1) * numerator / denominator * np.exp(1j * (u + m) * phi)

    qmax = min(n, v, (n + v + 1 - abs(m + u)) // 2)
    sum_term = 0
    for q in range(1, qmax + 1):
        p = n + v - 2 * q
        bq = b_func(m, n, u, v, p)
        A = 1j ** (p + 1) * (((p + 1) ** 2 - (n - v) ** 2) * ((n + v + 1) ** 2 - (p + 1) ** 2)) ** 0.5 * bq

        Pnm = miepy.vsh.associated_legendre(p + 1, u + m)
        sum_term += A * miepy.vsh.spherical_hn(p + 1, k * r) * Pnm(np.cos(theta))

    return factor * sum_term


A = A_translation(8, 10, -9, 12, 2, 0.5, 0.5, 1)
B = B_translation(8, 10, -9, 12, 2, 0.5, 0.5, 1)
print(f"A: {A:.10e}", f"B: {B:.10e}", "", sep="\n")

A = A_translation(0, 10, 0, 10, 2, 0.5, 0.5, 1)
B = B_translation(0, 10, 0, 10, 2, 0.5, 0.5, 1)
print(f"A: {A:.10e}", f"B: {B:.10e}", "", sep="\n")

A = A_translation(-2, 11, 3, 9, 2, 0.5, 0.5, 1)
B = B_translation(-2, 11, 3, 9, 2, 0.5, 0.5, 1)
print(f"A: {A:.10e}", f"B: {B:.10e}", "", sep="\n")

A = A_translation(-2, 6, -2, 10, 2, 0.5, 0.5, 1)
B = B_translation(-2, 6, -2, 10, 2, 0.5, 0.5, 1)
print(f"A: {A:.10e}", f"B: {B:.10e}", "", sep="\n")

# A = A_translation(0, 1, 1, 1, 2, 0.5, 0.5, 1)
# print(A)
# A = miepy.vsh.A_translation(0, 1, 1, 1, 2, 0.5, 0.5, 1)
# print(A)


from IPython import embed

embed()
