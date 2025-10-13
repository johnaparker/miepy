from math import factorial

import matplotlib.pyplot as plt
import numpy as np

import miepy

eps = 1e-3
theta = np.linspace(eps, np.pi - eps, 50)

n = 6
m = 5

plt.plot(theta, miepy.vsh.pi_func(n, m)(theta), label="pi", color="C0")
plt.plot(theta, miepy.vsh.tau_func(n, m)(theta), label="tau", color="C1")


def pi(theta, Nmax):
    values = np.zeros((Nmax + 1, 2 * Nmax + 1) + theta.shape)
    x = np.cos(theta)

    # n = 1, m = 1
    values[1, 2] = 1

    # n = 1, m = -1
    values[1, 0] = 2

    if Nmax > 1:
        # n = 2, m = 2,-2
        values[2, 4] = 6 * np.sqrt(1 - x**2)
        values[2, 0] = -1 / 4 * values[2, 4]

        # n = 2, m = 1,-1
        values[2, 3] = 3 * x
        values[2, 1] = 1 / 6 * values[2, 3]

    for n in range(3, Nmax + 1):
        nid = n

        # m = n, m = -n
        mid = n + n
        values[nid, mid] = (1 - x**2) ** 0.5 * n * (2 * n - 1) / (n - 1) * values[nid - 1, mid - 2]
        values[nid, n - n] = (-1) ** (n + 1) * 1 / factorial(n + n) * values[nid, mid]

        # m = -n+2, ..., n-2
        for m in range(0, n - 1):
            mid = n + m
            values[nid, mid] = (2 * n - 1) / (n - m) * x * values[nid - 1, mid - 1] - (n + m - 1) / (n - m) * values[
                nid - 2, mid - 2
            ]
            values[nid, n - m] = (-1) ** (m + 1) * factorial(n - m) / factorial(n + m) * values[nid, mid]

        # m = -n+1, m = n-1
        m = n - 1
        mid = n + m
        if m != 2:
            values[nid, mid] = (
                -m * (n + m - 1) * (n - m + 2) / (m - 2) * values[nid, mid - 1]
                + 2 * m * x / (1 - x**2) ** 0.5 * values[nid, mid - 2]
            )
        else:
            values[nid, mid] = (
                (values[nid, mid + 1] + (m + 1) * (n + m) * (n - m + 1) / (m - 1) * values[nid, mid - 1])
                * (1 - x**2) ** 0.5
                / (2 * (m + 1) * x)
            )

        values[nid, n - m] = (-1) ** (m + 1) * factorial(n - m) / factorial(n + m) * values[nid, mid]

    return values[1:]


values = pi(theta, n)

plt.plot(theta, values[n - 1, m + n], label="pi", color="C0", linestyle="--")
# plt.plot(theta, miepy.vsh.tau_func(n,m)(theta), label='tau', color='C1', linestyle='--')

plt.legend()
plt.show()
