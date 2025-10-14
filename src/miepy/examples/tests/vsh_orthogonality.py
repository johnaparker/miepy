import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import dblquad
from tqdm import tqdm

import miepy

nm = 1e-9

sampling = 50
THETA, PHI = miepy.vsh.sphere_mesh(sampling)
print(np.max(THETA))
tau = np.linspace(-1, 1, sampling)
phi = np.linspace(0, 2 * np.pi, 2 * sampling)

k = 2 * np.pi / (600 * nm)
N, M = miepy.vsh.VSH(1, 0)


### Test 1/r^2 dependence
def r_dependence():
    y = np.zeros(30)
    rvals = np.linspace(50 * nm, 1500 * nm, len(y))
    # rvals = np.linspace(5000*nm, 8000*nm, len(y))
    for i, r in enumerate(tqdm(rvals)):
        E1 = N(r, THETA, PHI, k)
        M(r, THETA, PHI, k)

        integrand = np.sum(E1 * np.conj(E1), axis=0)
        integral = miepy.vsh.misc.simps_2d(tau, phi, integrand).real
        y[i] = integral

    plt.plot(rvals, y)
    plt.plot(rvals, np.max(y) / (rvals / np.min(rvals)) ** 2)
    plt.show()


# r_dependence()

for n in range(1, 6):
    m = n
    N, M = miepy.vsh.VSH(n, m)
    r = 600 * nm
    E = N(r, THETA, PHI, k)

    integrand = np.sum(E * np.conj(E), axis=0)
    integral = miepy.vsh.misc.simps_2d(tau, phi, integrand).real
    err = 0

    def newF(theta, phi):
        E = N(r, theta, phi, k)
        return np.real(np.vdot(E, E)) * np.sin(theta)

    integral, err = dblquad(newF, 0, 2 * np.pi, lambda x: 0, lambda x: np.pi)

    Emn = miepy.vsh.Emn(m, n)

    # exact norm for M dot M
    zn = miepy.vsh.spherical_hn(n, k * r)
    radial_term = np.abs(zn) ** 2
    angular_term = 4 * np.pi * n * (n + 1) / np.abs(Emn)
    exact = angular_term * radial_term

    # exact norm for N dot N
    zn = miepy.vsh.spherical_hn(n, k * r)
    znp = miepy.vsh.spherical_hn(n, k * r, derivative=True)

    radial_term = (np.abs(zn + k * r * znp) ** 2 + n * (n + 1) * np.abs(zn) ** 2) / (k * r) ** 2
    exact = angular_term * radial_term

    print(f"n = {n}:  {integral:.5f} +/- {err:.5f} ({exact:.5f})")

# result = miepy.vsh.project_fields_onto(E, r, k, 'electric', 1, 0)
# print(result)
