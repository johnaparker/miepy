import numpy as np
import miepy
from my_pytools.my_numpy.integrate import simps_2d
import matplotlib.pyplot as plt
from tqdm import tqdm
from math import factorial
from scipy.integrate import dblquad

nm = 1e-9

sampling = 50
THETA, PHI = miepy.vsh.sphere_mesh(sampling)
print(np.max(THETA))
tau = np.linspace(-1, 1, sampling)
phi = np.linspace(0, 2*np.pi, 2*sampling)

k = 2*np.pi/(600*nm)
N, M = miepy.vsh.VSH(1,0)

### Test 1/r^2 dependence
def r_dependence():
    y = np.zeros(30)
    rvals = np.linspace(50*nm, 1500*nm, len(y))
    # rvals = np.linspace(5000*nm, 8000*nm, len(y))
    for i,r in enumerate(tqdm(rvals)):
        E1 = N(r, THETA, PHI, k)
        E2 = M(r, THETA, PHI, k)

        integrand  = np.sum(E1*np.conj(E1), axis=0)
        integral   = simps_2d(tau, phi, integrand).real
        y[i] = integral

    plt.plot(rvals, y)
    plt.plot(rvals, np.max(y)/(rvals/np.min(rvals))**2)
    plt.show()

r_dependence()

for n in range(1,6):
    m = 0
    N, M = miepy.vsh.VSH(n,m)
    r = 6000*nm
    E = N(r, THETA, PHI, k)

    integrand  = np.sum(E*np.conj(E), axis=0)
    integral   = simps_2d(tau, phi, integrand).real
    err = 0

    def newF(theta, phi):
        E = N(r, theta, phi, k)
        return np.vdot(E, E)*np.sin(theta)

    # integral,err = dblquad(newF, 0, 2*np.pi, lambda x: 0, lambda x: np.pi)

    Emn = miepy.vsh.Emn(m, n, 1)
    Emn = 1j**(n+2*m-1)/(2*np.pi**.5)*((2*n+1)*factorial(n-m)/factorial(n+m))**.5


    # exact = 4*np.pi*n*(n+1)/(np.abs(Emn)*k**2*r**2)
    # exact = 1/(np.abs(Emn)**2*k**2*r**2)
    exact = n*(n+1)/(np.abs(Emn)**2*k**2*r**2)
    print(f'n = {n}:  {integral:.5f} +/- {err:.5f} ({exact:.5f})')

# result = miepy.vsh.project_fields_onto(E, r, k, 'electric', 1, 0)
# print(result)