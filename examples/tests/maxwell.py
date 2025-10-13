"""Displaying the fields in an xy cross section of the sphere (x polarized light, z-propagating)."""

import numpy as np

from miepy import sphere
from miepy.materials import material

# wavelength from 400nm to 1000nm
wav = np.linspace(300, 1100, 1000)

# create a material with n = 3.7 (eps = n^2) at all wavelengths
eps = 3.7**2 * np.ones(1000)
mu = 1 * np.ones(1000)
dielectric = material(wav, eps, mu)  # material object

# calculate scattering coefficients
rad = 100  # 100 nm radius
Nmax = 1  # Use up to 10 multipoles
m = sphere(Nmax, dielectric, rad)

E_func = m.E_field(999)
H_func = m.H_field(999)

eps = 0.01
x0 = 0 * rad
y0 = 10000 * rad
z0 = 0 * rad

x = np.linspace(x0, x0 + eps, 2)
y = np.linspace(y0, y0 + eps, 2)
z = np.linspace(z0, z0 + eps, 2)

X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
# X = X.T
# Y = Y.T
# Z = Z.T

R = (X**2 + Y**2 + Z**2) ** 0.5
THETA = np.arccos(Z / R)
PHI = np.arctan2(Y, X)

E = E_func(R, THETA, PHI)
E = np.squeeze(E)
Ex = E[0] * np.sin(THETA) * np.cos(PHI) + E[1] * np.cos(THETA) * np.cos(PHI) - E[2] * np.sin(PHI)
Ey = E[0] * np.sin(THETA) * np.sin(PHI) + E[1] * np.cos(THETA) * np.sin(PHI) + E[2] * np.cos(PHI)
Ez = E[0] * np.cos(THETA) - E[1] * np.sin(THETA)


H = H_func(R, THETA, PHI)
H = np.squeeze(H)
Hx = H[0] * np.sin(THETA) * np.cos(PHI) + H[1] * np.cos(THETA) * np.cos(PHI) - H[2] * np.sin(PHI)
Hy = H[0] * np.sin(THETA) * np.sin(PHI) + H[1] * np.cos(THETA) * np.sin(PHI) + H[2] * np.cos(PHI)
Hz = H[0] * np.cos(THETA) - H[1] * np.sin(THETA)

E = np.array([Ex, Ey, Ez])
H = np.array([Hx, Hy, Hz])
k = 2 * np.pi / wav[999]
print(k * 10000 * rad % (2 * np.pi) - np.pi)


def div(A):
    C1 = np.average(np.gradient(A[0], eps, axis=0))
    C2 = np.average(np.gradient(A[1], eps, axis=1))
    C3 = np.average(np.gradient(A[2], eps, axis=2))
    return C1 + C2 + C3


def curl(A):
    Cx = np.average(np.gradient(A[2], eps, axis=1)) - np.average(np.gradient(A[1], eps, axis=2))
    Cy = np.average(np.gradient(A[0], eps, axis=2)) - np.average(np.gradient(A[2], eps, axis=0))
    Cz = np.average(np.gradient(A[1], eps, axis=0)) - np.average(np.gradient(A[0], eps, axis=1))
    return np.array([Cx, Cy, Cz]) / k


from itertools import chain


def get_args(A):
    return list(chain.from_iterable(zip(np.abs(A), np.angle(A), strict=False)))


curlE = curl(E)
args = get_args(curlE)
print("curl(E) = ({:.2e} exp({:.2f}), {:.2e} exp({:.2f}), {:.2e} exp({:.2f}) )".format(*args))

avgH = 1j * np.average(H, axis=(1, 2, 3))
args = get_args(avgH)
print("iH      = ({:.2e} exp({:.2f}), {:.2e} exp({:.2f}), {:.2e} exp({:.2f}) )".format(*args))

print("")

curlH = curl(H)
args = get_args(curlH)
print("curl(H) = ({:.2e} exp({:.2f}), {:.2e} exp({:.2f}), {:.2e} exp({:.2f}) )".format(*args))

avgE = -1j * np.average(E, axis=(1, 2, 3))
args = get_args(avgE)
print("-iE     = ({:.2e} exp({:.2f}), {:.2e} exp({:.2f}), {:.2e} exp({:.2f}) )".format(*args))

print("")

divE = div(E)
print(f"div(E) = {np.abs(divE):.2e} exp({np.angle(divE):.2f})")
divH = div(H)
print(f"div(H) = {np.abs(divH):.2e} exp({np.angle(divH):.2f})")
# plt.streamplot(np.squeeze(X), np.squeeze(Z), np.real(Ex), np.real(Ez), color='white', linewidth=2)k
# plt.show()
