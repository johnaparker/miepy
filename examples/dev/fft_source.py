import matplotlib.pyplot as plt
import numpy as np

import miepy

x = np.linspace(-2, 2, 70)
y = np.linspace(-2, 2, 70)
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)
k = 2 * np.pi / 1

source = miepy.sources.gaussian_beam(width=0.2, polarization=[1, 0])
E = source.E_field(X, Y, Z, k)
I = np.sum(np.abs(E[2:]) ** 2, axis=0)

fig, axes = plt.subplots(ncols=2, figsize=(12, 5))
ax = axes[0]
ax.pcolormesh(X, Y, I**0.5)
ax.set_aspect("equal")


kx = np.fft.fftfreq(len(x), d=x[1] - x[0])
ky = np.fft.fftfreq(len(y), d=y[1] - y[0])

res = 5
kmax = 2 * np.pi * k * res / 2
kx = np.fft.ifftshift(np.linspace(-kmax, kmax, 51 * res))
ky = np.fft.ifftshift(np.linspace(-kmax, kmax, 51 * res))
Kx, Ky = np.meshgrid(kx, ky)
Kz = np.sqrt((k**2 - Kx**2 - Ky**2).astype(complex))

_, THETA, PHI = miepy.coordinates.cart_to_sph(Kx, Ky, -Kz.real)
print(THETA[:, 0])
THETA[THETA < np.pi / 2 + 0.01] = np.pi / 2 + 0.01
print(np.min(THETA), np.max(THETA))

integrand = source.angular_spectrum(np.pi - THETA, PHI, k) / Kz
idx = (Kx**2 + Ky**2) > k**2
integrand[:, idx] = 0
integrand = np.insert(integrand, 0, 0, axis=0)
integrand = miepy.coordinates.vec_sph_to_cart(integrand, THETA, PHI)

E = np.array([np.fft.ifftshift(np.fft.ifft2(integrand[i])) for i in range(3)])
# E = np.array([np.fft.ifft2(integrand[i]) for i in range(3)])

x = np.fft.ifftshift(np.fft.fftfreq(len(kx), kx[1] - kx[0]) * 2 * np.pi)
y = np.fft.ifftshift(np.fft.fftfreq(len(ky), ky[1] - ky[0]) * 2 * np.pi)
X, Y = np.meshgrid(x, y)

I = np.sum(np.abs(E[2:]) ** 2, axis=0)
ax = axes[1]
ax.pcolormesh(X, Y, I**0.5)
ax.set_aspect("equal")
ax.set(xlim=[-2, 2], ylim=[-2, 2])


plt.show()
