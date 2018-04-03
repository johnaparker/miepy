"""
demonstration of a whispering gallery mode in a single dielectric sphere excited by a point-dipole near-field
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy

nm = 1e-9

f = miepy.sources.point_dipole([700*nm,1100*nm,0], polarization=[0,0,1])
x = np.linspace(-500*nm, 500*nm, 50)
y = np.linspace(200*nm, 1500*nm, 50)
X,Y = np.meshgrid(x, y, indexing='ij')
Z = np.zeros_like(X)
k = 2*np.pi/(400*nm)

E = f.E([X,Y,Z], k).squeeze()

I = np.sum(np.abs(E)**2, axis=0)

fig, ax = plt.subplots()
vmax = np.max(E[2].real)
ax.pcolormesh(X/nm, Y/nm, E[2].real, cmap='bwr', shading='gouraud', vmax=vmax, vmin=-vmax)
ax.set_aspect('equal')

plt.show()

