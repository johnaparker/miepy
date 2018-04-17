"""
Far-field analysis with the GMT
"""

import numpy as np
import matplotlib.pyplot as plt
from my_pytools.my_matplotlib.colors import cmap
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import miepy

nm = 1e-9

sol = miepy.cluster(position=[[-100*nm,0,0], [100*nm, 0, 0]],
                    radius=75*nm,
                    material=miepy.constant_material(3.6**2),
                    source=miepy.sources.y_polarized_plane_wave(),
                    wavelength=600*nm,
                    Lmax=2)

### xy plane far-field
fig, ax = plt.subplots(subplot_kw={'projection':'polar'})

r = 10000*nm
phi = np.linspace(0, 2*np.pi, 100)
theta = np.pi/2

THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
E = sol.E_angular(THETA, PHI).squeeze()
I = np.sum(np.abs(E)**2, axis=0)

ax.plot(phi,I)

### 3D far-field
fig, ax = plt.subplots(subplot_kw={'projection':'3d'})

phi = np.linspace(0, 2*np.pi, 50)
theta = np.linspace(0, np.pi, 50)
THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
E = sol.E_angular(THETA, PHI).squeeze()
I = np.sum(np.abs(E)**2, axis=0)
I /= np.max(I)
X, Y, Z = miepy.coordinates.sph_to_cart(I, THETA, PHI)
colors = mpl.cm.viridis(I)

surface = ax.plot_surface(X, Y, Z, facecolors=colors, cmap='viridis', edgecolors='#000000',
        cstride=1, rstride=1, linewidth=.1, shade=False)
surface.set_edgecolor('k')

a = 1.5
ax.set(xlim=[-a,a], ylim=[-a,a], zlim=[-a,a],
        xlabel='x', ylabel='y', zlabel='z')
ax.view_init(ax.elev, ax.azim+90)

ax.contourf(X,Y,Z, zdir='x', offset=-a, cmap='viridis')
ax.contourf(X,Y,Z, zdir='y', offset=-a, cmap='viridis')
ax.contourf(X,Y,Z, zdir='z', offset=-a, cmap='viridis')

plt.show()

