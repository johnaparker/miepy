import numpy as np
import matplotlib.pyplot as plt
from my_pytools.my_matplotlib.colors import cmap
import miepy

nm = 1e-9
radius = 20*nm



r = 10000*nm
phi = np.linspace(0,2*np.pi,100)

R,PHI = np.meshgrid(r,phi, indexing='ij')
X = R*np.cos(PHI)
Y = R*np.sin(PHI)
Z = np.zeros_like(X)

system = miepy.gmt(miepy.spheres([0,0,0], radius, miepy.constant_material(1.3)), 
            miepy.sources.y_polarized_plane_wave(),
            600*nm, 2, interactions=False)

E = np.squeeze(system.E_field(X,Y,Z,False))
I = np.sum(np.abs(E)**2, axis=0)



fig,ax = plt.subplots(subplot_kw={'projection':'polar'})

ax.plot(phi,I)
plt.show()

