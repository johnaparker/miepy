import numpy as np
import matplotlib.pyplot as plt
from my_pytools.my_matplotlib.colors import cmap
import miepy

nm = 1e-9
sep = 50*nm
r = 20*nm

fig,axes = plt.subplots(nrows=2, figsize=plt.figaspect(2))

Nx = 180
Ny = 181
x = np.linspace(-sep/2 - 2*r, sep/2 + 2*r, Nx)
y = np.linspace(-2*r, 2*r, Ny)
z = 0
X,Y,Z = np.meshgrid(x,y,z, indexing='ij') 

for i,sim in enumerate([True,False]):
    system = miepy.gmt(miepy.spheres([[-sep/2,0,0],[sep/2,0,0]], r, miepy.materials.predefined.Ag()), 
                miepy.sources.y_polarized_plane_wave(),
                600*nm, 2, interactions=sim)
    E = np.squeeze(system.E_field(X,Y,Z,True))
    I = np.sum(np.abs(E)**2, axis=0)
    print(I.shape)

    mask = np.zeros((Nx,Ny), dtype=bool)
    mask[(np.squeeze(X)-sep/2)**2 + np.squeeze(Y)**2 < r**2] = True
    mask[(np.squeeze(X)+sep/2)**2 + np.squeeze(Y)**2 < r**2] = True
    I[mask] = 0

    im = axes[i].pcolormesh(np.squeeze(X)/nm,np.squeeze(Y)/nm,I, shading='gouraud', cmap=cmap['parula'], rasterized=True)
    plt.colorbar(im, ax=axes[i], label='Intensity')
    axes[i].set(aspect='equal', xlabel='x (nm)', ylabel='y (nm)')


plt.show()