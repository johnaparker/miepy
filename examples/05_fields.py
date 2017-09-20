"""
Displaying the fields in an xy cross section of the sphere (x polarized light, z-propagating)
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

Ag = miepy.materials.predefined.Ag()

# calculate scattering coefficients, 800 nm illumination
radius = 200e-9      # 200 nm radius
Lmax = 5             # Use up to 5 multipoles
sphere = miepy.single_mie_sphere(radius, Ag, 800e-9, Lmax)

# create discretized xy plane
x = np.linspace(-2*radius,2*radius,100)
y = np.linspace(-2*radius,2*radius,100)
z = np.array([radius*0.0])

X,Y,Z = np.meshgrid(x,y,z, indexing='xy')
R = (X**2 + Y**2 + Z**2)**0.5
THETA = np.arccos(Z/R)
PHI = np.arctan2(Y,X)

# electric and magnetic field functions
E_func = sphere.E_field(index=0)
E = E_func(R,THETA,PHI).squeeze()
IE = np.sum(np.abs(E)**2, axis=0)
H_func = sphere.H_field(index=0)
H = H_func(R,THETA,PHI).squeeze()
IH = np.sum(np.abs(H)**2, axis=0)

# plot results
fig,axes = plt.subplots(ncols=2, figsize=plt.figaspect(1/2.7))
for i,ax in enumerate(axes):
    plt.subplot(ax)
    I = IE if i == 0 else IH
    plt.pcolormesh(np.squeeze(X)*1e9,np.squeeze(Y)*1e9, I, shading="gouraud", cmap=cm.viridis)
    plt.colorbar(label='field intensity')

THETA = np.squeeze(THETA)
PHI = np.squeeze(PHI)
for i,ax in enumerate(axes):
    F = E if i == 0 else H
    Fx = F[0]*np.sin(THETA)*np.cos(PHI) + F[1]*np.cos(THETA)*np.cos(PHI) - F[2]*np.sin(PHI)
    Fy = F[0]*np.sin(THETA)*np.sin(PHI) + F[1]*np.cos(THETA)*np.sin(PHI) + F[2]*np.cos(PHI)
    step=10
    ax.streamplot(np.squeeze(X)*1e9, np.squeeze(Y)*1e9, np.real(Fx), np.real(Fy), color='white', linewidth=1.0)

for ax in axes:
    ax.set(xlim=[-2*radius*1e9, 2*radius*1e9], ylim=[-2*radius*1e9, 2*radius*1e9],
            aspect='equal', xlabel="X (nm)", ylabel="Y (nm)")
axes[0].set_title("Electric Field")
axes[1].set_title("Magnetic Field")

plt.show()

# theta = np.linspace(0,np.pi,50)
# phi = np.linspace(0,2*np.pi,50)
# r = np.array([10000])

# R,THETA,PHI = np.meshgrid(r,theta,phi)
# X = R*np.sin(THETA)*np.cos(PHI)
# Y = R*np.sin(THETA)*np.sin(PHI)
# Z = R*np.cos(THETA)

# X = X.squeeze()
# Y = Y.squeeze()
# Z = Z.squeeze()

# E = E_func(R,THETA,PHI)
# I = np.sum(np.abs(E)**2, axis=0)
# I = np.squeeze(I)
# I -= np.min(I)
# I /= np.max(I)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# shape = X.shape
# C = np.zeros((shape[0], shape[1], 4))
# cmap_3d = cm.viridis
# for i in range(shape[0]):
    # for j in range(shape[1]):
        # C[i,j,:] = cmap_3d(I[i,j])
# surf = ax.plot_surface(X*1e9, Y*1e9, Z*1e9, rstride=1, cstride=1,shade=False, facecolors=C,linewidth=.0, edgecolors='#000000', antialiased=False)
# m = cm.ScalarMappable(cmap=cmap_3d)
# m.set_array(I)
# plt.colorbar(m)
# surf.set_edgecolor('k')
# ax.set_xlabel('X')
