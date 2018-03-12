"""
Far-field analysis with the GMT
"""

import numpy as np
import matplotlib.pyplot as plt
from my_pytools.my_matplotlib.colors import cmap
import miepy

nm = 1e-9
radius = 20*nm

r = 10000*nm
phi = np.linspace(0,2*np.pi,100)
theta = np.pi/2

R,PHI = np.meshgrid(r,phi, indexing='ij')
X = R*np.cos(PHI)
Y = R*np.sin(PHI)
Z = np.zeros_like(X)

system = miepy.gmt(miepy.spheres([0,0,0], radius, miepy.constant_material(1.3)), 
            miepy.sources.y_polarized_plane_wave(),
            600*nm, 2)

# E = np.squeeze(system.E_field(X,Y,Z,False))
# I = np.sum(np.abs(E)**2, axis=0)

Etheta = np.zeros(len(phi), dtype=complex)
Ephi = np.zeros(len(phi), dtype=complex)

k = system.material_data['k'][0]
factor = np.exp(1j*k*radius)/(-1j*k*radius)

system.solve_cluster_coefficients()
for r in range(system.rmax):
    n = system.n_indices[r]
    m = system.m_indices[r]
    Emn = miepy.vsh.Emn(m, n)

    tau = miepy.vsh.tau_func(n,m)(theta)
    pi = miepy.vsh.pi_func(n,m)(theta)

    Etheta += factor*Emn*(system.p_cluster[0,r]*tau + system.q_cluster[0,r]*pi)*np.exp(1j*m*phi)
    Ephi += 1j*factor*Emn*(system.p_cluster[0,r]*pi + system.q_cluster[0,r]*tau)*np.exp(1j*m*phi)

I = np.abs(Etheta)**2 + np.abs(Ephi)**2

fig,ax = plt.subplots(subplot_kw={'projection':'polar'})
ax.plot(phi,I)

plt.show()

