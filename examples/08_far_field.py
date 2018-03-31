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

x = r*np.cos(phi)
y = r*np.sin(phi)
z = np.zeros_like(x)

sol = miepy.cluster(position=[0,0,0],
                    radius=radius,
                    material=miepy.constant_material(1.3),
                    source=miepy.sources.y_polarized_plane_wave(),
                    wavelength=600*nm,
                    Lmax=1)

# E = sol.E_field(x,y,z,False)
# I = np.sum(np.abs(E)**2, axis=0)

Etheta = np.zeros(len(phi), dtype=complex)
Ephi = np.zeros(len(phi), dtype=complex)

k = sol.material_data.k
factor = np.exp(1j*k*r)/(-1j*k*r)

sol.solve_cluster_coefficients()
for r in range(sol.rmax):
    n = sol.n_indices[r]
    m = sol.m_indices[r]
    Emn = miepy.vsh.Emn(m, n)

    tau = miepy.vsh.tau_func(n,m)(theta)
    pi = miepy.vsh.pi_func(n,m)(theta)

    Etheta += factor*Emn*(sol.p_cluster[r]*tau + sol.q_cluster[r]*pi)*np.exp(1j*m*phi)
    Ephi += 1j*factor*Emn*(sol.p_cluster[r]*pi + sol.q_cluster[r]*tau)*np.exp(1j*m*phi)

I = np.abs(Etheta)**2 + np.abs(Ephi)**2

fig,ax = plt.subplots(subplot_kw={'projection':'polar'})
ax.plot(phi,I)

plt.show()

