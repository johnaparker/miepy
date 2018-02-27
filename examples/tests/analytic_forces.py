"""
GMT for dimer, varying the separation distance between the dimer pair
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
import miepy
from math import factorial
from tqdm import tqdm
from scipy import constants

nm = 1e-9

Ag = miepy.materials.predefined.Ag()
radius = 75*nm
source = miepy.sources.rhc_polarized_plane_wave(amplitude=1)
separations = np.linspace(2*radius + 10*nm, 2*radius + 700*nm, 50)

spheres = miepy.spheres([[separations[0]/2,0,0], [-separations[0]/2,0,0]], radius, Ag)
mie = miepy.gmt(spheres, source, 600*nm, 2)

force = np.zeros((3,) + separations.shape)
torque = np.zeros_like(force)

def get_fz(self,i):
    fid = 0
    k = self.material_data['k'][fid]
    radius = self.spheres.radius[fid]
    E0 = self.source.amplitude
    eps_b = self.material_data['eps_b'][fid]*constants.epsilon_0
    mu_b = self.material_data['mu_b'][fid]


    Fz = 0
    A = -k**2*E0**2/(8*np.pi)


    p = np.zeros_like(self.p)
    q = np.zeros_like(self.q)
    p_inc = np.zeros_like(self.p_inc)
    q_inc = np.zeros_like(self.q_inc)

    for n in range(1,self.Lmax+1):
        for m in range(-n,n+1):
            Emn = miepy.vsh.Emn(m,n,1)
            alpha = np.sqrt(4*np.pi*factorial(n+m)/(2*n+1)/factorial(n-m))
            r = n**2 + n - 1 + m

            p_inc[...,r] = Emn/(1j*k**2)*alpha*self.p_inc[...,r]
            q_inc[...,r] = -Emn/(k**2)*alpha*np.sqrt(eps_b/mu_b)*self.q_inc[...,r]
            p[...,r] = -Emn/(1j*k**2)*alpha*self.p[...,r]
            q[...,r] = Emn/(k**2)*alpha*np.sqrt(eps_b/mu_b)*self.q[...,r]

    for n in range(1,self.Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m
            factor = A*m*eps_b**0.5

            Fz += factor*(2*p[fid,i,r]*np.conj(q[fid,i,r]) \
                    + p[fid,i,r]*np.conj(q_inc[fid,i,r]) \
                    + p_inc[fid,i,r]*np.conj(q[fid,i,r]))

            if n < self.Lmax:
                factor = A*n*(n+2)*np.sqrt((n-m+1)*(n+m+1)/(2*n+3)/(2*n+1))
                r1 = (n+1)**2 + (n+1) - 1 + m

                Fz += factor*(2*eps_b*p[fid,i,r1]*np.conj(p[fid,i,r]) \
                        + eps_b*p[fid,i,r1]*np.conj(p_inc[fid,i,r]) \
                        + eps_b*p_inc[fid,i,r1]*np.conj(p[fid,i,r]) \
                        + 2*q[fid,i,r1]*np.conj(q[fid,i,r]) \
                        + q[fid,i,r1]*np.conj(q_inc[fid,i,r]) \
                        + q_inc[fid,i,r1]*np.conj(q[fid,i,r]))

    return 4*np.pi*np.imag(Fz)

for i, separation in enumerate(tqdm(separations)):
    mie.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))
    force[2,i] = get_fz(mie, 0)

fig, ax = plt.subplots()
with h5py.File('temp.h5', 'r') as f:
    F = f['F'][...]
    ax.plot(separations/nm, F, 'o', label='Numerical Stress Tensor')


ax.plot(separations/nm, force[2], label='Analytic Equation')
ax.set(xlabel='separation (nm)', ylabel='force')
plt.legend()
plt.show()

