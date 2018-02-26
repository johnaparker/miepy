"""
GMT for dimer, varying the separation distance between the dimer pair
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
from math import factorial
from tqdm import tqdm

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
    eps_b = self.material_data['eps_b'][fid]


    Fz = 0
    A = -1/k**2*E0**2/(8*np.pi)

    for n in range(1,self.Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m
            Emn = miepy.vsh.Emn(m,n,1)*np.sqrt(4*np.pi*factorial(n+m)/(2*n+1)/factorial(n-m))
            factor = A*np.abs(Emn)**2*m*eps_b**0.5

            Fz += factor*(2*self.p[fid,i,r]*np.conj(self.q[fid,i,r]) \
                    + self.p[fid,i,r]*np.conj(self.q_inc[fid,i,r]) \
                    + self.p_inc[fid,i,r]*np.conj(self.q[fid,i,r]))

            if n < self.Lmax:
                Emn1 = miepy.vsh.Emn(m,n+1,1)*np.sqrt(4*np.pi*factorial(n+1+m)/(2*(n+1)+1)/factorial(n+1-m))
                factor = A*Emn1*np.conj(Emn)*n*(n+2)*np.sqrt((n-m+1)*(n+m+1)/(2*n+3)/(2*n+1))
                r1 = (n+1)**2 + (n+1) - 1 + m

                Fz += factor*(2*eps_b*self.p[fid,i,r1]*np.conj(self.p[fid,i,r]) \
                        + eps_b*self.p[fid,i,r1]*np.conj(self.p_inc[fid,i,r]) \
                        + eps_b*self.p_inc[fid,i,r1]*np.conj(self.p[fid,i,r]) \
                        + 2*self.q[fid,i,r1]*np.conj(self.q[fid,i,r]) \
                        + self.q[fid,i,r1]*np.conj(self.q_inc[fid,i,r]) \
                        + self.q_inc[fid,i,r1]*np.conj(self.q[fid,i,r]))

    return np.imag(Fz)

for i, separation in enumerate(tqdm(separations)):
    mie.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))
    force[2,i] = get_fz(mie, 0)

plt.plot(separations/nm, force[2], color='C2')
plt.show()

