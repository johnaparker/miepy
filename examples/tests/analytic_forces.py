"""
Comapre analytic force and torque expressions to integrated Maxwell stress tensor
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

def get_force(self,i):
    fid = 0
    k = self.material_data['k'][fid]
    E0 = self.source.amplitude
    eps_b = self.material_data['eps_b'][fid]*constants.epsilon_0
    mu_b = self.material_data['mu_b'][fid]


    Fxy = 0
    Fz = 0
    Axy = np.pi*E0**2/k**2
    Az = -2*np.pi*E0**2/k**2

    for n in range(1,self.Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m

            # Fxy, term 1/3
            if m != n:
                factor = Axy*np.sqrt((n+m+1)*(n-m))*(eps_b/mu_b)/(n+1)
                r1 = n**2 + n - 1 + m + 1
                Fxy += factor*(2*self.p[fid,i,r]*np.conj(self.q[fid,i,r1]) \
                         - self.p[fid,i,r]*np.conj(self.q_inc[fid,i,r1]) \
                         - self.p_inc[fid,i,r]*np.conj(self.q[fid,i,r1]) \
                         + 2*self.q[fid,i,r]*np.conj(self.p[fid,i,r1]) \
                         - self.q[fid,i,r]*np.conj(self.p_inc[fid,i,r1]) \
                         - self.q_inc[fid,i,r]*np.conj(self.p[fid,i,r1]))

            # Fz, term 1/2
            factor = Az*m*(eps_b/mu_b)/(n*(n+1))
            Fz += factor*(2*self.p[fid,i,r]*np.conj(self.q[fid,i,r]) \
                    - self.p[fid,i,r]*np.conj(self.q_inc[fid,i,r]) \
                    - self.p_inc[fid,i,r]*np.conj(self.q[fid,i,r]))


            if n < self.Lmax:
                # Fxy, term 2/3
                factor = -Axy*np.sqrt((n+m+2)*(n+m+1)*n*(n+2)/((2*n+3)*(2*n+1)))/(n+1)
                r1 = (n+1)**2 + (n+1) - 1 + m + 1
                Fxy += factor*(2*eps_b*self.p[fid,i,r]*np.conj(self.p[fid,i,r1]) \
                         - eps_b*self.p[fid,i,r]*np.conj(self.p_inc[fid,i,r1]) \
                         - eps_b*self.p_inc[fid,i,r]*np.conj(self.p[fid,i,r1]) \
                         + 2*eps_b/mu_b*self.q[fid,i,r]*np.conj(self.q[fid,i,r1]) \
                         - eps_b/mu_b*self.q[fid,i,r]*np.conj(self.q_inc[fid,i,r1]) \
                         - eps_b/mu_b*self.q_inc[fid,i,r]*np.conj(self.q[fid,i,r1]))

                # Fxy, term 3/3
                factor = Axy*np.sqrt((n-m+1)*(n-m+2)*n*(n+2)/((2*n+3)*(2*n+1)))/(n+1)
                r1 = (n+1)**2 + (n+1) - 1 + m - 1
                Fxy += factor*(2*eps_b*self.p[fid,i,r1]*np.conj(self.p[fid,i,r]) \
                         - eps_b*self.p[fid,i,r1]*np.conj(self.p_inc[fid,i,r]) \
                         - eps_b*self.p_inc[fid,i,r1]*np.conj(self.p[fid,i,r]) \
                         + 2*eps_b/mu_b*self.q[fid,i,r1]*np.conj(self.q[fid,i,r]) \
                         - eps_b/mu_b*self.q[fid,i,r1]*np.conj(self.q_inc[fid,i,r]) \
                         - eps_b/mu_b*self.q_inc[fid,i,r1]*np.conj(self.q[fid,i,r]))

                # Fz, term 2/2
                factor = Az/(n+1)*np.sqrt((n-m+1)*(n+m+1)*n*(n+2)/(2*n+3)/(2*n+1))
                r1 = (n+1)**2 + (n+1) - 1 + m
                Fz += factor*(2*eps_b*self.p[fid,i,r1]*np.conj(self.p[fid,i,r]) \
                        - eps_b*self.p[fid,i,r1]*np.conj(self.p_inc[fid,i,r]) \
                        - eps_b*self.p_inc[fid,i,r1]*np.conj(self.p[fid,i,r]) \
                        + eps_b/mu_b*2*self.q[fid,i,r1]*np.conj(self.q[fid,i,r]) \
                        - eps_b/mu_b*self.q[fid,i,r1]*np.conj(self.q_inc[fid,i,r]) \
                        - eps_b/mu_b*self.q_inc[fid,i,r1]*np.conj(self.q[fid,i,r]))



    return np.array([np.real(Fxy), np.imag(Fxy), np.real(Fz)])

def get_torque(self,i):
    fid = 0
    k = self.material_data['k'][fid]
    E0 = self.source.amplitude
    eps_b = self.material_data['eps_b'][fid]
    mu_b = self.material_data['mu_b'][fid]

    Tx = 0
    Ty = 0
    Tz = 0
    A = -2*np.pi*E0**2/k**3

    for n in range(1,self.Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m

            if m != n:
                # Tx
                factor = -A*np.sqrt((n-m)*(n+m+1))
                r1 = n**2 + n - 1 + m + 1
                Tx += factor*np.real(eps_b*self.p[fid,i,r]*np.conj(self.p[fid,i,r1]) \
                        + mu_b*self.q[fid,i,r]*np.conj(self.q[fid,i,r1]) \
                        -0.5*(eps_b*self.p[fid,i,r1]*np.conj(self.p_inc[fid,i,r]) \
                        + eps_b*self.p[fid,i,r]*np.conj(self.p_inc[fid,i,r1]) \
                        + mu_b*self.q[fid,i,r1]*np.conj(self.q_inc[fid,i,r]) \
                        + mu_b*self.q[fid,i,r]*np.conj(self.q_inc[fid,i,r1])))

                # Ty
                Ty += factor*np.imag(eps_b*self.p[fid,i,r]*np.conj(self.p[fid,i,r1]) \
                        + mu_b*self.q[fid,i,r]*np.conj(self.q[fid,i,r1]) \
                        +0.5*(eps_b*self.p[fid,i,r1]*np.conj(self.p_inc[fid,i,r]) \
                        - eps_b*self.p[fid,i,r]*np.conj(self.p_inc[fid,i,r1]) \
                        + mu_b*self.q[fid,i,r1]*np.conj(self.q_inc[fid,i,r]) \
                        - mu_b*self.q[fid,i,r]*np.conj(self.q_inc[fid,i,r1])))

            # Tz
            factor = A*m/n
            Tz += factor* (eps_b*np.abs(self.p[fid,i,r])**2 + mu_b*np.abs(self.q[fid,i,r])**2 \
                    - np.real(eps_b*self.p[fid,i,r]*np.conj(self.p_inc[fid,i,r]) \
                    + mu_b*self.q[fid,i,r]*np.conj(self.q_inc[fid,i,r])))

    return np.array([Tx, Ty, Tz])*constants.epsilon_0

for i, separation in enumerate(tqdm(separations)):
    mie.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))
    force[:,i] = get_force(mie, 0)
    torque[:,i] = get_torque(mie, 0)

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=plt.figaspect(2/3)*2)

with h5py.File('temp.h5', 'r') as f:
    F = f['F'][...]
    T = f['T'][...]

    for i in range(3):
        axes[0,i].plot(separations/nm, F[i], 'o', color=f'C{i}', label='Numerical Stress Tensor')
        axes[1,i].plot(separations/nm, T[i], 'o', color=f'C{i}', label='Numerical Stress Tensor')

for i in range(3):
    comp = ['x', 'y', 'z'][i]

    axes[0,i].plot(separations/nm, force[i], color=f'C{i}', label='Analytic Equation')
    axes[0,i].set_title(label=f'F{comp}', weight='bold')

    axes[1,i].plot(separations/nm, torque[i], color=f'C{i}', label='Analytic Equation')
    axes[1,i].set_title(label=f'T{comp}', weight='bold')

for ax in axes.flatten():
    ax.axhline(y=0, color='black', linestyle='--')
    ax.legend()

for ax in axes[1]:
    ax.set_xlabel('separation (nm)')

axes[0,0].set_ylabel('force (N)')
axes[1,0].set_ylabel('torque (mN)')

plt.tight_layout()
plt.show()
