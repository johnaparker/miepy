"""
GMT for dimer, varying the separation distance between the dimer pair
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
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

for i, separation in enumerate(tqdm(separations)):
    mie.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))
    force[:,i] = mie.force_on_particle(0).squeeze()
    torque[:,i] = mie.torque_on_particle(0).squeeze()

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=plt.figaspect(2/3))

for i in range(3):
    comp = ['x', 'y', 'z'][i]

    axes[0,i].plot(separations/nm, force[i], color=f'C{i}')
    axes[0,i].set_title(label=f'F{comp}', weight='bold')

    axes[1,i].plot(separations/nm, torque[i], color=f'C{i}')
    axes[1,i].set_title(label=f'T{comp}', weight='bold')

for ax in axes.flatten():
    ax.axhline(y=0, color='black', linestyle='--')

for ax in axes[1]:
    ax.set_xlabel('separation (nm)')

axes[0,0].set_ylabel('force (N)')
axes[1,0].set_ylabel('torque (mN)')

plt.tight_layout()
plt.show()

