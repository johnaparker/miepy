"""
GMT for dimer, varying the separation distance between the dimer pair, with and without interactions
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
from tqdm import tqdm

Ag = miepy.materials.predefined.Ag()
radius = 75e-9
source = miepy.sources.x_polarized_plane_wave(amplitude=1)
separations = np.linspace(2*radius+10e-9,2*radius+700e-9, 50)

force1 = []
force2 = []

torque1 = []
torque2 = []

spheres = miepy.spheres([[separations[0]/2,0,0], [-separations[0]/2,0,0]], radius, Ag)
sol1 = miepy.gmt(spheres, source, 600e-9, 3, interactions=False)
sol2 = miepy.gmt(spheres, source, 600e-9, 1, interactions=True)

for separation in tqdm(separations):
    sol1.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))
    F,T = map(np.squeeze, sol1.force_on_particle(0))
    force1.append(F[1])
    torque1.append(T[2])

    sol2.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))
    F,T = map(np.squeeze, sol2.force_on_particle(0))
    force2.append(F[1])
    torque2.append(T[2])

plt.figure(1)
plt.plot(separations*1e9, force1, label="Fx, no interactions")
plt.plot(separations*1e9, force2, label="Fx, interactions")
plt.axhline(y=0, color='black')
plt.legend()

plt.figure(2)
plt.plot(separations*1e9, torque1, label="Tz, no interactions")
plt.plot(separations*1e9, torque2, label="Tz, interactions")
plt.legend()

plt.show()

