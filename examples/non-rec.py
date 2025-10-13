"""GMT for dimer, varying the separation distance between the dimer pair, with and without interactions."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm

import miepy

pdf = PdfPages("out.pdf")

Ag = miepy.materials.Ag()
radius = [100e-9, 75e-9]

sources = {
    "x-polarized": miepy.sources.plane_wave.from_string(polarization="x", amplitude=1e6),
    "y-polarized": miepy.sources.plane_wave.from_string(polarization="y", amplitude=1e6),
    "rhc-polarized": miepy.sources.plane_wave.from_string(polarization="rhc", amplitude=1e6),
}
for name, source in sources.items():
    plt.figure()
    separations = np.linspace(sum(radius) + 10e-9, sum(radius) + 700e-9, 70)

    force1 = []
    force2 = []

    sol1 = miepy.sphere_cluster(
        position=[[separations[0] / 2, 0, 0], [-separations[0] / 2, 0, 0]],
        radius=radius,
        material=Ag,
        source=source,
        wavelength=600e-9,
        lmax=2,
        interactions=False,
    )
    sol2 = miepy.sphere_cluster(
        position=[[separations[0] / 2, 0, 0], [-separations[0] / 2, 0, 0]],
        radius=radius,
        material=Ag,
        source=source,
        wavelength=600e-9,
        lmax=2,
        interactions=True,
    )

    for separation in tqdm(separations):
        sol1.update_position(np.array([[separation / 2, 0, 0], [-separation / 2, 0, 0]]))
        F, T = sol1.force()
        force1.append(F[0])

        sol2.update_position(np.array([[separation / 2, 0, 0], [-separation / 2, 0, 0]]))
        F, T = sol2.force()
        force2.append(F[0])

    force1 = np.asarray(force1)
    force2 = np.asarray(force2)
    plt.plot((separations - sum(radius)) * 1e9, force1 * 1e12, label="no interactions")
    plt.plot((separations - sum(radius)) * 1e9, force2 * 1e12, label="interactions")
    plt.axhline(y=0, color="black")
    plt.xlabel("surface separation (nm)")
    plt.ylabel("force on CoM (pN)")
    plt.title(name)
    plt.legend()
    pdf.savefig()

pdf.close()
plt.show()
