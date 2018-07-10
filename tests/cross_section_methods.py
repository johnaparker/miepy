"""
cluster cross-sections vs. particle cross-section when summing over particles
"""
import numpy as np
import miepy
import matplotlib.pyplot as plt
from tqdm import tqdm

nm = 1e-9


wavelengths = np.linspace(300*nm, 1000*nm, 30)

### particle cross-sections
C, A, E = (np.zeros_like(wavelengths) for _ in range(3))

### cluster cross-sections
Cp, Ap, Ep = (np.zeros_like(wavelengths) for _ in range(3))

for i, wavelength in enumerate(tqdm(wavelengths)):
    dimer = miepy.sphere_cluster(position=[[-100*nm,0,0], [100*nm, 0, 0]],
                                 radius=75*nm,
                                 material=miepy.constant_material(3.6**2 + 1j),
                                 source=miepy.sources.plane_wave.from_string(polarization='y'),
                                 wavelength=wavelength,
                                 lmax=2)
    C[i], A[i], E[i] = dimer.cross_sections()

    for j in range(dimer.Nparticles):
        ### absoprtion is calculated using p_inc
        Ap[i] += np.sum(miepy.flux.cluster_cross_sections(dimer.p_scat[j],
                        dimer.p_inc[j], dimer.material_data.k_b).absorption)

        ### extinction is calculated using p_src
        Ep[i] += np.sum(miepy.flux.cluster_cross_sections(dimer.p_scat[j],
                        dimer.p_src[j], dimer.material_data.k_b).extinction)

        ### scattering is the difference
        Cp[i] = Ep[i] - Ap[i]

fig, axes = plt.subplots(ncols=3, figsize=(11,4))

axes[0].plot(wavelengths/nm, C)
axes[0].plot(wavelengths/nm, Cp, 'o')

axes[1].plot(wavelengths/nm, A)
axes[1].plot(wavelengths/nm, Ap, 'o')

axes[2].plot(wavelengths/nm, E)
axes[2].plot(wavelengths/nm, Ep, 'o')

plt.show()
