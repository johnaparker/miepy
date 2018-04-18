"""
cluster cross-sections vs. particle cross-section when summing over particles
"""
import numpy as np
import miepy
import matplotlib.pyplot as plt
from tqdm import tqdm

nm = 1e-9


wavelengths = np.linspace(300*nm, 1000*nm, 100)

### particle cross-sections
C, A, E = (np.zeros_like(wavelengths) for _ in range(3))

### cluster cross-sections
Cp, Ap, Ep = (np.zeros_like(wavelengths) for _ in range(3))

for i, wavelength in enumerate(tqdm(wavelengths)):
    dimer = miepy.cluster(position=[[-100*nm,0,0], [100*nm, 0, 0]],
                          radius=75*nm,
                          material=miepy.constant_material(3.6**2),
                          source=miepy.sources.y_polarized_plane_wave(),
                          wavelength=wavelength,
                          Lmax=2)
    C[i], A[i], E[i] = dimer.cross_sections()

    for j in range(dimer.Nparticles):
        Cp[i] += np.sum(miepy.flux.cluster_cross_sections(dimer.p_scat[j],
                        dimer.p_src[j], dimer.material_data.k).scattering)
        Ap[i] += np.sum(miepy.flux.cluster_cross_sections(dimer.p_scat[j],
                        dimer.p_src[j], dimer.material_data.k).absorption)
        Ep[i] += np.sum(miepy.flux.cluster_cross_sections(dimer.p_scat[j],
                        dimer.p_src[j], dimer.material_data.k).extinction)

fig, axes = plt.subplots(ncols=3, figsize=(11,4))

### apparently...., scattering is Cp + Ap
axes[0].plot(wavelengths/nm, C)
axes[0].plot(wavelengths/nm, Cp + Ap)

### apparently...., get absorption from (extinction - scattering)
axes[1].plot(wavelengths/nm, A)
axes[1].plot(wavelengths/nm, Ep - (Cp + Ap))

### apparently...., exctinction is correct
axes[2].plot(wavelengths/nm, E)
axes[2].plot(wavelengths/nm, Ep)

plt.show()
