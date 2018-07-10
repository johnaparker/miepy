"""
test cluster scattering coefficients for a monomer and dimer
"""
import numpy as np
import miepy
from tqdm import tqdm

nm = 1e-9
wavelengths = np.linspace(300*nm, 1000*nm, 5)

def test_cross_section_methods_monomer(plot=False):
    """test two cross-section methods for a single particle (monomer)"""
    C1 = np.zeros_like(wavelengths)
    C2 = np.zeros([len(wavelengths), 2, 2])

    for i,wavelength in enumerate(tqdm(wavelengths)):
        cluster = miepy.sphere_cluster(position=[0,0,0],
                                       radius=75*nm,
                                       material=miepy.constant_material(3.7**2),
                                       source=miepy.sources.plane_wave.from_string(polarization='x'),
                                       lmax=2,
                                       wavelength=wavelength)

        C1[i] = cluster.cross_sections().scattering
        C2[i] = cluster.cross_sections_per_multipole().scattering
    C2_sum =np.sum(C2, axis=(1,2))

    if not plot:
        L2 = np.linalg.norm(C1 - C2_sum)/C1.shape[0]
        avg = np.average(C1 + C2_sum)/2
        print(L2, avg)
        assert L2 < 1e-4*avg
    else:
        plt.figure()
        plt.plot(wavelengths/nm, C1, label='total scattering')
        plt.plot(wavelengths/nm, C2_sum, 'o', label='total scattering')
        plt.plot(wavelengths/nm, C2[:,0,0], label='eD')
        plt.plot(wavelengths/nm, C2[:,1,0], label='mD')
        plt.plot(wavelengths/nm, C2[:,0,1], label='eQ')
        plt.plot(wavelengths/nm, C2[:,1,1], label='mQ')

        plt.legend()

def test_cross_section_methods_dimer(plot=False):
    """test two cross-section methods for a dimer"""

    C1 = np.zeros_like(wavelengths)
    C2 = np.zeros([len(wavelengths), 2, 4])
    separation = 200*nm

    for i,wavelength in enumerate(tqdm(wavelengths)):
        cluster = miepy.sphere_cluster(position=[[-separation/2, 0, 0], [separation/2, 0, 0]],
                                       radius=75*nm,
                                       material=miepy.constant_material(3.7**2),
                                       source=miepy.sources.plane_wave.from_string(polarization='x'),
                                       lmax=2,
                                       wavelength=wavelength)

        C1[i] = cluster.cross_sections().scattering
        C2[i] = cluster.cross_sections_per_multipole(lmax=4).scattering
    C2_sum =np.sum(C2, axis=(1,2))

    if not plot:
        L2 = np.linalg.norm(C1 - C2_sum)/C1.shape[0]
        avg = np.average(C1 + C2_sum)/2
        print(L2, avg)
        assert L2 < 1e-3*avg
    else:
        plt.figure()

        plt.plot(wavelengths/nm, C1, label='total scattering')
        plt.plot(wavelengths/nm, C2_sum, 'o', label='total scattering')
        plt.plot(wavelengths/nm, C2[:,0,0], label='eD')
        plt.plot(wavelengths/nm, C2[:,1,0], label='mD')
        plt.plot(wavelengths/nm, C2[:,0,1], label='eQ')
        plt.plot(wavelengths/nm, C2[:,1,1], label='mQ')

        plt.legend()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_cross_section_methods_monomer(plot=True)
    test_cross_section_methods_dimer(plot=True)
    plt.show()
