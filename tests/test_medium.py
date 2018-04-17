"""
Test usage of the medium material
"""

import numpy as np
import miepy

nm = 1e-9

def test_medium_scaling_force(plot=False):
    """ensure that the radiation pressure on a single sphere scales with 
       background index as n^5 for small, high-index sphere"""

    n_b = np.linspace(1,5,10)
    F = np.zeros_like(n_b)
    for i,n in enumerate(n_b):
        sphere = miepy.cluster(position=[0,0,0],
                               radius=2*nm,
                               material=miepy.constant_material(40**2),
                               medium=miepy.constant_material(n**2),
                               source=miepy.sources.x_polarized_plane_wave(),
                               wavelength=2800*nm,
                               Lmax=2)

        F[i] = sphere.force_on_particle(0)[2]

    F_fit = np.max(F)*(n_b/np.max(n_b))**5

    if not plot:
        L2 = np.linalg.norm(F - F_fit)/F.shape[0]
        avg = np.average(F + F_fit)/2
        print(L2, avg)
        assert L2 < 1e-2*avg
    else:
        plt.plot(n_b, F)
        plt.plot(n_b, F_fit)
        plt.show()

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_medium_scaling_force(plot=True)
