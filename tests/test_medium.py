"""
Test usage of the medium material
"""

import numpy as np
import miepy

nm = 1e-9

def test_medium_scaling_force(plot=False):
    """ensure that the radiation pressure on a single sphere scales with 
       background index as n^6 for small, high-index sphere"""

    n_b = np.linspace(1,5,10)
    F = np.zeros_like(n_b)
    for i,n in enumerate(n_b):
        sphere = miepy.sphere_cluster(position=[0,0,0],
                                      radius=2*nm,
                                      material=miepy.constant_material(40**2),
                                      medium=miepy.dielectric(n),
                                      source=miepy.sources.plane_wave.from_string(polarization='x'),
                                      wavelength=2800*nm,
                                      lmax=5)

        F[i] = sphere.force_on_particle(0)[2]

    F_fit = np.max(F)*(n_b/np.max(n_b))**6

    if not plot:
        L2 = np.linalg.norm(F - F_fit)/F.shape[0]
        avg = np.average(F + F_fit)/2
        print(L2, avg)
        assert L2 < 1e-2*avg
    else:
        plt.figure()
        plt.plot(n_b, F, label='exact force')
        plt.plot(n_b, F_fit, 'o', label='$n_b^6$ scaling')
        plt.xlabel('$n_b$')
        plt.ylabel('force')
        plt.legend()
        plt.title(test_medium_scaling_force.__name__, weight='bold')

def test_medium_cross_sections(plot=False):
    """verify cross-sections of an Au dimer in water by comparing with Poynting vector approach"""
    Nwav = 5
    Au = miepy.materials. Au()
    radius = 50*nm

    lmax = 2
    nb = 1.33
    medium = miepy.constant_material(nb**2)

    wavelengths = np.linspace(500*nm, 1000*nm, Nwav)
    separation = 2*radius + 40*nm
    source = miepy.sources.plane_wave.from_string(polarization='x')

    THETA, PHI = miepy.coordinates.sphere_mesh(20)
    R = (separation/2 + radius + 20*nm)*np.ones_like(THETA)
    X, Y, Z = miepy.coordinates.sph_to_cart(R, THETA, PHI)

    scat, absorb, extinct, C, A = (np.zeros_like(wavelengths) for i in range(5))
    for i,wavelength in enumerate(wavelengths):
        sol = miepy.sphere_cluster(position=[[separation/2,0,0], [-separation/2,0,0]],
                                   radius=radius,
                                   material=Au,
                                   source=source,
                                   medium=medium,
                                   wavelength=wavelength,
                                   lmax=lmax)

        scat[i], absorb[i], extinct[i] = sol.cross_sections()

        E = sol.E_field(X,Y,Z, source=False)
        H = sol.H_field(X,Y,Z, source=False)
        C[i] = miepy.flux.flux_from_poynting_sphere(E, H, R, eps=nb**2)

        E = sol.E_field(X,Y,Z, source=True)
        H = sol.H_field(X,Y,Z, source=True)
        A[i] = -miepy.flux.flux_from_poynting_sphere(E, H, R, eps=nb**2)

    if not plot:
        for a,b,tol in [(C,scat,1e-3), (A,absorb,1e-3), (C+A,extinct,1e-3)]:
            L2 = np.linalg.norm(a - b)/a.shape[0]
            avg = np.average(a + b)/2
            print(L2, avg)
            assert L2 < tol*avg
    else:
        plt.figure()
        line, = plt.plot(wavelengths/nm, C, 'o',   color='C0', label='scattering (Poynting)')
        line, = plt.plot(wavelengths/nm, A, 'o',   color='C1', label='absorbption (Poynting)')
        line, = plt.plot(wavelengths/nm, C + A, 'o',   color='C2', label='extinction (Poynting)')

        plt.axhline(color='black', linestyle='--')
        line, = plt.plot(wavelengths/nm, scat,    color='C0', label='scattering (analytic)')
        line, = plt.plot(wavelengths/nm, absorb,  color='C1', label='absorbption (analytic)')
        line, = plt.plot(wavelengths/nm, extinct, color='C2', label='extinction (analytic)')

        plt.xlabel('wavelength (nm)')
        plt.ylabel(r'cross-section ($\mu$m$^2$)')
        plt.legend()
        plt.title(test_medium_cross_sections.__name__, weight='bold')

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    test_medium_scaling_force(plot=True)
    test_medium_cross_sections(plot=True)
    plt.show()
