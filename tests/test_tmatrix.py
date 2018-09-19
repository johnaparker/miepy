import numpy as np
import miepy

nm = 1e-9

radius = 75*nm
material = miepy.materials.Ag()
source = miepy.sources.plane_wave([1,0])
wavelength = 600*nm
lmax = 2
eps_b = 1.5


def test_tmatrix_sphere_is_sphere():
    """tmatrix method with spheres should be equivalent to sphere cluster"""

    position = [[-300*nm, 0, 0], [300*nm, 0, 0]]

    spheres = miepy.sphere_cluster(position=position,
                                   radius=radius,
                                   material=material,
                                   source=source,
                                   wavelength=wavelength,
                                   lmax=lmax,
                                   medium=miepy.constant_material(eps_b))

    particles = [miepy.sphere(pos, radius, material) for pos in position]

    cluster = miepy.cluster(particles=particles,
                            source=source,
                            wavelength=wavelength,
                            lmax=lmax,
                            medium=miepy.constant_material(eps_b))

    print(np.max(np.abs(spheres.p_inc - cluster.p_inc)))
    assert np.allclose(spheres.p_inc, cluster.p_inc, rtol=0, atol=2e-16)

def test_tmatrix_spheroid_is_sphere():
    """tmatrix of spheroid with aspect ratio 1 is equal to tmatrix of sphere"""
    sphere   = miepy.sphere([0,0,0], radius, material)
    spheroid = miepy.spheroid([0,0,0], radius, radius, material, tmatrix_lmax=4)

    T1 = sphere.compute_tmatrix(lmax, wavelength, eps_b)
    T2 = spheroid.compute_tmatrix(lmax, wavelength, eps_b)

    print(np.max(np.abs(T1-T2)))
    assert np.allclose(T1, T2, rtol=0, atol=3e-4)
