import numpy as np
import miepy

nm = 1e-9

position = [[-300*nm, 0, 0], [300*nm, 0, 0]]
radius = 75*nm
material = miepy.materials.Ag()
source = miepy.sources.plane_wave([1,0])
wavelength = 600*nm
lmax = 2


spheres = miepy.sphere_cluster(position=position,
                               radius=radius,
                               material=material,
                               source=source,
                               wavelength=wavelength,
                               lmax=lmax)

sphere_cs = spheres.cross_sections()


def test_tmatrix_sphere_is_sphere():
    """tmatrix method with spheres should be equivalent to sphere cluster"""
    particles = [miepy.sphere(pos, radius, material) for pos in position]
    cluster = miepy.cluster(particles=particles,
                            source=source,
                            wavelength=wavelength,
                            lmax=lmax)

    assert np.allclose(spheres.p_inc, cluster.p_inc, rtol=0, atol=1e-17)

def test_tmatrix_spheroid_is_sphere():
    """tmatrix method with spheroids of aspect ratio = 1 should be equivalent to sphere cluster"""
    particles = [miepy.spheroid(pos, radius, radius, material, tmatrix_lmax=4) for pos in position]
    cluster = miepy.cluster(particles=particles,
                            source=source,
                            wavelength=wavelength,
                            lmax=lmax)

    assert np.allclose(spheres.p_inc, cluster.p_inc, rtol=0, atol=2e-5)
