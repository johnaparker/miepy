import numpy as np
import miepy
import pytest

nm = 1e-9

Ag = miepy.materials.Ag()
metal = miepy.materials.metal()
eps_b = 1.5
medium = miepy.constant_material(index=eps_b**2)

radius = 75*nm
source = miepy.sources.plane_wave([1,0])
wavelength = 600*nm
lmax = 2


@pytest.mark.parametrize("material,atol", [
    (Ag, 2e-16),
    (metal, 2e-17),
])
def test_tmatrix_sphere_is_sphere(material, atol):
    """tmatrix method with spheres should be equivalent to sphere cluster"""

    position = [[-300*nm, 0, 0], [300*nm, 0, 0]]

    spheres = miepy.sphere_cluster(position=position,
                                   radius=radius,
                                   material=material,
                                   source=source,
                                   wavelength=wavelength,
                                   lmax=lmax,
                                   medium=medium)

    particles = [miepy.sphere(pos, radius, material) for pos in position]

    cluster = miepy.cluster(particles=particles,
                            source=source,
                            wavelength=wavelength,
                            lmax=lmax,
                            medium=medium)

    print(np.max(np.abs(spheres.p_inc - cluster.p_inc)))
    assert np.allclose(spheres.p_inc, cluster.p_inc, rtol=0, atol=atol)


@pytest.mark.parametrize("material,atol", [
    (Ag, 3e-4),
    (metal, 5e-9),
])
def test_tmatrix_spheroid_is_sphere(material, atol):
    """tmatrix of spheroid with aspect ratio 1 is equal to tmatrix of sphere"""
    sphere   = miepy.sphere([0,0,0], radius, material)
    spheroid = miepy.spheroid([0,0,0], radius, radius, material, tmatrix_lmax=4)

    T1 = sphere.compute_tmatrix(lmax, wavelength, eps_b)
    T2 = spheroid.compute_tmatrix(lmax, wavelength, eps_b)

    print(np.max(np.abs(T1-T2)))
    assert np.allclose(T1, T2, rtol=0, atol=atol)


@pytest.mark.parametrize("material,atol", [
    (Ag, 8e-29),
    (metal, 5e-29),
])
def test_rotated_spheroid_equals_rotated_light(material, atol):
    """cross-sections are the same if a spheroid is rotated or if the light is rotated"""

    theta = 1.3
    phi = 0.7

    source = miepy.sources.plane_wave([1,0], theta=theta, phi=phi)
    z_oriented = miepy.cluster(particles=miepy.spheroid([0,0,0], radius, 2*radius, material),
                               source=source,
                               wavelength=wavelength,
                               lmax=lmax,
                               medium=medium)

    C1 = z_oriented.cross_sections()

    source = miepy.sources.plane_wave([1,0])
    q = miepy.quaternion.from_spherical_coords(theta, phi).inverse()
    oriented = miepy.cluster(particles=miepy.spheroid([0,0,0], radius, 2*radius, material, orientation=q),
                             source=source,
                             wavelength=wavelength,
                             lmax=lmax,
                             medium=medium)

    C2 = oriented.cross_sections()

    assert np.allclose(C1, C2, rtol=0, atol=atol)

def test_sphere_cluster_tmatrix_particle():
    """sphere_cluster_particle in miepy.cluster yields the same results as miepy.sphere_cluster"""
    L = 155*nm
    lmax = 6

    cluster = miepy.sphere_cluster(position=[[-L/2, 0,0], [L/2,0,0]],
                                   material=Ag,
                                   radius=75*nm,
                                   source=source,
                                   lmax=lmax,
                                   medium=medium,
                                   wavelength=800*nm)
    C1 = cluster.cross_sections()
    cluster.solve_cluster_coefficients()
    p1 = cluster.p_cluster

    particle = miepy.sphere_cluster_particle(cluster.position, cluster.radius, Ag, lmax=lmax)
    cluster = miepy.cluster(particles=particle,
                            source=source,
                            lmax=lmax,
                            medium=medium,
                            wavelength=800*nm)
    C2 = cluster.cross_sections()
    p2 = cluster.p_scat

    assert np.allclose(C1, C2, rtol=7e-4, atol=0), 'equal cross-sections'
    assert np.allclose(p1, p2, rtol=1e-2, atol=1e-12), 'equal cluster scattering coefficients'
