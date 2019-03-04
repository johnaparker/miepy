"""
Verify that the near fields, when evaluated far from the origin, approximate the expressions for the angular far fields
"""

import numpy as np
import miepy
import pytest

### parameters
nm = 1e-9
wav = 600*nm
k = 2*np.pi/wav
width = 100*nm
polarization = [1,1j]

### angular grid
radius = 150.3*wav
theta = np.linspace(0., np.pi, 4)
phi = np.linspace(0, 2*np.pi, 5)[:-1]
THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
X, Y, Z = miepy.coordinates.sph_to_cart(radius, THETA, PHI)

@pytest.mark.parametrize("source,atol,rtol", [
    (miepy.sources.gaussian_beam(width=width, polarization=polarization), 0, 2e-2),
    (miepy.sources.hermite_gaussian_beam(2, 0, width=width, polarization=polarization), 0, 2e-2),
    (miepy.sources.laguerre_gaussian_beam(1, 1, width=width, polarization=polarization), 1, 5e-2),
])
def test_source_electric_field_near_to_far(source, atol, rtol):
    """
    Compare E-field of source in far field using near and far field expressions
    Expressions are expected to converge in the limit r -> infinity
    """
    E1 = source.E_field(X, Y, Z, k, sampling=300)
    E1 = miepy.coordinates.vec_cart_to_sph(E1, THETA, PHI)[1:]
    E2 = source.E_angular(THETA, PHI, k, radius=radius)

    assert np.allclose(E1, E2, atol=atol, rtol=rtol)


@pytest.mark.parametrize("source,atol,rtol", [
    (miepy.sources.gaussian_beam(width=width, polarization=polarization), 0, 2e-2),
    (miepy.sources.hermite_gaussian_beam(2, 0, width=width, polarization=polarization), 0, 2e-2),
    (miepy.sources.laguerre_gaussian_beam(1, 1, width=width, polarization=polarization), 1, 5e-2),
])
def test_source_magnetic_field_near_to_far(source, atol, rtol):
    """
    Compare H-field of source in far field using near and far field expressions
    Expressions are expected to converge in the limit r -> infinity
    """
    H1 = source.H_field(X, Y, Z, k, sampling=300)
    H1 = miepy.coordinates.vec_cart_to_sph(H1, THETA, PHI)[1:]
    H2 = source.H_angular(THETA, PHI, k, radius=radius)

    assert np.allclose(H1, H2, atol=atol, rtol=rtol)

def test_cluster_field_near_to_far():
    """
    Compare scattered E/H-field of a cluster in far field using near and far field expressions
    Expressions are expected to converge in the limit r -> infinity
    """
    x = np.linspace(-600*nm, 600*nm, 3)
    y = np.linspace(-600*nm, 600*nm, 3)

    cluster = miepy.sphere_cluster(position=[[xv, yv, 0] for xv in x for yv in y],
                                   radius=100*nm,
                                   material=miepy.constant_material(index=2),
                                   wavelength=wav,
                                   source=miepy.sources.plane_wave([1,1]),
                                   lmax=3)

    theta = np.linspace(0, np.pi, 5)
    phi = np.linspace(0, 2*np.pi, 5)
    THETA, PHI = np.meshgrid(theta, phi)
    radius = np.ones_like(THETA)

    E1 = cluster.E_field(radius, THETA, PHI, spherical=True, source=False)
    E2 = cluster.E_angular(THETA, PHI, radius=radius, source=False)

    H1 = cluster.H_field(radius, THETA, PHI, spherical=True, source=False)
    H2 = cluster.H_angular(THETA, PHI, radius=radius, source=False)

    assert np.allclose(E1[0], 0, atol=1e-10), 'radial component of E goes to 0'
    assert np.allclose(E1[1:], E2, atol=0, rtol=1e-6), 'E converges'
    assert np.allclose(H1[0], 0, atol=1e-10), 'radial component of H goes to 0'
    assert np.allclose(H1[1:], H2, atol=0, rtol=1e-6), 'H converges'
