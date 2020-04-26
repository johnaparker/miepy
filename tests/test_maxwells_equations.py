"""
Test the validy of Maxwell's equations
"""

import numpy as np
import miepy
import pytest

def div(A, eps):
    C1 = np.average(np.gradient(A[0], eps, axis=0))
    C2 = np.average(np.gradient(A[1], eps, axis=1))
    C3 = np.average(np.gradient(A[2], eps, axis=2))
    return C1 + C2 + C3

def curl(A, eps):
    Cx = np.average(np.gradient(A[2], eps, axis=1)) - np.average(np.gradient(A[1], eps, axis=2))
    Cy = np.average(np.gradient(A[0], eps, axis=2)) - np.average(np.gradient(A[2], eps, axis=0))
    Cz = np.average(np.gradient(A[1], eps, axis=0)) - np.average(np.gradient(A[0], eps, axis=1))
    return np.array([Cx,Cy,Cz])

nm = 1e-9
wav = 600*nm
k = 2*np.pi/wav
width = 800*nm
polarization = [1,1j]
theta, phi = .6, .3

@pytest.mark.parametrize("source", [
    miepy.sources.gaussian_beam(width, polarization, theta=theta, phi=phi),
    miepy.sources.hermite_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
    miepy.sources.laguerre_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
])
def test_maxwells_equations_source_near_field(source):
    """
    Verify Maxwell's equations for a source at a point in the near field
    """
    x0, y0, z0 = 50*nm, 100*nm, 150*nm
    eps = .1*nm

    x = np.linspace(x0, x0 + eps, 2)
    y = np.linspace(y0, y0 + eps, 2)
    z = np.linspace(z0, z0 + eps, 2)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    E_grid = source.E_field(X, Y, Z, k)
    E = np.average(E_grid, axis=(1,2,3))
    divE = div(E_grid, eps)
    curlE = curl(E_grid, eps)

    H_grid = source.H_field(X, Y, Z, k)
    H = np.average(H_grid, axis=(1,2,3))
    divH = div(H_grid, eps)
    curlH = curl(H_grid, eps)

    assert np.abs(divE/(k*np.linalg.norm(E))) < 1e-5, 'div(E) = 0'
    assert np.abs(divH/(k*np.linalg.norm(H))) < 1e-5, 'div(H) = 0'
    assert np.allclose(curlE, 1j*k*H, atol=0, rtol=1e-3), 'curl(E) = ikH'
    assert np.allclose(curlH, -1j*k*E, atol=0, rtol=1e-3), 'curl(H) = -ikE'


#TODO: this test fails if x,y,z = (0,0,1), i.e. theta = 0. Is this a problem?
#The definition of polarization in terms of spherical coordinates is discontinuous here
@pytest.mark.parametrize("source", [
    miepy.sources.gaussian_beam(width, polarization, theta=theta, phi=phi),
    miepy.sources.hermite_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
    miepy.sources.laguerre_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
])
def test_maxwells_equations_source_far_field(source):
    """
    Verify Maxwell's equations for a source at a point in the far field
    """
    x0, y0, z0 = 0.1, 0.1, 1
    eps = 1*nm

    x = np.linspace(x0, x0 + eps, 2)
    y = np.linspace(y0, y0 + eps, 2)
    z = np.linspace(z0, z0 + eps, 2)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)

    E_grid = source.E_angular(THETA, PHI, k, radius=R)
    E_grid = np.insert(E_grid, 0, 0, axis=0)
    E_grid = miepy.coordinates.vec_sph_to_cart(E_grid, THETA, PHI)

    E = np.average(E_grid, axis=(1,2,3))
    divE = div(E_grid, eps)
    curlE = curl(E_grid, eps)

    H_grid = source.H_angular(THETA, PHI, k, radius=R)
    H_grid = np.insert(H_grid, 0, 0, axis=0)
    H_grid = miepy.coordinates.vec_sph_to_cart(H_grid, THETA, PHI)

    H = np.average(H_grid, axis=(1,2,3))
    divH = div(H_grid, eps)
    curlH = curl(H_grid, eps)

    assert np.abs(divE/(k*np.linalg.norm(E))) < 1e-5, 'div(E) = 0'
    assert np.abs(divH/(k*np.linalg.norm(H))) < 1e-5, 'div(H) = 0'
    assert np.allclose(curlE[:2], 1j*k*H[:2], atol=0, rtol=1e-5), 'curl(E) = ikH'
    assert np.allclose(curlH[:2], -1j*k*E[:2], atol=0, rtol=1e-5), 'curl(H) = -ikE'



@pytest.mark.parametrize("source", [
    miepy.sources.plane_wave(polarization, theta=theta, phi=phi),
    miepy.sources.gaussian_beam(width, polarization, theta=theta, phi=phi),
    miepy.sources.hermite_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
    miepy.sources.laguerre_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
])
def test_maxwells_equations_cluster_near_field(source):
    """
    Verify Maxwell's equations for a cluster at a point in the near field
    """
    cluster = miepy.sphere_cluster(position=[[-400*nm, -200*nm, 0], [200*nm, 200*nm, 100*nm]],
                                   radius=100*nm,
                                   material=miepy.constant_material(index=3.7),
                                   source=source,
                                   wavelength=wav,
                                   lmax=2)

    x0, y0, z0 = 0, 0, 0
    eps = .1*nm

    x = np.linspace(x0, x0 + eps, 2)
    y = np.linspace(y0, y0 + eps, 2)
    z = np.linspace(z0, z0 + eps, 2)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    E_grid = cluster.E_field(X, Y, Z)
    E = np.average(E_grid, axis=(1,2,3))
    divE = div(E_grid, eps)
    curlE = curl(E_grid, eps)

    H_grid = cluster.H_field(X, Y, Z, k)
    H = np.average(H_grid, axis=(1,2,3))
    divH = div(H_grid, eps)
    curlH = curl(H_grid, eps)

    assert np.abs(divE/(k*np.linalg.norm(E))) < 1e-6, 'div(E) = 0'
    assert np.abs(divH/(k*np.linalg.norm(H))) < 1e-6, 'div(H) = 0'
    assert np.allclose(curlE, 1j*k*H, atol=0, rtol=1e-6), 'curl(E) = ikH'
    assert np.allclose(curlH, -1j*k*E, atol=0, rtol=1e-6), 'curl(H) = -ikE'


@pytest.mark.parametrize("source", [
    miepy.sources.plane_wave(polarization, theta=theta, phi=phi),
    miepy.sources.gaussian_beam(width, polarization, theta=theta, phi=phi),
    miepy.sources.hermite_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
    miepy.sources.laguerre_gaussian_beam(1, 1, width, polarization, theta=theta, phi=phi),
])
def test_maxwells_equations_cluster_far_field(source):
    """
    Verify Maxwell's equations for a cluster at a point in the far field
    """
    cluster = miepy.sphere_cluster(position=[[-400*nm, -200*nm, 0], [200*nm, 200*nm, 100*nm]],
                                   radius=100*nm,
                                   material=miepy.constant_material(index=3.7),
                                   source=source,
                                   wavelength=wav,
                                   lmax=2)

    x0, y0, z0 = 0, 0, 1
    eps = 1*nm

    x = np.linspace(x0, x0 + eps, 2)
    y = np.linspace(y0, y0 + eps, 2)
    z = np.linspace(z0, z0 + eps, 2)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)

    E_grid = cluster.E_angular(THETA, PHI, radius=R)
    E_grid = np.insert(E_grid, 0, 0, axis=0)
    E_grid = miepy.coordinates.vec_sph_to_cart(E_grid, THETA, PHI)

    E = np.average(E_grid, axis=(1,2,3))
    divE = div(E_grid, eps)
    curlE = curl(E_grid, eps)

    H_grid = cluster.H_angular(THETA, PHI, radius=R)
    H_grid = np.insert(H_grid, 0, 0, axis=0)
    H_grid = miepy.coordinates.vec_sph_to_cart(H_grid, THETA, PHI)

    H = np.average(H_grid, axis=(1,2,3))
    divH = div(H_grid, eps)
    curlH = curl(H_grid, eps)

    assert np.abs(divE/(k*np.linalg.norm(E))) < 1e-6, 'div(E) = 0'
    assert np.abs(divH/(k*np.linalg.norm(H))) < 1e-6, 'div(H) = 0'
    assert np.allclose(curlE[:2], 1j*k*H[:2], atol=0, rtol=1e-5), 'curl(E) = ikH'
    assert np.allclose(curlH[:2], -1j*k*E[:2], atol=0, rtol=1e-5), 'curl(H) = -ikE'

@pytest.mark.skip(reason="not implemented")
def test_maxwells_equation_interior():
    """
    Verify Maxwell's equations in the interior of a particle
    """
    pass

@pytest.mark.skip(reason="not implemented")
def test_maxwells_equation_in_medium():
    """
    Verify Maxwell's equations in a non-vacuum medium
    """
    pass
