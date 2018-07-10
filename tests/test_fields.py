"""
Tests for E and H field calculations
"""
import numpy as np
import miepy

nm = 1e-9

def test_boundary_conditions():
    """verifies the continunity of tangential components of E and H at the surface of a particle"""
    radius = 50*nm

    cluster = miepy.sphere_cluster(position=[[0,0,0]],
                                   radius=radius,
                                   material=miepy.materials. Ag(),
                                   lmax=2,
                                   wavelength=600*nm,
                                   source=miepy.sources.plane_wave.from_string(polarization='y'),
                                   medium=miepy.constant_material(1.2**2))

    theta = 0.3
    phi = 0.3
    eps = .1*nm

    E_out = cluster.E_field(radius + eps, theta, phi, spherical=True)
    E_in = cluster.E_field(radius - eps, theta, phi, spherical=True)

    H_out = cluster.H_field(radius + eps, theta, phi, spherical=True)
    H_in = cluster.H_field(radius - eps, theta, phi, spherical=True)

    assert np.allclose(E_out[1:], E_in[1:], atol=4e-2, rtol=0)
    assert np.allclose(H_out[1:], H_in[1:], atol=4e-2, rtol=0)
