"""
Compare the T-matrics of spheroids to ellipsoids
"""

import miepy
import numpy as np

nm = 1e-9
rx = ry = 40*nm
rz = 70*nm
material = miepy.materials.Ag()
wavelength = 700*nm
lmax = 3

def test_ellipsoid_equal_spheroid_z_oriented():
    """An (rx, ry, rz) ellipsoid should equal a z-oriented (rxy, rz) spheroid"""
    e = miepy.spheroid([0,0,0], rx, rz, material)
    T1 = e.compute_tmatrix(lmax, wavelength, 1.0)

    e = miepy.ellipsoid([0,0,0], rx, ry, rz, material)
    T2 = e.compute_tmatrix(lmax, wavelength, 1.0)

    assert np.allclose(T1, T2, atol=2e-5)

def test_ellipsoid_equal_spheroid_x_oriented():
    """An (rz, rx, ry) ellipsoid should equal an x-oriented (rxy, rz) spheroid"""
    q = miepy.quaternion.from_spherical_coords(np.pi/2, 0)
    e = miepy.spheroid([0,0,0], rx, rz, material, orientation=q)
    T1 = e.compute_tmatrix(lmax, wavelength, 1.0)

    e = miepy.ellipsoid([0,0,0], rz, rx, ry, material)
    T2 = e.compute_tmatrix(lmax, wavelength, 1.0)

    assert np.allclose(T1, T2, atol=5e-5)

def test_ellipsoid_equal_spheroid_y_oriented():
    """An (rx, rz, ry) ellipsoid should equal a y-oriented (rxy, rz) spheroid"""
    q = miepy.quaternion.from_spherical_coords(np.pi/2, np.pi/2)
    e = miepy.spheroid([0,0,0], rx, rz, material, orientation=q)
    T1 = e.compute_tmatrix(lmax, wavelength, 1.0)

    e = miepy.ellipsoid([0,0,0], rx, rz, ry, material)
    T2 = e.compute_tmatrix(lmax, wavelength, 1.0)

    assert np.allclose(T1, T2, atol=5e-5)
