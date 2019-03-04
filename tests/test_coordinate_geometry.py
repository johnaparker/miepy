"""
Tests for coordinate transformations and related geometric functions
"""

import miepy
import numpy as np


class Test_spherical_coords_rotation:
    """
    Verify rotation of spherical coordinates (theta, phi) with `miepy.coordinates.rotate_sph`
    """

    def test_compare_to_cartestion_rotation(self):
        """
        Spherical rotation is identical to cartestion rotation after coordinate transformations
        """
        theta = .5
        phi = 1.2
        r = 1

        q_rot = miepy.quaternion.from_spherical_coords(.2, .1)

        x, y, z = miepy.coordinates.sph_to_cart(r, theta, phi)
        xr, yr, zr = miepy.coordinates.rotate(x, y, z, q_rot)
        _, theta_r1, phi_r1 = miepy.coordinates.cart_to_sph(xr, yr, zr)

        theta_r2, phi_r2 = miepy.coordinates.rotate_sph(theta, phi, q_rot)

        assert np.allclose(theta_r1, theta_r2), 'theta components are equal'
        assert np.allclose(phi_r1, phi_r2), 'phi components are equal'

    def test_theta_equal_to_zero(self):
        """
        When theta is equal to zero, rotate_sph should map phi -> phi + phi_rot correctly
        """

        theta = 0
        phi = 1.2
        r = 1
        phi_rot = .4

        q_rot = miepy.quaternion.from_spherical_coords(0, phi_rot)
        theta_r, phi_r = miepy.coordinates.rotate_sph(theta, phi, q_rot)

        assert np.allclose(theta, theta_r), 'theta components are equal'
        assert np.allclose(phi + phi_rot, phi_r), 'phi components are rotated correctly'
