"""
Tests for the power quantity of different incident beams
"""

import numpy as np
import miepy
from miepy.vsh.misc import trapz_2d
from miepy.constants import Z0
import pytest

nm = 1e-9
wav = 600*nm
k = 2*np.pi/wav

width = 800*nm
power = .2
polarization = [0,1]

xmax = 2000*nm
x = np.linspace(-xmax, xmax, 21)
y = np.linspace(-xmax, xmax, 21)
X, Y = np.meshgrid(x, y)
Z = np.zeros_like(X)

@pytest.mark.parametrize("source,rtol", [
    (miepy.sources.gaussian_beam(width, polarization, power=power), 2e-2),
    (miepy.sources.hermite_gaussian_beam(1, 0, width=width, polarization=polarization, power=power), 4e-4),
    (miepy.sources.laguerre_gaussian_beam(1, 1, width=width, polarization=polarization, power=power), 8e-3),
    (miepy.sources.azimuthal_beam(width=width, power=power), 3e-5),
    (miepy.sources.bigaussian_beam(width_x=width, width_y=width/2, polarization=polarization, power=power), 8e-3),
])
def test_power_by_near_field_poynting_vector(source, rtol):
    """Check that the power is correct by integrating the near field Poynting vector in the xy plane"""
    E = source.E_field(X, Y, Z, k)
    H = source.H_field(X, Y, Z, k)
    S = 0.5/Z0*np.cross(E, np.conjugate(H), axis=0)[2]
    P = trapz_2d(x, y, S).real

    assert np.allclose(P, power, rtol=rtol)

def test_power_by_far_field_poynting_vector_upper_hemisphere():
    """Check that the power is correct by integrating the far field Poynting vector in the forward hemisphere plane
    Expect the power to be positive when the unit vectors are parallel to the outgoing r_hat"""

    source = miepy.sources.gaussian_beam(width=100*nm, polarization=[1,0], power=power)

    theta = np.linspace(0, np.pi/2 - 1e-5, 30)
    phi = np.linspace(0, 2*np.pi, 60)
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    RAD = 1e6*wav*np.ones_like(THETA)

    # S = np.sum(np.abs(E)**2, axis=0)/(2*Z0)

    E = source.E_angular(THETA, PHI, k, RAD)
    E = np.insert(E, 0, 0, axis=0)
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
    H = source.H_angular(THETA, PHI, k, RAD)
    H = np.insert(H, 0, 0, axis=0)
    H = miepy.coordinates.vec_sph_to_cart(H, THETA, PHI)

    S = np.cross(E, np.conjugate(H), axis=0)/(2*Z0)

    rhat, *_ = miepy.coordinates.sph_basis_vectors(THETA, PHI)

    integrand = np.sum(rhat*S, axis=0)*np.sin(THETA)*RAD**2
    P = miepy.vsh.misc.trapz_2d(theta, phi, integrand).real

    assert np.allclose(P, power, rtol=8e-4)

def test_power_by_far_field_poynting_vector_lower_hemisphere():
    """Check that the power is correct by integrating the far field Poynting vector in the backward hemisphere.
    Expect the power to be negative when the unit vectors are parallel to the outgoing r_hat"""

    source = miepy.sources.gaussian_beam(width=100*nm, polarization=[1,0], power=power)

    theta = np.linspace(np.pi/2 + 1e-5, np.pi, 30)
    phi = np.linspace(0, 2*np.pi, 60)
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    RAD = 1e6*wav*np.ones_like(THETA)

    # S = np.sum(np.abs(E)**2, axis=0)/(2*Z0)

    E = source.E_angular(THETA, PHI, k, RAD)
    E = np.insert(E, 0, 0, axis=0)
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
    H = source.H_angular(THETA, PHI, k, RAD)
    H = np.insert(H, 0, 0, axis=0)
    H = miepy.coordinates.vec_sph_to_cart(H, THETA, PHI)

    S = np.cross(E, np.conjugate(H), axis=0)/(2*Z0)

    rhat, *_ = miepy.coordinates.sph_basis_vectors(THETA, PHI)

    integrand = np.sum(rhat*S, axis=0)*np.sin(THETA)*RAD**2
    P = miepy.vsh.misc.trapz_2d(theta, phi, integrand).real

    assert np.allclose(P, -power, rtol=8e-4)

class Test_gaussian_beam_numeric_power:
    """compute power of a Gaussian beam using various methods"""
    width = 400*nm
    power = 1
    polarization = [1,0]

    source = miepy.sources.gaussian_beam(width, polarization, power=power)

    def test_power_through_aperature(self):
        """power by using power_through_aperature function for a large aperature"""
        P = miepy.vsh.misc.power_through_aperature(self.source, [0,0,0], 3*self.width, k, sampling=20)
        assert np.allclose(P, self.power, rtol=.03)

    def test_power_spherically_ingoing_waves(self):
        """power by far-field integration of spherically ingoing waves"""
        lmax = 6
        radius = 1e6*(2*np.pi/k)

        p_src = self.source.structure([[0,0,0]], k, lmax)[0]
        theta = np.linspace(np.pi/2, np.pi, 60)
        phi = np.linspace(0, 2*np.pi, 60)
        THETA, PHI = np.meshgrid(theta, phi)
        R = radius*np.ones_like(THETA)

        Efunc = miepy.expand_E(p_src/2, k, miepy.vsh_mode.ingoing)
        E = Efunc(R, THETA, PHI)
        S = 0.5/Z0*np.sum(np.abs(E)**2, axis=0)*np.sin(THETA)
        P = radius**2*trapz_2d(theta, phi, S.T).real

        print(P)
        assert np.allclose(P, self.power, rtol=.04)

    def test_power_p_src_analytic(self):
        """power by analytic sum over p_src coefficients"""
        lmax = 6
        p_src = self.source.structure([[0,0,0]], k, lmax)[0]
        factor = 0.5/Z0*np.pi/k**2
        P = factor*np.sum(np.abs(p_src)**2)

        assert np.allclose(P, self.power, rtol=.04)

    def test_power_2d_integration(self):
        """power by  integrating Poynting vector in 2D plane"""
        lmax = 6
        p_src = self.source.structure([[0,0,0]], k, lmax)[0]
        a = 3000*nm
        x = np.linspace(-a, a, 20)
        y = np.linspace(-a, a, 20)
        X, Y = np.meshgrid(x, y)
        R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, 0)
        Efunc = miepy.expand_E(p_src, k, miepy.vsh_mode.incident)
        Hfunc = miepy.expand_H(p_src, k, miepy.vsh_mode.incident, 1, 1)
        E = Efunc(R, THETA, PHI)
        H = Hfunc(R, THETA, PHI)
        E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)
        H = miepy.coordinates.vec_sph_to_cart(H, THETA, PHI)/Z0

        S = 0.5*np.linalg.norm(np.cross(E, np.conjugate(H), axis=0), axis=0)
        S = 0.5*np.cross(E, np.conjugate(H), axis=0)[2]
        # S = 0.5*np.sum(np.abs(E)**2, axis=0)

        P = trapz_2d(x, y, S).real

        assert np.allclose(P, self.power, rtol=.04)
