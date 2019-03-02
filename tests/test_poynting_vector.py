"""
Verify the direction of the Poynting vector in the near and far field case for sources
"""

import numpy as np
import miepy
from miepy.constants import Z0
import pytest

nm = 1e-9
wav = 600*nm
k = 2*np.pi/wav
width = 400*nm


@pytest.mark.parametrize("source", [
    (miepy.sources.plane_wave(polarization=[1,0])),
    (miepy.sources.gaussian_beam(width=width, polarization=[1,0])),
])
def test_near_field_poytning_vector_direction_at_origin(source):
    """
    S should be in the +z direction at the origin
    """
    E = source.E_field(0, 0, 0, k)
    H = source.H_field(0, 0, 0, k)
    S = np.real(np.cross(E, np.conj(H))/(2*Z0))

    assert np.allclose(np.linalg.norm(S), S[2], rtol=1e-20)


def test_plane_wave_poynting_vector_non_normal_incidence():
    """
    S should be in the direction [1, 1, sqrt(2)] for (theta = pi/4, phi = pi/4) plane-wave
    """
    source = miepy.sources.plane_wave(polarization=[1,2], theta=np.pi/4, phi=np.pi/4)
    E = source.E_field(200*nm, 50*nm, 100*nm, k)
    H = source.H_field(200*nm, 50*nm, 100*nm, k)
    S = np.real(np.cross(E, np.conj(H))/(2*Z0))

    assert np.allclose(S[0], S[1], atol=0, rtol=1e-12), 'Sx should equal Sy'
    assert np.allclose(np.sqrt(2)*S[0], S[2], atol=0, rtol=1e-12), 'Sx should equal Sz/sqrt(2)'


def test_gaussian_poynting_vector_far_field():
    """
    S should be in the direction of ±r_hat (+ for upper hemisphere, - for lower hemisphere) in the angular far fields
    """
    source = miepy.sources.gaussian_beam(width=width, polarization=[1,2j])

    theta = np.linspace(0, np.pi, 6)
    phi = np.linspace(0, 2*np.pi, 6)
    THETA, PHI = np.meshgrid(theta, phi)

    E = source.E_angular(THETA, PHI, k)
    E = np.insert(E, 0, 0, axis=0)
    E = miepy.coordinates.vec_sph_to_cart(E, THETA, PHI)

    H = source.H_angular(THETA, PHI, k)
    H = np.insert(H, 0, 0, axis=0)
    H = miepy.coordinates.vec_sph_to_cart(H, THETA, PHI)

    r_hat, *_ = miepy.coordinates.sph_basis_vectors(THETA, PHI)
    r_hat[:,THETA > np.pi/2] *= -1

    S = np.real(np.cross(E, np.conj(H), axis=0)/(2*Z0))
    S = S/np.linalg.norm(S, axis=0)

    assert np.allclose(S, r_hat, atol=1e-24, rtol=1e-15)

