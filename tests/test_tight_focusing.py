"""
Test tighly focused beams by comparing the E/H field from two methods:
    (i) Directly integrating the far-field angular spectrum to obtain focal fields
    (ii) Using the expansion coefficients around the center of the beam
"""

import numpy as np
import miepy
import pytest

nm = 1e-9
wav = 600*nm
k = 2*np.pi/wav

### grid points to evaluate field over
x = np.linspace(-200*nm, 200*nm, 2)
y = np.linspace(-200*nm, 200*nm, 2)
z = np.linspace(-200*nm, 200*nm, 2)
X, Y, Z = np.meshgrid(x, y, z)

width = 200*nm
polarization = [1,0]

@pytest.mark.parametrize("source,rtol", [
    (miepy.sources.gaussian_beam(width=width, polarization=polarization), 2e-4),
    (miepy.sources.hermite_gaussian_beam(1, 0, width=width, polarization=polarization), 4e-4),
    (miepy.sources.laguerre_gaussian_beam(1, 1, width=width, polarization=polarization), 9e-3),
    (miepy.sources.azimuthal_beam(width=width), 3e-4),
    (miepy.sources.bigaussian_beam(width_x=1.5*width, width_y=width/1.5, polarization=polarization), 4e-4),
])
def test_E_field_tight_focusing(source, rtol):
    E1 = source.E_field(X, Y, Z, k)

    lmax = 8
    p_src = source.structure([0, 0, 0], k, lmax)
    Efunc = miepy.vsh.expand_E(p_src, k, miepy.vsh_mode.incident)
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)
    E2 = Efunc(R, THETA, PHI)
    E2 = miepy.coordinates.vec_sph_to_cart(E2, THETA, PHI)

    assert np.allclose(E1, E2, rtol=rtol, atol=1e-10)

@pytest.mark.parametrize("source,rtol", [
    (miepy.sources.gaussian_beam(width=width, polarization=polarization), 2e-4),
    (miepy.sources.hermite_gaussian_beam(1, 0, width=width, polarization=polarization), 4e-4),
    (miepy.sources.laguerre_gaussian_beam(1, 1, width=width, polarization=polarization), 9e-3),
    (miepy.sources.azimuthal_beam(width=width), 3e-4),
    (miepy.sources.bigaussian_beam(width_x=1.5*width, width_y=width/1.5, polarization=polarization), 4e-4),
])
def test_H_field_tight_focusing(source, rtol):
    H1 = source.H_field(X, Y, Z, k)

    lmax = 8
    p_src = source.structure([0, 0, 0], k, lmax)
    Hfunc = miepy.vsh.expand_H(p_src, k, miepy.vsh_mode.incident, 1, 1)
    R, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)
    H2 = Hfunc(R, THETA, PHI)
    H2 = miepy.coordinates.vec_sph_to_cart(H2, THETA, PHI)

    assert np.allclose(H1, H2, rtol=rtol, atol=1e-10)
