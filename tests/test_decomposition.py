import numpy as np

import miepy


def test_vsh_source_decomposition():
    """Verify the decomposition of the vsh_source."""
    x, y, z = (-0.1, 0.1, 0.1)
    k = 2 * np.pi
    lmax = 7
    p = np.array([[0.1, 0, -0.1]])
    R, THETA, PHI = miepy.coordinates.cart_to_sph(x - p[0, 0], y - p[0, 1], z - p[0, 2])

    source = miepy.sources.vsh_source(2, 2, ftype="electric", center=(0.3, 0.2, 0.1))
    E1 = source.E_field(x, y, z, k)
    p_src = source.structure(p, k, lmax)[0]
    Efunc = miepy.vsh.expand_E(p_src, k, miepy.vsh_mode.incident)
    E2 = Efunc(R, THETA, PHI)
    E2 = miepy.coordinates.vec_sph_to_cart(E2, THETA, PHI)

    assert np.allclose(E1[:2], E2[:2], atol=0, rtol=1e-7), "x,y components equal (electric mode)"
    assert np.allclose(E1[2], 0, atol=1e-15), "z components of E1 goes to zero (electric mode)"
    assert np.allclose(E2[2], 0, atol=1e-7), "z component of E2 goes to zero (electric mode)"

    source = miepy.sources.vsh_source(2, 2, ftype="magnetic", center=(0.3, 0.2, 0.1))
    E1 = source.E_field(x, y, z, k)
    p_src = source.structure(p, k, lmax)[0]
    Efunc = miepy.vsh.expand_E(p_src, k, miepy.vsh_mode.incident)
    E2 = Efunc(R, THETA, PHI)
    E2 = miepy.coordinates.vec_sph_to_cart(E2, THETA, PHI)

    assert np.allclose(E1[2], E2[2], atol=0, rtol=1e-7), "z components equal (magnetic mode)"
    assert np.allclose(E1[:2], 0, atol=1e-15), "x,yz components of E1 go to zero (magnetic mode)"
    assert np.allclose(E2[:2], 0, atol=1e-7), "x,y component of E2 go to zero (magnetic mode)"
