"""
tests for far-fields
"""
import numpy as np
import miepy

nm = 1e-9

# dimer = miepy.cluster(position=[[-100*nm,0,0], [100*nm, 0, 0]],
dimer = miepy.cluster(position=[0,0,0],
                      radius=75*nm,
                      material=miepy.constant_material(3.6**2),
                      source=miepy.sources.y_polarized_plane_wave(),
                      wavelength=600*nm,
                      Lmax=2)

theta = np.linspace(0, np.pi, 5)
phi = np.linspace(0, 2*np.pi, 5)
THETA, PHI = np.meshgrid(theta, phi)
R = 1e6
X, Y, Z = miepy.coordinates.sph_to_cart(R, THETA, PHI)
E_exact = dimer.E_field(R, THETA, PHI, source=False, interior=False, spherical=True)
# E_exact = dimer.E_field(X, Y, Z, source=False, interior=False)

def test_far_field_convergence():
    """far-field E and H field should agree with exact field in the large radius limit"""
    E_far = dimer.E_field(R, THETA, PHI, far=True, source=False, interior=False, spherical=True)
    # E_far = dimer.E_field(X, Y, Z, far=True, source=False, interior=False)

    print(E_exact[:,1,3])
    print(E_far[:,1,3])
    print(E_exact[:,4,3])
    print(E_far[:,4,3])

    assert False

def test_far_field_cluster_coefficient():
    """far-fields calculated from the cluster coefficients should be the same as the sum-over particle coefficients"""
    dimer.solve_cluster_coefficients(Lmax=4)
    E_func = miepy.vsh.expand_E_far(dimer.p_cluster, dimer.material_data.k)
    E_far = E_func(R, THETA, PHI)

    print(np.abs(E_exact[:,2,3]))
    print(np.abs(E_far[:,2,3]))

    assert False

if __name__ == "__main__":
    test_far_field_convergence()
    # test_far_field_cluster_coefficient()
