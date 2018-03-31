"""
Functions related to the flux: Poynting vector, cross-sections, etc.
"""

import numpy as np
from scipy import constants
import miepy
from my_pytools.my_numpy.integrate import simps_2d

#TODO eps/mu role here (related to our definition of the H field, eps/mu factor)
#TODO factor of 1/2 for complex fields not present???
def poynting_vector(E, H, eps=1, mu=1):
    """Compute the Poynting vector
    
       Arguments:
           E[3,...]   electric field data
           H[3,...]   magnetic field data
           eps        medium permitvitty (default: 1)
           mu         medium permeability (default: 1)

       Returns S[3,...]
    """

    S = np.cross(E, np.conj(H), axis=0)
    n_b = np.sqrt(eps*mu)

    return np.real(S)/n_b

#TODO is this right, can it be more useful?
def flux_from_poynting(E, H, Ahat, eps=1, mu=1):
    """Compute the flux from the E and H field over some area using the Poynting vector

       Arguments:
           E[3,...]             electric field values on some surface
           H[3,...]             magnetic field values on some surface
           Ahat[3,...]          normal vectors of the surface
           eps                  medium permitvitty (default: 1)
           mu                   medium permeability (default: 1)

       Returns flux (scalar)
    """
    S = poynting_vector(E, H, eps, mu)
    integrand = np.einsum('i...,i...->...', S, Ahat)

    return np.sum(integrand)

def flux_from_poynting_sphere(E, H, radius, eps=1, mu=1):
    """Compute the flux from the E and H field on the surface of a sphere using the Poynting vector

       Arguments:
           E[3,Ntheta,Nphi]     electric field values on the surface of a sphere
           H[3,Ntheta,Nphi]     magnetic field values on the surface of a sphere
           radius               radius of sphere 
           eps                  medium permitvitty (default: 1)
           mu                   medium permeability (default: 1)

       Returns flux (scalar)
    """
    S = poynting_vector(E, H, eps, mu)

    Ntheta, Nphi = E.shape[1:]
    THETA, PHI = miepy.coordinates.sphere_mesh(Ntheta)
    rhat,*_ = miepy.coordinates.sph_basis_vectors(THETA, PHI)

    tau = np.linspace(-1, 1, Ntheta)
    phi = np.linspace(0, 2*np.pi, Nphi)
    dA = radius**2

    integrand = np.einsum('ixy,ixy->xy', S, rhat)*dA
    flux = simps_2d(tau, phi, integrand)

    return flux

def _gmt_cross_sections_from_poynting(gmt, radius, sampling=30):
    """FOR TESTING ONLY!
    Given GMT object and particle number i, return cross-sections (C,A,E) from poynting vector
    """
    X,Y,Z,THETA,PHI,tau,phi = miepy.coordinates.cart_sphere_mesh(radius, gmt.origin, sampling)

    E_tot = gmt.E_field(X, Y, Z)
    H_tot = gmt.H_field(X, Y, Z)

    E_scat = gmt.E_field(X, Y, Z, source=False)
    H_scat = gmt.H_field(X, Y, Z, source=False)

    eps_b = gmt.material_data.eps_b
    mu_b = gmt.material_data.mu_b
    C = flux_from_poynting_sphere(E_scat, H_scat, radius, eps_b, mu_b)
    A = -flux_from_poynting_sphere(E_tot, H_tot, radius, eps_b, mu_b)

    return C, A, C+A
