"""
Defines functions used to construct basis vectors, convert between coordinate systems, build meshes, etc.
"""

import numpy as np

def sph_to_cart(r, theta, phi, origin=None):
    """convert spherical coordinates (r, theta, phi) centered at origin to cartesian coordinates (x, y, z)"""
    if origin is None:
        origin = np.zeros(3, dtype=float)

    x = origin[0] + r*np.sin(theta)*np.cos(phi)
    y = origin[1] + r*np.sin(theta)*np.sin(phi)
    z = origin[2] + r*np.cos(theta)

    return x,y,z

def cart_to_sph(x, y, z, origin=None):
    """convert cartesian coordinates (x, y, z) to spherical coordinates (r, theta, phi) centered at origin"""
    if origin is None:
        origin = np.zeros(3, dtype=float)

    x0,y0,z0 = origin
    r = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2)**0.5
    theta = np.arccos((z - z0)/r)
    phi = np.arctan2(y - y0, x - x0)

    return r, theta, phi

#TODO: implement origin
def sph_basis_vectors(theta, phi, origin=None):
    """obtain the spherical basis vectors (r_hat, theta_hat, phi_hat) for given theta, phi"""
    if origin is None:
        origin = np.zeros(3, dtype=float)

    r_hat = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
    theta_hat = np.array([np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -1*np.sin(theta)])
    phi_hat = np.array([-1*np.sin(phi), np.cos(phi), np.zeros_like(phi)])

    return r_hat, theta_hat, phi_hat

#TODO: implement origin
def vec_cart_to_sph(F, theta, phi, origin=None):
    """convert a vector field F from cartesian to spherical coordinates

    Arguments:
        F[3,...]     vector field values
        theta        theta coordinates
        phi          phi coordinates
    """
    if origin is None:
        origin = np.zeros(3, dtype=float)

    Fsph = np.zeros_like(F)
    r_hat, theta_hat, phi_hat = sph_basis_vectors(theta, phi)
    Fsph[0] = np.sum(F*r_hat, axis=0)
    Fsph[1] = np.sum(F*theta_hat, axis=0)
    Fsph[2] = np.sum(F*phi_hat, axis=0)

    return Fsph

#TODO: implement origin
def vec_sph_to_cart(F, theta, phi, origin=None):
    """convert a vector field F from spherical to cartesian coordinates

    Arguments:
        F[3,...]     vector field values
        theta        theta coordinates
        phi          phi coordinates
    """
    if origin is None:
        origin = np.zeros(3, dtype=float)

    Fcart = np.zeros_like(F)
    r_hat, theta_hat, phi_hat = sph_basis_vectors(theta, phi)
    for i in range(3):
        Fcart[i] = F[0]*r_hat[i] + F[1]*theta_hat[i] + F[2]*phi_hat[i]

    return Fcart

def sphere_mesh(sampling):
    """
    Obtain a THETA,PHI mesh for discretizing the surface of the sphere, consistent
    with the format required by the project and decompose functions
    Returns (THETA,PHI) meshgrids

    Arguments:
        sampling   number of points to sample between 0 and pi
    """

    phi = np.linspace(0, 2*np.pi, 2*sampling)
    tau = np.linspace(-1, 1, sampling)
    theta = np.arccos(tau)

    THETA,PHI = np.meshgrid(theta, phi, indexing='ij')
    return THETA, PHI

def cart_sphere_mesh(radius, origin, sampling):
    """Given a radius, origin and sampling, return the
    X,Y,Z,THETA,PHI,tau,phi coordinates of a discretized sphere""" 

    r = np.array([radius])
    tau = np.linspace(-1,1, sampling) 
    theta = np.arccos(tau)
    phi = np.linspace(0, 2*np.pi, 2*sampling)
    R, THETA, PHI = np.meshgrid(r,theta,phi, indexing='ij')

    X = origin[0] + R*np.sin(THETA)*np.cos(PHI)
    Y = origin[1] + R*np.sin(THETA)*np.sin(PHI) 
    Z = origin[2] + R*np.cos(THETA)

    return map(np.squeeze, (X,Y,Z,THETA,PHI,tau,phi))
