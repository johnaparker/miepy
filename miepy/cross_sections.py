"""
Expressions for the cross-sections given the expansion coefficients or the fields
"""

import numpy as np

#TODO move elsewhere
def flux_from_particle(self, i, source=False):
    """Determine the scattered flux from a single particle

        Arguments:
            i         Particle index
            source    Include the source field (bool, default=False)
        
        Returns: flux[M], M = number of wavelengths
    """
    r = self.spheres.radius[i]
    X,Y,Z,THETA,PHI,tau,phi = discrete_sphere(r, self.Ntheta, self.Nphi, self.spheres.position[i])
    rhat,*_ = sph_unit_vectors(THETA, PHI)

    E_all = self.E_field_from_particle(i, X, Y, Z, source)
    H_all = self.H_field_from_particle(i, X, Y, Z, source)

    flux = np.zeros(self.Nfreq, dtype=float)

    for k in range(self.Nfreq):
        E = E_all[:,k]
        H = H_all[:,k]

        S = np.real(np.einsum('ijk,jxy,kxy->ixy', levi, E, np.conj(H)))
        dA = r**2

        integrand = np.einsum('ixy,ixy->xy', S, rhat)*dA
        flux[k] = simps_2d(tau, phi, integrand)

    return flux

#TODO move elsewhere
def flux(self, source=False):
    """Determine the scattered flux from every particle

        Arguments:
            source    Include the source field (bool, default=False)
        
        Returns: flux[M,N], N = number of particle, M = number of wavelengths
    """
    flux_data = np.zeros([self.Nfreq, self.Nparticles], dtype=float)
    for i in range(self.Nparticles):
        flux_data[:,i] = self.flux_from_particle(i, source=source)

    return flux_data

def discrete_sphere(radius, Ntheta, Nphi, center=None):
    """Given a radius, center, and theta phi sampling, return the 
    X,Y,Z,THETA,PHI,tau,phi coordinates of a discretized sphere""" 
    if center is None:
        center = np.zeros(3)

    r = np.array([radius])
    tau = np.linspace(-1,1, Ntheta) 
    theta = np.pi - np.arccos(tau)
    phi = np.linspace(0, 2*np.pi, Nphi)
    R, THETA, PHI = np.meshgrid(r,theta,phi, indexing='ij')

    X = center[0] + R*np.sin(THETA)*np.cos(PHI)
    Y = center[1] + R*np.sin(THETA)*np.sin(PHI) 
    Z = center[2] + R*np.cos(THETA)

    return map(np.squeeze, (X,Y,Z,THETA,PHI,tau,phi))

