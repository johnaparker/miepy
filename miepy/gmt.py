"""
The Generalized Mie Theory (GMT) for a collection of spheres.
"""
import numpy as np
import miepy
from my_pytools.my_numpy.integrate import simps_2d
from my_pytools.my_numpy.indices import levi_civita
from my_pytools.my_numpy.array import atleast
from collections import namedtuple

levi = levi_civita()

def cart_to_sph(x, y, z, center=None):
    """convert the cartesian coordinates (x,y,z) to sphereical coordinats around center point"""
    if center is None:
        center = np.zeros(3)
    Xs = x - center[0]
    Ys = y - center[1]
    Zs = z - center[2]

    R = np.sqrt(Xs**2 + Ys**2 + Zs**2)
    THETA = np.arccos(Zs/R)
    PHI = np.arctan2(Ys,Xs)
    return R, THETA, PHI

def sph_unit_vectors(THETA, PHI):
    """return the spherical coordinate unit vectors (rhat,that,phat) given angular coordinats"""
    rhat = np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])
    that = np.array([np.cos(THETA)*np.cos(PHI), np.cos(THETA)*np.sin(PHI), -np.sin(THETA)])
    phat = np.array([-np.sin(PHI), np.cos(PHI), np.zeros_like(THETA)])
    return rhat, that, phat

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

#TODO swap all indices, so that [N,3] => [3,N]
#TODO make position a property so that it can be set properly (if input is a list)
class spheres:
    """A collection of N spheres"""
    sphere_type = namedtuple('sphere', ['position', 'radius', 'material'])

    def __init__(self, position, radius, material):
        """Arguments:
                position[N,3] or [3]        sphere positions
                radius[N] or scalar         sphere radii
                material[N] or scalar       sphere materials
        """
        self.position = np.asarray(np.atleast_2d(position), dtype=float)
        self.radius = atleast(radius, dim=1, length=self.position.shape[0], dtype=float)
        self.material = atleast(material, dim=1, length=self.position.shape[0], dtype=np.object)

        if (self.position.shape[0] != self.radius.shape[0] != self.material.shape[0]):
            raise ValueError("The shapes of position, radius, and material do not match")

    def __len__(self):
        return self.position.shape[0]

    def __iter__(self):
        return (self.sphere_type(position=self.position[i], radius=self.radius[i],
                            material=self.material[i]) for i in range(len(self)))

    def __getitem__(self, i):
        return self.sphere_type(position=self.position[i], radius=self.radius[i],
                            material=self.material[i])



class gmt:
    """Solve Generalized Mie Theory: N particles in an arbitray source profile"""
    def __init__(self, spheres, source, wavelength, Lmax, medium=None, interactions=True,
                   Ntheta=51, Nphi=31):
        """Arguments:
               spheres          spheres object specifying the positions, radii, and materials
               source           source object specifying the incident E and H functions
               wavelength[M]    wavelength(s) to solve the system at
               Lmax             maximum number of orders to use in angular momentum expansion (int)
               medium           (optional) material medium (must be non-absorbing; default=vacuum)
               interactions     (optional) If True, include particle interactions (bool, default=True) 
               Ntheta           (optional) number of points in theta to use in force/flux calculations (default = 51)
               Nphi             (optional) number of points in phi to use in force/flux calculations (default = 31)
        """
        self.spheres = spheres
        self.source = source
        self.wavelength = np.asarray(np.atleast_1d(wavelength), dtype=float)
        self.Lmax = Lmax
        self.interactions = interactions

        self.Ntheta = Ntheta
        self.Nphi = Nphi

        if medium is None:
            self.medium = miepy.constant_material(1.0, 1.0)
        else:
            self.medium = medium
            if (self.medium.eps(self.wavelength).imag != 0).any()  \
                         or (self.medium.mu(self.wavelength).imag != 0).any():
                raise ValueError('medium must be non-absorbing')

        self.Nfreq = len(self.wavelength)
        self.Nparticles = len(spheres)

        self.material_data = {}
        self.material_data['wavelength'] = self.wavelength
        self.material_data['eps']        = np.zeros([self.Nparticles, self.Nfreq], dtype=complex) 
        self.material_data['mu']         = np.zeros([self.Nparticles, self.Nfreq], dtype=complex) 
        self.material_data['n']          = np.zeros([self.Nparticles, self.Nfreq], dtype=complex) 
        self.material_data['eps_b']      = self.medium.eps(self.wavelength)
        self.material_data['mu_b']       = self.medium.mu(self.wavelength)
        self.material_data['n_b']        = np.sqrt(self.material_data['eps_b']*self.material_data['mu_b'])
        self.material_data['k']          = 2*np.pi*self.material_data['n_b']/self.wavelength

        self.a = np.zeros([self.Nparticles,self.Nfreq,self.Lmax], dtype=complex)
        self.b = np.zeros([self.Nparticles,self.Nfreq,self.Lmax], dtype=complex)
        self.p = np.zeros([2,self.Nparticles,self.Nfreq], dtype=complex)
        self.q = np.zeros([2,self.Nparticles,self.Nfreq], dtype=complex)

        for i in range(self.Nparticles):
            sphere = miepy.single_mie_sphere(self.spheres.radius[i], self.spheres.material[i],
                        self.wavelength, self.Lmax, self.medium)
            self.material_data['eps'][i] = sphere.material_data['eps']
            self.material_data['mu'][i] = sphere.material_data['mu']
            self.material_data['n'][i] = sphere.material_data['n']
            self.a[i], self.b[i] = sphere.solve_exterior()

        if (self.interactions):
            self._solve_interactions()
        else:
            self._set_without_interactions()
    
    def E_field(self, x, y, z, inc=True):
        """Compute the electric field due to all particles
             
            Arguments:
                x      x position (array-like) 
                y      y position (array-like) 
                z      z position (array-like) 
                inc    Include the incident field (bool, default=True)

            Returns: E[3,M], M = number of wavelengths
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        E = np.zeros((3,self.Nfreq) + x.shape, dtype=complex)
        for i in range(self.Nparticles):
            R, THETA, PHI = cart_to_sph(x, y, z, self.spheres.position[i])
            rhat, that, phat = sph_unit_vectors(THETA, PHI)

            for k in range(self.Nfreq):
                Ax, Ay = self.p[:,i,k]
                E_func = miepy.scattering.scattered_E(self.a[i,k], self.b[i,k], self.material_data['k'][k])
                E_sph = Ax*E_func(R,THETA,PHI) + Ay*E_func(R,THETA,PHI-np.pi/2)
                E[:,k] += E_sph[0]*rhat + E_sph[1]*that + E_sph[2]*phat      # convert to cartesian

        if inc:
            for k in range(self.Nfreq):
                E[:,k] += self.source.E(np.array([x,y,z]), self.material_data['k'][k])
        
        return E

    def H_field(self, x, y, z, inc=True):
        """Compute the magnetic field due to all particles
             
            Arguments:
                x      x position (array-like) 
                y      y position (array-like) 
                z      z position (array-like) 
                inc    Include the incident field (bool, default=True)

            Returns: H[3,M], M = number of wavelengths
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        H = np.zeros((3,self.Nfreq) + x.shape, dtype=complex)
        for i in range(self.Nparticles):
            R, THETA, PHI = cart_to_sph(x, y, z, self.spheres.position[i])
            rhat, that, phat = sph_unit_vectors(THETA, PHI)

            for k in range(self.Nfreq):
                Ax, Ay = self.p[:,i,k]
                H_func = miepy.scattering.scattered_H(self.a[i,k], self.b[i,k], self.material_data['k'][k],
                              self.material_data['n_b'][k], self.material_data['mu_b'][k])
                H_sph = Ax*H_func(R,THETA,PHI) + Ay*H_func(R,THETA,PHI-np.pi/2)
                H[:,k] += H_sph[0]*rhat + H_sph[1]*that + H_sph[2]*phat      # convert to cartesian

        if inc:
            for k in range(self.Nfreq):
                H[:,k] += self.source.H(np.array([x,y,z]), self.material_data['k'][k]) \
                            * (self.material_data['eps_b'][k]/self.material_data['mu_b'][k])**0.5
        return H

    def flux_from_particle(self, i, inc=False):
        """Determine the scattered flux from a single particle

            Arguments:
                i         Particle index
                inc       Include the incident field (bool, default=False)
            
            Returns: flux[M], M = number of wavelengths
        """
        r = self.spheres.radius[i]
        X,Y,Z,THETA,PHI,tau,phi = discrete_sphere(r, self.Ntheta, self.Nphi, self.spheres.position[i])
        rhat,*_ = sph_unit_vectors(THETA, PHI)

        E_all = self.E_field(X, Y, Z, inc)
        H_all = self.H_field(X, Y, Z, inc)

        flux = np.zeros(self.Nfreq, dtype=float)

        for k in range(self.Nfreq):
            E = E_all[:,k]
            H = H_all[:,k]

            S = np.real(np.einsum('ijk,jxy,kxy->ixy', levi, E, np.conj(H)))
            dA = r**2

            integrand = np.einsum('ixy,ixy->xy', S, rhat)*dA
            flux[k] = simps_2d(tau, phi, integrand)

        return flux

    def force_on_particle(self, i, inc=True):
        """Determine the force on a single particle

            Arguments:
                i         Particle index
                inc       Include the incident field (bool, default=False)
            
            Returns: (F[3,M],T[3,M]), M = number of wavelengths
        """
        r = self.spheres.radius[i]
        X,Y,Z,THETA,PHI,tau,phi = discrete_sphere(r, self.Ntheta, self.Nphi, self.spheres.position[i])
        rhat,*_ = sph_unit_vectors(THETA, PHI)

        E_all = self.E_field(X, Y, Z, inc)
        H_all = self.H_field(X, Y, Z, inc)

        F = np.zeros([3, self.Nfreq], dtype=float)
        T = np.zeros([3, self.Nfreq], dtype=float)

        for k in range(self.Nfreq):
            E = E_all[:,k]
            H = H_all[:,k]

            eps_b = self.material_data['eps_b'][k]
            mu_b = self.material_data['mu_b'][k]
            sigma = eps_b*np.einsum('ixy,jxy->ijxy', E, np.conj(E)) \
                    + mu_b*np.einsum('ixy,jxy->ijxy', H, np.conj(H)) \
                    - 0.5*np.einsum('ij,xy->ijxy', np.identity(3), eps_b*np.sum(np.abs(E)**2, axis=0)) \
                    - 0.5*np.einsum('ij,xy->ijxy', np.identity(3), mu_b*np.sum(np.abs(H)**2, axis=0))

            # compute F
            dA = r**2
            integrand = np.einsum('ijxy,jxy->ixy', sigma, rhat)*dA
            F[:,k] = np.array([simps_2d(tau, phi, integrand[x].real) for x in range(3)])

            # compute T
            integrand = np.einsum('imn,mxy,njxy,jxy->ixy', levi, r*rhat, sigma, rhat)*dA
            T[:,k] = np.array([simps_2d(tau, phi, integrand[x].real) for x in range(3)])

        return F,T

    def flux(self, inc=False):
        """Determine the scattered flux from every particle

            Arguments:
                inc       Include the incident field (bool, default=False)
            
            Returns: flux[M,N], N = number of particle, M = number of wavelengths
        """
        flux_data = np.zeros([self.Nfreq, self.Nparticles], dtype=float)
        for i in range(self.Nparticles):
            flux_data[:,i] = self.flux_from_particle(i, inc=inc)

        return flux_data

    def force(self, inc=True):
        """Determine the force on every particle

            Arguments:
                inc       Include the incident field (bool, default=False)
            
            Returns: (F[3,M,N],T[3,M,N]), 
                     N = number of particles, M = number of wavelengths
        """
        force_data = np.zeros([3, self.Nfreq, self.Nparticles], dtype=float)
        torque_data = np.zeros([3, self.Nfreq, self.Nparticles], dtype=float)
        for i in range(self.Nparticles):
            force_data[...,i], torque_data[...,i] = self.force_on_particle(i, inc=inc)

        return force_data, torque_data

    def update_position(self, position):
        """Update the positions of the spheres

            Arguments
                position[N,3]       new particle positions
        """
        self.spheres.position = position

        if (self.interactions):
            self._solve_interactions()
        else:
            self._set_without_interactions()

    def _set_without_interactions(self):
        pos = self.spheres.position.T
        for k in range(self.Nfreq):
            Einc = self.source.E(pos,self.material_data['k'][k])
            self.p[...,k] = Einc[:2,:]

    #TODO vectorize for loops. Avoid transpose of position->pass x,y,z to source instead...?
    def _solve_interactions(self):
        pos = self.spheres.position.T
        
        identity = np.zeros(shape = (2, self.Nparticles, 2, self.Nparticles), dtype=np.complex)
        np.einsum('xixi->xi', identity)[...] = 1
        
        for k in range(self.Nfreq):
            MieMatrix = np.zeros(shape = (2, self.Nparticles, 2, self.Nparticles), dtype=np.complex)
            Einc = self.source.E(pos,self.material_data['k'][k])
            Einc = Einc[:2,:]

            for i in range(self.Nparticles):
                for j in range(self.Nparticles):
                    if i == j: continue
                    pi = self.spheres.position[i]
                    pj = self.spheres.position[j]
                    dji = pi -  pj
                    r_ji = np.linalg.norm(dji)
                    theta_ji = np.arccos(dji[2]/r_ji)
                    phi_ji = np.arctan2(dji[1], dji[0])
                    
                    rhat = np.array([np.sin(theta_ji)*np.cos(phi_ji), np.sin(theta_ji)*np.sin(phi_ji), np.cos(theta_ji)])
                    that = np.array([np.cos(theta_ji)*np.cos(phi_ji), np.cos(theta_ji)*np.sin(phi_ji), -np.sin(theta_ji)])
                    phat = np.array([-np.sin(phi_ji), np.cos(phi_ji), np.zeros_like(theta_ji)])
                    
                    E_func = miepy.scattering.scattered_E(self.a[j,k], self.b[j,k], self.material_data['k'][k])
                    xsol = E_func(r_ji, theta_ji, phi_ji)
                    ysol = E_func(r_ji, theta_ji, phi_ji - np.pi/2)
                    xsol = xsol[0]*rhat + xsol[1]*that + xsol[2]*phat
                    ysol = ysol[0]*rhat + ysol[1]*that + ysol[2]*phat
                    
                    MieMatrix[:,i,0,j] = xsol[:2]
                    MieMatrix[:,i,1,j] = ysol[:2]

            A = identity - MieMatrix
            sol = np.linalg.tensorsolve(A, Einc)
            self.p[...,k] = sol
        
if __name__ == "__main__":
    system = gmt(spheres([[0,-50e-9,0],[0,50e-9,0]], 20e-9, miepy.constant_material(2)), 
                  miepy.sources.x_polarized_plane_wave(),
                  600e-9, 2)
    from IPython import embed
    embed()
