"""
The Generalized Mie Theory (GMT) for a collection of particles
"""
import numpy as np
import miepy
from miepy.special_functions import riccati_1,riccati_2,vector_spherical_harmonics
from miepy.interface import atleast

#TODO: make several properties... such as wavelength, source, position, etc.
#TODO: prefer manual solve calls ratehr than auto solve calls (or have an auto_solve option)
#TODO: swap position indices, so that [N,3] => [3,N]
class cluster:
    """Solve Generalized Mie Theory: N particle cluster in an arbitray source profile"""
    def __init__(self, *, particles, lmax,
                 source=None, wavelength=None, medium=None, origin=None,
                 interactions=True):
        """Arguments:
               particles[N]  list of particle objects
               lmax          maximum number of orders to use in angular momentum expansion (int)
               source        (optional) source object specifying the incident E and H functions (deault: defer source until later)
               wavelength    (optional) wavelength to solve the system at (default: defer wavelength until later)
               medium        (optional) material medium (must be non-absorbing; default=vacuum)
               origin        (optional) system origin around which to compute cluster quantities (default = [0,0,0]). Choose 'auto' to automatically choose origin as center of geometry.
               interactions  (optional) If True, include particle interactions (bool, default=True) 
        """
        ### system properties
        self.source = source
        self.wavelength = wavelength
        self.lmax = lmax
        self.rmax = miepy.vsh.lmax_to_rmax(lmax)
        self.interactions = interactions

        ### set the medium
        if medium is None:
            self.medium = miepy.materials.vacuum()
        else:
            self.medium = medium
            if (self.medium.eps(self.wavelength).imag != 0)  \
                    or (self.medium.mu(self.wavelength).imag != 0):
                raise ValueError('medium must be non-absorbing')

        ### particle properties
        self.particles  = particles if isinstance(particles, list) else [particles]
        self.Nparticles = len(self.particles)
        self.position = np.empty([self.Nparticles,3], dtype=float)
        self.material = np.empty([self.Nparticles], dtype=object)
        self.tmatrix  = np.empty([self.Nparticles, 2, self.rmax, 2, self.rmax], dtype=complex)

        ### calculate T-matrices 
        tmatrices = {}
        for i in range(self.Nparticles):
            self.position[i] = self.particles[i].position
            self.material[i] = self.particles[i].material
            
            key = self.particles[i]._dict_key(wavelength)
            if key in tmatrices:
                self.particles[i].tmatrix_fixed = tmatrices[key]
                self.particles[i]._rotate_fixed_tmatrix()
                self.tmatrix[i] = self.particles[i].tmatrix
            else:
                self.tmatrix[i]  = self.particles[i].compute_tmatrix(lmax, wavelength, self.medium.eps(wavelength))
                tmatrices[key] = self.particles[i].tmatrix_fixed

        ### set the origin
        self.auto_origin = False    
        if origin is None:
            self.origin = np.zeros(3)
        elif origin is 'auto':
            self.auto_origin = True
            self.origin = np.average(self.position, axis=0)
        else:
            self.origin = np.asarray(origin)

        ### build material data of particles
        self.material_data = miepy.materials.material_struct(self.material, self.medium, wavelength=self.wavelength)

        ### expansion coefficients
        self.p_inc  = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)
        self.p_scat = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)
        self.p_int  = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)
        self.p_src  = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)

        ### cluster expansion coefficients
        self.p_cluster = None

        ### solve the interactions
        self.solve()

    #TODO: interface more like E_field
    def E_field_from_particle(self, i, x, y, z, source=True):
        """Compute the electric field around particle i
             
            Arguments:
                i        particle number
                x        x position (array-like) 
                y        y position (array-like) 
                z        z position (array-like) 
                source   Include the source field (bool, default=True)

            Returns: E[3,...]
        """
        rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])

        E_sph = miepy.expand_E(self.p_scat[i], self.material_data.k_b,
                     mode=miepy.vsh_mode.outgoing)(rad,theta,phi)
        Escat = miepy.coordinates.vec_sph_to_cart(E_sph, theta, phi)

        p = self.p_inc[i]
        if not source:
            p -= self.p_src[i]
        E_sph = miepy.expand_E(p, self.material_data.k_b,
                     mode=miepy.vsh_mode.incident)(rad,theta,phi)
        Einc = miepy.coordinates.vec_sph_to_cart(E_sph, theta, phi)

        return Escat + Einc

    def H_field_from_particle(self, i, x, y, z, source=True):
        """Compute the magnetic field around particle i
             
            Arguments:
                i        particle number
                x        x position (array-like) 
                y        y position (array-like) 
                z        z position (array-like) 
                source   Include the source field (bool, default=True)

            Returns: H[3,...]
        """
        rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])

        H_sph = miepy.expand_H(self.p_scat[i], self.material_data.k_b,
                  mode=miepy.vsh_mode.outgoing, eps=self.material_data.eps_b,
                  mu=self.material_data.mu_b)(rad,theta,phi)
        Hscat = miepy.coordinates.vec_sph_to_cart(H_sph, theta, phi)

        p = self.p_inc[i]
        if not source:
            p -= self.p_src[i]

        H_sph = miepy.expand_H(p, self.material_data.k_b,
                  mode=miepy.vsh_mode.incident, eps=self.material_data.eps_b,
                  mu=self.material_data.mu_b)(rad,theta,phi)
        Hinc = miepy.coordinates.vec_sph_to_cart(H_sph, theta, phi)

        return Hscat + Hinc
    
    def E_field(self, x1, x2, x3, interior=True, source=True, mask=False, far=False, spherical=False):
        """Compute the electric field due to all particles
             
            Arguments:
                x1        x/r position (array-like) 
                x2        y/theta position (array-like) 
                x3        z/phi position (array-like) 
                interior  (optional) compute interior fields (bool, default=True)
                source    (optional) include the source field (bool, default=True)
                mask      (optional) set interior fields to 0 (bool, default=False)
                far       (optional) use expressions valid only for far-field (bool, default=False)
                spherical (optional) input/output in spherical coordinates (bool, default=False)

            Returns: E[3,...]
        """
        x1, x2, x3 = (np.asarray(x) for x in (x1, x2, x3))
        shape = max(*[x.shape for x in (x1, x2, x3)], key=len)
        E = np.zeros((3,) + shape, dtype=complex)

        if spherical:
            (x, y, z) = miepy.coordinates.sph_to_cart(x1, x2, x3, origin=self.origin)
        else:
            (x, y, z) = (x1, x2, x3)

        if far:
            expand = miepy.expand_E_far
        else:
            expand = lambda p,k: miepy.expand_E(p, k, mode=miepy.vsh_mode.outgoing)

        for i in range(self.Nparticles):
            rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
            E_sph = expand(self.p_scat[i], self.material_data.k_b)(rad,theta,phi)
            E += miepy.coordinates.vec_sph_to_cart(E_sph, theta, phi)

        #TODO: [x,y,z] to x,y,z
        if source:
            E += self.source.E_field(x, y, z, self.material_data.k_b)

        #TODO: what if x is scalar...
        if interior and not mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.particles[i].enclosed_radius()**2)
                k_int = 2*np.pi*self.material_data.n[i]/self.wavelength

                rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
                E_sph = miepy.expand_E(self.p_int[i], k_int, 
                                mode=miepy.vsh_mode.interior)(rad[idx], theta[idx], phi[idx])
                E[:,idx] = miepy.coordinates.vec_sph_to_cart(E_sph, theta[idx], phi[idx])

        if mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.particles[i].enclosed_radius()**2)
                E[:,idx] = 0

        #TODO: does this depend on the origin?
        if spherical:
            E = miepy.coordinates.vec_cart_to_sph(E, theta=x2, phi=x3)
        
        return E

    def H_field(self, x1, x2, x3, interior=True, source=True, mask=False, far=False, spherical=False):
        """Compute the magnetic field due to all particles
             
            Arguments:
                x1        x/r position (array-like) 
                x2        y/theta position (array-like) 
                x3        z/phi position (array-like) 
                interior  (optional) compute interior fields (bool, default=True)
                source    (optional) include the source field (bool, default=True)
                mask      (optional) set interior fields to 0 (bool, default=False)
                far       (optional) use expressions valid only for far-field (bool, default=False)
                spherical (optional) input/output in spherical coordinates (bool, default=False)

            Returns: H[3,...]
        """
        x1, x2, x3 = (np.asarray(x) for x in (x1, x2, x3))
        shape = max(*[x.shape for x in (x1, x2, x3)], key=len)
        H = np.zeros((3,) + shape, dtype=complex)

        if spherical:
            (x, y, z) = miepy.coordinates.sph_to_cart(x1, x2, x3, origin=self.origin)
        else:
            (x, y, z) = (x1, x2, x3)

        if far:
            expand = miepy.expand_H_far
        else:
            expand = lambda p,k,eps,mu: miepy.expand_H(p, k, mode=miepy.vsh_mode.outgoing, eps=eps, mu=mu)

        for i in range(self.Nparticles):
            rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
            H_sph = expand(self.p_scat[i], self.material_data.k_b, eps=self.material_data.eps_b,
                               mu=self.material_data.mu_b)(rad,theta,phi)
            H += miepy.coordinates.vec_sph_to_cart(H_sph, theta, phi)

        #TODO: [x,y,z] to x,y,z
        if source:
            factor = (self.material_data.eps_b/self.material_data.mu_b)**0.5
            H += factor*self.source.H_field(x, y, z, self.material_data.k_b)

        #TODO: what if x is scalar...
        if interior and not mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.particles[i].enclosed_radius()**2)
                k_int = 2*np.pi*self.material_data.n[i]/self.wavelength

                rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
                H_sph = miepy.expand_H(self.p_int[i], k_int, 
                            eps=self.material_data.eps[i], mu=self.material_data.mu[i],
                            mode=miepy.vsh_mode.interior)(rad[idx], theta[idx], phi[idx])
                H[:,idx] = miepy.coordinates.vec_sph_to_cart(H_sph, theta[idx], phi[idx])

        if mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.particles[i].enclosed_radius()**2)
                H[:,idx] = 0

        #TODO: does this depend on the origin?
        if spherical:
            H = miepy.coordinates.vec_cart_to_sph(H, theta=x2, phi=x3)
        
        return H

    def E_angular(self, theta, phi, radius=None):
        """Compute the electric field due to all particles in the far-field in spherical coordinates
             
            Arguments:
                theta    theta position (array-like) 
                phi      phi position (array-like) 
                radius   r position (default: large value)
        """
        #TODO better expression far default far-radius
        if radius is None:
            radius = 1e10*self.wavelength

        return self.E_field(radius, theta, phi, interior=False, source=False, far=True, spherical=True)

    def H_angular(self, theta, phi, radius=None):
        """Compute the magnetic field due to all particles in the far-field in spherical coordinates
             
            Arguments:
                theta    theta position (array-like) 
                phi      phi position (array-like) 
                radius   r position (default: large value)
        """
        #TODO better expression far default far-radius
        if radius is None:
            radius = 1e10*self.wavelength

        return self.H_field(radius, theta, phi, interior=False, source=False, far=True, spherical=True)

    def cross_sections_per_multipole(self, lmax=None):
        """Compute the scattering, absorption, and extinction cross-section of the cluster per multipole

        Arguments:
            lmax    (optional) compute scattering for up to lmax terms (defult: self.lmax)
        """

        if lmax is None:
            lmax = self.lmax

        self.solve_cluster_coefficients(lmax)
        p0 = self.source.structure(self.origin, self.material_data.k_b, lmax, radius=lmax/self.material_data.k_b)

        return miepy.flux.cluster_cross_sections(self.p_cluster, p0, self.material_data.k_b)

    def cross_sections(self):
        """Compute the scattering, absorption, and extinction cross-section of the cluster"""

        Cscat = 0
        Cabs = 0
        Cext = 0

        for i in range(self.Nparticles):
            C = self.cross_sections_of_particle(i)
            Cscat += C.scattering
            Cabs += C.absorption
            Cext += C.extinction

        return miepy.flux.cross_sections(Cscat, Cabs, Cext)

    def cross_sections_per_multipole_of_particle(self, i):
        """Compute the scattering, absorption, and extinction cross-section per multipole of a single particle

        Arguments:
            i    particle index
        """
        return miepy.flux.particle_cross_sections(self.p_scat[i], self.p_inc[i], 
                    self.p_src[i], self.material_data.k_b)

    def cross_sections_of_particle(self, i):
        """Compute the scattering, absorption, and extinction cross-section of a single particle

        Arguments:
            i    particle index
        """

        C = self.cross_sections_per_multipole_of_particle(i)
        return miepy.flux.cross_sections(*[np.sum(A) for A in C])

    def force_on_particle(self, i, source=True):
        """Determine the force on a single particle

            Arguments:
                i         Particle index
                source    Include the source field (bool, default=True)
            
            Returns: F[3]
        """

        if source:
            p_inc = self.p_inc
        else:
            p_inc = self.p_inc - self.p_src

        F = miepy.forces.force(self.p_scat[i].reshape(-1), p_inc[i].reshape(-1),
                self.material_data.k_b, self.material_data.eps_b, self.material_data.mu_b)

        return F

    def torque_on_particle(self, i, source=True):
        """Determine the torque on a single particle

            Arguments:
                i         Particle index
                source    Include the source field (bool, default=True)
            
            Returns: T[3]
        """

        if source:
            p_inc = self.p_inc
        else:
            p_inc = self.p_inc - self.p_src

        T = miepy.forces.torque(self.p_scat[i].reshape(-1), p_inc[i].reshape(-1),
                self.material_data.k_b, self.material_data.eps_b, self.material_data.mu_b)

        return T

    def force(self, source=True):
        """Determine the force on every particle

            Arguments:
                source    Include the source field (bool, default=False)
            
            Returns: F[3,Nparticles]
        """
        F = np.zeros([3, self.Nparticles], dtype=float)
        for i in range(self.Nparticles):
            F[:,i] = self.force_on_particle(i, source=source)

        return F

    def torque(self, source=True):
        """Determine the torque on every particle

            Arguments:
                source    Include the source field (bool, default=False)
            
            Returns: T[3,Nparticles]
        """
        T = np.zeros([3, self.Nparticles], dtype=float)
        for i in range(self.Nparticles):
            T[:,i] = self.torque_on_particle(i, source=source)

        return T

    def update(self, position=None, orientation=None):
        """Update properties of the particles

            Arguments
                position[N,3]     new particle positions
                orientation[N]    new particle orientations (array of quaternions)
        """

        if position is not None:
            self.position = np.asarray(np.atleast_2d(position), dtype=float)

            if self.auto_origin:
                self.origin = np.average(self.position, axis=0)

        if orientation is not None:
            for i in raneg(self.Nparticles):
                self.particles[i].orientation = orientation[i]
                self.tmatrix[i] = self.particles[i].tmatrix

        self._reset_cluster_coefficients()

        if position is not None or orientation is None:
            self.solve()
        elif position is None and orientation is not None:
            if self.interactions:
                self._solve_interactions()
            else:
                self._solve_without_interactions()

    def solve_cluster_coefficients(self, lmax=None):
        """Solve for the p,q coefficients of the entire cluster around the origin

        Arguments:
            lmax    (optional) compute scattering for up to lmax terms (default: self.lmax)
        """

        if lmax is None:
            lmax = self.lmax

        self.p_cluster = miepy.cluster_coefficients(self.position, 
                self.p_scat, self.material_data.k_b, origin=self.origin, lmax=lmax)

    def solve(self, wavelength=None, source=None):
        """solve for the p,q incident and scattering coefficients

           Arguments:
               wavelength   wavelength to solve at (default: current wavelength)
               source       source to use (default: current source). If current source is also None, solve the particle's T-matrix instead
        """
        self._solve_source_decomposition()
        if self.interactions:
            self._solve_interactions()
        else:
            self._solve_without_interactions()

    def _reset_cluster_coefficients(self):
        self.p_cluster = None

    def _solve_source_decomposition(self):
        pos = self.position.T
        for i in range(self.Nparticles):
            pos = self.position[i]
            self.p_src[i] = self.source.structure(pos, self.material_data.k_b, self.lmax, self.particles[i].enclosed_radius())

    def _solve_without_interactions(self):
        self.p_inc[...] = self.p_src
        self.p_scat[...] = np.einsum('naibj,nbj->nai', self.tmatrix, self.p_inc)

    def _solve_interactions(self):
        agg_tmatrix = miepy.interactions.particle_aggregate_tmatrix(self.position, self.tmatrix,
                                  self.material_data.k_b)
        self.p_inc[...] = miepy.interactions.solve_linear_system(agg_tmatrix, self.p_src, method=miepy.solver.bicgstab)

        self.p_scat[...] = np.einsum('naibj,nbj->nai', self.tmatrix, self.p_inc)
