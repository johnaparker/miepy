"""
The Generalized Mie Theory (GMT) for a collection of spheres.
"""
import numpy as np
import miepy
from miepy.special_functions import riccati_1,riccati_2,vector_spherical_harmonics
from miepy.utils import atleast
from functools import partial

#TODO: make several properties... such as wavelength, source, position, etc.
#TODO: prefer manual solve calls ratehr than auto solve calls (or have an auto_solve option)
#TODO: swap position indices, so that [N,3] => [3,N]
class sphere_cluster:
    """Solve Generalized Mie Theory for an N particle sphere cluster in an arbitray source profile"""
    def __init__(self, *, position, radius, material, source, wavelength,
                 lmax, medium=None, origin=None, symmetry=None, interface=None,
                 interactions=True):
        """Arguments:
               position[N,3] or [3]    sphere positions
               radius[N] or scalar     sphere radii
               material[N] or scalar   sphere materials
               source        source object specifying the incident E and H functions
               wavelength    wavelength to solve the system at
               lmax          maximum number of orders to use in angular momentum expansion (int)
               medium        (optional) material medium (must be non-absorbing; default=vacuum)
               origin        (optional) system origin around which to compute cluster quantities (default = [0,0,0]). Choose 'auto' to automatically choose origin as center of geometry.
               symmetry      (optional) specify system symmetries (default: no symmetries)
               interface     (optional) include an infinite interface (default: no interface)
               interactions  (optional) If True, include particle interactions (bool, default=True) 
        """
        ### sphere properties
        self.position = np.asarray(np.atleast_2d(position), dtype=float)
        self.radius = atleast(radius, dim=1, length=self.position.shape[0], dtype=float)
        self.material = atleast(material, dim=1, length=self.position.shape[0], dtype=np.object)
        if (self.position.shape[0] != self.radius.shape[0] != self.material.shape[0]):
            raise ValueError("The shapes of position, radius, and material do not match")
        self.Nparticles = self.radius.shape[0]
        self.symmetry = symmetry
        self.interface = interface

        ### system properties
        self.source = source
        self.wavelength = wavelength
        self.lmax = lmax
        self.rmax = lmax*(lmax + 2)
        self.interactions = interactions

        ### set the origin
        self.auto_origin = False    
        if origin is None:
            self.origin = np.zeros(3)
        elif origin is 'auto':
            self.auto_origin = True
            self.origin = np.average(self.position, axis=0)
        else:
            self.origin = np.asarray(origin)

        ### set the medium
        if medium is None:
            self.medium = miepy.constant_material(eps=1.0, mu=1.0)
        else:
            self.medium = medium
            if (self.medium.eps(self.wavelength).imag != 0)  \
                    or (self.medium.mu(self.wavelength).imag != 0):
                raise ValueError('medium must be non-absorbing')

        ### build material data of particles
        self.material_data = miepy.material_functions.material_struct(self.material, self.medium, wavelength=self.wavelength)

        ### mie coefficients
        self.mie_scat = np.zeros([self.Nparticles, 2, self.lmax], dtype=complex)
        self.mie_int = np.zeros([self.Nparticles, 2, self.lmax], dtype=complex)

        for i in range(self.Nparticles):
            conducting = (self.material[i].name == 'metal')
            for n in range(1, self.lmax+1):
                self.mie_scat[i,:,n-1] = \
                    miepy.mie_single.mie_sphere_scattering_coefficients(self.radius[i],
                    n, self.material_data.eps[i], self.material_data.mu[i],
                    self.material_data.eps_b, self.material_data.mu_b, self.material_data.k_b,
                    conducting=conducting)

                self.mie_int[i,:,n-1] = \
                    miepy.mie_single.mie_sphere_interior_coefficients(self.radius[i],
                    n, self.material_data.eps[i], self.material_data.mu[i],
                    self.material_data.eps_b, self.material_data.mu_b, self.material_data.k_b,
                    conducting=conducting)

        ### modified coefficients
        self.p_inc  = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)
        self.p_scat = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)
        self.p_int  = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)
        self.p_src  = np.zeros([self.Nparticles, 2, self.rmax], dtype=complex)

        ### cluster coefficients
        self.p_cluster = None

        ### solve the interactions
        self.solve()

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    Nparticles = {self.Nparticles}
    source = {self.source}
    wavelength = {self.wavelength:.2e}
    medium = {self.medium}
    lmax = {self.lmax}
    origin = {self.origin}
    '''

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
    
    def E_source(self, x1, x2, x3, far=False, spherical=False):
        """Compute the electric field from the source

        Arguments:
            x1        x/r position (array-like) 
            x2        y/theta position (array-like) 
            x3        z/phi position (array-like) 
            far       (optional) use expressions valid only for far-field (bool, default=False)
            spherical (optional) input/output in spherical coordinates (bool, default=False)

        Returns: E[3,...]
        """
        E = self.source.E_field(x1, x2, x3, self.material_data.k_b, far=far, spherical=spherical)
        return E

    def H_source(self, x1, x2, x3, far=False, spherical=False):
        """Compute the magnetic field from the source

        Arguments:
            x1        x/r position (array-like) 
            x2        y/theta position (array-like) 
            x3        z/phi position (array-like) 
            far       (optional) use expressions valid only for far-field (bool, default=False)
            spherical (optional) input/output in spherical coordinates (bool, default=False)

        Returns: H[3,...]
        """
        factor = (self.material_data.eps_b/self.material_data.mu_b)**0.5
        H = self.source.H_field(x1, x2, x3, self.material_data.k_b, far=far, spherical=spherical)
        return factor*H

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
            expand = partial(miepy.expand_E, mode=miepy.vsh_mode.outgoing)

        for i in range(self.Nparticles):
            rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
            E_sph = expand(self.p_scat[i], self.material_data.k_b)(rad,theta,phi)
            E += miepy.coordinates.vec_sph_to_cart(E_sph, theta, phi)

        if source:
            E += self.E_source(x, y, z, far=far, spherical=False)

            if self.interface is not None:
                idx = z <= self.interface.z
                reflected = self.source.reflect(self.interface, self.medium, self.wavelength)
                E += reflected.E_field(x[idx], y[idx], z[idx], self.material_data.k_b, far=far, spherical=False)

                transmitted = self.source.transmit(self.interface, self.medium, self.wavelength)
                E += transmitted.E_field(x[~idx], y[~idx], z[~idx], self.material_data.k_b, far=far, spherical=False)

        #TODO: what if x is scalar...
        if interior and not mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.radius[i]**2)
                k_int = 2*np.pi*self.material_data.n[i]/self.wavelength

                rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
                E_sph = miepy.expand_E(self.p_int[i], k_int, 
                                mode=miepy.vsh_mode.interior)(rad[idx], theta[idx], phi[idx])
                E[:,idx] = miepy.coordinates.vec_sph_to_cart(E_sph, theta[idx], phi[idx])

        if mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.radius[i]**2)
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
            expand = partial(miepy.expand_H, mode=miepy.vsh_mode.outgoing)

        for i in range(self.Nparticles):
            rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
            H_sph = expand(self.p_scat[i], self.material_data.k_b, eps=self.material_data.eps_b,
                               mu=self.material_data.mu_b)(rad,theta,phi)
            H += miepy.coordinates.vec_sph_to_cart(H_sph, theta, phi)

        if source:
            H += self.H_source(x, y, z, far=far, spherical=False)

            if self.interface is not None:
                idx = z <= self.interface.z
                reflected = self.source.reflect(self.interface, self.medium, self.wavelength)
                H += reflected.H_field(x[idx], y[idx], z[idx], self.material_data.k_b, far=far, spherical=False)

                transmitted = self.source.transmit(self.interface, self.medium, self.wavelength)
                H += transmitted.H_field(x[~idx], y[~idx], z[~idx], self.material_data.k_b, far=far, spherical=False)

        #TODO: what if x is scalar...
        if interior and not mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.radius[i]**2)
                k_int = 2*np.pi*self.material_data.n[i]/self.wavelength

                rad, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position[i])
                H_sph = miepy.expand_H(self.p_int[i], k_int, 
                            eps=self.material_data.eps[i], mu=self.material_data.mu[i],
                            mode=miepy.vsh_mode.interior)(rad[idx], theta[idx], phi[idx])
                H[:,idx] = miepy.coordinates.vec_sph_to_cart(H_sph, theta[idx], phi[idx])

        if mask and not far:
            for i in range(self.Nparticles):
                x0, y0, z0 = self.position[i]
                idx = ((x - x0)**2 + (y - y0)**2 + (z - z0)**2 < self.radius[i]**2)
                H[:,idx] = 0

        #TODO: does this depend on the origin?
        if spherical:
            H = miepy.coordinates.vec_cart_to_sph(H, theta=x2, phi=x3)

        return H

    def E_angular(self, theta, phi, radius=None, source=False):
        """Compute the electric field due to all particles in the far-field in spherical coordinates
             
            Arguments:
                theta    theta position (array-like) 
                phi      phi position (array-like) 
                radius   r position (default: large value)
                source   (bool) include the angular source fields (default: False)
        """
        #TODO better expression far default far-radius
        if radius is None:
            radius = 1e6*2*np.pi/self.material_data.k_b

        return self.E_field(radius, theta, phi, interior=False, source=source, far=True, spherical=True)[1:]

    def H_angular(self, theta, phi, radius=None, source=False):
        """Compute the magnetic field due to all particles in the far-field in spherical coordinates
             
            Arguments:
                theta    theta position (array-like) 
                phi      phi position (array-like) 
                radius   r position (default: large value)
                source   (bool) include the angular source fields (default: False)
        """
        #TODO better expression far default far-radius
        if radius is None:
            radius = 1e6*2*np.pi/self.material_data.k_b

        return self.H_field(radius, theta, phi, interior=False, source=source, far=True, spherical=True)[1:]

    def cross_sections_per_multipole(self, lmax=None):
        """Compute the scattering, absorption, and extinction cross-section of the cluster per multipole

        Arguments:
            lmax    (optional) compute scattering for up to lmax terms (defult: self.lmax)
        """

        if lmax is None:
            lmax = self.lmax

        self.solve_cluster_coefficients(lmax)
        p0 = self.source.structure(self.origin, self.material_data.k_b, lmax)

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

    def local_density_of_states(self, enhancement=True):
        """Compute the local density of states (LDOS)

        Arguments:
            enhancement      (bool) if True, return the relative enhancement LDOS (default: True)
        """
        if type(self.source) != miepy.sources.point_dipole:
            raise ValueError("The source must be a single point dipole to compute the local density of states, not of type '{}'".format(
                            type(self.source)))

        factor = -2/np.pi*self.material_data.eps_b
        pos = self.source.position + np.array([1e-12, 0, 0])
        E = self.E_field(*pos)
        p = 1j*self.source.E_field(*pos, self.material_data.k_b)

        projection = self.source.direction
        E_comp = np.dot(E, projection)
        p_comp = np.dot(p, projection)

        if enhancement:
            return np.real(1j*E_comp)/np.real(p_comp)
        else:
            return factor*np.real(E_comp*np.conj(p_comp))/np.abs(p_comp)**2

    def update_position(self, position):
        """Update the positions of the spheres

            Arguments
                position[N,3]       new particle positions
        """
        self.position = np.asarray(np.atleast_2d(position), dtype=float)
        self._reset_cluster_coefficients()

        if self.auto_origin:
            self.origin = np.average(self.position, axis=0)

        self.solve()

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
        for i in range(self.Nparticles):
            self.p_src[i] = self.source.structure(self.position[i], self.material_data.k_b, self.lmax)

        if self.interface is not None:
            reflected = self.source.reflect(self.interface, self.medium, self.wavelength)

            for i in range(self.Nparticles):
                self.p_src[i] += reflected.structure(self.position[i], self.material_data.k_b, self.lmax)

    def _solve_without_interactions(self):
        self.p_inc[...] = self.p_src

        for r,n,m in miepy.mode_indices(self.lmax):
            self.p_scat[...,r] = self.p_inc[...,r]*self.mie_scat[...,n-1]
            self.p_int[...,r] = self.p_inc[...,r]*self.mie_int[:,::-1,n-1]

    def _solve_interactions(self):
        if self.symmetry is None:
            agg_tmatrix = miepy.interactions.sphere_aggregate_tmatrix(self.position, self.mie_scat,
                                      self.material_data.k_b)
        else:
            agg_tmatrix = miepy.interactions.sphere_aggregate_tmatrix_periodic(self.position, self.mie_scat,
                                      self.material_data.k_b, self.symmetry, self.source.k_hat)

        if self.interface is not None:
            r0 = self.interface.reflection_coefficients(theta=0, wavelength=self.wavelength, medium=self.medium)[0]
            z = self.interface.z
            R_matrix = miepy.interactions.reflection_matrix_nia(self.position, self.mie_scat, self.material_data.k_b, r0, z)
            agg_tmatrix -= R_matrix


        self.p_inc[...] = miepy.interactions.solve_linear_system(agg_tmatrix, self.p_src, method=miepy.solver.bicgstab)

        for r,n,m in miepy.mode_indices(self.lmax):
            self.p_scat[...,r] = self.p_inc[...,r]*self.mie_scat[...,n-1]
            self.p_int[...,r] = self.p_inc[...,r]*self.mie_int[:,::-1,n-1]
