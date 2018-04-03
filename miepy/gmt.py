"""
The Generalized Mie Theory (GMT) for a collection of spheres.
"""
import numpy as np
import miepy
from my_pytools.my_numpy.array import atleast
from miepy.special_functions import riccati_1,riccati_2,vector_spherical_harmonics

#TODO: make several properties... such as wavelength, source, position, etc.
#TODO: prefer manual solve calls ratehr than auto solve calls (or have an auto_solve option)
#TODO: swap position indices, so that [N,3] => [3,N]
class cluster:
    """Solve Generalized Mie Theory: N particle cluster in an arbitray source profile"""
    def __init__(self, position, radius, material, Lmax,
                 source=None, wavelength=None, medium=None, origin=None,
                 interactions=True):
        """Arguments:
               position[N,3] or [3]    sphere positions
               radius[N] or scalar     sphere radii
               material[N] or scalar   sphere materials
               Lmax          maximum number of orders to use in angular momentum expansion (int)
               source        (optional) source object specifying the incident E and H functions (deault: defer source until later)
               wavelength    (optional) wavelength to solve the system at (default: defer wavelength until later)
               medium        (optional) material medium (must be non-absorbing; default=vacuum)
               origin        (optional) system origin around which to compute cluster quantities (default = [0,0,0]). Choose 'auto' to automatically choose origin as center of geometry.
               interactions  (optional) If True, include particle interactions (bool, default=True) 
        """
        ### sphere properties
        self.position = np.asarray(np.atleast_2d(position), dtype=float)
        self.radius = atleast(radius, dim=1, length=self.position.shape[0], dtype=float)
        self.material = atleast(material, dim=1, length=self.position.shape[0], dtype=np.object)
        if (self.position.shape[0] != self.radius.shape[0] != self.material.shape[0]):
            raise ValueError("The shapes of position, radius, and material do not match")
        self.Nparticles = self.radius.shape[0]

        ### system properties
        self.source = source
        self.wavelength = wavelength
        self.Lmax = Lmax
        self.rmax = Lmax*(Lmax + 2)
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
        self.material_data = miepy.materials.material_struct(self.material, self.medium, wavelength=self.wavelength)

        ### mie coefficients
        self.a = np.zeros([self.Nparticles, self.Lmax], dtype=complex)
        self.b = np.zeros([self.Nparticles, self.Lmax], dtype=complex)

        #TODO: a,b function calls instead of class
        for i in range(self.Nparticles):
            sphere = miepy.single_mie_sphere(self.radius[i], self.material[i],
                        self.wavelength, self.Lmax, self.medium)
            self.a[i], self.b[i] = sphere.solve_exterior()

        ### modified coefficients
        self.p_scat = np.zeros([self.Nparticles, self.rmax], dtype=complex)
        self.q_scat = np.zeros([self.Nparticles, self.rmax], dtype=complex)
        self.p_inc  = np.zeros([self.Nparticles, self.rmax], dtype=complex)
        self.q_inc  = np.zeros([self.Nparticles, self.rmax], dtype=complex)
        self.p_src  = np.zeros([self.Nparticles, self.rmax], dtype=complex)
        self.q_src  = np.zeros([self.Nparticles, self.rmax], dtype=complex)

        ### cluster coefficients
        self.p_cluster = None
        self.q_cluster = None

        ### n and m indices for iteration
        self.n_indices, self.m_indices = miepy.vsh.get_indices(self.Lmax)

        ### solve the interactions
        self.solve()

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
        Escat = miepy.vsh.expand_E(self.p_scat[i], self.q_scat[i], self.material_data.k,
                  mode=miepy.vsh.VSH_mode.outgoing, origin=self.position[i])(x,y,z)

        p = self.p_inc[i]
        q = self.q_inc[i]
        if not source:
            p -= self.p_src[i]
            q -= self.q_src[i]

        Einc = miepy.vsh.expand_E(p, q, self.material_data.k,
                  mode=miepy.vsh.VSH_mode.ingoing, origin=self.position[i])(x,y,z)

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
        Hscat = miepy.vsh.expand_H(self.p_scat[i], self.q_scat[i], self.material_data.k,
                  mode=miepy.vsh.VSH_mode.outgoing, eps_b=self.material_data.eps_b,
                  mu_b=self.material_data.mu_b, origin=self.position[i])(x,y,z)

        p = self.p_inc[i]
        q = self.q_inc[i]
        if not source:
            p -= self.p_src[i]
            q -= self.q_src[i]

        Hinc = miepy.vsh.expand_H(p, q, self.material_data.k,
                  mode=miepy.vsh.VSH_mode.ingoing, eps_b=self.material_data.eps_b,
                  mu_b=self.material_data.mu_b, origin=self.position[i])(x,y,z)

        return Hscat + Hinc
    
    def E_field(self, x, y, z, source=True):
        """Compute the electric field due to all particles
             
            Arguments:
                x         x position (array-like) 
                y         y position (array-like) 
                z         z position (array-like) 
                source    Include the source field (bool, default=True)

            Returns: E[3,...]
        """
        E = sum((miepy.vsh.expand_E(self.p_scat[i], self.q_scat[i], self.material_data.k,
                  mode=miepy.vsh.VSH_mode.outgoing, origin=self.position[i])(x,y,z)
                for i in range(self.Nparticles)))

        if source:
            E += self.source.E(np.array([x,y,z]), self.material_data.k)
        
        return E

    def H_field(self, x, y, z, source=True):
        """Compute the magnetic field due to all particles
             
            Arguments:
                x         x position (array-like) 
                y         y position (array-like) 
                z         z position (array-like) 
                source    Include the source field (bool, default=True)

            Returns: H[3,...]
        """
        H = sum((miepy.vsh.expand_H(self.p_scat[i], self.q_scat[i], self.material_data.k,
                  mode=miepy.vsh.VSH_mode.outgoing, eps_b=self.material_data.eps_b,
                  mu_b=self.material_data.mu_b, origin=self.position[i])(x,y,z)
                for i in range(self.Nparticles)))

        if source:
            factor = (self.material_data.eps_b/self.material_data.mu_b)**0.5
            H += factor*self.source.H(np.array([x,y,z]), self.material_data.k)
        
        return H

    def cross_sections_per_multipole(self, Lmax=None):
        """Compute the scattering, absorption, and extinction cross-section of the cluster per multipole

        Arguments:
            Lmax    (optional) compute scattering for up to Lmax terms (defult: self.Lmax)
        """

        if Lmax is None:
            Lmax = self.Lmax
        n_indices, m_indices = miepy.vsh.get_indices(Lmax)
        rmax = n_indices.shape[0]

        p0 =  np.zeros(rmax, dtype=complex)
        q0 =  np.zeros(rmax, dtype=complex)
        self.solve_cluster_coefficients(Lmax)
        for r in range(rmax):
            n = n_indices[r]
            m = m_indices[r]
            p0[r], q0[r] = self.source.structure_of_mode(n, m, self.origin, self.material_data.k)

        return miepy.flux.cluster_cross_sections(self.p_cluster, self.q_cluster,
                  p0, q0, self.material_data.k)


    def cross_sections(self):
        """Compute the scattering, absorption, and extinction cross-section of the cluster"""

        Cscat = 0
        Cabs = 0
        Cext = 0

        for i in range(self.Nparticles):
            C,A,E = self.cross_sections_of_particle(i)
            Cscat += C
            Cabs += A
            Cext += E

        return Cscat, Cabs, Cext

    def cross_sections_per_multipole_of_particle(self, i):
        """Compute the scattering, absorption, and extinction cross-section per multipole of a single particle

        Arguments:
            i    particle index
        """

        return miepy.flux.particle_cross_sections(self.p_scat[i], self.q_scat[i], self.p_src[i], self.q_src[i],
                    self.radius[i], self.material_data.k, self.material_data.n[i], self.material_data.mu[i],
                    self.material_data.n_b, self.material_data.mu_b)

    def cross_sections_of_particle(self, i):
        """Compute the scattering, absorption, and extinction cross-section of a single particle

        Arguments:
            i    particle index
        """

        Cscat, Cabs, Cext = self.cross_sections_per_multipole_of_particle(i)
        return map(lambda C: np.sum(C), [Cscat, Cabs, Cext])

    def force_on_particle(self, i, source=True):
        """Determine the force on a single particle

            Arguments:
                i         Particle index
                source    Include the source field (bool, default=True)
            
            Returns: F[3]
        """

        if source:
            p_inc = self.p_inc
            q_inc = self.q_inc
        else:
            p_inc = self.p_inc - self.p_src
            q_inc = self.q_inc - self.q_src

        F = miepy.forces.force(self.p_scat[i], self.q_scat[i], p_inc[i], q_inc[i],
                self.material_data.k, self.material_data.eps_b,
                self.material_data.mu_b)

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
            q_inc = self.q_inc
        else:
            p_inc = self.p_inc - self.p_src
            q_inc = self.q_inc - self.q_src

        T = miepy.forces.torque(self.p_scat[i], self.q_scat[i], p_inc[i], q_inc[i],
                self.material_data.k, self.material_data.eps_b,
                self.material_data.mu_b)

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

    def solve_cluster_coefficients(self, Lmax=None):
        """Solve for the p,q coefficients of the entire cluster around the origin

        Arguments:
            Lmax    (optional) compute scattering for up to Lmax terms (default: self.Lmax)
        """

        if Lmax is None:
            Lmax = self.Lmax

        self.p_cluster, self.q_cluster = miepy.vsh.cluster_coefficients(self.position, 
                self.p_scat, self.q_scat, self.material_data.k, origin=self.origin, Lmax=Lmax)

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
        self.q_cluster = None

    def _solve_source_decomposition(self):
        pos = self.position.T
        for r in range(self.rmax):
            for n in range(self.Nparticles):
                pos = self.position[n]
                self.p_src[n,r], self.q_src[n,r] = \
                    self.source.structure_of_mode(self.n_indices[r], self.m_indices[r], pos, self.material_data.k)

    def _solve_without_interactions(self):
        self.p_inc[...] = self.p_src
        self.q_inc[...] = self.q_src

        for r in range(self.rmax):
            n = self.n_indices[r]
            self.p_scat[...,r] = self.p_inc[...,r]*self.a[...,n-1]
            self.q_scat[...,r] = self.q_inc[...,r]*self.b[...,n-1]

    def _solve_interactions(self):
        sol = miepy.interactions.solve_sphere_cluster(self.position, self.a, self.b, 
                self.p_src, self.q_src, self.material_data.k)
        self.p_inc[...] = sol[0]
        self.q_inc[...] = sol[1]

        for r in range(self.rmax):
            n = self.n_indices[r]
            self.p_scat[:,r] = self.p_inc[:,r]*self.a[:,n-1]
            self.q_scat[:,r] = self.q_inc[:,r]*self.b[:,n-1]
