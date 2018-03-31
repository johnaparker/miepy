"""
The Generalized Mie Theory (GMT) for a collection of spheres.
"""
import numpy as np
import miepy
from my_pytools.my_numpy.array import atleast
from miepy.special_functions import riccati_1,riccati_2,vector_spherical_harmonics

def get_indices(Lmax):
    """return n_indices, m_indices arrays for a given Lmax"""
    rmax = Lmax*(Lmax + 2)
    n_indices = np.zeros(rmax, dtype=int)
    m_indices = np.zeros(rmax, dtype=int)

    counter = 0
    for n in range(1,Lmax+1):
        for m in range(-n,n+1):
            n_indices[counter] = n
            m_indices[counter] = m
            counter += 1

    return n_indices, m_indices

class material_struct:
    """struct to hold material data, and update it with changing wavelength"""
    def __init__(self, materials, medium, wavelength=None):
        """Arguments:
               materials    list of materials
               medium       medium material
               wavelength   wavelength (default: None)
        """
        Nparticles = len(materials)

        self.materials = materials
        self.medium = medium
        self._wavelength = wavelength

        self.eps   = np.zeros(Nparticles, dtype=complex) 
        self.mu    = np.zeros(Nparticles, dtype=complex) 
        self.n     = np.zeros(Nparticles, dtype=complex) 
        self.eps_b = None
        self.mu_b  = None
        self.n_b   = None
        self.k     = None

        self.wavelength = wavelength

    @property
    def wavelength(self):
        """get the wavelength"""
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        """set the wavelength to some value, changing all eps and mu data with it"""
        self._wavelength = value

        if value is not None:
            self.eps_b = self.medium.eps(self._wavelength)
            self.mu_b  = self.medium.mu(self._wavelength)
            self.n_b = (self.eps_b*self.mu_b)**0.5
            self.k = 2*np.pi*self.n_b/self._wavelength

            for i,material in enumerate(self.materials):
                self.eps[i] = material.eps(self._wavelength)
                self.mu[i]  = material.mu(self._wavelength)

            self.n[...]   = (self.eps*self.mu)**0.5

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
        self.material_data = material_struct(self.material, self.medium, wavelength=self.wavelength)

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
        self.n_indices, self.m_indices = get_indices(self.Lmax)

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
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        E = np.zeros((3,) + x.shape, dtype=complex)

        R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.position[i])
        rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

        E_sph = np.zeros((3,) + x.shape, dtype=complex)
        for r in range(self.rmax):
            n = self.n_indices[r]
            m = self.m_indices[r]
            factor = 1j*miepy.vsh.Emn(m, n)

            N,M = miepy.vsh.VSH(n,m)
            E_sph += factor*self.p_scat[i,r]*N(R, THETA, PHI, self.material_data.k)
            E_sph += factor*self.q_scat[i,r]*M(R, THETA, PHI, self.material_data.k)

            N,M = miepy.vsh.VSH(n, m, miepy.vsh.VSH_mode.incident)
            p = self.p_inc[i,r]
            q = self.q_inc[i,r]
            if not source:
                p -= self.p_src[i,r]
                q -= self.q_src[i,r]

            E_sph += -factor*p*N(R, THETA, PHI, self.material_data.k)
            E_sph += -factor*q*M(R, THETA, PHI, self.material_data.k)

        E += E_sph[0]*rhat + E_sph[1]*that + E_sph[2]*phat      # convert to cartesian
        
        return E

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
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        H = np.zeros((3,) + x.shape, dtype=complex)

        R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.position[i])
        rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

        H_sph = np.zeros((3,) + x.shape, dtype=complex)
        for r in range(self.rmax):
            n = self.n_indices[r]
            m = self.m_indices[r]
            factor = miepy.vsh.Emn(m, n)
            N,M = miepy.vsh.VSH(n,m)
            H_sph += factor*self.q_scat[i,r]*N(R, THETA, PHI, self.material_data.k)
            H_sph += factor*self.p_scat[i,r]*M(R, THETA, PHI, self.material_data.k)

            N,M = miepy.vsh.VSH(n, m, miepy.vsh.VSH_mode.incident)
            p = self.p_inc[i,r]
            q = self.q_inc[i,r]
            if not source:
                p -= self.p_src[i,r]
                q -= self.q_src[i,r]

            H_sph += -factor*q*N(R,THETA,PHI,self.material_data.k)
            H_sph += -factor*p*M(R,THETA,PHI,self.material_data.k)

        H += H_sph[0]*rhat + H_sph[1]*that + H_sph[2]*phat      # convert to cartesian
        
        return H*(self.material_data.eps_b/self.material_data.mu_b)**0.5
    
    #TODO: should source be computed through E func or through p/q src coefficients?
    def E_field(self, x, y, z, source=True):
        """Compute the electric field due to all particles
             
            Arguments:
                x         x position (array-like) 
                y         y position (array-like) 
                z         z position (array-like) 
                source    Include the source field (bool, default=True)

            Returns: E[3,...]
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        E = np.zeros((3,) + x.shape, dtype=complex)
        for i in range(self.Nparticles):
            R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.position[i])
            rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

            E_sph = np.zeros((3,) + x.shape, dtype=complex)
            for r in range(self.rmax):
                n = self.n_indices[r]
                m = self.m_indices[r]
                factor = 1j*miepy.vsh.Emn(m, n)
                N,M = miepy.vsh.VSH(n,m)
                E_sph += factor*self.p_scat[i,r]*N(R, THETA, PHI, self.material_data.k)
                E_sph += factor*self.q_scat[i,r]*M(R, THETA, PHI, self.material_data.k)

            E += E_sph[0]*rhat + E_sph[1]*that + E_sph[2]*phat      # convert to cartesian

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

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        H = np.zeros((3,) + x.shape, dtype=complex)
        for i in range(self.Nparticles):
            R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.position[i])
            rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

            H_sph = np.zeros((3,) + x.shape, dtype=complex)
            for r in range(self.rmax):
                n = self.n_indices[r]
                m = self.m_indices[r]
                factor = miepy.vsh.Emn(m, n)
                N,M = miepy.vsh.VSH(n,m)
                H_sph += factor*self.q_scat[i,r]*N(R,THETA,PHI,self.material_data.k)
                H_sph += factor*self.p_scat[i,r]*M(R,THETA,PHI,self.material_data.k)

            H += H_sph[0]*rhat + H_sph[1]*that + H_sph[2]*phat      # convert to cartesian

        if source:
            H += self.source.H(np.array([x,y,z]), self.material_data.k)

        return H*(self.material_data.eps_b/self.material_data.mu_b)**0.5

    def cross_sections_per_multipole(self, Lmax=None):
        """Compute the scattering, absorption, and extinction cross-section of the cluster per multipole

        Arguments:
            Lmax    (optional) compute scattering for up to Lmax terms (defult: self.Lmax)
        """

        if Lmax is None:
            Lmax = self.Lmax

        n_indices, m_indices = get_indices(Lmax)
        rmax = n_indices.shape[0]

        p0 =  np.zeros(rmax, dtype=complex)
        q0 =  np.zeros(rmax, dtype=complex)
        self.solve_cluster_coefficients(Lmax)

        Cscat = np.zeros([2, Lmax], dtype=float)
        Cext  = np.zeros([2, Lmax], dtype=float)

        factor = 4*np.pi/self.material_data.k**2
        for r in range(rmax):
            n = n_indices[r]
            m = m_indices[r]

            p0[r], q0[r] = self.source.structure_of_mode(n, m, self.origin, self.material_data.k)

            Cscat[0,n-1] += factor*np.abs(self.p_cluster[r])**2
            Cscat[1,n-1] += factor*np.abs(self.q_cluster[r])**2

            Cext[0,n-1] += factor*np.real(np.conj(p0[r])*self.p_cluster[r])
            Cext[1,n-1] += factor*np.real(np.conj(q0[r])*self.q_cluster[r])

        Cabs = Cext - Cscat
        return Cscat, Cabs, Cext

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

        Cscat = np.zeros([2, self.Lmax], dtype=float)
        Cext  = np.zeros([2, self.Lmax], dtype=float)
        Cabs  = np.zeros([2, self.Lmax], dtype=float)

        riccati = miepy.special_functions.riccati_1_single

        factor = 4*np.pi/self.material_data.k**2

        xj = self.material_data.k*self.radius[i]
        mj = self.material_data.n[i]/self.material_data.n_b
        yj = xj*mj
        mu_b = self.material_data.mu_b
        mu = self.material_data.mu[i]

        for r in range(self.rmax):
            n = self.n_indices[r]

            # Cscat[0,n-1] += factor*np.abs(self.p_scat[i,r])**2
            # Cscat[1,n-1] += factor*np.abs(self.q_scat[i,r])**2

            psi_x, psi_xp = riccati(n, xj)
            psi_y, psi_yp = riccati(n, yj)

            
            Dn = -np.divide(np.real(1j*mj*mu_b*mu*psi_y*np.conj(psi_yp)),
                            np.abs(mu_b*mj*psi_y*psi_xp - mu*psi_x*psi_yp)**2)
            Cn = -np.divide(np.real(1j*np.conj(mj)*mu_b*mu*psi_y*np.conj(psi_yp)),
                            np.abs(mu*psi_y*psi_xp - mu_b*mj*psi_x*psi_yp)**2)

            Cabs[0,n-1] += Dn*factor*np.abs(self.p_scat[i,r])**2
            Cabs[1,n-1] += Cn*factor*np.abs(self.q_scat[i,r])**2

            #TODO should this be p_src or p_inc?
            Cext[0,n-1] += factor*np.real(np.conj(self.p_src[i,r])*self.p_scat[i,r])
            Cext[1,n-1] += factor*np.real(np.conj(self.q_src[i,r])*self.q_scat[i,r])

        Cscat[...] = Cext - Cabs
        return Cscat, Cabs, Cext

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
                self.material_data.mu_b, self.Lmax)

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
                self.material_data.mu_b, self.Lmax)

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

        n_indices, m_indices = get_indices(Lmax)
        rmax = n_indices.shape[0]

        self.p_cluster = np.zeros(rmax, dtype=complex)
        self.q_cluster = np.zeros(rmax, dtype=complex)

        for i in range(self.Nparticles):
            if np.all(self.position[i] == self.origin):
                self.p_cluster[...] = self.p_scat[:,i]
                self.q_cluster[...] = self.q_scat[:,i]
                continue

            rij = self.origin - self.position[i]
            rad, theta, phi = miepy.coordinates.cart_to_sph(*rij)
            
            for r in range(rmax):
                n = n_indices[r]
                m = m_indices[r]

                for rp in range(self.rmax):
                    v = self.n_indices[rp]
                    u = self.m_indices[rp]

                    a = self.p_scat[i,rp]
                    b = self.q_scat[i,rp]

                    A = miepy.vsh.A_translation(m, n, u, v, rad, theta, phi, self.material_data.k, miepy.vsh.VSH_mode.incident)
                    B = miepy.vsh.B_translation(m, n, u, v, rad, theta, phi, self.material_data.k, miepy.vsh.VSH_mode.incident)

                    self.p_cluster[r] += a*A + b*B
                    self.q_cluster[r] += a*B + b*A

    def solve(self, wavelength=None, source=None):
        """solve for the p,q incident and scattering coefficients

           Arguments:
               wavelength   wavelength to solve at (default: current wavelength)
               source       source to use (default: current source). If current source is also None, solve the particle's T-matrix instead
        """
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
        self._solve_source_decomposition()
        self.p_inc[...] = self.p_src
        self.q_inc[...] = self.q_src

        for r in range(self.rmax):
            n = self.n_indices[r]
            self.p_scat[...,r] = self.p_inc[...,r]*self.a[...,n-1]
            self.q_scat[...,r] = self.q_inc[...,r]*self.b[...,n-1]

    #TODO vectorize for loops. Avoid transpose of position->pass x,y,z to source instead...?
    def _solve_interactions(self):
        self._solve_source_decomposition()
        identity = np.zeros(shape = (2, self.Nparticles, self.rmax, 2, self.Nparticles, self.rmax), dtype=np.complex)
        np.einsum('airair->air', identity)[...] = 1
        
        interaction_matrix = np.zeros(shape = (2, self.Nparticles, self.rmax, 2, self.Nparticles, self.rmax), dtype=np.complex)

        for i in range(self.Nparticles):
            for j in range(self.Nparticles):
                if i == j: continue

                pi = self.position[i]
                pj = self.position[j]
                dji = pi -  pj
                r_ji = np.linalg.norm(dji)
                theta_ji = np.arccos(dji[2]/r_ji)
                phi_ji = np.arctan2(dji[1], dji[0])

                for r in range(self.rmax):
                    n = self.n_indices[r]
                    m = self.m_indices[r]
                    for s in range(self.rmax):
                        v = self.n_indices[s]
                        u = self.m_indices[s]

                        A_transfer = miepy.vsh.A_translation(m,n,u,v,r_ji,theta_ji,phi_ji,self.material_data.k, miepy.vsh.VSH_mode.outgoing)
                        B_transfer = miepy.vsh.B_translation(m,n,u,v,r_ji,theta_ji,phi_ji,self.material_data.k, miepy.vsh.VSH_mode.outgoing)

                        interaction_matrix[0,i,r,0,j,s] = A_transfer*self.a[j,v-1]
                        interaction_matrix[0,i,r,1,j,s] = B_transfer*self.b[j,v-1]
                        interaction_matrix[1,i,r,0,j,s] = B_transfer*self.a[j,v-1]
                        interaction_matrix[1,i,r,1,j,s] = A_transfer*self.b[j,v-1]

        A = identity + interaction_matrix
        b = np.array([self.p_src,self.q_src])
        sol = np.linalg.tensorsolve(A, b)
        self.p_inc[...] = sol[0]
        self.q_inc[...] = sol[1]

        for r in range(self.rmax):
            n = self.n_indices[r]
            self.p_scat[:,r] = self.p_inc[:,r]*self.a[:,n-1]
            self.q_scat[:,r] = self.q_inc[:,r]*self.b[:,n-1]
