"""
The Generalized Mie Theory (GMT) for a collection of spheres.
"""
import numpy as np
import miepy
from my_pytools.my_numpy.array import atleast
from collections import namedtuple
from miepy.special_functions import riccati_1,riccati_2,vector_spherical_harmonics

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

class material_struct:
    def __init__(self, materials, medium, wavelength=None):
        """struct to hold material data, and update it with changing wavelength
            
           Arguments:
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

#TODO: make several properties... such as wavelength, source, position, etc.
#TODO: prefer manual solve calls ratehr than auto solve calls (or have an auto_solve option)
class gmt:
    """Solve Generalized Mie Theory: N particles in an arbitray source profile"""
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
        self.spheres = spheres
        self.source = source
        self.wavelength = np.asarray(np.atleast_1d(wavelength), dtype=float)
        self.Lmax = Lmax
        self.rmax = Lmax*(Lmax + 2)
        self.interactions = interactions

        self.auto_origin = False    
        if origin is None:
            self.origin = np.zeros(3)
        elif origin is 'auto':
            self.auto_origin = True
            self.origin = np.average(self.spheres.position, axis=0)
        else:
            self.origin = np.asarray(origin)

        if medium is None:
            self.medium = miepy.constant_material(1.0, 1.0)
        else:
            self.medium = medium
            if (self.medium.eps(self.wavelength).imag != 0).any()  \
                         or (self.medium.mu(self.wavelength).imag != 0).any():
                raise ValueError('medium must be non-absorbing')

        self.Nfreq = len(self.wavelength)
        self.Nparticles = len(spheres)

        #TODO this should not be a dictionary -- a struct is better since the number of items is known
        self.material_data = {}
        self.material_data['wavelength'] = self.wavelength
        self.material_data['eps']        = np.zeros([self.Nparticles, self.Nfreq], dtype=complex) 
        self.material_data['mu']         = np.zeros([self.Nparticles, self.Nfreq], dtype=complex) 
        self.material_data['n']          = np.zeros([self.Nparticles, self.Nfreq], dtype=complex) 
        self.material_data['eps_b']      = self.medium.eps(self.wavelength)
        self.material_data['mu_b']       = self.medium.mu(self.wavelength)
        self.material_data['n_b']        = np.sqrt(self.material_data['eps_b']*self.material_data['mu_b'])
        self.material_data['k']          = 2*np.pi*self.material_data['n_b']/self.wavelength

        self.a = np.zeros([self.Nfreq,self.Nparticles,self.Lmax], dtype=complex)
        self.b = np.zeros([self.Nfreq,self.Nparticles,self.Lmax], dtype=complex)

        self.p_scat = np.zeros([self.Nfreq,self.Nparticles,self.rmax], dtype=complex)
        self.q_scat = np.zeros([self.Nfreq,self.Nparticles,self.rmax], dtype=complex)
        self.p_inc  = np.zeros([self.Nfreq,self.Nparticles,self.rmax], dtype=complex)
        self.q_inc  = np.zeros([self.Nfreq,self.Nparticles,self.rmax], dtype=complex)
        self.p_src  = np.zeros([self.Nfreq,self.Nparticles,self.rmax], dtype=complex)
        self.q_src  = np.zeros([self.Nfreq,self.Nparticles,self.rmax], dtype=complex)

        self.p_cluster = None
        self.q_cluster = None

        for i in range(self.Nparticles):
            sphere = miepy.single_mie_sphere(self.spheres.radius[i], self.spheres.material[i],
                        self.wavelength, self.Lmax, self.medium)
            self.material_data['eps'][i] = sphere.material_data['eps']
            self.material_data['mu'][i] = sphere.material_data['mu']
            self.material_data['n'][i] = sphere.material_data['n']
            self.a[:,i], self.b[:,i] = sphere.solve_exterior()

        self.n_indices, self.m_indices = get_indices(self.Lmax)

        self.solve()

    def E_field_from_particle(self, i, x, y, z, source=True):
        """Compute the electric field around particle i
             
            Arguments:
                i        particle number
                x        x position (array-like) 
                y        y position (array-like) 
                z        z position (array-like) 
                source   Include the source field (bool, default=True)

            Returns: E[3,M], M = number of wavelengths
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        E = np.zeros((3,self.Nfreq) + x.shape, dtype=complex)

        R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.spheres.position[i])
        rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

        for k in range(self.Nfreq):
            E_sph = np.zeros((3,) + x.shape, dtype=complex)
            for r in range(self.rmax):
                n = self.n_indices[r]
                m = self.m_indices[r]
                factor = 1j*miepy.vsh.Emn(m, n)

                N,M = miepy.vsh.VSH(n,m)
                E_sph += factor*self.p_scat[k,i,r]*N(R,THETA,PHI,self.material_data['k'][k])
                E_sph += factor*self.q_scat[k,i,r]*M(R,THETA,PHI,self.material_data['k'][k])

                N,M = miepy.vsh.VSH(n, m, miepy.vsh.VSH_mode.incident)
                p = self.p_inc[k,i,r]
                q = self.q_inc[k,i,r]
                if not source:
                    p -= self.p_src[k,i,r]
                    q -= self.q_src[k,i,r]

                E_sph += -factor*p*N(R,THETA,PHI,self.material_data['k'][k])
                E_sph += -factor*q*M(R,THETA,PHI,self.material_data['k'][k])

            E[:,k] += E_sph[0]*rhat + E_sph[1]*that + E_sph[2]*phat      # convert to cartesian
        
        return E

    def H_field_from_particle(self, i, x, y, z, source=True):
        """Compute the magnetic field around particle i
             
            Arguments:
                i        particle number
                x        x position (array-like) 
                y        y position (array-like) 
                z        z position (array-like) 
                source   Include the source field (bool, default=True)

            Returns: H[3,M], M = number of wavelengths
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        H = np.zeros((3,self.Nfreq) + x.shape, dtype=complex)

        R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.spheres.position[i])
        rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

        for k in range(self.Nfreq):
            H_sph = np.zeros((3,) + x.shape, dtype=complex)
            for r in range(self.rmax):
                n = self.n_indices[r]
                m = self.m_indices[r]
                factor = miepy.vsh.Emn(m, n)
                N,M = miepy.vsh.VSH(n,m)
                H_sph += factor*self.q_scat[k,i,r]*N(R,THETA,PHI,self.material_data['k'][k])
                H_sph += factor*self.p_scat[k,i,r]*M(R,THETA,PHI,self.material_data['k'][k])

                N,M = miepy.vsh.VSH(n, m, miepy.vsh.VSH_mode.incident)
                p = self.p_inc[k,i,r]
                q = self.q_inc[k,i,r]
                if not source:
                    p -= self.p_src[k,i,r]
                    q -= self.q_src[k,i,r]

                H_sph += -factor*q*N(R,THETA,PHI,self.material_data['k'][k])
                H_sph += -factor*p*M(R,THETA,PHI,self.material_data['k'][k])

            H[:,k] += H_sph[0]*rhat + H_sph[1]*that + H_sph[2]*phat      # convert to cartesian
        
        return H*(self.material_data['eps_b'][k]/self.material_data['mu_b'][k])**0.5
    
    #TODO: should source be computed through E func or through p/q src coefficients?
    def E_field(self, x, y, z, source=True):
        """Compute the electric field due to all particles
             
            Arguments:
                x         x position (array-like) 
                y         y position (array-like) 
                z         z position (array-like) 
                source    Include the source field (bool, default=True)

            Returns: E[3,M], M = number of wavelengths
        """
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        E = np.zeros((3,self.Nfreq) + x.shape, dtype=complex)
        for i in range(self.Nparticles):
            R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.spheres.position[i])
            rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

            for k in range(self.Nfreq):
                E_sph = np.zeros((3,) + x.shape, dtype=complex)
                for r in range(self.rmax):
                    n = self.n_indices[r]
                    m = self.m_indices[r]
                    factor = 1j*miepy.vsh.Emn(m, n)
                    N,M = miepy.vsh.VSH(n,m)
                    E_sph += factor*self.p_scat[k,i,r]*N(R,THETA,PHI,self.material_data['k'][k])
                    E_sph += factor*self.q_scat[k,i,r]*M(R,THETA,PHI,self.material_data['k'][k])

                E[:,k] += E_sph[0]*rhat + E_sph[1]*that + E_sph[2]*phat      # convert to cartesian

        if source:
            for k in range(self.Nfreq):
                E[:,k] += self.source.E(np.array([x,y,z]), self.material_data['k'][k])
        
        return E

    def H_field(self, x, y, z, source=True):
        """Compute the magnetic field due to all particles
             
            Arguments:
                x         x position (array-like) 
                y         y position (array-like) 
                z         z position (array-like) 
                source    Include the source field (bool, default=True)

            Returns: H[3,M], M = number of wavelengths
        """

        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)

        H = np.zeros((3,self.Nfreq) + x.shape, dtype=complex)
        for i in range(self.Nparticles):
            R, THETA, PHI = miepy.coordinates.cart_to_sph(x, y, z, self.spheres.position[i])
            rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)

            for k in range(self.Nfreq):
                H_sph = np.zeros((3,) + x.shape, dtype=complex)
                for r in range(self.rmax):
                    n = self.n_indices[r]
                    m = self.m_indices[r]
                    factor = miepy.vsh.Emn(m, n)
                    N,M = miepy.vsh.VSH(n,m)
                    H_sph += factor*self.q_scat[k,i,r]*N(R,THETA,PHI,self.material_data['k'][k])
                    H_sph += factor*self.p_scat[k,i,r]*M(R,THETA,PHI,self.material_data['k'][k])

                H[:,k] += H_sph[0]*rhat + H_sph[1]*that + H_sph[2]*phat      # convert to cartesian

        if source:
            for k in range(self.Nfreq):
                H[:,k] += self.source.H(np.array([x,y,z]), self.material_data['k'][k])

        return H*(self.material_data['eps_b'][k]/self.material_data['mu_b'][k])**0.5

    def cross_sections_per_multipole(self, Lmax=None):
        """Compute the scattering, absorption, and extinction cross-section of the cluster per multipole

        Arguments:
            Lmax    (optional) compute scattering for up to Lmax terms (defult: self.Lmax)
        """

        if Lmax is None:
            Lmax = self.Lmax

        n_indices, m_indices = get_indices(Lmax)
        rmax = n_indices.shape[0]

        p0 =  np.zeros([rmax, self.Nfreq], dtype=complex)
        q0 =  np.zeros([rmax, self.Nfreq], dtype=complex)
        self.solve_cluster_coefficients(Lmax)

        Cscat = np.zeros([2, Lmax, self.Nfreq], dtype=float)
        Cext  = np.zeros([2, Lmax, self.Nfreq], dtype=float)

        for k in range(self.Nfreq):
            factor = 4*np.pi/self.material_data['k'][k]**2
            for r in range(rmax):
                n = n_indices[r]
                m = m_indices[r]

                p0[r,k], q0[r,k] = self.source.structure_of_mode(n, m, self.origin, self.material_data['k'][k])

                Cscat[0,n-1,k] += factor*np.abs(self.p_cluster[k,r])**2
                Cscat[1,n-1,k] += factor*np.abs(self.q_cluster[k,r])**2

                Cext[0,n-1,k] += factor*np.real(np.conj(p0[r,k])*self.p_cluster[k,r])
                Cext[1,n-1,k] += factor*np.real(np.conj(q0[r,k])*self.q_cluster[k,r])

        Cabs = Cext - Cscat
        return Cscat, Cabs, Cext

    def cross_sections(self):
        """Compute the scattering, absorption, and extinction cross-section of the cluster"""

        Cscat = np.zeros(self.Nfreq, dtype=float)
        Cabs = np.zeros(self.Nfreq, dtype=float)
        Cext = np.zeros(self.Nfreq, dtype=float)

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

        Cscat = np.zeros([2, self.Lmax, self.Nfreq], dtype=float)
        Cext  = np.zeros([2, self.Lmax, self.Nfreq], dtype=float)
        Cabs  = np.zeros([2, self.Lmax, self.Nfreq], dtype=float)

        riccati = miepy.special_functions.riccati_1_single

        for k in range(self.Nfreq):
            factor = 4*np.pi/self.material_data['k'][k]**2

            xj = self.material_data['k'][k]*self.spheres.radius[i]
            mj = self.material_data['n'][i,k]/self.material_data['n_b'][k]
            yj = xj*mj
            mu_b = self.material_data['mu_b'][k]
            mu = self.material_data['mu'][i,k]

            for r in range(self.rmax):
                n = self.n_indices[r]

                # Cscat[0,n-1,k] += factor*np.abs(self.p_scat[k,i,r])**2
                # Cscat[1,n-1,k] += factor*np.abs(self.q_scat[k,i,r])**2

                psi_x, psi_xp = riccati(n, xj)
                psi_y, psi_yp = riccati(n, yj)

                
                Dn = -np.divide(np.real(1j*mj*mu_b*mu*psi_y*np.conj(psi_yp)),
                                np.abs(mu_b*mj*psi_y*psi_xp - mu*psi_x*psi_yp)**2)
                Cn = -np.divide(np.real(1j*np.conj(mj)*mu_b*mu*psi_y*np.conj(psi_yp)),
                                np.abs(mu*psi_y*psi_xp - mu_b*mj*psi_x*psi_yp)**2)

                Cabs[0,n-1,k] += Dn*factor*np.abs(self.p_scat[k,i,r])**2
                Cabs[1,n-1,k] += Cn*factor*np.abs(self.q_scat[k,i,r])**2

                #TODO should this be p_src or p_inc?
                Cext[0,n-1,k] += factor*np.real(np.conj(self.p_src[k,i,r])*self.p_scat[k,i,r])
                Cext[1,n-1,k] += factor*np.real(np.conj(self.q_src[k,i,r])*self.q_scat[k,i,r])

        Cscat[...] = Cext - Cabs
        return Cscat, Cabs, Cext

    def cross_sections_of_particle(self, i):
        """Compute the scattering, absorption, and extinction cross-section of a single particle

        Arguments:
            i    particle index
        """

        Cscat, Cabs, Cext = self.cross_sections_per_multipole_of_particle(i)
        return map(lambda C: np.sum(C, axis=(0,1)), [Cscat, Cabs, Cext])

    def force_on_particle(self, i, source=True):
        """Determine the force on a single particle

            Arguments:
                i         Particle index
                source    Include the source field (bool, default=True)
            
            Returns: F[3,M], M = number of wavelengths
        """

        F = np.zeros((3,self.Nfreq))

        if source:
            p_inc = self.p_inc
            q_inc = self.q_inc
        else:
            p_inc = self.p_inc - self.p_src
            q_inc = self.q_inc - self.q_src

        for k in range(self.Nfreq):
            F[:,k] = miepy.forces.force(self.p_scat[k,i], self.q_scat[k,i], p_inc[k,i], q_inc[k,i],
                        self.material_data['k'][k], self.material_data['eps_b'][k],
                        self.material_data['mu_b'][k], self.Lmax)

        return F

    def torque_on_particle(self, i, source=True):
        """Determine the torque on a single particle

            Arguments:
                i         Particle index
                source    Include the source field (bool, default=True)
            
            Returns: T[3,M], M = number of wavelengths
        """

        T = np.zeros((3,self.Nfreq))

        if source:
            p_inc = self.p_inc
            q_inc = self.q_inc
        else:
            p_inc = self.p_inc - self.p_src
            q_inc = self.q_inc - self.q_src

        for k in range(self.Nfreq):
            T[:,k] = miepy.forces.torque(self.p_scat[k,i], self.q_scat[k,i], p_inc[k,i], q_inc[k,i],
                        self.material_data['k'][k], self.material_data['eps_b'][k],
                        self.material_data['mu_b'][k], self.Lmax)

        return T

    def force(self, source=True):
        """Determine the force on every particle

            Arguments:
                source    Include the source field (bool, default=False)
            
            Returns: F[3,M,N]
                     N = number of particles, M = number of wavelengths
        """
        F = np.zeros([3, self.Nfreq, self.Nparticles], dtype=float)
        for i in range(self.Nparticles):
            F[...,i] = self.force_on_particle(i, source=source)

        return F

    def torque(self, source=True):
        """Determine the torque on every particle

            Arguments:
                source    Include the source field (bool, default=False)
            
            Returns: T[3,M,N]
                     N = number of particles, M = number of wavelengths
        """
        T = np.zeros([3, self.Nfreq, self.Nparticles], dtype=float)
        for i in range(self.Nparticles):
            T[...,i] = self.torque_on_particle(i, source=source)

        return T

    def update_position(self, position):
        """Update the positions of the spheres

            Arguments
                position[N,3]       new particle positions
        """
        self.spheres.position = np.asarray(np.atleast_2d(position), dtype=float)
        self._reset_cluster_coefficients()

        if self.auto_origin:
            self.origin = np.average(self.spheres.position, axis=0)

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

        self.p_cluster = np.zeros([self.Nfreq, rmax], dtype=complex)
        self.q_cluster = np.zeros([self.Nfreq, rmax], dtype=complex)

        for i in range(self.Nparticles):
            if np.all(self.spheres.position[i] == self.origin):
                self.p_cluster[...] = self.p_scat[:,i]
                self.q_cluster[...] = self.q_scat[:,i]
                continue

            rij = self.origin - self.spheres.position[i]
            rad, theta, phi = miepy.coordinates.cart_to_sph(*rij)
            
            for k in range(self.Nfreq):
                for r in range(rmax):
                    n = n_indices[r]
                    m = m_indices[r]

                    for rp in range(self.rmax):
                        v = self.n_indices[rp]
                        u = self.m_indices[rp]

                        a = self.p_scat[k,i,rp]
                        b = self.q_scat[k,i,rp]

                        A = miepy.vsh.A_translation(m, n, u, v, rad, theta, phi, self.material_data['k'][k], miepy.vsh.VSH_mode.incident)
                        B = miepy.vsh.B_translation(m, n, u, v, rad, theta, phi, self.material_data['k'][k], miepy.vsh.VSH_mode.incident)

                        self.p_cluster[k,r] += a*A + b*B
                        self.q_cluster[k,r] += a*B + b*A

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
        pos = self.spheres.position.T
        for k in range(self.Nfreq):
            for r in range(self.rmax):
                for n in range(self.Nparticles):
                    pos = self.spheres.position[n]
                    self.p_src[k,n,r], self.q_src[k,n,r] = \
                        self.source.structure_of_mode(self.n_indices[r], self.m_indices[r], pos, self.material_data['k'][k])

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
        
        for k in range(self.Nfreq):
            interaction_matrix = np.zeros(shape = (2, self.Nparticles, self.rmax, 2, self.Nparticles, self.rmax), dtype=np.complex)

            for i in range(self.Nparticles):
                for j in range(self.Nparticles):
                    if i == j: continue

                    pi = self.spheres.position[i]
                    pj = self.spheres.position[j]
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

                            A_transfer = miepy.vsh.A_translation(m,n,u,v,r_ji,theta_ji,phi_ji,self.material_data['k'][k], miepy.vsh.VSH_mode.outgoing)
                            B_transfer = miepy.vsh.B_translation(m,n,u,v,r_ji,theta_ji,phi_ji,self.material_data['k'][k], miepy.vsh.VSH_mode.outgoing)

                            interaction_matrix[0,i,r,0,j,s] = A_transfer*self.a[k,j,v-1]
                            interaction_matrix[0,i,r,1,j,s] = B_transfer*self.b[k,j,v-1]
                            interaction_matrix[1,i,r,0,j,s] = B_transfer*self.a[k,j,v-1]
                            interaction_matrix[1,i,r,1,j,s] = A_transfer*self.b[k,j,v-1]

            A = identity + interaction_matrix
            b = np.array([self.p_src[k],self.q_src[k]])
            sol = np.linalg.tensorsolve(A, b)
            self.p_inc[k] = sol[0]
            self.q_inc[k] = sol[1]

            for r in range(self.rmax):
                n = self.n_indices[r]
                self.p_scat[k,:,r] = self.p_inc[k,:,r]*self.a[k,:,n-1]
                self.q_scat[k,:,r] = self.q_inc[k,:,r]*self.b[k,:,n-1]
