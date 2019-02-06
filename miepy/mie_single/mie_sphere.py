"""
mie_sphere calculates the scattering coefficients of a sphere using Mie theory
"""
import numpy as np
import miepy
from miepy.special_functions import (riccati_1, riccati_2, riccati_1_single,
        riccati_2_single, riccati_3_single, vector_spherical_harmonics)
from miepy.mie_single.scattering import scattered_E,scattered_H,interior_E,interior_H
from scipy import constants

def mie_sphere_scattering_coefficients(radius, n, eps, mu, eps_b, mu_b, k, conducting=False):
    """Solve for the exterior of the sphere, the an and bn coefficients
        
       Arguments:
           radius      sphere radius
           n           coefficient order (n=1,2,...)
           eps         sphere permitivitty
           mu          sphere permeability
           eps_b       medium permitivitty
           mu_b        medium permeability
           k           medium wavenumber
           conducting  if True, calculate for conducting sphere (default: False)
    """
    xvals = k*radius
    m = (eps/eps_b)**0.5
    jn = riccati_1_single(n, xvals)
    yn = riccati_3_single(n, xvals)

    if conducting:
        a = jn[1]/yn[1]
        b = jn[0]/yn[0]
    else:
        jnm = riccati_1_single(n, m*xvals)
        m = (eps/eps_b)**0.5
        mt = m*mu_b/mu

        a = (mt*jnm[0]*jn[1] - jn[0]*jnm[1])/(mt*jnm[0]*yn[1] - yn[0]*jnm[1])
        b = (jnm[0]*jn[1] - mt*jn[0]*jnm[1])/(jnm[0]*yn[1] - mt*yn[0]*jnm[1])

    return a, b

def mie_sphere_interior_coefficients(radius, n, eps, mu, eps_b, mu_b, k, conducting=False):
    """Solve for the exterior of the sphere, the an and bn coefficients
        
       Arguments:
           radius      sphere radius
           n           coefficient order (n=1,2,...)
           eps         sphere permitivitty
           mu          sphere permeability
           eps_b       medium permitivitty
           mu_b        medium permeability
           k           medium wavenumber
           conducting  if True, calculate for conducting sphere (default: False)
    """
    xvals = k*radius
    m = (eps/eps_b)**0.5
    mt = m*mu_b/mu

    jn = riccati_1_single(n, xvals)
    yn = riccati_3_single(n, xvals)
    jnm = riccati_1_single(n, m*xvals)

    if conducting:
        c = 0.0
        d = 0.0
    else:
        c = (m*jn[0]*yn[1] - m*yn[0]*jn[1])/(jnm[0]*yn[1] - mt*yn[0]*jnm[1])
        d = (m*jn[0]*yn[1] - m*yn[0]*jn[1])/(mt*jnm[0]*yn[1] - yn[0]*jnm[1])

    return c, d

class single_mie_sphere:
    def __init__(self, radius, material, wavelength, lmax, medium=None):
        """Solve traditional Mie theory: a single sphere in x-polarized plane wave illumination

               radius           particle radius
               material         particle material
               wavelength[N]    wavelength(s) to solve the system at
               lmax             maximum number of orders to use in angular momentum expansion
               medium           material medium (must be non-absorbing; defaults to vacuum)
        """

        self.radius = radius
        self.material = material
        self.wavelength = np.asarray(np.atleast_1d(wavelength), dtype=float)
        self.lmax = lmax
        if medium is None:
            self.medium = miepy.constant_material(1.0, 1.0)
        else:
            self.medium = medium
            if (self.medium.eps(self.wavelength).imag != 0).any()  \
                         or (self.medium.mu(self.wavelength).imag != 0).any():
                raise ValueError('medium must be non-absorbing')

        self.Nfreq = len(self.wavelength)

        self.material_data = {}
        self.material_data['wavelength'] = self.wavelength
        self.material_data['eps']        = self.material.eps(self.wavelength)
        self.material_data['mu']         = self.material.mu(self.wavelength)
        self.material_data['n']          = np.sqrt(self.material_data['eps']*self.material_data['mu'])
        self.material_data['eps_b']      = self.medium.eps(self.wavelength)
        self.material_data['mu_b']       = self.medium.mu(self.wavelength)
        self.material_data['n_b']        = np.sqrt(self.material_data['eps_b']*self.material_data['mu_b'])
        self.material_data['k']          = 2*np.pi*self.material_data['n_b']/self.wavelength
               
        #TODO (performance) swap an/bn shape for lmax by Nfreq, remove transpose
        self.an = np.zeros((self.Nfreq, self.lmax), dtype=np.complex)
        self.bn = np.zeros((self.Nfreq, self.lmax), dtype=np.complex)
        self.cn = np.zeros((self.Nfreq, self.lmax), dtype=np.complex)
        self.dn = np.zeros((self.Nfreq, self.lmax), dtype=np.complex)

        self.scattering_properties = (self.an, self.bn, self.material_data['k'])

        self.exterior_computed = False
        self.interior_computed = False
    
    def solve_exterior(self):
        """solve for the exterior of the sphere, the an and bn coefficients"""
        xvals = self.material_data['k']*self.radius
        m = (self.material_data['eps']/self.material_data['eps_b'])**.5
        mt = m*self.material_data['mu_b']/self.material_data['mu']

        jn = riccati_1(self.lmax,xvals)
        yn = riccati_2(self.lmax,xvals)

        if self.material.name == 'metal':
            a = jn[1]/yn[1]
            b = jn[0]/yn[0]
        else:
            jnm = riccati_1(self.lmax,m*xvals)
            a = (mt*jnm[0]*jn[1] - jn[0]*jnm[1])/(mt*jnm[0]*yn[1] - yn[0]*jnm[1])
            b = (jnm[0]*jn[1] - mt*jn[0]*jnm[1])/(jnm[0]*yn[1] - mt*yn[0]*jnm[1])

        self.an[...] = np.nan_to_num(a.T)
        self.bn[...] = np.nan_to_num(b.T)

        self.exterior_computed = True
        return self.an, self.bn

    def solve_interior(self):
        """solve for the interior of the sphere, the cn and dn coefficients"""
        #TODO shouldn't this be multiplied by self.radius?
        xvals = self.material_data['k']*self.radius
        m = (self.material_data['eps']/self.material_data['eps_b'])**.5
        mt = m*self.material_data['mu_b']/self.material_data['mu']

        jn = riccati_1(self.lmax,xvals)
        yn = riccati_2(self.lmax,xvals)

        if self.material.name == 'metal':
            f = self.material_data['mu_b']/self.material_data['mu']
            c = (jn[0]*yn[1] - yn[0]*jn[1])/(f*yn[0]*jnm[1])
            d = (jn[0]*yn[1] - yn[0]*jn[1])/(f*jnm[0]*yn[1])
        else:
            jnm = riccati_1(self.lmax,m*xvals)
            c = (m*jn[0]*yn[1] - m*yn[0]*jn[1])/(jnm[0]*yn[1] - mt*yn[0]*jnm[1])
            d = (m*jn[0]*yn[1] - m*yn[0]*jn[1])/(mt*jnm[0]*yn[1] - yn[0]*jnm[1])

        self.cn[...] = np.nan_to_num(c.T)
        self.dn[...] = np.nan_to_num(d.T)

        self.interior_computed = True
        return self.cn, self.dn

    def cross_sections(self):
        """Return the 3 cross-sections: (Scattering, Absorbption, Extinction)"""
        if not self.exterior_computed: self.solve_exterior()
        return miepy.cross_sections(*self.scattering_properties)

    def cross_sections_per_multipole(self):
        """Return the 3 cross-sections per multipole: (Scattering, Absorbption, Extinction)"""
        if not self.exterior_computed: self.solve_exterior()
        return miepy.flux.cross_sections(miepy.scattering_per_multipole(*self.scattering_properties),
                miepy.absorbption_per_multipole(*self.scattering_properties),
                miepy.extinction_per_multipole(*self.scattering_properties))

    #TODO: moved to functions script
    def radiation_force(self):
        """Return the radiation force (Fz)"""

        Fz = np.zeros(self.Nfreq)
        nvals = np.arange(1, self.lmax+1)

        k = self.material_data['k']
        for i in range(self.lmax):
            n = 1 + i

            factor = 2*np.pi*(2*n+1)/k**2
            Fz += factor*np.real(self.an[:,i] + self.bn[:,i])

            factor = 4*np.pi*(2*n+1)/(n*(n+1))/k**2
            Fz -= factor*np.real(self.an[:,i]*np.conj(self.bn[:,i]))

            if n < self.lmax:
                factor = 4*np.pi*(n*(n+2))/(n+1)/k**2
                Fz -= factor*np.real(self.an[:,i]*np.conj(self.an[:,i+1])
                                   + self.bn[:,i]*np.conj(self.bn[:,i+1]))

        return Fz * constants.epsilon_0*self.material_data['n_b']/2

    def E_field(self, index=None, lmax=None):
        """Return an electric field function E(r,theta,phi) for a given wavenumber index"""
        if lmax is None: lmax = self.lmax
        if index is None:
            index = np.s_[:]

        if not self.interior_computed: self.solve_interior()
        if not self.exterior_computed: self.solve_exterior()

        an = self.an[index, :lmax+1]
        bn = self.bn[index, :lmax+1]
        cn = self.cn[index, :lmax+1]
        dn = self.dn[index, :lmax+1]
        mat = self.material_data

        def E_func(r, theta, phi):
            E = np.zeros(shape = [3] + list(r.shape), dtype=np.complex)
            id_inside = r <= self.radius

            k = mat['k'][index]
            E[:,~id_inside] = scattered_E(an, bn, k)(r[~id_inside], theta[~id_inside], phi[~id_inside])

            k = 2*np.pi*mat['n'][index]/self.wavelength[index]
            E[:,id_inside] = interior_E(cn, dn, k)(r[id_inside], theta[id_inside], phi[id_inside])
            return E

        return E_func

    def H_field(self, index=None, lmax=None):
        """Return a magnetic field function H(r,theta,phi) for a given wavenumber index"""

        if lmax is None: lmax = self.lmax
        if index is None:
            index = np.s_[:]

        if not self.interior_computed: self.solve_interior()
        if not self.exterior_computed: self.solve_exterior()

        an = self.an[index, :lmax+1]
        bn = self.bn[index, :lmax+1]
        cn = self.cn[index, :lmax+1]
        dn = self.dn[index, :lmax+1]
        mat = self.material_data

        def H_func(r, theta, phi):
            H = np.zeros(shape = [3] + list(r.shape), dtype=np.complex)
            id_inside = r <= self.radius

            k = mat['k'][index]
            n = mat['n_b'][index]
            mu = mat['mu_b'][index]
            H[:,~id_inside] = scattered_H(an, bn, k, n, mu)(r[~id_inside], theta[~id_inside], phi[~id_inside])

            k = 2*np.pi*mat['n'][index]/self.wavelength[index]
            n = mat['n'][index]
            mu = mat['mu'][index]
            H[:,id_inside] = interior_H(cn, dn, k, n, mu)(r[id_inside], theta[id_inside], phi[id_inside])
            return H

        return H_func
