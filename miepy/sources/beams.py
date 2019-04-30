"""
Base classes for beams. Defines:

    beam________________base clase for beams
    polarized_beam______beam with a global polarization state (TE, TM pair)
    reflected_beam______beam reflected by an interface
    transmitted_beam____beam transmitted by an interface
"""

import numpy as np
from abc import ABCMeta, abstractmethod
import miepy
from miepy.sources import propagating_source, polarized_propagating_source
from functools import partial
from miepy.constants import Z0

class beam(propagating_source):
    """abstract base class for beam sources"""
    __metaclass__ = ABCMeta

    def __init__(self, power=1, theta_max=np.pi/2, phase=0, center=None, theta=0, phi=0, standing=False):
        propagating_source.__init__(self, amplitude=1, phase=phase, origin=center, theta=theta, phi=phi, standing=standing)
        self.power = power
        self.theta_max = theta_max
        self.center = self.origin

        self.k_stored = None
        self.p_src_func = None

    def E0(self, k):
        """
        Compute the amplitude constant of the beam
        
        Arguments:
            k    medium wavenumber
        """
        theta_c = self.theta_cutoff(k)

        theta = np.linspace(0, theta_c, 20)
        phi = np.linspace(0, 2*np.pi, 40)
        THETA, PHI = np.meshgrid(theta, phi)

        E = self.angular_spectrum(THETA, PHI, k)
        S = 0.5/Z0*np.linalg.norm(E, axis=0)**2*np.sin(THETA)
        P = miepy.vsh.misc.trapz_2d(theta, phi, S.T).real
        return np.sqrt(self.power/P)

    def E_field(self, x1, x2, x3, k, far=False, spherical=False, sampling=20, origin=None):
        """
        Compute the electric field

        Arguments:
            x1          x/r position (array-like) 
            x2          y/theta position (array-like) 
            x3          z/phi position (array-like) 
            k           medium wavenumber
            far         use expressions valid only for far-field (bool, default=False)
            spherical   input/output in spherical coordinates (bool, default=False)
            sampling    sampling to use in theta/phi integrals (default: 20)
            origin      origin for when spherical=True (default: self.center)
        """
        return self._field(self.angular_spectrum, self.E_angular, x1, x2, x3, k,
                far=far, spherical=spherical, sampling=sampling, origin=origin)

    def H_field(self, x1, x2, x3, k, far=False, spherical=False, sampling=20, origin=None):
        """
        Compute the magnetic field

        Arguments:
            x1          x/r position (array-like) 
            x2          y/theta position (array-like) 
            x3          z/phi position (array-like) 
            k           medium wavenumber
            far         use expressions valid only for far-field (bool, default=False)
            spherical   input/output in spherical coordinates (bool, default=False)
            sampling    sampling to use in theta/phi integrals (default: 20)
            origin      origin for when spherical=True (default: self.center)
        """
        return -1*self._field(self.H_angular_spectrum, self.H_angular, x1, x2, x3, k,
                    far=far, spherical=spherical, sampling=sampling, origin=origin)

    def E_angular(self, theta, phi, k, radius=None, origin=None):
        """
        Obtain the far-field electric field in spherical coordinates

        Arguments:
            theta    theta coordinates (array-like)
            phi      phi coordinates (array-like)
            k        medium wavenumber
            radius   radius coordinate (scalar or array-like)
            origin   origin around which to compute angular fields (default: self.center)
        """
        return self._angular(self.angular_spectrum, theta, phi, k, radius=radius, origin=origin)

    def H_angular(self, theta, phi, k, radius=None, origin=None):
        """
        Obtain the far-field magnetic field in spherical coordinates

        Arguments:
            theta    theta coordinates (array-like)
            phi      phi coordinates (array-like)
            k        medium wavenumber
            radius   radius coordinate (scalar or array-like)
            origin   origin around which to compute angular fields (default: self.center)
        """
        return -1*self._angular(self.H_angular_spectrum, theta, phi, k, radius=radius, origin=origin)

    def theta_cutoff(self, k, cutoff=1e-6, tol=1e-9):
        """
        Cutoff angle required for numerical intergration of the angular spectrum

        Arguments:
            k        medium wavenumber
            cutoff   fraction to determine where to cutoff theta_max
            tol      tolerance to obtain given cutoff
        """
        Nphi = 60
        theta = np.linspace(0, self.theta_max, Nphi)
        phi = np.linspace(0, 2*np.pi, Nphi)
        THETA, PHI = np.meshgrid(theta, phi)
        E = self.angular_spectrum(THETA, PHI, k)
        I = np.sum(np.abs(E)**2, axis=0)
        Imax = np.max(I)

        theta = self.theta_max
        dtheta = 0.01*theta

        err = np.abs(I/Imax - cutoff)

        while (err > tol).all():
            I = np.zeros_like(err)
            while (I/Imax < cutoff).all():
                theta -= dtheta
                E = self.angular_spectrum(theta, phi, k)
                I = np.sum(np.abs(E)**2, axis=0)

            err = np.abs(I/Imax - cutoff)
            theta += dtheta
            dtheta /= 2

        return theta

    def structure(self, position, k, lmax):
        """
        Obtain the structure coefficients of the beam

        Arguments:
            position     (x, y, z) position
            k            medium wavenumber
            lmax         maximum expansion order
        """
        theta_c = self.theta_cutoff(k)

        if k != self.k_stored:
            self.k_stored = k
            self.p_src_func = miepy.vsh.decomposition.integral_project_source_far(self, 
                               k, lmax, theta_0=np.pi - theta_c)

        pos = miepy.coordinates.rotate(*(position - self.center), self.orientation.inverse())
        p_src = self.p_src_func(pos)

        if self.orientation != miepy.quaternion.one:
            p_src = miepy.vsh.rotate_expansion_coefficients(p_src, self.orientation)

        return p_src

    def reflect(self, interface, medium, wavelength):
        return reflected_beam(self, interface, wavelength, medium)

    def transmit(self, interface, medium, wavelength):
        return transmitted_beam(self, interface, wavelength, medium)

    #TODO: Test dependence on center and orientation
    def _angular(self, angular_func, theta, phi, k, radius=None, origin=None):
        """
        Method to obtain E or H far field angular fields

        Arguments:
            angular_func    function for the angular spectrum
            theta    theta coordinates (array-like)
            phi      phi coordinates (array-like)
            k        medium wavenumber
            radius   radius coordinate (scalar or array-like)
            origin   origin around which to compute angular fields (default: self.center)
        """
        if radius is None:
            radius = 1e6*(2*np.pi/k)

        theta_r, phi_r = miepy.coordinates.rotate_sph(theta, phi, self.orientation.inverse())

        #TODO: can the lines below be reduced to some rotate_vec_sph function?
        E_inf = angular_func(theta_r, phi_r, k)
        E_inf = np.insert(E_inf, 0, 0, axis=0)
        E_inf = miepy.coordinates.vec_sph_to_cart(E_inf, theta_r, phi_r)
        E_inf = miepy.coordinates.rotate_vec(E_inf, self.orientation)
        E_inf = miepy.coordinates.vec_cart_to_sph(E_inf, theta, phi)[1:]

        sign = np.sign(np.pi/2 - theta_r)
        E0 = 1j*self.E0(k)*np.exp(1j*self.phase)*np.exp(sign*1j*k*radius)/radius

        factor = -(2*((theta_r < np.pi/2)) - 1)
        E_inf[1] *= factor

        if origin is not None:
            rhat, *_ = miepy.coordinates.sph_basis_vectors(theta, phi)
            phase = k*np.einsum('i...,i', rhat, origin - self.center)
            E_inf *= np.exp(1j*phase)

        return E0*E_inf


    def _field(self, angular_func, far_func, x1, x2, x3, k, far=False, spherical=False, sampling=20, origin=None):
        """
        Method to compute E or H field

        Arguments:
            angular_func    function for the angular spectrum
            far_func        function for the far-fields (used when far=True)
            x1          x/r position (array-like) 
            x2          y/theta position (array-like) 
            x3          z/phi position (array-like) 
            far         use expressions valid only for far-field (bool, default=False)
            spherical   input/output in spherical coordinates (bool, default=False)
            sampling    sampling to use in theta/phi integrals (default: 20)
            origin      origin for when spherical=True (default: self.center)
        """
        if origin is None:
            origin = self.center

        if far:
            if not spherical:
                x1, x2, x3 = miepy.coordinates.cart_to_sph(x1, x2, x3, origin=origin)
            E = far_func(x2, x3, k, radius=x1)

            if not spherical:
                E = np.insert(E, 0, 0, axis=0)
                E = miepy.coordinates.vec_sph_to_cart(E, x2, x3)
            return E

        if spherical:
            x1, x2, x3 = miepy.coordinates.sph_to_cart(x1, x2, x3, origin=origin)

        x1r, x2r, x3r = miepy.coordinates.translate(x1, x2, x3, -self.center)
        x1r, x2r, x3r = miepy.coordinates.rotate(x1r, x2r, x3r, self.orientation.inverse())
        rho, angle, z = miepy.coordinates.cart_to_cyl(x1r, x2r, x3r)

        theta_c = self.theta_cutoff(k)
        theta = np.linspace(np.pi - theta_c, np.pi, sampling)
        phi = np.linspace(0, 2*np.pi, 2*sampling)
        THETA, PHI = np.meshgrid(theta, phi, indexing='ij')

        E_inf = angular_func(THETA, PHI, k)
        E_inf = np.insert(E_inf, 0, 0, axis=0)
        E_inf = miepy.coordinates.vec_sph_to_cart(E_inf, THETA, PHI)

        @partial(np.vectorize, signature='(),(),()->(n)')
        def far_to_near(rho, angle, z):
            integrand = np.exp(1j*k*(-z*np.cos(THETA) - rho*np.sin(THETA)*np.cos(PHI - angle))) \
                        * E_inf*np.sin(THETA)
            return np.array([miepy.vsh.misc.trapz_2d(theta, phi, integrand[i]) for i in range(3)])

        E = far_to_near(rho, angle, z)
        E = np.moveaxis(E, source=-1, destination=0)
        E = miepy.coordinates.rotate_vec(E, self.orientation)

        A = k*self.E0(k)*np.exp(1j*self.phase)/(2*np.pi)
        E *= A

        if spherical:
            E = miepy.coordinates.vec_cart_to_sph(E, x2, x3)
        return E

class polarized_beam(beam, polarized_propagating_source):
    """abstract base class for polarized beam sources"""
    __metaclass__ = ABCMeta

    def __init__(self, polarization, power=1, theta_max=np.pi/2, phase=0, center=None, theta=0, phi=0, standing=False):
        polarized_propagating_source.__init__(self, polarization=polarization)
        beam.__init__(self, power=power, theta_max=theta_max, phase=phase, center=center,
            theta=theta, phi=phi, standing=standing)

class reflected_beam(beam):
    def __init__(self, incident_beam, interface, wavelength, medium):
        center = np.copy(incident_beam.center)
        center[2] += 2*(interface.z - center[2])
        theta = np.pi - incident_beam.theta
        phi = incident_beam.phi
        beam.__init__(self, power=incident_beam.power, theta_max=incident_beam.theta_max,
                phase=incident_beam.phase, center=center, theta=theta, phi=phi)

        self.incident_beam = incident_beam
        self.interface = interface
        self.wavelength = wavelength
        self.medium = medium

    def angular_spectrum(self, theta, phi, k):
        q = miepy.quaternion.from_spherical_coords(self.theta, self.phi)
        theta, phi = miepy.coordinates.rotate_sph(theta - np.pi, phi, q)
        theta = np.pi - theta

        r_parallel, r_perp = self.interface.reflection_coefficients(theta, self.wavelength, self.medium)
        q = miepy.quaternion.from_spherical_coords(self.incident_beam.theta, self.incident_beam.phi)
        theta, phi = miepy.coordinates.rotate_sph(theta - np.pi, phi, q.inverse())

        U = self.incident_beam.angular_spectrum(theta, phi, k)
        U[0] *= r_parallel
        U[1] *= r_perp

        return -U

    def E0(self, k):
        return self.incident_beam.E0(k)

    def theta_cutoff(self, k, cutoff=1e-6, tol=1e-9):
        return self.incident_beam.theta_cutoff(k, cutoff=cutoff, tol=tol)

class transmitted_beam(beam):
    def __init__(self, incident_beam, interface, wavelength, medium):
        theta = interface.transmission_angle(incident_beam.theta, wavelength, medium)
        phi = incident_beam.phi
        beam.__init__(self, power=incident_beam.power, theta_max=incident_beam.theta_max,
                phase=incident_beam.phase, center=incident_beam.center, theta=theta, phi=phi)

        self.incident_beam = incident_beam
        self.interface = interface
        self.wavelength = wavelength
        self.medium = medium

    def angular_spectrum(self, theta, phi, k):
        q = miepy.quaternion.from_spherical_coords(self.theta, self.phi)
        theta, phi = miepy.coordinates.rotate_sph(theta - np.pi, phi, q)
        m = self.interface.get_relative_index(self.wavelength, self.medium)
        theta = np.arcsin(m*np.sin(theta)).real

        t_parallel, t_perp = self.interface.transmission_coefficients(theta, self.wavelength, self.medium)
        q = miepy.quaternion.from_spherical_coords(self.incident_beam.theta, self.incident_beam.phi)
        theta, phi = miepy.coordinates.rotate_sph(theta - np.pi, phi, q.inverse())

        U = self.incident_beam.angular_spectrum(theta, phi, k)
        U[0] *= t_parallel
        U[1] *= t_perp

        return U

    def E0(self, k):
        return self.incident_beam.E0(k)

    def theta_cutoff(self, k, cutoff=1e-6, tol=1e-9):
        return self.incident_beam.theta_cutoff(k, cutoff=cutoff, tol=tol)
