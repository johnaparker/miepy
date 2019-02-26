"""
Abstract base classes for beams. Defines:

    beam______________base clase for beams
    polarized_beam____beam with a global polarization state (TE, TM pair)
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
        radius = 1e6*2*np.pi/k
        theta_c = self.theta_cutoff(k)

        theta = np.linspace(theta_c, np.pi, 30)
        phi = np.linspace(0, 2*np.pi, 20)
        THETA, PHI = np.meshgrid(theta, phi)
        R = radius*np.ones_like(THETA)

        E = self.angular_spectrum(THETA, PHI, k)
        S = 0.5/Z0*np.abs(E)**2*np.sin(THETA)
        P = radius**2*miepy.vsh.misc.trapz_2d(theta, phi, S.T).real/4

        return np.sqrt(P)

    def E_field(self, x1, x2, x3, k, far=False, spherical=False, sampling=20):
        x1r, x2r, x3r = miepy.coordinates.translate(x1, x2, x3, -self.center)
        x1r, x2r, x3r = miepy.coordinates.rotate(x1r, x2r, x3r, self.orientation.inverse())
        rho, angle, z = miepy.coordinates.cart_to_cyl(x1r, x2r, x3r)

        theta_c = self.theta_cutoff(k)
        theta = np.linspace(np.pi - theta_c, np.pi, sampling)
        phi = np.linspace(0, 2*np.pi, 2*sampling)
        THETA, PHI = np.meshgrid(theta, phi)

        E_inf = np.zeros((3,) + THETA.shape, dtype=complex)
        E_inf[1:] = self.angular_spectrum(THETA, PHI, k)
        E_inf = miepy.coordinates.vec_sph_to_cart(E_inf, THETA, PHI)

        @partial(np.vectorize, signature='(),(),()->(n)')
        def far_to_near(rho, angle, z):
            integrand = np.exp(1j*k*(z*np.cos(THETA) + rho*np.sin(THETA)*np.cos(PHI - angle))) \
                        * E_inf*np.sin(THETA)
            return np.array([miepy.vsh.misc.trapz_2d(theta, phi, integrand[i].T) for i in range(3)])

        E = far_to_near(rho, angle, z)
        E = np.moveaxis(E, source=-1, destination=0)
        E = miepy.coordinates.rotate_vec(E, self.orientation)

        E0 = self.E0(k)*np.exp(1j*self.phase)
        return E0*E

    def H_field(self, x, y, z, k):
        x1r, x2r, x3r = miepy.coordinates.translate(x1, x2, x3, -self.center)
        x1r, x2r, x3r = miepy.coordinates.rotate(x1r, x2r, x3r, self.orientation.inverse())
        rho, angle, z = miepy.coordinates.cart_to_cyl(x1r, x2r, x3r)

        theta_c = self.theta_cutoff(k)
        theta = np.linspace(0, theta_c, sampling)
        phi = np.linspace(0, 2*np.pi, 2*sampling)
        THETA, PHI = np.meshgrid(theta, phi)

        H_inf = np.zeros((3,) + THETA.shape, dtype=complex)
        H_inf[1:] = self.angular_spectrum(THETA, PHI, k)[::-1]
        H_inf = miepy.coordinates.vec_sph_to_cart(H_inf, THETA, PHI)

        @partial(np.vectorize, signature='(),(),()->(n)')
        def far_to_near(rho, angle, z):
            integrand = np.exp(1j*k*(z*np.cos(THETA) + rho*np.sin(THETA)*np.cos(PHI - angle))) \
                        * H_inf*np.sin(THETA)
            return np.array([miepy.vsh.misc.trapz_2d(theta, phi, integrand[i].T) for i in range(3)])

        H = far_to_near(rho, angle, z)
        H = np.moveaxis(J, source=-1, destination=0)
        H = miepy.coordinates.rotate(*H, self.orientation)

        E0 = self.E0(k)*np.exp(1j*self.phase)
        return E0*H

    def E_angular(self, theta, phi, k, radius=None):
        if radius is None:
            radius = 1e6*(2*np.pi/k)

        E0 = self.E0(k)*np.exp(1j*self.phase)*np.exp(1j*k*radius)/(k*radius)
        E = E0*self.angular_spectrum(theta - self.theta, phi - self.phi, k)
        return E

    def H_angular(self, theta, phi, k, radius=None):
        if radius is None:
            radius = 1e6*(2*np.pi/k)

        E0 = self.E0(k)*np.exp(1j*self.phase)*np.exp(1j*k*radius)/(k*radius)
        E = E0*self.angular_spectrum(theta - self.theta, phi - self.phi, k)
        H = E[::-1]
        return H

    def theta_cutoff(self, k, eps=1e-3):
        return min(self.theta_max, np.pi/2)

    def structure(self, position, k, lmax):
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


class polarized_beam(beam, polarized_propagating_source):
    """abstract base class for polarized beam sources"""
    __metaclass__ = ABCMeta

    def __init__(self, polarization, power=1, theta_max=np.pi/2, phase=0, center=None, theta=0, phi=0, standing=False):
        polarized_propagating_source.__init__(self, polarization=polarization)
        beam.__init__(self, power=power, theta_max=theta_max, phase=phase, center=center,
            theta=theta, phi=phi, standing=standing)
