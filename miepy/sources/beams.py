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
        theta_c = self.theta_cutoff(k)

        theta = np.linspace(0, theta_c, 20)
        phi = np.linspace(0, 2*np.pi, 40)
        THETA, PHI = np.meshgrid(theta, phi)

        E = self.angular_spectrum(THETA, PHI, k)
        S = 0.5/Z0*np.linalg.norm(E, axis=0)**2*np.sin(THETA)
        P = miepy.vsh.misc.trapz_2d(theta, phi, S.T).real
        return np.sqrt(self.power/P)

    def E_field(self, x1, x2, x3, k, far=False, spherical=False, sampling=20):
        x1r, x2r, x3r = miepy.coordinates.translate(x1, x2, x3, -self.center)
        x1r, x2r, x3r = miepy.coordinates.rotate(x1r, x2r, x3r, self.orientation.inverse())
        rho, angle, z = miepy.coordinates.cart_to_cyl(x1r, x2r, x3r)

        theta_c = self.theta_cutoff(k)
        theta = np.linspace(np.pi - theta_c, np.pi, sampling)
        phi = np.linspace(0, 2*np.pi, 2*sampling)
        THETA, PHI = np.meshgrid(theta, phi)

        E_inf = self.angular_spectrum(THETA, PHI, k)
        E_inf = np.insert(E_inf, 0, 0, axis=0)
        E_inf = miepy.coordinates.vec_sph_to_cart(E_inf, THETA, PHI)

        @partial(np.vectorize, signature='(),(),()->(n)')
        def far_to_near(rho, angle, z):
            integrand = np.exp(-1j*k*(z*np.cos(THETA) + rho*np.sin(THETA)*np.cos(PHI - angle))) \
                        * E_inf*np.sin(THETA)
            return np.array([miepy.vsh.misc.trapz_2d(theta, phi, integrand[i].T) for i in range(3)])

        E = far_to_near(rho, angle, z)
        E = np.moveaxis(E, source=-1, destination=0)
        E = miepy.coordinates.rotate_vec(E, self.orientation)

        A = k*self.E0(k)*np.exp(1j*self.phase)/(2*np.pi)
        return A*E

    def H_field(self, x1, x2, x3, k, far=False, spherical=False, sampling=20):
        x1r, x2r, x3r = miepy.coordinates.translate(x1, x2, x3, -self.center)
        x1r, x2r, x3r = miepy.coordinates.rotate(x1r, x2r, x3r, self.orientation.inverse())
        rho, angle, z = miepy.coordinates.cart_to_cyl(x1r, x2r, x3r)

        theta_c = self.theta_cutoff(k)
        theta = np.linspace(np.pi - theta_c, np.pi, sampling)
        phi = np.linspace(0, 2*np.pi, 2*sampling)
        THETA, PHI = np.meshgrid(theta, phi)

        H_inf = self.angular_spectrum(THETA, PHI, k)[::-1]
        H_inf[0] *= -1
        H_inf = np.insert(H_inf, 0, 0, axis=0)
        H_inf = miepy.coordinates.vec_sph_to_cart(H_inf, THETA, PHI)

        @partial(np.vectorize, signature='(),(),()->(n)')
        def far_to_near(rho, angle, z):
            integrand = np.exp(-1j*k*(z*np.cos(THETA) + rho*np.sin(THETA)*np.cos(PHI - angle))) \
                        * H_inf*np.sin(THETA)
            return np.array([miepy.vsh.misc.trapz_2d(theta, phi, integrand[i].T) for i in range(3)])

        H = far_to_near(rho, angle, z)
        H = np.moveaxis(H, source=-1, destination=0)
        H = miepy.coordinates.rotate_vec(H, self.orientation)

        A = k*self.E0(k)*np.exp(1j*self.phase)/(2*np.pi)
        return -A*H

    def E_angular(self, theta, phi, k, radius=None):
        if radius is None:
            radius = 1e6*(2*np.pi/k)

        x, y, z = miepy.coordinates.sph_to_cart(radius, theta, phi)
        xr, yr, zr = miepy.coordinates.rotate(x, y, z, self.orientation.inverse())
        _, theta_r, phi_r = miepy.coordinates.cart_to_sph(xr, yr, zr)

        E_inf = self.angular_spectrum(np.pi - theta_r, phi_r, k)
        E_inf = np.insert(E_inf, 0, 0, axis=0)
        E_inf = miepy.coordinates.vec_sph_to_cart(E_inf, theta_r, phi_r)
        E_inf = miepy.coordinates.rotate_vec(E_inf, self.orientation)
        E_inf = miepy.coordinates.vec_cart_to_sph(E_inf, theta, phi)

        sign = np.sign(zr)
        E0 = self.E0(k)*np.exp(1j*self.phase)*np.exp(sign*1j*k*radius)/radius

        return E0*E_inf[1:]

    def H_angular(self, theta, phi, k, radius=None):
        if radius is None:
            radius = 1e6*(2*np.pi/k)

        x, y, z = miepy.coordinates.sph_to_cart(radius, theta, phi)
        xr, yr, zr = miepy.coordinates.rotate(x, y, z, self.orientation.inverse())
        _, theta_r, phi_r = miepy.coordinates.cart_to_sph(xr, yr, zr)

        H_inf = self.angular_spectrum(np.pi - theta_r, phi_r, k)[::-1]
        H_inf[0] *= -1
        H_inf = np.insert(H_inf, 0, 0, axis=0)
        H_inf = miepy.coordinates.vec_sph_to_cart(H_inf, theta_r, phi_r)
        H_inf = miepy.coordinates.rotate_vec(H_inf, self.orientation)
        H_inf = miepy.coordinates.vec_cart_to_sph(H_inf, theta, phi)

        sign = np.sign(zr)
        E0 = self.E0(k)*np.exp(1j*self.phase)*np.exp(sign*1j*k*radius)/radius

        return sign*E0*H_inf[1:]

    def theta_cutoff(self, k, cutoff=1e-6, tol=1e-9):
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
