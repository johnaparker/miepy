"""
plane wave sources
"""

import numpy as np
import miepy
from miepy.vsh.special import pi_func, tau_func
from miepy.sources import polarized_propagating_source

class plane_wave(polarized_propagating_source):
    def __init__(self, polarization, amplitude=1, phase=0, theta=0, phi=0, standing=False):
        """
        Create a plane-wave source. Default arguments provide a unit-amplitude, zero-propagating wave

        Arguments:
            polarization[2]      (TM, TE) values representing the polarization
            theta                theta spherical angle of k-vector
            phi                  phi spherical angle of k-vector
            amplitude            electric field amplitude E0
            phase                phase factor
        """
        polarized_propagating_source.__init__(self, polarization=polarization,
                amplitude=amplitude, phase=phase, origin=None, theta=theta, phi=phi, standing=standing)

    def __repr__(self):
        return f'plane_wave(polarization={self.polarization}, amplitude={self.amplitude}, theta={self.theta}, phi={self.phi})'

    @classmethod
    def from_string(cls, polarization, direction='z', amplitude=1, phase=0, standing=False):
        """Create a plane wave from string values for the polarization and direction
        
        Arguments:
            polarization     x, y, or z
            direction        x, y, z, -x, -y, or -z
            amplitude        electric field amplitude E0
            phase            phase factor
        """
        if direction in ['z', '+z']:
            theta = 0
            phi = 0
        elif direction == '-z':
            theta = np.pi
            phi = 0
        elif direction in ['x', '+x']:
            theta = np.pi/2
            phi = 0
        elif direction == '-x':
            theta = np.pi/2
            phi = np.pi
        elif direction in ['y', '+y']:
            theta = np.pi/2
            phi = np.pi/2
        elif direction == '-y':
            theta = np.pi/2
            phi = 3*np.pi/2
        else:
            raise ValueError("'{direction}' is not a valid direction of propagation. Use one of ['x', 'y', 'z', '-x', '-y', '-z']".format(direction=direction))

        if polarization == direction[-1]:
            raise ValueError('polarization cannot be the same as the direction of propagation')

        if polarization == 'x':
            if direction[-1] == 'z':
                pol = [1, 0]
            else:
                pol = [0, -1]
        elif polarization == 'y':
            if direction[-1] == 'x':
                pol = [0, 1]
            else:
                pol = [0, 1]
        elif polarization == 'z':
            pol = [-1, 0]
        elif polarization == 'rhc':
            pol = [1, 1j]
        elif polarization == 'lhc':
            pol = [1, -1j]
        else:
            raise ValueError("'{polarization}' is not a valid polarization. Use one of ['x', 'y', 'z', 'rhc', 'lhc']".format(polarization=polarization))

        return cls(polarization=pol, theta=theta, phi=phi, amplitude=amplitude, phase=phase, standing=standing)
    
    def E_field(self, x1, x2, x3, k, far=False, spherical=False):
        if spherical:
            x1, x2, x3 = miepy.coordinates.sph_to_cart(x1, x2, x3)

        amp = self.amplitude*np.exp(1j*k*(self.k_hat[0]*x1 + self.k_hat[1]*x2 + self.k_hat[2]*x3))*np.exp(1j*self.phase)
        pol = self.n_tm*self.polarization[0] + self.n_te*self.polarization[1]

        E = np.einsum('i...,...->i...', pol, amp)
        if spherical:
            E = miepy.coordinates.cart_to_sph(*E)

        return E

    def H_field(self, x1, x2, x3, k, far=False, spherical=False):
        if spherical:
            x1, x2, x3 = miepy.coordinates.sph_to_cart(x1, x2, x3)

        amp = self.amplitude*np.exp(1j*k*(self.k_hat[0]*x1 + self.k_hat[1]*x2 + self.k_hat[2]*x3))*np.exp(1j*self.phase)
        pol = self.n_te*self.polarization[0] - self.n_tm*self.polarization[1]

        H = np.einsum('i...,...->i...', pol, amp)
        if spherical:
            H = miepy.coordinates.cart_to_sph(*H)

        return H

    def structure(self, position, k, lmax):
        rmax = miepy.vsh.lmax_to_rmax(lmax)
        p_src = np.zeros([2, rmax], dtype=complex)
        phase = k*(self.k_hat[0]*position[0] + self.k_hat[1]*position[1] + self.k_hat[2]*position[2]) + self.phase

        for i,n,m in miepy.mode_indices(lmax):
            pi_value = pi_func(n, m, self.theta)
            tau_value = tau_func(n, m, self.theta)
            Emn = np.abs(miepy.vsh.Emn(m, n))
            factor = self.amplitude*np.exp(1j*(phase - m*self.phi))*Emn

            p_src[0,i] = factor*(tau_value*self.polarization[0] - 1j*pi_value*self.polarization[1])
            p_src[1,i] = factor*(pi_value*self.polarization[0]  - 1j*tau_value*self.polarization[1])

        return p_src

    def angular_spectrum(self, theta, phi, k):
        if theta == self.theta and phi == self.phi:
            return np.inf
        else:
            return 0

    def E_angular(self, theta, phi, k, radius=None):
        return self.angular_spectrum(theta, phi, j)

    def H_angular(self, theta, phi, k, radius=None):
        return self.angular_spectrum(theta, phi, j)

    def reflect(self, interface, medium, wavelength):
        theta = np.pi - self.theta
        phi = self.phi
        k = 2*np.pi*medium.index(wavelength)/wavelength
        phase = self.phase - 2*k*self.k_hat[2]*interface.z

        r_parallel, r_perp = interface.reflection_coefficients(self.theta, wavelength, medium)

        a_theta = r_parallel*self.polarization[0]
        a_phi = r_perp*self.polarization[1]
        polarization = [a_theta, a_phi]
        amplitude = np.linalg.norm(polarization)*self.amplitude

        if amplitude == 0:
            polarization = self.polarization

        return plane_wave(polarization=polarization, theta=theta, phi=phi, amplitude=amplitude, phase=phase)

    def transmit(self, interface, medium, wavelength):
        m = interface.get_relative_index(wavelength, medium)

        theta = np.arcsin(np.sin(self.theta)/m)
        phi = self.phi
        phase = self.phase

        t_parallel, t_perp = interface.transmission_coefficients(self.theta, wavelength, medium)

        a_theta = t_parallel*self.polarization[0]
        a_phi = t_perp*self.polarization[1]
        polarization = [a_theta, a_phi]
        amplitude = np.linalg.norm(polarization)*self.amplitude

        if amplitude == 0:
            polarization = self.polarization

        return plane_wave(polarization=polarization, theta=theta, phi=phi, amplitude=amplitude, phase=phase)
