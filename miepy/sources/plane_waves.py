"""
plane wave sources
"""

import numpy as np
import miepy
from miepy.vsh.special import pi_func, tau_func
from miepy.sources.source_base import source

class plane_wave(source):
    def __init__(self, polarization, theta=0, phi=0, amplitude=1, phase=0):
        """
        Create a plane-wave source. Default arguments provide a unit-amplitude, zero-propagating wave

        Arguments:
            polarization[2]      (TE, TM) values representing the polarization
            theta                theta spherical angle of k-vector
            phi                  phi spherical angle of k-vector
            amplitude            electric field amplitude E0
            phase                phase factor
        """
        super().__init__(amplitude, phase)
        self.polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization /= np.linalg.norm(self.polarization)

        self.theta = theta
        self.phi   = phi

        ### TE and TM vectors
        self.k_hat, self.n_te, self.n_tm = miepy.coordinates.sph_basis_vectors(theta, phi)

    @classmethod
    def from_string(cls, polarization, direction='z', amplitude=1, phase=0):
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

        return cls(pol, theta, phi, amplitude, phase)
    
    def E_field(self, x, y, z, k):
        amp = self.amplitude*np.exp(1j*k*(self.k_hat[0]*x + self.k_hat[1]*y + self.k_hat[2]*z))*np.exp(1j*self.phase)
        pol = self.n_te*self.polarization[0] + self.n_tm*self.polarization[1]
        return np.einsum('i...,...->i...', pol, amp)

    def H_field(self, x, y, z, k):
        amp = self.amplitude*np.exp(1j*k*(self.k_hat[0]*x + self.k_hat[1]*y + self.k_hat[2]*z))*np.exp(1j*self.phase)
        pol = self.n_tm*self.polarization[0] - self.n_te*self.polarization[1]
        return np.einsum('i...,...->i...', pol, amp)

    #TODO: use analytic expressions (commented) instead of using vsh_rotation 
    def structure(self, position, k, lmax, radius=None):
        rmax = miepy.vsh.lmax_to_rmax(lmax)
        p_src = np.zeros([2, rmax], dtype=complex)
        phase = k*(self.k_hat[0]*position[0] + self.k_hat[1]*position[1] + self.k_hat[2]*position[2]) + self.phase

        for i,n,m in miepy.mode_indices(lmax):
            pi_value = pi_func(n, m, 0)
            tau_value = tau_func(n, m, 0)

            p_src[0,i] = self.amplitude*np.exp(1j*phase)*np.sqrt(2*n+1)*tau_value*(self.polarization[0] - 1j*m*self.polarization[1])
            p_src[1,i] = self.amplitude*np.exp(1j*phase)*np.sqrt(2*n+1)*pi_value*(self.polarization[0] - 1j*m*self.polarization[1])

            if m == 1:
                p_src[:,i] /= (n*(n+1))

            # pi_value = pi_func(n, m)(self.theta)
            # tau_value = tau_func(n, m)(self.theta)
            # Emn = np.abs(miepy.vsh.Emn(m, n))

            # p_src[0,i] = -1*self.amplitude*np.exp(1j*phase)*Emn*(
                      # (tau_value*(self.polarization[0]*np.cos(self.phi) + self.polarization[1]*np.cos(self.phi - np.pi/2))) \
                  # + 1j*(pi_value*(self.polarization[0]*np.sin(self.phi) + self.polarization[1]*np.sin(self.phi - np.pi/2))))

            # p_src[1,i] = -1*self.amplitude*np.exp(1j*phase)*Emn*(
                      # (pi_value*(self.polarization[0]*np.cos(self.phi) + self.polarization[1]*np.cos(self.phi - np.pi/2))) \
                # + 1j*(tau_value*(self.polarization[0]*np.sin(self.phi) + self.polarization[1]*np.sin(self.phi - np.pi/2))))

        if self.theta != 0 or self.phi != 0:
            quat = miepy.quaternion.from_spherical_coords(self.theta, self.phi)
            p_src = miepy.vsh.rotate_expansion_coefficients(p_src, quat)

        return p_src

    def is_paraxial(self, k):
        return True
