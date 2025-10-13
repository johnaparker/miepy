"""VSH sources."""

import numpy as np

import miepy
from miepy.sources import source


class vsh_source(source):
    def __init__(self, n, m, ftype="electric", center=None, mode=None, amplitude=1, phase=0):
        """Arguments:
        n            n value of VSH mode
        m            m value of VSH mode
        ftype        'electric' or 'magnetic' dipole (default: electric)
        center       center position of vsh_source
        mode         type of vsh_mode (default: incident)
        amplitude    amplitude of the source (default: 1)
        phase        additional phase factor (default: 0).
        """
        self.n = n
        self.m = m
        self.ftype = ftype

        self.mode = mode
        if mode is None:
            self.mode = miepy.vsh_mode.incident

        if self.ftype == "electric":
            self.N, self.M = miepy.vsh.VSH(self.n, self.m, self.mode)
            self.N_far, self.M_far = miepy.vsh.VSH_far(self.n, self.m, self.mode)
        elif self.ftype == "magnetic":
            self.M, self.N = miepy.vsh.VSH(self.n, self.m, self.mode)
            self.M_far, self.N_far = miepy.vsh.VSH_far(self.n, self.m, self.mode)
        else:
            raise ValueError("ftype must be either 'electric' or 'magnetic'")

        if center is None:
            self.center = np.array([0, 0, 0], dtype=float)
        else:
            self.center = np.asarray(center, dtype=float)

        self.amplitude = amplitude
        self.phase = phase

    def E_field(self, x1, x2, x3, k, far=False, spherical=False):
        factor = self.amplitude * np.exp(1j * self.phase)

        if not spherical:
            x1, x2, x3 = miepy.coordinates.cart_to_sph(x1, x2, x3, origin=self.center)

        if far:
            E = self.N_far(x1, x2, x3, k)
        else:
            E = self.N(x1, x2, x3, k)

        if not spherical:
            return factor * miepy.coordinates.vec_sph_to_cart(E, x2, x3)
        else:
            return factor * E

    def H_field(self, x1, x2, x3, k, far=False, spherical=False):
        factor = self.amplitude * np.exp(1j * self.phase)

        if not spherical:
            x1, x2, x3 = miepy.coordinates.cart_to_sph(x1, x2, x3, origin=self.center)

        if far:
            E = self.M_far(x1, x2, x3, k)
        else:
            E = self.M(x1, x2, x3, k)

        if not spherical:
            return factor * miepy.coordinates.vec_sph_to_cart(E, x2, x3)
        else:
            return factor * E

    def structure(self, position, k, lmax):
        position = np.asarray(position)
        Nparticles = len(position)

        rmax = miepy.vsh.lmax_to_rmax(lmax)
        p_src = np.zeros([Nparticles, 2, rmax], dtype=complex)
        factor = self.amplitude * np.exp(1j * self.phase)

        for i in range(Nparticles):
            dr = position[i] - self.center

            if not np.any(dr):
                for r, n, m in miepy.mode_indices(lmax):
                    if n == self.n and m == self.m:
                        Emn = miepy.vsh.Emn(m, n)
                        p_src[i, 0, r] = 1 / (-1j * Emn)
            else:
                rad, theta, phi = miepy.coordinates.cart_to_sph(*dr)
                for r, n, m in miepy.mode_indices(lmax):
                    Emn = miepy.vsh.Emn(self.m, self.n)
                    miepy.vsh.Emn(m, n)
                    A, B = miepy.cpp.vsh_translation.vsh_translation(
                        m, n, self.m, self.n, rad, theta, phi, k, self.mode
                    )
                    p_src[i, 0, r] = A / (-1j * Emn)
                    p_src[i, 1, r] = B / (-1j * Emn)

        if self.ftype == "magnetic":
            p_src = p_src[:, ::-1]

        return factor * p_src

    def H_angular(self, theta, phi, k, radius=None, origin=None):
        if radius is None:
            radius = 1e6 * 2 * np.pi / k

        return self.H_field(radius, radius, theta, phi, far=True, spherical=True)[1:]

    def E_angular(self, theta, phi, k, radius=None, origin=None):
        if radius is None:
            radius = 1e6 * 2 * np.pi / k

        return self.E_field(radius, radius, theta, phi, far=True, spherical=True)[1:]

    def angular_spectrum(self, theta, phi, k):
        return self.E_angular(theta, phi, k)
