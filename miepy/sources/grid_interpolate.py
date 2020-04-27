"""
Sources that use interpolation on a grid to evaluate structure coefficients
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from miepy.sources import source
import miepy
from miepy.cpp.decomposition import grid_interpolate
from functools import partial

class grid_interpolate_source(source):
    def __init__(self, source, grid):
        """
        Arguments:
            source     miepy.source object
            grid       [x, y, z] arrays representing the grid
        """
        self.source = source
        self.grid = grid
        self.k_stored = None
        self.lmax_stored = None

        if len(grid) != 3:
            raise ValueError('grid must specify x, y, and z values')
        if np.isscalar(grid[2]):
            self.ndim = 2
        else:
            self.ndim = 3

    def compute_functions(self):
        rmax = miepy.vsh.lmax_to_rmax(self.lmax_stored)

        if self.ndim == 3:
            data = np.empty((len(self.grid[0]), len(self.grid[1]), len(self.grid[2]), 2, rmax), dtype=complex)

            for i, x in enumerate(self.grid[0]):
                for j, y in enumerate(self.grid[1]):
                    for k, z in enumerate(self.grid[2]):
                        pos = (x, y, z)
                        data[i,j,k] = self.source.structure([pos], self.k_stored, self.lmax_stored)

            self.f_interp = RegularGridInterpolator(self.grid, data)

        elif self.ndim == 2:
            data = np.empty((len(self.grid[0]), len(self.grid[1]), 2, rmax), dtype=complex)

            for i, x in enumerate(self.grid[0]):
                for j, y in enumerate(self.grid[1]):
                    pos = (x, y, self.grid[2])
                    data[i,j] = self.source.structure([pos], self.k_stored, self.lmax_stored)

            # self.f_interp = RegularGridInterpolator(self.grid[:2], data)
            data=data.reshape([len(self.grid[0]), len(self.grid[1]), -1])
            def f_interp(pts):
                vals = grid_interpolate(grid=self.grid[:2], data=data, pts=pts)
                return vals.reshape([-1, 2, rmax])

            self.f_interp = f_interp

    def angular_spectrum(self, theta, phi, k):
        return self.source.angular_spectrum(theta, phi, k)

    def E_field(self, x1, x2, x3, k, far=False, spherical=False):
        return self.source.E_field(x1, x2, x3, k, far, spherical)

    def H_field(self, x1, x2, x3, k, far=False, spherical=False):
        return self.source.H_field(x1, x2, x3, k, far, spherical)

    def E_angular(self, theta, phi, k, radius=None, origin=None):
        return self.source.E_angular(theta, phi, k, radius, origin)

    def H_angular(self, theta, phi, k, radius=None, origin=None):
        return self.source.H_angular(theta, phi, k, radius, origin)

    def reflect(self, interface, medium, wavelength):
        return self.source.reflect(interface, medium, wavelength)

    def transmit(self, interface, medium, wavelength):
        return self.source.transmit(interface, medium, wavelength)

    def structure(self, position, k, lmax):
        if k != self.k_stored or lmax != self.lmax_stored:
            self.k_stored = k
            self.lmax_stored = lmax

            self.compute_functions()

        return self.f_interp(position[:, :self.ndim])
