import miepy
import numpy as np
from .particle_base import particle

class cylinder(particle):
    def __init__(self, position, radius, height, material, orientation=None, rounded=False, extended_precision=False, Nint=200, tmatrix_lmax=0):
        """A cylinder object

        Arguments:
            position[3]   x,y,z position of particle
            radius        cylinder radius
            height        cylinder height
            material      particle material (miepy.material object)
            orientation   particle orientation
        """
        super().__init__(position, orientation, material)
        self.radius = radius
        self.height = height
        self.extended_precision = extended_precision
        self.Nint = 200
        self.rounded = rounded

        self.tmatrix_lmax = tmatrix_lmax

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    position = {self.position} m
    orientation = {self.orientation}
    radius = {self.radius:.2e} m
    height = {self.height:.2e} m
    material = {self.material}
    rounded = {self.rounded}'''

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        calc_lmax = max(lmax+2, self.tmatrix_lmax)

        self.tmatrix_fixed = miepy.tmatrix.tmatrix_cylinder(self.radius, self.height, wavelength, 
                self.material.eps(wavelength), eps_m, calc_lmax, rounded=self.rounded, extended_precision=self.extended_precision, Nint=self.Nint,
                conducting=self.conducting)

        if lmax < calc_lmax:
            self.tmatrix_fixed = miepy.tmatrix.tmatrix_reduce_lmax(self.tmatrix_fixed, lmax)

        self._rotate_fixed_tmatrix()

        return self.tmatrix

    def enclosed_radius(self):
        return np.sqrt((self.height/2)**2 + self.radius**2)

    def _dict_key(self, wavelength):
        return (cylinder, self.radius, self.height, self.rounded, self.material.eps(wavelength).item(), self.material.mu(wavelength).item())
