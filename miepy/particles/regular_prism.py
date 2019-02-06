import miepy
from .particle_base import particle
import numpy as np

class regular_prism(particle):
    def __init__(self, position, N, height, material, width=None, radius=None, orientation=None, tmatrix_lmax=0):
        """A regular prism object (an extruded regular polygon). By convention, the initial orientation is such that
        the bottom edge is parallel to the x-axis.

        Arguments:
            position[3]   x,y,z position of particle
            N             number of vertices
            height        height of prism
            material      particle material (miepy.material object)
            width         width (side length) of prism (specifiy width or radius)
            radius        radius of prism, from center to vertex (specifiy radius or width)
            orientation   particle orientation
        """
        if width is None and radius is None:
            raise ValueError('A width or a radius must be specified')

        super().__init__(position, orientation, material)
        self.N = N
        self.height = height

        factor = np.sqrt(2*(1 - np.cos(2*np.pi/N)))
        if width is None:
            self.radius = radius
            self.width = radius*factor
        else:
            self.radius = width/factor
            self.width = width

        self.tmatrix_lmax = tmatrix_lmax

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    position = {self.position} m
    orientation = {self.orientation}
    vertices = {self.N}
    width = {self.width:.2e} m
    height = {self.height:.2e} m
    material = {self.material}'''

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        calc_lmax = max(lmax+2, self.tmatrix_lmax)

        self.tmatrix_fixed = miepy.tmatrix.tmatrix_regular_prism(self.N, self.width, self.height, wavelength, 
                self.material.eps(wavelength), eps_m, calc_lmax, extended_precision=False, conducting=self.conducting)

        if lmax < calc_lmax:
            self.tmatrix_fixed = miepy.tmatrix.tmatrix_reduce_lmax(self.tmatrix_fixed, lmax)

        self._rotate_fixed_tmatrix()
        return self.tmatrix

    def enclosed_radius(self):
        return np.sqrt(self.radius**2 + (self.height/2)**2)

    def _dict_key(self, wavelength):
        return (regular_prism, self.N, self.width, self.height, self.material.eps(wavelength).item(), self.material.mu(wavelength).item())


def cube(position, width, material, orientation=None, tmatrix_lmax=0):
    """A cube is a type of regular prism

    Arguments:
        position[3]   x,y,z position of particle
        width         width (side length) of cube
        material      particle material (miepy.material object)
        orientation   particle orientation
    """
    return regular_prism(position, N=4, width=width, height=width, material=material,
            orientation=orientation, tmatrix_lmax=tmatrix_lmax)
