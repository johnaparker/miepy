import miepy
from .particle_base import particle

class sphere(particle):
    def __init__(self, position, radius, material):
        """A sphere object

        Arguments:
            position[3]   x,y,z position of particle
            radius        sphere radius
            material      particle material (miepy.material object)
        """
        super().__init__(position, None, material)
        self.radius = radius

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    position = {self.position} m
    radius = {self.radius:.2e} m
    material = {self.material}'''

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        print(self.conducting)
        self.tmatrix = miepy.tmatrix.tmatrix_sphere(self.radius, wavelength, 
                self.material.eps(wavelength), eps_m, lmax, conducting=self.conducting)
        self.tmatrix_fixed = self.tmatrix

        return self.tmatrix

    def enclosed_radius(self):
        return self.radius

    def _dict_key(self, wavelength):
        return (sphere, self.radius, self.material.eps(wavelength).item(), self.material.mu(wavelength).item())
