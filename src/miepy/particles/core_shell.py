import miepy
from .particle_base import particle

class core_shell(particle):
    def __init__(self, position, core_radius, shell_thickness, core_material, shell_material):
        """A sphere object

        Arguments:
            position[3]   x,y,z position of particle
            radius        sphere radius
            material      particle material (miepy.material object)
        """
        super().__init__(position, None, core_material)
        self.core_radius = core_radius
        self.shell_thickness = shell_thickness
        self.core_material = core_material
        self.shell_material = shell_material

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    position = {self.position} m
    radius = {self.radius:.2e} m
    material = {self.material}'''

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        self.tmatrix = miepy.tmatrix.tmatrix_core_shell(self.core_radius, self.shell_thickness, wavelength, 
                self.core_material.eps(wavelength), self.shell_material.eps(wavelength), eps_m, lmax)
        self.tmatrix_fixed = self.tmatrix

        return self.tmatrix

    def enclosed_radius(self):
        return self.core_radius + self.shell_thickness

    def _dict_key(self, wavelength):
        return (core_shell, self.core_radius, self.shell_thickness, 
                self.core_material.eps(wavelength).item(), self.core_material.mu(wavelength).item(),
                self.shell_material.eps(wavelength).item(), self.shell_material.mu(wavelength).item())
