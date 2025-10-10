import miepy
from .particle_base import particle

class spheroid(particle):
    def __init__(self, position, axis_xy, axis_z, material, orientation=None, tmatrix_lmax=0):
        """A spheroid object

        Arguments:
            position[3]   x,y,z position of particle
            axis_xy       length of semiaxes perpendicular to the axis of symmetry
            axis_z        length of semiaxis along axis of symmetry
            material      particle material (miepy.material object)
            orientation   particle orientation
        """
        super().__init__(position, orientation, material)
        self.axis_xy = axis_xy
        self.axis_z  = axis_z

        self.tmatrix_lmax = tmatrix_lmax

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    position = {self.position} m
    orientation = {self.orientation}
    axis_xy = {self.axis_xy:.2e} m
    axis_z = {self.axis_z:.2e} m
    material = {self.material}'''

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        calc_lmax = max(lmax+2, self.tmatrix_lmax)

        self.tmatrix_fixed = miepy.tmatrix.tmatrix_spheroid(self.axis_xy, self.axis_z, wavelength, 
                self.material.eps(wavelength), eps_m, calc_lmax, extended_precision=False, conducting=self.conducting)

        if lmax < calc_lmax:
            self.tmatrix_fixed = miepy.tmatrix.tmatrix_reduce_lmax(self.tmatrix_fixed, lmax)

        self._rotate_fixed_tmatrix()
        return self.tmatrix

    def enclosed_radius(self):
        return max(self.axis_xy, self.axis_z)

    def _dict_key(self, wavelength):
        return (spheroid, self.axis_xy, self.axis_z, self.material.eps(wavelength).item(), self.material.mu(wavelength).item())
