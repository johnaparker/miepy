import miepy
from .particle_base import particle

class spheroid(particle):
    def __init__(position, axis_xy, axis_z, material, orientation=None):
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

    def is_inside(pos):
        pass

    def compute_tmatrix(lmax, wavelength, eps_m, **kwargs):
        self.tmatrix = miepy.tmatrix.tmatrix_spheroid(self.axis_xy, self.axis_z, wavelength, 
                self.material.eps(wavelength), eps_m, lmax, extended_precision=False)
        return self.tmatrix
