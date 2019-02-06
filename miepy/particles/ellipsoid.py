import miepy
from .particle_base import particle

class ellipsoid(particle):
    def __init__(self, position, rx, ry, rz, material, orientation=None, tmatrix_lmax=0):
        """An ellipsoid object

        Arguments:
            position[3]   x,y,z position of particle
            rx,ry,rz      radii of the 3 axes
            material      particle material (miepy.material object)
            orientation   particle orientation
        """
        super().__init__(position, orientation, material)
        self.rx  = rx
        self.ry  = ry
        self.rz  = rz

        self.tmatrix_lmax = tmatrix_lmax

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    position = {self.position} m
    orientation = {self.orientation}
    rx = {self.rx:.2e} m
    ry = {self.rx:.2e} m
    rz = {self.rx:.2e} m
    material = {self.material}'''

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        calc_lmax = max(lmax+2, self.tmatrix_lmax)

        self.tmatrix_fixed = miepy.tmatrix.tmatrix_ellipsoid(self.rx, self.ry, self.rz, wavelength, 
                self.material.eps(wavelength), eps_m, calc_lmax, extended_precision=False,
                conducting=self.conducting)

        if lmax < calc_lmax:
            self.tmatrix_fixed = miepy.tmatrix.tmatrix_reduce_lmax(self.tmatrix_fixed, lmax)

        self._rotate_fixed_tmatrix()
        return self.tmatrix

    def enclosed_radius(self):
        return max(self.rx, self.ry, self.rz)

    def _dict_key(self, wavelength):
        return (ellipsoid, self.rx, self.ry, self.rz, self.material.eps(wavelength).item(), self.material.mu(wavelength).item())
