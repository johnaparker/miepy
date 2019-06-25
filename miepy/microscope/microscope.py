import numpy as np
from tqdm import tqdm
import miepy
from functools import partial

def cluster_microscope(cluster, medium=None, orientation=None, focal_img=100, focal_obj=1, theta_obj=np.pi/2, sampling=30, source=False):
    """
    Arguments:
        cluster        miepy cluster
        medium         the outer medium of the microscope (default: air)
        orientation    orientation of the microscope (as a quaternion; default +z)
        focal_img      focal length of the imaging lens (default: 100)
        focal_obj      focal length of the objective lens (default: 100)
        theta_obj      maximum collection angle of the objective lens
        sampling       far-field sampling
        source         (bool) include the angular source fields (default: False)
    """
    if medium is None:
        medium = miepy.materials.air()

    E_angular = partial(cluster.E_angular, source=source)
    wavelength = cluster.wavelength
    n1 = cluster.medium.index(wavelength)
    n2 = medium.index(wavelength)

    return microscope(E_angular, wavelength, n1, n2=n2, orientation=orientation, focal_img=focal_img, focal_obj=focal_obj, theta_obj=theta_obj, sampling=sampling)


class microscope:
    """A microscope to produce images"""
    def __init__(self, E_angular, wavelength, n1, n2=1, orientation=None, focal_img=100, focal_obj=1, theta_obj=np.pi/2, sampling=30):
        """
        Arguments:
            E_angular      far-field angular function for E(theta, phi)
            wavelength     wavelength of light
            n1             refractive index at focal plane
            n2             refractive index at image plane (default: 1)
            orientation    orientation of the microscope (as a quaternion; default +z)
            focal_img      focal length of the imaging lens (default: 100)
            focal_obj      focal length of the objective lens (default: 100)
            theta_obj      maximum collection angle of the objective lens
            sampling       far-field sampling
            source         (bool) include the angular source fields (default: False)
        """
        self.focal_img = focal_img
        self.focal_obj = focal_obj
        self.theta_obj = theta_obj

        self.theta = np.linspace(0, self.theta_obj, sampling)
        self.phi   = np.linspace(0, 2*np.pi, 2*sampling)
        self.THETA, self.PHI = np.meshgrid(self.theta, self.phi, indexing='ij')

        self.n1 = n1
        self.n2 = n2
        self.k1 = 2*np.pi*self.n1/wavelength
        self.k2 = 2*np.pi*self.n2/wavelength

        if orientation is not None:
            THETA, PHI = miepy.coordinates.rotate_sph(self.THETA, self.PHI, orientation)
        else:
            THETA, PHI = self.THETA, self.PHI

        self.E_far = E_angular(THETA, PHI)
        self.E_far = np.insert(self.E_far, 0, 0, axis=0)
        self.E_far = miepy.coordinates.vec_sph_to_cart(self.E_far, THETA, PHI)
        if orientation is not None:
            self.E_far = miepy.coordinates.rotate_vec(self.E_far, orientation)

        self.magnification = self.n1*self.focal_img/(self.n2*self.focal_obj)
        self.numerical_aperature = self.n1*np.sin(theta_obj)

    def image(self, x, y, z_val=0, magnify=False):
        """
        Create an image

        Arguments:
            x_array    image x-values (array-like)
            y_array    image y-values (array-like)
            z_val      z-value of the image (relative to the focus)
            magnify    If True, magnify the input by the microscope's magnification
        """
        if magnify:
            M = self.magnification
            x, y = M*x, M*y
        k = self.k2
        f1 = self.focal_obj
        f2 = self.focal_img
        factor = 1j*k*f2*np.exp(-1j*k*(f2+z_val))/(2*np.pi)*np.sqrt(self.n1/self.n2)*(f1/f2)**2 \
                    *np.sin(self.THETA)*np.sqrt(np.cos(self.THETA)) \
                    *np.exp(0.5j*k*z_val*(f1/f2)**2*np.sin(self.THETA)**2)

        integrand = np.empty((2,) + self.THETA.shape, dtype=complex)

        @partial(np.vectorize, signature='(),()->(n)')
        def compute_pixel(x, y):
            rho = np.sqrt(x**2 + y**2)
            angle = np.arctan2(y, x)
            integrand[...] = factor*self.E_far[:2]*np.exp(1j*k*f1/f2*rho*np.sin(self.THETA)*np.cos(self.PHI-angle))
            val = np.zeros(2, dtype=complex)
            for p in range(2):
                val[p] = miepy.vsh.misc.trapz_2d(self.theta, self.phi, integrand[p])

            return val

        image = compute_pixel(x, y)
        image = np.moveaxis(image, -1, 0)
        
        return image
