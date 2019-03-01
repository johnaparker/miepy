import numpy as np
from tqdm import tqdm
import miepy
from functools import partial

class microscope:
    """A microscope to produce images of a cluster"""
    def __init__(self, cluster, medium=None, focal_img=100, focal_obj=1, theta_obj=np.pi/2, sampling=30, source=False):
        """
        Arguments:
            cluster        miepy cluster
            medium         the outer medium of the microscope (default: air)
            focal_img      focal length of the imaging lens (default: 100)
            focal_obj      focal length of the objective lens (default: 100)
            theta_obj      maximum collection angle of the objective lens
            sampling       far-field sampling
            source         (bool) include the angular source fields (default: False)
        """
        self.cluster = cluster
        self.focal_img = focal_img
        self.focal_obj = focal_obj
        self.theta_obj = theta_obj

        if medium is None:
            self.medium = miepy.materials.air()
        else:
            self.medium = medium

        self.theta = np.linspace(0, self.theta_obj, sampling)
        self.phi   = np.linspace(0, 2*np.pi, 2*sampling)
        self.THETA, self.PHI = np.meshgrid(self.theta, self.phi, indexing='ij')

        self.n1 = self.cluster.medium.eps(cluster.wavelength)**0.5
        self.n2 = self.medium.eps(cluster.wavelength)**0.5
        self.k1 = 2*np.pi*self.n1/cluster.wavelength
        self.k2 = 2*np.pi*self.n2/cluster.wavelength

        self.E_far = cluster.E_angular(self.THETA, self.PHI, source=source)
        self.E_far = np.insert(self.E_far, 0, 0, axis=0)
        self.E_far = miepy.coordinates.vec_sph_to_cart(self.E_far, self.THETA, self.PHI)

        self.magnification = self.n1*self.focal_img/(self.n2*self.focal_obj)
        self.numerical_aperature = self.n1*np.sin(theta_obj)

    def image(self, x, y, z_val=0):
        """
        Create an image of the cluster

        Arguments:
            x_array    1D array of image x-values
            y_array    1D array of image y-values
            z_val      z-value of the image (relative to the focus)
        """
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
