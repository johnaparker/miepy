import numpy as np
from tqdm import tqdm
import miepy

class microscope:
    """A microscope to produce images of a cluster"""
    def __init__(self, cluster, focal_length, theta_obj=np.pi/2, sampling=30):
        """
        Arguments:
            cluster        miepy cluster
            focal_length   focal length of the microscope
            theta_obj      maximum collection angle of the collection objective
            sampling       far-field sampling
        """
        self.cluster = cluster
        self.focal_length = focal_length
        self.theta_obj = theta_obj
        self.k = cluster.material_data.k_b

        self.theta = np.linspace(0, self.theta_obj, sampling)
        self.phi   = np.linspace(0, 2*np.pi, 2*sampling)
        self.THETA, self.PHI = np.meshgrid(self.theta, self.phi, indexing='ij')

        n1 = n2 = self.cluster.medium.eps(2*np.pi/self.k)**0.5

        self.E_far = cluster.E_angular(self.THETA, self.PHI)
        self.E_far *= np.sqrt(n1/n2)*np.sqrt(np.cos(self.THETA))
        # self.E_far[1] *= np.cos(self.THETA)**2

    def image(self, x_array, y_array, z_val=0):
        """
        Create an image of the cluster

        Arguments:
            x_array    1D array of image x-values
            y_array    1D array of image y-values
            z_val      z-value of the image (relative to the focus)
        """
        k = self.k
        f = self.focal_length
        factor = 1j*k*f*np.exp(-1j*k*f)/(2*np.pi)

        image = np.zeros((3, len(x_array), len(y_array)), dtype=complex)

        for i,x in enumerate(tqdm(x_array)):
            for j,y in enumerate(y_array):
                rho = np.sqrt(x**2 + y**2)
                angle = np.arctan2(y, x)
                angle = np.arctan2(x, y)
                integrand = factor*self.E_far*np.exp(1j*k*z_val*np.cos(self.THETA))*np.exp(1j*k*rho*np.sin(self.THETA)*np.cos(self.PHI-angle))*np.sin(self.THETA)
                integrand = miepy.coordinates.vec_sph_to_cart(integrand, self.THETA, self.PHI)

                for p in range(3):
                    image[p,i,j] = miepy.vsh.misc.trapz_2d(self.theta, self.phi, integrand[p])
        
        return image
