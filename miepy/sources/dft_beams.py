import numpy as np
import miepy
from miepy.sources import beam
from scipy.interpolate import RectBivariateSpline

def scalar_dft_beam(Ufunc, polarization, xmax, sampling=60, power=1, theta_max=np.pi/2, phase=0, center=None,
                    theta=0, phi=0, standing=False):
    """
    Scalar-function version of the DFT beam

    Arguments:
        Ufunc    (x,y)->scalar function for the complex scalar field
        polarization[2]   (TM, TE) values representing the polarization
        xmax     maximum x-value to integrate Ufunc over
    """
    polarization = np.asarray(polarization, dtype=complex)
    polarization /= np.linalg.norm(polarization)

    def Efunc(x, y):
        U = Ufunc(x,y)
        return np.array([polarization[0]*U, polarization[1]*U])

    return dft_beam(Efunc, xmax, sampling=sampling, power=power, theta_max=theta_max, phase=phase,
             center=center, theta=theta, phi=phi, standing=standing)

class dft_beam(beam):
    """
    Use DFT on fields in an xy-plane to produce a beam
    """
    def __init__(self, Efunc, xmax, sampling=60, power=1, theta_max=np.pi/2, phase=0, center=None,
                theta=0, phi=0, standing=False):
        """
        Arguments:
            Efunc    (x,y)->[2] function for the complex (Ex,Ey) field
            polarization[2]   (TM, TE) values representing the polarization
            xmax     maximum x-value to integrate Ufunc over
        """
        super().__init__(power=power, theta_max=theta_max,
                phase=phase, center=center, theta=theta, phi=phi, standing=standing)

        self.Efunc = Efunc
        self.xmax = xmax
        self.sampling = sampling
        self.Ex_func = None
        self.Ey_func = None

    def compute_spectrum(self, k):
        theta = np.linspace(0, self.theta_max, self.sampling)
        phi = np.linspace(0, 2*np.pi, self.sampling)
        THETA, PHI = np.meshgrid(theta, phi, indexing='ij')

        kz = k*np.cos(THETA)
        kx = k*np.sin(THETA)*np.cos(PHI)
        ky = k*np.sin(THETA)*np.sin(PHI)
        
        amp = np.zeros((2,) + kx.shape, dtype=complex)
        x = np.linspace(-self.xmax, self.xmax, self.sampling)
        y = np.linspace(-self.xmax, self.xmax, self.sampling)
        X, Y = np.meshgrid(x, y, indexing='ij')
        E = self.Efunc(X, Y)
        for i in range(kx.shape[0]):
            for j in range(kx.shape[1]):
                exp_term = np.exp(-1j*(kx[i,j]*X + ky[i,j]*Y))
                amp[0,i,j] = miepy.vsh.misc.trapz_2d(x, y, E[0]*exp_term)
                amp[1,i,j] = miepy.vsh.misc.trapz_2d(x, y, E[1]*exp_term)

        Exr = RectBivariateSpline(theta, phi, amp[0].real)
        Exi = RectBivariateSpline(theta, phi, amp[0].imag)
        Eyr = RectBivariateSpline(theta, phi, amp[1].real)
        Eyi = RectBivariateSpline(theta, phi, amp[1].imag)

        self.Ex_func = lambda theta, phi: Exr.ev(theta, phi) + 1j*Exi.ev(theta, phi)
        self.Ey_func = lambda theta, phi: Eyr.ev(theta, phi) + 1j*Eyi.ev(theta, phi)

    def __repr__(self):
        return f'dft_beam(width={self.width}, polarization={self.polarization}, power={self.power}, ' \
               f'center={self.center}, theta={self.theta}, phi={self.phi})'

    def angular_spectrum(self, theta, phi, k):
        if self.Ex_func is None:
            self.compute_spectrum(k)

        theta = np.minimum(theta, np.pi - theta)  # functions are valid for theta < pi/2

        Ex = self.Ex_func(theta, phi)
        Ey = self.Ey_func(theta, phi)

        idx = theta > self.theta_max  #if theta > theta_max, angular_spectrum should vanish
        Ex[idx] = 0
        Ey[idx] = 0

        Esph = np.array([-Ex*np.cos(phi) - Ey*np.sin(phi),
                         -Ex*np.sin(phi) + Ey*np.cos(phi)], dtype=complex)

        return Esph
