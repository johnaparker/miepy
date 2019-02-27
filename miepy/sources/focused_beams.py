"""
Focused vector beam sources that represent exact solutions to Maxwell's equations.
Includes lowest order Hermite-Gaussian and Laguerre-Gaussian modes.
For paraxial approximate versions of these beams, see sources/paraxial_beams.
"""
import numpy as np
import scipy.special as special
import scipy.integrate as integrate
from scipy.optimize import curve_fit
import miepy
from miepy.sources import source, combined_source
j0 = special.j0
j1 = special.j1
jv = special.jv

def filling_factor(f, w0, theta, theta_max):
    f0 = w0/(f*np.sin(theta_max))
    arg = -1/f0**2 * (np.sin(theta)/np.sin(theta_max))**2
    return np.exp(arg)

def I_generic(func, k, theta_max, f, w0):
    @np.vectorize
    def compute_integral(rho, z):
        def integral(theta):
            integrand = filling_factor(f, w0, theta, theta_max)*np.cos(theta)**.5*func(theta,rho)*np.e**(1j*k*z*np.cos(theta))
            return integrand

        return integrate.quad(lambda theta: integral(theta).real, 0, theta_max)[0]  \
                      + 1j*integrate.quad(lambda theta: integral(theta).imag, 0, theta_max)[0]

    return compute_integral

def I00(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)*(1+np.cos(theta))*j0(k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def I01(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)**2*j1(k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def I02(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)*(1-np.cos(theta))*jv(2,k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def I10(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)**3*j0(k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def I11(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)**2*(1+3*np.cos(theta))*j1(k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def I12(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)**2*(1-np.cos(theta))*j1(k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def I13(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)**3*jv(2,k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def I14(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)**2*(1-np.cos(theta))*jv(3,k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def Irad(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.cos(theta)*np.sin(theta)**2*j1(k*rho*np.sin(theta)), 
            k, theta_max, f, w0)

def Iazi(k, theta_max, f, w0):
    return I_generic(
            lambda theta,rho: np.sin(theta)**2*j1(k*rho*np.sin(theta)), 
            k, theta_max, f, w0)


class HG_00(source):
    def __init__(self, width, theta_max, focal_length):
        self.width = width
        self.theta_max = theta_max
        self.focal_length = focal_length

    def E_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g00 = I00(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g01 = I01(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g02 = I02(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)
        Ex = g00 + g02*np.cos(2*phi)
        Ey = g02*np.sin(2*phi)
        Ez = -2j*g01*np.cos(phi)
        E = np.array([Ex,Ey,Ez])*factor

        return E

    def H_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g00 = I00(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g01 = I01(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g02 = I02(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)

        Hx = g02*np.sin(2*phi)
        Hy = g00 - g02*np.cos(2*phi)
        Hz = -2j*g01*np.sin(phi)
        H = np.array([Hx,Hy,Hz])*factor
        
        return H

class HG_10(source):
    def __init__(self, width, theta_max, focal_length):
        self.width = width
        self.theta_max = theta_max
        self.focal_length = focal_length

    def E_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g10 = I10(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g11 = I11(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g12 = I12(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g13 = I13(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g14 = I14(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)
    
        Ex = 1j*g11*np.cos(phi) + 1j*g14*np.cos(3*phi)
        Ey = -1j*g12*np.sin(phi) + 1j*g14*np.sin(3*phi)
        Ez = -2*g10 + 2*g13*np.cos(2*phi)
        E = np.array([Ex,Ey,Ez])*factor

        return E

    def H_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g10 = I10(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g11 = I11(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g12 = I12(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g13 = I13(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g14 = I14(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)
    
        Hx = -1j*g12*np.sin(phi) + 1j*g14*np.sin(3*phi)
        Hy = 1j*(g11+2*g12)*np.cos(phi) - 1j*g14*np.cos(3*phi)
        Hz = 2*g13*np.sin(2*phi)
        H = np.array([Hx,Hy,Hz])*factor
    
        return H

class HG_01(source):
    def __init__(self, width, theta_max, focal_length):
        self.width = width
        self.theta_max = theta_max
        self.focal_length = focal_length

    def E_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g10 = I10(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g11 = I11(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g12 = I12(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g13 = I13(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g14 = I14(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)
    
        Ex = 1j*(g11+2*g12)*np.sin(phi) + 1j*g14*np.sin(3*phi)
        Ey = -1j*g12*np.cos(phi) - 1j*g14*np.cos(3*phi)
        Ez = 2*g13*np.sin(2*phi)
        E = np.array([Ex,Ey,Ez])*factor

        return E

    def H_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g10 = I10(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g11 = I11(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g12 = I12(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g13 = I13(k, self.theta_max, self.focal_length, self.width)(rho, z)
        g14 = I14(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)
    
        Hx = -1j*g12*np.cos(phi) - 1j*g14*np.cos(3*phi)
        Hy = 1j*g11*np.sin(phi) - 1j*g14*np.sin(3*phi)
        Hz = -2*g10 - 2*g13*np.cos(2*phi)
        H = np.array([Hx,Hy,Hz])*factor
    
        return H

class azimuthal_beam(source):
    def __init__(self, width, theta_max, focal_length):
        self.width = width
        self.theta_max = theta_max
        self.focal_length = focal_length

    def E_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        gazi = Iazi(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)

        Ex = 1j*gazi*np.sin(phi)
        Ey = -1j*gazi*np.cos(phi)
        Ez = np.zeros_like(Ex)
        E = np.array([Ex,Ey,Ez])*factor

        return E

    def H_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g10  = I10(k , self.theta_max, self.focal_length, self.width)(rho, z)
        grad = Irad(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)

        Hx = 1j*grad*np.cos(phi)
        Hy = 1j*grad*np.sin(phi)
        Hz = -4*g10
        H = np.array([Hx,Hy,Hz])*factor

        return H

class radial_beam(source):
    def __init__(self, width, theta_max, focal_length):
        self.width = width
        self.theta_max = theta_max
        self.focal_length = focal_length

    def E_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g10  = I10(k,  self.theta_max, self.focal_length, self.width)(rho, z)
        grad = Irad(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)

        Ex = 1j*grad*np.cos(phi)
        Ey = 1j*grad*np.sin(phi)
        Ez = -4*g10
        E = np.array([Ex,Ey,Ez])*factor
        E = np.array([Ex,Ey,Ez])*factor

        return E

    def H_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        gazi = Iazi(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)

        Hx = 1j*gazi*np.sin(phi)
        Hy = 1j*gazi*np.cos(phi)
        Hz = np.zeros_like(Hx)
        H = np.array([Hx,Hy,Hz])*factor

        return H

#TODO: fix
class shear_beam(source):
    def __init__(self, width, theta_max, focal_length):
        self.width = width
        self.theta_max = theta_max
        self.focal_length = focal_length

    def E_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        gazi = Iazi(k, self.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)

        Ex = 1j*gazi*np.sin(phi)
        Ey = 1j*gazi*np.cos(phi)
        Ez = 0*phi
        E = np.array([Ex,Ey,Ez])*factor

        return E

    def H_field(self, x, y, z, k):
        rho, phi, z = miepy.coordinates.cart_to_cyl(x, y, z)

        g10  = I10(k,  sel.theta_max, self.focal_length, self.width)(rho, z)
        grad = Irad(k, sel.theta_max, self.focal_length, self.width)(rho, z)

        factor = 1j*k*self.focal_length/2*np.exp(-1j*k*self.focal_length)

        Hx = 1j*grad*np.cos(phi)
        Hy = -1j*grad*np.sin(phi)
        Hz = -4*g10
        H = np.array([Hx,Hy,Hz])*factor

        return H
