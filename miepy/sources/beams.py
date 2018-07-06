"""
beam sources
"""

import numpy as np
import miepy
from miepy.sources.source_base import source, combined_source
from math import factorial
from scipy.special import eval_genlaguerre, eval_hermite, erfc
from scipy.constants import physical_constants

Z0 = physical_constants['characteristic impedance of vacuum'][0]

def zr(w0, wav):
    return np.pi*w0**2/wav

def w(z, w0, wav):
    return w0*np.sqrt(1 + (z/zr(w0,wav))**2)
    
def Rinv(z, w0, wav):
    return z/(z**2 + zr(w0,wav)**2)

def gouy(z, w0, wav):
    return np.arctan2(z, zr(w0,wav))

def find_cutoff(f, cutoff, tol=1e-7):
    theta = np.pi/2
    val = f(theta)
    dtheta = 0.01*np.pi
    err = np.abs(val - cutoff)
    
    while err > tol:
        val = 0
        while val < cutoff:
            theta += dtheta
            val = f(theta)

        err = abs(val - cutoff)
        theta -= dtheta
        dtheta /= 2

    return theta

class beam(source):
    def __init__(self, polarization, power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(amplitude)
        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
        self.center = np.asarray(center)
        self.phase = phase

        self.amplitude = amplitude
        self.power = power
        if power is None and amplitude is None:
            raise ValueError('either power or amplitude must be specified')
        elif power is not None and amplitude is not None:
            raise ValueError('cannot specify power and amplitude simultaneously')

    def scalar_potenital(self, x, y, z, k): pass

    def scalar_potenital_ingoing(self, theta, phi, k): pass

    def is_paraxial(self, x, y, z, k): pass

    def E_field(self, x, y, z, k):
        U = self.scalar_potenital(x, y, z, k)
        E = np.zeros((3,) + U.shape, dtype=complex)
        E[0] = U*self.polarization[0]
        E[1] = U*self.polarization[1]
        E *= np.exp(1j*self.phase)

        return E

    def H_field(self, x, y, z, k):
        U = self.scalar_potenital(x, y, z, k)
        H = np.zeros((3,) + U.shape, dtype=complex)
        H[0] = -U*self.polarization[1]
        H[1] = U*self.polarization[0]
        H *= np.exp(1j*self.phase)

        return H

    def spherical_ingoing(self, theta, phi, k):
        """Determine far-field spherical fields from the scalar potential.
        Returns E[2,...], the theta & phi components"""
        U = self.scalar_potenital_ingoing(theta, phi, k)
        Esph = np.zeros((2,) + U.shape, dtype=complex)

        Ex = U*self.polarization[0]
        Ey = U*self.polarization[1]
        Esph[0] = -Ex*np.cos(phi) - Ey*np.sin(phi)
        Esph[1] = -Ex*np.sin(phi) + Ey*np.cos(phi)
        Esph *= np.exp(1j*self.phase)

        return Esph

    def structure(self, position, k, lmax, radius):
        if self.is_paraxial(k):
            sampling = miepy.vsh.decomposition.sampling_from_lmax(lmax, method='near')
            return miepy.vsh.decomposition.near_field_point_matching(self, 
                              position, 2*radius, k, lmax, sampling)
        else:
            #TODO implement a better way of finding the maximum value... per source object
            f = lambda theta: np.linalg.norm(self.spherical_ingoing(theta, 0, k))
            g = lambda theta: np.linalg.norm(self.spherical_ingoing(theta, np.pi/2, k))
            theta = np.linspace(np.pi/2, np.pi, 500)
            f_max = np.max(f(theta))
            g_max = np.max(g(theta))

            if g_max > f_max:
                cutoff = find_cutoff(lambda theta: g(theta)/g_max, 1e-6, tol=1e-9)
            else:
                cutoff = find_cutoff(lambda theta: f(theta)/f_max, 1e-6, tol=1e-9)

            return miepy.vsh.decomposition.integral_project_source_far(self, 
                              k, lmax, origin=position, theta_0=cutoff)

class paraxial_beam(beam):
    def __init__(self, Ufunc, polarization, power=None, amplitude=1, phase=0, center=np.zeros(3)):
        super().__init__(polarization, power, amplitude, phase, center)
        self.Ufunc = Ufunc

    def scalar_potenital(self, x, y, z, k):
        return self.Ufunc(x, y, z, k)

    def is_paraxial(self, k):
        return True

class gaussian_beam(beam):
    def __init__(self, width, polarization, power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, power, amplitude, phase, center)
        self.width = width
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            E0 = 2/self.width*np.sqrt(Z0*self.power/np.pi)
        else:
            E0 = self.amplitude

        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k
        amp = E0*self.width/w(rp[2], self.width, wav) * np.exp(-rho_sq/w(rp[2],self.width,wav)**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - gouy(rp[2],self.width,wav)

        return amp*np.exp(1j*phase)

    def scalar_potenital_ingoing(self, theta, phi, k):
        r = 1e6*(2*np.pi/k)
        if self.amplitude is None:
            c = 0.5*(k*self.width)**2
            if self.is_paraxial(k):
                U0 = 2*k*self.width*np.sqrt(Z0*self.power/np.pi)/r
            else:
                U0 = 2*np.sqrt(Z0*self.power/(np.pi*(1 - np.sqrt(np.pi*c)*np.exp(c)*erfc(np.sqrt(c)))))/r
        else:
            U0 = k*self.width**2*self.amplitude/r/2

        U = U0*np.exp(-(k*self.width*np.tan(theta)/2)**2)
        return U

    def is_paraxial(self, k):
        wav = 2*np.pi/k
        return self.width > 4*wav

class bigaussian_beam(beam):
    def __init__(self, width_x, width_y, polarization, power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, power, amplitude, phase, center)
        self.width_x = width_x
        self.width_y = width_y
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            E0 = 2/np.sqrt(self.width_x*self.width_y)*np.sqrt(Z0*self.power/np.pi)
        else:
            E0 = self.amplitude

        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2/self.width_x**2 + rp[1]**2/self.width_y**2
        wav = 2*np.pi/k
        amp = E0*np.exp(-rho_sq)
        phase = k*rp[2]

        return amp*np.exp(1j*phase)

    def scalar_potenital_ingoing(self, theta, phi, k):
        pass

    def is_paraxial(self, k):
        return True
        # return 2*np.pi/k < min(self.width_x, self.width_y)

class hermite_gaussian_beam(beam):
    def __init__(self, l, m, width, polarization, power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, power, amplitude, phase, center)
        self.width = width
        self.l = l
        self.m = m
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            factor = 1/self.width*np.sqrt(2/(np.pi*2**self.l*2**self.m*factorial(self.l)*factorial(self.m)))
            E0 = factor*np.sqrt((2*Z0*self.power))
        else:
            E0 = self.amplitude

        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        wav = 2*np.pi/k

        wz = w(rp[2], self.width, wav)
        HG_l = eval_hermite(self.l, np.sqrt(2)*rp[0]/wz)
        HG_m = eval_hermite(self.m, np.sqrt(2)*rp[1]/wz)
        N = self.l + self.m

        amp = E0*self.width/wz * HG_l * HG_m * np.exp(-rho_sq/wz**2)
        phase = k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - (N+1)*gouy(rp[2],self.width,wav)

        return amp*np.exp(1j*phase)

    def scalar_potenital_ingoing(self, theta, phi, k):
        wav = 2*np.pi/k
        r = 1e6*wav
        x, y, z = miepy.coordinates.sph_to_cart(r, theta, phi)
        
        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2

        wz = w(rp[2], self.width, wav)
        HG_l = eval_hermite(self.l, np.sqrt(2)*rp[0]/wz)
        HG_m = eval_hermite(self.m, np.sqrt(2)*rp[1]/wz)
        N = self.l + self.m

        amp = self.amplitude*self.width/wz * HG_l * HG_m * np.exp(-rho_sq/wz**2)

        return amp

    def is_paraxial(self, k):
        wav = 2*np.pi/k
        return self.width > 4*wav

class laguerre_gaussian_beam(beam):
    def __init__(self, p, l, width, polarization, power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, power, amplitude, phase, center)
        self.width = width
        self.p = p
        self.l = l
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            E0 = np.sqrt((2*Z0*self.power))
        else:
            E0 = self.amplitude

        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        phi = np.arctan2(rp[1], rp[0])
        wav = 2*np.pi/k

        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        wz = w(rp[2], self.width, wav)

        Lpl = eval_genlaguerre(self.p, abs(self.l), 2*rho_sq/wz**2)
        N = abs(self.l) + 2*self.p

        amp = E0*C/wz * np.exp(-rho_sq/wz**2) * ((2*rho_sq)**0.5/wz)**abs(self.l) * Lpl
        phase = self.l*phi + k*rp[2] + k*rho_sq*Rinv(rp[2],self.width,wav)/2 - (N+1)*gouy(rp[2],self.width,wav)

        return amp*np.exp(1j*phase)

    def scalar_potenital_ingoing(self, theta, phi, k):
        if self.amplitude is None:
            E0 = np.sqrt((2*Z0*self.power))
        else:
            E0 = self.amplitude

        wav = 2*np.pi/k
        r = 1e6*wav
        x, y, z = miepy.coordinates.sph_to_cart(r, theta, phi)

        rp = np.array([x - self.center[0], y - self.center[1], z - self.center[2]])
        rho_sq = rp[0]**2 + rp[1]**2
        phi = np.arctan2(rp[1], rp[0])

        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        wz = w(rp[2], self.width, wav)

        Lpl = eval_genlaguerre(self.p, abs(self.l), 2*rho_sq/wz**2)
        N = abs(self.l) + 2*self.p

        amp = E0*C/wz * np.exp(-rho_sq/wz**2) * ((2*rho_sq)**0.5/wz)**abs(self.l) * Lpl
        phase = self.l*phi

        return amp*np.exp(1j*phase)

    def is_paraxial(self, k):
        wav = 2*np.pi/k
        return self.width > 4*wav

def azimuthal_beam(width, amplitude=None, power=None, phase=0, center=np.zeros(3)):
    """azimuthally polarized beam"""
    if power is not None:
        power /= 2
    HG_1 = hermite_gaussian_beam(1, 0, width, [0,1],  amplitude=amplitude, power=power, phase=phase, center=center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [-1,0], amplitude=amplitude, power=power, phase=phase, center=center)
    return HG_1 + HG_2

def radial_beam(width, amplitude=None, power=None, phase=0, center=np.zeros(3)):
    """radially polarized beam"""
    if power is not None:
        power /= 2
    HG_1 = hermite_gaussian_beam(1, 0, width, [1,0], amplitude=amplitude, power=power, phase=phase, center=center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [0,1], amplitude=amplitude, power=power, phase=phase, center=center)
    return HG_1 + HG_2

def shear_beam(width, amplitude=None, power=None, phase=0, center=np.zeros(3)):
    """shear polarized beam"""
    if power is not None:
        power /= 2
    HG_1 = hermite_gaussian_beam(1, 0, width, [0,1], amplitude=amplitude, power=power, phase=phase, center=center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [1,0], amplitude=amplitude, power=power, phase=phase, center=center)
    return HG_1 + HG_2
