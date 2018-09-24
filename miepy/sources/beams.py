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

def power_numeric(beam, k):
    """Calculate the total power in a beam via numeric integration"""
    radius = 1e6*2*np.pi/k
    theta = np.linspace(np.pi/2, np.pi, 30)
    phi = np.linspace(0, 2*np.pi, 20)
    THETA, PHI = np.meshgrid(theta, phi)
    R = radius*np.ones_like(THETA)

    E = beam.scalar_potenital_ingoing(THETA, PHI, k, norm=False)
    S = 0.5/Z0*np.abs(E)**2*np.sin(THETA)
    P = radius**2*miepy.vsh.misc.trapz_2d(theta, phi, S.T).real/4

    return P

def structure_function(src, k, lmax):
    #TODO implement a better way of finding the maximum value... per source object
    f = lambda theta: np.linalg.norm(src.spherical_ingoing(theta, 0, k))
    g = lambda theta: np.linalg.norm(src.spherical_ingoing(theta, np.pi/2, k))
    theta = np.linspace(np.pi/2, np.pi, 500)
    f_max = np.max(f(theta))
    g_max = np.max(g(theta))

    if g_max > f_max:
        cutoff = find_cutoff(lambda theta: g(theta)/g_max, 1e-6, tol=1e-9)
    else:
        cutoff = find_cutoff(lambda theta: f(theta)/f_max, 1e-6, tol=1e-9)

    return miepy.vsh.decomposition.integral_project_source_far(src, 
                           k, lmax, theta_0=cutoff)

#TODO: implement .from_string constructors
class beam(source):
    def __init__(self, polarization, theta=0, phi=0, power=None, 
                    amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(amplitude, phase)
        self.polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization /= np.linalg.norm(self.polarization)
        self.center = np.asarray(center)

        self.theta = theta
        self.phi   = phi
        self.orientation = miepy.quaternion.from_spherical_coords(self.theta, self.phi)

        ### TE and TM vectors
        self.k_hat, self.n_te, self.n_tm = miepy.coordinates.sph_basis_vectors(theta, phi)

        self.amplitude = amplitude
        self.power = power
        if power is None and amplitude is None:
            raise ValueError('either power or amplitude must be specified')
        elif power is not None and amplitude is not None:
            raise ValueError('cannot specify power and amplitude simultaneously')

        self.current_k = None

    def scalar_potenital(self, x, y, z, k): pass

    def scalar_potenital_ingoing(self, theta, phi, k): pass

    def is_paraxial(self, k): pass

    def E_field(self, x, y, z, k):
        xp, yp, zp = miepy.coordinates.translate(x, y, z, -self.center)
        xp, yp, zp = miepy.coordinates.rotate(xp, yp, zp, self.orientation.inverse())

        U = self.scalar_potenital(xp, yp, zp, k)
        E = np.zeros((3,) + U.shape, dtype=complex)
        amp = U*np.exp(1j*self.phase)

        pol = self.n_te*self.polarization[0] + self.n_tm*self.polarization[1]
        return np.einsum('i...,...->i...', pol, amp)

    def H_field(self, x, y, z, k):
        xp, yp, zp = miepy.coordinates.translate(x, y, z, -self.center)
        xp, yp, zp = miepy.coordinates.rotate(xp, yp, zp, -self.orientation)

        U = self.scalar_potenital(xp, yp, zp, k)
        H = np.zeros((3,) + U.shape, dtype=complex)
        amp = U*np.exp(1j*self.phase)

        pol = self.n_tm*self.polarization[0] - self.n_te*self.polarization[1]
        return np.einsum('i...,...->i...', pol, amp)

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
        #TODO: fix orientation_copy hack
        orientation_copy = self.orientation
        self.orientation = miepy.quaternion.one
        if self.is_paraxial(k):
            sampling = miepy.vsh.decomposition.sampling_from_lmax(lmax, method='near')
            p_src = miepy.vsh.decomposition.near_field_point_matching(self, 
                              position, 2*radius, k, lmax, sampling)
        else:
            if k != self.current_k:
                self.current_k = k
                self.p_src_func = structure_function(self, k, lmax)

            p_src = self.p_src_func(position)

        self.orientation = orientation_copy

        if self.theta != 0 or self.phi != 0:
            p_src = miepy.vsh.rotate_expansion_coefficients(p_src, self.orientation)

        return p_src

class paraxial_beam(beam):
    def __init__(self, Ufunc, polarization, theta=0, phi=0, 
                   power=None, amplitude=1, phase=0, center=np.zeros(3)):
        super().__init__(polarization, theta, phi, power, amplitude, phase, center)
        self.Ufunc = Ufunc

    def scalar_potenital(self, x, y, z, k):
        return self.Ufunc(x, y, z, k)

    def is_paraxial(self, k):
        return True

class gaussian_beam(beam):
    def __init__(self, width, polarization, theta=0, phi=0, 
                  power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, theta, phi, power, amplitude, phase, center)
        self.width = width
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            E0 = 2/self.width*np.sqrt(Z0*self.power/np.pi)
        else:
            E0 = self.amplitude

        rho_sq = x**2 + y**2
        wav = 2*np.pi/k
        amp = E0*self.width/w(z, self.width, wav) * np.exp(-rho_sq/w(z,self.width,wav)**2)
        phase = k*z[2] + k*rho_sq*Rinv(z[2],self.width,wav)/2 - gouy(z[2],self.width,wav)

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
        # return True
        wav = 2*np.pi/k
        return self.width > 4*wav

class bigaussian_beam(beam):
    def __init__(self, width_x, width_y, polarization, theta=0, phi=0,
                    power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, theta, phi, power, amplitude, phase, center)
        self.width_x = width_x
        self.width_y = width_y
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            E0 = 2/np.sqrt(self.width_x*self.width_y)*np.sqrt(Z0*self.power/np.pi)
        else:
            E0 = self.amplitude

        rho_sq = x**2/self.width_x**2 + y**2/self.width_y**2
        wav = 2*np.pi/k
        amp = E0*np.exp(-rho_sq)
        phase = k*z

        return amp*np.exp(1j*phase)

    def scalar_potenital_ingoing(self, theta, phi, k):
        pass

    def is_paraxial(self, k):
        return True
        # return 2*np.pi/k < min(self.width_x, self.width_y)

class hermite_gaussian_beam(beam):
    def __init__(self, l, m, width, polarization, theta=0, phi=0, 
                    power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, theta, phi, power, amplitude, phase, center)
        self.width = width
        self.l = l
        self.m = m
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            factor = 1/self.width*np.sqrt(2/(np.pi*2**self.l*2**self.m*factorial(self.l)*factorial(self.m)))
            E0 = factor*np.sqrt((2*Z0*self.power))
        else:
            E0 = self.amplitude

        rho_sq = x**2 + y**2
        wav = 2*np.pi/k

        wz = w(z, self.width, wav)
        HG_l = eval_hermite(self.l, np.sqrt(2)*x/wz)
        HG_m = eval_hermite(self.m, np.sqrt(2)*y/wz)
        N = self.l + self.m

        amp = E0*self.width/wz * HG_l * HG_m * np.exp(-rho_sq/wz**2)
        phase = k*z + k*rho_sq*Rinv(z,self.width,wav)/2 - (N+1)*gouy(z,self.width,wav)

        return amp*np.exp(1j*phase)

    def scalar_potenital_ingoing(self, theta, phi, k, norm=True):
        #TODO: remove the norm hack
        if norm:
            if self.amplitude is None:
                U0 = np.sqrt(self.power/power_numeric(self, k))
            else:
                U0 = self.amplitude
        else:
            U0 = 1

        wav = 2*np.pi/k
        r = 1e6*wav
        x, y, z = miepy.coordinates.sph_to_cart(r, theta, phi)
        
        rho_sq = x**2 + y**2

        wz = w(z, self.width, wav)
        HG_l = eval_hermite(self.l, np.sqrt(2)*x/wz)
        HG_m = eval_hermite(self.m, np.sqrt(2)*y/wz)
        N = self.l + self.m

        amp = U0*self.width/wz * HG_l * HG_m * np.exp(-rho_sq/wz**2)

        return amp

    def is_paraxial(self, k):
        wav = 2*np.pi/k
        return self.width > 4*wav

class laguerre_gaussian_beam(beam):
    def __init__(self, p, l, width, polarization, theta=0, phi=0,
                    power=None, amplitude=None, phase=0, center=np.zeros(3)):
        super().__init__(polarization, theta, phi, power, amplitude, phase, center)
        self.width = width
        self.p = p
        self.l = l
    
    def scalar_potenital(self, x, y, z, k):
        if self.amplitude is None:
            E0 = np.sqrt((2*Z0*self.power))
        else:
            E0 = self.amplitude

        rho_sq = x**2 + y**2
        phi = np.arctan2(y, x)
        wav = 2*np.pi/k

        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        wz = w(z, self.width, wav)

        Lpl = eval_genlaguerre(self.p, abs(self.l), 2*rho_sq/wz**2)
        N = abs(self.l) + 2*self.p

        amp = E0*C/wz * np.exp(-rho_sq/wz**2) * ((2*rho_sq)**0.5/wz)**abs(self.l) * Lpl
        phase = self.l*phi + k*z + k*rho_sq*Rinv(z,self.width,wav)/2 - (N+1)*gouy(z,self.width,wav)

        return amp*np.exp(1j*phase)

    def scalar_potenital_ingoing(self, theta, phi, k):
        if self.amplitude is None:
            E0 = np.sqrt((2*Z0*self.power))
        else:
            E0 = self.amplitude

        wav = 2*np.pi/k
        r = 1e6*wav
        x, y, z = miepy.coordinates.sph_to_cart(r, theta, phi)

        rho_sq = x**2 + y**2
        phi = np.arctan2(y, x)

        C = np.sqrt(2*factorial(self.p)/(np.pi*factorial(self.p + abs(self.l))))
        wz = w(z, self.width, wav)

        Lpl = eval_genlaguerre(self.p, abs(self.l), 2*rho_sq/wz**2)
        N = abs(self.l) + 2*self.p

        amp = E0*C/wz * np.exp(-rho_sq/wz**2) * ((2*rho_sq)**0.5/wz)**abs(self.l) * Lpl
        phase = self.l*phi

        return amp*np.exp(1j*phase)

    def is_paraxial(self, k):
        wav = 2*np.pi/k
        return self.width > 4*wav

def azimuthal_beam(width, theta=0, phi=0, amplitude=None, power=None, phase=0, center=np.zeros(3)):
    """azimuthally polarized beam"""
    if power is not None:
        power /= 2
    HG_1 = hermite_gaussian_beam(1, 0, width, [0,1],  theta=theta, phi=phi, 
                  amplitude=amplitude, power=power, phase=phase, center=center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [-1,0], theta=theta, phi=phi,
                  amplitude=amplitude, power=power, phase=phase, center=center)
    return HG_1 + HG_2

def radial_beam(width, theta=0, phi=0, amplitude=None, power=None, phase=0, center=np.zeros(3)):
    """radially polarized beam"""
    if power is not None:
        power /= 2
    HG_1 = hermite_gaussian_beam(1, 0, width, [1,0], theta=theta, phi=phi,
                  amplitude=amplitude, power=power, phase=phase, center=center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [0,1], theta=theta, phi=phi,
                  amplitude=amplitude, power=power, phase=phase, center=center)
    return HG_1 + HG_2

def shear_beam(width, theta=0, phi=0, amplitude=None, power=None, phase=0, center=np.zeros(3)):
    """shear polarized beam"""
    if power is not None:
        power /= 2
    HG_1 = hermite_gaussian_beam(1, 0, width, [0,1], theta=theta, phi=phi,
                  amplitude=amplitude, power=power, phase=phase, center=center)
    HG_2 = hermite_gaussian_beam(0, 1, width, [1,0], theta=theta, phi=phi,
                  amplitude=amplitude, power=power, phase=phase, center=center)
    return HG_1 + HG_2
