import miepy
import numpy as np

class interface:
    def __init__(self, material, z=0):
        self.material = material
        self.z = z

    def get_relative_index(self, wavelength, medium):
        eps_m = medium.eps(wavelength)
        eps = self.material.eps(wavelength)
        return np.sqrt(complex(eps/eps_m))

    def reflection_coefficients(self, theta, wavelength, medium):
        m = self.get_relative_index(wavelength, medium)
        theta_t = self.transmission_angle(theta, wavelength, medium)
        r_parallel = (m*np.cos(theta) - np.cos(theta_t))/(m*np.cos(theta) + np.cos(theta_t))
        r_perp     = (np.cos(theta) - m*np.cos(theta_t))/(np.cos(theta) + m*np.cos(theta_t))

        return r_parallel, r_perp

    #TODO: implement
    def transmission_kz(self, theta, wavelength, medium):
        pass

    def transmission_angle(self, theta, wavelength, medium):
        m = self.get_relative_index(wavelength, medium)
        theta_t = np.arcsin(np.sin(theta)/m)
        return theta_t

    def transmission_coefficients(self, theta, wavelength, medium):
        m = self.get_relative_index(wavelength, medium)
        theta_t = self.transmission_angle(theta, wavelength, medium)
        t_parallel = 2*np.cos(theta)/(m*np.cos(theta) + np.cos(theta_t))
        t_perp     = 2*np.cos(theta)/(np.cos(theta) + m*np.cos(theta_t))

        return t_parallel, t_perp

if __name__ == '__main__':
    n1 = 1
    n2 = 500 + 2j
    wavelength = 1
    k1 = n1*2*np.pi/wavelength
    k2 = n2*2*np.pi/wavelength
    medium = miepy.constant_material(index=n1)
    interface = miepy.interface(material=miepy.constant_material(index=n2))

    incident = miepy.sources.plane_wave([1,0], 0)
    reflected = interface.reflected_plane_wave(incident, wavelength, medium)
    transmitted = interface.transmitted_plane_wave(incident, wavelength, medium)

    x = np.linspace(-3, 3, 300)
    z = np.linspace(-3, 3, 300)
    X, Z = np.meshgrid(x, z)
    Y = np.zeros_like(X)

    E = np.zeros((3,) + X.shape, dtype=complex)

    idx = Z < 0
    E[:,idx] += incident.E_field(X[idx], Y[idx], Z[idx], k1)
    E[:,idx] += reflected.E_field(X[idx], Y[idx], Z[idx], k1)

    idx = Z > 0
    # E[:,idx] = transmitted.E_field(X[idx], Y[idx], Z[idx], k2)

    I = np.sum(np.abs(E)**2, axis=0)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.pcolormesh(X, Z, I, vmin=0)
    ax.set_aspect('equal')

    fig, ax = plt.subplots()
    ax.plot(x, I[:,0])
    ax.plot(x, 4*np.sin(k1*x)**2)

    fig, ax = plt.subplots()
    pos = [-1.5,0,0]
    lmax = 9
    p1 = incident.structure(pos, k1, lmax) 
    p2 = reflected.structure(pos, k1, lmax) 
    
    x = np.linspace(-1, 1, 100)
    z = np.linspace(-1, 1, 100)
    X, Z = np.meshgrid(x, z)
    Y = .1*np.ones_like(X)
    RAD, THETA, PHI = miepy.coordinates.cart_to_sph(X, Y, Z)

    E = miepy.expand_E(p1 + p2, k1, miepy.vsh_mode.incident)(RAD, THETA, PHI)
    I = np.sum(np.abs(E)**2, axis=0)
    ax.pcolormesh(X, Z, I, vmin=0)
    ax.set_aspect('equal')

    plt.show()

