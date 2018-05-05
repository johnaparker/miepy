import numpy as np
import matplotlib.pyplot as plt
from my_pytools.my_numpy.integrate import simps_2d
from tqdm import tqdm
import h5py
import miepy

def sph_to_cart(r, theta, phi, origin=[0,0,0]):
    """convert spherical coordinates (r, theta, phi) centered at origin to cartesian coordinates (x, y, z)"""
    x = origin[0] + r*np.sin(theta)*np.cos(phi)
    y = origin[1] + r*np.sin(theta)*np.sin(phi)
    z = origin[2] + r*np.cos(theta)

    return x,y,z

def sphere_mesh(sampling):
    """
    Obtain a THETA,PHI mesh for discretizing the surface of the sphere, consistent
    with the format required by the project and decompose functions
    Returns (THETA,PHI) meshgrids

    Arguments:
        sampling   number of points to sample between 0 and pi
    """

    phi = np.linspace(0, 2*np.pi, 2*sampling)
    tau = np.linspace(-1, 1, sampling)
    theta = np.arccos(tau)

    THETA,PHI = np.meshgrid(theta, phi, indexing='ij')
    return THETA, PHI

nm = 1e-9
sampling = 40
Nwav = 200

Ag = miepy.materials. Ag()
radius = 75*nm
source = miepy.sources.x_polarized_plane_wave(amplitude=1)
# separations = np.linspace(2*radius+10e-9,2*radius+700e-9, 50)
separation = 153*nm
wavelengths = np.linspace(330*nm, 1000*nm, Nwav)

separations = np.linspace(153*nm, 300*nm, 10)

plt.figure(figsize=(5,10))
for sep_idx,separation in enumerate(tqdm(separations)):
    spheres = miepy.spheres([[separation/2,0,0], [-separation/2,0,0]], radius, Ag)
    sol = miepy.gmt(spheres, source, wavelengths, 3, interactions=True)

    # sol.update_position(np.array([[separation/2,0,0], [-separation/2,0,0]]))

    r = 10000*nm
    THETA,PHI = sphere_mesh(sampling)
    X,Y,Z = sph_to_cart(r, THETA, PHI)

    E = sol.E_field(X,Y,Z, inc=False)
    # H = sol.H_field(X,Y,Z, inc=False)
    I = np.sum(np.abs(E)**2, axis=0)

    dA = r**2
    tau = np.linspace(-1, 1, sampling)
    phi = np.linspace(0, 2*np.pi, 2*sampling)
    flux = np.zeros_like(wavelengths)

    for i in range(Nwav):
        flux[i] = simps_2d(tau, phi, I[i])

    plt.plot(wavelengths/nm, flux - .0004*sep_idx, color='C0')

plt.xlabel('wavelength (nm)')
plt.ylabel('scattering')
plt.savefig('out.pdf', bbox_inches='tight')
plt.show()

