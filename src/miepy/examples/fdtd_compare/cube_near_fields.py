"""Compare the near fields of a cube in FDTD and GMMT."""

import matplotlib.pyplot as plt
import meep
import meep_ext
import numpy as np
from numpipe import scheduler

import miepy

job = scheduler()
nm = 1e-9
um = 1e-6

### Parameters
W = 250 * nm
H = 120 * nm

material = meep_ext.material.Au()
material = meep.Medium(index=3)

geometry = []
q = miepy.quaternion.from_spherical_coords(0, 0)
R = miepy.quaternion.as_rotation_matrix(q)

X = W / np.sqrt(2 * (1 - np.cos(2 * np.pi / 3)))
vertices = [
    meep.Vector3(X, 0, 0),
    meep.Vector3(-X / 2, np.sqrt(3) / 2 * X, 0),
    meep.Vector3(-X / 2, -np.sqrt(3) / 2 * X, 0),
]

vertices = [
    meep.Vector3(W / 2, -W / 2, 0),
    meep.Vector3(W / 2, +W / 2, 0),
    meep.Vector3(-W / 2, +W / 2, 0),
    meep.Vector3(-W / 2, -W / 2, 0),
]

geometry.append(
    meep.Prism(center=meep.Vector3(0, 0, 0), height=H, vertices=vertices, axis=meep.Vector3(0, 0, 1), material=material)
)


L = W + 400 * nm
box = [L, L, L]
medium = meep.Medium(index=1)
resolution = 1 / (8 * nm)

fcen, df = meep_ext.freq_data(1 / (400 * nm), 1 / (1000 * nm))
polarization = "x"

src_time = meep.GaussianSource(frequency=fcen, fwidth=df)
def source(sim):
    return meep_ext.rhc_polarized_plane_wave(sim, src_time)
decay = meep.Ex

### monitor info
gap = 80 * nm

### grid
pml = meep.PML(15 / resolution)
lx = box[0] + 2 * gap
ly = box[1] + 2 * gap
lz = box[2] + 2 * gap
cell = meep.Vector3(lx, ly, lz)
Nx, Ny, Nz = map(round, cell * resolution)


@job.cache
def sim():
    """Perform scattering simulation."""
    sim = meep.Simulation(
        cell_size=cell, boundary_layers=[pml], geometry=geometry, default_material=medium, resolution=resolution
    )
    sim.init_fields()
    sim.init_sim()
    source(sim)

    freq = 1 / (600 * nm)
    # dft = sim.add_dft_fields([meep.Ex, meep.Ey, meep.Ez], freq, freq, 1, center=meep.Vector3(0,0,0),
    # size=meep.Vector3(L,L,0))

    vol = meep.volume(meep.vec(-L / 2, -L / 2, 0), meep.vec(L / 2, L / 2, 0))
    dft = sim.fields.add_dft_fields([meep.Ex, meep.Ey, meep.Ez], vol, freq, freq, 1)

    sim.run(
        until_after_sources=meep.stop_when_fields_decayed(
            0.5 * um, decay, pt=meep.Vector3(0, 0, box[2] / 2), decay_by=1e-5
        )
    )
    # for i in range(100):
    # sim.fields.step()

    Ex = sim.get_dft_array(dft, meep.Ex, 0)
    Ey = sim.get_dft_array(dft, meep.Ey, 0)
    Ez = sim.get_dft_array(dft, meep.Ez, 0)
    E = np.array([Ex, Ey, Ez])

    return dict(E=E)


@job.cache
def gmt_sim():
    material = miepy.materials.constant_material(3**2)

    q = miepy.quaternion.from_spherical_coords(0, 0)
    particles = [miepy.regular_prism([0, 0, 0], 4, width=W, height=H, material=material, orientation=q)]

    c = miepy.cluster(particles=particles, source=miepy.sources.plane_wave([1, 1j]), wavelength=600 * nm, lmax=8)

    x = np.linspace(-1.5 * W, 1.5 * W, 150)
    y = np.linspace(-1.5 * W, 1.5 * W, 150)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    E = c.E_field(X, Y, Z, mask=True)
    Ixy = np.sum(np.abs(E) ** 2, axis=0)

    R = particles[0].enclosed_radius() + 20 * nm
    phi = np.linspace(0, 2 * np.pi, 50)
    x = R * np.cos(phi)
    y = R * np.sin(phi)
    z = np.zeros_like(x)
    E = c.E_field(x, y, z)
    Iphi = np.sum(np.abs(E) ** 2, axis=0)

    return dict(phi=phi, X=X, Y=Y, Ixy=Ixy, Iphi=Iphi)


@job.plots
def vis():
    fdtd = job.load(sim)
    gmmt = job.load(gmt_sim)

    fig, axes = plt.subplots(ncols=2)
    I = np.sum(np.abs(fdtd.E) ** 2, axis=0)

    x = np.linspace(-L / 2, L / 2, fdtd.E.shape[1])
    y = np.linspace(-L / 2, L / 2, fdtd.E.shape[2])
    X, Y = np.meshgrid(x, y, indexing="ij")
    axes[0].pcolormesh(X / nm, Y / nm, I, shading="gouraud", vmax=np.max(I) / 10)
    axes[0].set_title("FDTD")
    axes[1].pcolormesh(gmmt.X / nm, gmmt.Y / nm, gmmt.Ixy, shading="gouraud")
    axes[1].set_title("GMMT")

    for ax in axes:
        ax.set_aspect("equal")

    from scipy.interpolate import interp2d

    f = interp2d(X, Y, I)

    R = 206 * nm
    phi = np.linspace(0, 2 * np.pi, 50)
    x = R * np.cos(phi)
    y = R * np.sin(phi)
    I = np.zeros_like(x)
    for i in range(len(phi)):
        I[i] = f(x[i], y[i])

    fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
    ax.plot(gmmt.phi, gmmt.Iphi / np.max(gmmt.Iphi), label="GMMT")
    ax.plot(phi, I / np.max(I), "o", label="FDTD")

    ax.legend()

    plt.show()


job.run()
