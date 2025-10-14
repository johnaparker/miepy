import matplotlib.pyplot as plt
import meep
import meep_ext
import numpy as np
from numpipe import pbar, scheduler

import miepy

job = scheduler()
nm = 1e-9
um = 1e-6

### Parameters
radius = 75 * nm
height = 110 * nm
material = meep_ext.material.Au()

x = np.linspace(-220 * nm, 220 * nm, 3)
y = np.linspace(-220 * nm, 220 * nm, 3)
z = np.linspace(-40 * nm, 40 * nm, 9)
X, Y = np.meshgrid(x, y)
x = X.flatten()
y = Y.flatten()

theta = np.linspace(0, np.pi / 2, 9)
phi = np.linspace(0, 2 * np.pi, 9)

geometry = []
for i in range(9):
    axis = meep.Vector3(np.sin(theta[i]) * np.cos(phi[i]), np.sin(theta[i]) * np.sin(phi[i]), np.cos(theta[i]))
    geometry.append(
        meep.Cylinder(center=meep.Vector3(x[i], y[i], z[i]), radius=radius, height=height, material=material, axis=axis)
    )

box = [2 * radius + 2 * np.max(x)] * 2 + [2 * radius + 2 * np.max(z)]
resolution = 1 / (8 * nm)
medium = meep.Medium(index=1)

fcen, df = meep_ext.freq_data(1 / (400 * nm), 1 / (1000 * nm))
nfreq = 40
polarization = "x"

src_time = meep.GaussianSource(frequency=1.3 / um, fwidth=4.0 / um)
if polarization == "x":
    def source(sim):
        return meep_ext.x_polarized_plane_wave(sim, src_time)
    decay = meep.Ex
else:
    def source(sim):
        return meep_ext.y_polarized_plane_wave(sim, src_time)
    decay = meep.Ey

### monitor info
particle_monitor_gap = 50 * nm
pml_monitor_gap = 50 * nm
norm_file_ext = "norm_box"
monitor_size = [box[0] + 2 * particle_monitor_gap, box[1] + 2 * particle_monitor_gap, box[2] + 2 * particle_monitor_gap]

### grid
pml = meep.PML(15 / resolution)
lx = box[0] + 2 * particle_monitor_gap + 2 * pml_monitor_gap + 2 * pml.thickness
ly = box[1] + 2 * particle_monitor_gap + 2 * pml_monitor_gap + 2 * pml.thickness
lz = box[2] + 2 * particle_monitor_gap + 2 * pml_monitor_gap + 2 * pml.thickness
cell = meep.Vector3(lx, ly, lz)
Nx, Ny, Nz = map(round, cell * resolution)


@job.cache
def norm_sim():
    """Perform normalization simulation."""
    norm = meep.Simulation(
        cell_size=cell, boundary_layers=[pml], geometry=[], default_material=medium, resolution=resolution
    )
    norm.init_fields()
    source(norm)

    flux_box_inc = meep_ext.add_flux_box(norm, fcen, df, nfreq, [0, 0, 0], monitor_size)
    flux_inc = meep_ext.add_flux_plane(norm, fcen, df, nfreq, [0, 0, 0], [box[0], box[1], 0])

    norm.run(
        until_after_sources=meep.stop_when_fields_decayed(
            0.5 * um, decay, pt=meep.Vector3(0, 0, monitor_size[2] / 2), decay_by=1e-3
        )
    )

    norm.save_flux(norm_file_ext, flux_box_inc)

    return {
        "frequency": np.array(meep.get_flux_freqs(flux_inc)),
        "area": box[0] * box[1],
        "incident": np.asarray(meep.get_fluxes(flux_inc)),
    }


@job.cache
def scat_sim():
    """Perform scattering simulation."""
    scat = meep.Simulation(
        cell_size=cell, boundary_layers=[pml], geometry=geometry, default_material=medium, resolution=resolution
    )
    scat.init_fields()
    source(scat)

    flux_box_absorb = meep_ext.add_flux_box(scat, fcen, df, nfreq, [0, 0, 0], monitor_size)
    flux_box_scat = meep_ext.add_flux_box(scat, fcen, df, nfreq, [0, 0, 0], monitor_size)
    scat.load_minus_flux(norm_file_ext, flux_box_scat)

    scat.run(
        until_after_sources=meep.stop_when_fields_decayed(
            0.5 * um, decay, pt=meep.Vector3(0, 0, monitor_size[2] / 2), decay_by=1e-5
        )
    )

    return {
        "scattering": np.array(meep.get_fluxes(flux_box_scat)),
        "absorption": -np.array(meep.get_fluxes(flux_box_absorb)),
        "frequency": np.array(meep.get_flux_freqs(flux_box_scat)),
    }


@job.plots
def vis():
    fig, ax = plt.subplots()

    norm = job.load(norm_sim)
    scat = job.load(scat_sim)

    ax.plot(
        (1 / nm) / norm.frequency,
        scat.scattering / norm.incident * norm.area,
        "o",
        color="C0",
        label="scattering (FDTD)",
    )
    ax.plot(
        (1 / nm) / norm.frequency,
        scat.absorption / norm.incident * norm.area,
        "o",
        color="C1",
        label="absorption (FDTD)",
    )
    ax.plot(
        (1 / nm) / norm.frequency,
        (scat.scattering + scat.absorption) / norm.incident * norm.area,
        "o",
        color="C2",
        label="extinction (FDTD)",
    )

    wavelengths = np.linspace(400 * nm, 1000 * nm, 100)
    eps = meep_ext.get_eps(material)(wavelengths)
    Au = miepy.data_material(wavelengths, eps)
    # Au = miepy.constant_material(3.5**2)

    C, A, E = [np.zeros_like(wavelengths) for i in range(3)]

    particles = []
    for i in range(9):
        orientation = miepy.quaternion.from_spherical_coords(theta[i], phi[i])
        particles.append(
            miepy.cylinder(
                position=[x[i], y[i], z[i]], radius=radius, height=height, material=Au, orientation=orientation
            )
        )

    for i, wavelength in enumerate(pbar(wavelengths)):
        sol = miepy.cluster(particles=particles, source=miepy.sources.plane_wave([1, 0]), wavelength=wavelength, lmax=3)

        C[i], A[i], E[i] = sol.cross_sections()

    ax.axhline(0, linestyle="--", color="black")
    ax.plot(wavelengths / nm, C, color="C0", label="scattering (GMT)")
    ax.plot(wavelengths / nm, A, color="C1", label="absorption (GMT)")
    ax.plot(wavelengths / nm, E, color="C2", label="extinction (GMT)")

    ax.legend()
    ax.set(xlabel="wavelength (nm)", ylabel="cross-section")

    plt.show()


job.run()
