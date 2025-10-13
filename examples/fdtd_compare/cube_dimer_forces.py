import matplotlib.pyplot as plt
import meep
import meep_ext
import numpy as np
from numpipe import pbar, scheduler
from scipy import constants

import miepy

job = scheduler()
nm = 1e-9
um = 1e-6

### Parameters
W = 200 * nm
sep = 330 * nm

material = meep_ext.material.Au()
material = meep.Medium(index=3.5)

geometry = []
q = miepy.quaternion.from_spherical_coords(0, np.pi / 4)
R = miepy.quaternion.as_rotation_matrix(q)

geometry.append(
    meep.Block(
        center=meep.Vector3(-sep / 2, 0, 0),
        size=meep.Vector3(W, W, W),
        e1=meep.Vector3(*R[:, 0]),
        e2=meep.Vector3(*R[:, 1]),
        e3=meep.Vector3(*R[:, 2]),
        material=material,
    )
)
geometry.append(
    meep.Block(
        center=meep.Vector3(sep / 2, 0, 0),
        size=meep.Vector3(W, W, W),
        e1=meep.Vector3(*R[:, 0]),
        e2=meep.Vector3(*R[:, 1]),
        e3=meep.Vector3(*R[:, 2]),
        material=material,
    )
)


box = [sep + 2**0.5 * W, 2**0.5 * W, 2**0.5 * W]
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
particle_monitor_gap = 30 * nm
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
def force_norm():
    """Perform normalization simulation."""
    norm = meep.Simulation(cell_size=cell, boundary_layers=[pml], geometry=[], resolution=resolution)
    norm.init_fields()
    source(norm)

    flux_inc = meep_ext.add_flux_plane(norm, fcen, df, nfreq, [0, 0, 0], [W, W, 0])

    norm.run(
        until_after_sources=meep.stop_when_fields_decayed(0.5 * um, decay, pt=meep.Vector3(0, 0, 0), decay_by=1e-3)
    )

    return {
        "frequency": np.array(meep.get_flux_freqs(flux_inc)),
        "area": (W) ** 2,
        "incident": np.asarray(meep.get_fluxes(flux_inc)),
    }


@job.cache
def force_sim():
    """Perform scattering simulation."""
    sim = meep.Simulation(
        cell_size=cell, boundary_layers=[pml], geometry=geometry, default_material=medium, resolution=resolution
    )
    sim.init_fields()
    source(sim)

    L = 2**0.5 * W + 2 * particle_monitor_gap
    p = meep.Vector3(sep / 2, 0, 0)
    Fx_mon = meep_ext.add_force_box(sim, fcen, df, nfreq, p, [L, L, L], meep.X)
    Fy_mon = meep_ext.add_force_box(sim, fcen, df, nfreq, p, [L, L, L], meep.Y)
    Fz_mon = meep_ext.add_force_box(sim, fcen, df, nfreq, p, [L, L, L], meep.Z)

    sim.run(
        until_after_sources=meep.stop_when_fields_decayed(
            0.5 * um, decay, pt=meep.Vector3(0, 0, monitor_size[2] / 2), decay_by=1e-5
        )
    )

    Fx = np.array(meep.get_forces(Fx_mon))
    Fy = np.array(meep.get_forces(Fy_mon))
    Fz = np.array(meep.get_forces(Fz_mon))
    frequency = np.array(meep.get_force_freqs(Fx_mon))
    return dict(frequency=frequency, Fx=Fx, Fy=Fy, Fz=Fz)


@job.cache
def force_gmmt():
    wavelengths = np.linspace(400 * nm, 1000 * nm, 100)
    # eps = meep_ext.get_eps(material)(wavelengths)
    # Au = miepy.data_material(wavelengths, eps)
    material = miepy.constant_material(3.5**2)

    particles = [miepy.cube([-sep / 2, 0, 0], 200 * nm, material), miepy.cube([sep / 2, 0, 0], 200 * nm, material)]

    F1 = np.zeros((3,) + wavelengths.shape, dtype=float)
    F2 = np.zeros((3,) + wavelengths.shape, dtype=float)

    for i, wavelength in enumerate(pbar(wavelengths)):
        c = miepy.cluster(particles=particles, source=miepy.sources.plane_wave([1, 0]), wavelength=wavelength, lmax=4)

        q = miepy.quaternion.from_spherical_coords(0, 0)
        c.update(orientation=[q, q])
        F1[:, i] = c.force_on_particle(1)

        q = miepy.quaternion.from_spherical_coords(0, np.pi / 4)
        c.update(orientation=[q, q])
        F2[:, i] = c.force_on_particle(1)

    return dict(wavelengths=wavelengths, F1=F1, F2=F2)


@job.plots
def vis():
    norm = job.load(force_norm)
    scat = job.load(force_sim)
    gmmt = job.load(force_gmmt)

    fig, axes = plt.subplots(ncols=3, figsize=(16, 5))
    for j, c in enumerate(["x", "y", "z"]):
        ax = axes[j]
        key = f"F{c}"
        ax.plot(
            (1 / nm) / norm.frequency,
            scat[key] / norm.incident * norm.area * constants.epsilon_0 / 2 * 1e25,
            "o",
            color=f"C{j}",
            label=f"F{c} (FDTD)",
        )
        ax.plot(gmmt.wavelengths / nm, gmmt.F1[j] * 1e25, "-", color=f"C{j}", label=f"F{c}, 1 (GMMT)")
        ax.plot(gmmt.wavelengths / nm, gmmt.F2[j] * 1e25, "--", color=f"C{j}", label=f"F{c}, 2 (GMMT)")
        ax.set(xlabel="wavelength (nm)", ylabel="force")
        ax.set_title(f"F{c}", weight="bold")
        ax.legend()

    plt.show()


job.run()
