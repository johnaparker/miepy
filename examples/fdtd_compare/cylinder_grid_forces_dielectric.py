import matplotlib.pyplot as plt
import meep
import meep_ext
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from numpipe import pbar, scheduler
from scipy import constants

import miepy

job = scheduler()
nm = 1e-9
um = 1e-6

### Parameters
radius = 75 * nm
height = 110 * nm
# material = meep_ext.material.Au()
material = meep.Medium(index=3)

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
resolution = 1 / (4 * nm)
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
def force_norm():
    """Perform normalization simulation."""
    norm = meep.Simulation(cell_size=cell, boundary_layers=[pml], geometry=[], resolution=resolution)
    norm.init_fields()
    source(norm)

    flux_inc = meep_ext.add_flux_plane(norm, fcen, df, nfreq, [0, 0, 0], [2 * radius, 2 * radius, 0])

    norm.run(
        until_after_sources=meep.stop_when_fields_decayed(0.5 * um, decay, pt=meep.Vector3(0, 0, 0), decay_by=1e-5)
    )

    return {
        "frequency": np.array(meep.get_flux_freqs(flux_inc)),
        "area": (2 * radius) ** 2,
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

    L = 2 * radius + 2 * particle_monitor_gap
    forces = {}
    for i in range(9):
        p = meep.Vector3(x[i], y[i], z[i])
        key = "N{i}_F{c}"
        forces[key.format(i=i, c="x")] = meep_ext.add_force_box(sim, fcen, df, nfreq, p, [L, L, L], meep.X)
        forces[key.format(i=i, c="y")] = meep_ext.add_force_box(sim, fcen, df, nfreq, p, [L, L, L], meep.Y)
        forces[key.format(i=i, c="z")] = meep_ext.add_force_box(sim, fcen, df, nfreq, p, [L, L, L], meep.Z)

    sim.run(
        until_after_sources=meep.stop_when_fields_decayed(
            0.5 * um, decay, pt=meep.Vector3(0, 0, monitor_size[2] / 2), decay_by=1e-8
        )
    )

    ret = {}
    for i in range(9):
        for c in ["x", "y", "z"]:
            key = f"N{i}_F{c}"
            ret[key] = np.array(meep.get_forces(forces[key]))

    ret["frequency"] = np.array(meep.get_force_freqs(forces["N1_Fx"]))
    return ret


@job.cache
def force_gmmt():
    wavelengths = np.linspace(400 * nm, 1000 * nm, 100)
    eps = meep_ext.get_eps(material)(wavelengths)
    Au = miepy.data_material(wavelengths, eps)
    # Au = miepy.constant_material(3.5**2)

    particles = []
    for i in range(9):
        orientation = miepy.quaternion.from_spherical_coords(theta[i], phi[i])
        particles.append(
            miepy.cylinder(
                position=[x[i], y[i], z[i]], radius=radius, height=height, material=Au, orientation=orientation
            )
        )

    F = np.zeros([3, 9, len(wavelengths)], dtype=float)

    for i, wavelength in enumerate(pbar(wavelengths)):
        sol = miepy.cluster(particles=particles, source=miepy.sources.plane_wave([1, 0]), wavelength=wavelength, lmax=4)

        F[..., i] = sol.force()
    miepy.visualize(sol)

    return dict(wavelengths=wavelengths, force=F)


@job.plots
def vis():
    norm = job.load(force_norm)
    scat = job.load(force_sim)
    gmmt = job.load(force_gmmt)

    with PdfPages("forces.pdf") as pdf:
        for i in range(9):
            fig, axes = plt.subplots(ncols=3, figsize=(16, 5))
            for j, c in enumerate(["x", "y", "z"]):
                ax = axes[j]
                key = f"N{i}_F{c}"
                ax.plot(
                    (1 / nm) / norm.frequency,
                    scat[key] / norm.incident * norm.area * constants.epsilon_0 / 2 * 1e25,
                    "o",
                    color=f"C{j}",
                    label=f"F{c} (FDTD)",
                )
                ax.plot(gmmt.wavelengths / nm, gmmt.force[j, i] * 1e25, "-", color=f"C{j}", label=f"F{c} (GMMT)")
                ax.set(xlabel="wavelength (nm)", ylabel="force")
                ax.set_title(f"F{c}", weight="bold")
                ax.legend()
            fig.suptitle(f"Particle {i + 1}", fontsize=18)
            pdf.savefig(fig)

    # plt.show()


job.run()
